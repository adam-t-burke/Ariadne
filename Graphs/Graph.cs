using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;
using Rhino.Input.Custom;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Text.Json;
using System.Text.Json.Serialization;
using System.Threading;
using System.Threading.Tasks;

namespace Ariadne.Graphs
{
    /// <summary>
    /// Topology and geometry of a network: nodes (points) and edges (curves) built from input curves
    /// with optional tolerance-based vertex merging. Used by <see cref="FDM.FDM_Network"/>.
    /// </summary>
    public class Graph
    {
        /// <summary>Distance tolerance for merging coincident or near-coincident vertices when building the graph.</summary>
        [JsonIgnore]
        public double Tolerance { get; set; }

        /// <summary>All nodes (vertices) in the graph.</summary>
        [JsonIgnore]
        public List<Node> Nodes { get; set; }

        /// <summary>All edges (branches) in the graph, each connecting two nodes.</summary>
        [JsonIgnore]
        public List<Edge> Edges { get; set; }

        /// <summary>Grasshopper tree mapping to branch/item indices for graph elements.</summary>
        [JsonIgnore]
        public GH_Structure<GH_Number> IndicesTree { get; set; }

        /// <summary>Grasshopper tree representing adjacency structure.</summary>
        [JsonIgnore]
        public GH_Structure<GH_Number> AdjacencyTree { get; set; }

        /// <summary>Maps each edge index to (branchIndex, itemIndex) in the original input tree.</summary>
        [JsonIgnore]
        public List<(int branchIndex, int itemIndex)> EdgeInputMap { get; set; }

        /// <summary>Output edge tree matching Grasshopper structure for downstream components.</summary>
        [JsonIgnore]
        public GH_Structure<Edge> OutputEdgeTree { get; set; }

        /// <summary>Parameterless constructor; initialize properties before use.</summary>
        public Graph()
        {
            Nodes = new List<Node>();
            Edges = new List<Edge>();
            IndicesTree = new GH_Structure<GH_Number>();
            AdjacencyTree = new GH_Structure<GH_Number>();
            EdgeInputMap = new List<(int branchIndex, int itemIndex)>();
            OutputEdgeTree = new GH_Structure<Edge>();
        }

        /// <summary>Builds a graph from a flat list of curves using the given merge tolerance.</summary>
        /// <param name="_InputCurves">Curves to convert into edges; endpoints become nodes.</param>
        /// <param name="_Tolerance">Distance tolerance for merging coincident vertices.</param>
        public Graph(List<GH_Curve> _InputCurves, double _Tolerance)
        {
            Nodes = new List<Node>();
            Edges = new List<Edge>();
            IndicesTree = new GH_Structure<GH_Number>();
            AdjacencyTree = new GH_Structure<GH_Number>();
            EdgeInputMap = new List<(int branchIndex, int itemIndex)>();
            OutputEdgeTree = new GH_Structure<Edge>();
            ConstructGraph(_InputCurves, _Tolerance);
        }

        /// <summary>Builds a graph from a Grasshopper tree of curves.</summary>
        /// <param name="inputTree">Tree of curves; structure is preserved in output.</param>
        /// <param name="tol">Distance tolerance for merging coincident vertices.</param>
        public Graph(GH_Structure<GH_Curve> inputTree, double tol)
        {
            Nodes = new List<Node>();
            Edges = new List<Edge>();
            IndicesTree = new GH_Structure<GH_Number>();
            AdjacencyTree = new GH_Structure<GH_Number>();
            EdgeInputMap = new List<(int branchIndex, int itemIndex)>();
            OutputEdgeTree = new GH_Structure<Edge>();
            ConstructGraphFromTree(inputTree, tol);
        }

        /// <summary>Copy constructor; shallow copy of node and edge lists.</summary>
        /// <param name="other">Graph to copy from.</param>
        public Graph(Graph other)
        {
            Tolerance = other.Tolerance;
            Nodes = new List<Node>(other.Nodes);
            Edges = new List<Edge>(other.Edges);
            IndicesTree = other.IndicesTree ?? new GH_Structure<GH_Number>();
            AdjacencyTree = other.AdjacencyTree ?? new GH_Structure<GH_Number>();
            EdgeInputMap = other.EdgeInputMap != null ? new List<(int, int)>(other.EdgeInputMap) : new List<(int, int)>();
            OutputEdgeTree = other.OutputEdgeTree ?? new GH_Structure<Edge>();
        }

        // Internal construction mode (kept private per request)
        private enum GraphConstructionMode
        {
            SequentialHashed,
            ParallelFast
        }

        private const int ParallelThreshold = 256; // default threshold for enabling parallel build

        private void ConstructGraph(List<GH_Curve> curves, double tol)
        {
            Tolerance = tol;
            if (curves == null || curves.Count == 0)
            {
                Nodes = new List<Node>();
                Edges = new List<Edge>();
                return;
            }

            // Decide mode (fast parallel if large enough)
            var mode =
                (tol > 0 &&
                 curves.Count >= ParallelThreshold &&
                 Environment.ProcessorCount > 1)
                 ? GraphConstructionMode.ParallelFast
                 : GraphConstructionMode.SequentialHashed;

            if (mode == GraphConstructionMode.ParallelFast)
                ConstructGraphParallelFast(curves, tol);
            else
                ConstructGraphSequentialHashed(curves, tol);

            // Assign stable indices (current order accepted as "fast mode")
            for (int i = 0; i < Nodes.Count; i++)
                Nodes[i].Index = i;
        }

        #region Sequential (Hashed) Construction
        private void ConstructGraphSequentialHashed(List<GH_Curve> curves, double tol)
        {
            double cell = tol;
            double tol2 = tol * tol;

            // Spatial hash: cell key -> list of node indices
            Dictionary<(long, long, long), List<int>> cells = new();
            Nodes = new List<Node>(curves.Count * 2);
            Edges = new List<Edge>(curves.Count);

            static (long, long, long) Key(Point3d p, double c) =>
                ((long)Math.Floor(p.X / c),
                 (long)Math.Floor(p.Y / c),
                 (long)Math.Floor(p.Z / c));

            static IEnumerable<(long, long, long)> NeighborKeys((long, long, long) k)
            {
                for (long dx = -1; dx <= 1; dx++)
                    for (long dy = -1; dy <= 1; dy++)
                        for (long dz = -1; dz <= 1; dz++)
                            yield return (k.Item1 + dx, k.Item2 + dy, k.Item3 + dz);
            }

            int GetOrCreate(Point3d p)
            {
                // With non-positive tolerance just treat each unique coordinate separately
                if (Tolerance <= 0)
                {
                    // Linear search fallback (very rare path)
                    for (int i = 0; i < Nodes.Count; i++)
                    {
                        if (Nodes[i].Value.DistanceToSquared(p) <= 0)
                            return i;
                    }
                    var n = new Node { Value = p };
                    Nodes.Add(n);
                    return Nodes.Count - 1;
                }

                var key = Key(p, cell);
                int found = -1;

                // Search existing neighbors
                foreach (var nk in NeighborKeys(key))
                {
                    if (!cells.TryGetValue(nk, out var list)) continue;
                    for (int i = 0; i < list.Count; i++)
                    {
                        var candidate = Nodes[list[i]];
                        double dx = candidate.Value.X - p.X;
                        double dy = candidate.Value.Y - p.Y;
                        double dz = candidate.Value.Z - p.Z;
                        if (dx * dx + dy * dy + dz * dz <= tol2)
                        {
                            found = list[i];
                            break;
                        }
                    }
                    if (found >= 0) break;
                }

                if (found >= 0) return found;

                // Create new node
                var node = new Node { Value = p };
                Nodes.Add(node);
                int id = Nodes.Count - 1;

                if (!cells.TryGetValue(key, out var bucket))
                {
                    bucket = new List<int>(4);
                    cells[key] = bucket;
                }
                bucket.Add(id);
                return id;
            }

            // Build arrays of node ids first
            int edgeCount = curves.Count;
            int[] startIds = new int[edgeCount];
            int[] endIds = new int[edgeCount];

            for (int i = 0; i < edgeCount; i++)
            {
                var c = curves[i];
                var s = c.Value.PointAtStart;
                var e = c.Value.PointAtEnd;
                int u = GetOrCreate(s);
                int v = GetOrCreate(e);
                startIds[i] = u;
                endIds[i] = v;
            }

            // Create edges & neighbor lists
            for (int i = 0; i < edgeCount; i++)
            {
                var ghCurve = curves[i];
                var edge = new Edge
                {
                    Start = Nodes[startIds[i]],
                    End = Nodes[endIds[i]],
                    Value = ghCurve.Value
                };
                if (ghCurve.IsReferencedGeometry)
                    edge.ReferenceID = ghCurve.ReferenceID;

                // Preserve original behavior (including self-loop duplication)
                edge.Start.Neighbors.Add(edge.End);
                edge.End.Neighbors.Add(edge.Start);

                Edges.Add(edge);
            }
        }
        #endregion

        #region Parallel Fast Construction
        // Optimized parallel construction with region-based locking to prevent duplicates
        private void ConstructGraphParallelFast(List<GH_Curve> curves, double tol)
        {
            double cell = tol;
            double tol2 = tol * tol;

            Nodes = new List<Node>();
            Edges = new List<Edge>(curves.Count);

            // Region-based locking: lock groups of adjacent cells to reduce contention
            var regionLocks = new ConcurrentDictionary<(long, long, long), object>();
            var nodesById = new ConcurrentDictionary<int, Node>();
            var cellBuckets = new ConcurrentDictionary<(long, long, long), List<int>>();
            int nextId = 0;

            static (long, long, long) Key(Point3d p, double c) =>
                ((long)Math.Floor(p.X / c),
                 (long)Math.Floor(p.Y / c),
                 (long)Math.Floor(p.Z / c));

            static (long, long, long) RegionKey((long, long, long) cellKey) =>
                (cellKey.Item1 >> 2, cellKey.Item2 >> 2, cellKey.Item3 >> 2); // 4x4x4 cell regions

            static IEnumerable<(long, long, long)> NeighborKeys((long, long, long) k)
            {
                for (long dx = -1; dx <= 1; dx++)
                    for (long dy = -1; dy <= 1; dy++)
                        for (long dz = -1; dz <= 1; dz++)
                            yield return (k.Item1 + dx, k.Item2 + dy, k.Item3 + dz);
            }

            int GetOrCreate(Point3d p)
            {
                if (Tolerance <= 0)
                {
                    // Exact matching fallback
                    foreach (var kv in nodesById)
                    {
                        if (kv.Value.Value.DistanceToSquared(p) <= 0) return kv.Key;
                    }
                    int id = Interlocked.Increment(ref nextId) - 1;
                    nodesById[id] = new Node { Value = p };
                    return id;
                }

                var cellKey = Key(p, cell);
                var regionKey = RegionKey(cellKey);
                var regionLock = regionLocks.GetOrAdd(regionKey, _ => new object());

                lock (regionLock)
                {
                    // Search for existing node within tolerance
                    foreach (var nk in NeighborKeys(cellKey))
                    {
                        if (!cellBuckets.TryGetValue(nk, out var bucket)) continue;
                        foreach (var nid in bucket)
                        {
                            if (!nodesById.TryGetValue(nid, out var node)) continue;
                            double dx = node.Value.X - p.X;
                            double dy = node.Value.Y - p.Y;
                            double dz = node.Value.Z - p.Z;
                            if (dx * dx + dy * dy + dz * dz <= tol2)
                                return nid;
                        }
                    }

                    // Create new node
                    int newId = Interlocked.Increment(ref nextId) - 1;
                    var newNode = new Node { Value = p };
                    nodesById[newId] = newNode;

                    var cellBucket = cellBuckets.GetOrAdd(cellKey, _ => new List<int>());
                    cellBucket.Add(newId);

                    return newId;
                }
            }

            int edgeCount = curves.Count;
            int[] startIds = new int[edgeCount];
            int[] endIds = new int[edgeCount];

            // Parallel endpoint processing with region locking
            Parallel.For(0, edgeCount, i =>
            {
                var c = curves[i];
                var s = c.Value.PointAtStart;
                var e = c.Value.PointAtEnd;
                int u = GetOrCreate(s);
                int v = GetOrCreate(e);
                startIds[i] = u;
                endIds[i] = v;
            });

            // Build final node list
            int nodeCount = nextId;
            var nodeArray = new Node[nodeCount];
            foreach (var kv in nodesById)
                nodeArray[kv.Key] = kv.Value;

            Nodes = nodeArray.ToList();

            // Sequential edge creation and neighbor linking
            for (int i = 0; i < edgeCount; i++)
            {
                var ghCurve = curves[i];
                var edge = new Edge
                {
                    Start = nodeArray[startIds[i]],
                    End = nodeArray[endIds[i]],
                    Value = ghCurve.Value
                };
                if (ghCurve.IsReferencedGeometry)
                    edge.ReferenceID = ghCurve.ReferenceID;

                edge.Start.Neighbors.Add(edge.End);
                edge.End.Neighbors.Add(edge.Start);

                Edges.Add(edge);
            }
        }
        #endregion

        // Update EdgeIndicesToTree for nested structure
        public void EdgeIndicesToTree()
        {
            IndicesTree = new GH_Structure<GH_Number>();

            if (EdgeInputMap == null || EdgeInputMap.Count != Edges.Count)
            {
                // Fallback to existing behavior if no mapping exists
                int counter = 0;
                foreach (var edge in Edges)
                {
                    IndicesTree.Append(new GH_Number(edge.Start.Index), new GH_Path(counter));
                    IndicesTree.Append(new GH_Number(edge.End.Index), new GH_Path(counter));
                    counter++;
                }
                return;
            }

            for (int edgeIdx = 0; edgeIdx < Edges.Count; edgeIdx++)
            {
                var (branchIdx, itemIdx) = EdgeInputMap[edgeIdx];
                var path = new GH_Path(branchIdx, itemIdx);

                IndicesTree.Append(new GH_Number(Edges[edgeIdx].Start.Index), path);
                IndicesTree.Append(new GH_Number(Edges[edgeIdx].End.Index), path);
            }
        }


        public void AdjacencyListToTree()
        {
            AdjacencyTree = new GH_Structure<GH_Number>();
            int counter = 0;
            foreach (var node in Nodes)
            {
                foreach (var nbr in node.Neighbors)
                {
                    AdjacencyTree.Append(new GH_Number(nbr.Index), new GH_Path(counter));
                }
                counter++;
            }
        }

        public Graph Duplicate() => new Graph(this);

        public int Ne => Edges.Count;
        public int Nn => Nodes.Count;

        public void UpdateNodeIndices()
        {
            for (int i = 0; i < Nodes.Count; i++)
            {
                Nodes[i].Index = i;
            }
        }

        private void ConstructGraphFromTree(GH_Structure<GH_Curve> inputTree, double tol)
        {
            // Flatten curves while preserving mapping
            var allCurves = new List<GH_Curve>();
            EdgeInputMap = new List<(int branchIndex, int itemIndex)>();
            
            for (int branchIdx = 0; branchIdx < inputTree.PathCount; branchIdx++)
            {
                var path = inputTree.Paths[branchIdx];
                var branch = inputTree.get_Branch(path);
                
                for (int itemIdx = 0; itemIdx < branch.Count; itemIdx++)
                {
                    allCurves.Add((GH_Curve)branch[itemIdx]!);
                    EdgeInputMap.Add((branchIdx, itemIdx));
                }
            }
            
            // Use existing construction logic
            ConstructGraph(allCurves, tol);
        }

        // New method to rebuild tree structure for output
        public void BuildOutputEdgeTree()
        {

            if (EdgeInputMap == null || EdgeInputMap.Count != Edges.Count)
            {
                // Fallback to simple list structure if no mapping exists
                OutputEdgeTree = new GH_Structure<Edge>();
                for (int i = 0; i < Edges.Count; i++)
                {
                    OutputEdgeTree.Append(Edges[i], new GH_Path(0));
                }
                return;
            }

            OutputEdgeTree = new GH_Structure<Edge>();
    
            for (int edgeIdx = 0; edgeIdx < Edges.Count; edgeIdx++)
            {
                var (branchIdx, itemIdx) = EdgeInputMap[edgeIdx];
                var path = new GH_Path(branchIdx);
                OutputEdgeTree.Append(Edges[edgeIdx], path);
            }
        }

        // ── Face-finding from 3D graph topology ──────────────────

        /// <summary>
        /// Result of the face-finding algorithm.
        /// </summary>
        public class FaceFindResult
        {
            /// <summary>Bounded face cycles as ordered lists of node indices.</summary>
            public List<List<int>> Faces { get; init; } = [];
            /// <summary>Euler characteristic (V - E + F, including exterior face).</summary>
            public int EulerCharacteristic { get; init; }
            /// <summary>True if the graph appears to be a valid 2-manifold.</summary>
            public bool IsManifold { get; init; }
            /// <summary>Edges shared by more than 2 faces (non-manifold).</summary>
            public List<(int, int)> NonManifoldEdges { get; init; } = [];
            /// <summary>Warning messages for the user.</summary>
            public List<string> Warnings { get; init; } = [];
        }

        /// <summary>
        /// Find all bounded face cycles in a topologically planar graph embedded in 3D.
        /// Uses local tangent frames at each vertex for angular sorting, then walks
        /// half-edges to extract faces. No global projection is needed.
        /// </summary>
        public static FaceFindResult FindFaces(Graph graph)
        {
            int V = graph.Nn;
            int E = graph.Ne;
            var warnings = new List<string>();

            // Build adjacency: for each node, list of (neighbor_index, edge_index)
            var adj = new List<List<(int neighbor, int edgeIdx)>>(V);
            for (int i = 0; i < V; i++)
                adj.Add(new List<(int, int)>());

            for (int e = 0; e < E; e++)
            {
                int s = graph.Edges[e].Start.Index;
                int t = graph.Edges[e].End.Index;
                adj[s].Add((t, e));
                adj[t].Add((s, e));
            }

            // Compute local normal at each vertex from cross products of successive edge pairs
            var localNormals = new Vector3d[V];
            for (int v = 0; v < V; v++)
            {
                var neighbors = adj[v];
                var pos = graph.Nodes[v].Value;
                var normal = Vector3d.Zero;

                if (neighbors.Count >= 2)
                {
                    for (int i = 0; i < neighbors.Count; i++)
                    {
                        int j = (i + 1) % neighbors.Count;
                        var u1 = graph.Nodes[neighbors[i].neighbor].Value - pos;
                        var u2 = graph.Nodes[neighbors[j].neighbor].Value - pos;
                        normal += Vector3d.CrossProduct(u1, u2);
                    }
                }

                if (normal.Length < 1e-12)
                {
                    // Fallback: pick an arbitrary normal perpendicular to the first edge
                    if (neighbors.Count > 0)
                    {
                        var edge = graph.Nodes[neighbors[0].neighbor].Value - pos;
                        normal = PerpVector(edge);
                    }
                    else
                    {
                        normal = Vector3d.ZAxis;
                    }
                }

                normal.Unitize();
                localNormals[v] = normal;
            }

            // Build local tangent frame at each vertex and sort edges by angle
            // sortedAdj[v] = list of neighbor indices sorted by angle in local tangent plane
            var sortedAdj = new List<int>[V];
            for (int v = 0; v < V; v++)
            {
                var neighbors = adj[v];
                var pos = graph.Nodes[v].Value;
                var n = localNormals[v];

                // Build orthonormal tangent frame
                Vector3d t1, t2;
                if (neighbors.Count > 0)
                {
                    var firstEdge = graph.Nodes[neighbors[0].neighbor].Value - pos;
                    t1 = firstEdge - (firstEdge * n) * n;
                    if (t1.Length < 1e-12)
                        t1 = PerpVector(n);
                    t1.Unitize();
                }
                else
                {
                    t1 = PerpVector(n);
                    t1.Unitize();
                }
                t2 = Vector3d.CrossProduct(n, t1);
                t2.Unitize();

                // Sort neighbors by angle in (t1, t2) plane
                var angles = new List<(double angle, int neighbor)>(neighbors.Count);
                foreach (var (nb, _) in neighbors)
                {
                    var dir = graph.Nodes[nb].Value - pos;
                    double x = dir * t1;
                    double y = dir * t2;
                    angles.Add((Math.Atan2(y, x), nb));
                }
                angles.Sort((a, b) => a.angle.CompareTo(b.angle));
                sortedAdj[v] = angles.Select(a => a.neighbor).ToList();
            }

            // Walk half-edges to find face cycles.
            // A half-edge is (from, to). For each unvisited half-edge, walk the face.
            var visited = new HashSet<(int, int)>();
            var allFaces = new List<List<int>>();

            for (int v = 0; v < V; v++)
            {
                foreach (int nb in sortedAdj[v])
                {
                    if (visited.Contains((v, nb)))
                        continue;

                    var face = WalkFace(v, nb, sortedAdj, visited);
                    if (face != null)
                        allFaces.Add(face);
                }
            }

            // Count half-edge usage per undirected edge to detect non-manifold
            var edgeUsage = new Dictionary<(int, int), int>();
            foreach (var face in allFaces)
            {
                for (int i = 0; i < face.Count; i++)
                {
                    int a = face[i];
                    int b = face[(i + 1) % face.Count];
                    var key = a < b ? (a, b) : (b, a);
                    edgeUsage.TryGetValue(key, out int count);
                    edgeUsage[key] = count + 1;
                }
            }

            var nonManifold = new List<(int, int)>();
            foreach (var (key, count) in edgeUsage)
            {
                if (count > 2)
                    nonManifold.Add(key);
            }

            if (nonManifold.Count > 0)
                warnings.Add($"{nonManifold.Count} non-manifold edge(s) found (shared by >2 faces).");

            // Identify and discard the exterior face (largest face by vertex count)
            int totalFaces = allFaces.Count;
            if (totalFaces > 0)
            {
                int maxIdx = 0;
                int maxLen = allFaces[0].Count;
                for (int i = 1; i < totalFaces; i++)
                {
                    if (allFaces[i].Count > maxLen)
                    {
                        maxLen = allFaces[i].Count;
                        maxIdx = i;
                    }
                }
                allFaces.RemoveAt(maxIdx);
            }

            // Euler formula check: V - E + F = chi
            // totalFaces includes the exterior face we just removed, so F_total = allFaces.Count + 1
            int F_total = allFaces.Count + 1;
            int chi = V - E + F_total;
            bool isManifold = nonManifold.Count == 0 && (chi == 1 || chi == 2);

            if (chi != 1 && chi != 2)
                warnings.Add($"Euler characteristic χ = {chi} (expected 1 for open surface, 2 for closed). Graph may be non-planar or have higher genus.");

            return new FaceFindResult
            {
                Faces = allFaces,
                EulerCharacteristic = chi,
                IsManifold = isManifold,
                NonManifoldEdges = nonManifold,
                Warnings = warnings,
            };
        }

        /// <summary>
        /// Walk a face cycle starting from half-edge (start -> current).
        /// At each vertex, find the predecessor of the reverse half-edge in the
        /// sorted adjacency (turn right relative to local normal).
        /// </summary>
        private static List<int>? WalkFace(
            int start, int next,
            List<int>[] sortedAdj,
            HashSet<(int, int)> visited)
        {
            var face = new List<int>();
            int from = start;
            int to = next;
            int maxSteps = 0;
            int totalEdges = visited.Count + sortedAdj.Sum(a => a.Count);

            while (true)
            {
                if (visited.Contains((from, to)))
                    break;

                visited.Add((from, to));
                face.Add(from);

                // At vertex 'to', find reverse half-edge (to -> from) in sorted adjacency
                var adjList = sortedAdj[to];
                int reverseIdx = adjList.IndexOf(from);
                if (reverseIdx < 0)
                    return null; // broken topology

                // Next half-edge: predecessor in cyclic order (turn right)
                int prevIdx = (reverseIdx - 1 + adjList.Count) % adjList.Count;
                int nextNode = adjList[prevIdx];

                from = to;
                to = nextNode;

                maxSteps++;
                if (maxSteps > totalEdges)
                    return null; // infinite loop guard

                if (from == start && to == next)
                    break;
            }

            return face.Count >= 3 ? face : null;
        }

        /// <summary>
        /// Returns an arbitrary vector perpendicular to the given vector.
        /// </summary>
        private static Vector3d PerpVector(Vector3d v)
        {
            if (Math.Abs(v.X) < 0.9)
                return Vector3d.CrossProduct(v, Vector3d.XAxis);
            return Vector3d.CrossProduct(v, Vector3d.YAxis);
        }
    }

    /// <summary>
    /// A vertex in a <see cref="Graph"/>, with position, neighbor list, and optional anchor flag.
    /// </summary>
    public class Node : GH_Point
    {
        /// <summary>When true, this node is treated as a fixed support (anchor).</summary>
        [JsonIgnore]
        public bool Anchor { get; set; }

        /// <summary>Nodes connected to this node by an edge.</summary>
        [JsonIgnore]
        public List<Node> Neighbors { get; set; }

        /// <summary>Stable index of this node in the graph's Nodes list (0-based).</summary>
        [JsonIgnore]
        public int Index { get; set; }

        /// <summary>Default constructor; non-anchor, empty neighbors, index -1.</summary>
        public Node()
        {
            Anchor = false;
            Neighbors = new List<Node>();
            Value = new Point3d();
            Index = -1;
        }

        /// <summary>Copy constructor.</summary>
        /// <param name="other">Node to copy from.</param>
        public Node(Node other)
        {
            Anchor = other.Anchor;
            Neighbors = other.Neighbors;
            Value = other.Value;
            Index = other.Index;
        }

        /// <summary>Constructs a node with the given anchor flag and neighbor list.</summary>
        /// <param name="anchor">Whether this node is a fixed support.</param>
        /// <param name="neighbors">Initial neighbor list.</param>
        public Node(bool anchor, List<Node> neighbors)
        {
            Anchor = anchor;
            Neighbors = neighbors;
            Index = -1;
        }
    }

    /// <summary>
    /// A branch in a <see cref="Graph"/> connecting two <see cref="Node"/>s, with optional force density Q and geometry.
    /// </summary>
    public class Edge : GH_Curve
    {
        /// <summary>Start node of this edge.</summary>
        [JsonIgnore]
        public Node Start { get; set; }

        /// <summary>End node of this edge.</summary>
        [JsonIgnore]
        public Node End { get; set; }

        /// <summary>Force density (q) for this member; used by the FDM solver.</summary>
        [JsonIgnore]
        public double Q { get; set; }

        /// <summary>Optional reference ID for Grasshopper geometry linking.</summary>
        [JsonIgnore]
        public new Guid ReferenceID { get; set; }

        /// <summary>Default constructor; empty start/end nodes, Q = 0.</summary>
        public Edge()
        {
            Start = new Node();
            End = new Node();
            Q = 0;
            ReferenceID = Guid.Empty;
        }

        /// <summary>Constructs an edge from two nodes, force density, and curve geometry.</summary>
        /// <param name="start">Start node.</param>
        /// <param name="end">End node.</param>
        /// <param name="q">Force density for this member.</param>
        /// <param name="curve">Curve geometry (length and shape).</param>
        public Edge(Node start, Node end, double q, GH_Curve curve)
        {
            Start = start;
            End = end;
            Q = q;
            Value = curve.Value;
            ReferenceID = Guid.Empty;
        }

        /// <summary>Copy constructor.</summary>
        /// <param name="other">Edge to copy from.</param>
        public Edge(Edge other)
        {
            Start = other.Start;
            End = other.End;
            Q = other.Q;
            Value = other.Value;
            ReferenceID = other.ReferenceID;
        }
    }
}
