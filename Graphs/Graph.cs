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
        public List<Node> Nodes { get; set; } = null!;

        /// <summary>All edges (branches) in the graph, each connecting two nodes.</summary>
        [JsonIgnore]
        public List<Edge> Edges { get; set; } = null!;

        /// <summary>Grasshopper tree mapping to branch/item indices for graph elements.</summary>
        [JsonIgnore]
        public GH_Structure<GH_Number> IndicesTree { get; set; } = null!;

        /// <summary>Grasshopper tree representing adjacency structure.</summary>
        [JsonIgnore]
        public GH_Structure<GH_Number> AdjacencyTree { get; set; } = null!;

        /// <summary>Maps each edge index to (branchIndex, itemIndex) in the original input tree.</summary>
        [JsonIgnore]
        public List<(int branchIndex, int itemIndex)>? EdgeInputMap { get; set; }

        /// <summary>Output edge tree matching Grasshopper structure for downstream components.</summary>
        [JsonIgnore]
        public GH_Structure<Edge> OutputEdgeTree { get; set; } = null!;

        /// <summary>Parameterless constructor; initialize properties before use.</summary>
        public Graph() { }

        /// <summary>Builds a graph from a flat list of curves using the given merge tolerance.</summary>
        /// <param name="_InputCurves">Curves to convert into edges; endpoints become nodes.</param>
        /// <param name="_Tolerance">Distance tolerance for merging coincident vertices.</param>
        public Graph(List<GH_Curve> _InputCurves, double _Tolerance)
        {
            ConstructGraph(_InputCurves, _Tolerance);
        }

        /// <summary>Builds a graph from a Grasshopper tree of curves.</summary>
        /// <param name="inputTree">Tree of curves; structure is preserved in output.</param>
        /// <param name="tol">Distance tolerance for merging coincident vertices.</param>
        public Graph(GH_Structure<GH_Curve> inputTree, double tol)
        {
            ConstructGraphFromTree(inputTree, tol);
        }

        /// <summary>Copy constructor; shallow copy of node and edge lists.</summary>
        /// <param name="other">Graph to copy from.</param>
        public Graph(Graph other)
        {
            Tolerance = other.Tolerance;
            Nodes = new List<Node>(other.Nodes);
            Edges = new List<Edge>(other.Edges);
            IndicesTree = other.IndicesTree;
            AdjacencyTree = other.AdjacencyTree;
            EdgeInputMap = other.EdgeInputMap != null ? new List<(int, int)>(other.EdgeInputMap) : null;
            OutputEdgeTree = other.OutputEdgeTree;
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
                    if (branch[itemIdx] is GH_Curve curve)
                    {
                        allCurves.Add(curve);
                        EdgeInputMap.Add((branchIdx, itemIdx));
                    }
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
