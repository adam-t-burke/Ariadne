using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.Json;
using System.Text.Json.Serialization;
using System.Collections.Concurrent;
using System.Threading;
using System.Threading.Tasks;
using Rhino.Geometry;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;

namespace Ariadne.Graphs
{
    internal class Graph
    {
        [JsonIgnore]
        public double Tolerance { get; set; }

        [JsonIgnore]
        public List<Node> Nodes { get; set; }

        [JsonIgnore]
        public List<Edge> Edges { get; set; }

        [JsonIgnore]
        public GH_Structure<GH_Number> IndicesTree { get; set; }

        [JsonIgnore]
        public GH_Structure<GH_Number> AdjacencyTree { get; set; }

        public Graph() { }

        public Graph(List<GH_Curve> _InputCurves, double _Tolerance)
        {
            ConstructGraph(_InputCurves, _Tolerance);
        }

        public Graph(Graph other)
        {
            Tolerance = other.Tolerance;
            Nodes = other.Nodes;
            Edges = other.Edges;
            IndicesTree = other.IndicesTree;
            AdjacencyTree = other.AdjacencyTree;
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

        public void EdgeIndicesToTree()
        {
            IndicesTree = new GH_Structure<GH_Number>();
            int counter = 0;
            foreach (var edge in Edges)
            {
                // Using pre-assigned Node indices (avoids O(N) Node lookup)
                IndicesTree.Append(new GH_Number(edge.Start.Index), new GH_Path(counter));
                IndicesTree.Append(new GH_Number(edge.End.Index), new GH_Path(counter));
                counter++;
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
    }

    internal class Node : GH_Point
    {
        [JsonIgnore]
        public bool Anchor { get; set; }

        [JsonIgnore]
        public List<Node> Neighbors { get; set; }

        [JsonIgnore]
        public int Index { get; internal set; }

        public Node()
        {
            Anchor = false;
            Neighbors = new List<Node>();
            Value = new Point3d();
            Index = -1;
        }

        public Node(Node other)
        {
            Anchor = other.Anchor;
            Neighbors = other.Neighbors;
            Value = other.Value;
            Index = other.Index;
        }

        public Node(bool anchor, List<Node> neighbors)
        {
            Anchor = anchor;
            Neighbors = neighbors;
            Index = -1;
        }
    }

    internal class Edge : GH_Curve
    {
        [JsonIgnore]
        public Node Start { get; set; }

        [JsonIgnore]
        public Node End { get; set; }

        [JsonIgnore]
        public double Q { get; set; }

        public Edge()
        {
            Start = new();
            End = new();
            Q = 0;
        }

        public Edge(Node start, Node end, double q, GH_Curve curve)
        {
            Start = start;
            End = end;
            Q = q;
            Value = curve.Value;
        }

        public Edge(Edge other)
        {
            Start = other.Start;
            End = other.End;
            Q = other.Q;
            Value = other.Value;
        }
    }
}
