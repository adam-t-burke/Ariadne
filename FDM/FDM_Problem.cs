using Ariadne.Optimization;
using Rhino.Collections;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.Json.Serialization;
using Ariadne.Graphs;
using Ariadne.Objectives;

namespace Ariadne.FDM
{
    internal class FDM_Problem
    {
        public List<double> Q { get; set; }
        public List<int> I { get; set; }
        public List<int> J { get; set; }
        public List<int> V { get; set; }
        public List<double[]> P { get; set; }
        public List<double[]> XYZf { get; set; }
        public FDM_Network Network { get; set; }

        [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)]
        public List<int> LoadNodes { get; set; }

        [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)]
        public OBJParameters Parameters { get; set; }

        [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)]
        public List<Anchor> VariableAnchors { get; set; }

        [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)]
        public List<int> VariableAnchorIndices { get; set; }

        [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingNull)]
        public List<int> FixedAnchorIndices { get; set; }

        public FDM_Problem() { }

        // Centralized constructor logic
        private void Initialize(
            FDM_Network network,
            List<double> q,
            List<Vector3d> p,
            OBJParameters parameters = null,
            List<int> loadNodes = null,
            List<Anchor> variableAnchors = null)
        {
            Q = q;
            Network = network;
            Parameters = parameters;
            LoadNodes = loadNodes;
            VariableAnchors = variableAnchors;

            ExtractPxyz(p);
            IJV(network);
            ExtractXYZf();

            if (Parameters != null)
                GetParameterIndices(Parameters);
        }

        // Refactored constructors
        public FDM_Problem(FDM_Network network, List<double> q, List<Vector3d> p)
            => Initialize(network, q, p);

        public FDM_Problem(FDM_Network network, List<double> q, List<Vector3d> p, List<int> loadNodes)
            => Initialize(network, q, p, null, loadNodes);

        public FDM_Problem(FDM_Network network, List<double> q, List<Vector3d> p, OBJParameters parameters)
            => Initialize(network, q, p, parameters);

        public FDM_Problem(FDM_Network network, List<double> q, List<Vector3d> p, OBJParameters parameters, List<int> loadNodes)
            => Initialize(network, q, p, parameters, loadNodes);

        public FDM_Problem(FDM_Network network, List<double> q, List<Vector3d> p, OBJParameters parameters, List<Anchor> variableAnchors)
            => Initialize(network, q, p, parameters, null, variableAnchors);

        public FDM_Problem(FDM_Network network, List<double> q, List<Vector3d> p, OBJParameters parameters, List<int> loadNodes, List<Anchor> variableAnchors)
            => Initialize(network, q, p, parameters, loadNodes, variableAnchors);

        private void GetParameterIndices(OBJParameters parameters)
        {
            var objs = parameters.Objectives;
            var edgeType = typeof(OBJEdges);
            var nodeType = typeof(OBJNodes);

            foreach (var obj in objs)
            {
                if (obj.GetType().IsSubclassOf(edgeType))
                    ((OBJEdges)obj).SetIndices(Network, ((OBJEdges)obj).Edges);
                if (obj.GetType().IsSubclassOf(nodeType))
                    ((OBJNodes)obj).SetIndices(Network, ((OBJNodes)obj).Nodes);
            }
        }

        private void ExtractPxyz(List<Vector3d> pxyz)
        {
            P = new List<double[]>(pxyz.Count);
            for (int i = 0; i < pxyz.Count; i++)
                P.Add(new[] { pxyz[i].X, pxyz[i].Y, pxyz[i].Z });
        }

        private void ExtractXYZf()
        {
            var nodes = Network.Fixed;
            XYZf = new List<double[]>(nodes.Count);
            for (int i = 0; i < nodes.Count; i++)
            {
                var p = nodes[i].Value;
                XYZf.Add(new[] { p.X, p.Y, p.Z });
            }
        }

        private void IJV(FDM_Network network)
        {
            int edgeCount = network.Graph.Ne;
            int nodeCount = network.Graph.Nodes.Count;

            I = new List<int>(edgeCount * 2);
            J = new List<int>(edgeCount * 2);
            V = new List<int>(edgeCount * 2);

            // Build node index map once for O(1) lookup
            var nodeIndexMap = new Dictionary<Node, int>(nodeCount);
            for (int idx = 0; idx < nodeCount; idx++)
                nodeIndexMap[network.Graph.Nodes[idx]] = idx;

            for (int i = 0; i < edgeCount; i++)
            {
                var edge = network.Graph.Edges[i];
                int start = nodeIndexMap[edge.Start];
                int end = nodeIndexMap[edge.End];

                I.Add(i); J.Add(start); V.Add(-1);
                I.Add(i); J.Add(end);   V.Add(1);
            }
        }
    }
}
