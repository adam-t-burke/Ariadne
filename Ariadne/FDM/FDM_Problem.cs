using Ariadne.Optimization;
using Rhino.Collections;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Security.Authentication.ExtendedProtection;
using System.Text;
using System.Threading.Tasks;
using System.Text.Json;
using System.Text.Json.Serialization;
using Ariadne.Graphs;
using Eto.Forms;
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

        public FDM_Problem()
        {
        }

        public FDM_Problem(FDM_Network _network, List<double> _Q, List<Vector3d> _P)
        {
            Q = _Q; // force density vector
            Network = _network;
            ExtractPxyz(_P); //force vectors for each axis
            IJV(Network); // CSC format of connectivity matrix
            ExtractXYZf();

        }
        public FDM_Problem(FDM_Network _network, List<double> _Q, List<Vector3d> _P, List<int> _LoadNodes)
        {
            Q = _Q; // force density vector
            Network = _network;
            ExtractPxyz(_P); //force vectors for each axis
            IJV(Network); // CSC format of connectivity matrix
            LoadNodes = _LoadNodes;
            ExtractXYZf();

        }

        public FDM_Problem(FDM_Network _Network, List<double> _Q, List<Vector3d> _P, OBJParameters _Parameters)
        {
            Q = _Q; // force density vector
            Network = _Network;
            Parameters = _Parameters; // optimization parameters
            GetParameterIndices(Parameters);
            ExtractPxyz(_P); //force vectors for each axis
            IJV(Network); // CSC format of connectivity matrix
            ExtractXYZf();

        }
        public FDM_Problem(FDM_Network _Network, List<double> _Q, List<Vector3d> _P, OBJParameters _Parameters, List<int> _LoadNodes)
        {
            Q = _Q; // force density vector
            Network = _Network;
            Parameters = _Parameters; // optimization parameters
            GetParameterIndices(Parameters);
            ExtractPxyz(_P); //force vectors for each axis
            IJV(Network); // CSC format of connectivity matrix
            LoadNodes = _LoadNodes;
            ExtractXYZf();

        }

        public FDM_Problem(FDM_Network _Network, List<double> _Q, List<Vector3d> _P, OBJParameters _Parameters, List<Anchor> _Anchors)
        {
            Q = _Q; // force density vector
            Network = _Network;
            Parameters = _Parameters; // optimization parameters
            GetParameterIndices(Parameters);
            ExtractPxyz(_P); //force vectors for each axis
            IJV(Network); // CSC format of connectivity matrix
            ExtractXYZf();
            VariableAnchors = _Anchors;
        }

        public FDM_Problem(FDM_Network _Network, List<double> _Q, List<Vector3d> _P, OBJParameters _Parameters, List<int> _LoadNodes, List<Anchor> _Anchors)
        {
            Q = _Q; // force density vector
            Network = _Network;
            Parameters = _Parameters; // optimization parameters
            GetParameterIndices(Parameters);
            ExtractPxyz(_P); //force vectors for each axis
            IJV(Network); // CSC format of connectivity matrix
            LoadNodes = _LoadNodes;
            ExtractXYZf();
            VariableAnchors = _Anchors;
        }

        private void GetParameterIndices(OBJParameters _parameters)
        {
            List<OBJ> objs = _parameters.Objectives;
            Type edge_type = typeof(OBJEdges);
            Type node_type = typeof(OBJNodes);

            foreach (OBJ obj in objs)
            {
                if(obj.GetType().IsSubclassOf(edge_type))
                {
                    OBJEdges edges = (OBJEdges)obj;
                    edges.SetIndices(Network, edges.Edges);
                }
                if(obj.GetType().IsSubclassOf(node_type))
                {
                    OBJNodes nodes = (OBJNodes)obj;
                    nodes.SetIndices(Network, nodes.Nodes);
                }
            }
        }


        

        private void ExtractPxyz(List<Vector3d> PXYZ)
        {
            P = new List<double[]>();

            foreach (Vector3d p in PXYZ)
            {
                double[] _P = new double[] { p.X, p.Y, p.Z };
                P.Add(_P);
            }           
        }

        private void ExtractXYZf()
        {
            List<Node> nodes = Network.Fixed;

            XYZf = new List<double[]>();

            foreach (Node node in nodes)
            {
                Point3d p = node.Value;
                double[] _XYZ = new double[] { p.X, p.Y, p.Z };
                
                XYZf.Add(_XYZ);
            }            
        }

        private void IJV(FDM_Network network)
        {
            I = new List<int>();
            J = new List<int>();
            V = new List<int>();

            for(int i = 0; i < network.Graph.Ne; i++)
            {
                int start = network.Graph.Nodes.IndexOf(network.Graph.Edges[i].Start);
                int end = network.Graph.Nodes.IndexOf(network.Graph.Edges[i].End);

                I.Add(i);
                J.Add(start);
                V.Add(-1);


                I.Add(i);
                J.Add(end);
                V.Add(1);
            }
        }
      
    }
}
