using Rhino.Collections;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Text.Json;
using System.Text.Json.Serialization;
using System.Windows.Markup;
using Ariadne.FDM;
using Ariadne.Graphs;
using Rhino.Geometry;
using Ariadne.Utilities;

namespace Ariadne.Objectives
{
    /// <summary>
    /// Abstract class for objective functions
    /// Must include ID field and Weight field at a minimum
    /// </summary>
    internal abstract class OBJ
    {
        public int OBJID { get; set; }

        public double Weight { get; set; }

        public List<int> Indices { get; set; }

        [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingDefault)]
        public List<double[]> Points { get; set; }

        [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingDefault)]
        public List<double[]> Vectors { get; set; }

        [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingDefault)]
        public List<double> Values { get; set; }

        [JsonIgnore(Condition = JsonIgnoreCondition.Always)]
        public bool IsValid { get; set; }

        public OBJ() { }
    }

    internal class OBJEdges : OBJ
    {
        [JsonIgnore(Condition = JsonIgnoreCondition.Always)]
        public List<Edge> Edges { get; set; }
        public void SetIndices(FDM_Network network, List<Edge> edges)
        {
            if (Indices is null)
            {
                Indices = new();
                foreach (Edge edge in edges)
                {
                    Indices.Add(network.Graph.Edges.IndexOf(edge));
                }
            }
            else {
                
            }
        }
    }

    internal class OBJNodes : OBJ
    {
        [JsonIgnore(Condition = JsonIgnoreCondition.Always)]
        public List<Node> Nodes { get; set; }
        public void SetIndices(FDM_Network network, List<Node> nodes)
        {
            if (Indices is null)
            {
                Indices = new();
                foreach (Node node in nodes)
                {
                    Indices.Add(network.Graph.Nodes.IndexOf(node));
                }
                
            }
            else
            {
                if (Points is null)
                {
                    Points = UtilityFunctions.PointsToArray(network.Free.Select(x => x.Value).ToList());
                }
            }
         

        }
    }

    internal class OBJNull : OBJ
    {        
        public OBJNull()
        {
            OBJID = -1;
            Weight = 0;
            IsValid = false;
        }
    }

    internal class OBJTarget : OBJNodes
    {
        /// <summary>
        /// Empty construtor
        /// </summary>
        public OBJTarget()
        {
            IsValid = false;
        }

        /// <summary>
        /// Constructor for targeting a geometry that uses all free nodes in the network as targets.
        /// </summary> 
        /// <param name="_weight">Weight of objective function</param>

        public OBJTarget(double _weight)
        {
            OBJID = 1;
            Weight = _weight;
            Indices = new List<int>() { -1 };
            Points = null;
            IsValid = true;
        }

        /// <summary>
        /// Constructor for targeting a geometry that uses a subset of nodes in the network as targets.
        /// </summary>
        /// <param name="_weight"></param>
        /// <param name="_nodes"></param>
        /// <param name="_points"></param>
        public OBJTarget(double _weight, List<Node> _nodes)
        {
            OBJID = 1;
            Weight = _weight;
            Nodes = _nodes;
            Indices = null;
            Points = UtilityFunctions.PointsToArray(_nodes.Select(x => x.Value).ToList());
            IsValid = true;
        }
    }

    internal class OBJTarget_XY : OBJNodes
    {
        /// <summary>
        /// Empty construtor
        /// </summary>
        public OBJTarget_XY()
        {
            IsValid = false;
        }

        /// <summary>
        /// Constructor for targeting a geometry that uses all free nodes in the network as targets.
        /// </summary> 
        /// <param name="_weight">Weight of objective function</param>

        public OBJTarget_XY(double _weight)
        {
            OBJID = 10;
            Weight = _weight;
            Indices = new List<int>() { -1 };
            Points = null;
            IsValid = true;
        }

        /// <summary>
        /// Constructor for targeting a geometry that uses a subset of nodes in the network as targets.
        /// </summary>
        /// <param name="_weight"></param>
        /// <param name="_nodes"></param>
        /// <param name="_points"></param>
        public OBJTarget_XY(double _weight, List<Node> _nodes)
        {
            OBJID = 10;
            Weight = _weight;
            Nodes = _nodes;
            Indices = null;
            Points = UtilityFunctions.PointsToArray(_nodes.Select(x => x.Value).ToList());
            IsValid = true;
        }
    }

    internal class OBJlengthvariation : OBJEdges
    {
        /// <summary>
        /// Empty construtor
        /// </summary>
        public OBJlengthvariation()
        {
            IsValid = false;
        }

        /// <summary>
        /// Constructor with only weight to apply to all edges.
        /// </summary>
        /// <param name="_weight"></param>
        public OBJlengthvariation(double _weight)
        {
            OBJID = 2;
            Weight = _weight;
            Indices = new List<int> { -1 };
            IsValid = true;
        }

        /// <summary>
        /// Constructor with weight and list of edges to apply to.
        /// </summary>
        /// <param name="_weight"></param>
        /// <param name="_edges"></param>
        public OBJlengthvariation(double _weight, List<Edge> _edges)
        {
            OBJID = 2;
            Weight = _weight;
            Edges = _edges;
            Indices = null;
            IsValid = true;
        }

    }

    internal class OBJforcevariation : OBJEdges
    {
        /// <summary>
        /// Empty construtor
        /// </summary>
        public OBJforcevariation()
        {
            IsValid = false;
        }

        /// <summary>
        /// Constructor with only weight to apply to all edges.
        /// </summary>
        /// <param name="_weight"></param>
        public OBJforcevariation(double _weight)
        {
            OBJID = 3;
            Weight = _weight;
            Indices = new List<int>() { -1 };
            IsValid = true;
        }

        /// <summary>
        /// Constructor with weight and list of edges to apply to.
        /// </summary>
        /// <param name="_weight"></param>
        /// <param name="_edges"></param>
        public OBJforcevariation(double _weight, List<Edge> _edges)
        {
            OBJID = 3;
            Weight = _weight;
            Edges = _edges;
            Indices = null;
            IsValid = true;
        }

    }

    internal class OBJPerformance : OBJEdges
    {
        public OBJPerformance()
        {
            IsValid = false;
        }
        public OBJPerformance(double _weight)
        {
            OBJID = 4;
            Weight = _weight;
            Indices = new List<int>() { -1 };
            IsValid = true;
        }
    }

    internal class OBJMinlength : OBJEdges
    {
        /// <summary>
        /// Empty Constructor
        /// </summary>
        public OBJMinlength()
        {
            IsValid = false;
        }

        /// <summary>
        /// Minimum length objective function from a single value to apply to all edges.
        /// </summary>
        /// <param name="_weight"></param>
        /// <param name="_value"></param>
        public OBJMinlength(double _weight, List<double> _value)
        {
            OBJID = 5;
            Weight = _weight;
            Values = _value;
            Indices = new List<int>() { -1 };
            IsValid = true;
        }

        /// <summary>
        /// Minimum length objective function from a single value and list of edges to apply to.
        /// </summary>
        /// <param name="_weight"></param>
        /// <param name="_value"></param>
        /// <param name="_edges"></param>
        public OBJMinlength(double _weight, List<double> _value, List<Edge> _edges)
        {
            OBJID = 5;
            Weight = _weight;

            if (_edges.Count == 1 && _value.Count == 1)
            {
                Values = _value;
                Edges = _edges;
                Indices = null;
                IsValid = true;
            }
            else if (_value.Count == 1 && _edges.Count > 1)
            {
                Values = Enumerable.Repeat(_value[0], _edges.Count).ToList();
                Edges = _edges;
                Indices = null;
                IsValid = true;
            }
            else if (_edges.Count > 1 && _value.Count > 1)
            {
                if (_edges.Count == _value.Count)
                {
                    Values = _value;
                    Edges = _edges;
                    Indices = null;
                    IsValid = true;
                }
                else
                {
                    throw new Exception($"Number of values must match number of edges. Current count of values is {_value.Count}");
                }
            }     
            
            else
            {
                throw new Exception("Invalid objective function.");
            }
        }
    }

    internal class OBJMaxlength : OBJEdges
    {
        /// <summary>
        /// Empty constructor
        /// </summary>
        public OBJMaxlength()
        {
            IsValid = false;
        }

        /// <summary>
        /// Maximum length objective function from a single value to apply to all edges.
        /// </summary>
        /// <param name="_weight"></param>
        /// <param name="_value"></param>
        public OBJMaxlength(double _weight, List<double> _value)
        {
            OBJID = 6;
            Weight = _weight;
            Values = _value;
            Indices = new List<int>() { -1 };
            IsValid = true;
        }


        /// <summary>
        /// Maximum length objective function from a list of values and list of edges to apply to.
        /// </summary>
        /// 
        public OBJMaxlength(double _weight, List<double> _value, List<Edge> _edges)
        {
            OBJID = 6;
            Weight = _weight;

            if (_edges.Count == 1 && _value.Count == 1)
            {
                Values = _value;
                Edges = _edges;
                Indices = null;
                IsValid = true;
            }
            else if (_value.Count == 1 && _edges.Count > 1)
            {
                Values = Enumerable.Repeat(_value[0], _edges.Count).ToList();
                Edges = _edges;
                Indices = null;
                IsValid = true;
            }
            else if (_edges.Count > 1 && _value.Count > 1)
            {
                if (_edges.Count == _value.Count)
                {
                    Values = _value;
                    Edges = _edges;
                    Indices = null;
                    IsValid = true;
                }
                else
                {
                    throw new Exception($"Number of values must match number of edges. Current count of values is {_value.Count}");
                }
            }

            else
            {
                throw new Exception("Invalid objective function.");
            }
        }
    }

    internal class OBJMinforce : OBJEdges
    {
        /// <summary>
        /// Empty constructor
        /// </summary>
        public OBJMinforce()
        {
            IsValid = false;
        }

        /// <summary>
        /// Minimum force objective function from a single value to apply to all edges.
        /// </summary>
        /// <param name="_weight"></param>
        /// <param name="_values"></param>
        public OBJMinforce(double _weight, List<double> _values)
        {
            OBJID = 7;
            Weight = _weight;
            Values = _values;
            Indices = new List<int>() { -1 };
            IsValid = true;
        }

        /// <summary>
        /// Minimum force objective function from a list of values and list of edges to apply to.
        /// </summary>
        /// <param name="_weight"></param>
        /// <param name="_values"></param>
        /// <param name="_edges"></param>
        public OBJMinforce(double _weight, List<double> _value, List<Edge> _edges)
        {
            OBJID = 7;
            Weight = _weight;

            if (_edges.Count == 1 && _value.Count == 1)
            {
                Values = _value;
                Edges = _edges;
                Indices = null;
                IsValid = true;
            }
            else if (_value.Count == 1 && _edges.Count > 1)
            {
                Values = Enumerable.Repeat(_value[0], _edges.Count).ToList();
                Edges = _edges;
                Indices = null;
                IsValid = true;
            }
            else if (_edges.Count > 1 && _value.Count > 1)
            {
                if (_edges.Count == _value.Count)
                {
                    Values = _value;
                    Edges = _edges;
                    Indices = null;
                    IsValid = true;
                }
                else
                {
                    throw new Exception($"Number of values must match number of edges. Current count of values is {_value.Count}");
                }
            }

            else
            {
                throw new Exception("Invalid objective function.");
            }
        }
    }

    internal class OBJMaxforce : OBJEdges
    {
        /// <summary>
        /// Empty constructor
        /// </summary>
        public OBJMaxforce()
        {
            IsValid = false;
        }

        /// <summary>
        /// Minimum force objective function from a single value to apply to all edges.
        /// </summary>
        /// <param name="_weight"></param>
        /// <param name="_values"></param>
        public OBJMaxforce(double _weight, List<double> _values)
        {
            OBJID = 8;
            Weight = _weight;
            Values = _values;
            Indices = new List<int>() { -1 };
            IsValid = true;
        }

        /// <summary>
        /// Minimum force objective function from a list of values and list of edges to apply to.
        /// </summary>
        /// <param name="_weight"></param>
        /// <param name="_values"></param>
        /// <param name="_edges"></param>
        public OBJMaxforce(double _weight, List<double> _value, List<Edge> _edges)
        {
            OBJID = 8;
            Weight = _weight;

            if (_edges.Count == 1 && _value.Count == 1)
            {
                Values = _value;
                Edges = _edges;
                Indices = null;
                IsValid = true;
            }
            else if (_value.Count == 1 && _edges.Count > 1)
            {
                Values = Enumerable.Repeat(_value[0], _edges.Count).ToList();
                Edges = _edges;
                Indices = null;
                IsValid = true;
            }
            else if (_edges.Count > 1 && _value.Count > 1)
            {
                if (_edges.Count == _value.Count)
                {
                    Values = _value;
                    Edges = _edges;
                    Indices = null;
                    IsValid = true;
                }
                else
                {
                    throw new Exception($"Number of values must match number of edges. Current count of values is {_value.Count}");
                }
            }

            else
            {
                throw new Exception("Invalid objective function.");
            }
        }
    }

    internal class OBJPointSet : OBJNodes
    {
        /// <summary>
        /// Empty construtor
        /// </summary>
        public OBJPointSet()
        {
            IsValid = false;
        }

        /// <summary>
        /// Constructor for targeting a geometry that uses all free nodes in the network as targets.
        /// </summary> 
        /// <param name="_weight">Weight of objective function</param>

        public OBJPointSet(double _weight)
        {
            OBJID = 11;
            Weight = _weight;
            Indices = new List<int>() { -1 };
            Points = null;
            IsValid = true;
        }

        /// <summary>
        /// Constructor for targeting a geometry that uses a subset of nodes in the network as targets.
        /// </summary>
        /// <param name="_weight"></param>
        /// <param name="_nodes"></param>
        /// <param name="_points"></param>
        public OBJPointSet(double _weight, List<Node> _nodes)
        {
            OBJID = 11;
            Weight = _weight;
            Nodes = _nodes;
            Indices = null;
            Points = UtilityFunctions.PointsToArray(_nodes.Select(x => x.Value).ToList());
            IsValid = true;
        }

    }

    internal class OBJTargetLength : OBJEdges
    {
        /// <summary>
        /// Empty constructor
        /// </summary>
        public OBJTargetLength()
        {
            IsValid = false;
        }

        /// <summary>
        /// Minimum force objective function from a single value to apply to all edges.
        /// </summary>
        /// <param name="_weight"></param>
        /// <param name="_values"></param>
        public OBJTargetLength(double _weight, List<double> _values)
        {
            OBJID = 9;
            Weight = _weight;
            Values = _values;
            Indices = new List<int>() { -1 };
            IsValid = true;
        }

        /// <summary>
        /// Minimum force objective function from a single value and list of edges to apply to.
        /// </summary>
        /// <param name="_weight"></param>
        /// <param name="_value"></param>
        /// <param name="_edges"></param>
        public OBJTargetLength(double _weight, double _value, List<Edge> _edges)
        {
            OBJID = 9;
            Weight = _weight;
            Values = new List<double>() { _value };
            Edges = _edges;
            Indices = null;
            IsValid = true;
        }

        /// <summary>
        /// Minimum force objective function from a list of values and list of edges to apply to.
        /// </summary>
        /// <param name="_weight"></param>
        /// <param name="_values"></param>
        /// <param name="_edges"></param>
        public OBJTargetLength(double _weight, List<double> _values, List<Edge> _edges)
        {
            OBJID = 9;
            Weight = _weight;
            Values = _values;
            Edges = _edges;
            Indices = null;
            if (_values.Count == _edges.Count)
            {
                IsValid = true;
            }
            else
            {
                IsValid = false;
            }
        }
    }

    internal class OBJReaction : OBJNodes
    {
        public OBJReaction(double _weight, List<Node> _nodes, List<Vector3d> _vectors, bool _useMag)
        {
            Weight = _weight;
            Nodes = _nodes;
            Vectors = UtilityFunctions.VectorsToArray(_vectors);
            Indices = [.. Nodes.Select(x => x.Index)];

            if (!_useMag)
            {
                OBJID = 12;
            }
            else
            {
                OBJID = 13;
            }
            foreach (Node node in _nodes)
            {
                if (node.Anchor)
                {
                    IsValid = true;
                }
                else
                {
                    IsValid = false;
                    throw new Exception("All nodes must be fixed for a reaction objective.");
                }
            }
        }
    }
}
