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
        public string OBJID { get; set; }

        [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingDefault)]
        public double Weight { get; set; }

        [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingDefault)]
        public List<int> Indices { get; set; }

        [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingDefault)]
        public List<double[]> Points { get; set; }


        [JsonIgnore(Condition = JsonIgnoreCondition.WhenWritingDefault)]
        public List<double> Values { get; set; }

        [JsonIgnore(Condition = JsonIgnoreCondition.Always)]
        public bool IsValid { get; set; }

        public OBJ() { }
    }

    internal class OBJNull : OBJ
    {
        
        public OBJNull()
        {
            OBJID = "None";
            Weight = 0;
            IsValid = false;
        }
    }

    internal class OBJTarget : OBJ
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
        /// <param name="_points">List of points to apply to</param>

        public OBJTarget(double _weight, List<Point3d> _points)
        {
            OBJID = "Target";
            Weight = _weight;
            Indices = new List<int>() { -1 };
            Points = UtilityFunctions.PointsToArray(_points);
            IsValid = true;
        }

        /// <summary>
        /// Constructor for targeting a geometry that uses a subset of nodes in the network as targets.
        /// </summary>
        /// <param name="_weight"></param>
        /// <param name="_nodes"></param>
        /// <param name="_points"></param>
        public OBJTarget(double _weight, List<int> _nodes, List<Point3d> _points)
        {
            OBJID = "Target";
            Weight = _weight;
            Indices = _nodes;
            Points = UtilityFunctions.PointsToArray(_points);
            if (_nodes.Count == _points.Count)
            {
                IsValid = true;
            }
            else
            {
                IsValid = false;
            }
        }
    }

    internal class OBJlengthvariation : OBJ
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
            OBJID = "LengthVar";
            Weight = _weight;
            Indices = new List<int> { -1 };
            IsValid = true;
        }

        /// <summary>
        /// Constructor with weight and list of edges to apply to.
        /// </summary>
        /// <param name="_weight"></param>
        /// <param name="_edges"></param>
        public OBJlengthvariation(double _weight, List<int> _edges)
        {
            OBJID = "LengthVar";
            Weight = _weight;
            Indices = _edges;
            IsValid = true;
        }

    }

    internal class OBJforcevariation : OBJ
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
            OBJID = "ForceVar";
            Weight = _weight;
            Indices = new List<int>() { -1 };
            IsValid = true;
        }

        /// <summary>
        /// Constructor with weight and list of edges to apply to.
        /// </summary>
        /// <param name="_weight"></param>
        /// <param name="_edges"></param>
        public OBJforcevariation(double _weight, List<int> _edges)
        {
            OBJID = "ForceVar";
            Weight = _weight;
            Indices = _edges;
            IsValid = true;
        }

    }

    internal class OBJPerformance : OBJ
    {
        public OBJPerformance()
        {
            IsValid = false;
        }
        public OBJPerformance(double _weight)
        {
            OBJID = "Performance";
            Weight = _weight;
            IsValid = true;
        }
    }

    internal class OBJMinlength : OBJ
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
        public OBJMinlength(double _weight, double _value)
        {
            OBJID = "MinLength";
            Weight = _weight;
            Values = new List<double>() { _value };
            Indices = new List<int>() { -1 };
            IsValid = true;
        }

        /// <summary>
        /// Minimum length objective function from a single value and list of edges to apply to.
        /// </summary>
        /// <param name="_weight"></param>
        /// <param name="_value"></param>
        /// <param name="_edges"></param>
        public OBJMinlength(double _weight, double _value, List<int> _edges)
        {
            OBJID = "MinLength";
            Weight = _weight;
            Values = new List<double>() { _value };
            Indices = _edges;
            IsValid = true;
        }

        /// <summary>
        /// Minimum length objective function from a list of values and list of edges to apply to.
        /// </summary>
        /// 
        public OBJMinlength(double _weight, List<double> _values, List<int> _edges)
        {
            OBJID = "MinLength";
            Weight = _weight;
            Values = _values;
            Indices = _edges;
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

    internal class OBJMaxlength : OBJ
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
        public OBJMaxlength(double _weight, double _value)
        {
            OBJID = "MaxLength";
            Weight = _weight;
            Values = new List<double>() { _value };
            Indices = new List<int>() { -1 };
            IsValid = true;
        }

        /// <summary>
        /// Maximum length objective function from a single value and list of edges to apply to.
        /// </summary>
        /// <param name="_weight"></param>
        /// <param name="_value"></param>
        /// <param name="_edges"></param>
        public OBJMaxlength(double _weight, double _value, List<int> _edges)
        {
            OBJID = "MaxLength";
            Weight = _weight;
            Values = new List<double>() { _value };
            Indices = _edges;
            IsValid = true;
        }

        /// <summary>
        /// Maximum length objective function from a list of values and list of edges to apply to.
        /// </summary>
        /// 
        public OBJMaxlength(double _weight, List<double> _value, List<int> _edges)
        {
            OBJID = "MaxLength";
            Weight = _weight;
            Values = _value;
            Indices = _edges;
            if (_value.Count == _edges.Count)
            {
                IsValid = true;
            }
            else
            {
                IsValid = false;
            }
        }
    }

    internal class OBJMinforce : OBJ
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
            OBJID = "MinForce";
            Weight = _weight;
            Values = _values;
            Indices.Add(-1);
            IsValid = true;
        }

        /// <summary>
        /// Minimum force objective function from a single value and list of edges to apply to.
        /// </summary>
        /// <param name="_weight"></param>
        /// <param name="_value"></param>
        /// <param name="_edges"></param>
        public OBJMinforce(double _weight, double _value, List<int> _edges)
        {
            OBJID = "MinForce";
            Weight = _weight;
            Values = new List<double>() { _value };
            Indices = _edges;
            IsValid = true;
        }

        /// <summary>
        /// Minimum force objective function from a list of values and list of edges to apply to.
        /// </summary>
        /// <param name="_weight"></param>
        /// <param name="_values"></param>
        /// <param name="_edges"></param>
        public OBJMinforce(double _weight, List<double> _values, List<int> _edges)
        {
            OBJID = "MinForce";
            Weight = _weight;
            Values = _values;
            Indices = _edges;
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

    internal class OBJMaxforce : OBJ
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
            OBJID = "MaxForce";
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
        public OBJMaxforce(double _weight, double _value, List<int> _edges)
        {
            OBJID = "MaxForce";
            Weight = _weight;
            Values = new List<double>() { _value };
            Indices = new List<int>() { -1 };
            IsValid = true;
        }

        /// <summary>
        /// Minimum force objective function from a list of values and list of edges to apply to.
        /// </summary>
        /// <param name="_weight"></param>
        /// <param name="_values"></param>
        /// <param name="_edges"></param>
        public OBJMaxforce(double _weight, List<double> _values, List<int> _edges)
        {
            OBJID = "MaxForce";
            Weight = _weight;
            Values = _values;
            Indices = _edges;
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

    internal class OBJTargetLength : OBJ
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
            OBJID = "TargetLen";
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
        public OBJTargetLength(double _weight, double _value, List<int> _edges)
        {
            OBJID = "TargetLen";
            Weight = _weight;
            Values = new List<double>() { _value };
            Indices = _edges;
            IsValid = true;
        }

        /// <summary>
        /// Minimum force objective function from a list of values and list of edges to apply to.
        /// </summary>
        /// <param name="_weight"></param>
        /// <param name="_values"></param>
        /// <param name="_edges"></param>
        public OBJTargetLength(double _weight, List<double> _values, List<int> _edges)
        {
            OBJID = "TargetLen";
            Weight = _weight;
            Values = _values;
            Indices = _edges;
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
}
