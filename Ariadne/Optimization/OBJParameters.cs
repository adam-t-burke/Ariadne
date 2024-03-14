using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Ariadne.Objectives;
using System.Text.Json;
using System.Text.Json.Serialization;
using Eto.Forms;

namespace Ariadne.Optimization
{

    /// <summary>
    /// Contains all parameters necessary for an optimization run
    /// </summary>
    internal class OBJParameters
    {
        public double LowerBound { get; set; }
        public double UpperBound { get; set; }
        public double AbsTol { get; set; }
        public double RelTol { get; set; }
        public List<OBJ> Objectives { get; set; }
        public bool ShowIterations { get; set; }
        public int UpdateFrequency { get; set; }
        public int MaxIterations { get; set; }

        public bool NodeTrace { get; set; }

        public OBJParameters()
        {

        }

        public OBJParameters(double lb, double ub, double abstol, double reltol, List<OBJ> objs, bool showiter, int updatefreq, int maxiters, bool nodeTrace)
        {
            LowerBound = lb;
            UpperBound = ub;
            AbsTol = abstol;
            RelTol = reltol;
            Objectives = objs;
            ShowIterations = showiter;
            UpdateFrequency = updatefreq;
            MaxIterations = maxiters;
            NodeTrace = nodeTrace;
        }

        public OBJParameters(double lb, double ub, double abstol, double reltol, OBJ objs, bool showiter, int updatefreq, int maxiters, bool nodeTrace)
        {
            LowerBound = lb;
            UpperBound = ub;
            AbsTol = abstol;
            RelTol = reltol;
            ShowIterations = showiter;
            UpdateFrequency = updatefreq;

            List<OBJ> objectiveList = new List<OBJ> { objs };
            Objectives = objectiveList;

            MaxIterations = maxiters;
            NodeTrace = nodeTrace;
        }
    }
}
