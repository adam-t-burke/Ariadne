using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Text.Json;
using System.Text.Json.Serialization;
using Rhino.Collections;
using Rhino.Geometry;
using Ariadne.Graphs;
using Ariadne.Utilities;

namespace Ariadne.FDM
{
    internal class FDM_Network
    {
        /// <summary>
        /// Pass along the graph
        /// </summary>
        public Graph Graph { get; set; }

        /// <summary>
        /// List of anchor points
        /// </summary>
        [JsonIgnore(Condition = JsonIgnoreCondition.Always)]
        public List<Point3d> Anchors { get; set; }

        /// <summary>
        /// Tolerance for grabbing anchor points
        /// </summary>
        [JsonIgnore(Condition = JsonIgnoreCondition.Always)]
        public double ATol { get; set; }

        /// <summary>
        /// Tolerance end-to-end overlap tolerance for merging curves into a graph
        /// </summary>
        [JsonIgnore(Condition = JsonIgnoreCondition.Always)]
        public double ETol { get; set; }

        /// <summary>
        /// List of free nodes
        /// </summary>
        public List<Node> Free { get; set; }

        /// <summary>
        /// List of fixed nodes
        /// </summary>
        public List<Node> Fixed { get; set; }

        /// <summary>
        /// Boolean flag indicating if the network is valid
        /// </summary>
        [JsonIgnore(Condition = JsonIgnoreCondition.Always)]
        public bool Valid { get; set; }

        /// <summary>
        /// Empty constructor
        /// </summary>
        public FDM_Network()
        {
            Valid = false;
            Free = new List<Node>();
            Fixed = new List<Node>();
        }

        /// <summary>
        /// Constructor from a list of curves
        /// </summary>
        /// <param name="_Edges"></param>
        /// <param name="_ETol"></param>
        /// <param name="_Anchors"></param>
        /// <param name="_ATol"></param>
        public FDM_Network(List<Curve> _Edges, double _ETol, List<Point3d> _Anchors, double _ATol)
        {
            //Set fields
            Graph = new Graph(_Edges, _ETol);
            Anchors = _Anchors;
            ATol = _ATol;
            Free = new List<Node>();
            Fixed = new List<Node>();
            //Identify fixed and free nodes
            FixedFree();
            ValidCheck();
        }

        /// <summary>
        /// Constructor from another graph
        /// </summary>
        /// <param name="graph"></param>
        /// <param name="anchors"></param>
        /// <param name="_ForceDensities"></param>
        /// <param name="_Tolerance"></param>
        public FDM_Network(Graph _Graph, List<Point3d> _Anchors, double _ATol)
        {
            //Set fields
            Graph = _Graph;
            Anchors = _Anchors;
            ATol = _ATol;
            Free = new List<Node>();
            Fixed = new List<Node>();

            //Identify fixed and free nodes
            FixedFree();
            ValidCheck();
        }

        /// <summary>
        /// Construct a copy of an existing network
        /// </summary>
        /// <param name = "other" ></ param >
        public FDM_Network(FDM_Network other)
        {
            Graph = other.Graph;
            Anchors = other.Anchors;
            ATol = other.ATol;
            Free = other.Free;
            Fixed = other.Fixed;

            FixedFree();
            ValidCheck();
        }

        /// <summary>
        /// Identifies fixed and free nodes. Updates the order of Graph.Nodes to place fixed nodes at the end
        /// </summary>
        private void FixedFree()
        {
            List<Node> nodes = Graph.Nodes;

            Parallel.ForEach(Anchors, (anchor, state) =>
            {
                (bool anchorFound, int ianchor) = UtilityFunctions.WithinTolerance(nodes, anchor, ATol);
                if (anchorFound)
                {
                    nodes[ianchor].Anchor = true;
                    Fixed.Add(nodes[ianchor]);
                }
            });

            Free = nodes.Except(Fixed).ToList();

            // update the order of nodes to place fixed nodes at the end
            Graph.Nodes = Free.Concat(Fixed).ToList();
        }

        /// <summary>
        /// Checks if the network is valid
        /// </summary>
        private void ValidCheck()
        {
            // must have fixed nodes
            if (Anchors.Count < 2) Valid = false;

            // check N/F was properly captured
            else if (Fixed.Count != Anchors.Count || Fixed.Count + Free.Count != Graph.Nn) Valid = false;

            // all conditions pass
            else Valid = true;
        }

        /// <summary>
        /// checks for sufficient number of anchors
        /// </summary>
        /// <returns></returns>
        public bool AnchorCheck()
        {
            if (Anchors.Count < 2) return false;
            else return true;
        }

        /// <summary>
        /// Checks that N and F were properly generated
        /// </summary>
        /// <returns></returns>
        public bool NFCheck()
        {
            if (Fixed.Count != Anchors.Count || Fixed.Count + Free.Count != Graph.Nn) return false;
            else return true;
        }
      
    }
}
