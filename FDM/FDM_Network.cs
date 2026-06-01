using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;
using System.Text.Json;
using System.Text.Json.Serialization;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino.Collections;
using Rhino.Geometry;
using Ariadne.Graphs;
using Ariadne.Utilities;

namespace Ariadne.FDM
{
    /// <summary>
    /// Force density method network: topology (<see cref="Graph"/>), fixed/free node partition,
    /// and anchor positions. Consumed by the solver for forward and optimization solves.
    /// </summary>
    public class FDM_Network
    {
        /// <summary>
        /// Underlying graph (nodes and edges) for this network.
        /// </summary>
        public Graph Graph { get; set; }

        /// <summary>
        /// List of anchor points
        /// </summary>
        [JsonIgnore]
        public List<Point3d> Anchors { get; set; }

        /// <summary>
        /// Tolerance for grabbing anchor points
        /// </summary>
        [JsonIgnore]
        public double ATol { get; set; }

        /// <summary>
        /// Tolerance end-to-end overlap tolerance for merging curves into a graph
        /// </summary>
        [JsonIgnore]
        public double ETol { get; set; }

        /// <summary>
        /// List of free nodes
        /// </summary>
        [JsonIgnore]
        public List<Node> Free { get; set; }

        /// <summary>Indices of free nodes in <see cref="Graph.Nodes"/> (nodes that are not anchors).</summary>
        public List<int> FreeNodes { get; set; }

        /// <summary>
        /// List of fixed nodes
        /// </summary>
        [JsonIgnore]
        public List<Node> Fixed { get; set; }

        /// <summary>Indices of fixed (anchor) nodes in <see cref="Graph.Nodes"/>.</summary>
        public List<int> FixedNodes { get; set; }

        /// <summary>
        /// Boolean flag indicating if the network is valid
        /// </summary>
        [JsonIgnore]
        public bool Valid { get; set; }

        /// <summary>When true, the network is being updated (e.g. during solver callbacks); used to avoid re-entrancy.</summary>
        public bool IsUpdating { get; set; }

        /// <summary>
        /// Empty constructor
        /// </summary>
        public FDM_Network()
        {
            Valid = false;
            Free = new List<Node>();
            Fixed = new List<Node>();
            FixedNodes = new List<int>();
            FreeNodes = new List<int>();
            Graph = new Graph();
            Anchors = new List<Point3d>();
            ATol = 0.001;
            ETol = 0.001;
            IsUpdating = false;
        }

        /// <summary>
        /// Builds an FDM network from curves, merging endpoints within tolerance and classifying nodes by anchor positions.
        /// </summary>
        /// <param name="_Edges">Input curves; endpoints become nodes, curves become edges.</param>
        /// <param name="_ETol">End-to-end tolerance for building the graph (vertex merging).</param>
        /// <param name="_Anchors">Fixed support positions; nodes within <paramref name="_ATol"/> of these become fixed.</param>
        /// <param name="_ATol">Distance tolerance for matching nodes to anchor positions.</param>
        public FDM_Network(List<GH_Curve> _Edges, double _ETol, List<Point3d> _Anchors, double _ATol)
        {
            //Set fields
            Graph = new Graph(_Edges, _ETol);
            Anchors = _Anchors;
            ATol = _ATol;
            Free = new List<Node>();
            Fixed = new List<Node>();
            FixedNodes = new List<int>();
            FreeNodes = new List<int>();
            //Identify fixed and free nodes
            FixedFree();
            ValidCheck();

            IsUpdating = false;
        }

        /// <summary>
        /// Builds an FDM network from an existing graph and anchor positions.
        /// </summary>
        /// <param name="_Graph">Pre-built graph (nodes and edges).</param>
        /// <param name="_Anchors">Fixed support positions; nodes within <paramref name="_ATol"/> become fixed.</param>
        /// <param name="_ATol">Distance tolerance for matching nodes to anchor positions.</param>
        public FDM_Network(Graph _Graph, List<Point3d> _Anchors, double _ATol)
        {
            //Set fields
            Graph = _Graph;
            Anchors = _Anchors;
            ATol = _ATol;
            Free = new List<Node>();
            Fixed = new List<Node>();
            FixedNodes = new List<int>();
            FreeNodes = new List<int>();

            //Identify fixed and free nodes
            FixedFree();
            ValidCheck();

            IsUpdating = false;
        }

        /// <summary>
        /// Construct a copy of an existing network.
        /// Copies lists to avoid shared-mutation issues.
        /// Does NOT re-run FixedFree since the source is already partitioned.
        /// </summary>
        /// <param name="other">Network to copy from.</param>
        public FDM_Network(FDM_Network other)
        {
            Graph = new Graph(other.Graph);
            Anchors = new List<Point3d>(other.Anchors);
            ATol = other.ATol;
            ETol = other.ETol;
            Free = other.FreeNodes.ConvertAll(i => Graph.Nodes[i]);
            Fixed = other.FixedNodes.ConvertAll(i => Graph.Nodes[i]);
            FixedNodes = new List<int>(other.FixedNodes);
            FreeNodes = new List<int>(other.FreeNodes);
            Valid = other.Valid;
            IsUpdating = false;
        }

        /// <summary>
        /// Identifies fixed and free nodes. Updates the order of Graph.Nodes to place fixed nodes at the end
        /// </summary>
        private void FixedFree()
        {
            List<Node> nodes = Graph.Nodes;
            var fixedNodes = new HashSet<Node>();

            foreach (Point3d anchor in Anchors)
            {
                int ianchor = AnchorMatcher.FindNodeIndex(nodes, anchor, ATol);
                if (ianchor >= 0)
                {
                    nodes[ianchor].Anchor = true;
                    fixedNodes.Add(nodes[ianchor]);
                }
            }

            Fixed = fixedNodes.ToList();
            Free = nodes.Except(Fixed).ToList();

            // update the order of nodes to place fixed nodes at the end
            Graph.Nodes = Free.Concat(Fixed).ToList();

            // update the indices of fixed and free nodes
            Graph.UpdateNodeIndices();


            foreach (Node node in Graph.Nodes)
                {
                    if (node.Anchor == false)
                    {
                        FreeNodes.Add(node.Index);
                    }
                    else
                    {
                        FixedNodes.Add(node.Index);
                    }
                };
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
        /// Checks that the network has at least two anchor positions (required for a valid FDM setup).
        /// </summary>
        /// <returns>True if Anchors has at least two elements.</returns>
        public bool AnchorCheck()
        {
            if (Anchors.Count < 2) return false;
            else return true;
        }

        /// <summary>
        /// Checks that fixed and free node counts match anchors and total node count.
        /// </summary>
        /// <returns>True if Fixed.Count equals Anchors.Count and Fixed + Free equals Graph.Nn.</returns>
        public bool NFCheck()
        {
            if (Fixed.Count != Anchors.Count || Fixed.Count + Free.Count != Graph.Nn) return false;
            else return true;
        }

        /// <summary>
        /// Content-based hash of network topology and current geometry for solve input-change detection.
        /// Excludes solved edge q values so optimizer output does not become implicit input state.
        /// </summary>
        public int GetTopologyHashCode()
        {
            var h = new HashCode();
            h.Add(ATol);
            h.Add(ETol);
            h.Add(Graph.Nn);
            h.Add(Graph.Ne);
            h.Add(Graph.Tolerance);

            foreach (var node in Graph.Nodes)
            {
                h.Add(node.Index);
                h.Add(node.Anchor);
                h.Add(node.Value.X);
                h.Add(node.Value.Y);
                h.Add(node.Value.Z);
            }

            foreach (var edge in Graph.Edges)
            {
                h.Add(edge.Start.Index);
                h.Add(edge.End.Index);
                h.Add(edge.ReferenceID);
                h.Add(edge.UserKey);
            }

            foreach (var anchor in Anchors)
            {
                h.Add(anchor.X);
                h.Add(anchor.Y);
                h.Add(anchor.Z);
            }

            return h.ToHashCode();
        }

        /// <summary>
        /// Content-based hash including mutable solved edge q state for diagnostics.
        /// Solver recompute logic should prefer <see cref="GetTopologyHashCode"/>.
        /// </summary>
        public int GetContentHashCode()
        {
            var h = new HashCode();
            h.Add(GetTopologyHashCode());
            foreach (var edge in Graph.Edges)
                h.Add(edge.Q);
            return h.ToHashCode();
        }
      
    }
}
