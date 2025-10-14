using Ariadne.FDM;
using Rhino.Collections;
using Rhino.Display;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Ariadne.Graphs;
using System.Linq.Expressions;
using Ed.Eto;

namespace Ariadne.Utilities
{
    internal class UtilityFunctions
    {
        public static bool WithinTolerance(Point3dList points, Point3d point, double tolerance)
        {
            try
            {
                bool closeEnough = Point3dList.ClosestPointInList(points, point).EpsilonEquals(point, tolerance);
                return closeEnough;
            }
            catch
            {
                return false;
            }
        }

        public static (bool, int) WithinTolerance(List<Node> nodes, Point3d point, double tolerance)
        {
            bool match = false;
            int index = -1;

            Parallel.ForEach(nodes, (node, state) =>
                {
                    bool closeEnough = node.Value.EpsilonEquals(point, tolerance);
                    if (closeEnough)
                    {
                        match = true;
                        index = nodes.IndexOf(node);
                        state.Break();
                    }
                });

            return (match, index);


        }

        public static List<double> GetLengths(CurveList curves)
        {
            List<double> lengths = new List<double>();
            foreach (Curve curve in curves)
            {
                lengths.Add(curve.GetLength());
            }
            return lengths;
        }

        public static List<double> GetLengths(List<Edge> edges)
        {

            List<double> lengths = new List<double>();
            foreach (Edge edge in edges)
            {
                lengths.Add(edge.Value.GetLength());
            }
            return lengths;
        }

        public static List<double> GetForces(List<double> lengths, List<Edge> forceDensities)
        {
            List<double> forces = new();
            for (int i = 0; i < lengths.Count; i++)
            {
                forces.Add(lengths[i] * forceDensities[i].Q);
            }
            return forces;
        }

        public static List<double[]> PointsToArray(List<Point3d> points)
        {
            List<double[]> xyz = new List<double[]>();

            foreach (Point3d p in points)
            {
                double[] _XYZ = new double[] { p.X, p.Y, p.Z };

                xyz.Add(_XYZ);
            }
            return xyz;
        }

        public static List<double[]> VectorsToArray(List<Vector3d> vectors)
        {
            List<double[]> vecs = new List<double[]>();
            foreach (Vector3d v in vectors)
            {
                double[] _XYZ = new double[] { v.X, v.Y, v.Z };
                vecs.Add(_XYZ);
            }
            return vecs;
        }

        public static List<Vector3d> GetReactions(List<double> forces, FDM_Network network)
        {
            List<Vector3d> reactions = new();
            foreach(Node anchor in network.Fixed)
            {
                Vector3d reaction = new(0,0,0);
                foreach(Node neighbor in anchor.Neighbors)
                {
                    foreach(Edge edge in network.Graph.Edges)
                    {
                        if (edge.Start == anchor && edge.End == neighbor || edge.Start == neighbor && edge.End == anchor)
                        {
                            double force = forces[network.Graph.Edges.IndexOf(edge)];
                            Vector3d direction = anchor.Value - neighbor.Value;
                            direction.Unitize();
                            reaction += direction * force;
                        }
                    }
                }
                reactions.Add(reaction);
            }
            return reactions;
        }


    }
}
