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
                    bool closeEnough = node.Position.EpsilonEquals(point, tolerance);
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
                lengths.Add(edge.Curve.GetLength());
            }
            return lengths;
        }

        public static List<double> GetForces(List<double> lengths, List<double> forceDensities)
        {
            List<double> forces = new List<double>();
            for (int i = 0; i < lengths.Count; i++)
            {
                forces.Add(lengths[i] * forceDensities[i]);
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


    }
}
