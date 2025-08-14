using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Ariadne.Optimization
{
    internal class Anchor
    {
        public int NodeIndex { get; set; }

        public double InitialX { get; set; }
        public double InitialY { get; set; }
        public double InitialZ { get; set; }



        public Anchor()
        {
        }

        public Anchor(int _nodeIndex, double _InitialX, double _InitialY, double _InitialZ)
        {
            NodeIndex = _nodeIndex;
            InitialX = _InitialX;
            InitialY = _InitialY;
            InitialZ = _InitialZ;
        }
    }
}
