using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Ariadne.Optimization
{
    /// <summary>
    /// Class to directly parse output JSON messages from server
    /// </summary>
    internal class Receiver
    {
        public bool Finished { get; set; }
        public int Iter { get; set; }
        public double Loss { get; set; }
        public List<double> Q { get; set; }
        public List<double> X { get; set; }
        public List<double> Y { get; set; }
        public List<double> Z { get; set; }
        public List<double> Losstrace { get; set; }
    }
}
