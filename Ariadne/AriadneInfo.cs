using Grasshopper;
using Grasshopper.Kernel;
using System;
using System.Drawing;

namespace Ariadne
{
    public class AriadneInfo : GH_AssemblyInfo
    {
        public override string Name => "Ariadne";

        //Return a 24x24 pixel bitmap to represent this GHA library.
        public override Bitmap Icon => null;

        //Return a short string describing the purpose of this GHA library.
        public override string Description => "Grasshopper plugin for the inverse design of form found structures in Rhino 8.";

        public override Guid Id => new Guid("86cc04a3-959e-4c13-8595-a4603ef43c7e");

        //Return a string identifying you or your company.
        public override string AuthorName => "Adam Burke";

        //Return a string representing your preferred contact details.
        public override string AuthorContact => "aburke3@mit.edu";
    }
}