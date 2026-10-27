using Grasshopper;
using Grasshopper.Kernel;
using System;
using System.Drawing;

namespace Ariadne
{
    /// <summary>
    /// Assembly metadata for the Ariadne Grasshopper plugin (name, description, author).
    /// </summary>
    public class AriadneInfo : GH_AssemblyInfo
    {
        /// <summary>
        /// Preview / breaking-change channel published on Yak as <c>ariadne-dev</c>.
        /// Distinct from the stable <c>Ariadne</c> assembly GUID so Package Restore
        /// and side-by-side installs do not collide with the public release.
        /// </summary>
        public override string Name => "Ariadne Dev";

        public override Bitmap? Icon => null;

        //Return a short string describing the purpose of this GHA library.
        public override string Description =>
            "PREVIEW / BREAKING: development build of Ariadne (Yak: ariadne-dev). "
            + "Not API-compatible with the stable Ariadne package — do not install "
            + "both at once. Grasshopper plugin for inverse design of form-found "
            + "structures in Rhino 8.";

        public override Guid Id => new Guid("d4f8a2e1-9c3b-4e7f-a1b2-6c8d0e2f4a5b");

        //Return a string identifying you or your company.
        public override string AuthorName => "Adam Burke";

        //Return a string representing your preferred contact details.
        public override string AuthorContact => "aburke3@mit.edu";
    }
}