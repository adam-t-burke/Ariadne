# Native artifact provenance

The tracked `theseus.dll` and `libtheseus.dylib` are source-build bootstrap
artifacts for IDE usability. They are not authoritative release artifacts and
their repository timestamps or hashes must not be presented as reproducible
provenance.

Release binaries are built from the checked-out Rust workspace:

- Windows: `./build.ps1 -Configuration Release`
- macOS universal binary: `./build.sh release`

Each script replaces the matching bootstrap binary and copies the complete
license/attribution bundle beside it. CI independently rebuilds each platform,
checks the native file type and notice bundle, and uploads the rebuilt binary
and notices as a platform-specific artifact. Releases should use those CI
artifacts, not a pre-existing binary from the source tree.

Rebuilt native libraries use Basin for L-BFGS and L-BFGS-B. Distribute
`basin-LICENSE-MIT.txt`, `basin-LICENSE-APACHE.txt`, and `basin-COPYRIGHT.txt`
with the binary. The source copies and their version are recorded in
[`third-party/basin`](../third-party/basin/README.md).

The tracked bootstrap binaries predate the Basin migration and may include
`ariadne-lbfgsb`. Its `ariadne-lbfgsb-*` notices remain in the distribution bundle
to cover those binaries. The standalone crate also retains its own notices.
The surrounding Ariadne project is MIT licensed; see `Ariadne-LICENSE.txt`.
