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

The native library includes `ariadne-lbfgsb`, a BSD-3-Clause solver. Distribute
all files named `ariadne-lbfgsb-*` with the native binary. The surrounding
Ariadne project is MIT licensed; see `Ariadne-LICENSE.txt`.
