#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
RUST_DIR="$SCRIPT_DIR/Theseus/rust"
OUTPUT_DIR="$SCRIPT_DIR/Theseus"

CONFIGURATION="${1:-release}"
PROFILE_FLAG=""
TARGET_SUBDIR="debug"
if [ "$CONFIGURATION" = "release" ]; then
    PROFILE_FLAG="--release"
    TARGET_SUBDIR="release"
fi

echo "Building theseus ($CONFIGURATION) universal binary..."

rustup target add aarch64-apple-darwin x86_64-apple-darwin

cargo build $PROFILE_FLAG --target aarch64-apple-darwin --manifest-path "$RUST_DIR/Cargo.toml"
cargo build $PROFILE_FLAG --target x86_64-apple-darwin  --manifest-path "$RUST_DIR/Cargo.toml"

ARM64="$RUST_DIR/target/aarch64-apple-darwin/$TARGET_SUBDIR/libtheseus.dylib"
X86_64="$RUST_DIR/target/x86_64-apple-darwin/$TARGET_SUBDIR/libtheseus.dylib"

lipo -create -output "$OUTPUT_DIR/libtheseus.dylib" "$ARM64" "$X86_64"

echo "Created universal binary at $OUTPUT_DIR/libtheseus.dylib"
file "$OUTPUT_DIR/libtheseus.dylib"
