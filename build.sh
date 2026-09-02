#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
WORKSPACE_DIR="$SCRIPT_DIR"
RUST_DIR="$SCRIPT_DIR/crates/theseus"
OUTPUT_DIR="$SCRIPT_DIR/Theseus"

CONFIGURATION="release"
SKIP_TESTS=0
for argument in "$@"; do
    case "$argument" in
        release|debug)
            CONFIGURATION="$argument"
            ;;
        --skip-tests|--fast)
            SKIP_TESTS=1
            ;;
        *)
            echo "Usage: $0 [release|debug] [--skip-tests|--fast]" >&2
            exit 2
            ;;
    esac
done

if [ ! -d "$RUST_DIR" ]; then
    echo "Rust source not found at $RUST_DIR." >&2
    exit 1
fi
for command in cargo rustup lipo file; do
    if ! command -v "$command" >/dev/null 2>&1; then
        echo "Required command '$command' was not found on PATH." >&2
        exit 1
    fi
done

PROFILE_FLAGS=()
TARGET_SUBDIR="debug"
if [ "$CONFIGURATION" = "release" ]; then
    PROFILE_FLAGS=(--release)
    TARGET_SUBDIR="release"
fi

if [ "$SKIP_TESTS" -eq 0 ]; then
    echo "Testing Rust workspace..."
    cargo test --workspace --manifest-path "$WORKSPACE_DIR/Cargo.toml"
else
    echo "Skipping Rust tests (--skip-tests)."
fi

echo "Building theseus ($CONFIGURATION) universal binary..."

rustup target add aarch64-apple-darwin x86_64-apple-darwin

cargo build -p theseus "${PROFILE_FLAGS[@]}" --target aarch64-apple-darwin --manifest-path "$WORKSPACE_DIR/Cargo.toml"
cargo build -p theseus "${PROFILE_FLAGS[@]}" --target x86_64-apple-darwin --manifest-path "$WORKSPACE_DIR/Cargo.toml"

ARM64="$WORKSPACE_DIR/target/aarch64-apple-darwin/$TARGET_SUBDIR/libtheseus.dylib"
X86_64="$WORKSPACE_DIR/target/x86_64-apple-darwin/$TARGET_SUBDIR/libtheseus.dylib"

if [ ! -f "$ARM64" ] || [ ! -f "$X86_64" ]; then
    echo "Expected architecture-specific dylibs were not produced." >&2
    exit 1
fi

lipo -create -output "$OUTPUT_DIR/libtheseus.dylib" "$ARM64" "$X86_64"

NOTICE_SOURCES=(
    "$WORKSPACE_DIR/LICENSE.txt"
    "$WORKSPACE_DIR/crates/lbfgsb/LICENSE"
    "$WORKSPACE_DIR/crates/lbfgsb/UPSTREAM_LICENSE.txt"
    "$WORKSPACE_DIR/crates/lbfgsb/THIRD_PARTY_NOTICES.md"
)
NOTICE_DESTINATIONS=(
    "$OUTPUT_DIR/Ariadne-LICENSE.txt"
    "$OUTPUT_DIR/ariadne-lbfgsb-LICENSE.txt"
    "$OUTPUT_DIR/ariadne-lbfgsb-UPSTREAM_LICENSE.txt"
    "$OUTPUT_DIR/ariadne-lbfgsb-THIRD_PARTY_NOTICES.md"
)
for index in "${!NOTICE_SOURCES[@]}"; do
    if [ ! -f "${NOTICE_SOURCES[$index]}" ]; then
        echo "Required notice '${NOTICE_SOURCES[$index]}' was not found." >&2
        exit 1
    fi
    cp "${NOTICE_SOURCES[$index]}" "${NOTICE_DESTINATIONS[$index]}"
done

echo "Created universal binary at $OUTPUT_DIR/libtheseus.dylib"
echo "Copied native license and provenance bundle into $OUTPUT_DIR"
file "$OUTPUT_DIR/libtheseus.dylib"
