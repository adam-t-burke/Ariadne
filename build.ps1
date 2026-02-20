<#
.SYNOPSIS
    Builds the Rust theseus native library and copies it into the Theseus/ folder.

.DESCRIPTION
    Runs 'cargo build --release' in Theseus/rust/ and copies the resulting
    theseus.dll to Theseus/theseus.dll so the .NET build can pick it up.

    Requires the Rust toolchain (rustup / cargo) to be installed and on PATH.

.EXAMPLE
    .\build.ps1
    .\build.ps1 -Configuration Debug
#>
param(
    [ValidateSet("Release", "Debug")]
    [string]$Configuration = "Release"
)

$ErrorActionPreference = "Stop"

$rustDir = Join-Path (Join-Path $PSScriptRoot "Theseus") "rust"
$outputDir = Join-Path $PSScriptRoot "Theseus"

if (-not (Test-Path $rustDir)) {
    Write-Error "Rust source not found at $rustDir. Copy your Rust project into Theseus/rust/ first."
    exit 1
}

Write-Host "Building theseus ($Configuration)..." -ForegroundColor Cyan

$profile = if ($Configuration -eq "Release") { "--release" } else { "" }
$targetSubdir = if ($Configuration -eq "Release") { "release" } else { "debug" }

Push-Location $rustDir
try {
    if ($profile) {
        cargo build $profile
    } else {
        cargo build
    }
    if ($LASTEXITCODE -ne 0) {
        Write-Error "Cargo build failed with exit code $LASTEXITCODE"
        exit $LASTEXITCODE
    }
} finally {
    Pop-Location
}

$dllSource = Join-Path (Join-Path (Join-Path $rustDir "target") $targetSubdir) "theseus.dll"
$dllDest = Join-Path $outputDir "theseus.dll"

if (-not (Test-Path $dllSource)) {
    Write-Error "Expected DLL not found at $dllSource"
    exit 1
}

Copy-Item $dllSource $dllDest -Force
Write-Host "Copied $dllSource -> $dllDest" -ForegroundColor Green
Write-Host "Done." -ForegroundColor Green
