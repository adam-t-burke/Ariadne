<#
.SYNOPSIS
    Builds the Rust native library and its redistribution notice bundle.

.DESCRIPTION
    Runs the Rust workspace tests, builds the selected theseus profile, and
    replaces Theseus/theseus.dll so the .NET build can pick it up. Copies the
    MIT, BSD-3-Clause, verbatim upstream, and attribution notices beside it.

    Requires the Rust toolchain (rustup / cargo) to be installed and on PATH.

.EXAMPLE
    .\build.ps1
    .\build.ps1 -Configuration Debug
    .\build.ps1 -SkipTests
#>
param(
    [ValidateSet("Release", "Debug")]
    [string]$Configuration = "Release",

    [Alias("Fast")]
    [switch]$SkipTests
)

$ErrorActionPreference = "Stop"

$workspaceDir = $PSScriptRoot
$rustDir = Join-Path (Join-Path $PSScriptRoot "crates") "theseus"
$outputDir = Join-Path $PSScriptRoot "Theseus"

if (-not (Test-Path $rustDir)) {
    throw "Rust source not found at $rustDir."
}

if (-not (Get-Command cargo -ErrorAction SilentlyContinue)) {
    throw "Cargo was not found on PATH. Install a Rust toolchain before running this script."
}

$profile = if ($Configuration -eq "Release") { "--release" } else { "" }
$targetSubdir = if ($Configuration -eq "Release") { "release" } else { "debug" }

Push-Location $workspaceDir
try {
    if (-not $SkipTests) {
        Write-Host "Testing Rust workspace..." -ForegroundColor Cyan
        cargo test --workspace
        if ($LASTEXITCODE -ne 0) {
            throw "Cargo tests failed with exit code $LASTEXITCODE."
        }
    } else {
        Write-Host "Skipping Rust tests (-SkipTests)." -ForegroundColor Yellow
    }

    Write-Host "Building theseus ($Configuration)..." -ForegroundColor Cyan
    if ($profile) {
        cargo build -p theseus $profile
    } else {
        cargo build -p theseus
    }
    if ($LASTEXITCODE -ne 0) {
        throw "Cargo build failed with exit code $LASTEXITCODE."
    }
} finally {
    Pop-Location
}

$dllSource = Join-Path (Join-Path (Join-Path $workspaceDir "target") $targetSubdir) "theseus.dll"
$dllDest = Join-Path $outputDir "theseus.dll"

if (-not (Test-Path $dllSource)) {
    throw "Expected DLL not found at $dllSource."
}

Copy-Item $dllSource $dllDest -Force
Write-Host "Copied $dllSource -> $dllDest" -ForegroundColor Green
$noticeFiles = @(
    @{ Source = (Join-Path $workspaceDir "LICENSE.txt"); Destination = "Ariadne-LICENSE.txt" },
    @{ Source = (Join-Path $workspaceDir "crates/lbfgsb/LICENSE"); Destination = "ariadne-lbfgsb-LICENSE.txt" },
    @{ Source = (Join-Path $workspaceDir "crates/lbfgsb/UPSTREAM_LICENSE.txt"); Destination = "ariadne-lbfgsb-UPSTREAM_LICENSE.txt" },
    @{ Source = (Join-Path $workspaceDir "crates/lbfgsb/THIRD_PARTY_NOTICES.md"); Destination = "ariadne-lbfgsb-THIRD_PARTY_NOTICES.md" }
)
foreach ($notice in $noticeFiles) {
    if (-not (Test-Path $notice.Source)) {
        throw "Required notice not found at $($notice.Source)."
    }
    Copy-Item $notice.Source (Join-Path $outputDir $notice.Destination) -Force
}
Write-Host "Copied native license and provenance bundle into $outputDir" -ForegroundColor Green
Write-Host "Done." -ForegroundColor Green
