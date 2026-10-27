<#
.SYNOPSIS
    Build an ariadne-dev .yak package from the Release .gha output.

.DESCRIPTION
    Copies the Grasshopper distribution bundle into packaging/ariadne-dev/stage,
    runs Yak build, and optionally pushes to the Package Manager.

    Prerequisites:
      - .\build.ps1 (native theseus.dll + notices)
      - dotnet build Ariadne.csproj -c Release
      - Rhino 8 Yak CLI at the default install path (or set $env:YAK_EXE)

.EXAMPLE
    .\packaging\Pack-Yak.ps1
    .\packaging\Pack-Yak.ps1 -Push
#>
param(
    [switch]$Push,
    [switch]$SkipBuild
)

$ErrorActionPreference = "Stop"

$repoRoot = Split-Path $PSScriptRoot -Parent
if (-not (Test-Path (Join-Path $repoRoot "Ariadne.csproj"))) {
    $repoRoot = $PSScriptRoot
    if (-not (Test-Path (Join-Path $repoRoot "Ariadne.csproj"))) {
        throw "Could not locate repo root containing Ariadne.csproj."
    }
}

$yakExe = if ($env:YAK_EXE) {
    $env:YAK_EXE
} else {
    "C:\Program Files\Rhino 8\System\Yak.exe"
}
if (-not (Test-Path $yakExe)) {
    throw "Yak CLI not found at '$yakExe'. Install Rhino 8 or set YAK_EXE."
}

# Dedicated output avoids file locks when Rhino has bin\Release\...\Ariadne.gha open.
$output = Join-Path $PSScriptRoot "ariadne-dev\build-out"

if (-not $SkipBuild) {
    Write-Host "Building native theseus..." -ForegroundColor Cyan
    & (Join-Path $repoRoot "build.ps1") -SkipTests
    if ($LASTEXITCODE -ne 0) { throw "build.ps1 failed." }

    Write-Host "Building Ariadne (Release) into $output ..." -ForegroundColor Cyan
    if (Test-Path $output) {
        Remove-Item $output -Recurse -Force
    }
    Push-Location $repoRoot
    try {
        dotnet build Ariadne.csproj -c Release -o $output
        if ($LASTEXITCODE -ne 0) { throw "dotnet build failed." }
    } finally {
        Pop-Location
    }
}

$gha = Join-Path $output "Ariadne.gha"
if (-not (Test-Path $gha)) {
    throw "Missing $gha — build Release first (or omit -SkipBuild)."
}

$stage = Join-Path $PSScriptRoot "ariadne-dev\stage"
if (Test-Path $stage) {
    Remove-Item $stage -Recurse -Force
}
New-Item -ItemType Directory -Path $stage | Out-Null

$required = @(
    "Ariadne.gha",
    "theseus.dll",
    "Ariadne-LICENSE.txt",
    "ariadne-lbfgsb-LICENSE.txt",
    "ariadne-lbfgsb-UPSTREAM_LICENSE.txt",
    "ariadne-lbfgsb-THIRD_PARTY_NOTICES.md",
    "NATIVE_ARTIFACTS.md"
)

foreach ($name in $required) {
    $src = Join-Path $output $name
    if (-not (Test-Path $src)) {
        throw "Missing distribution file: $src"
    }
    Copy-Item $src (Join-Path $stage $name) -Force
}

$dylib = Join-Path $output "libtheseus.dylib"
if (Test-Path $dylib) {
    Copy-Item $dylib (Join-Path $stage "libtheseus.dylib") -Force
}

Copy-Item (Join-Path $PSScriptRoot "ariadne-dev\manifest.yml") (Join-Path $stage "manifest.yml") -Force

Write-Host "Building yak package from $stage ..." -ForegroundColor Cyan
Push-Location $stage
try {
    Get-ChildItem -Filter "*.yak" | Remove-Item -Force
    & $yakExe build
    if ($LASTEXITCODE -ne 0) { throw "yak build failed." }

    $yak = Get-ChildItem -Filter "ariadne-dev-*.yak" | Select-Object -First 1
    if (-not $yak) {
        throw "yak build did not produce an ariadne-dev-*.yak file."
    }

    $destDir = Join-Path $PSScriptRoot "ariadne-dev"
    $dest = Join-Path $destDir $yak.Name
    Copy-Item $yak.FullName $dest -Force
    Write-Host "Created $dest" -ForegroundColor Green

    if ($Push) {
        Write-Host "Pushing $($yak.Name) to yak.rhino3d.com ..." -ForegroundColor Cyan
        & $yakExe push $yak.FullName
        if ($LASTEXITCODE -ne 0) { throw "yak push failed (are you logged in? run: yak login)." }
        Write-Host "Pushed. Search Package Manager for 'ariadne-dev'." -ForegroundColor Green
        & $yakExe search --all --prerelease ariadne-dev
    } else {
        Write-Host "Package ready. Push with: .\packaging\Pack-Yak.ps1 -SkipBuild -Push" -ForegroundColor Yellow
        Write-Host "Or: & '$yakExe' push '$dest'" -ForegroundColor Yellow
    }
} finally {
    Pop-Location
}
