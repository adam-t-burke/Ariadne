<#
.SYNOPSIS
    Builds tetgen_wrapper.dll from TetGen source + C-compatible shim using MSVC.

.DESCRIPTION
    Locates Visual Studio via vswhere, invokes vcvars64.bat to set up the MSVC
    environment, then compiles tetgen_wrapper.cpp + tetgen.cxx + predicates.cxx
    into a shared library (DLL).

.EXAMPLE
    .\build_tetgen.ps1
#>

$ErrorActionPreference = "Stop"

$scriptDir = $PSScriptRoot

# Locate Visual Studio
$vsWhere = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe"
if (-not (Test-Path $vsWhere)) {
    Write-Error "vswhere.exe not found. Install Visual Studio with C++ build tools."
    exit 1
}

$vsPath = & $vsWhere -latest -property installationPath
if (-not $vsPath) {
    Write-Error "No Visual Studio installation found."
    exit 1
}

$vcvars = Join-Path $vsPath "VC\Auxiliary\Build\vcvars64.bat"
if (-not (Test-Path $vcvars)) {
    Write-Error "vcvars64.bat not found at $vcvars. Install the C++ desktop workload."
    exit 1
}

Write-Host "Building tetgen_wrapper.dll..." -ForegroundColor Cyan
Write-Host "Using Visual Studio at: $vsPath" -ForegroundColor Gray

# Build via cmd.exe so vcvars64.bat can set environment variables
$buildCmd = @"
call "$vcvars" >nul 2>&1 && ^
cl /O2 /DTETLIBRARY /DNDEBUG /EHsc /LD /nologo ^
   /Fe:"$scriptDir\tetgen_wrapper.dll" ^
   /Fo:"$scriptDir\\" ^
   "$scriptDir\tetgen_wrapper.cpp" ^
   "$scriptDir\tetgen.cxx" ^
   "$scriptDir\predicates.cxx"
"@

cmd /c $buildCmd
if ($LASTEXITCODE -ne 0) {
    Write-Error "Compilation failed with exit code $LASTEXITCODE"
    exit $LASTEXITCODE
}

# Clean up intermediate files
Remove-Item "$scriptDir\*.obj" -ErrorAction SilentlyContinue
Remove-Item "$scriptDir\*.exp" -ErrorAction SilentlyContinue
Remove-Item "$scriptDir\tetgen_wrapper.lib" -ErrorAction SilentlyContinue

if (Test-Path "$scriptDir\tetgen_wrapper.dll") {
    $size = (Get-Item "$scriptDir\tetgen_wrapper.dll").Length / 1KB
    Write-Host "Built tetgen_wrapper.dll ($([math]::Round($size)) KB)" -ForegroundColor Green
} else {
    Write-Error "tetgen_wrapper.dll was not produced"
    exit 1
}

Write-Host "Done." -ForegroundColor Green
