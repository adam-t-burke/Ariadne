# SPDX-License-Identifier: BSD-3-Clause
[CmdletBinding()]
param(
    [int]$Repeats = 20
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$sourceUrl = 'http://users.iems.northwestern.edu/~nocedal/Software/Lbfgsb.3.0.tar.gz'
$expectedSha256 = 'f5b9a1c8c30ff6bcc8df9b5d5738145f4cbe4c7eadec629220e808dcf0e54720'
$workDir = Join-Path ([System.IO.Path]::GetTempPath()) "ariadne-lbfgsb-bench-$([guid]::NewGuid())"
$archive = Join-Path $workDir 'Lbfgsb.3.0.tar.gz'
$sourceDir = Join-Path $workDir 'source'
$executable = Join-Path $workDir 'benchmark_driver.exe'

try {
    New-Item -ItemType Directory -Path $workDir, $sourceDir -Force | Out-Null
    Invoke-WebRequest -Uri $sourceUrl -OutFile $archive
    $actual = (Get-FileHash -Path $archive -Algorithm SHA256).Hash.ToLowerInvariant()
    if ($actual -ne $expectedSha256) {
        throw "Archive SHA-256 mismatch: expected $expectedSha256, got $actual"
    }
    tar -xzf $archive -C $sourceDir
    if ($LASTEXITCODE -ne 0) { throw "tar failed with exit code $LASTEXITCODE" }

    $upstream = Join-Path $sourceDir 'Lbfgsb.3.0'
    $sources = @(
        (Join-Path $upstream 'lbfgsb.f'),
        (Join-Path $upstream 'linpack.f'),
        (Join-Path $upstream 'blas.f'),
        (Join-Path $upstream 'timer.f'),
        (Join-Path $PSScriptRoot 'benchmark_driver.f90')
    )
    & gfortran -O3 -std=legacy -ffixed-line-length-none -ffree-line-length-none -fno-fast-math @sources -o $executable
    if ($LASTEXITCODE -ne 0) { throw "gfortran failed with exit code $LASTEXITCODE" }

    Write-Host (& gfortran --version | Select-Object -First 1)
    foreach ($case in @(@(25, 5), @(1000, 10))) {
        & $executable $case[0] $case[1] $Repeats
        if ($LASTEXITCODE -ne 0) { throw "benchmark failed with exit code $LASTEXITCODE" }
    }
}
finally {
    if (Test-Path $workDir) {
        Remove-Item -Path $workDir -Recurse -Force
    }
}
