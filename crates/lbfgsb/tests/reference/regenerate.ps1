# SPDX-License-Identifier: BSD-3-Clause
[CmdletBinding()]
param()

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$sourceUrl = 'http://users.iems.northwestern.edu/~nocedal/Software/Lbfgsb.3.0.tar.gz'
$expectedSha256 = 'f5b9a1c8c30ff6bcc8df9b5d5738145f4cbe4c7eadec629220e808dcf0e54720'
$referenceDir = $PSScriptRoot
$fixturesDir = Join-Path $referenceDir 'fixtures'
$workDir = Join-Path ([System.IO.Path]::GetTempPath()) "ariadne-lbfgsb-$([guid]::NewGuid())"
$archive = Join-Path $workDir 'Lbfgsb.3.0.tar.gz'
$extractDir = Join-Path $workDir 'source'
$generatedDir = Join-Path $workDir 'generated'
$backupDir = Join-Path $workDir 'backup'
$executable = Join-Path $workDir 'trace_generator.exe'
$expectedFixtures = @(
    'driver1.csv',
    'driver2.csv',
    'driver3-large-n.csv',
    'edge-mixed-bounds.csv',
    'audit-mixed-rollover.csv',
    'dcsrch-cases.csv'
)
$utf8NoBom = [System.Text.UTF8Encoding]::new($false)

try {
    New-Item -ItemType Directory -Path $workDir, $extractDir, $generatedDir, $backupDir -Force |
        Out-Null

    Write-Host "Downloading authoritative L-BFGS-B 3.0 source..."
    Invoke-WebRequest -Uri $sourceUrl -OutFile $archive
    $actualSha256 = (Get-FileHash -Path $archive -Algorithm SHA256).Hash.ToLowerInvariant()
    if ($actualSha256 -ne $expectedSha256) {
        throw "Archive SHA-256 mismatch: expected $expectedSha256, got $actualSha256"
    }

    tar -xzf $archive -C $extractDir
    if ($LASTEXITCODE -ne 0) {
        throw "tar failed with exit code $LASTEXITCODE"
    }

    $upstream = Join-Path $extractDir 'Lbfgsb.3.0'
    $compilerVersion = (& gfortran --version | Select-Object -First 1)
    if ($LASTEXITCODE -ne 0) {
        throw 'gfortran is required to regenerate the reference fixtures'
    }
    Write-Host "Compiler: $compilerVersion"

    $flags = @(
        '-O0',
        '-std=legacy',
        '-ffixed-line-length-none',
        '-ffree-line-length-none',
        '-fno-fast-math'
    )
    $sources = @(
        (Join-Path $upstream 'lbfgsb.f'),
        (Join-Path $upstream 'linpack.f'),
        (Join-Path $upstream 'blas.f'),
        (Join-Path $upstream 'timer.f'),
        (Join-Path $referenceDir 'trace_generator.f90')
    )
    & gfortran @flags @sources '-o' $executable
    if ($LASTEXITCODE -ne 0) {
        throw "gfortran failed with exit code $LASTEXITCODE"
    }

    & $executable $generatedDir
    if ($LASTEXITCODE -ne 0) {
        throw "fixture generator failed with exit code $LASTEXITCODE"
    }

    $actualFixtures = @(Get-ChildItem -Path $generatedDir -File -Filter '*.csv' |
        Select-Object -ExpandProperty Name |
        Sort-Object)
    if (Compare-Object ($expectedFixtures | Sort-Object) $actualFixtures) {
        throw "fixture generator produced an unexpected file set: $($actualFixtures -join ', ')"
    }

    $hashes = @{}
    foreach ($name in $expectedFixtures) {
        $path = Join-Path $generatedDir $name
        $text = [System.IO.File]::ReadAllText($path)
        $canonical = $text.Replace("`r`n", "`n").Replace("`r", "`n")
        if (-not $canonical.EndsWith("`n")) {
            $canonical += "`n"
        }
        [System.IO.File]::WriteAllText($path, $canonical, $utf8NoBom)

        $firstLine = ([System.IO.File]::ReadAllLines($path))[0]
        $expectedHeader = if ($name -eq 'dcsrch-cases.csv') {
            'case,step,task'
        }
        else {
            'record,iteration,evaluations,f,projected_gradient,x,g,task'
        }
        if ($firstLine -ne $expectedHeader) {
            throw "fixture $name has an unexpected header"
        }
        $hashes[$name] = (Get-FileHash -Path $path -Algorithm SHA256).Hash.ToLowerInvariant()
    }

    # All generated output is validated before replacing any committed fixture.
    # Backups permit rollback if a replacement fails partway through.
    New-Item -ItemType Directory -Path $fixturesDir -Force | Out-Null
    foreach ($name in $expectedFixtures) {
        $destination = Join-Path $fixturesDir $name
        if (Test-Path $destination) {
            Copy-Item $destination (Join-Path $backupDir $name)
        }
    }
    try {
        foreach ($name in $expectedFixtures) {
            Copy-Item (Join-Path $generatedDir $name) (Join-Path $fixturesDir $name) -Force
        }
    }
    catch {
        foreach ($name in $expectedFixtures) {
            $backup = Join-Path $backupDir $name
            if (Test-Path $backup) {
                Copy-Item $backup (Join-Path $fixturesDir $name) -Force
            }
        }
        throw
    }

    Write-Host "Replaced validated LF fixtures in $fixturesDir"
    foreach ($name in ($expectedFixtures | Sort-Object)) {
        Write-Host "$name  $($hashes[$name])"
    }
}
finally {
    if (Test-Path $workDir) {
        Remove-Item -Path $workDir -Recurse -Force
    }
}
