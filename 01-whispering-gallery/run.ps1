#!/usr/bin/env pwsh
# ── Run this Kaimon Slate notebook live (Windows / PowerShell) ───────────────────────────────
# Auto-generated. Runs the sibling run.jl, which installs Kaimon + KaimonSlate into a
# dedicated environment, fetches this notebook's reproducible bundle, and serves it — the
# notebook's exact environment reconstructs on open. Prerequisite: Julia 1.10+ (juliaup /
# https://julialang.org/downloads). Double-click run.bat, or right-click this file →
# "Run with PowerShell".
$ErrorActionPreference = "Stop"
Set-Location -Path $PSScriptRoot
# Keep the window open on exit so a message is readable — UNLESS run.bat launched us (it pauses
# itself, so we'd otherwise double-prompt). The bat sets SLATE_LAUNCHED_BY_BAT before calling.
function Hold { if (-not $env:SLATE_LAUNCHED_BY_BAT) { Read-Host "`nPress Enter to close" } }
if (-not (Get-Command julia -ErrorAction SilentlyContinue)) {
    Write-Host "Julia was not found on your PATH." -ForegroundColor Yellow
    Write-Host "Install it from https://julialang.org/downloads/  (or run:  winget install julia )"
    Hold; exit 1
}
julia (Join-Path $PSScriptRoot "run.jl")
Hold
