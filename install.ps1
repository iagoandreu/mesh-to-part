# Push this folder into Fusion's AddIns directory.
# MeshToPart.py is the shell - changing it needs the add-in stopped and re-run.
# mtp_core.py hot-reloads, so pushing only that file takes effect on the next command.
$src = Split-Path -Parent $MyInvocation.MyCommand.Definition
$dst = Join-Path $env:APPDATA 'Autodesk\Autodesk Fusion 360\API\AddIns\MeshToPart'

New-Item -ItemType Directory -Force -Path $dst | Out-Null

$cache = Join-Path $dst '__pycache__'
if (Test-Path -LiteralPath $cache) { Remove-Item -LiteralPath $cache -Recurse -Force }

foreach ($f in @('MeshToPart.py', 'mtp_core.py', 'mtp_organic.py', 'MeshToPart.manifest')) {
    Copy-Item -LiteralPath (Join-Path $src $f) -Destination $dst -Force
}
# Copy-Item -Recurse into a directory that ALREADY holds one of the same name
# nests it (resources/Wrap/Wrap) instead of merging, and the old icons then stay
# where Fusion is looking. Clear the destination first so the copy is a replace.
$resDst = Join-Path $dst 'resources'
if (Test-Path -LiteralPath $resDst) { Remove-Item -LiteralPath $resDst -Recurse -Force }
Copy-Item -LiteralPath (Join-Path $src 'resources') -Destination $dst -Recurse -Force

# the dev bridge is opt-in: the marker file is only copied if it exists here
$flag = Join-Path $src 'devbridge.on'
$dstFlag = Join-Path $dst 'devbridge.on'
if (Test-Path -LiteralPath $flag) {
    Copy-Item -LiteralPath $flag -Destination $dst -Force
} elseif (Test-Path -LiteralPath $dstFlag) {
    Remove-Item -LiteralPath $dstFlag -Force
}

Write-Host "Installed to $dst"
if (Test-Path -LiteralPath $dstFlag) { Write-Host "dev bridge: ON" } else { Write-Host "dev bridge: off" }
