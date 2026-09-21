param(
    [string]$Config = "configs/bler_curves_result028_a100_transparent_cdd_small_delay.yaml",
    [int]$MaxInvocations = 500
)

$ErrorActionPreference = "Continue"
$repo = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$python = "D:\venvs\cdd-s102\Scripts\python.exe"
$runner = Join-Path $repo "tools\run_bler_curves.py"
$configPath = Join-Path $repo $Config
$output = Join-Path $repo "outputs\experiment028_csi_curves\20260803_main\transparent_cdd_small_delay"
$progress = Join-Path $output "a100\run_progress.jsonl"
$log = Join-Path $output "background_runner.log"
$status = Join-Path $output "background_runner_status.json"
$env:MPLCONFIGDIR = Join-Path $repo ".mplconfig"
New-Item -ItemType Directory -Force -Path $output | Out-Null
$initialProgress = if (Test-Path -LiteralPath $progress) { (Get-Item -LiteralPath $progress).Length } else { 0 }
@{status="running"; invocation=0; progress_bytes=$initialProgress} |
    ConvertTo-Json | Set-Content -LiteralPath $status -Encoding UTF8

for ($invocation = 1; $invocation -le $MaxInvocations; $invocation++) {
    $before = if (Test-Path -LiteralPath $progress) { (Get-Item -LiteralPath $progress).Length } else { 0 }
    & $python $runner --config $configPath --stage run *>> $log
    if ($LASTEXITCODE -ne 0) {
        @{status="failed"; invocation=$invocation; exit_code=$LASTEXITCODE} |
            ConvertTo-Json | Set-Content -LiteralPath $status -Encoding UTF8
        exit $LASTEXITCODE
    }
    $after = if (Test-Path -LiteralPath $progress) { (Get-Item -LiteralPath $progress).Length } else { 0 }
    @{status="running"; invocation=$invocation; progress_bytes=$after} |
        ConvertTo-Json | Set-Content -LiteralPath $status -Encoding UTF8
    if ($after -eq $before) {
        @{status="complete"; invocation=$invocation; progress_bytes=$after} |
            ConvertTo-Json | Set-Content -LiteralPath $status -Encoding UTF8
        exit 0
    }
}

@{status="max_invocations_reached"; invocation=$MaxInvocations} |
    ConvertTo-Json | Set-Content -LiteralPath $status -Encoding UTF8
exit 2
