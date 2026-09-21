param(
    [int]$StartIndex = 0
)

$ErrorActionPreference = "Continue"
$repo = Split-Path -Parent (Split-Path -Parent $MyInvocation.MyCommand.Path)
$python = "D:\venvs\cdd-s102\Scripts\python.exe"
$runner = Join-Path $repo "tools\run_bler_curves.py"
$output = Join-Path $repo "outputs\experiment028_csi_curves\20260811_4rx_subset\formal"
$log = Join-Path $output "background_runner.log"
$status = Join-Path $output "background_runner_status.json"
$env:MPLCONFIGDIR = Join-Path $repo ".mplconfig"
$configs = @(
    "configs\bler_curves_result028_a100_4rx_subset_formal_cdd.yaml",
    "configs\bler_curves_result028_a100_4rx_subset_formal_prg6.yaml",
    "configs\bler_curves_result028_a100_4rx_subset_formal_small_delay_estimated.yaml",
    "configs\bler_curves_result028_a100_4rx_subset_formal_small_delay_ideal.yaml"
)

New-Item -ItemType Directory -Force -Path $output | Out-Null
for ($index = $StartIndex; $index -lt $configs.Count; $index++) {
    $config = $configs[$index]
    @{
        status = "running"
        config_index = $index
        config = $config
        updated_at = (Get-Date).ToString("o")
    } | ConvertTo-Json | Set-Content -LiteralPath $status -Encoding UTF8
    & $python $runner --config (Join-Path $repo $config) --stage all *>> $log
    if ($LASTEXITCODE -ne 0) {
        @{
            status = "failed"
            config_index = $index
            config = $config
            exit_code = $LASTEXITCODE
            updated_at = (Get-Date).ToString("o")
        } | ConvertTo-Json | Set-Content -LiteralPath $status -Encoding UTF8
        exit $LASTEXITCODE
    }
}

@{
    status = "complete"
    config_index = $configs.Count
    updated_at = (Get-Date).ToString("o")
} | ConvertTo-Json | Set-Content -LiteralPath $status -Encoding UTF8
exit 0
