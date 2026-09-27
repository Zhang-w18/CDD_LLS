# Plan-037 第 13 节执行手册

本手册只覆盖 `E100_NT32_NR2_V60_BEAM8`。长期统计是一次性波束选择标定；冻结
manifest 后，smoke、prescan 和 formal 只读取该 manifest，不重新选择波束。

## 1. 环境变量

以下 PowerShell 变量仅用于缩短命令，不写入产物：

```powershell
$cfg = 'configs/plan037_cdl_e_32tx_2rx_beam8_statistics.yaml'
$root = 'outputs/experiment037_cdl_e_32tx_beam8'
$gateA = "$root/calibration/gate_a"
$stats = "$root/calibration/statistics"
$freeze = "$root/calibration/freeze"
```

## 2. Gate A：几何、profile 与解析协方差

```powershell
python tools/run_plan037_cdl_e_beam8.py --config $cfg --stage gate-a --output $gateA
```

只有 `$gateA/gate_a_report.json` 的 `status` 为 `PASS` 才进入长期统计。

Gate A 和当前累计统计完成后，可生成角度功率谱及解析/Monte Carlo 波束功率热力图：

```powershell
python tools/analyze_plan037_cdl_e_beam8.py --stage calibration --gate-a-dir $gateA --manifest "$freeze/beam8_frozen_manifest.json" --output "$root/calibration/plots"
```

图中的 `*` 表示当前 manifest 选中的 top-8；manifest 尚为
`NEEDS_MORE_STATISTICS` 时只表示当前累计样本下的临时选择。

## 3. Gate B：一次性长期统计与冻结

初始 2000 个 realization 分四批运行；每批使用绝对且不重叠的区间：

```powershell
python tools/run_plan037_cdl_e_beam8.py --config $cfg --stage statistics --gate-a-dir $gateA --output $stats --start 0 --count 500
python tools/run_plan037_cdl_e_beam8.py --config $cfg --stage statistics --gate-a-dir $gateA --output $stats --start 500 --count 500
python tools/run_plan037_cdl_e_beam8.py --config $cfg --stage statistics --gate-a-dir $gateA --output $stats --start 1000 --count 500
python tools/run_plan037_cdl_e_beam8.py --config $cfg --stage statistics --gate-a-dir $gateA --output $stats --start 1500 --count 500
python tools/run_plan037_cdl_e_beam8.py --config $cfg --stage freeze --gate-a-dir $gateA --statistics-dir $stats --output $freeze --bootstrap-repeats 2000
```

若状态为 `NEEDS_MORE_STATISTICS`，继续以 500 为批次追加，例如：

```powershell
python tools/run_plan037_cdl_e_beam8.py --config $cfg --stage statistics --gate-a-dir $gateA --output $stats --start 2000 --count 500
python tools/run_plan037_cdl_e_beam8.py --config $cfg --stage freeze --gate-a-dir $gateA --statistics-dir $stats --output $freeze --bootstrap-repeats 2000
```

每次追加的 `--start` 等于当前累计 $D$。最大允许 $D=10000$。出现
`TOP8_NOT_IDENTIFIABLE` 或 `CDL_ANGLE_POWER_MISMATCH` 时停止，不进入 BLER。

`beam8_frozen_manifest.json` 的 `status` 为 `FROZEN` 后，不重新运行波束扫描；直接复用该
top-8 manifest 完成第 13.5 节 PDP 提取。此步骤会生成包含 reference PDP/hash 的新 manifest：

```powershell
$pdp = "$root/calibration/pdp"
python tools/run_plan037_cdl_e_beam8.py --config $cfg --stage pdp --manifest "$freeze/beam8_frozen_manifest.json" --output $pdp
python tools/analyze_plan037_cdl_e_beam8.py --stage pdp --pdp-report "$pdp/pdp_report.json" --output "$pdp/plots"
```

只有 `$pdp/pdp_report.json` 为 `PASS` 才从 PDP 增强后的 manifest 生成链路配置：

```powershell
python tools/run_plan037_cdl_e_beam8.py --config $cfg --stage render-configs --manifest "$pdp/beam8_frozen_manifest.json" --output configs
```

生成的 smoke/prescan YAML 含 manifest 的绝对路径和 SHA-256。只要冻结场景定义不变，
这份 manifest 可以复用；修改天面、映射、profile、频率、速度、姿态、资源采样或统计 seed
后必须重新执行 Gate A/B。

## 4. Gate C：链路 validate 与 smoke

```powershell
$smokeCfg = 'configs/plan037_cdl_e_32tx_2rx_beam8_smoke.yaml'
python tools/run_plan037_cdl_e_beam8.py --config $smokeCfg --stage validate
python tools/run_plan037_cdl_e_beam8.py --config $smokeCfg --stage run
```

检查 smoke 输出中三方案各 2 trials，`noise_var` 完全相同，absolute trial、
channel realization、payload/noise seed 配对一致，且 precoder 功率误差不超过 $10^{-12}$。

## 5. Prescan 与 formal 网格

```powershell
$prescanCfg = 'configs/plan037_cdl_e_32tx_2rx_beam8_prescan.yaml'
python tools/run_plan037_cdl_e_beam8.py --config $prescanCfg --stage validate
python tools/run_plan037_cdl_e_beam8.py --config $prescanCfg --stage run
python tools/analyze_plan037_cdl_e_beam8.py --stage prescan --run-dir "$root/prescan/prescan" --output "$root/prescan_analysis" --base-config $prescanCfg --formal-config 'configs/plan037_cdl_e_32tx_2rx_beam8_formal.yaml'
```

若分析报告某个目标缺双侧 bracket，只沿报告方向按 2 dB 扩展，并保持
`[-10,22] dB` 硬边界。向高 SNR 扩展时使用自动生成的
`configs/plan037_cdl_e_32tx_2rx_beam8_prescan_extend22.yaml`，它只运行
`[16,18,20,22] dB`；分析时同时传入基础 prescan 和扩展目录。不得外推。得到 formal YAML
后先人工核对网格和预算，再运行首批：

```powershell
$formalCfg = 'configs/plan037_cdl_e_32tx_2rx_beam8_formal.yaml'
python tools/run_plan037_cdl_e_beam8.py --config $formalCfg --stage validate
python tools/run_plan037_cdl_e_beam8.py --config $formalCfg --stage run
```

## 6. Formal 追加与分析

需要追加时，从当前 formal YAML 生成新的、不重叠的 1000-trial batch。示例为追加
absolute trial `[1000,2000)`：

```powershell
python tools/run_plan037_cdl_e_beam8.py --config $formalCfg --stage render-batch --start 1000 --count 1000 --output 'configs/plan037_cdl_e_32tx_2rx_beam8_formal_01000_02000.yaml'
python tools/run_plan037_cdl_e_beam8.py --config 'configs/plan037_cdl_e_32tx_2rx_beam8_formal_01000_02000.yaml' --stage validate
python tools/run_plan037_cdl_e_beam8.py --config 'configs/plan037_cdl_e_32tx_2rx_beam8_formal_01000_02000.yaml' --stage run
```

分析时可重复给出 `--run-dir`，入口会拒绝 absolute trial 重复、缺口和三方案配对失配：

```powershell
python tools/analyze_plan037_cdl_e_beam8.py --stage formal `
  --run-dir "$root/formal/formal_00000_01000" `
  --run-dir "$root/formal/formal_01000_02000" `
  --output "$root/formal_analysis" --bootstrap-repeats 2000
```

正式停止规则仍以 plan 第 13.9 节为准：10% 与 1% bracket 端点满足相应 BLER 区间时
至少 200 errors，单点最多 50000 trials；crossing 只用真实相邻点的 log-BLER 线性插值。
