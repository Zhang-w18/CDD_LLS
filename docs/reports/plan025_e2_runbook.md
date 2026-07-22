# plan-025 E2 异机执行说明

## 1. 执行范围

只执行 QC 与 Sidon 的 known-delay matched/oracle 二维 RMMSE 链路对比。不运行 unknown-delay BLER。

固定条件：Sionna 1.0.2 TDL-A、RMS delay spread 5 ns、UE 速度 0 km/h、载频 3.5 GHz、8 Tx / 1 Rx、48 RB × 10 OFDM symbols、DMRS symbols `[2,7]`、DMRS comb 24、16QAM、MCS 8、LDPC 8 次迭代。

## 2. 环境

使用 Python 3.11，并在短路径建立虚拟环境。Windows 示例：

```powershell
& "$env:LOCALAPPDATA\Programs\Python\Python311\python.exe" -m venv "D:\venvs\cdd-s102"
& "D:\venvs\cdd-s102\Scripts\python.exe" -m pip install --upgrade pip
& "D:\venvs\cdd-s102\Scripts\python.exe" -m pip install -r requirements-sionna102.txt
```

验证：

```powershell
& "D:\venvs\cdd-s102\Scripts\python.exe" -m pytest -q
```

预期为 20 项测试通过。为了与本轮 CPU 基线保持一致，正式对比建议使用 CPU；CPU/GPU 不要求逐比特一致。

## 3. 400-trial 粗扫

```powershell
$py = "D:\venvs\cdd-s102\Scripts\python.exe"
& $py tools\run_plan024_sionna_tdl.py `
  --snrs 13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18 `
  --trials 400 --batch-size 20 --branches known `
  --run-id 20260722_main/prescan
```

程序在每个完整 SNR 点后写入 CSV。相同命令重复执行时会跳过已经完成且 trial 数相同的 QC/Sidon 点，可用于断点续跑。

当前工作区的 `prescan` 目录只完成了 13.0 和 13.5 dB，属于中断后的部分数据。若把整个 `outputs/experiment025_sionna_tdl_rmmse/20260722_main/prescan/` 一并复制到另一台电脑，以上命令会从 14.0 dB 继续；若不复制 outputs，则从头运行。

粗扫完成后检查 `prescan/sidon_qc_bler.csv`。只有原精化区间仍覆盖 10% 和 1% 交叉时，才执行下面两步。

## 4. 3000-trial 精化

10% BLER 区间：

```powershell
& $py tools\run_plan024_sionna_tdl.py `
  --snrs 14.25,14.5,14.75 `
  --trials 3000 --batch-size 20 --branches known `
  --run-id 20260722_main/refine_10pct
```

1% BLER 区间：

```powershell
& $py tools\run_plan024_sionna_tdl.py `
  --snrs 16,16.25,16.5,16.75,17,17.25,17.5 `
  --trials 3000 --batch-size 20 --branches known `
  --run-id 20260722_main/refine_1pct
```

若粗扫显示交叉已经移出上述区间，不运行不覆盖交叉点的精化命令，先根据粗扫重新确定 SNR 点。

## 5. 合并分析

```powershell
& $py tools\analyze_plan025_sionna_e2.py `
  --prescan outputs/experiment025_sionna_tdl_rmmse/20260722_main/prescan `
  --refine10 outputs/experiment025_sionna_tdl_rmmse/20260722_main/refine_10pct `
  --refine1 outputs/experiment025_sionna_tdl_rmmse/20260722_main/refine_1pct `
  --out outputs/experiment025_sionna_tdl_rmmse/20260722_main/final
```

主要输出：

- 每阶段 `resolved_experiment.json`、`environment.json`；
- 每阶段 `sidon_qc_bler.csv`、`paired_error_counts.csv`；
- `final/final_summary.json`：10%/1% 目标 SNR、保守 95% 区间、Sidon 相对 QC 改善及与 result-024 的差值；
- `final/sidon_qc_refined_combined.csv`：精化区原始数字；
- `final/sidon_qc_tdl_bler.png`：最终曲线。

1% 目标附近任一候选累计错误块少于 30 时，只报告先导结果，不做确定性 1% 结论。
