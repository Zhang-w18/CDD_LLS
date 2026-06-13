# CDD LLS Simulation Platform

这是基于 Sionna LDPC 的 QC 显式 CDD 链路级仿真平台第一版。平台复用 `lls_platform_slm` 的工程风格：YAML 配置驱动、`run.py` 统一入口、CSV/JSON/PNG 输出。

## 环境

优先复用本机已有 Sionna 环境：

```bash
/Users/zhangwei/Downloads/lls_platform_sc_mimo/.venv-sionna1/bin/python run.py --config configs/smoke.yaml
```

该环境已验证包含：

- TensorFlow 2.19.1
- Sionna 1.2.0
- NumPy / PyYAML / Matplotlib

## 快速运行

最小闭环 smoke：

```bash
/Users/zhangwei/Downloads/lls_platform_sc_mimo/.venv-sionna1/bin/python run.py --config configs/smoke.yaml
```

覆盖 PRG、ideal CSI、RMMSE、pairwise reconstruction、basis reconstruction 的小测试：

```bash
/Users/zhangwei/Downloads/lls_platform_sc_mimo/.venv-sionna1/bin/python run.py --config configs/smoke_variants.yaml
```

QC 静态 TDL 趋势复现配置：

```bash
/Users/zhangwei/Downloads/lls_platform_sc_mimo/.venv-sionna1/bin/python run.py --config configs/config_qc_static_tdl.yaml
```

RMMSE vs reconstruction 对比配置：

```bash
/Users/zhangwei/Downloads/lls_platform_sc_mimo/.venv-sionna1/bin/python run.py --config configs/config_recon_vs_rmmse.yaml
```

## 输出

每次运行会生成一个时间戳目录，例如：

```text
outputs/smoke/sim_YYYYMMDD_HHMMSS/
```

主要文件：

- `resolved_config.yaml`：本次运行展开后的完整配置。
- `summary.csv` / `summary.json`：每个 scenario/variant/SNR 的 BLER、NMSE、goodput 和配置元数据。
- `summary_10pct_bler_snr.csv`：按目标 BLER 插值得到的 SNR 和相对 PRG baseline 增益。
- `trial_metrics.csv`：仅当 `simulation.save_trial_metrics: true` 时输出逐 trial 记录。
- `bler_curves.png`：BLER vs SNR 曲线。
- `ce_nmse_curves.png`：等效信道 NMSE vs SNR 曲线。

## 测试

不依赖 pytest，直接用标准库 unittest：

```bash
/Users/zhangwei/Downloads/lls_platform_sc_mimo/.venv-sionna1/bin/python -m unittest discover -s tests
```

## 当前范围

第一版实现 static TDL、单层 PDSCH、2Tx/4Tx 到 4Rx、CDD/PRG/NO_CDD、IDEAL/LS/RMMSE/reconstruction 信道估计，以及 Sionna LDPC bit-level BLER。UE mobility、CDL、宽带 massive-MIMO 预编码暂未纳入第一版闭环。
