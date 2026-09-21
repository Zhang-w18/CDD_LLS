# PDCCH 仿真使用说明

平台现支持单层 PDCCH DCI BLER：标准 CRC24C（全 1 初始化）与 RNTI 掩码、Polar 编码/速率匹配、PDCCH 扰码、QPSK、CCE–REG 映射、DMRS、Sionna TDL-A～E、AWGN、MRC，以及 ideal CSI、REG-bundle 内 LS 线性估计或频域 LMMSE 信道估计。PDCCH 没有 PDSCH 式 MCS 索引，调制固定为 QPSK；`E=108×aggregation_level`。

YAML 可配置：payload/RNTI/scrambling ID、CORESET RB 与 1～3 symbols、首 CCE、CCE–REG mapping、REG bundle、DFT 轮询、天线、TDL 参数、SNR 和 trial 预算。预编码在同一 bundle 的 data/DMRS 上不变。端到端支持 AL1/2/4/8：AL4 由 Sionna 完成 puncturing，AL8 以 `N=512` 为母码，在 wrapper 中标准 repetition 至 `E=864`，接收端合并重复 LLR；AL16 尚未开放。不含候选散列和 CFO/ICI。

`pdcch-cdd-bler-v2` 源自 plan-031 的频域 CDD/相位方案比较，现支持通过 `antenna.n_rx` 配置正整数接收天线数。所有 CDD 和循环 DFT 预编码向量均归一化为每个 RE 总发射功率 1，噪声方差为 `1/SNR_linear`，且不随 Tx 或 Rx 数缩放。默认各 Rx 分支独立且无空间相关；UE 对各 Rx 分支独立执行相同的频域或时频 LMMSE，再在 Rx 维执行相干 MRC。CDD 候选按候选带宽 `K=72×AL` 构造，循环 DFT 候选按 6RB PRG 独立估计。候选可以通过 `--candidate` 独立运行；此时工具从基础种子和候选 ID 派生互不重叠的随机流，适用于并行执行和断点续跑。既有 plan-031 配置仍显式固定为 1Rx，其历史结果口径不变。

AL8 示例：`configs/pdcch_result030_al8_formal.yaml`。先校验，再运行：

```powershell
python tools/run_pdcch_bler_curves.py --config configs/pdcch_result030_al8_formal.yaml --stage validate
python tools/run_pdcch_bler_curves.py --config configs/pdcch_result030_al8_formal.yaml --stage run
```

输出目录由 `output_dir` 指定，包含展开配置、逐 SNR BLER/Wilson 95% 区间/CE NMSE CSV、逐 trial 错误标志与 CE 数组、日志、元数据和曲线 PNG。

plan-031 的独立候选运行示例：

```powershell
python tools/run_pdcch_bler_curves.py --config configs/pdcch_result031_al4_formal.yaml --stage validate
python tools/run_pdcch_bler_curves.py --config configs/pdcch_result031_al4_formal.yaml --stage run --candidate B0_QC
```
