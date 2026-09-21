# result-031-text：PDCCH CDD 时延 BLER（无图版）

> 对应 `research/plan-031-PDCCH-CDD时延-BLER.md`。实现、正式运行和结果审计已完成，结论待研究者确认。本文件不依赖图片，可由表格和精确输出独立核验。

## 1. 结论

事实：AL2/4/8 共完成 389 个正式 SNR 点、8,845,000 trials 和 206,190 个 DCI block errors，其中新增透明小时延 CDD 为 25 点、379,500 trials、14,962 errors。所有逐 trial error flags 的长度与错误和均和 CSV 一致，CE 汇总复算误差为 0。目标 SNR 使用 Jeffreys 平滑、按 trials 加权的单调递减 isotonic 曲线和局部 `log10(BLER)` 插值；95% 区间来自 2,000 次独立 Bernoulli bootstrap。两个没有双侧实测 bracket 的旧候选目标按 plan 不外推；新增透明小时延的六个目标均以 0.25 dB 原始网格双侧夹住，2,000 次 bootstrap 全部有效。

2026-09-07 补充的 TDL-C 300 ns、4Tx、2-symbol、AL1/2/4 场景另完成 333 个 SNR 点、5,628,400 trials 和 190,071 errors；56/58 个目标获得不宽于 0.25 dB 的原始 bracket。AL1/AL2 中 `AP_TEPS_T1` 与 `S0_SIDON` 最有竞争力，AL4 则是 `B0_QC` 最低；AL1 两条 GEO soft-control 出现高 SNR error floor，1%不外推。新增结果和二维 LMMSE 实现详见第 7 节。

2026-09-10 严格 Sidon 搜索补充先对 48 个候选预扫，按预注册 CE/BLER floor 联合判据停止 5 个明确失败候选，再由 3,000-trial 粗确认选出 18 个可能最优候选。最终细扫描与边界补点合计 172 个新 SNR 点、2,644,000 trials、92,766 errors，36/36 个 10%/1% 目标均由不宽于 0.25 dB 的原始双侧点夹住。1% 主判据下，A100 AL2/AL4/AL8 的点估计领先者分别为候选 02/01/01，C300 AL1/AL2/AL4 分别为候选 04/04/06；A100 AL4 的 01/06 与 C300 AL2 的 04/01/06 在 95% 水平下不能区分，其余四场景为唯一统计领先者。详见第 9 节。

2026-09-14 的 C300、4Tx/2Rx、2-symbol 正式补充目前已完成 AL1/AL2，合计 571 个 SNR 点、10,718,050 trials 和 377,050 errors；AL4 仍在运行，未纳入阶段性结论。AL1/AL2 的 72 个目标中 67 个具备不宽于 0.25 dB 的原始双侧 bracket。2Rx 相对对应 1Rx 的可比较目标全部降低所需 SNR，但该差值同时包含独立接收分支的阵列/分集收益和 BLER 尾部变化，不能解释为信道估计器本身获得同等 dB 改善。详见第 10 节。

推断：AL2 的 `B0_QC`、`AP_TEPS_T1`、`GEO_T1_CTRL` 在 10%/1% 处未显示可分辨差异。AL4 的 `S0_SIDON` 相对 `B0_QC` 在 10%/1% 改善 0.177 dB `[0.116,0.252]`/0.217 dB `[0.096,0.320]`。AL8 的 cycling 相对 `B0_QC` 改善 0.676 dB `[0.607,0.747]`/0.334 dB `[0.187,0.514]`；CDD 中 `AP_RMS_T1` 最好，改善 0.567 dB `[0.500,0.618]`/0.462 dB `[0.373,0.560]`。新增透明小时延 CDD 在 AL2/4/8 的 1% 目标分别为 9.899/3.688/-0.899 dB，相对透明 PRG 分别差 1.103/1.557/1.943 dB；带宽越大，忽略 CDD delay 的全带协方差失配越明显。这些结论只在同一 AL 内成立。

## 2. 配置、实现与统计口径

- 固定 A=41、QPSK、8 Tx/1 Rx、48-RB/1-symbol CORESET、30 kHz、FFT 4096、CP 288、non-interleaved mapping、first CCE 0、L=6、static Sionna TDL-A 100 ns、3.5 GHz、0 km/h。AL2/4/8 分别实占 12/24/48 RB，`K=144/288/576`、`E=216/432/864`、DMRS 为 36/72/144 RE。
- CDD 的 data/DMRS 共用单位范数频域 precoder。研究者于 2026-09-07 确认的最终接收机口径如下；“透明”表示 UE 不知道具体预编码向量或 CDD delay，也不使用候选特定的真实等效协方差。

| 方案 | 透明性 | UE 已知量 | 信道估计 |
|---|---|---|---|
| `A100_PRG_DFT8_6RB` | 透明 | 物理 PDP、6-RB PRG 边界；不知道各 PRG 的 DFT 向量索引 | 每个 PRG 独立使用 `R_phy` 做 frequency LMMSE，不跨 PRG |
| 小时延 CDD（修订 ID：`A100_SMALL_CDD_QSTEP0P25_TRANSPARENT_CDD`） | 透明 | 物理 PDP；不知道 CDD delay/`V` | 整个候选占用带宽使用 `R_phy` 做 frequency LMMSE |
| `A100_SMALL_CDD_QSTEP0P25_MATCHED_CDD`（原曲线） | 非透明 | 小时延 CDD 的真实 `V`、物理 PDP | 整个候选占用带宽使用真实 `R_g = R_phy ⊙ (V V^H)` 做 matched frequency LMMSE |
| `B0_QC`、`S0_SIDON`、`AP_RMS_T1`、`AP_TEPS_T1`、`GEO_T1_CTRL` | 非透明 | 各自的 `V`、物理 PDP | 整个候选占用带宽使用真实 `R_g = R_phy ⊙ (V V^H)` 做 matched frequency LMMSE |

`A100_SMALL_CDD_QSTEP0P25_MATCHED_CDD` 的原正式数据已原样保留，并在图例中明确标为 non-transparent；新增 `A100_SMALL_CDD_QSTEP0P25_TRANSPARENT_CDD` 使用相同的 data/DMRS 发射 precoder 和物理时延，只把估计器协方差改为全带 `R_phy`。全局 `channel_estimation: frequency_lmmse` 只指定算法类别；候选字段 `receiver_covariance_mode` 分别为 `matched_effective`、`physical_fullband` 和 cycling 的 `physical_prg`，不再依靠 candidate 名称推断口径。

这与 result-028 的 small-delay transparent 定义一致：028 在总发射功率 8 的归一化下使用 `8 R_phy`，本轮采用单位总发射功率，因而使用 `R_phy`；协方差和噪声同比缩放后 LMMSE 权重等价。

- cycling 在 AL2/4/8 分别轮询 DFT8 `[0,1]`、`[0,1,2,3]`、`[0,1,2,3,4,5,6,7]`。
- SNR 为单位总发射功率、单位能量 active PDCCH RE 的 Es/N0。每点至少 10,000 trials、目标 200 errors、上限 50,000。
- 根据研究者许可，各方案使用独立随机流；seed 由全局 seed `20260903`、candidate ID 和固定标签唯一派生。增益为“基线目标 SNR − candidate 目标 SNR”，正值表示 candidate 所需 SNR 更低；区间使用独立而非 paired bootstrap。
- 正式运行只生成候选 active subcarrier 上的 TDL 频响；同 seed 与 4096-bin 全 FFT 版本在 active bins 上通过 `rtol=atol=1e-12` 等价测试。首次 full-FFT/batch=200 尝试在产生任何新正式点前因内存压力终止，不计入本结果。

包含两条小时延曲线的当前配置 SHA-256 分别为 AL2 `c251c72f27182b7aa826b3ab612e11339d2148f2e8c19aa554d9592e3228f63f`、AL4 `1125f94cd5a59b9419f4488e9dab001bc9cca8624902ecd85b2cb828874de9b4`、AL8 `79d4b6b5d6112a3b09bbd08720713c432eb6569aec647fa670ef42b8e85ce3c1`。基准提交 `83b25bc2ea77324ae98f7dc67136fc1e93b1f32c`；运行环境为 Python 3.11.9、NumPy 1.26.4、SciPy 1.15.3、TensorFlow 2.15.1、Sionna 1.0.2，CPU-only。

## 3. 目标 SNR 与相对增益

单位均为 dB，格式为“点估计 `[95% CI]`”；`—` 表示未形成双侧实测 bracket，目标及相关增益均不估计。

### AL2

| Candidate | SNR@10% | gain vs B0 | gain vs PRG | SNR@1% | gain vs B0 | gain vs PRG |
|---|---:|---:|---:|---:|---:|---:|
| `B0_QC` | 2.831 [2.738, 2.909] | 0.000 [0.000, 0.000] | 1.226 [1.082, 1.374] | 5.018 [4.934, 5.067] | 0.000 [0.000, 0.000] | 3.778 [3.543, 3.965] |
| `S0_SIDON` | 3.036 [2.974, 3.098] | -0.205 [-0.322, -0.103] | 1.021 [0.891, 1.159] | 5.326 [5.228, 5.390] | -0.308 [-0.417, -0.198] | 3.470 [3.233, 3.664] |
| `AP_RMS_T1` | 3.079 [2.982, 3.160] | -0.247 [-0.376, -0.125] | 0.979 [0.835, 1.124] | 6.593 [6.429, 6.694] | -1.575 [-1.710, -1.404] | 2.203 [1.961, 2.428] |
| `AP_TEPS_T1` | 2.784 [2.722, 2.839] | 0.047 [-0.066, 0.154] | 1.273 [1.141, 1.406] | 4.991 [4.917, 5.048] | 0.026 [-0.079, 0.117] | 3.804 [3.563, 3.978] |
| `GEO_T1_CTRL` | 2.819 [2.756, 2.865] | 0.012 [-0.093, 0.113] | 1.238 [1.113, 1.359] | 4.885 [4.742, 5.011] | 0.133 [-0.013, 0.285] | 3.911 [3.655, 4.114] |
| `A100_PRG_DFT8_6RB` | 4.057 [3.948, 4.168] | -1.226 [-1.374, -1.082] | 0.000 [0.000, 0.000] | 8.796 [8.562, 8.952] | -3.778 [-3.965, -3.543] | 0.000 [0.000, 0.000] |
| `A100_SMALL_CDD_QSTEP0P25_MATCHED_CDD` | 5.028 [4.822, 5.079] | -2.196 [-2.303, -1.979] | -0.970 [-1.093, -0.744] | 9.480 [9.300, 9.760] | -4.462 [-4.748, -4.271] | -0.684 [-1.042, -0.440] |
| `A100_SMALL_CDD_QSTEP0P25_TRANSPARENT_CDD` | 4.880 [4.798, 4.980] | -2.049 [-2.180, -1.934] | -0.823 [-0.968, -0.680] | 9.899 [9.645, 10.040] | -4.881 [-5.051, -4.630] | -1.103 [-1.378, -0.780] |

### AL4

| Candidate | SNR@10% | gain vs B0 | gain vs PRG | SNR@1% | gain vs B0 | gain vs PRG |
|---|---:|---:|---:|---:|---:|---:|
| `B0_QC` | -1.062 [-1.101, -1.009] | 0.000 [0.000, 0.000] | 0.067 [-0.027, 0.142] | 0.881 [0.813, 0.935] | 0.000 [0.000, 0.000] | 1.250 [1.076, 1.424] |
| `S0_SIDON` | -1.240 [-1.292, -1.197] | 0.177 [0.116, 0.252] | 0.245 [0.152, 0.334] | 0.665 [0.573, 0.766] | 0.217 [0.096, 0.320] | 1.466 [1.275, 1.646] |
| `AP_RMS_T1` | — | — | — | 0.863 [0.736, 0.946] | 0.018 [-0.090, 0.154] | 1.268 [1.075, 1.463] |
| `AP_TEPS_T1` | -1.088 [-1.145, -1.013] | 0.026 [-0.055, 0.103] | 0.093 [-0.013, 0.178] | 0.779 [0.697, 0.849] | 0.102 [0.005, 0.198] | 1.352 [1.172, 1.533] |
| `GEO_T1_CTRL` | -0.983 [-1.036, -0.945] | -0.079 [-0.135, -0.007] | -0.012 [-0.098, 0.074] | 0.925 [0.801, 1.024] | -0.044 [-0.165, 0.095] | 1.206 [1.015, 1.401] |
| `A100_PRG_DFT8_6RB` | -0.995 [-1.077, -0.932] | -0.067 [-0.142, 0.027] | 0.000 [0.000, 0.000] | 2.131 [1.969, 2.285] | -1.250 [-1.424, -1.076] | 0.000 [0.000, 0.000] |
| `A100_SMALL_CDD_QSTEP0P25_MATCHED_CDD` | -0.077 [-0.161, 0.026] | -0.986 [-1.096, -0.888] | -0.918 [-1.054, -0.811] | 3.559 [3.416, 3.637] | -2.678 [-2.782, -2.525] | -1.428 [-1.613, -1.211] |
| `A100_SMALL_CDD_QSTEP0P25_TRANSPARENT_CDD` | 0.014 [-0.057, 0.084] | -1.076 [-1.156, -0.988] | -1.009 [-1.120, -0.910] | 3.688 [3.615, 3.787] | -2.807 [-2.933, -2.717] | -1.557 [-1.757, -1.382] |

### AL8

| Candidate | SNR@10% | gain vs B0 | gain vs PRG | SNR@1% | gain vs B0 | gain vs PRG |
|---|---:|---:|---:|---:|---:|---:|
| `B0_QC` | -4.255 [-4.309, -4.217] | 0.000 [0.000, 0.000] | -0.676 [-0.747, -0.607] | -2.508 [-2.577, -2.449] | 0.000 [0.000, 0.000] | -0.334 [-0.514, -0.187] |
| `S0_SIDON` | -4.552 [-4.596, -4.495] | 0.297 [0.218, 0.357] | -0.379 [-0.460, -0.316] | -2.792 [-2.866, -2.721] | 0.284 [0.178, 0.380] | -0.049 [-0.243, 0.104] |
| `AP_RMS_T1` | -4.822 [-4.854, -4.780] | 0.567 [0.500, 0.618] | -0.109 [-0.179, -0.054] | -2.970 [-3.046, -2.915] | 0.462 [0.373, 0.560] | 0.128 [-0.057, 0.282] |
| `AP_TEPS_T1` | -4.330 [-4.364, -4.287] | 0.075 [0.008, 0.127] | -0.601 [-0.677, -0.546] | -2.596 [-2.663, -2.527] | 0.087 [-0.012, 0.175] | -0.246 [-0.427, -0.107] |
| `GEO_T1_CTRL` | -4.165 [-4.213, -4.129] | -0.090 [-0.156, -0.028] | -0.766 [-0.832, -0.704] | -2.376 [-2.448, -2.307] | -0.132 [-0.228, -0.044] | -0.466 [-0.648, -0.323] |
| `A100_PRG_DFT8_6RB` | -4.931 [-4.992, -4.888] | 0.676 [0.607, 0.747] | 0.000 [0.000, 0.000] | -2.842 [-3.013, -2.725] | 0.334 [0.187, 0.514] | 0.000 [0.000, 0.000] |
| `A100_SMALL_CDD_QSTEP0P25_MATCHED_CDD` | -4.283 [-4.334, -4.214] | 0.028 [-0.060, 0.095] | -0.648 [-0.736, -0.578] | — | — | — |
| `A100_SMALL_CDD_QSTEP0P25_TRANSPARENT_CDD` | -3.423 [-3.480, -3.380] | -0.832 [-0.904, -0.762] | -1.508 [-1.580, -1.437] | -0.899 [-1.033, -0.781] | -1.610 [-1.744, -1.465] | -1.943 [-2.138, -1.768] |

新增透明小时延候选的原始 bracket 如下；每格格式为 `SNR dB: errors/trials=BLER`，直接来自各 AL candidate 目录的 `bler_points.csv`：

| AL | 10% 左/右原始点 | 1% 左/右原始点 |
|---:|---|---|
| 2 | `4.75: 1076/10000=0.1076`; `5.00: 934/10000=0.0934` | `9.75: 201/18300=0.010984`; `10.00: 200/21400=0.009346` |
| 4 | `0.00: 1008/10000=0.1008`; `0.25: 864/10000=0.0864` | `3.50: 202/15900=0.012704`; `3.75: 200/21700=0.009217` |
| 8 | `-3.50: 1072/10000=0.1072`; `-3.25: 854/10000=0.0854` | `-1.00: 200/18200=0.010989`; `-0.75: 200/23100=0.008658` |

## 4. 可辨识性与 CE floor

`-300 dB` 是数值零的报告下限。cycling 在 AL2/4 的全候选 rank 2/4 是 DFT bundle 轮询导致，其接收机按 bundle 单独估计，因此全候选 condition=`inf` 不表示单 bundle 不可估。

| AL | Candidate | rank | condition | zero-noise CE floor / dB |
|---:|---|---:|---:|---:|
| 2 | `B0_QC` | 8 | 1.223 | -67.38 |
| 2 | `S0_SIDON` | 8 | 1.000 | -14.96 |
| 2 | `AP_RMS_T1` | 8 | 861.1 | -300.00 |
| 2 | `AP_TEPS_T1` | 8 | 1.131 | -76.04 |
| 2 | `GEO_T1_CTRL` | 8 | 1.000 | -25.45 |
| 2 | `A100_PRG_DFT8_6RB` | 2 | inf | -300.00 |
| 2 | `A100_SMALL_CDD_QSTEP0P25_MATCHED_CDD` | 8 | 1.554e9 | -300.00 |
| 2 | `A100_SMALL_CDD_QSTEP0P25_TRANSPARENT_CDD` | 8 | 1.554e9 | -56.72 |
| 4 | `B0_QC` | 8 | 1.118 | -65.82 |
| 4 | `S0_SIDON` | 8 | 1.000 | -39.58 |
| 4 | `AP_RMS_T1` | 8 | 2.379 | -79.32 |
| 4 | `AP_TEPS_T1` | 8 | 1.104 | -79.01 |
| 4 | `GEO_T1_CTRL` | 8 | 1.000 | -15.50 |
| 4 | `A100_PRG_DFT8_6RB` | 4 | inf | -300.00 |
| 4 | `A100_SMALL_CDD_QSTEP0P25_MATCHED_CDD` | 8 | 1.048e7 | -300.00 |
| 4 | `A100_SMALL_CDD_QSTEP0P25_TRANSPARENT_CDD` | 8 | 1.048e7 | -40.15 |
| 8 | `B0_QC` | 8 | 1.000 | -300.00 |
| 8 | `S0_SIDON` | 8 | 1.000 | -300.00 |
| 8 | `AP_RMS_T1` | 8 | 1.407 | -85.00 |
| 8 | `AP_TEPS_T1` | 8 | 1.060 | -300.00 |
| 8 | `GEO_T1_CTRL` | 8 | 1.000 | -15.59 |
| 8 | `A100_PRG_DFT8_6RB` | 8 | 1.000 | -300.00 |
| 8 | `A100_SMALL_CDD_QSTEP0P25_MATCHED_CDD` | 8 | 5.970e4 | -69.03 |
| 8 | `A100_SMALL_CDD_QSTEP0P25_TRANSPARENT_CDD` | 8 | 5.970e4 | -34.90 |

## 5. 两组排序变化的机理分析

### 5.1 AL8 PDCCH 中 `AP_RMS_T1` 相对 cycling 的优势为何缩小

先校正数字：result-028 的 `AP_RMS_T1` 1% estimated-CSI BLER 目标为 16.039 dB；6-RB cycling 的实测 17/17.25 dB 点为 `0.010984/0.009038`，按本轮相同的局部 `log10(BLER)` 插值复算 crossing 为 17.120 dB，故可复现差值约 **1.081 dB**，不是严格 2 dB。result-031 AL8 的对应差值为 0.128 dB，独立 bootstrap 95%区间 `[-0.057,0.282]`，在 1% 处尚不可分辨。

| 项目 | result-028 PDSCH | result-031 PDCCH AL8 |
|---|---|---|
| 时频资源 | 10 symbols，5568 data RE | 1 symbol，432 data RE |
| DMRS | 2 symbols、comb-6，共 192 RE | 单 symbol、每 RB 3 RE，共 144 RE |
| 调制/编码 | 16QAM、LDPC、码率 553/1024 | QPSK、Polar，`K=65,N=512,E=864`，`K/E=0.0752` 且 repetition |
| CDD / cycling CE | 全带 matched / 每 6-RB PRG physical-covariance LMMSE | 全带 matched / 每 6-RB PRG physical-covariance LMMSE |

两轮的功率归一化等价：result-028 的向量范数平方 8 配合 `noise=8/SNR`，result-031 的单位范数配合 `noise=1/SNR`。同 SNR 事实也排除了“平均 CE NMSE 单独解释差值”：result-028 @16 dB 的 AP/PRG CE NMSE 为 `-24.497/-24.265 dB`、BLER 为 `0.010412/0.034661`；result-031 @-3 dB 的 CE NMSE 为 `-7.012/-7.718 dB`，PRG 估计反而好 0.706 dB，而 BLER 仍为接近的 `0.010526/0.011167`。

机理推断：单位范数下，两方案主要改变跨 RE 的联合相关性，而不是单 RE 边缘功率。译码错误可近似理解为 `Pr[mean_i I_M(γ_i)<R_eff]`；候选改变相关的 `γ_i`，调制改变 `I_M` 的非线性，码率/码长/交织和 rate matching 改变门限及各 RE 权重。result-028 的较高码率 16QAM LDPC 长码字跨 5568 RE，static 信道又让相同频率状态在多个 symbol 重复，cycling 的 6-RB 深衰落会成簇影响大量 coded bits；AP 的全带相关结构因此可在低尾部体现约 1.1 dB 优势。result-031 的极低码率 QPSK Polar 能容忍更多弱 RE，且 `E>N` 会重复部分 mother-code bits，16QAM 位信道不均衡和较高译码门限所放大的相关性差异可能被压缩；短 Polar 的 interleaving/error events 也不同。

结论：编码调制很可能参与优势压缩，但两轮还同时改变了 DMRS、码长、data RE 数、时域复用和估计观测，不能把全部变化归因于编码调制。定量拆分需要同批信道上的 ideal-CSI QPSK/16QAM BICM outage，以及固定 DMRS/资源、逐项替换调制、码率和编码器的 factorial 实验。

result-028 精确证据来自 `outputs/experiment028_csi_curves/20260803_main/section37_receiver_comparison/a100/final/a100_section37_<estimated_csi_bler|ce_nmse>_points.csv` 和 `outputs/experiment028_csi_curves/20260803_main/transparent_prg_6rb_1pct_extension/a100/final/bler_bracket_audit.csv`；result-031 使用本轮 `analysis/formal_points.csv`、`target_snr.csv`、`target_gains.csv`。

### 5.2 AL2 的 GEO 与 AL4 的 S0 为什么互换排序

两者都有 36 个互异的模 K 无序二元和，均满足严格 Sidon；但 GEO 不是硬厚 Sidon。S0 是跨 K 固定的 `[0,1,3,7,12,20,30,65]`，只保证精确唯一；GEO 按每个 K 重新生成，以 `T_0.01` 对 pair-sum/fold gap 做软 Pareto 平衡。区别不只是是否加余量，还包括整套坐标、物理时延及 PDP 加权相关结构。

| AL | candidate | pair/fold gap | delay span | CE floor | 共同 SNR 的 CE NMSE | SNR@1% |
|---:|---|---:|---:|---:|---:|---:|
| 2 | `S0_SIDON` | 1q / 1q | 65q = 15.046 µs | -14.96 dB | -7.612 dB @ 5 dB | 5.326 dB |
| 2 | `GEO_T1_CTRL` | 1q / 3q | 133q = 30.787 µs | -25.45 dB | -8.067 dB @ 5 dB | 4.885 dB |
| 4 | `S0_SIDON` | 1q / 1q | 65q = 7.523 µs | -39.58 dB | -6.909 dB @ 0.75 dB | 0.665 dB |
| 4 | `GEO_T1_CTRL` | 3q / 3q | 271q = 31.366 µs | -15.50 dB | -6.188 dB @ 0.75 dB | 0.925 dB |

事实：GEO 在 AL4 的最小 pair/fold gap 仍更大却更差，说明最小余量不是单调 BLER 指标。四种组合的 pilot rank/condition 都是 `8/≈1`，但 CE floor、有限 SNR CE NMSE 和 1% BLER 的排序同步反转；问题不在 pilot 列可辨识性，而在实际 `R_PP/R_DP` 的数据可预测性。

机理推断：严格 Sidon 只消除纯平坦 DFT 模型的精确四阶共振，不控制 PDP 加权近碰撞、其他阶相关矩或 matched LMMSE 的 `R_DP`。K 从 144 变为 288 后，`q_K` 减半、DMRS 模环从 36 变为 72；固定 S0 的物理 span 减半，而 GEO 被换成全新的坐标并保持约 31 µs span。人工相位变化、fold 中心与物理 PDP 的重叠和 data-to-pilot 预测关系因此都变了。当前最强解释是 AL2 GEO 的 layout 更利于 CE，而 AL4 GEO 的新 layout 使 `R_DP` 和 CE floor 变差，S0 因而反超。要形成因果结论，还需 ideal-CSI PDCCH、实际 PDP 下的 `J_CDD` 分阶/conditional-NMSE/MI 尾部，以及“固定物理 delay”对“固定 j”的跨 AL 消融。

### 5.3 小时延 CDD 的 transparent 与 non-transparent 曲线

这两条曲线的真实 CDD delay、逐 RE `V`、物理信道分布、data/DMRS 映射和功率归一化完全相同；唯一有意的配置差异是估计器使用 `R_g`（non-transparent）还是只使用 `R_phy`（transparent）。两者使用独立随机流而非相同 realization，因此点估计差值还包含有限样本波动；其系统性差异对应接收机协方差知识，而不是发射波形变化。

| AL | transparent SNR@10% | non-transparent SNR@10% | 点估计差值（transparent − non-transparent） | transparent SNR@1% | non-transparent SNR@1% | 点估计差值 |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 4.880 | 5.028 | -0.148 | 9.899 | 9.480 | +0.419 |
| 4 | 0.014 | -0.077 | +0.091 | 3.688 | 3.559 | +0.129 |
| 8 | -3.423 | -4.283 | +0.860 | -0.899 | — | — |

事实：AL2/4 的差值较小且两条目标 SNR 的各自 95%区间有重叠，不应仅凭点估计宣称稳定优劣；AL8 在 10% 已出现明显右移。旧 non-transparent AL8 1% 没有下侧实测点，不能外推或定量比较该目标。transparent 的 zero-noise CE floor 随 AL 从 -56.72、-40.15 升到 -34.90 dB，而 matched 为 -300、-300、-69.03 dB（报告下限 -300 dB）。

机理推断：物理小时延固定为 0–101.273 ns，AL 增大时占用带宽从 12/24/48 RB 扩展，CDD 相位在整个候选带内的累计旋转变大，真实 `R_g` 与假设 `R_phy` 的差异也变大。全带 transparent LMMSE 把 CDD 引入的额外频率结构误当成纯物理信道相关性，因此 data-from-pilot 的预测失配随带宽加重。这也解释了为何 transparent 小时延在 AL8 比逐 6-RB transparent PRG 差 1.943 dB：PRG 内 precoder 向量恒定，接收机只在 6-RB 内使用 `R_phy`，不会跨向量切换边界外推；小时延方案则在全 48 RB 上持续相位旋转，却要求接收机完全不知道该旋转。

## 6. 证据、验收与限制

| AL | points | trials | errors | adjacent BLER increases | 50k 且 errors<200 的点 |
|---:|---:|---:|---:|---:|---:|
| 2 | 164 | 3,920,900 | 82,716 | 5 | 44 |
| 4 | 126 | 3,047,700 | 56,649 | 2 | 36 |
| 8 | 99 | 1,876,400 | 66,825 | 0 | 10 |

事实：原始相邻反向均为有限样本尾部波动；原始点未被改写，目标估计使用规定的单调 isotonic 曲线。所有点满足最少 trials；未达到 200 errors 的点均达到 50,000 上限。

未完全满足的计划条件：AL4 `AP_RMS_T1` 在最低 -1.5 dB 已为 `977/10000=0.0977`，10% 未获上侧点；AL8 **non-transparent** 小时延 CDD 在最高 -1.75 dB 仍为 `201/17000=0.011824`，1% 未获下侧点。另有八个旧候选的可估计目标实测 bracket 超过 0.25 dB：AL2 `S0_SIDON`/`AP_RMS_T1` 的 10% 为 0.75 dB，AL2 cycling/non-transparent 小时延 CDD 的 1% 为 1.25 dB；AL4 `B0_QC`、`AP_RMS_T1`、`AP_TEPS_T1`、`GEO_T1_CTRL` 的 1% 为 0.5 dB。上述旧目标的网格分辨率未达到 plan 要求，bootstrap 区间不能替代该缺口；新增 transparent 小时延的六个目标均满足 0.25 dB。

精确证据位于 `outputs/experiment031_pdcch_cdd/20260903_main/`。各 candidate 目录含展开配置、`bler_points.csv`、`run_metadata.json`、`run.log` 和逐 trial flags；分析目录含 `formal_points.csv`、`diagnostics.csv`、`target_snr.csv`、`target_gains.csv`、`analysis_metadata.json`。这组表格足以在不读取图片的条件下复核全部结论。

复现：

```powershell
python tools/run_pdcch_bler_curves.py --config configs/pdcch_result031_al2_formal.yaml --stage validate
python tools/run_plan031_candidates.py --config configs/pdcch_result031_al2_formal.yaml --config configs/pdcch_result031_al4_formal.yaml --config configs/pdcch_result031_al8_formal.yaml --stage run --max-workers 4
python tools/analyze_plan031.py --bootstrap-repeats 2000
python -m pytest tests/test_pdcch.py tests/test_plan031_analysis.py -q
```

最终测试为 `19 passed`；警告仅来自 Matplotlib/PyParsing deprecation 和只读 font cache。当前模型不含候选散列、CFO/ICI、同步误差、跨时隙信道变化；独立随机流不支持 paired-sample 方差缩减。结果尚未确认，未更新 `KNOWLEDGE.md`/`GOALS.md`，未创建 Git checkpoint。

## 7. 新增场景：TDL-C 300 ns、4Tx、2-symbol CORESET

本节对应 `plan-031` 的 2026-09-07 补充计划，是独立于前述 TDL-A/8Tx/1-symbol 结果的新场景。两组场景同时改变了 TDL profile/DS、天线数、CORESET duration、AL、载频和移动性，因此只在本节各 AL 内比较 candidate，不把新旧场景差值归因于单一因素。

### 7.1 配置与接收机口径

固定 A=41、CRC24C/RNTI `0xFFFF`、QPSK、4 Tx/1 Rx、48-RB CORESET、2 symbols、30 kHz、FFT 4096、CP 288、non-interleaved、first CCE 0、REG bundle L=6、Sionna TDL-C 300 ns、3 km/h、4 GHz、20 sinusoids。速度为 0.833333 m/s，最大 Doppler 为 11.1188 Hz；实际信道 realization 在两个 symbol 间变化，不做静态复制或 DMRS 平均。

为避免与 Polar 编码中的信息长度混淆，本节把占用的唯一频率子载波数记作 $K_f$。每个 CCE 的 6 REG 在 2-symbol 映射中覆盖 `3 RB × 2 symbols`；两个 symbol 各有 comb-4 DMRS。

| AL | 占用 RB | $K_f$ | 每 symbol 唯一 DMRS SC | 总 DMRS RE | data RE | Polar $K_{info+CRC}$ | $E$ |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 3 | 36 | 9 | 18 | 54 | 65 | 108 |
| 2 | 6 | 72 | 18 | 36 | 108 | 65 | 216 |
| 4 | 12 | 144 | 36 | 72 | 216 | 65 | 432 |

这里 $E$ 是 rate matching 后映射到 QPSK 的编码比特数；`A=41` 加 24-bit CRC 后 Polar 输入为 $K_{info+CRC}=65$，并非表中的频率维数 $K_f$。

| candidate | 透明性 | 二维 LMMSE 使用的协方差及范围 |
|---|---|---|
| `C300_PRG_DFT4_6REG` | 透明 | UE 不知道各 bundle 的 DFT4 索引；每个 `3 RB × 2 symbols` bundle 独立使用物理 TDL-C 时频协方差，不跨 bundle |
| `C300_SMALL_CDD_QSTEP0P25_TRANSPARENT_CDD` | 透明 | UE 不知道 CDD delay/`V`；全候选带宽只使用物理 TDL-C 时频协方差 |
| `C300_SMALL_CDD_QSTEP0P25_MATCHED_CDD` | 非透明 | 发射波形与上一行完全相同；全带使用由真实 `V` 构造的等效时频协方差 |
| `C300_B0_QC`、`C300_S0_SIDON`、`C300_AP_RMS_T1`、`C300_AP_TEPS_T1`、全部 `C300_GEO_T1_CTRL_*` | 非透明 | 全带使用各自真实 `V` 构造的等效时频协方差 |

cycling 的 DFT4 顺序为 AL1 `[0]`、AL2 `[0,1]`、AL4 `[0,1,2,3]`。AL1 只有一个 bundle，因而该候选在 AL1 实际没有跨 bundle cycling；AL2/4 分别使用 2/4 个正交空间向量。

CDD 的物理 delay（ns）为：B0 `[0,520.833,1041.667,1562.500]`，AP-RMS `[0,300,600,900]`，AP-$T_{0.01}$ `[0,1891.950,3783.900,5675.850]`，small CDD `[0,14.468,28.935,43.403]`。`S0_SIDON` 在各 AL 均使用整数坐标 `[0,1,3,7]`；GEO 的 AL1 balance/pair、fold 为 `[0,2,6,14]`、`[0,2,5,16]`，AL2 balance/pair/fold 为 `[0,4,12,33]`、`[0,5,16,49]`、`[0,4,10,32]`，AL4 为 `[0,10,30,97]`、`[0,11,44,66]`、`[0,9,18,27]`。坐标单位 $q=1/(K_f\Delta f)$。穷举回执按含自和的 10 个 $i\le j$ pair-sum 计算，三组结果与 plan 冻结值完全一致。TDL-C 的 24-tap PDP 给出 $T_{0.01}=1891.950$ ns，区间内功率 0.99334723。

seed 为 `20260907`，candidate 间使用独立派生随机流。先以 300 trials/点、1 dB 预扫，正式点至少 10,000 trials、目标 200 errors、上限 50,000；目标采用 Jeffreys 平滑、按 trials 加权的递减 isotonic 拟合和局部 `log10(BLER)` 插值，95%区间用 2,000 次独立逐点 Bernoulli bootstrap。正式网格偏离预扫后，另追加 31 个可审计 bracket-completion 点；该后验补点不覆盖旧数据，完整选择见 `bracket_completion_receipt.json`。

### 7.2 用 `(symbol, subcarrier)` 坐标构造二维 LMMSE

实现没有先把两个 DMRS symbol 平均，也没有把二维网格误当作一条仅频率序列。每个 RE 由真实坐标 $c=(\ell,k)$ 标识，其中 $\ell\in\{0,1\}$ 是 CORESET 内 symbol index，$k$ 是实际 active-subcarrier index。按资源映射得到有序集合

$$
\mathcal P=\{p_i=(\ell_i,k_i)\}_{i=1}^{N_P},\qquad
\mathcal D=\{d_u=(\ell_u,k_u)\}_{u=1}^{N_D},
$$

其中 $(N_P,N_D)=(18,54),(36,108),(72,216)$ 分别对应 AL1/2/4。同一频率在两个 symbol 上出现时是两个不同 pilot 坐标和两次不同观测。

底层 TDL-C 采用可分离时频协方差。物理频率协方差由 24-tap PDP 直接在 active subcarrier 上计算，

$$
R_{f,phy}[k,k']
=\sum_a P_a\exp\{-j2\pi(f_k-f_{k'})\tau_a\},
$$

时间协方差 $R_t[\ell,\ell']$ 使用 Sionna TDL time-correlation，输入 3 km/h、4 GHz 和实际 OFDM symbol duration

$$
T_{sym}=\frac{1+288/4096}{30\,000}=35.677083\ \mu\mathrm{s}.
$$

因此透明接收机假设的任意两坐标间协方差为

$$
R_{phy}\big[(\ell,k),(\ell',k')\big]
=R_t[\ell,\ell']R_{f,phy}[k,k'].
$$

对已知 CDD/precoder 的 non-transparent 接收机，先令 $V$ 的第 $k$ 行为空间向量 $v_k$，再构造

$$
R_{f,g}=R_{f,phy}\odot(VV^H),\qquad
R_g\big[(\ell,k),(\ell',k')\big]
=R_t[\ell,\ell']R_{f,g}[k,k'].
$$

随后直接按二维坐标抽取 pilot-to-pilot 和 data-to-pilot 子矩阵：

$$
[R_{PP}]_{i,j}=R[p_i,p_j],\qquad
[R_{DP}]_{u,i}=R[d_u,p_i].
$$

DMRS 为单位模。接收值除以已知 DMRS 后得到

$$
z_i=\frac{y(p_i)}{x_{DMRS}(p_i)}=g(p_i)+n_i,
\qquad E[nn^H]=\sigma^2I.
$$

最终权重和 data-RE 信道估计为

$$
W=R_{DP}(R_{PP}+\sigma^2I+\epsilon I)^{-1},\qquad
\widehat{\boldsymbol g}_{D}=W\boldsymbol z_P,
$$

其中 $W$ 的形状是 $N_D\times N_P$，$\epsilon$ 仅为 Cholesky 所需的数值 diagonal loading。代码用 Cholesky 线性求解实现右侧逆，不显式求矩阵逆。这样 $W_{u,i}$ 同时利用 pilot $i$ 与 data $u$ 的 symbol 间时间相关和 subcarrier 间频率相关，正是 `(symbol, subcarrier)` 二维 LMMSE，而非逐 symbol frequency-LMMSE。

两种透明模式只改变 $R$ 的选择与作用范围：transparent small CDD 在全带使用 $R_{phy}$；cycling 为每个 36-subcarrier bundle 单独构造上述 $R_{PP}/R_{DP}$，把跨 bundle 的 $W$ 元素严格置零。non-transparent CDD 全带使用 $R_g$。报告的 zero-noise CE floor 则令 $\sigma^2=0$，用接收机所假设的 $W$ 作用于真实等效信道协方差；因此 transparent floor 会包含“不知道 CDD/DFT 向量”造成的模型失配。

### 7.3 BLER 目标与相对增益

单位均为 dB，格式为“点估计 `[95% CI]`”。gain 定义为“baseline SNR − candidate SNR”，正值表示 candidate 更好；`—` 表示未形成合格的 0.25 dB 双侧实测 bracket，不外推。

#### AL1

| candidate | SNR@10% | gain vs B0 | gain vs PRG | SNR@1% | gain vs B0 | gain vs PRG |
|---|---:|---:|---:|---:|---:|---:|
| `C300_B0_QC` | 8.974 [8.861,9.064] | +0.000 [+0.000,+0.000] | +1.872 [+1.662,+2.100] | 13.655 [13.525,13.794] | +0.000 [+0.000,+0.000] | +3.269 [+3.109,+3.428] |
| `C300_S0_SIDON` | 8.504 [8.418,8.619] | +0.470 [+0.292,+0.604] | +2.342 [+2.122,+2.560] | 12.789 [12.659,12.924] | +0.866 [+0.677,+1.064] | +4.135 [+3.973,+4.308] |
| `C300_AP_RMS_T1` | 9.551 [9.394,9.751] | -0.577 [-0.787,-0.388] | +1.295 [+1.048,+1.548] | 14.893 [14.699,15.054] | -1.237 [-1.432,-0.993] | +2.032 [+1.844,+2.261] |
| `C300_AP_TEPS_T1` | 8.441 [8.381,8.526] | +0.532 [+0.391,+0.648] | +2.404 [+2.211,+2.606] | 12.124 [11.953,12.279] | +1.531 [+1.346,+1.747] | +4.800 [+4.625,+4.999] |
| `C300_GEO_T1_CTRL_BAL_PAIR` | 9.536 [9.434,9.603] | -0.562 [-0.698,-0.424] | +1.310 [+1.116,+1.525] | — | — | — |
| `C300_GEO_T1_CTRL_FOLD` | 9.139 [9.044,9.268] | -0.166 [-0.340,-0.036] | +1.707 [+1.489,+1.928] | — | — | — |
| `C300_PRG_DFT4_6REG` | 10.846 [10.675,11.042] | -1.872 [-2.100,-1.662] | +0.000 [+0.000,+0.000] | 16.924 [16.842,17.036] | -3.269 [-3.428,-3.109] | +0.000 [+0.000,+0.000] |
| `C300_SMALL_CDD_QSTEP0P25_MATCHED_CDD` | 10.819 [10.689,10.946] | -1.846 [-2.024,-1.677] | +0.026 [-0.191,+0.250] | 16.859 [16.704,17.040] | -3.203 [-3.419,-3.001] | +0.065 [-0.138,+0.258] |
| `C300_SMALL_CDD_QSTEP0P25_TRANSPARENT_CDD` | 10.870 [10.769,10.956] | -1.896 [-2.035,-1.761] | -0.024 [-0.216,+0.193] | 16.992 [16.895,17.293] | -3.336 [-3.661,-3.174] | -0.068 [-0.385,+0.077] |

##### AL1 仿真结果分析

AL1 在 $10^{-2}$ 附近出现的局部非单调首先是有限样本现象。本轮各 candidate、各 SNR 点独立生成信道、噪声与码字样本，并在累计约 200 个误块后停止；若某点恰有 $e\approx200$ 个误块，则 BLER 的相对标准误差量级约为 $1/\sqrt{e}=7.1\%$，95% 相对波动量级约为 $\pm14\%$。例如，`S0_SIDON` 在 13.50/13.75 dB 的原始 BLER 分别为 200/32400=0.00617 和 200/32000=0.00625，对应 Wilson 区间 `[0.00538,0.00709]` 与 `[0.00544,0.00717]`；`AP_RMS_T1` 在 14.50/14.75 dB 为 200/18700=0.01070 和 200/17400=0.01149，区间 `[0.00932,0.01227]` 与 `[0.01001,0.01319]`；small CDD matched 在 16.25/16.50 dB 为 202/16100=0.01255 和 200/14600=0.01370，区间 `[0.01094,0.01439]` 与 `[0.01194,0.01572]`。这些相邻点的区间高度重叠，因而不能据此判定提高 SNR 后真实性能变差；目标 SNR 提取中的 isotonic 处理只用于抑制这种抽样起伏，原始数据未被修改。

GEO 曲线还叠加了另一种、具有物理意义的非单调来源。`GEO_BAL_PAIR` 与 `GEO_FOLD` 的 zero-noise CE floor 分别仅为 -12.69 dB 和 -13.67 dB，高 SNR 区 CE NMSE 也分别停留在约 -11.25～-11.49 dB 和 -12.09～-12.30 dB；对应 18～20 dB 的 BLER 分别约为 1.60%～2.30% 和 0.90%～1.48%。尽管二者 pilot precoder 均满秩且 condition 约为 1，完整的 data-from-pilot 协方差 $R_{DP}$ 仍存在不可消除的预测误差。机理上，当前 coherent equalizer/demapper 把估计信道当作精确信道，LLR 噪声方差计入 AWGN、但未显式计入残余信道估计误差；当 AWGN 随 SNR 降低时，残余误差可能形成过置信的错误 LLR，从而产生 error floor 及轻微反弹。这是与现有 CE/BLER 证据一致的机理解释，而非已证明的唯一原因。相比之下，`AP_RMS_T1` 和 small CDD matched 的 CE floor 分别低至 -70.77 dB 和 -93.00 dB，其上述小幅反弹应主要归因于 Monte Carlo 波动；`S0_SIDON` 的 -16.90 dB floor 可能影响更高 SNR 的渐近行为，但当前两个相邻点仍无法在统计上区分。

#### AL2

| candidate | SNR@10% | gain vs B0 | gain vs PRG | SNR@1% | gain vs B0 | gain vs PRG |
|---|---:|---:|---:|---:|---:|---:|
| `C300_B0_QC` | 2.905 [2.841,2.994] | +0.000 [+0.000,+0.000] | +1.255 [+1.139,+1.573] | 6.029 [5.849,6.208] | +0.000 [+0.000,+0.000] | +2.695 [+2.290,+2.898] |
| `C300_S0_SIDON` | 2.744 [2.673,2.813] | +0.161 [+0.064,+0.278] | +1.416 [+1.308,+1.727] | 5.585 [5.474,5.673] | +0.444 [+0.251,+0.653] | +3.140 [+2.803,+3.281] |
| `C300_AP_RMS_T1` | 3.445 [3.345,3.545] | -0.540 [-0.667,-0.403] | +0.714 [+0.573,+1.007] | 7.020 [6.844,7.281] | -0.990 [-1.322,-0.738] | +1.705 [+1.229,+1.904] |
| `C300_AP_TEPS_T1` | 2.881 [2.813,2.956] | +0.024 [-0.079,+0.130] | +1.279 [+1.159,+1.595] | 5.516 [5.424,5.587] | +0.513 [+0.314,+0.716] | +3.208 [+2.877,+3.341] |
| `C300_GEO_T1_CTRL_BAL` | 3.195 [3.147,3.262] | -0.290 [-0.383,-0.191] | +0.965 [+0.852,+1.294] | 6.254 [6.135,6.389] | -0.224 [-0.452,-0.000] | +2.471 [+2.097,+2.628] |
| `C300_GEO_T1_CTRL_PAIR` | 3.253 [3.174,3.332] | -0.348 [-0.451,-0.233] | +0.907 [+0.791,+1.224] | 6.511 [6.164,6.617] | -0.481 [-0.687,-0.059] | +2.214 [+1.862,+2.586] |
| `C300_GEO_T1_CTRL_FOLD` | 2.988 [2.936,3.062] | -0.083 [-0.184,+0.016] | +1.172 [+1.051,+1.482] | 6.028 [5.942,6.090] | +0.002 [-0.193,+0.197] | +2.697 [+2.370,+2.819] |
| `C300_PRG_DFT4_6REG` | 4.160 [4.068,4.502] | -1.255 [-1.573,-1.139] | +0.000 [+0.000,+0.000] | 8.725 [8.414,8.807] | -2.695 [-2.898,-2.290] | +0.000 [+0.000,+0.000] |
| `C300_SMALL_CDD_QSTEP0P25_MATCHED_CDD` | 5.166 [5.091,5.268] | -2.261 [-2.384,-2.153] | -1.006 [-1.134,-0.702] | 9.869 [9.689,10.050] | -3.840 [-4.095,-3.574] | -1.144 [-1.541,-0.939] |
| `C300_SMALL_CDD_QSTEP0P25_TRANSPARENT_CDD` | 5.137 [4.984,5.294] | -2.232 [-2.407,-2.057] | -0.977 [-1.154,-0.667] | 9.748 [9.663,9.862] | -3.719 [-3.939,-3.525] | -1.023 [-1.380,-0.895] |

#### AL4

| candidate | SNR@10% | gain vs B0 | gain vs PRG | SNR@1% | gain vs B0 | gain vs PRG |
|---|---:|---:|---:|---:|---:|---:|
| `C300_B0_QC` | -1.450 [-1.506,-1.411] | +0.000 [+0.000,+0.000] | +0.474 [+0.376,+0.588] | 0.805 [0.711,0.878] | +0.000 [+0.000,+0.000] | +1.311 [+1.156,+1.460] |
| `C300_S0_SIDON` | -1.389 [-1.464,-1.314] | -0.061 [-0.152,+0.022] | +0.413 [+0.297,+0.539] | 1.041 [0.914,1.111] | -0.236 [-0.358,-0.086] | +1.075 [+0.922,+1.251] |
| `C300_AP_RMS_T1` | -1.329 [-1.371,-1.277] | -0.122 [-0.196,-0.064] | +0.352 [+0.251,+0.462] | 1.248 [1.158,1.325] | -0.442 [-0.571,-0.323] | +0.868 [+0.713,+1.023] |
| `C300_AP_TEPS_T1` | -1.237 [-1.297,-1.193] | -0.214 [-0.282,-0.139] | +0.260 [+0.163,+0.384] | 0.882 [0.820,0.940] | -0.077 [-0.186,+0.020] | +1.234 [+1.090,+1.376] |
| `C300_GEO_T1_CTRL_BAL` | -1.140 [-1.187,-1.099] | -0.310 [-0.377,-0.248] | +0.163 [+0.069,+0.278] | 1.076 [0.972,1.156] | -0.271 [-0.397,-0.141] | +1.040 [+0.895,+1.203] |
| `C300_GEO_T1_CTRL_PAIR` | -1.070 [-1.108,-1.021] | -0.380 [-0.448,-0.322] | +0.094 [-0.007,+0.199] | 1.406 [1.318,1.508] | -0.601 [-0.735,-0.483] | +0.710 [+0.555,+0.864] |
| `C300_GEO_T1_CTRL_FOLD` | -1.241 [-1.296,-1.181] | -0.209 [-0.290,-0.141] | +0.265 [+0.163,+0.380] | 1.016 [0.939,1.063] | -0.210 [-0.314,-0.106] | +1.101 [+0.955,+1.241] |
| `C300_PRG_DFT4_6REG` | -0.977 [-1.060,-0.878] | -0.474 [-0.588,-0.376] | +0.000 [+0.000,+0.000] | 2.116 [1.979,2.251] | -1.311 [-1.460,-1.156] | +0.000 [+0.000,+0.000] |
| `C300_SMALL_CDD_QSTEP0P25_MATCHED_CDD` | 0.231 [0.093,0.343] | -1.682 [-1.806,-1.534] | -1.208 [-1.357,-1.035] | 3.951 [3.834,4.128] | -3.146 [-3.351,-3.011] | -1.835 [-2.055,-1.661] |
| `C300_SMALL_CDD_QSTEP0P25_TRANSPARENT_CDD` | 0.256 [0.090,0.328] | -1.706 [-1.797,-1.537] | -1.233 [-1.348,-1.041] | 3.854 [3.627,4.084] | -3.049 [-3.324,-2.809] | -1.738 [-2.020,-1.472] |

##### AL4 仿真结果分析

在 1% BLER 处，small CDD non-transparent（matched $R_g$）与 transparent（仅使用 $R_{phy}$）的门限点估计分别为 3.951 dB `[3.834,4.128]` 和 3.854 dB `[3.627,4.084]`，因此结果中 transparent 表观领先 0.097 dB。然而，对两条独立样本曲线直接 bootstrap 得到 $\Delta=\mathrm{SNR}_{matched}-\mathrm{SNR}_{transparent}=+0.097$ dB，95% 区间为 `[-0.178,+0.406]` dB；在 10% BLER 处同一定义的差值为 -0.025 dB，区间 `[-0.183,+0.187]` dB。两个区间都跨越 0，故本轮数据不支持“transparent 真实性能优于 non-transparent”的结论。

原始点说明了这一表观排序的来源：3.75 dB 时，matched 与 transparent 的 BLER 分别为 200/16900=1.183% `[1.031%,1.358%]` 和 201/19900=1.010% `[0.880%,1.159%]`，transparent 的点估计较低；到 4.00 dB 时，两者分别为 201/21000=0.957% `[0.834%,1.098%]` 和 200/20400=0.980% `[0.854%,1.125%]`，次序已经反转，而且各点 Wilson 区间均重叠。isotonic 与局部插值把这一独立抽样造成的交叉反映为约 0.1 dB 的门限次序，不应解释为透明接收带来的增益。

同 SNR 的信道估计结果也不支持 transparent 具有估计优势：在 3.75 dB，matched/transparent 的 CE NMSE 为 -12.001/-11.985 dB；在 4.00 dB 为 -12.184/-12.112 dB，matched 均略好。二者 zero-noise CE floor 更分别为 -122.45 dB 与 -48.14 dB，说明使用真实 $R_g$ 的 matched LMMSE 在模型层面保留了更低的渐近估计误差。由于小时延引入的 $R_g$ 与 $R_{phy}$ 差异在当前有限 SNR 区仍较小，BLER 又由独立样本估计，约 0.1 dB 的差别被统计不确定度与门限插值主导。若要分辨这一量级的差异，应在 3.75～4.25 dB 附近使用 matched/transparent 共用的信道、噪声和码字样本，或把每点目标误块数由约 200 提高到约 1000；在当前证据下只能结论为二者 1% BLER 性能不可分辨。

### 7.4 Pilot 可辨识性、二维 CE NMSE 与 zero-noise floor

表中 `M/P/T` 分别表示 `matched_effective`、transparent `physical_prg` 和 transparent `physical_fullband`。`rank/cond` 是 pilot 上 precoder 矩阵的列秩/条件数；cycling 的 AL1/2 只使用 1/2 个 DFT4 向量，故 rank 为 1/2、condition 为 `inf`。`CE@ref` 取最接近本 candidate 1%目标的实测 SNR 点；两条未获得 1%目标的 AL1 GEO 取原始 BLER 最接近 1%的点。由于各行 reference SNR 不同，CE 数字用于核对每条曲线的量级，不能直接作为同 SNR candidate 排名。

| AL | candidate | mode | rank/cond | CE floor / dB | ref SNR / BLER | CE@ref / dB |
|---:|---|---:|---:|---:|---:|---:|
| 1 | `B0_QC` | M | 4/6.47 | -59.40 | 13.75/0.0091 | -17.09 |
| 1 | `S0_SIDON` | M | 4/1 | -16.90 | 12.75/0.0103 | -12.68 |
| 1 | `AP_RMS_T1` | M | 4/47.6 | -70.77 | 15.00/0.0092 | -18.06 |
| 1 | `AP_TEPS_T1` | M | 4/1.07 | -21.31 | 12.00/0.0109 | -13.64 |
| 1 | `GEO_BAL_PAIR` | M | 4/1 | -12.69 | 18.25/0.0160 | -11.27 |
| 1 | `GEO_FOLD` | M | 4/1 | -13.67 | 18.25/0.0100 | -12.12 |
| 1 | `PRG_DFT4_6REG` | P | 1/inf | -96.13 | 17.00/0.0092 | -19.23 |
| 1 | `SMALL_MATCHED` | M | 4/5.47e5 | -93.00 | 16.75/0.0108 | -18.91 |
| 1 | `SMALL_TRANSPARENT` | T | 4/5.47e5 | -87.67 | 17.00/0.0099 | -19.16 |
| 2 | `B0_QC` | M | 4/1.22 | -90.87 | 6.00/0.0102 | -12.03 |
| 2 | `S0_SIDON` | M | 4/1 | -55.19 | 5.50/0.0110 | -10.78 |
| 2 | `AP_RMS_T1` | M | 4/3.67 | -101.18 | 7.00/0.0101 | -13.14 |
| 2 | `AP_TEPS_T1` | M | 4/1.07 | -25.94 | 5.50/0.0102 | -9.92 |
| 2 | `GEO_BAL` | M | 4/1 | -13.56 | 6.25/0.0100 | -8.85 |
| 2 | `GEO_PAIR` | M | 4/1 | -12.95 | 6.50/0.0101 | -8.71 |
| 2 | `GEO_FOLD` | M | 4/1 | -15.65 | 6.00/0.0104 | -9.19 |
| 2 | `PRG_DFT4_6REG` | P | 2/inf | -96.13 | 8.75/0.0100 | -14.18 |
| 2 | `SMALL_MATCHED` | M | 4/6.37e4 | -113.27 | 9.75/0.0108 | -14.88 |
| 2 | `SMALL_TRANSPARENT` | T | 4/6.37e4 | -79.57 | 9.75/0.0100 | -14.89 |
| 4 | `B0_QC` | M | 4/1.18 | -102.71 | 0.75/0.0107 | -8.69 |
| 4 | `S0_SIDON` | M | 4/1 | -104.95 | 1.00/0.0105 | -9.07 |
| 4 | `AP_RMS_T1` | M | 4/1.37 | -111.75 | 1.25/0.0100 | -9.72 |
| 4 | `AP_TEPS_T1` | M | 4/1.05 | -31.15 | 1.00/0.0080 | -7.83 |
| 4 | `GEO_BAL` | M | 4/1 | -14.68 | 1.00/0.0109 | -7.18 |
| 4 | `GEO_PAIR` | M | 4/1 | -13.19 | 1.50/0.0089 | -7.14 |
| 4 | `GEO_FOLD` | M | 4/1 | -23.82 | 1.00/0.0103 | -7.67 |
| 4 | `PRG_DFT4_6REG` | P | 4/1 | -96.13 | 2.00/0.0112 | -10.24 |
| 4 | `SMALL_MATCHED` | M | 4/7.78e3 | -122.45 | 4.00/0.0096 | -12.18 |
| 4 | `SMALL_TRANSPARENT` | T | 4/7.78e3 | -48.14 | 3.75/0.0101 | -11.98 |

#### B0_QC 与 precoder cycling 的 CE 差异解释

先澄清处理范围：只有 AL1 时整个 candidate 恰好只有一个 `6 REG = 3 RB × 2 symbols` bundle，因此 `B0_QC` 与 `PRG_DFT4_6REG` 都使用该 candidate 内的全部 18 个 DMRS 观测估计 54 个 data RE。AL2/4 分别有 2/4 个 bundle；`B0_QC` 使用全 candidate 的 matched 二维时频 LMMSE，而 cycling 的权重矩阵按 36-subcarrier bundle 分块、每块联合两个 symbol，但不跨 bundle。因而“二者都在相同全带宽做频域 MMSE”只对 AL1 成立；本实现实际是二维时频 LMMSE，不是把两个 DMRS symbol 平均后只做一维频域估计。

原始 CE 字段采用 $10\log_{10}E[\|e_b\|^2/\|g_b\|^2]$，即先按每个 trial 的真实 data-channel 能量归一化，再对 trial 求均值。预扫中选取共同 SNR 后，`B0_QC/PRG_DFT4_6REG` 的 CE NMSE 分别为：AL1、10 dB 时 `-14.390/-13.755 dB`；AL2、4 dB 时 `-10.726/-10.737 dB`；AL4、2 dB 时 `-9.746/-10.220 dB`。这些预扫点每点只有 300 trials，只用于解释曲线形状；正式网格中两者唯一的共同点是 AL4、-1.25 dB，对应 `-7.237/-8.105 dB`，方向与 AL4 预扫一致。因此 AL1 的“小幅领先”和 AL2 的“接近”不应从各 candidate 不同 reference SNR 的 7.4 表格直接推出。

用户关于协方差的判断方向是对的，但还需加上上述归一化统计量。AL1 的 cycling 实际没有发生向量轮换：整个 candidate 使用一个固定 DFT4 向量，等效信道保持 TDL-C 在 3 RB 内较强的频率相关，整块同时落入深衰落时 $\|g_b\|^2$ 较小，per-trial ratio 会放大该 trial。`B0_QC` 的频率相位斜率引入人工频率选择性，使 54 个 data RE 的块能量更分散，减弱这种小分母效应；这是与 AL1 预扫现象一致的机理推断，不是“B0 的无条件 LMMSE MSE 更低”的证明。按冻结的 $R_{PP}$、$R_{DP}$ 直接计算 ratio-of-expectations $10\log_{10}(E\|e\|^2/E\|g\|^2)$，在上述 AL1/AL2 SNR 以及 AL4 的正式共同 SNR 处，`B0_QC/PRG_DFT4_6REG` 分别为 `-15.599/-17.169 dB`、`-11.335/-12.565 dB`、`-7.727/-8.980 dB`，反而均是 cycling 较低。这说明 AL1 图上的轻微反序主要涉及不同协方差下的衰落能量分布与 per-trial 归一化，而不是全带观测优势。

AL 增大后，cycling 出现 2/4 个采用不同 DFT4 向量的 bundle，归一化分母汇总多个 bundle 的信道能量，小分母效应随之减弱；每个 bundle 内仍是较平滑的物理 TDL-C 协方差，并有固定 18 个 DMRS 观测预测 54 个 data RE。这解释了 cycling 的 zero-noise floor 在三个 AL 基本固定为 `-96.13 dB`，以及它相对 AL1 改善后在 AL2 接近、AL4 低于 B0 的趋势。表中的 cycling pilot rank 1/2/4 只是 candidate 出现了多少个 DFT 向量；接收机估计的是每 bundle 的标量等效信道，所以 AL1 的 rank 1 本身不是 CE 缺陷。

最后，全带 LMMSE 的增益只对“同一个真实协方差、同一个发射波形下，全带权重与强制分块权重”的反事实比较有必然性；当前 `B0_QC` 与 cycling 同时改变了发射预编码和真实等效协方差。`B0_QC` 的人工时延降低远距离频点的条件相关性，使跨 bundle DMRS 的边际信息有限；cycling 虽然分块，却在每个较平滑的局部信道上做匹配估计。因此 B0 的全带处理增益可以被其更难预测的等效协方差抵消。若要单独量化“跨 bundle 联合处理增益”，应保持 `B0_QC` 的同一真实信道与 paired trials，只把接收机从 matched fullband 改成 forced-PRG；本轮没有做该消融，不能由当前两条不同波形的曲线反推该增益。

### 7.5 结果分析

事实：AL1 在可估计候选中 `AP_TEPS_T1` 最好，10%/1% 相对 B0 分别改善 0.532 dB `[0.391,0.648]` 和 1.531 dB `[1.346,1.747]`；`S0_SIDON` 次之，改善 0.470/0.866 dB。AL1 在 1% 附近的统计抖动与 GEO error floor 已在本节 AL1 分析中结合原始计数和 CE 结果说明。

事实：AL2 的 `S0_SIDON` 与 `AP_TEPS_T1` 最有竞争力。S0 在 10%/1% 相对 B0 改善 0.161 dB `[0.064,0.278]`/0.444 dB `[0.251,0.653]`；AP-$T_{0.01}$ 的 10% 差异不可分辨，但 1% 改善 0.513 dB `[0.314,0.716]`。两者各自区间重叠，当前表不提供二者直接差值的 paired 结论。

事实：AL4 的最低点估计是 B0。`AP_TEPS_T1` 在 1% 只差 0.077 dB，区间 `[-0.186,0.020]` 跨 0，尚不可分辨；S0 在 1% 比 B0 差 0.236 dB `[-0.358,-0.086]`。这再次说明“严格 Sidon”或较大 pair/fold gap 不是 BLER 的单调充分条件。

三个 AL 中 transparent cycling 均比 B0 差；AL1 由于只有一个 bundle，实际没有 cycling，且每 bundle 估计不能从候选外借用 pilot。小时延 CDD 的 matched/transparent 点估计彼此很近：AL1 对 PRG 的 10%/1% 差值区间均跨 0；AL2 两种接收口径的差异也很小；AL4 的 0.097 dB 表观反序已在本节 AL4 分析中通过直接差值区间、原始计数和同 SNR CE 结果说明。总体上，本轮有限 SNR 工作点仍由热噪声、短码译码事件和有限样本不确定度主导，不能仅凭 zero-noise CE floor 预测 BLER 排名。

机理推断：候选排序由完整的 $R_{PP}$ 与 $R_{DP}$、资源映射和 Polar rate matching 共同决定。GEO 只优化整数 pair/fold 的最小圆周间距，不优化 TDL-C PDP 加权的 data-from-pilot 条件协方差；所以 rank 均为 4、condition 约 1 也不保证低 CE floor。AL1 GEO 的 floor 直接证明这种几何摘要遗漏了关键预测结构。AL1/2/4 同时改变 $K_f$、$E$、码率、bundle 数和 DFT 覆盖，故不能把排序变化单独归因于时延规则或移动性。

### 7.6 审计、证据与限制

| AL | SNR points | trials | errors | adjacent raw BLER increases | 50k 且 errors<200 |
|---:|---:|---:|---:|---:|---:|
| 1 | 109 | 1,612,200 | 64,748 | 17 | 0 |
| 2 | 115 | 2,030,600 | 63,145 | 0 | 4 |
| 4 | 109 | 1,985,600 | 62,178 | 0 | 5 |
| 合计 | 333 | 5,628,400 | 190,071 | 17 | 9 |

最终分析文件以表中 333 点为准。所有逐 trial flags 的长度与错误和都和 CSV 一致，CE 线性均值复算到 dB 的最大误差为 0；56 个已估计目标的原始 bracket 最大为 0.25 dB，bootstrap 有效 replicate 为 1,974–2,000。AL1 的 17 个相邻 BLER increases 均与 GEO 高 SNR floor 或有限样本波动有关，原始数据保留，目标拟合使用预定递减 isotonic。

正式前冻结的主网格满足每点预算，但预扫位置误差使部分实际 target 落在两段局部网格之间；为满足 0.25 dB bracket，又按固定缺口规则后验追加 31 个点。该补点选择以及 prescan/formal 在相同 SNR 下复用了相同 seed 流的前 300 trials，会产生未计入普通 Bernoulli bootstrap 的网格选择依赖；因此 CI 只描述冻结后各点的有限 Bernoulli 样本误差，不包含预扫或补点选择不确定性。这是本轮明确的统计限制。

最终仍未完成的两个目标仅为 AL1 `GEO_T1_CTRL_BAL_PAIR` 和 `GEO_T1_CTRL_FOLD` 的 1%；二者不报告目标 SNR或相对增益。当前模型只覆盖两相邻 CORESET symbols 内的 3 km/h 时间变化，不含跨 slot 预测、CFO/ICI、同步误差、候选散列或空间相关性。

精确证据位于 `outputs/experiment031_pdcch_cdd/20260907_c300_4tx_2sym/`：`candidate_receipt.json`、`geo_exhaustive_receipt.json`、`pdp_t0p01_receipt.json`、`prescan_freeze_receipt.json`、`bracket_completion_receipt.json` 和 `environment_receipt.json`；每个 `formal/al*/<candidate>/` 目录含展开配置、`bler_points.csv`、`run_metadata.json`、日志及逐 trial flags/CE；最终汇总为 `analysis/formal_points.csv`、`diagnostics.csv`、`ce_near_1pct.csv`、`target_snr.csv`、`target_gains.csv` 与 `analysis_metadata.json`。最终补点配置 SHA-256 为 AL1 `6c26ba6c4f942ef0000e2637198a523d3bb7bd702aa6ad5c181d12bd23e8a988`、AL2 `53046887dafa83c00bbe6c40ed00b3ee838cae009648fc91d8e914ab62589483`、AL4 `d47435e2f3b7b80456db8b796b9227add7eb1c8a75088e8a0bf0232f2b035d23`；基准提交仍为 `83b25bc2ea77324ae98f7dc67136fc1e93b1f32c`，实现和结果是未提交工作区增量。

运行环境为 Windows、Python 3.11.9、NumPy 1.26.4、SciPy 1.15.3、Matplotlib 3.10.3、TensorFlow 2.15.1、Sionna 1.0.2。复现与验证命令为：

```powershell
python tools/prepare_plan031_c300_formal.py
python tools/verify_plan031_c300_geometry.py
python tools/run_plan031_candidates.py --config configs/pdcch_result031_c300_al1_formal.yaml --config configs/pdcch_result031_c300_al2_formal.yaml --config configs/pdcch_result031_c300_al4_formal.yaml --stage run --max-workers 3
python tools/prepare_plan031_c300_bracket_completion.py
python tools/run_plan031_candidates.py --config configs/pdcch_result031_c300_al1_bracket_completion.yaml --config configs/pdcch_result031_c300_al2_bracket_completion.yaml --config configs/pdcch_result031_c300_al4_bracket_completion.yaml --stage run --max-workers 3
python tools/analyze_plan031_c300.py --bootstrap-repeats 2000
python -m pytest tests/test_pdcch.py tests/test_rmmse_time_frequency.py tests/test_plan031_analysis.py -q
```

最终为 `31 passed`；14 条 warning 均来自 Matplotlib/PyParsing deprecated API，不影响数值结果。结果尚待研究者确认，因此未更新 `KNOWLEDGE.md`/`GOALS.md`，也未创建 Git checkpoint。

## 8. C300 2-symbol ideal-CSI BLER 补充（2026-09-10）

### 8.1 配置与判据

本补充保持第 7 节的 TDL-C 300 ns、4Tx/1Rx、3 km/h、两个相邻 CORESET symbols、资源映射、Polar 编译码与候选集合不变，只把接收机切换为 `channel_estimation: ideal`。解调直接使用每个 data RE 的真实等效信道做 coherent MRC；DMRS 仍发送但不参与估计，不构造 LMMSE filter，记录的 CE NMSE 精确为 0。小时延 CDD 的 matched/transparent 发射波形在此接收机下完全相同，因此合并为 `C300_SMALL_CDD_QSTEP0P25_IDEAL_CSI` 一条曲线。

预扫采用 -12:2:18 dB、每点 300 trials；正式网格为 0.25 dB，停止规则为至少 10,000 trials 且累计 200 errors，最多 50,000 trials。正式后仅向缺失方向追加了 AL2 `AP_RMS_T1` 的 5.75 dB、AL2 `PRG_DFT4_6REG` 的 6.5/6.75/7.0 dB，以及 AL4 small-CDD 的 2.25/2.5 dB；既有点由绝对 trial 流断点续跑机制跳过，未重算或相加重复 trial。

### 8.2 目标 SNR 与相对增益

下表增益定义为基准目标 SNR 减候选目标 SNR，正值表示候选更好。全部 52 个 10%/1% 目标都有不宽于 0.25 dB 的相邻原始括点；目标 SNR 由预定递减 isotonic 曲线在括点内插值得到，完整 95% bootstrap CI 见 `analysis/target_snr.csv` 与 `analysis/target_gains.csv`。

| AL | candidate | 10% SNR | gain vs B0 | gain vs PRG | 1% SNR | gain vs B0 | gain vs PRG |
|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | `B0_QC` | 7.575 | 0.000 | 2.223 | 12.002 | 0.000 | 3.505 |
| 1 | `S0_SIDON` | 6.235 | 1.340 | 3.563 | 9.234 | 2.768 | 6.274 |
| 1 | `AP_RMS_T1` | 8.317 | -0.741 | 1.481 | 13.280 | -1.278 | 2.228 |
| 1 | `AP_TEPS_T1` | 6.390 | 1.185 | 3.408 | 9.518 | 2.485 | 5.990 |
| 1 | `GEO_T1_CTRL_BAL_PAIR` | 6.124 | 1.451 | 3.674 | 9.098 | 2.904 | 6.410 |
| 1 | `GEO_T1_CTRL_FOLD` | 6.172 | 1.404 | 3.626 | 9.146 | 2.856 | 6.362 |
| 1 | `PRG_DFT4_6REG` | 9.798 | -2.223 | 0.000 | 15.508 | -3.505 | 0.000 |
| 1 | `SMALL_CDD_IDEAL_CSI` | 9.976 | -2.401 | -0.179 | 15.545 | -3.543 | -0.038 |
| 2 | `B0_QC` | 1.266 | 0.000 | 1.589 | 4.356 | 0.000 | 2.434 |
| 2 | `S0_SIDON` | 1.001 | 0.265 | 1.854 | 3.609 | 0.748 | 3.182 |
| 2 | `AP_RMS_T1` | 1.966 | -0.700 | 0.889 | 5.443 | -1.087 | 1.347 |
| 2 | `AP_TEPS_T1` | 0.814 | 0.452 | 2.042 | 3.196 | 1.160 | 3.594 |
| 2 | `GEO_T1_CTRL_BAL` | 0.736 | 0.530 | 2.119 | 3.165 | 1.191 | 3.625 |
| 2 | `GEO_T1_CTRL_PAIR` | 0.753 | 0.513 | 2.102 | 3.281 | 1.075 | 3.509 |
| 2 | `GEO_T1_CTRL_FOLD` | 0.822 | 0.443 | 2.033 | 3.221 | 1.135 | 3.569 |
| 2 | `PRG_DFT4_6REG` | 2.855 | -1.589 | 0.000 | 6.791 | -2.434 | 0.000 |
| 2 | `SMALL_CDD_IDEAL_CSI` | 4.089 | -2.823 | -1.234 | 8.378 | -4.022 | -1.588 |
| 4 | `B0_QC` | -3.274 | 0.000 | 0.742 | -1.140 | 0.000 | 1.355 |
| 4 | `S0_SIDON` | -3.166 | -0.108 | 0.634 | -0.930 | -0.210 | 1.144 |
| 4 | `AP_RMS_T1` | -3.003 | -0.271 | 0.471 | -0.555 | -0.586 | 0.769 |
| 4 | `AP_TEPS_T1` | -3.446 | 0.172 | 0.914 | -1.502 | 0.362 | 1.716 |
| 4 | `GEO_T1_CTRL_BAL` | -3.456 | 0.182 | 0.923 | -1.533 | 0.392 | 1.747 |
| 4 | `GEO_T1_CTRL_PAIR` | -3.483 | 0.209 | 0.951 | -1.545 | 0.405 | 1.759 |
| 4 | `GEO_T1_CTRL_FOLD` | -3.537 | 0.263 | 1.005 | -1.611 | 0.471 | 1.825 |
| 4 | `PRG_DFT4_6REG` | -2.532 | -0.742 | 0.000 | 0.214 | -1.355 | 0.000 |
| 4 | `SMALL_CDD_IDEAL_CSI` | -1.019 | -2.255 | -1.513 | 2.332 | -3.472 | -2.118 |

### 8.3 信道估计代价与结论边界

定义 $\Delta_{\rm CE}=\mathrm{SNR}_{\rm estimated}-\mathrm{SNR}_{\rm ideal}$，正值表示第 7 节的估计 CSI 接收机需要更多 SNR。下表给出 10%/1% 的点估计；small-CDD 依次列 matched/transparent 两种历史接收口径。AL1 两个 GEO 的历史 estimated-CSI 1% 目标不合格，故对应 $\Delta_{\rm CE}$ 不报告。

| AL | candidate | ΔCE @10% (dB) | ΔCE @1% (dB) |
|---:|---|---:|---:|
| 1 | `B0_QC` | 1.398 | 1.653 |
| 1 | `S0_SIDON` | 2.269 | 3.555 |
| 1 | `AP_RMS_T1` | 1.234 | 1.613 |
| 1 | `AP_TEPS_T1` | 2.051 | 2.606 |
| 1 | `GEO_BAL_PAIR` | 3.411 | — |
| 1 | `GEO_FOLD` | 2.967 | — |
| 1 | `PRG_DFT4_6REG` | 1.048 | 1.416 |
| 1 | `SMALL_CDD` matched / transparent | 0.843 / 0.893 | 1.314 / 1.446 |
| 2 | `B0_QC` | 1.639 | 1.673 |
| 2 | `S0_SIDON` | 1.744 | 1.976 |
| 2 | `AP_RMS_T1` | 1.480 | 1.576 |
| 2 | `AP_TEPS_T1` | 2.068 | 2.320 |
| 2 | `GEO_BAL` | 2.459 | 3.088 |
| 2 | `GEO_PAIR` | 2.499 | 3.230 |
| 2 | `GEO_FOLD` | 2.166 | 2.806 |
| 2 | `PRG_DFT4_6REG` | 1.305 | 1.934 |
| 2 | `SMALL_CDD` matched / transparent | 1.077 / 1.048 | 1.491 / 1.370 |
| 4 | `B0_QC` | 1.824 | 1.946 |
| 4 | `S0_SIDON` | 1.776 | 1.971 |
| 4 | `AP_RMS_T1` | 1.675 | 1.802 |
| 4 | `AP_TEPS_T1` | 2.209 | 2.384 |
| 4 | `GEO_BAL` | 2.316 | 2.609 |
| 4 | `GEO_PAIR` | 2.413 | 2.951 |
| 4 | `GEO_FOLD` | 2.295 | 2.626 |
| 4 | `PRG_DFT4_6REG` | 1.555 | 1.902 |
| 4 | `SMALL_CDD` matched / transparent | 1.250 / 1.275 | 1.619 / 1.522 |

事实：去除信道估计误差后，AL1/2 的最佳点估计仍来自 GEO/AP/Sidon 类 CDD，而 AL4 的最佳点估计变为 `GEO_FOLD`；`PRG_DFT4_6REG` 与 small-CDD 在三个 AL 均未超过 B0。所有可比较的 $\Delta_{\rm CE}$ 都为正，范围为 0.843–3.555 dB，说明第 7 节的 estimated-CSI 链路确有显著估计代价。

解释边界：estimated 与 ideal 的差值同时包含 noisy DMRS、二维插值/外推、候选或 bundle 处理范围及其与译码事件的相互作用；本轮不能把 $\Delta_{\rm CE}$ 单独归因于频率估计、时间估计或 PRG 分块。候选在 ideal CSI 下仍有明显排序，证明排序不只来自 CE NMSE，但这也不是任一时延设计规则的普适最优性证明。

### 8.4 审计、证据与复现

| AL | formal points | trials | errors | adjacent raw BLER increases | 50k 且 errors<200 |
|---:|---:|---:|---:|---:|---:|
| 1 | 128 | 2,405,900 | 78,133 | 4 | 8 |
| 2 | 148 | 2,639,600 | 95,600 | 2 | 13 |
| 4 | 146 | 2,851,000 | 102,263 | 0 | 18 |
| 合计 | 422 | 7,896,500 | 275,996 | 6 | 39 |

所有逐 trial flags 的长度和求和均与 CSV 一致；CE 线性均值复算误差为 0，ideal CE NMSE 最大绝对误差为 0，52/52 个目标合格。验证覆盖 ideal CE 逐元素精确为零、无 LMMSE filter、无噪声 DCI 解码、coherent-MRC 分母/噪声、非零 Doppler 下两 symbol 信道不同、small-CDD 两接收模式等价、candidate/AL/stage 独立随机流，以及 estimated-CSI 回归；定向测试结果为 `37 passed`，14 条 warning 均来自 Matplotlib/PyParsing deprecated API。

精确证据位于 `outputs/experiment031_pdcch_cdd/20260910_c300_4tx_2sym_ideal_csi/`。每个 `formal/al*/<candidate>/` 目录保存展开配置、逐点 CSV、日志和 trial flags；汇总文件为 `analysis/formal_points.csv`、`diagnostics.csv`、`target_snr.csv`、`target_gains.csv`、`delta_ce.csv` 与 `analysis_metadata.json`。最终配置 SHA-256 为 AL1 `a45923640db3a3a7f311895efd9f88f354a77dcfba68310f2eabbd071dda5130`、AL2 `b678efa03d2de16362451a19bb782d2c47b68f538ef8e857516f3c6d62a1e546`、AL4 `ce23c706445d109403b4516537b7237d4cb9b5414e8b077ae05bff1e0650df23`；基准提交为 `83b25bc2ea77324ae98f7dc67136fc1e93b1f32c`，实现与结果仍是未提交工作区增量。

复现命令为：

```powershell
python tools/prepare_plan031_c300_ideal.py --stage prescan
python tools/run_plan031_candidates.py --config configs/pdcch_result031_c300_ideal_al1_prescan.yaml --config configs/pdcch_result031_c300_ideal_al2_prescan.yaml --config configs/pdcch_result031_c300_ideal_al4_prescan.yaml --stage run --max-workers 6
python tools/prepare_plan031_c300_ideal.py --stage formal
python tools/run_plan031_candidates.py --config configs/pdcch_result031_c300_ideal_al1_formal.yaml --config configs/pdcch_result031_c300_ideal_al2_formal.yaml --config configs/pdcch_result031_c300_ideal_al4_formal.yaml --stage run --max-workers 6
python tools/analyze_plan031_c300_ideal.py --bootstrap-repeats 2000
python -m pytest tests/test_pdcch.py tests/test_rmmse_time_frequency.py tests/test_plan031_analysis.py -q
```

本补充结果尚待研究者确认，因此未更新 `KNOWLEDGE.md`/`GOALS.md`，也未创建 Git checkpoint。

## 9. 严格 Sidon 候选搜索补充（2026-09-10）

### 9.1 范围、筛选和统计口径

本补充在既有六个 `(Tx, AL)` 场景中搜索表列严格 Sidon 候选：A100 为 static TDL-A 100 ns、8Tx/1Rx、1 symbol、AL2/4/8；C300 为 TDL-C 300 ns、4Tx/1Rx、2 symbols、3 km/h、AL1/2/4。其余 PDCCH、功率、DMRS、资源映射和 matched effective LMMSE 定义与第 2、7 节一致。所有候选都满足完整 pilot rank、condition 约为 1、模 $K$ 无序二元和唯一和 DMRS fold residue 唯一；完整整数坐标与 ns 换算见 `research/plan-031-PDCCH-CDD时延-BLER.md`。

执行链为 48 候选、300 trials/点预扫；5 个候选按预注册的 BLER bracket、Wilson 下界、CE floor 距离和高端 BLER 平台四项联合判据停止；其余 43 个候选做 3,000 trials/点粗确认；18 个仍可能最优者进入正式细扫描。正式扫描和边界补点合计 172 个 candidate/SNR 点，每点至少 10,000 trials、达到 200 errors 后停止、上限 50,000。全局 seed 为 `20260908`，正式细扫描随机流命名空间为 `plan031_strict_sidon_final_fine_v1`。运行环境为 Python 3.11.9、NumPy 1.26.4、SciPy 1.15.3、TensorFlow 2.15.1、Sionna 1.0.2，CPU-only；基准提交为 `83b25bc2ea77324ae98f7dc67136fc1e93b1f32c`，本补充仍是未提交工作区增量。目标 SNR 使用 Jeffreys 平滑、trial 数加权的递减 isotonic fit 和局部 `log10(BLER)` 插值；95% 区间来自 10,000 次独立 Bernoulli bootstrap。36/36 个目标均有宽度不超过 0.25 dB 的原始双侧 bracket，且至少 95% bootstrap 重采样保持 bracket。

### 9.2 目标 SNR 与相对 S0 增益

单位均为 dB；格式为点估计 `[95% CI]`。gain 定义为 `S0 目标 SNR − candidate 目标 SNR`，正值表示优于同场景 S0。A100 AL2 和 C300 AL1 的 S0 来自既有正式曲线，其余场景的 candidate 01 即 S0；均按独立 bootstrap 比较。

| scene | candidate | SNR@10% | gain@10% vs S0 | SNR@1% | gain@1% vs S0 |
|---|---|---:|---:|---:|---:|
| A100 AL2 | 02 | 2.841 `[2.788,2.885]` | 0.196 `[0.118,0.278]` | 4.868 `[4.775,4.942]` | 0.458 `[0.337,0.569]` |
| A100 AL4 | 01 (S0) | -1.250 `[-1.299,-1.193]` | 0 | 0.709 `[0.606,0.789]` | 0 |
| A100 AL4 | 02 | -1.025 `[-1.076,-0.965]` | -0.225 `[-0.304,-0.148]` | 0.998 `[0.914,1.070]` | -0.289 `[-0.418,-0.165]` |
| A100 AL4 | 03 | -1.066 `[-1.108,-1.013]` | -0.183 `[-0.255,-0.112]` | 0.840 `[0.765,0.890]` | -0.130 `[-0.250,-0.021]` |
| A100 AL4 | 06 | -1.130 `[-1.171,-1.088]` | -0.119 `[-0.185,-0.049]` | 0.807 `[0.701,0.896]` | -0.098 `[-0.239,0.044]` |
| A100 AL8 | 01 (S0) | -4.493 `[-4.550,-4.455]` | 0 | -2.720 `[-2.811,-2.674]` | 0 |
| A100 AL8 | 02 | -4.286 `[-4.332,-4.233]` | -0.207 `[-0.283,-0.147]` | -2.570 `[-2.639,-2.486]` | -0.150 `[-0.273,-0.066]` |
| A100 AL8 | 05 | -4.237 `[-4.286,-4.191]` | -0.257 `[-0.330,-0.192]` | -2.431 `[-2.524,-2.369]` | -0.289 `[-0.396,-0.184]` |
| A100 AL8 | 06 | -4.419 `[-4.451,-4.392]` | -0.074 `[-0.138,-0.024]` | -2.592 `[-2.698,-2.486]` | -0.129 `[-0.262,-0.016]` |
| C300 AL1 | 04 | 8.197 `[8.134,8.281]` | 0.307 `[0.181,0.444]` | 11.805 `[11.694,11.912]` | 0.984 `[0.812,1.172]` |
| C300 AL2 | 01 (S0) | 2.755 `[2.667,2.817]` | 0 | 5.641 `[5.498,5.766]` | 0 |
| C300 AL2 | 03 | 3.066 `[2.976,3.143]` | -0.311 `[-0.431,-0.198]` | 5.925 `[5.790,6.075]` | -0.285 `[-0.478,-0.105]` |
| C300 AL2 | 04 | 2.700 `[2.655,2.769]` | 0.055 `[-0.053,0.134]` | 5.548 `[5.386,5.767]` | 0.093 `[-0.165,0.297]` |
| C300 AL2 | 06 | 2.930 `[2.881,2.995]` | -0.175 `[-0.283,-0.094]` | 5.838 `[5.702,6.039]` | -0.197 `[-0.450,0.008]` |
| C300 AL4 | 01 (S0) | -1.395 `[-1.471,-1.320]` | 0 | 0.983 `[0.901,1.061]` | 0 |
| C300 AL4 | 02 | -1.030 `[-1.075,-0.972]` | -0.365 `[-0.458,-0.279]` | 1.374 `[1.243,1.488]` | -0.392 `[-0.525,-0.249]` |
| C300 AL4 | 03 | -1.168 `[-1.217,-1.127]` | -0.227 `[-0.311,-0.138]` | 0.976 `[0.889,1.075]` | 0.006 `[-0.129,0.128]` |
| C300 AL4 | 06 | -1.405 `[-1.458,-1.357]` | 0.010 `[-0.079,0.100]` | 0.817 `[0.730,0.872]` | 0.165 `[0.065,0.284]` |

### 9.3 最优判定与 CE 关系

1% 主判据的点估计领先者及 95% 判定如下。这里“并列”只表示 `candidate SNR − leader SNR` 的 95% bootstrap 区间包含 0，并非证明两者等效。

| scene | 点估计领先者；整数坐标 | 95% 判定 |
|---|---|---|
| A100 AL2 | 02；`[0,3,7,25,46,57,124,139]` | 分阶段筛选后的唯一领先者；相对 S0 改善 0.458 dB `[0.337,0.569]` |
| A100 AL4 | 01；`[0,1,3,7,12,20,30,65]` | 01 与 06 并列；06 比 01 多需 0.098 dB `[-0.044,0.239]` |
| A100 AL8 | 01；`[0,1,3,7,12,20,30,65]` | 唯一领先；最近的 06 多需 0.129 dB `[0.016,0.262]` |
| C300 AL1 | 04；`[0,1,4,6]` | 分阶段筛选后的唯一领先者；相对 S0 改善 0.984 dB `[0.812,1.172]` |
| C300 AL2 | 04；`[0,1,4,6]` | 04、01、06 并列；01/06 分别多需 0.093 `[-0.165,0.297]`/0.290 `[-0.019,0.580]` dB |
| C300 AL4 | 06；`[0,3,7,136]` | 唯一领先；相对 S0 改善 0.165 dB `[0.065,0.284]` |

同一 SNR 的 CE NMSE 抽样如下，mean 为先在线性域求均值再转 dB；median、q10/q90 为逐 trial 线性 NMSE 分位数转 dB。完整 q2.5/q97.5 也保存在原始汇总中。

| scene/SNR | candidate | mean | median | q10/q90 |
|---|---|---:|---:|---:|
| A100 AL4 / 0.75 | 01 | -6.904 | -7.160 | -8.936 / -5.304 |
| A100 AL4 / 0.75 | 06 | -6.502 | -6.739 | -8.516 / -4.935 |
| A100 AL8 / -2.50 | 01 | -6.114 | -6.322 | -7.970 / -4.632 |
| A100 AL8 / -2.50 | 06 | -5.665 | -5.860 | -7.426 / -4.222 |
| C300 AL2 / 5.75 | 01 | -10.965 | -11.615 | -14.389 / -8.654 |
| C300 AL2 / 5.75 | 04 | -11.074 | -11.735 | -14.542 / -8.751 |
| C300 AL2 / 5.75 | 06 | -9.460 | -10.029 | -12.679 / -7.229 |
| C300 AL4 / 0.75 | 01 | -8.882 | -9.414 | -11.813 / -6.775 |
| C300 AL4 / 0.75 | 06 | -7.989 | -8.333 | -10.543 / -6.139 |

事实：A100 AL4/AL8 和 C300 AL2 中，较低 CE NMSE 与较好的 BLER 大体一致；但 C300 AL4 的 candidate 01 在 0.75 dB 的平均 CE NMSE 比 06 低约 0.89 dB，1% BLER 却由 06 显著领先 0.165 dB。因此平均 CE NMSE 可解释部分排序，但不是 BLER 排序的充分统计量。严格 Sidon、较大 pair/fold gap 或较低 CE floor 也都不是跨场景单调最优条件。

### 9.4 审计、证据与限制

| scene | candidates | points | trials | errors | adjacent BLER increases | 50k 且 errors<200 |
|---|---:|---:|---:|---:|---:|---:|
| A100 AL2 | 1 | 6 | 85,300 | 3,847 | 0 | 0 |
| A100 AL4 | 4 | 36 | 528,300 | 21,250 | 0 | 1 |
| A100 AL8 | 4 | 31 | 558,700 | 16,712 | 0 | 2 |
| C300 AL1 | 1 | 16 | 189,100 | 7,573 | 0 | 0 |
| C300 AL2 | 4 | 44 | 659,500 | 23,815 | 0 | 0 |
| C300 AL4 | 4 | 39 | 623,100 | 19,569 | 0 | 0 |
| 合计 | 18 | 172 | 2,644,000 | 92,766 | 0 | 3 |

所有逐 trial error flags 的长度与求和均和 CSV 一致，CE 线性均值复算最大误差为 0 dB；原始曲线没有相邻 BLER 回升。完整证据位于 `outputs/experiment031_pdcch_cdd/20260908_strict_sidon_search/`：预扫筛除记录在 `screening/`，粗筛记录在 `coarse_analysis/`，正式原始数据在 `fine/` 与 `fine_supplement/`，最终汇总为 `fine_analysis/formal_points.csv`、`nmse_trial_summary.csv`、`target_snr.csv`、`target_bracket_points.csv`、`winner_summary.csv`、`pairwise_1pct.csv`、`gain_vs_s0.csv`、`diagnostics.csv`、`curve_styles.json` 和 `analysis_metadata.json`。后者保存全部配置与源码 SHA-256、Git HEAD 和统计参数。

限制：A100 AL2 和 C300 AL1 的正式细网格只保留了粗确认后唯一可能领先的候选，因此其“唯一领先”同时依赖前两阶段的排除证据；5 个 CE/BLER floor 候选没有被强行外推到 1%。结论仅适用于当前信道、资源、接收机、功率和独立候选随机流口径，不证明这些整数集合在其他 PDP、速度、接收机或 AL 下仍最优。本补充尚待研究者确认，未更新 `KNOWLEDGE.md`/`GOALS.md`。

复现命令为：

```powershell
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan031_strict_sidon_prescan.py --phase prescan
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan031_strict_sidon_prescan.py --phase extension
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan031_strict_sidon_prescan.py --phase coarse-confirmation
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan031_strict_sidon_prescan.py --phase fine
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan031_strict_sidon_prescan.py --phase fine-supplement
& D:\venvs\cdd-s102\Scripts\python.exe tools\analyze_plan031_strict_sidon_fine.py --bootstrap-repeats 10000
```

## 10. C300 4Tx/2Rx、2-symbol 补充：AL1/AL2 阶段性结果（2026-09-14）

### 10.1 范围、口径与完整性

本节只冻结已经完整结束的 AL1/AL2；AL4 尚在正式运行，不读取其中间曲线，也不据此排序。相对第 7、8 节的 1Rx C300 配置，本补充保持 4Tx、TDL-C 300 ns、2-symbol CORESET、单位总发射功率、每接收分支噪声定义、候选 precoder 和接收机协方差口径不变，只把 `n_rx` 改为 2；两个接收分支独立生成并在检测器中联合使用。estimated/ideal CSI 各自使用独立候选随机流，故跨模式和跨 Rx 比较采用独立样本 bootstrap，而不是逐 trial 配对。

目标 SNR 沿用 Jeffreys 平滑、按 trials 加权的单调递减 isotonic 曲线和局部 `log10(BLER)` 插值；只有目标被相邻且间距不超过 0.25 dB 的原始点双侧夹住时才报告。95% 区间来自 2,000 次 Bernoulli bootstrap。

| AL / CSI | candidates | points | trials | errors | adjacent BLER increases | 50k 且 errors<200 |
|---|---:|---:|---:|---:|---:|---:|
| AL1 / estimated | 9 | 144 | 2,568,800 | 92,429 | 3 | 9 |
| AL1 / ideal | 8 | 127 | 2,686,150 | 81,759 | 1 | 21 |
| AL2 / estimated | 10 | 158 | 2,775,950 | 106,785 | 2 | 13 |
| AL2 / ideal | 9 | 142 | 2,687,150 | 96,077 | 0 | 19 |
| 合计 | 36 | 571 | 10,718,050 | 377,050 | 6 | 62 |

逐 trial error flags 的长度与求和均和 CSV 一致，CE NMSE 线性均值复算最大误差为 0 dB；ideal-CSI 的 CE NMSE 全为 0。72 个 10%/1% 目标中 67 个合格。以下 5 个 1% 目标因正式网格在目标附近留有大于 0.25 dB 的间隙而不外推：AL1 ideal 的 `AP_RMS_T1`、`GEO_T1_CTRL_FOLD`、`PRG_DFT4_6REG`，AL2 estimated 的 `GEO_T1_CTRL_FOLD`，以及 AL2 ideal 的 `SMALL_CDD_IDEAL_CSI`。

### 10.2 BLER 与目标 SNR

表中数值为所需 SNR（dB）；`—` 表示没有合格原始 bracket。

| AL | CSI | candidate | BLER=10% | BLER=1% |
|---|---|---|---:|---:|
| 1 | estimated | B0 QC | 3.525 | 6.560 |
| 1 | estimated | S0 Sidon | 3.504 | 5.904 |
| 1 | estimated | AP RMS T1 | 3.817 | 7.087 |
| 1 | estimated | AP Teps T1 | 3.464 | 5.883 |
| 1 | estimated | GEO balanced/pair | 3.821 | 6.489 |
| 1 | estimated | GEO fold | 3.752 | 6.392 |
| 1 | estimated | PRG DFT4 6REG | 4.372 | 8.852 |
| 1 | estimated | small CDD matched | 4.465 | 8.805 |
| 1 | estimated | small CDD transparent | 4.383 | 8.875 |
| 1 | ideal | B0 QC | 1.944 | 4.658 |
| 1 | ideal | S0 Sidon | 1.311 | 3.483 |
| 1 | ideal | AP RMS T1 | 2.425 | — |
| 1 | ideal | AP Teps T1 | 1.324 | 3.416 |
| 1 | ideal | GEO balanced/pair | 1.262 | 3.354 |
| 1 | ideal | GEO fold | 1.269 | — |
| 1 | ideal | PRG DFT4 6REG | 3.413 | — |
| 1 | ideal | small CDD ideal | 3.299 | 7.451 |
| 2 | estimated | B0 QC | -1.255 | 0.861 |
| 2 | estimated | S0 Sidon | -1.172 | 0.892 |
| 2 | estimated | AP RMS T1 | -1.047 | 1.559 |
| 2 | estimated | AP Teps T1 | -1.087 | 0.865 |
| 2 | estimated | GEO balanced | -0.827 | 1.172 |
| 2 | estimated | GEO pair | -0.823 | 1.347 |
| 2 | estimated | GEO fold | -0.902 | — |
| 2 | estimated | PRG DFT4 6REG | -0.772 | 2.361 |
| 2 | estimated | small CDD matched | -0.132 | 3.452 |
| 2 | estimated | small CDD transparent | -0.173 | 3.531 |
| 2 | ideal | B0 QC | -3.163 | -1.151 |
| 2 | ideal | S0 Sidon | -3.241 | -1.344 |
| 2 | ideal | AP RMS T1 | -2.710 | -0.375 |
| 2 | ideal | AP Teps T1 | -3.349 | -1.499 |
| 2 | ideal | GEO balanced | -3.376 | -1.533 |
| 2 | ideal | GEO pair | -3.380 | -1.502 |
| 2 | ideal | GEO fold | -3.359 | -1.567 |
| 2 | ideal | PRG DFT4 6REG | -2.262 | 0.598 |
| 2 | ideal | small CDD ideal | -1.466 | — |

阶段性事实与推断：

- AL1 estimated 下，`AP_TEPS_T1` 的 10%/1% 为 3.464 `[3.416,3.520]` / 5.883 `[5.817,5.945]` dB；相对 B0 的增益分别为 0.061 `[-0.027,0.122]` / 0.677 `[0.526,0.767]` dB。因此 10% 尚不可分辨，1% 明确优于 B0。`S0_SIDON` 与它很接近。
- AL1 ideal 下，`GEO_BAL_PAIR` 的 10%/1% 为 1.262 `[1.204,1.311]` / 3.354 `[3.298,3.402]` dB，相对 B0 改善 0.682 `[0.598,0.785]` / 1.304 `[1.175,1.436]` dB。
- AL2 estimated 下，B0 的 10%/1% 为 -1.255 `[-1.304,-1.200]` / 0.861 `[0.772,0.924]` dB；1% 与 `AP_TEPS_T1` 的 0.865 `[0.789,0.923]` dB不可分辨，但相对 PRG 分别改善 0.483 `[0.381,0.588]` / 1.500 `[1.379,1.621]` dB。
- AL2 ideal 下，GEO 家族最低；10% 点估计最低的是 pair，1% 最低的是 fold。`GEO_FOLD` 相对 B0 在 10%/1% 改善 0.196 `[0.134,0.252]` / 0.416 `[0.288,0.519]` dB。

### 10.3 CE NMSE、估计代价与 1Rx→2Rx 收益

定义 `ΔCE = SNR_estimated,2Rx - SNR_ideal,2Rx`。全部可比较目标的 `ΔCE` 均为正：10% 处范围为 0.959–2.559 dB，1% 处为 1.354–3.135 dB。代表值如下：

| AL | candidate | ΔCE@10% (95% CI) | ΔCE@1% (95% CI) |
|---|---|---:|---:|
| 1 | B0 QC | 1.581 [1.480,1.659] | 1.902 [1.736,2.042] |
| 1 | AP Teps T1 | 2.140 [2.069,2.231] | 2.467 [2.376,2.544] |
| 1 | GEO balanced/pair | 2.559 [2.464,2.655] | 3.135 [3.040,3.250] |
| 1 | PRG DFT4 6REG | 0.959 [0.792,1.149] | — |
| 2 | B0 QC | 1.909 [1.846,1.982] | 2.012 [1.897,2.130] |
| 2 | AP Teps T1 | 2.262 [2.192,2.330] | 2.364 [2.278,2.464] |
| 2 | GEO pair | 2.557 [2.499,2.629] | 2.849 [2.752,2.962] |
| 2 | PRG DFT4 6REG | 1.490 [1.377,1.598] | 1.763 [1.619,1.923] |

定义 `ΔRx = SNR_1Rx - SNR_2Rx`，正值表示 2Rx 所需 SNR 更低。全部 65 个具备两侧合格 bracket 的比较均为正。estimated CSI 下，AL1 的 10%/1% `ΔRx` 范围为 4.978–6.487 / 6.241–8.117 dB，AL2 为 3.890–5.310 / 4.651–6.418 dB；ideal CSI 下，AL1 为 4.862–6.678 / 5.744–8.095 dB，AL2 为 4.112–5.555 / 4.695–6.193 dB。B0 的代表值为：

| AL | CSI | ΔRx@10% (95% CI) | ΔRx@1% (95% CI) |
|---|---|---:|---:|
| 1 | estimated | 5.449 [5.330,5.566] | 7.095 [6.958,7.277] |
| 1 | ideal | 5.631 [5.447,5.783] | 7.345 [7.131,7.490] |
| 2 | estimated | 4.160 [4.071,4.256] | 5.169 [4.986,5.373] |
| 2 | ideal | 4.429 [4.335,4.505] | 5.507 [5.369,5.645] |

上述 `ΔRx` 大于单纯 3 dB 阵列增益并不矛盾：这里改变的是随机信道的接收分支数，差值同时包含 MRC/联合检测的阵列增益、独立衰落分集以及 BLER 曲线尾部变化。它不是“2Rx 使单分支 CE NMSE 改善了同样 dB”的证据；完整 NMSE 数值在 `formal_points.csv`，端到端结论仍以 BLER 为准。

### 10.4 证据、复现与未完成项

阶段性分析目录为 `outputs/experiment031_pdcch_cdd/20260911_c300_4tx_2rx_2sym/analysis/partial_al1_al2/`。其中 `formal_points.csv`、`target_snr.csv`、`target_gains.csv`、`delta_ce.csv`、`delta_rx.csv`、`diagnostics.csv` 和 `analysis_metadata.json` 分别保存逐点汇总、完整目标及区间、候选相对增益、估计代价、1Rx/2Rx 差值、审计结果和配置/源码哈希。estimated AL1/ideal AL1/estimated AL2/ideal AL2 的展开配置 SHA-256 分别为 `8184d23f...9fe7`、`f9d98d79...4253`、`f9a72bd4...34a2`、`5e7ce5ad...d45f`；完整值见 metadata。分析时 Git HEAD 为 `83b25bc2ea77324ae98f7dc67136fc1e93b1f32c`。

复现阶段性分析：

```powershell
python tools/analyze_plan031_c300_2rx_partial.py --bootstrap-repeats 2000
```

限制：AL4 尚未完成，故本节不是 2Rx 补充的最终验收；5 个无合格 bracket 的 1% 目标需要后续边界补点才能报告。当前结论只适用于本配置，不外推到其他 PDP、速度、接收机、功率口径或 AL。阶段性结果尚待研究者确认，未更新 `KNOWLEDGE.md`/`GOALS.md`。

## 11. CDD911 与 transparent CDD AL1 BLER 补充（2026-09-14）

### 11.1 配置、数据复用与统计口径

本补充只覆盖 TDL-C 300 ns、4Tx、2-symbol CORESET、AL1。新增 `CDD911` 的物理人工时延为 `[0,0,911,911] ns`；AL1 的 $K=36$、$q_K=925.925926$ ns，因此配置坐标为 `[0,0,0.98388,0.98388]`，phase denominator 为 36。两个重复时延是预定波形的一部分；配置以 `allow_duplicate_delays: true` 显式放宽默认互异校验。预编码保持单位总发射功率。

系统其余字段冻结为 48-RB CORESET、每 symbol 实占 3 RB、18 个 DMRS RE、54 个 data RE、$E=108$、30 kHz、FFT 4096、CP 288、non-interleaved、first CCE 0、L=6 REG、A=41、CRC24C/RNTI `0xFFFF`、QPSK、Polar list size 8、TDL-C 300 ns、3 km/h、4 GHz 和 20 sinusoids。estimated-CSI 使用 non-transparent `matched_effective` 二维时频 LMMSE；ideal-CSI 直接使用真实 data-RE 等效信道。2Rx 使用独立同分布分支和 coherent MRC，每个 Rx 的噪声方差不随 Rx 数缩放。

只新运行 `CDD911` 的 1Rx/2Rx × estimated/ideal 四组。1Rx Sidon 和 precoder cycling 分别复用第 7、8 节的正式数据，不重跑，也不纳入 2Rx Sidon/cycling。CDD911 不做预扫，直接使用对应历史 1Rx/2Rx 配置中 Sidon 与 cycling 正式 SNR 点的并集；固定网格详见四份 `configs/pdcch_result031_cdd911_*_formal.yaml`。每点至少 10,000 trials、累计 200 errors 后停止、上限 50,000；目标 SNR 仅在相邻原始 bracket 不宽于 0.25 dB 且至少 95% 的 4,000 次独立 Bernoulli bootstrap 仍有 bracket 时报告。

随后新增 transparent estimated-CSI 补充：`CDD911_transparent` 复用 `[0,0,911,911] ns` 发射波形但接收机改用不知道 CDD 的 `physical_fullband` 协方差；`CDD130_transparent` 使用 `[0,0,130,130] ns`、坐标 `[0,0,0.1404,0.1404]` 和同一接收机。两者均运行 1Rx/2Rx，共 44 个固定 SNR 点；不新增 ideal-CSI 曲线。

第 11.2 节还复用第 10 节 AL1 的 `C300_S0_SIDON` 2Rx 正式数据：estimated-CSI 为 `matched_effective` non-transparent 接收，ideal-CSI 使用真实 data-RE 信道。原始 CSV 分别为 `outputs/experiment031_pdcch_cdd/20260911_c300_4tx_2rx_2sym/formal/{estimated,ideal}/al1/C300_S0_SIDON/bler_points.csv`，不重跑、不与本节 trial 合并。

2026-09-15 再新增 Sidon transparent 1Rx estimated-CSI 和 CDD130 1Rx/2Rx ideal-CSI。前者复用 Sidon `[0,1,3,7]` 发射波形并使用 `physical_fullband` 接收；后者复用 `[0,0,130,130] ns` 波形并直接使用真实 data-RE 信道。三条曲线共 29 个固定点，使用独立 seed/namespace。

同日补入 2Rx `DFTcodebook`/`LTEcodebook` 的 estimated/ideal-CSI 曲线。AL1 只有一个 6-REG bundle，`DFTcodebook` 固定使用 4 维 DFT 索引 0，即 `[1,1,1,1]/2`；`LTEcodebook` 固定使用数值等价于 DFT 索引 2 的 `[1,-1,1,-1]/2`。estimated-CSI 均采用 transparent `physical_prg` LMMSE，ideal-CSI 使用真实 data-RE 信道。目标网格均为 `5,6,7,7.5,8,8.5,9,10 dB`；DFT 已有的 estimated `5,8.5,9,10` 和 ideal `8,8.5,9` 从 2026-09-11 正式数据复用，本次只运行其余 9 点及 LTE 的全部 16 点。

### 11.2 BLER 原始结果与目标

图对应的数据完整保存在 `outputs/experiment031_pdcch_cdd/20260914_c300_al1_cdd911/analysis/formal_points.csv`。该文件可独立重建十二曲线 estimated-CSI 图和九曲线 ideal-CSI 图。estimated 图例对每条 estimated-CSI 曲线统一显式标注 `transparent` 或 `non-transparent`；ideal 图例标注 `ideal-CSI`，因为 ideal 接收不存在信道估计透明性差别。绘图仅使用正式仿真点；跨采样空档的连线只作视觉引导，不用于目标判定。

下表中方括号对正式目标表示 bootstrap 95% CI；`—` 表示固定网格没有不宽于 0.25 dB 的合格 bracket，随后列出的范围只是跨越目标的相邻原始采样点，禁止在其间形成正式插值。

| CSI | curve | SNR@10% / dB | SNR@1% / dB |
|---|---|---:|---:|
| estimated | CDD911, 1Rx | 9.236 `[8.975,9.324]` | 14.316 `[14.169,14.443]` |
| estimated | CDD911, 2Rx | 3.515 `[3.433,3.591]` | 6.831 `[6.690,7.035]` |
| estimated | CDD911 transparent, 1Rx | —；`<10.00` | —；raw `[14.00,16.00]` |
| estimated | CDD911 transparent, 2Rx | —；raw `[4.00,5.00]` | —；raw `[7.00,7.50]` |
| estimated | CDD130 transparent, 1Rx | —；raw `[10.00,11.00]` | —；raw `[16.50,17.00]` |
| estimated | CDD130 transparent, 2Rx | —；raw `[4.00,5.00]` | —；raw `[8.50,9.00]` |
| estimated | Sidon, 1Rx | 8.504 `[8.419,8.618]` | 12.789 `[12.657,12.921]` |
| estimated | Sidon non-transparent, 2Rx | 3.504 `[3.442,3.557]` | 5.904 `[5.851,5.964]` |
| estimated | Sidon transparent, 1Rx | —；`>13.00` | —；`>13.00` |
| estimated | cycling, 1Rx | 10.846 `[10.675,11.033]` | 16.924 `[16.845,17.035]` |
| estimated | DFTcodebook transparent, 2Rx | —；`<5.00` | —；raw `[8.50,9.00]` |
| estimated | LTEcodebook transparent, 2Rx | —；`<5.00` | —；raw `[8.50,9.00]` |
| ideal | CDD911, 1Rx | 7.952 `[7.862,8.042]` | 12.789 `[12.641,12.903]` |
| ideal | CDD911, 2Rx | 2.251 `[2.182,2.325]` | 5.328 `[5.189,5.482]` |
| ideal | Sidon, 1Rx | 6.235 `[6.165,6.309]` | 9.234 `[9.153,9.519]` |
| ideal | Sidon, 2Rx | 1.311 `[1.270,1.342]` | 3.483 `[3.396,3.565]` |
| ideal | cycling, 1Rx | 9.798 `[9.677,9.898]` | 15.508 `[15.349,15.968]` |
| ideal | CDD130, 1Rx | —；`<10.00` | —；raw `[14.00,16.00]` |
| ideal | CDD130, 2Rx | —；`<4.00` | —；raw `[7.00,7.50]` |
| ideal | DFTcodebook, 2Rx | —；`<5.00` | —；raw `[7.50,8.00]` |
| ideal | LTEcodebook, 2Rx | —；`<5.00` | —；raw `[7.50,8.00]` |

CDD911 定向补点后八个目标均合格。estimated-CSI 的 1Rx→2Rx 收益在 10%/1% 为 5.722 `[5.432,5.849]` / 7.485 `[7.214,7.693]` dB；ideal-CSI 为 5.701 `[5.583,5.823]` / 7.461 `[7.256,7.651]` dB。estimated 相对 ideal 的代价在 1Rx 为 1.285 `[0.983,1.415]` / 1.527 `[1.339,1.724]` dB，在 2Rx 为 1.264 `[1.148,1.370]` / 1.503 `[1.284,1.785]` dB。四条 transparent 曲线均达到原始 BLER 不高于 1%，但其粗网格没有合格目标，因此不报告 matched-transparent、911-130 或 1Rx-2Rx 的正式门限差。

#### CDD911 从 1Rx 到 2Rx 的增益解释

上述 5.722/7.485 dB 已由新增 0.25 dB 局部点和 bootstrap 正式闭合，不再是跨空档示意插值。

该增益在当前模型下是合理的。这里的 SNR 定义为单个 Rx 分支上的 $E_s/N_0$；2Rx 采用两个独立、无空间相关的衰落分支并进行相干最大比合并（MRC），合并后瞬时 SNR 为

$$
\gamma_{\mathrm{MRC}}=\gamma_1+\gamma_2.
$$

因此，2Rx 除了平均合并功率翻倍所对应的约 3 dB 增益，还把接收分集阶数由一阶提高到二阶，显著降低深衰落事件概率。在 TDL 衰落的低 BLER 区，分集改变了曲线斜率，固定 BLER 处的水平间距可以明显超过 3 dB；“只增加 3 dB”更适用于固定信道、AWGN 或只比较平均合并后 SNR 的情形。

该额外收益不应归因于信道估计器改善。ideal-CSI 曲线仍显示相近的 1Rx/2Rx 分离趋势，且两种 Rx 的 estimated-ideal 差值均已正式闭合。实现中各 Rx 分支独立估计后再做 MRC，且每个分支的噪声方差不随 $N_{\mathrm{Rx}}$ 缩放。因此，当前观察到的主要是接收分集收益。需要注意，本轮假设两个 Rx 分支完全不相关；真实天线存在空间相关性时，实际增益会减小。

### 11.3 同 SNR 原始比较与 CE 诊断

以下代表点直接来自正式原始点；BLER 后的方括号为逐点 Wilson 95% 区间。它们不依赖跨网格插值。

| CSI | SNR / dB | CDD911 1Rx BLER | comparison 1Rx BLER | comparison |
|---|---:|---:|---:|---|
| estimated | 9.00 | 0.1028 `[0.0970,0.1089]` | 0.0784 `[0.0733,0.0838]` | Sidon |
| estimated | 14.00 | 0.01264 `[0.01102,0.01450]` | 0.005714 `[0.004977,0.006560]` | Sidon |
| estimated | 10.75 | 0.0546 `[0.0503,0.0592]` | 0.1025 `[0.0967,0.1086]` | cycling |
| estimated | 16.50 | 0.00320 `[0.002742,0.003735]` | 0.01163 `[0.01013,0.01334]` | cycling |
| ideal | 6.25 | 0.1864 `[0.1789,0.1942]` | 0.0990 `[0.0933,0.1050]` | Sidon |
| ideal | 10.75 | 0.0284 `[0.0253,0.0318]` | 0.002680 `[0.002263,0.003173]` | Sidon |
| ideal | 9.75 | 0.0448 `[0.0409,0.0490]` | 0.1020 `[0.0962,0.1081]` | cycling |
| ideal | 15.50 | 0.00160 `[0.001286,0.001991]` | 0.0100 `[0.008715,0.01147]` | cycling |

事实：所有共同 SNR 点的方向一致。estimated-CSI 下，CDD911/Sidon 的 BLER 比为 1.16–2.80，CDD911/cycling 为 0.218–0.566；ideal-CSI 下相应范围为 1.58–10.60 和 0.160–0.500。即本场景的 1Rx 排序稳定为 Sidon 优于 CDD911、CDD911 优于 precoder cycling。ideal-CSI 下仍保持同一方向，说明该排序不只由信道估计误差造成；这是本冻结场景的事实，不是对其他 PDP、AL 或时延的普适排序。

CDD911 的 4Tx pilot 矩阵秩为 2、四列条件数为无穷，这是 `[0,0,911,911] ns` 只有两个不同时延的直接结果。matched 接收机估计的是标量等效信道而不是分别恢复四个底层 Tx 分支；解析 zero-noise CE floor 为 -72.196 dB，未显示由 rank 2 引起的有效信道估计 floor。代表点的 CE NMSE 为：1Rx 在 9/14 dB 时 -13.298/-17.329 dB，2Rx 在 3.5/6.5 dB 时 -10.116/-12.412 dB。

同 SNR 下 CDD911 的 CE NMSE 还低于 Sidon，例如 estimated 1Rx、9 dB 时为 -13.298 dB，而 Sidon 为 -10.760 dB；但 CDD911 的 BLER 更高。这是“较低 CE NMSE 不保证更低 BLER”的直接证据，波形的频率分集、有限码长译码事件及其联合分布仍影响排序。相反，CDD911 相对 cycling 同时具有较低 BLER和略低 CE NMSE；当前结果不能把两者差值只归因于其中一个因素。

transparent 新曲线均通过 1% 原始点验收：CDD911 1Rx/2Rx 的最低点分别为 19.5 dB 时 0.00080 `[0.000588,0.001089]`（40/50,000）和 10 dB 时 0.00084 `[0.000622,0.001135]`（42/50,000）；CDD130 1Rx/2Rx 分别为 0.00260 `[0.002190,0.003086]`（130/50,000）和 0.00356 `[0.003075,0.004122]`（178/50,000）。共同 SNR 点上，CDD911 transparent/matched 的 BLER 比在 1Rx 为 1.10–1.85、2Rx 为 1.40–1.69，透明协方差失配一致地增加 BLER。CDD911/CDD130 transparent 的 BLER 比在 1Rx 为 0.31–0.64、2Rx 为 0.24–0.95；当前冻结网格上 911 ns 波形优于 130 ns，但该排序不能外推到其他 PDP 或资源配置。CDD911/CDD130 transparent 的解析 zero-noise CE floor 分别为 -33.398/-73.622 dB；更低的 CDD130 CE floor 并未对应更低 BLER，再次说明 CE NMSE 或 floor 不能单独决定链路排序。

Sidon transparent 是明确负结果：7 个点各运行 10,000 trials，BLER 从 8 dB 的 0.8692 仅下降到 13 dB 的 0.7514 `[0.7428,0.7598]`，远未达到 10% 或 1%；解析 zero-noise CE floor 为 +6.054 dB。与同波形 non-transparent Sidon 的共同点相比，transparent/non-transparent BLER 比由 8 dB 的 6.47 增至 13 dB 的 87.54。这表明当前大时延 Sidon 波形与只使用物理 PDP 的透明估计器严重失配；它不否定该波形在 matched 或 ideal 接收下的性能。

CDD130 ideal 的 1Rx/2Rx 均通过 1% 原始点验收：最低点分别为 19.5 dB 的 0.00112 `[0.000863,0.001454]`（56/50,000）和 10 dB 的 0.00124 `[0.000967,0.001589]`（62/50,000）。10% 目标已落在采样下界之外；1% raw bracket 分别为 `[14.00,16.00]` 和 `[7.00,7.50]`，均宽于 0.25 dB，因此不报告正式门限或 estimated-ideal 差值。

两类 codebook 的 2Rx 曲线均通过 1% 原始点验收。10 dB 时，estimated-CSI 的 DFT/LTE BLER 分别为 0.004077 `[0.003551,0.004682]`（200/49,050）和 0.004362 `[0.003799,0.005008]`（200/45,850）；ideal-CSI 分别为 0.00154 `[0.001232,0.001924]`（77/50,000）和 0.00130 `[0.001020,0.001656]`（65/50,000）。共同 SNR 点的 DFT/LTE BLER 比在 estimated 下为 0.810–1.088、ideal 下为 0.887–1.185，逐点 Wilson 区间均重叠；estimated CE NMSE 的两曲线最大绝对差为 0.0437 dB。事实是当前样本未显示可统计区分的差异；在独立同分布 Tx 分支下，两个单位范数固定向量具有相同边缘分布，这与观察一致，但不表示逐 trial 信道或判决相同。两组 1% raw bracket 宽度均为 0.5 dB，故不报告正式 codebook 门限差。

### 11.4 样本审计、证据与限制

| CSI | Rx | points | trials | errors | raw BLER increases | 50k 且 errors<200 |
|---|---:|---:|---:|---:|---:|---:|
| estimated CDD911 matched | 1 | 26 | 512,000 | 14,330 | 0 | 5 |
| estimated CDD911 matched | 2 | 27 | 636,400 | 11,520 | 0 | 8 |
| ideal CDD911 | 1 | 33 | 694,400 | 22,998 | 2 | 7 |
| ideal CDD911 | 2 | 32 | 686,750 | 25,422 | 2 | 8 |
| estimated CDD911 transparent | 1 | 12 | 424,600 | 2,995 | 0 | 6 |
| estimated CDD911 transparent | 2 | 10 | 298,100 | 3,144 | 0 | 4 |
| estimated CDD130 transparent | 1 | 12 | 307,800 | 4,820 | 0 | 2 |
| estimated CDD130 transparent | 2 | 10 | 193,850 | 3,979 | 0 | 1 |
| 本节新增 CDD 数据合计 | — | 162 | 3,753,900 | 89,208 | 4 | 41 |
| estimated Sidon transparent | 1 | 7 | 70,000 | 55,851 | 0 | 0 |
| ideal CDD130 | 1 | 12 | 406,600 | 3,290 | 0 | 6 |
| ideal CDD130 | 2 | 10 | 294,200 | 2,519 | 0 | 3 |
| 2026-09-15 补充合计 | — | 29 | 770,800 | 61,660 | 0 | 9 |
| estimated DFTcodebook（含 4 个历史点） | 2 | 8 | 139,000 | 2,569 | 0 | 0 |
| estimated LTEcodebook | 2 | 8 | 135,400 | 2,577 | 0 | 0 |
| ideal DFTcodebook（含 3 个历史点） | 2 | 8 | 223,700 | 1,730 | 0 | 2 |
| ideal LTEcodebook | 2 | 8 | 220,500 | 1,729 | 0 | 2 |
| codebook 新运行合计 | — | 25 | 504,200 | 6,648 | 0 | 4 |

十一份配置的全部 SNR 点均完成。216 个新增曲线点的 trial flags 长度、错误和、Wilson 区间及 CE 线性均值复算误差均为 0；ideal CE NMSE 逐点为 0。54 个点达到 50,000 trials 但不足 200 errors，符合停止规则。分析器已完整执行 4,000 次 bootstrap；TensorFlow 仅报告既有 deprecated API warning。

分析目录为 `outputs/experiment031_pdcch_cdd/20260914_c300_al1_cdd911/analysis/`。`formal_points.csv` 含 216 个新增曲线点和 94 个历史复用点，共 310 点；其余目标、差值、同 SNR 对照、审计、来源、样式及十一份配置 SHA-256 分别保存在同目录的 `target_snr.csv`、`target_differences.csv`、`same_snr_comparisons.csv`、`same_snr_cdd911_comparisons.csv`、`diagnostics.csv`、`source_receipt.json`、`curve_styles.json` 和 `analysis_metadata.json`。约 13 cm 宽的本地预览保存在 `analysis/previews_13cm/`，待研究者人工检查图例与曲线可辨识性；本分析未读取图片。分析时 Git HEAD 为 `83b25bc2ea77324ae98f7dc67136fc1e93b1f32c`，本补充仍是未提交工作区增量。

复现命令为：

```powershell
python tools/run_plan031_candidates.py --config configs/pdcch_result031_cdd911_1rx_estimated_formal.yaml --config configs/pdcch_result031_cdd911_2rx_estimated_formal.yaml --config configs/pdcch_result031_cdd911_1rx_ideal_formal.yaml --config configs/pdcch_result031_cdd911_2rx_ideal_formal.yaml --stage run --max-workers 4
python tools/run_plan031_candidates.py --config configs/pdcch_result031_cdd_transparent_1rx_formal.yaml --config configs/pdcch_result031_cdd_transparent_2rx_formal.yaml --stage run --max-workers 4
python tools/run_plan031_candidates.py --config configs/pdcch_result031_sidon_transparent_1rx_formal.yaml --config configs/pdcch_result031_cdd130_1rx_ideal_formal.yaml --config configs/pdcch_result031_cdd130_2rx_ideal_formal.yaml --stage run --max-workers 3
python tools/run_plan031_candidates.py --config configs/pdcch_result031_codebooks_2rx_estimated_formal.yaml --config configs/pdcch_result031_codebooks_2rx_ideal_formal.yaml --stage run --max-workers 4
python tools/analyze_plan031_cdd911.py --bootstrap-repeats 4000
```

限制：CDD transparent、CDD130 ideal 与两类 codebook 的粗网格不能形成 0.25 dB 正式目标；Sidon transparent 在 13 dB 内没有接近 10%/1%，当前只能记录严重失配负结果，不能给出门限。codebook 比较只覆盖 AL1 单 bundle 和独立同分布 Tx/Rx 分支，不能外推到多 bundle cycling、相关天线或其他 codebook。结论只适用于当前 TDL-C 300 ns、4Tx、AL1、2-symbol 与指定 matched/transparent/ideal 接收口径。结果尚待研究者确认，未更新 `KNOWLEDGE.md`/`GOALS.md`，也未创建 Git checkpoint。
