# plan-034-RSRP统计-PDCCH 4T4R：无噪声块平均接收功率 CDF 与固定 precoder 的 4Tx/4Rx PDCCH BLER

> 状态：原始计划和第 10 节逐 Rx RSRP 补充实验已执行，结果待研究者确认。原始计划只新增机理诊断和一条固定 precoder BLER 曲线；补充实验只扩展无噪声 RSRP 诊断，不新增 BLER 曲线、不搜索新预编码、不追加其他 AL/Rx/信道场景，也不以 RSRP-CDF 替代链路 BLER。

## 1. 目的与研究问题

本轮在 `research/plan-031-PDCCH-CDD时延-BLER.md` 已定义的 C300、PDCCH AL1 单 bundle 场景上，把接收天线数固定为 4，并回答两个相互区分的问题：

1. 在同一批小尺度衰落信道实现下，频变预编码是否会降低无噪声真实等效信道 $\mathbf H\mathbf W$ 的块平均接收功率波动，使其经验累积分布函数（CDF）更集中，并改善 1% 和 10% 低尾分位数；
2. 在 4 GHz、30 kHz SCS、4Tx/4Rx、TDL-C 300 ns、3 km/h、2-symbol CORESET、AL1、40-bit payload 加 24-bit CRC 条件下，一个仅含单个 6-REG bundle 的固定 precoder，其 estimated-CSI PDCCH BLER 在 `[-3,-2,-1,0,1,2,3,4,5] dB` 九个采样点上是多少。

第一个问题是分集机理诊断：独立同分布 Tx 分支下，不同单位范数预编码向量通常不改变单 RE 的边缘功率分布，但频变 $\mathbf W$ 会改变跨频率相关性，因而可能改变块平均功率的方差和低尾。该判断与 `DESIGN.md` 第 3.2、3.6 节关于“单点边缘均值不变、块统计由跨频率相关性决定”的框架一致。

第二个问题只报告指定 SNR 点的原始 BLER 和不确定性，不根据本轮数据搜索 precoder，也不插值或外推 10%/1% BLER 门限。RSRP-CDF 是解释量，不是链路验收量；即使其更集中，也不能单独推出 BLER 一定改善。

## 2. 冻结系统条件

| 参数 | 冻结取值 |
|---|---|
| 载频 / SCS | 4 GHz / 30 kHz |
| FFT / CP | 4096 / 288 samples |
| CORESET | 48 RB，2 个相邻 OFDM symbols，起始 symbol 0 |
| PDCCH candidate | AL1，first CCE 0，non-interleaved，`L=6 REG`；每 symbol 占 3 个连续 RB、36 个子载波；两个 symbols 共 72 个候选 RE |
| RE 计数 | 18 DMRS RE、54 data RE、QPSK，编码后比特数 $E=108$ |
| DCI | payload $A=40$ bit；CRC24C 共 24 bit；RNTI `0xFFFF`；Polar list size 8，`hybSCL`，CPU-only |
| 天线 / 层数 | 4Tx / 4Rx / 1 layer |
| 信道 | Sionna 1.0.2 TDL-C，RMS delay spread 300 ns，3 km/h，20 sinusoids；16 个 Tx--Rx 分支独立同分布、无空间相关 |
| 同步假设 | 完美定时和载波同步；不模拟 CFO、ICI、ISI、天线互耦、相关性或分支增益失配 |
| 信道归一化 | 只使用 TDL 模型的平均功率归一化；禁止逐 realization 归一化，以保留自然的小尺度功率起伏 |
| 接收机 | 每根 Rx 独立做二维时频 LMMSE；四根 Rx 在 data RE 上 coherent MRC；LLR 不额外加入 CE-error-aware 方差 |
| SNR 定义 | 单位总发射功率相对单根 Rx 分支噪声功率；每根 Rx 的噪声方差均为 $10^{-\mathrm{SNR}_{\rm dB}/10}$，不得因 4Rx 再除以 4 |

速度 3 km/h 在 4 GHz 下只用于两个相邻 CORESET symbols 的时间相关。不同 trial 是独立 TDL realization；同一 trial 内两个 symbols 来自同一连续 realization，不得分别生成。

## 3. 预编码、功率归一化与公平性

### 3.1 固定基线 `FIXED_DFT0`

AL1 只有一个 6-REG bundle。固定 precoder 在该 bundle 的两个 symbols、全部 data/DMRS RE 上均使用同一个 4 维 DFT index 0 向量：

$$
\mathbf W_q^{\rm fixed}
=\frac{1}{2}[1,1,1,1]^T,
\qquad q\in\mathcal Q.
$$

因此每个 RE 都满足 $\lVert\mathbf W_q^{\rm fixed}\rVert_2^2=1$。配置实现复用 `reg_bundle_dft_cycling`，但 `cycling_order=[0]` 且只有一个 bundle，result 中必须称为 fixed DFT0，不得称为多 bundle cycling。

### 3.2 频变诊断候选 `FREQ_SIDON_0137`

RSRP-CDF 的频变对照复用 plan-031 C300 AL1 的四路严格 Sidon 坐标

$$
\mathbf j=[0,1,3,7],
$$

并在每个 symbol 的 36 个局部子载波上构造

$$
W_{k,n}^{\rm freq}
=\frac{1}{2}\exp\!\left(-j2\pi k j_n/36\right),
\qquad k=0,\ldots,35.
$$

phase denominator 固定为 36；对应 delay-grid 单位为

$$
q_{36}=\frac{1}{36\times30\ \mathrm{kHz}}
=925.925926\ \mathrm{ns},
$$

物理人工时延为 `[0,925.925926,2777.777778,6481.481481] ns`。人工时延只进入数字相位，不叠加到真实传播 tap 或 CP。两个 CORESET symbols 在相同子载波上使用相同 $\mathbf W_k^{\rm freq}$。每行同样满足单位总发射功率。

该候选只参与 RSRP-CDF 机理诊断，不在 034 中新增 BLER 曲线。选择它是因为它是当前 C300 AL1 已定义、结构明确的频变预编码；本轮不利用 034 数据再挑选其他 delay 集。

### 3.3 配对与禁止项

RSRP-CDF 中两个预编码严格共享每个 absolute trial 的底层 $\mathbf H$，仅替换 $\mathbf W$。两者使用相同 4Tx/4Rx、同一 72 RE 集、相同 $P_{\rm tx}=1$ 和相同 TDL 参数。禁止：

- 对任一预编码单独归一化每个 trial 的 $\mathbf H\mathbf W$；
- 按观察到的 CDF 重新选择 delay、相位、RE 子集或 trial；
- 将四根 Rx 的合计功率当成改变后的发射功率；
- 把 RSRP 统计混入 AWGN、DMRS 信道估计、MRC 权重或译码结果。

BLER 只运行 `FIXED_DFT0`，不因中途曲线位置新增候选、SNR 点或改变接收机口径。

## 4. 无噪声真实 $\mathbf H\mathbf W$ 的 RSRP-CDF 诊断

### 4.1 统计量定义

本 plan 的“RSRP”是研究者指定的归一化块平均无噪声接收功率，不宣称等同于 3GPP UE RSRP 测量流程。令

$$
\mathcal Q=\{(\ell,k):\ell\in\{0,1\},\ k=0,\ldots,35\},
\qquad |\mathcal Q|=72,
$$

即 AL1 candidate 的全部两个 symbols 和 36 个占用子载波；$q=(\ell,k)$ 把 symbol 索引并入 RE 索引。对第 $t$ 个独立信道 realization，定义

$$
\overline P_t(\mathbf W)=
\frac{1}{N_r|\mathcal Q|P_{\rm tx}}
\sum_{q\in\mathcal Q}\sum_{r=0}^{N_r-1}
\left|
\sum_{n=0}^{N_t-1}H_{t,r,n,q}W_{q,n}
\right|^2,
$$

其中 $N_t=N_r=4$、$P_{\rm tx}=1$。保存线性值和

$$
X_t(\mathbf W)=10\log_{10}\overline P_t(\mathbf W)\quad[\mathrm{dB}],
$$

经验 CDF 为

$$
\widehat F_{\mathbf W}(x)=
\frac{1}{T}\sum_{t=1}^{T}
\mathbf 1\{X_t(\mathbf W)\le x\}.
$$

主统计使用 `T=100000` 个独立 realization、base seed `20260916`、随机流 namespace `plan034_rsrp_cdf_v1`。absolute trial 从 1 开始并可断点续跑；两个候选的每个 trial key 必须相同。保存逐 trial $\overline P_t$，不得只保存直方图。

经验分位数使用左连续逆 CDF

$$
Q_p(\mathbf W)=\inf\{x:\widehat F_{\mathbf W}(x)\ge p\},
$$

重点报告：

- $Q_{0.01}$：1% 低尾分位数，越高越好；
- $Q_{0.10}$：10% 低尾分位数，越高越好；
- $Q_{0.90}-Q_{0.10}$：中央 80% 宽度，越小表示越集中；
- 同时列出均值、标准差、$Q_{0.50}$ 和 $Q_{0.90}$，用于审计平均功率是否发生非预期偏移。

### 4.2 单 RE 边缘分布控制

预先固定 `q0=(symbol 0, local subcarrier 18)`，额外保存两个预编码在该单 RE 上按 4Rx 平均的无噪声功率。该控制只检查“边缘分布应基本不变、块平均分布可能改变”的机理：报告单 RE 的 $Q_{0.01}$、$Q_{0.10}$、$Q_{0.50}$ 及两条经验 CDF 的最大绝对差。不得在看到结果后改选更接近或更分离的 RE。

### 4.3 不确定性与判定

CDF 图直接由排序后的逐 trial 原始值绘制，不做 KDE、分布拟合或平滑。每条经验 CDF 同时给出 95% Dvoretzky--Kiefer--Wolfowitz 带；两个预编码共享同一 trial，因此分位数和宽度差使用 1000 次 paired nonparametric bootstrap，固定 bootstrap seed `20260916`，报告 percentile 95% 区间。

定义

$$
\Delta Q_p=Q_p(\mathbf W^{\rm freq})-Q_p(\mathbf W^{\rm fixed}),
$$

以及

$$
\Delta B_{80}=
\bigl[Q_{0.90}-Q_{0.10}\bigr]_{\rm freq}
-\bigl[Q_{0.90}-Q_{0.10}\bigr]_{\rm fixed}.
$$

判定冻结为：

- **支持“块平均功率更集中”**：$\Delta Q_{0.01}>0$、$\Delta Q_{0.10}>0$、$\Delta B_{80}<0$，且三者的 paired-bootstrap 95% 区间均不跨 0；
- **部分支持**：低尾和宽度方向不全一致，或至少一个区间跨 0；逐项报告，不合并成单一胜负；
- **不支持**：两个低尾均未改善，或中央 80% 宽度显著增大；
- **实现/公平性失败**：任一预编码存在行功率偏离 1、候选未共享信道、逐 realization 被归一化、存在噪声，或单 RE 边缘出现无法由有限样本解释的明显系统偏移。此时停止机理结论，先修复实现并从 trial 1 重跑。

单 RE 控制的“明显系统偏移”预注册为任一 $Q_{0.01}$、$Q_{0.10}$、$Q_{0.50}$ 差的绝对值超过 `0.10 dB`，或两条单 RE 经验 CDF 的最大绝对差超过 `0.01`。该门槛只用于触发实现审计，不把通过控制当成分集增益证据。

## 5. 固定 precoder 的 PDCCH BLER

### 5.1 冻结曲线

只运行一条曲线：

| curve ID | 发射端 | 接收端 |
|---|---|---|
| `C300_AL1_FIXED_DFT0_NR4_EST` | AL1 单 bundle，全部 data/DMRS RE 固定 `DFT4 index 0` | transparent `physical_prg` 二维时频 LMMSE；4Rx coherent MRC |

`physical_prg` 只使用 TDL-C 的物理时频协方差，不依赖频变 CDD。由于单 bundle 内 precoder 固定且 Tx 分支独立同分布，该接收口径与现有 fixed-codebook 基线一致。

### 5.2 SNR 点与动态 trial 预算

正式 SNR 点固定为

$$
[-3,-2,-1,0,1,2,3,4,5]\ \mathrm{dB}.
$$

不做 prescan，不追加中间点，不因结果偏离预期而移动 SNR 网格。为避免低 SNR 高 BLER 区无效堆样本，同时保证预计接近 1% 的高 SNR 点有足够错误块，按下表分三档运行。每档使用独立 YAML，但三份 YAML 除 SNR、trial 预算、输出子目录和随机流 namespace 外，物理、编码、接收机和 candidate 字段必须逐项一致，并生成递归 diff 回执。

| SNR / dB | 最少 trials/点 | 累计错误停止目标 | 最多 trials/点 | 预期用途 |
|---|---:|---:|---:|---|
| `-3,-2,-1,0` | 2,000 | 200 | 5,000 | 预计 BLER $>0.1$；达到最少 trial 且已有 200 errors 即停止 |
| `1,2` | 5,000 | 200 | 20,000 | 10% 到数个百分点的过渡区 |
| `3,4,5` | 10,000 | 200 | 50,000 | 预计接近 1%；优先保证低 BLER 统计量 |

batch size 固定为 50。每点在同时满足“最少 trials”和“至少 200 errors”后停止；若达到该档上限仍不足 200 errors，则保留原始点、Wilson 95% 区间和上限标记，不补伪计数、不后验增加上限。上述预算是正式规则，不根据中间 BLER 手工调整。

每个 SNR 使用由 base seed `20260916`、namespace、SNR 和 absolute trial 唯一确定的信道、payload、data AWGN 与 DMRS AWGN。三档的 namespace 分别为 `plan034_fixed_dft0_low_v1`、`plan034_fixed_dft0_mid_v1`、`plan034_fixed_dft0_high_v1`；不同 SNR 不要求 paired。保存逐 trial error flag 和 CE NMSE。

### 5.3 BLER 输出与解释

逐点报告：trials、TB errors、BLER、Wilson 95%区间、线性域平均后转 dB 的 data-RE CE NMSE、平均 MRC denominator 和实际停止原因。只画原始九点连线，连线只作视觉引导；不做 isotonic、PAVA、PCHIP、log-BLER 门限插值或外推。

用户关于“0 dB 前 BLER 大于 0.1、3--5 dB 接近 0.01”的判断只用于预分配 trial，不是验收门槛。result 必须逐项标注哪些采样点支持或不支持该预期，不因预期不成立而视为实现失败。

## 6. 实现范围、定向测试与 smoke

### 6.1 预计代码与配置修改

1. 在 `cdd_lls/` 增加可复用的无噪声等效信道块平均功率计算，输入显式为 `H`、`W`、RE 索引、$N_r$ 和 $P_{\rm tx}$，不把统计公式复制到绘图脚本；
2. 新增 `tools/run_plan034_rsrp_cdf.py`，复用现有 active-only Sionna TDL 生成、PDCCH 网格和 precoder 构造，支持 absolute-trial 重放、分批落盘和 resume；
3. 最小扩展 `pdcch-cdd-bler-v2` 校验，使冻结 C300 AL1 场景允许 `A=40` 和任意正整数 Rx，且现有 `A=41`、1Rx/2Rx 配置行为不变；把当前只在 `n_rx==2` 触发的接收分支 validation receipt 泛化为 `n_rx>1`，本轮必须实际记录 4Rx；
4. 新增三份固定 DFT0 BLER YAML 和一份 RSRP YAML；不复制链路主循环；
5. 新增 `tools/analyze_plan034_rsrp_pdcch_4t4r.py`，只读保存的 CSV/NPY 生成数值汇总、paired bootstrap、本地图和审计回执。

计划文件不授权重构无关 PDSCH runner、修改历史 result、更新 `KNOWLEDGE.md`/`GOALS.md` 或重跑 result-031。

### 6.2 只覆盖本轮增量的测试

不运行整个历史测试集。新增测试仅覆盖：

1. 小型手算复数数组的 $\overline P_t$ 与公式逐元素一致，$\lVert\mathbf W_q\rVert_2^2=1$ 时 $P_{\rm tx}=1$，且函数不做逐 trial 归一化；
2. AL1 网格确有一个 6-REG bundle、72 个候选 RE，`FIXED_DFT0` 在全部 RE 完全相同，`FREQ_SIDON_0137` 的 phase denominator、坐标和行功率正确；
3. `A=40`、4Tx/4Rx 配置可通过校验；一个 4Rx estimated-CSI batch 的信道、CE、MRC denominator、effective noise 和 error flag 维度正确，手算 MRC 与实现一致；
4. 相同 absolute trial 重跑得到逐元素相同的 RSRP；两个预编码共享相同底层信道 key。

定向回归只运行现有 PDCCH 无噪声编解码、C300 2-symbol 网格、fixed DFT codebook、正整数 `n_rx`、多 Rx MRC 和 active-only TDL 测试。除非上述定向测试暴露公共模块回归，不运行全仓测试。

### 6.3 smoke 与停止

正式运行前只做两个轻量 smoke：

- RSRP：100 个 paired trials，检查所有功率有限且为正、两候选 trial key 相同、逐 trial 文件可 resume、单 RE 与块平均字段均落盘；
- BLER：`SNR=2 dB`、20 trials，检查 A=40+CRC24C 无噪声编解码、4Rx 信道形状、estimated-CSI CE、MRC、CRC/error flag、噪声缩放和落盘字段。

smoke 不用于估计 CDF、分位数或正式 BLER。任一功率/维度/seed/无噪声译码检查失败时停止，不启动正式 trial；修复后 smoke 从头运行。

## 7. 执行顺序

1. 创建并校验冻结 YAML、candidate/资源/功率 manifest 和三档 BLER config-diff 回执；
2. 实现本轮最小代码增量并运行第 6.2 节定向测试；
3. 依次执行 RSRP 与 BLER smoke，记录 CPU、Python/NumPy/SciPy/TensorFlow/Sionna 版本、Git HEAD、工作区变更标识、墙钟和进程峰值 RSS；
4. 运行 100,000 个 paired RSRP trials；完成后先审计 trial 连续性、candidate 配对、功率和逐 trial 数量，再分析 CDF；
5. 按低、中、高三档运行九个 BLER 点。已有点续跑必须从已保存的下一个 absolute trial 继续，禁止从 trial 1 重跑后与旧数据相加；
6. 复算全部 error count、Wilson 区间、CE 线性均值、RSRP 分位数和 bootstrap；审计通过后才生成图表；
7. 成对生成 `research/result-034-RSRP统计-PDCCH 4T4R.md` 与 `research/result-034-RSRP统计-PDCCH 4T4R-text.md`。研究者确认前不更新全局结论，不创建 Git checkpoint。

## 8. 输出、复现与证据

正式产物写入

`outputs/experiment034_rsrp_pdcch_4t4r/<run_id>/`

并至少包含：

- `configs/`：原 YAML、展开 YAML、SHA-256、三档递归 diff 回执；
- `validation/`：资源计数、precoder 行功率、phase denominator、4Rx 分支形状/相关性、无噪声编解码、噪声口径和 smoke 回执；
- `rsrp/`：每个 absolute trial 的线性/dB 块平均功率、单 RE 功率、trial key、候选 ID；CDF 点、分位数、DKW 带、paired-bootstrap 差值与 bootstrap seed；
- `bler/{low,mid,high}/`：逐点 CSV、error flags、CE NMSE 数组、展开配置、运行日志和停止原因；
- `analysis/`：合并九点 BLER CSV、RSRP 汇总表、审计 JSON、本地图的精确输入 CSV、样式 JSON、主图和约 13 cm 宽预览；
- `environment_receipt.json`：运行环境、代码版本/变更标识、设备、耗时和峰值 RSS。

计划复现命令为：

```powershell
python tools/run_plan034_rsrp_cdf.py --config configs/pdcch_plan034_rsrp_cdf.yaml --stage validate
python tools/run_pdcch_bler_curves.py --config configs/pdcch_plan034_fixed_dft0_4rx_bler_low.yaml --stage validate
python tools/run_plan034_rsrp_cdf.py --config configs/pdcch_plan034_rsrp_cdf.yaml --stage run
python tools/run_pdcch_bler_curves.py --config configs/pdcch_plan034_fixed_dft0_4rx_bler_low.yaml --stage run
python tools/run_pdcch_bler_curves.py --config configs/pdcch_plan034_fixed_dft0_4rx_bler_mid.yaml --stage run
python tools/run_pdcch_bler_curves.py --config configs/pdcch_plan034_fixed_dft0_4rx_bler_high.yaml --stage run
python tools/analyze_plan034_rsrp_pdcch_4t4r.py --bootstrap-repeats 1000
```

实际执行时必须把 `<run_id>`、配置 SHA-256 和完整命令写入 result。图只能由本地脚本读取已保存原始数据生成；不得把原始仿真数据交给云端绘图或图像生成工具。

## 9. result-034 必须回答的问题

两版 result 的配置、数字、判定和限制必须一致，并逐项回答：

1. `FIXED_DFT0` 与 `FREQ_SIDON_0137` 的 $Q_{0.01}$、$Q_{0.10}$、$Q_{0.50}$、$Q_{0.90}$、$Q_{0.90}-Q_{0.10}$、均值和标准差分别是多少；
2. $\Delta Q_{0.01}$、$\Delta Q_{0.10}$、$\Delta B_{80}$ 及 paired-bootstrap 95%区间是否满足第 4.3 节的集中性判据；
3. 单 RE 边缘控制是否通过；若不通过，是否已停止机理解释并定位实现或模型原因；
4. 九个固定 SNR 点各自的 BLER、errors/trials、Wilson 95%区间、CE NMSE 和停止原因是什么；
5. `0 dB` 及以下是否确实高于 0.1，`3--5 dB` 是否接近 0.01；不符合预期时必须直接报告；
6. 所有功率、4Rx 噪声、MRC、payload/CRC、资源计数、seed/absolute-trial 和 resume 审计是否通过；
7. 结论只适用于当前独立 4Tx/4Rx、TDL-C 300 ns、3 km/h、AL1 单 bundle 和指定接收机；不得外推到相关天线、多 bundle、其他 AL、其他 PDP、实际 RSRP 测量或标准 PMI 行为。

图文版至少包含一张 RSRP-CDF 主图和一张固定 precoder BLER 图；无图版提供完整关键表、图输入 CSV 路径和审计路径，不加载或描述图片像素。若 100,000 个 RSRP trials 或任一九点 BLER 达到预定上限但统计仍不足，result 标记“样本上限/部分不可判定”，不后验改判据。

## 10. 2026-09-16 补充规划：逐 Rx、四种预编码的无噪声块 RSRP

### 10.1 目的、继承范围与候选

第 4.1 节原统计量的分母含 $N_r$，并对 $r=0,\ldots,3$ 求和，因此它是四根 Rx 所见真实等效信道功率的算术平均。本补充保留该已执行统计及其结果不变，新增“不跨 Rx 平均”的诊断：分别对 Rx0、Rx1、Rx2、Rx3 计算经验 CDF，用于确认四个独立同分布 Rx 分支上的结论是否一致，并比较四种冻结发射波形。

除本节明确修改的统计口径、候选集合、输出和分析外，逐字段继承第 2--4 节的 4Tx/4Rx、1 layer、TDL-C 300 ns、3 km/h、4 GHz、20 sinusoids、AL1、2-symbol、每 symbol 36 个占用子载波、72 个候选 RE、单位总发射功率、无 AWGN、无逐 realization 归一化和相同预编码用于两个 symbols 等条件。四种预编码为：

| candidate ID | 类型 | `delay_grid_coordinates` | 物理人工时延 / ns | `phase_denominator` |
|---|---|---|---|---:|
| `FIXED_DFT0` | 固定 DFT index 0，$[1,1,1,1]^T/2$ | 不适用 | 无人工时延（等效相位可写为 `[0,0,0,0] ns`） | 不适用 |
| `FREQ_SIDON_0137` | plan-034 当前严格 Sidon | `[0,1,3,7]` | `[0,925.925926,2777.777778,6481.481481]` | 36 |
| `CDD911` | plan-031 的重复时延对 | `[0,0,0.98388,0.98388]` | `[0,0,911,911]` | 36 |
| `CDD130` | plan-031 的重复时延对 | `[0,0,0.1404,0.1404]` | `[0,0,130,130]` | 36 |

表中 $q_{36}=925.925926\ \mathrm{ns}$。三种频变候选均使用

$$
W_{k,n}=\frac12\exp\!\left(-j2\pi k j_n/36\right),
\qquad k=0,\ldots,35.
$$

`CDD911` 和 `CDD130` 的两个 0 ns 与两个非零时延均是预定波形的一部分，必须保留重复值，不得去重；实现配置显式设置 `allow_duplicate_delays: true`。本补充只使用上述发射波形计算真实 $\mathbf H\mathbf W$，不涉及 plan-031 中 matched/transparent 接收机协方差的区别。

### 10.2 逐 Rx 统计量与配对

对第 $t$ 个 realization、第 $r$ 根接收天线和预编码 $\mathbf W$，定义

$$
\overline P_{t,r}(\mathbf W)=
\frac{1}{|\mathcal Q|P_{\rm tx}}
\sum_{q\in\mathcal Q}
\left|
\sum_{n=0}^{N_t-1}H_{t,r,n,q}W_{q,n}
\right|^2,
$$

以及

$$
X_{t,r}(\mathbf W)=10\log_{10}\overline P_{t,r}(\mathbf W)\quad[\mathrm{dB}].
$$

该定义不含对 $r$ 的求和或 $1/N_r$；Rx0--Rx3 分别形成四组统计，禁止把四根 Rx 拼接为 400,000 个样本后只报告一条 CDF。每个 Rx、每个 candidate 均单独报告线性域均值/标准差、dB 域均值/标准差、$Q_{0.01}$、$Q_{0.10}$、$Q_{0.50}$、$Q_{0.90}$ 和 $Q_{0.90}-Q_{0.10}$。

正式样本仍为 100,000 个 absolute trials。为使新增结果可与原始两候选统计逐 trial 核对，信道生成必须重放原始 base seed `20260916`、channel namespace `plan034_rsrp_cdf_v1` 和 absolute trial 1--100000；新结果写入独立补充目录，不得覆盖原始输出。每个 $(t,r)$ 上四个 candidate 必须使用完全相同的底层 $H_{t,r,n,q}$，四根 Rx 也必须来自同一个 4Rx realization，不得为每个 candidate 或每根 Rx 另行生成信道。

### 10.3 不确定性、Rx 一致性与判定边界

每根 Rx 上，以 `FIXED_DFT0` 为基线，分别报告另外三种波形的 $\Delta Q_{0.01}$、$\Delta Q_{0.10}$ 和 $\Delta B_{80}$。同一 Rx 内的 candidate 差值使用 1000 次 paired nonparametric bootstrap，按 absolute trial 成组重采样，bootstrap seed 固定为 `20260916`，报告 percentile 95% 区间。不得把一个 Rx 上的结果复制到其他 Rx，也不得先跨 Rx 平均再做 bootstrap。

当前模型的 16 个 Tx--Rx 分支独立同分布，因此事实层面的理论预期是：对任一固定 candidate，Rx0--Rx3 的总体分布相同；有限样本的经验分位数和 CDF 不要求逐元素相等。result 应分别展示四根 Rx，并报告同一 candidate 在四根 Rx 间各指标的最大值、最小值和极差。另以 16 条经验 CDF（4 candidates $\times$ 4 Rx）的 family-wise 95% Bonferroni-DKW 带作一致性审计：令 $M=16$、$T=100000$，每条带的半宽为

$$
\epsilon_{\rm FWER}=
\sqrt{\frac{\ln(2M/0.05)}{2T}}.
$$

对同一 candidate，若任意两个 Rx 的经验 CDF 最大绝对差不超过 $2\epsilon_{\rm FWER}$，则记为“与独立同分布 Rx 假设一致”；超过时只记为需要审计的异常，不据此直接断言总体分布不同，应先检查 seed、分支索引、归一化和信道生成。该审计不用于选择最佳 Rx、丢弃 Rx 或把四根 Rx 合并。

### 10.4 实现、验证、输出与 result 要求

复用 `tools/run_plan034_rsrp_cdf.py` 和 `tools/analyze_plan034_rsrp_pdcch_4t4r.py`，只做支持四候选和逐 Rx 原始功率落盘的最小扩展，不复制信道生成或预编码主循环。新增配置使用独立输出子目录 `outputs/experiment034_rsrp_pdcch_4t4r/20260916_plan034/per_rx_4waveform/`。原始 trial 文件至少显式记录 `absolute_trial`、`rx_index`、`candidate_id`、线性/dB 功率、精确 delay 数组、grid coordinates、phase denominator、channel key 和 $P_{\rm tx}=1$ 标签。

正式运行前先执行 validate 和 100-trial smoke。validate/smoke 必须核对：四个 candidate 白名单；表中精确时延和坐标；CDD 重复时延未被去重；每个 RE 的预编码行功率为 1；输出恰有四个 Rx 索引；同一 $(t,r)$ 的四候选 channel key 相同；逐 Rx 统计不含 $1/N_r$；由四个逐 Rx 线性功率的算术平均可逐 trial 重构原第 4.1 节的 $\overline P_t$，且 `FIXED_DFT0`、`FREQ_SIDON_0137` 的重构值与原始保存结果逐元素一致。任一核对失败则停止正式分析并修复，不接受只在 dB 域近似相等。

补充分析至少生成：逐 Rx 的四候选 CDF 图（四个 Rx panel，不合并样本）、完整统计表、相对 `FIXED_DFT0` 的 paired-bootstrap 差值表、Rx 间极差表、Bonferroni-DKW 一致性审计和包含输入 SHA-256 的分析回执。图和表必须在图例或列中同时写明 candidate、Rx 索引及人工时延；其中 `CDD911` 明示 `[0,0,911,911] ns`，`CDD130` 明示 `[0,0,130,130] ns`。

补充实验的复现命令为：

```powershell
python tools/run_plan034_rsrp_cdf.py --config configs/pdcch_plan034_rsrp_per_rx_4waveform.yaml --stage validate
python tools/run_plan034_rsrp_cdf.py --config configs/pdcch_plan034_rsrp_per_rx_4waveform.yaml --stage smoke
python tools/run_plan034_rsrp_cdf.py --config configs/pdcch_plan034_rsrp_per_rx_4waveform.yaml --stage run
python tools/analyze_plan034_rsrp_pdcch_4t4r.py --mode per-rx --input-root outputs/experiment034_rsrp_pdcch_4t4r/20260916_plan034/per_rx_4waveform --bootstrap-repeats 1000
```

执行完成后，在 `research/result-034-RSRP统计-PDCCH 4T4R.md` 与 `research/result-034-RSRP统计-PDCCH 4T4R-text.md` 增加同名配对章节，并回答：四种预编码在 Rx0--Rx3 上的低尾、宽度、均值和标准差；每根 Rx 上三种频变方案相对固定预编码的差值及 95%区间；四根 Rx 的经验分布是否通过一致性审计；逐 Rx 线性功率能否精确重构原四 Rx 平均统计。该补充仍是无噪声机理诊断，不新增或推断任何 candidate 的 BLER 结论；研究者确认前不更新 `KNOWLEDGE.md` 或 `GOALS.md`，也不创建 Git checkpoint。
