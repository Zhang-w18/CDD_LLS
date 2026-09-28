# plan-039：CDL-C 角域全覆盖超宽波束 CDD 与 cycling BLER（草案）

> 状态：草案。场景已收敛为仅 `CDL-C + 双极化 2R UE`，AoD 目标 ASD 冻结为
> $25^\circ$；正式 trial 1 前仍需由研究者确认第 9 节统计预算和第 12 节冻结区。
> 本次改版改变了 CDL 角度、波束定义、Sidon 选择和接收机曲线集合；此前
> `outputs/experiment039_cdl_beam_bler/prepare/` 下使用旧配置生成的产物不得作为本 plan 的
> validate、smoke 或性能证据，也不得复用其 beam/statistics cache。
>
> 2026-09-28 补充冻结：增加一条 `60 km/h + 40 ms outdated CSI + PRG6 MRT` 的
> estimated-CSI BLER 曲线。不做 prescan，直接在 `7.5:0.5:20 dB` 运行正式仿真；该补充曲线
> 改变速度和发射端知识，必须与原 3 km/h 三曲线分组报告，不作配对增益归因。

## 1. 目的、研究问题与解释边界

本轮使用 `fixed_cdl_statistics` 后端，在经角度平移和缩放的 CDL-C 上构造一个
“角域全覆盖超宽参考波束”和 8 个正交水平 DFT 窄波束，并在相同 8 维波束子空间内比较：

1. 等差时延 `B0_QC` CDD；
2. 从严格 Sidon 集中按 DMRS 折叠距离和独立 ideal-CSI 预扫描选择的 `SIDON_SELECTED` CDD；
3. 8 个 6-RB PRG 依次使用 8 个窄波束的 `PRG6_CYCLING`。

补充实验另运行 `PLAN039_AGED_MRT_PRG6`：UE 速度为 60 km/h，发射端只使用当前 PDSCH
slot 起点之前 40 ms 的同一 realization 完美物理 CSI，在每个 6-RB PRG 内由两根 Rx 和该 PRG
全部子载波的 Gram 矩阵最大特征向量形成未量化 rank-1 MRT。该曲线用于给出明确反馈时延下的
闭环参考，不改变前三条 3 km/h 曲线的冻结定义。

实验先完成波束/信道审计、严格 Sidon 搜索和阶段 1A ideal-CSI/NMSE 预扫描。
根据 2026-09-26 的研究者决策，当前不再运行原定阶段 1B ideal-CSI 确认，而是直接在独立
evaluation seeds 上运行仅使用方法一公共 SSB 参考 PDP 的 estimated-CSI BLER。根据同日后续决策，
CDD 人工时延透明的全带 MMSE 对照也暂缓；已完成并落盘的透明数据保留为补充数据，但不再追加。
ideal-CSI 正式确认和方法二同样暂缓，不作为当前运行的前置 gate。

主要回答：

- AoD 平移与目标 ASD 缩放后，CDL-C 的角度功率是否主要落在
  AoD $[-60^\circ,60^\circ]$、ZoD $[90^\circ,110^\circ]$；
- 新参考超宽波束能否连续覆盖该角域，8 个窄波束是否正交、共同覆盖该角域，并使波束域
  等效信道相关性低于同均值、目标 ASD $10^\circ$ 的窄角诊断对照；
- 在严格 Sidon 集中优先最大化 DMRS fold 最小圆周距离后，能否选出 ideal-CSI BLER 同时优于
  `B0_QC` 和 `PRG6_CYCLING`、且 CE NMSE 不明显失稳的候选；
- 入围 Sidon 在两种预声明的频域协方差近似下，相对 B0 和 cycling 是否仍保留
  estimated-CSI 10%/1% BLER 优势。

本轮属于 `DESIGN.md` 的物理展宽 C1 场景。ideal CSI、角度覆盖、波束相关矩阵和 CE NMSE
只用于分阶段筛选或解释，不能替代第二阶段 estimated-CSI BLER 结论。当前已实现平台对
`wide_beam_split` 强制使用 profile-native angles，且旧波束是“先选最强 SSB、再在其局部角域
切 8 个窄波束”；二者均不满足本 plan，必须实现新 codebook type 和新的 cache key。

## 2. 冻结的资源、波形和链路条件

| 参数 | 冻结取值 |
|---|---|
| 链路 | 下行、rank 1、双极化 2Rx MRC |
| 资源 | 48 PRB，576 个连续 active subcarriers，30 kHz SCS，FFT 4096，CP 288 samples |
| PDSCH | 10 个 OFDM symbols；假设定时与载波同步，不模拟 CFO、ICI、ISI 或路径损耗/阴影衰落 |
| DMRS | symbols `[2,7]`，comb-6，offset 0；每 symbol 96 个 pilot RE，共 192 个 LS 观测；5568 data RE |
| PRG | 6 RB/PRG，共 8 个 PRG；从首个 active RB 起依次编号 |
| 调制编码 | 16QAM；NR 256QAM MCS table 的 MCS 8；目标码率 553/1024；Sionna LDPC 最多 8 次迭代 |
| 信道 | Sionna 1.0.2 CDL-C，RMS delay spread 100 ns，`fixed_cdl_statistics`，下行 |
| 载频与速度 | 4 GHz，3 km/h；最大经典 Doppler 约 11.12 Hz |
| 时间演化 | 同一 realization 的 10 个 OFDM symbols 按固定 ray Doppler 连续演化 |
| 发射知识 | 只使用 profile、阵列和独立长期统计；不使用待评估 trial 的瞬时或过时 CSI |
| 接收合并 | 两根 Rx 分别估计等效 rank-1 信道，data RE 上做 MRC；不做跨 Rx 空间 LMMSE |
| LLR | 沿用 plan-035，不额外加入 channel-estimation-error-aware 噪声项 |
| 执行设备 | CPU-only；保存 placement、batch size、墙钟时间和峰值 RSS |

补充 MRT 曲线除以下覆盖项外沿用本表：速度改为 60 km/h，最大经典 Doppler 约 222.38 Hz；
发射端使用 40 ms outdated CSI；MRT 权值在每个 PRG 内和整个 10-symbol slot 内固定。接收端仍使用
两个 DMRS symbol 的 `PLAN039_PRG_COMMON_REFERENCE_PDP` estimated-CSI LMMSE 和 MRC，不使用
outdated CSI。原 3 km/h 曲线及其 reference SNR 定义保持不变；补充曲线仍以同一超宽参考波束
长期接收功率定义噪声，因此横轴口径一致，但速度和发射知识不同，不能据此把差值归因为 MRT 本身。

人工 CDD delay 是频域线性相位对应的数字循环移位，不占用真实传播时延或 CP 预算。
CDL 路径仍按连续物理时延直接求频响；本轮不把物理路径量化到 FFT tap。

## 3. CDL-C 角度变换与唯一 UE 场景

### 3.1 阵列与场景

BS 使用 128 个 antenna elements（AE）和 32 个 TXRU：

$$
(M,N,P,M_g,N_g;M_p,N_p)=(8,8,2,1,1;2,8).
$$

水平、垂直 AE 间距分别为 $(d_H,d_V)=(0.5,0.8)\lambda$，交叉双极化，使用
38.901 阵元方向图。每极化有 $2\times8$ 个 TXRU；每个 TXRU 以等相位、单位范数权值连接
同一水平列内连续 4 个垂直 AE。

本轮只保留一个正式场景：

| 场景 ID | CDL | UE 阵列 | 极化 | UE 间距 |
|---|---|---|---|---|
| `C_UE2R_DUAL_ASD25` | CDL-C，100 ns | $(1,1,2,1,1;1,1)$ | cross dual | $(0.5,0.5)\lambda$ |

必须保存并审计 AE/TXRU/极化 flatten 顺序、映射矩阵、orientation、阵元方向图、极化类型和
所有间距，不得仅靠场景名推断。单极化 UE 和 CDL-E 均不属于本 plan。

### 3.2 AoD translation/scaling 的定义和原理

3GPP TR 38.901 的 CDL 角度可以按目标均值平移并按目标角扩展缩放。本轮对每个展开 ray 的
AoD 使用

$$
\phi_{n,\mathrm{scaled}}
=\mathrm{Wrap}_{[-180^\circ,180^\circ)}
\left[
s\,\mathrm{Wrap}_{[-180^\circ,180^\circ)}
(\phi_{n,\mathrm{model}}-\mu_{\phi,\mathrm{model}})
+\mu_{\phi,\mathrm{desired}}
\right].
$$

这里的 ASD 明确定义为全部 CDL-C 展开 ray 的功率加权 RMS azimuth spread of departure：

$$
\mathrm{ASD}
=\sqrt{
\frac{\sum_n p_n\,\Delta\phi_n^2}{\sum_n p_n}
},
\qquad
\Delta\phi_n=\mathrm{Wrap}(\phi_n-\mu_{\rm AoD}),
$$

其中 $p_n$ 是未经过发射波束方向图加权的标准 ray 功率，$\mu_{\rm AoD}$ 用功率加权圆均值
计算。冻结目标为

$$
\mu_{\rm AoD}=0^\circ,
\qquad
\mathrm{ASD}_{\rm target}=25^\circ,
\qquad
s=\frac{25^\circ}{\mathrm{ASD}_{\rm model}}.
$$

因此配置不得把 `aod_scale: 25` 误当作目标 ASD；validate 必须先从实际 Sionna CDL-C
展开数组计算 $\mathrm{ASD}_{\rm model}$，再得到无量纲 $s$。AoA、ZoD、ZoA 保持 profile-native
角度，只有 AoD 做上述变换。变换在长期统计固结前执行，对全部方案、CSI 模式和 trial 完全相同。

物理意图是让主要出射功率覆盖更宽的水平角域，从而降低 BS 阵面和 8 个波束域等效端口的
空间相关性；“相关性降低”是待验证目标，不由 ASD 数值自动推出。为明确“调大”的比较口径，
另生成同为 $\mu_{\rm AoD}=0^\circ$、目标 ASD $10^\circ$ 的 `ASD10_DIAGNOSTIC`，只用于角度谱和
相关矩阵对照，不运行 BLER、NMSE 或候选选择。作为配置自检，当前 CDL-C
表预期在上述变换后约 $94\%$ 的总 ray 功率落入 $[-60^\circ,60^\circ]$，主要强 cluster
中心约覆盖 $-30^\circ$ 到 $+59^\circ$。这些数值必须由本轮保存的 ray 表重算，不得直接抄作结果。

### 3.3 角度功率谱与相关性审计

在任何 BLER 预扫描前，本地生成并保存：

1. profile-native、ASD10 诊断对照和 ASD25 正式场景的 power-weighted AoD marginal spectrum；
2. 变换后 AoD--ZoD 二维角度功率谱，叠加 cluster center、展开 ray 和目标矩形
   $[-60^\circ,60^\circ]\times[90^\circ,110^\circ]$；
3. 每个 cluster/ray 的原始角度、变换后角度、delay、线性功率和累计功率表；
4. ASD10/ASD25 的 TXRU 域相关矩阵摘要，以及第 4 节 8 个窄波束的 $8\times8$ 等效信道相关热力图。

角度 gate 为：功率加权圆均值与 $0^\circ$ 的偏差不超过 $0.1^\circ$，ASD 与
$25^\circ$ 的偏差不超过 $0.1^\circ$，AoD 区间功率覆盖率不低于 0.93，二维目标矩形内功率
覆盖率不低于 0.92。未通过时停止，不允许通过改图轴或舍弃弱 ray 追认通过。

令第 $d$ 个独立小尺度 realization、接收端口 $r$、子载波 $k$ 上的 8 波束等效信道向量为

$$
\widetilde{\mathbf h}_{d,r}[k]
=\left[
\mathbf e_r^T\mathbf H_d[k]\mathbf w_0,\ldots,
\mathbf e_r^T\mathbf H_d[k]\mathbf w_7
\right]^T.
$$

以子载波为样本做平均，并以独立 realization 提供统计重复：

$$
\widehat{\mathbf R}_{b}
=\frac{1}{D N_r K}
\sum_{d=1}^{D}\sum_{r=1}^{N_r}\sum_{k=0}^{K-1}
\widetilde{\mathbf h}_{d,r}[k]\widetilde{\mathbf h}_{d,r}^{H}[k],
\qquad
C_{ij}=\frac{|[\widehat{\mathbf R}_{b}]_{ij}|}
{\sqrt{[\widehat{\mathbf R}_{b}]_{ii}[\widehat{\mathbf R}_{b}]_{jj}}}\in[0,1].
$$

CDL 小尺度系数理论均值为零，因此主定义不减样本均值。频域样本相关，禁止把 $K=576$
当作 576 个独立样本计算置信区间；bootstrap 以 realization 为单位。报告完整热力图、非对角
均值/中位数/95 分位数/最大值和有效秩。ASD10 与 ASD25 使用相同 seeds 和除 AoD scale 外
完全相同的配置；若 ASD25 的非对角均值未低于 ASD10，标记
`ANGULAR_SCALING_NO_CORRELATION_REDUCTION`，
停止 BLER，先复核角度、阵列与波束定义。

## 4. 角域全覆盖超宽参考波束与 8 个窄波束

### 4.1 双极化和 Kronecker 顺序

每极化的 TXRU 网格为垂直 2、水平 8。冻结双极化共同权重

$$
\mathbf p=\frac{1}{\sqrt2}[1,1]^T.
$$

以下 Kronecker 积按实际 manifest 的 `(polarization, vertical, horizontal)` 顺序展开；若代码
内部 flatten 顺序不同，只允许显式 permutation，不得改变物理权重。所有波束必须满足单位范数。

### 4.2 8 个正交水平 DFT 窄波束

定义 8 个水平空间频率

$$
\xi_m=\frac{-7+2m}{8},
\qquad m=0,\ldots,7,
$$

以及长度 8 的 shifted-DFT vectors

$$
u_m[n]=\frac{1}{\sqrt8}\exp(j\pi n\xi_m),
\qquad n=0,\ldots,7.
$$

相邻 $\xi_m$ 相差 $1/4$，所以 $\mathbf u_i^H\mathbf u_j=\delta_{ij}$；在 $d_H=0.5\lambda$
且 ZoD 接近 $90^\circ$ 时，其中心约从 AoD $-61^\circ$ 排到 $+61^\circ$，用于共同覆盖
$[-60^\circ,60^\circ]$。垂直方向只设计一个覆盖 ZoD $[90^\circ,110^\circ]$ 的公共宽波束
$\mathbf v_V^{\rm wide}\in\mathbb C^2$。8 个 32-TXRU 窄波束为

$$
\mathbf w_m
=\mathbf p\otimes\mathbf v_V^{\rm wide}\otimes\mathbf u_m,
\qquad m=0,\ldots,7.
$$

由于三个因子单位范数且 $\mathbf u_m$ 正交，必须数值验证
$\mathbf W^H\mathbf W=\mathbf I_8$。波束编号按中心 AoD 从负到正固定为 0--7；CDD delay 与
PRG cycling 均按此顺序映射，不按后验波束功率重排。

### 4.3 角域全覆盖超宽参考波束

另行设计水平超宽权重 $\mathbf v_H^{\rm wide}\in\mathbb C^8$，使归一化方向图尽量平坦覆盖
AoD $[-60^\circ,60^\circ]$；垂直权重 $\mathbf v_V^{\rm wide}$ 覆盖
ZoD $[90^\circ,110^\circ]$。参考波束定义为

$$
\mathbf w_{\rm SSB}
=\mathbf p\otimes\mathbf v_V^{\rm wide}\otimes\mathbf v_H^{\rm wide}.
$$

该参考波束命名为 `ANGULAR_FULL_COVERAGE_ULTRAWIDE`。允许牺牲峰值阵列增益和覆盖距离；
优化目标的优先级固定为：最大化目标二维角域内的最小增益、最小化角域内
95%--5% 分位增益纹波、最大化角域内平均增益，最后才最小化角域外泄漏。优化只读取阵列几何、
阵元方向图和预冻结角域网格，不读取 CDL ray 功率、BLER、NMSE 或正式 trial。候选算法、
正则化网格和 tie-break 必须在 validate 配置中冻结，不能按性能结果换波束。

“包住 8 个窄波束”指参考波束在整个目标二维矩形连续无未声明深零点，并覆盖 8 个窄波束的
中心方向；不要求参考波束的绝对增益逐点高于窄波束。参考波束与 8 个窄波束同图时必须同时
给出各自峰值归一化图和统一绝对增益图，避免用分别归一化掩盖参考波束的增益代价。

### 4.4 波束图和验收

绘图与数据格式参考 `tools/run_plan037_cdl_platform.py` 和
`outputs/experiment037_cdl_platform/beam_patterns_e_c/cdl_c/`，但不得复用其中的旧波束权重。
本轮至少生成：

- 变换后 CDL-C AoD--ZoD 角度功率谱与目标矩形；
- 参考超宽波束和 8 个窄波束的二维方向图；
- 固定 ZoD 的水平切面、固定 AoD 的垂直切面，以及统一绝对增益叠图；
- 参考波束 in-band minimum/mean/percentiles/ripple、峰值和角域外泄漏表；
- 8 个窄波束的 Gram matrix、目标角域采样覆盖和第 3.3 节信道相关热力图。

数值 gate 为：每个波束范数误差不超过 $10^{-12}$，
$\|\mathbf W^H\mathbf W-\mathbf I_8\|_{\max}\le10^{-12}$，8 个中心全部位于目标角域边界
$1.5^\circ$ 容差内，参考波束在预冻结二维采样网格上无非有限值且无零增益采样点。
平坦度作为本轮需要最优化并完整报告的设计指标，不在看到图后追加阈值。

## 5. 发射方案和严格 Sidon 选择

### 5.1 CDD 相位和 B0 基线

CDD 相位以首个 active subcarrier 起算的局部索引 $k=0,\ldots,575$ 定义为

$$
e^{-j2\pi k j_m/576},
\qquad
d_m=\frac{j_m}{576\times30\ \mathrm{kHz}}.
$$

一个 delay-grid 单位为 57.870370370 ns。因为 8 个窄波束严格正交，CDD 合成固定使用
$1/\sqrt8$，不得沿用旧 `wide_beam_split` 的逐子载波 $\alpha[k]$。基线为：

| 方案 ID | delay-grid index $\mathbf j$ |
|---|---|
| `B0_QC` | `[0,9,18,27,36,45,54,63]` |

### 5.2 严格 Sidon 搜索、预筛选和冻结

在 $\mathbb Z_{576}$ 上生成严格 Sidon 候选：36 个无序二元和
$(j_a+j_b)\bmod576$ 必须全部不同。comb-6 下每个 DMRS symbol 有
$N_p=576/6=96$ 个频域观测位置，定义

$$
g_{\rm fold}(\mathbf j)
=\min_{a\ne b}d_{\mathbb Z_{96}}(j_a\bmod96,j_b\bmod96).
$$

搜索和排序规则在读取性能数据前冻结：

1. 规范化 $j_0=0$，按公共移位、反射和分支置换去重；
2. 首要目标是在严格 Sidon 集中最大化 $g_{\rm fold}$；理论上限
   $\lfloor96/8\rfloor=12$ 必须写入审计，但不得预设一定可达到；
3. 同一最大 $g_{\rm fold}$ 下，依次最大化模 576 的 pair-sum 最小圆周距离、最小化
   $\max_m j_m$，最后按 delay vector 字典序稳定排序；
4. 保存排序前 8 个候选；按研究者 2026-09-26 的执行指令，本次先只对几何 rank 1
   `[0,11,28,148,170,233,277,351]` 运行 ideal-CSI/NMSE 预扫描，其余 7 个候选标记为
   `DEFERRED_NOT_RUN`；后续追加时必须保留本次数据和 selection seed，不得根据 rank 1 的结果
   只挑选部分候选；
5. 当前单候选 prescan 自动把该候选冻结为 `SIDON_SELECTED`；若后续按预声明顺序追加多个候选，
   再按 10% ideal-CSI crossing 最低、1% crossing 最低、方法二 NMSE penalty 最小、几何排序最
   靠前的顺序自动重算唯一 `SIDON_SELECTED`，不得人工挑选；
6. 候选一经选择即冻结，其正式确认和第二阶段使用与 selection seeds 不相交的 evaluation seeds。

必须为每个候选保存精确有理数、秒、FFT-sample 等效值、comb-6 residue、36 个无序二元和、
`pair_gap`、`fold_gap`、pilot rank、condition number 和近零噪声 CE 检查。严格 Sidon 只是一项
候选生成条件，不自动保证物理展宽 CDL-C 中的 outage、BLER 或 CE 优势。

### 5.3 Precoder cycling

`PRG6_CYCLING` 不叠加人工 delay。8 个 6-RB PRG 依次使用窄波束 index
`[0,1,2,3,4,5,6,7]`，每个 PRG 内权值恒定，trial 间不随机重排。manifest 保存每个
RB/subcarrier 到 PRG 和波束的完整映射。

## 6. SNR、ideal CSI 和两种 estimated-CSI 协方差

### 6.1 公共参考 SNR

所有方案在每个 active subcarrier 上满足 $\|\mathbf w[k]\|_2^2=1$。固定参考功率为

$$
P_{\rm ref}=\mathbf w_{\rm SSB}^{H}\mathbf R_t\mathbf w_{\rm SSB},
\qquad
N_0=P_{\rm ref}10^{-\gamma_{\rm ref,dB}/10},
$$

其中 $\mathbf R_t=\mathbb E[\mathbf H^H\mathbf H]$ 包含两根 Rx 的总接收功率。
同一场景所有方案、CSI 模式和协方差方法共用同一个 $P_{\rm ref}$ 与每根 Rx 的 $N_0$；不得按
方案重新归一化接收功率。参考超宽波束允许牺牲绝对增益，因此横轴只表示该参考口径下的相对
性能，不得外推为小区覆盖距离。

### 6.2 曲线集合

当前阶段 1B 冻结为 3 条 estimated-CSI BLER，并在同一 trial 中记录 3 条 CE NMSE：

| 发射方案 | 方法一，知道 CDD delay | 透明全带 MMSE | PRG 内方法一 |
|---|---:|---:|---:|
| `B0_QC` | 1 | 0 | 0 |
| `SIDON_SELECTED` | 1 | 0 | 0 |
| `PRG6_CYCLING` | 0 | 0 | 1 |

阶段 1A 已有的 ideal CSI 只用于趋势和候选选择，不与当前 evaluation 数据合并。原计划的
三条 ideal-CSI 正式曲线、CDD 方法二和透明全带 MMSE 曲线均标记为 `DEFERRED_NOT_RUN`；
本决策前已落盘的透明点例外标记为 `RETAINED_SUPPLEMENT_NOT_EXTENDED`。

### 6.3 共同二维时频 LMMSE 和两种频域协方差

estimated CSI 使用两个 DMRS symbol 的独立 LS 观测，不做 DMRS symbol 平均。时间协方差
沿用 plan-035/Sionna 的 3 km/h、4 GHz Jakes 模型，并与频率协方差作可分离组合；这一模型
不等同于包含每条 CDL ray 方向 Doppler 的完整非可分离统计。

方法一、方法二的等效 PDP 和频域协方差实现以仓库根目录
`CDL模式下PDP计算方法.md` 第 1--3、5--6 节为规范来源；本节把其中一般 $K$ 分支公式具体化为
$K=8$。若本节与实现注释或旧 plan 冲突，以本节和该规范文档的公共时延轴、绝对路径功率、
PDP 平移后等功率合成定义为准。

以下公式对每根 Rx 分别成立，为简洁省略 Rx 上标。第 $m$ 个窄波束的物理频响写为

$$
h_m[k]=\sum_\ell a_{m,\ell}e^{-j2\pi f_k\tau_\ell},
$$

CDD 合成后的等效信道为

$$
g_{\mathbf j}[k]
=\sum_{m=0}^{7}c_m[k]h_m[k],
\qquad
c_m[k]=\frac{1}{\sqrt8}e^{-j2\pi k j_m/576}.
$$

令

$$
\mathbf a_\ell=[a_{0,\ell},\ldots,a_{7,\ell}]^T,
\qquad
\mathbf C_\ell=\mathbb E[\mathbf a_\ell\mathbf a_\ell^H],
\qquad
\mathbf c[k]=[c_0[k],\ldots,c_7[k]]^T.
$$

在不同物理 delay 分量互不相关的统计模型下，CDD 等效频域协方差的一般形式是

$$
R_f[k,k']
=\mathbb E[g_{\mathbf j}[k]g_{\mathbf j}^*[k']]
=\sum_\ell e^{-j2\pi(f_k-f_{k'})\tau_\ell}
\mathbf c^T[k]\mathbf C_\ell\mathbf c^*[k'].
$$

两种方法的区别不是 LMMSE 公式不同，而是接收机如何近似未知的
$\mathbf C_\ell$。令 $p_{{\rm ref},\ell}$ 为参考超宽波束在公共物理 delay
$\tau_\ell$ 上的绝对线性路径功率，$p_{m,\ell}$ 为第 $m$ 个窄波束的绝对线性路径功率。
若内部保存归一化 PDP shape，构造协方差前必须乘回对应总接收功率；不得只使用单位和为 1 的
shape 而丢失与噪声方差有关的绝对尺度。

#### 方法一：公共参考 PDP 复制

`common_reference_pdp` 只要求接收机知道参考超宽波束测得或由长期统计得到的一条 PDP，以及
已知的 CDD delay vector。它假设 8 个窄波束分支相互独立、具有完全相同的路径功率：

$$
\mathbf C_\ell^{(1)}=p_{{\rm ref},\ell}\mathbf I_8.
$$

其实现主定义是把同一条参考 PDP 复制 8 份，并按各分支人工时延整体平移后等功率相加：

$$
P_{\rm ref}(\tau)
=\sum_\ell p_{{\rm ref},\ell}\delta(\tau-\tau_\ell),
$$

$$
\boxed{
P_{\rm eq}^{(1)}(\tau)
=\frac{1}{8}\sum_{m=0}^{7}P_{\rm ref}(\tau-d_m)
=\frac{1}{8}\sum_{m=0}^{7}\sum_\ell p_{{\rm ref},\ell}
\delta(\tau-\tau_\ell-d_m)
}.
$$

从该等效 PDP 直接得到频域协方差

$$
R_f^{(1)}[k,k']
=\frac{1}{8}\sum_{m=0}^{7}\sum_\ell p_{{\rm ref},\ell}
e^{-j2\pi(f_k-f_{k'})(\tau_\ell+d_m)}.
$$

由于 $(f_k-f_{k'})d_m=(k-k')j_m/576$，该式与第 5.1 节局部 active-subcarrier
相位定义严格一致。实现必须用这一等价关系做逐元素测试，不允许把相位分母误换为 FFT 4096。

因此方法一并非忽略 CDD：它保留接收端已知的 8 组人工线性相位，但把每个窄波束真实不同的
角度选择性、总功率和 PDP 形状全部替换为同一条参考 PDP，并把波束间交叉相关设为零。这是
只维护一条参考 PDP 的低复杂度、预声明失配接收机，也对应本轮三种发射方案能够共用的接收
知识口径。

#### 方法二：逐窄波束 PDP、分支独立

`beam_specific_pdp_independent` 要求接收机知道 8 条窄波束各自的长期 PDP。它保留每个波束
对不同 angle--delay ray 的选择作用和总接收功率差异，但仍假设不同波束分支互不相关：

$$
\mathbf C_\ell^{(2)}
=\operatorname{diag}
\left(p_{0,\ell},p_{1,\ell},\ldots,p_{7,\ell}\right),
$$

其实现主定义是每条窄波束 PDP 只平移本分支的人工时延，再做等功率合成：

$$
P_m(\tau)=\sum_\ell p_{m,\ell}\delta(\tau-\tau_\ell),
$$

$$
\boxed{
P_{\rm eq}^{(2)}(\tau)
=\frac{1}{8}\sum_{m=0}^{7}P_m(\tau-d_m)
=\frac{1}{8}\sum_{m=0}^{7}\sum_\ell p_{m,\ell}
\delta(\tau-\tau_\ell-d_m)
}.
$$

对应频域协方差为

$$
R_f^{(2)}[k,k']
=\frac{1}{8}\sum_{m=0}^{7}\sum_\ell p_{m,\ell}
e^{-j2\pi(f_k-f_{k'})(\tau_\ell+d_m)}.
$$

#### 频域协方差矩阵的直接数值实现

连续形式的 $P_{\rm eq}^{(1)}(\tau)$ 和 $P_{\rm eq}^{(2)}(\tau)$ 只用于解释和绘图；LMMSE
实现不得先把它们采样到 delay bin、FFT tap 或直方图后再反推协方差。代码直接读取公共物理
delay 数组 $\{\tau_\ell\}_{\ell=1}^{L}$、人工 delay 数组
$\{d_m\}_{m=1}^{K}$ 和对应线性路径功率，并直接按下列双重求和计算每个矩阵元素。

方法一的实现公式为

$$
\boxed{
R_f^{(1)}[k,k']
=
\frac{1}{K}
\sum_{m=1}^{K}
\sum_{\ell=1}^{L}
p_{{\rm ref},\ell}
e^{-j2\pi(f_k-f_{k'})(\tau_\ell+d_m)}
},
\qquad K=8.
$$

方法二只把公共路径功率替换为波束特定路径功率：

$$
\boxed{
R_f^{(2)}[k,k']
=
\frac{1}{K}
\sum_{m=1}^{K}
\sum_{\ell=1}^{L}
p_{m,\ell}
e^{-j2\pi(f_k-f_{k'})(\tau_\ell+d_m)}
},
\qquad K=8.
$$

实现时先构造 $576\times576$ 的差频矩阵

$$
\Delta F_{kk'}=f_k-f_{k'}=(k-k')\Delta f,
\qquad \Delta f=30\ \mathrm{kHz},
$$

再对每个 $(m,\ell)$ 计算
$\exp[-j2\pi\Delta F(\tau_\ell+d_m)]$ 并按
$p_{{\rm ref},\ell}/K$ 或 $p_{m,\ell}/K$ 加权累加。可以对 $(m,\ell)$ 分块或广播以控制内存，
但主计算必须保持上述双重求和语义；不能截断弱路径、重采样 delay、量化到 FFT tap，或把
分母 $K=8$ 换成 FFT 4096/active-subcarrier 576。

矩阵形式 $\mathbf A\operatorname{diag}(\mathbf Q)\mathbf A^H$ 只作为单元测试中的等价计算，
不作为本轮主实现路径。方法一、方法二都必须检查：双重求和与该等价式的最大逐元素误差、
Hermitian 误差、最小特征值，以及

$$
R_f^{(1)}[k,k]=\sum_\ell p_{{\rm ref},\ell},
\qquad
R_f^{(2)}[k,k]=\frac{1}{K}\sum_{m=1}^{K}\sum_{\ell=1}^{L}p_{m,\ell}.
$$

方法二比方法一多使用“每个窄波束自己的边缘二阶统计”，但没有使用
$\mathbb E[a_{m,\ell}a_{n,\ell}^*]$（$m\ne n$）交叉项，因此即使
$\mathbf W^H\mathbf W=\mathbf I_8$，它也不能被称为 matched receiver：正交发射权重不保证
经过有限角扩展 CDL 信道后的 8 个等效信道统计独立。方法二相对方法一的 NMSE/BLER 差只衡量
“逐波束边缘 PDP 信息”的价值，不能解释为完整波束联合协方差的价值。

两种方法的实现顺序固定为：先在统一物理时延轴上得到
$p_{{\rm ref},\ell}$ 或 $p_{m,\ell}$；再形成 $\Delta F$ 和每个 $\tau_\ell+d_m$；最后按上述
双重求和直接计算 $\mathbf R_f^{(1)}$ 或 $\mathbf R_f^{(2)}$。数值相同的等效 delay 原子可以在
保存表格或绘图时合并，但协方差实现不依赖这一步。方法二不得把每条 $P_m$ 单独重新归一化为
相同总功率，否则会删除窄波束间真实总功率差异；仅用于画 shape 的归一化副本必须使用不同
字段名，不能作为 LMMSE 输入。

两种方法共享同一下行物理时延参考；任何波束 PDP 都不得逐波束把首径重新置零。两根 Rx
分别使用各自的长期路径功率构造 $R_f^{(1)}$ 或 $R_f^{(2)}$；若实现利用对称性复用同一矩阵，
必须先数值验证两根 Rx 的协方差在预声明容差内一致。二维 LMMSE 对 pilot 集 $P$ 和 data 集
$D$ 使用

$$
\widehat{\mathbf g}_D
=\mathbf R_{DP}
\left(\mathbf R_{PP}+\sigma_{\rm LS}^2\mathbf I\right)^{-1}
\widehat{\mathbf g}_{P}^{\rm LS},
$$

其中 $\mathbf R_{PP}$ 和 $\mathbf R_{DP}$ 从对应方法的频率协方差与公共时间协方差生成，
$\sigma_{\rm LS}^2$ 使用实际 DMRS 能量和噪声方差。两根 Rx 分别滤波，再在 data RE 上做
estimated MRC；记录 loading、condition number、求解残差和有限性检查。

#### Cycling 中的“方法一”

cycling 的每个 72-subcarrier PRG 只激活一个窄波束，因此 PRG 内不存在 8 分支 CDD 合成，
也没有人工相位因子 $c_m[k]$。这里称为“方法一”仅表示 8 个 PRG 都使用同一条参考超宽波束
PDP 构造各自的局部 LMMSE；它不是把 8 条参考 PDP 在一个 PRG 内相加。估计器不得跨 PRG
读取 pilot 或做协方差耦合。实际激活窄波束的 $p_{m,\ell}$ 与
$p_{{\rm ref},\ell}$ 之间的差异是预声明失配；本轮不为 cycling 增加逐波束 PDP 方法二曲线。

#### CDD 的透明全带 MMSE

`transparent_common_reference_pdp` 作为已实现但当前暂缓的补充接收机，与方法一使用完全相同的公共 SSB
参考 PDP，但接收机不知道也不
补偿 CDD delay vector；构造频域协方差时将人工时延统一取为 0：

$$
R_{f,\mathrm{transparent}}[k,k']
=\sum_\ell p_{{\rm ref},\ell}
e^{-j2\pi(f_k-f_{k'})\tau_\ell}.
$$

该协方差用于整个 576-subcarrier active band 的二维时频 MMSE，不在 PRG 边界分块。因此它与
cycling 的共同点是都不使用 CDD delay 调整，差异是 cycling 只在每个 6-RB PRG 内独立做 MMSE。

原草案的第三种 `beam_joint_covariance` 不实现、不运行、不计入曲线数或完成条件，标记为
`DEFERRED_RESEARCH`。后续若研究完整波束联合协方差，必须另建或更新 plan，明确非可分离时频
统计和与方法一/二的公平比较，不能把本轮数据后验补成“方法三”。

### 6.4 CE NMSE

仅对 estimated 接收机保存 trial 级线性 NMSE：先在两根 Rx 和全部 5568 个 data RE 上合并
误差能量与真实信道能量，再跨 trial 在线性域求均值，最后转 dB。保存 trial 值、sum、sum of
squares、均值和 95% Monte Carlo 区间；禁止先对 Rx、PRG、RE 或 trial 转 dB 后平均。

## 7. 分阶段实验与阶段门

### 7.1 阶段 0：实现、validate 和无结论 smoke

完成 CDL-C ASD25 角度变换、新波束、相关矩阵、严格 Sidon 搜索、两种 LMMSE 和 runner 后，
依次执行单元测试、`--stage validate` 和每个 SNR 20 个 paired trials 的 smoke。阶段 0 必须
通过第 3、4 节的角度/波束/相关 gate、单位功率、seed replay、直接 TXRU 与波束域合成一致、
DMRS/PRG 坐标和 LMMSE 有限性检查。smoke 不作性能结论。

### 7.2 阶段 1A：候选 selection prescan

先对几何 rank 1 严格 Sidon 候选 `[0,11,28,148,170,233,277,351]` 以及 B0/cycling，在公共网格

`[0,2,4,6,8,10,12,14,16,18,20] dB`

运行 400 个 selection paired trials/点。Sidon 候选只运行 ideal-CSI BLER；B0 和 cycling 的
ideal 曲线共享相同 realization/payload/noise。另对 B0、当前纳入 prescan 的 Sidon 候选的两种估计器以及
cycling 方法一运行 CE-only NMSE，不做 estimated-CSI LDPC 解码。若 10% crossing 无真实 bracket，
仅向缺失方向按 2 dB 扩展，硬边界暂定 `[-6,26] dB`。

按第 5.2 节规则自动冻结唯一 `SIDON_SELECTED`，保存 top-8 清单、当前运行候选及其结果，并把
未运行候选明确标记为 deferred。selection 数据不得与后续 evaluation 数据合并计算正式置信区间。

### 7.3 阶段 1B：独立 estimated-CSI BLER 确认

只对 `B0_QC`、`SIDON_SELECTED` 和 `PRG6_CYCLING` 使用独立 evaluation seeds 运行第 6.2 节的
3 条 estimated-CSI BLER/CE NMSE 曲线。CDD 的方法一在 7.5--20 dB 保留
0.5 dB 基础网格，只在 10--15 dB 加密到 0.25 dB；cycling 在同一总区间保留 0.5 dB
基础网格，只在 15--20 dB 加密到 0.25 dB。所有曲线执行 estimated-CSI LDPC 解码并同时保存 CE NMSE。

本阶段不运行 ideal CSI、方法二、CDD 透明接收或原计划的近零噪声方法二 gate。
本决策前已完整落盘的透明数据必须保留，但不参与后续自适应追加。阶段 1A 结果只用于
说明这一范围选择，不并入正式置信区间。

执行顺序可在 B0/Sidon 初始 1,000 trials 及第一轮配对追加落盘后暂停 CDD 自适应，
先通过 `--stage cycling_initial` 完成 cycling 的全网格初始 1,000 trials，不做 cycling 自适应追加。
研究者比较三条初始曲线后，再决定是否恢复自适应。后续追加必须从已落盘的
`absolute_trial_stop` 继续，按 variant、SNR 和不重叠的 absolute trial 区间合并。

### 7.4 暂缓项

三条独立 evaluation ideal-CSI 曲线、CDD 方法二/透明 estimated-CSI 曲线和原方法二 NMSE gate 暂缓。
后续若恢复，必须继续使用已冻结候选和 evaluation seed namespace，并通过新的 absolute trial 区间追加；
不得将阶段 1A selection trial 合并为正式数据。

2026-09-28 研究者恢复其中的 `SIDON_SELECTED` CDD 透明 estimated-CSI 曲线。不做
prescan；使用 `PLAN039_TRANSPARENT_COMMON_REFERENCE_PDP`、已冻结 Sidon delay、evaluation seeds
和方法一 Sidon 曲线已完成的相同 SNR 点及 absolute trial 区间直接正式扫描。runner
只从同时存在 `summary.csv`、`trial_metrics.csv` 和 `batch_receipt.json` 的方法一批次冻结
覆盖，不读取阶段 1A trial，不自适应增加新点或超过方法一的 trial 数。其余方法二、
ideal-CSI 和 NMSE gate 仍暂缓。

### 7.5 补充正式曲线：60 km/h、40 ms outdated-CSI MRT

本补充不做 selection/prescan，也不读取阶段 1A 性能来调整 SNR。唯一曲线为
`aged_mrt_prg6__60kmh__csi_age40ms`，发射方案 `PLAN039_AGED_MRT_PRG6`，接收机为
`PLAN039_PRG_COMMON_REFERENCE_PDP`。每个 absolute trial 先以同一组 CDL 随机初相位生成当前
slot 的 10 个 OFDM-symbol 信道，再以反向速度在 $+40$ ms 取样，等价得到同一 realization
在 $-40$ ms 的物理频域 CSI；必须核对两次生成的 $t=0$ 信道最大逐元素误差不超过
$10^{-10}$。不得把独立 realization 当作 outdated CSI。

对每个 6-RB PRG，令 $\mathbf H_{\rm old}[k]\in\mathbb C^{2\times32}$，冻结

$$
\mathbf G_{\rm old}^{(p)}
=\sum_{k\in p}\mathbf H_{\rm old}^{H}[k]\mathbf H_{\rm old}[k],
\qquad
\mathbf w_p=\operatorname{eigvec}_{\max}(\mathbf G_{\rm old}^{(p)}),
\qquad \|\mathbf w_p\|_2^2=1.
$$

不量化、不做波束码本投影，8 个 PRG 分别求权值。保存每 trial 的权值功率上下界、CSI age、
current/stale replay 误差、信道/载荷/噪声 seed 和 CE NMSE。该曲线不与 3 km/h 三曲线共享
realization，也不执行 paired bootstrap；只报告自身 BLER、Wilson 95% 区间和 CE NMSE。

## 8. 公平性、配对和随机种子

每个 `stage + SNR + absolute trial` 内，同阶段的全部方案共享 fixed-CDL realization、
transport block、编码比特和单位方差 data AWGN；estimated 接收机另共享单位方差 DMRS AWGN。
同一 CDD 方案的两种估计器复用完全相同的 DMRS LS 观测。不同方案只改变预编码；ideal 分支
不得另生成信道或 data noise。

阶段 1A selection 与阶段 1B evaluation 使用不相交 seed namespace。阶段 1B 某 SNR 点只要
任一 CDD 曲线需要追加 trial，该点 B0/Sidon 两条 CDD BLER/CE 记录均按相同 absolute-trial 区间追加；
cycling 独立按同样规则追加。两组在共同 SNR 点仍使用同一 deterministic seed 派生，保持可配对。
base seed 暂定 `20260939`；channel、payload、data noise、DMRS noise、statistics、
selection、evaluation 和 bootstrap 使用显式子命名空间，保存 seed derivation 与 replay audit。

同一场景的 gain、receiver penalty 和 ideal/estimated gap 使用 paired bootstrap。频域相关矩阵
的不确定性以独立 realization 为 bootstrap 单位，不以子载波为独立单位。

## 9. SNR 网格、统计预算和停止条件（待确认）

阶段 1B 的总 SNR 范围冻结为 7.5--20 dB，基础间隔为 0.5 dB。CDD 的方法一曲线
在 10--15 dB 加密到 0.25 dB；cycling 在 15--20 dB 加密到 0.25 dB。runner 必须把展开后的
两组精确网格及其 SHA-256 写入阶段报告。crossing 只使用真实双侧 bracket，不用 PAVA、
PCHIP、单调修正或外推制造。

每个 formal SNR 点采用：

1. 初始 1,000 个 paired trials；每批追加 1,000；单点上限 50,000；
2. 10% bracket 端点 BLER 位于 5%--20% 时，每端点至少 200 个 block errors；
3. 1% bracket 端点 BLER 位于 0.5%--2% 时，每端点至少 200 个 block errors；
4. crossing 使用相邻真实采样点的 log-BLER 线性插值；端点必须有正错误数；
5. 未被任何曲线命中为 bracket 的点完成初始 1,000 trials 后停止；
6. 达到 50,000 trials 仍不足时保留真实数据并标记 `capped`，不加伪计数；
7. 到 20 dB 仍无 bracket 时标记 `target_not_reached_at_snr_cap`，不扩展、不外推；
8. CE NMSE 与相应 BLER 使用完全相同的 evaluation trial，不运行用于正式结论的 NMSE-only
   追加 trial；阶段 1A 的 CE-only 数据只用于候选选择，不能并入正式 NMSE。

逐点 BLER 报告 Wilson 95% 区间；crossing 和方案 gain
使用 1,000 次预声明 paired bootstrap，并报告有效重复数。若曲线非单调导致多个真实
bracket，保存全部 bracket，并以最低 SNR 的首次下降 crossing 为主，同时报告敏感性表。

对未取得真实 crossing 的曲线，只报告最高已测 SNR、BLER/errors/trials/Wilson 区间和对应
NMSE；相关 gain 标记 `not_comparable_no_bracket`。不得根据结果后验修改预算、NMSE gate、
候选数、SNR cap 或主比较。

补充 aged-MRT 曲线冻结为 `7.5:0.5:20 dB` 共 26 点，每点固定 1,000 trials，总计
26,000 trials；不做 prescan、不加密、不自适应追加、不以错误数提前停止。若 10% 或 1%
crossing 有相邻真实双侧 bracket，可按本节相同 log-BLER 线性插值作描述性报告；否则只报告
实测点，不外推。正式运行前必须通过配置校验、aged/current 同 realization replay、40 ms 时间
索引、PRG 常值和单位功率定向测试；无需重复原 stage-0 角度/波束设计搜索。每个 SNR 点独立
保存展开配置、trial 数据和完成回执，最多并行运行 3 个互不重叠的 SNR 点；并行只改变调度，
不得改变逐点 seed 派生或 absolute trial 区间。

## 10. 实现范围、测试与 smoke

### 10.1 预计代码和配置范围

不得复制 PDSCH/CDL 主循环。预计工作为：

1. 扩展 `fixed_cdl_statistics` 配置和 cache，使新 codebook type 支持“仅 AoD 目标均值/ASD
   变换”，并把原始/变换角度、派生 scale 和目标 ASD 纳入 hash；
2. 在 `cdd_lls/phy/cdl_beam_platform.py` 中增加 shifted horizontal DFT、可分离水平/垂直
   超宽波束合成、角域覆盖指标和 8 波束相关矩阵；保留 037 原行为不变；
3. 在可复用模块实现严格 Sidon 最大 fold-gap 搜索、稳定去重/排序和候选 audit；
4. 扩展 rank-1 下行链路支持 ideal、两种 CDD LMMSE、CDD 透明全带 LMMSE、PRG-local LMMSE、trial 级 CE 指标和
   可恢复 interval 汇总；不实现第三种联合波束协方差；
5. 更新 `tools/run_plan039_cdl_beam_bler.py` 和 `configs/plan039_cdl_wide_beam_bler.yaml`，支持
   `validate/smoke/select/estimated_confirm/cycling_initial` 阶段；两个 evaluation stage 只能在
   selection freeze 存在后运行；
6. 新增或更新本地分析脚本，只读取保存的 CSV/NPY/JSON，完成角度谱、方向图、相关热力图、
   bracket、Wilson 区间、paired bootstrap、crossing、gain 和 NMSE；
7. 回归 plan-037 CDL 平台、`fixed_cdl_statistics` smoke 和现有通用 ideal 链路。
8. 增加 `--stage aged_mrt_formal`：从已通过的 ASD25 几何/波束 manifest 派生 60 km/h 时间
   协方差 manifest，只执行第 7.5 节固定 26 点正式批次，支持完整批次校验后幂等复用。

### 10.2 定向测试

至少覆盖：

1. 唯一双极化 UE 场景、BS 128 AE→32 TXRU 映射、极化/vertical/horizontal 顺序；
2. CDL-C 100 ns ray 表、AoD 圆均值、ASD scale、wrap、只变 AoD、目标角度和 seed replay；
3. 角度覆盖率、二维角度功率谱输入和 native/ASD10/ASD25 cache 完全隔离；
4. shifted-DFT 解析正交性、中心方向、共同垂直宽波束、双极化复制和波束范数；
5. 超宽参考波束优化目标、稳定 tie-break、方向图数值与全阵列直接求和一致；
6. 8 波束相关矩阵的手算小例子、归一化范围、realization 级 bootstrap 和热力图输入；
7. 严格 Sidon 二元和唯一性、最大 fold-gap、对称去重、稳定排序和 selection/evaluation seed 隔离；
8. B0/Sidon 相位、ns/FFT-sample 换算、分母 576、CDD 全子载波单位功率和直接合成一致；
9. 按 `CDL模式下PDP计算方法.md` 构造的方法一/二等效 delay--power 原子，与直接频域求和
   逐元素一致；同时检查 Hermitian/PSD/对角功率、公共 delay 轴，并确认第三种方法不可配置；
10. PRG cycling 映射准确，每个滤波器只读取本 PRG pilot；
11. estimated MRC、3 条 CE/BLER 记录、B0/Sidon 配对和 interval resume/merge；
12. plan-037/fixed-CDL/通用 PDSCH 定向回归不改变现有行为。
13. aged-MRT 的 40 ms 反向时间取样与同一 realization 重放一致、PRG Gram 最大特征向量使用
    两根 Rx、每子载波单位功率、CSI age 字段及 26 点固定正式网格。

### 10.3 Smoke

在两个远离/接近预期 waterfall 的 SNR 点各运行 20 个 paired trials。阶段 0 保留执行方案收缩前的
8 条功能性曲线覆盖，只用于检查 5 条 estimated 接收分支的
有限 CE NMSE、角度/波束 cache、单位功率、公共 $P_{\rm ref}$、噪声、DMRS/PRG 坐标、MRC
denominator、滤波器残差、seed 配对、CPU placement、RSS、耗时和数组引用。任一角度、几何、
相关性、功率、数值或配对审计失败即停止，不得进入 selection prescan。

## 11. 输出、图表、比较和 result 要求

正式产物写入：

`outputs/experiment039_cdl_beam_bler/<run_id>/<stage>/C_UE2R_DUAL_ASD25/`

至少交付：

- 原始/展开 YAML、SHA-256、代码版本/工作区变更标识和复现命令；
- CDL-C 标准 profile、原始/变换 ray 表、ASD/覆盖率 audit、阵列/TXRU manifest；
- 角度功率谱输入、超宽/窄波束权值、方向图输入、覆盖/纹波/绝对增益表和 cache manifest；
- ASD10/ASD25 的 TXRU 与 8 波束相关矩阵、热力图输入、汇总和 realization 级区间；
- 全部严格 Sidon 搜索空间摘要、top-8 候选、selection 指标、唯一入围候选和未入围负结果；
- B0/Sidon delay audit、CDD 单位功率、PRG cycling 映射和 receiver manifest；
- 参考波束 PDP、8 个窄波束 PDP、方法一与透明接收协方差；方法二实现审计保留但不生成本阶段曲线；
  condition number、loading 和 solve audit；
- selection 与 evaluation 完全分离的 interval CSV、trial error flags、CE trial arrays、seed/pairing
  audit、合并校验和与 adaptive history；
- BLER/NMSE 汇总、errors/trials/Wilson 区间、真实 bracket、crossing、paired-bootstrap gain/penalty 和
  所有未达到目标状态；
- CPU placement、batch/RSS/墙钟日志、失败或作废运行清单。
- 补充 aged-MRT 的 60 km/h 派生 manifest、展开 YAML、26 点 summary、26,000 行 trial metrics、
  40 ms replay audit 和不运行 prescan/自适应追加的冻结回执。

本地至少生成：角度功率谱、参考/8 窄波束方向图、统一增益切面、8×8 相关热力图、阶段 1A
ideal BLER 趋势、阶段 1B estimated BLER/CE NMSE 和 crossing/gain 汇总表。BLER 图使用
对数纵轴；零错误点在 CSV 保留 0，绘图时才使用 `0.5/trials` 下界并明确标记。所有图由本地
脚本读取已保存数据生成，绘图输入另存 CSV/NPZ；Agent 不读取仿真结果图片。

正式完成后成对生成：

- `research/result-039-CDL宽波束CDD与cycling-BLER.md`；
- `research/result-039-CDL宽波束CDD与cycling-BLER-text.md`。

两版必须逐项回答第 1 节问题，区分阶段 1A 筛选证据与阶段 1B evaluation 证据，报告角度变换、波束
平坦度/增益代价、相关性、Sidon 搜索、10%/1% crossing、errors/trials、Wilson/bootstrap 区间、
CE NMSE、异常、capped 点、证据路径和复现命令。未经研究者确认不得更新 `KNOWLEDGE.md` 或 `GOALS.md`，不得创建 Git
checkpoint。

## 12. 执行顺序与正式冻结区

执行顺序：

1. 研究者确认 ASD25、第 6.2/7.3/9 节的当前曲线集、SNR 网格、预算和停止条件；
2. 实现角度变换、新波束、相关矩阵、严格 Sidon 搜索、两种 estimated-CSI 接收机和分析入口；
3. 运行新增单元测试及 plan-037/fixed-CDL/通用 PDSCH 回归；失败即停止；
4. 运行阶段 0 validate/smoke，生成角度谱、波束图输入、相关热力图输入和全部审计；
5. 运行阶段 1A selection prescan，冻结唯一 `SIDON_SELECTED`；
6. 把候选 hash、evaluation seeds 和预算写入冻结清单；候选不再人工挑选；
7. 研究者手动启动 `--stage estimated_confirm`，runner 将展开后的 CDD/cycling SNR 网格与 hash
   写入阶段 1B 报告，并按第 9 节自适应追加完成配对/恢复审计；
8. 本地生成分析表、图和两版 result，交研究者确认；
9. ideal-CSI 正式确认和方法二只在研究者后续明确恢复时执行。
10. 研究者 2026-09-28 已授权直接运行 `--stage aged_mrt_formal`；该阶段不经过 prescan，完成
    固定 26 点后把结果追加到两版 result，保持与原 3 km/h 结果的比较边界。

正式 trial 1 前填写：

- `C_UE2R_DUAL_ASD25` 阶段 1B 展开 YAML 与 SHA-256；
- CDL/ray、阵列、ASD、波束、PDP/协方差 manifest hash；
- 原始 ASD、派生 $s$、变换后 ASD、AoD 和二维角域覆盖率；
- 参考波束平坦度/增益指标、8 波束 Gram/相关指标；
- 严格 Sidon 搜索版本、top-8 候选、`SIDON_SELECTED` delay 和 manifest hash；
- selection/evaluation/bootstrap seeds 与不相交审计；
- 阶段 1B CDD/cycling SNR 网格、batch size、初始和最大 trial 预算；
- validate/smoke/selection audit 路径与 SHA-256；
- 研究者确认日期。
