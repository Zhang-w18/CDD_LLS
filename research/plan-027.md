# plan-027：协方差加权有效矩、Sidon 适用条件与门控 estimated-CSI 验证

> 状态：已由研究者于 2026-07-26 通过“根据 plan27 执行仿真实验任务”确认执行。主基线、有效矩容差 `[0.10,0.25,0.50] dB`、每场景 200,000 个唯一搜索状态、200,000 个共同 outage 样本、1,000 次成对 bootstrap 和 E3 `0.5 dB` 触发门槛按本文冻结；E4 仍需前序结果后的第二次 manifest 确认。
>
> E4 执行补充确认（2026-07-27）：研究者确认 E4 只运行本轮五个 TDL-A 场景；每场景包含 `B0_QC`、`AP_RMS_T1`、`AP_TEPS_T1`、`S0_SIDON`，以及 `AP_T2_CTRL`、`GEO_T1_CTRL`、`MEFF_T2_CAND` 在 10%/1% ideal-CSI outage 上各自胜出代表的去重并集。若同一 family 的两个目标由不同物理候选胜出，两者均运行。五场景分别为 10/9/10/9/9 条物理曲线，共 47 条。链路固定 comb 24，不运行 comb 12/6；候选原样使用 E1 保存的 `delay_grid_coordinates`，连续时延不得量化。场景按 `A5 → A1 → A10 → A30 → A100` 逐个运行，每个场景生成结果并经研究者确认后才运行下一个场景。原 63 条正增益汇总 manifest 保留为历史门控产物，不作为本次链路执行清单。
>
> 系统几何基线增补确认（2026-07-27）：研究者要求加入 `AP_TU_NT`（等差公差 $T_u/N_t$）与 `AP_TALIAS_NT`（等差公差 $\tau_{\rm alias}/N_t$）两条基线，对五个场景补算 10%/1% ideal-CSI outage 和各自工作点 comb-24 CE，并对所有冻结基线/候选补充等效圆周 RMS delay spread 与 99%能量最短圆周支撑宽度。estimated-CSI BLER 只在 A5 增补这两条曲线；原47条 manifest 留存，执行 manifest 增补为49条，其中 A5由9条变为11条。A5沿用既有集成 preflight，不重复 smoke；新增曲线按相同400-trial自适应粗扫、3,000-trial细扫和2 dB扩展上限执行。原 `AP_T1_BASE` 定义不变；另报告五条基线及五者最优的补充比较。

> A10/A30 链路增补确认（2026-07-27）：研究者确认 A5 场景报告并要求接下来运行 A10、A30；两场景均把 `AP_TU_NT` 和 `AP_TALIAS_NT` 纳入五基线集合。执行 manifest 在49条基础上新增4条，变为53条，场景条目数为 `A5/A1/A10/A30/A100 = 11/10/12/11/9`；保留原47条和49条 manifest。此次授权执行顺序为 `A5（已确认）→ A10 → A30`，A10报告完成后才运行A30；A1/A100不在此次授权范围。comb-24、原始时延坐标、400-trial粗扫、3,000-trial细扫、2 dB扩展上限和不重复 smoke 的规则不变。

> A30 高导频密度链路增补确认（2026-07-27）：研究者要求本增补继续归入 `plan-027` 和 `result-027`，不新建实验编号。分别在 comb 12 和 comb 6 下重新生成依赖当前导频周期的 A30 候选，重算各自 ideal-CSI outage 与 matched CE，并运行 estimated-CSI BLER。固定物理构造 `B0_QC`、`S0_SIDON`、`AP_RMS_T1`、`AP_TEPS_T1`、`AP_TU_NT` 保持原时延；`AP_TALIAS_NT` 按当前 $\tau_{\rm alias}/N_t$ 更新；`AP_T2_CTRL`、`GEO_T1_CTRL`、`MEFF_T2_CAND` 按当前 pilot period 重新搜索和冻结。每个导频密度的链路 manifest 取五条基线、`S0_SIDON`，以及三个搜索 family 在10%/1% ideal-outage 上各自胜出代表的去重并集。搜索、outage、trial、扩展和统计口径沿用本 plan 已确认预算；comb-12/6 分别生成独立 manifest、hash、输出目录和结果，禁止覆盖 comb-24 产物。result 除 BLER 曲线外，必须生成10%/1% 的“matched CE NMSE at target SNR–target SNR”散点图；空心/实心含义、颜色、线型和 marker 与 comb-24 图完全一致，并逐项说明所有时延变化。
>
> A100 comb-6 链路增补确认（2026-07-30）：研究者要求在100 ns TDL-A下增加comb-6链路仿真，并交付10%/1%扫描曲线和二维散点图。本增补继续归入 `plan-027` 和 `result-027`，不新建实验编号。候选、基线、公平性、搜索与outage预算、链路trial、扩展、闭合判据和散点语义全部沿用A30 comb-6增补，只把物理RMS delay spread改为100 ns。所有依赖PDP或当前pilot period的量必须按A100/comb-6重新计算并冻结，禁止复用A30搜索候选或覆盖A30产物。
>
> A100 基线增补确认（2026-07-30）：研究者追加一条公差为$T_u/(N_t-1)$的系统几何等差基线。为保持既有`AP_TU_NT`历史结果可比，本轮不改写原$T_u/N_t$基线，而新增`AP_TU_NTM1`；A100 comb-6的链路基线集合由五条增为六条。该构造的首末时延相差一个完整$T_u$，在当前频域循环相位下两路相位完全等价，因此预期pilot rank不超过7；仍须原样运行并报告该退化，不得静默去重或把末端量化为0。

## 1. 研究背景

result-024 在 48 PRB、8 Tx / 1 Rx、平坦独立分支信道、DMRS comb 24 和 known-$\mathbf V$ matched LMMSE 条件下，得到严格 Sidon 时延相对等差 QC 的 10%/1% estimated-CSI BLER 增益 `0.33/0.85 dB`。result-025 和 result-026 表明：

1. TDL-A 5 ns 下严格 Sidon 只保留约 `0.16–0.17 dB` 的 10% BLER 增益，1% 排序不能确定；
2. 只使用 $T_\epsilon$ 和系统几何生成的硬厚 Sidon，在 TDL-A 10 ns、TDL-C 5 ns、TDL-C 10 ns 下相对 T1 等差基线得到 `0.160/0.202/0.174 dB` 的 10% ideal-CSI outage 增益，以及 `0.234/0.291/0.248 dB` 的 1%增益；
3. TDL-A 5 ns 下存在相对本轮等差 B1 outage 更好的非等差候选，但 result-026 没有在相同候选空间中隔离“完整协方差知识”和“解除等差约束”的独立作用；
4. result-026 的当前 $J_{\rm CDD}$ 数值核存在 Gauss–Hermite 欠分辨和非负裁剪造成的假零地板，不能继续用于候选排序；Monte Carlo outage、matched CE 和直接计算的 $M_2^{\rm eff}/M_4^{\rm eff}$ 不受该问题影响；
5. ideal-CSI outage 只能评价分集潜力，不能替代真实 DMRS、信道估计、LLR 和 LDPC 条件下的 estimated-CSI BLER。

本轮因此直接使用完整物理协方差计算有效二阶、四阶相关矩，研究其对 CDD 时延选择的价值；同时系统整理严格 Sidon、厚 Sidon的构造变量和适用条件。最后仅在前序结果完成并经研究者确认后，对入围候选运行 estimated-CSI 链路。

证据：`research/result-024-text.md`、`research/result-025-text.md`、`research/result-026-text.md`、`KNOWLEDGE.md` K15–K19。

## 2. 研究问题

### Q1：完整协方差是否有助于选择 CDD 时延

在 TDL-A RMS delay spread 为 `[1,5,10,30,100] ns` 时，发射端使用完整频域协方差计算 $M_2^{\rm eff}/M_4^{\rm eff}$ 并选择时延，能否找到相对“分别直接以 $T_{\rm RMS}$ 和 $T_\epsilon$ 为等差公差，并在目标上取较优者”的主基线具有更低 10%或1% ideal-CSI outage 的候选？

本轮区分以下比较：

- T2 非等差候选对 T1 展宽感知等差主基线：完整协方差和候选空间放宽的联合价值；
- T2 网格优化等差候选对两条 T1 直接公差等差基线：协方差感知等差设计的总体价值，其中同时包含知识增加和公差优化；
- T2 非等差候选对 T1 厚 Sidon/几何候选：近似相同非等差空间中的知识价值；
- T2 非等差候选对 T2 等差候选：相同完整协方差下解除等差约束的价值。

因此，“优于 T1 等差主基线”可以表述为完整协方差指导设计的实用价值；只有匹配候选空间的 T1/T2 比较也改善时，才进一步表述为协方差知识的独立价值。

### Q2：严格 Sidon和厚 Sidon的适用条件是什么

分析 Sidon/厚 Sidon的构造可行性和性能与以下变量的关系：

- 有效支撑 $T_\epsilon$ 及 $\epsilon$；
- DMRS频域间隔 $S_f$ 和导频折叠周期 $N_p=K/S_f$；
- 发射天线数 $N_t$；
- active subcarrier数 $K$；
- 子载波间隔 $\Delta f$ 和 delay-index物理分辨率；
- 实际 PDP/频域协方差对 $M_2^{\rm eff}/M_4^{\rm eff}$ 和 outage排序的影响。

在相同 TDL-A `[1,5,10,30,100] ns` 场景中，同时评价固定严格 Sidon、硬/软厚 Sidon、T1 等差、T2 等差和 T2 非等差候选。

### Q3：提高导频密度能否修复有 outage 潜力候选的 CE

默认使用 comb 24。只有候选的 ideal-CSI outage 优于主基线、但目标 SNR 附近 matched CE 明显劣化时，才追加 comb 12；comb 12 后仍明显劣化时才追加 comb 6。该阶段只诊断 CE，不把增加导频后的开销变化混入默认 ideal-CSI outage结论。

### Q4：入围候选的 estimated-CSI BLER 是否保持增益

前述筛选、outage和CE全部完成后，先汇总：

1. 本轮具有正 outage增益的候选；
2. result-026 中具有正 outage增益、且不与本轮候选重复的候选；
3. 每个候选的来源、适用场景、主基线、10%/1% outage增益、CE和几何指标。

只有研究者确认该去重清单后，才允许启动 estimated-CSI BLER。未确认时 E4 必须停止，不得自动运行。

## 3. 固定系统、知识条件与物理定义

### 3.1 固定系统

| 项目 | 取值 |
|---|---|
| 有效带宽 | 48 PRB，$K=576$ active subcarriers |
| FFT / CP | 4096 / 288 samples |
| 子载波间隔 | $\Delta f=30$ kHz |
| 时域资源 | 10 OFDM symbols |
| Tx / Rx / layer | 8 / 1 / 1 |
| 默认 DMRS | symbols `[2,7]`，comb 24，offset 0 |
| 默认 pilot / data RE | 48 / 5712 |
| 速度 / 载频 | 0 km/h / 3.5 GHz |
| 调制编码 | 16QAM，MCS 8，码率 553/1024，LDPC 最多 8 次迭代 |
| 信道 | Sionna 1.0.2 TDL-A，RMS delay spread `[1,5,10,30,100] ns` |
| 接收机 | known-$\mathbf V$、known-PDP matched 全带频域 LMMSE |
| DMRS处理 | 两个静态 DMRS等效平均，平均后 LS噪声方差 $N_0/2$ |
| CDD功率 | 每分支恒模 1；数据噪声 $N_0=8/\mathrm{SNR}$ |
| 共同随机数 | 同场景、SNR、trial内候选共享信道、payload、LS noise和data noise |

本轮主实验只使用 NLoS TDL-A。result-026 的 TDL-C 入围候选只在门控 E4 中按其原场景复核，不把 TDL-D/E 共享 specular 实现纳入本轮主结论。

### 3.2 CDD delay定义

$$
V_{k,n}
=
\exp\left(-j2\pi k\Delta f\tau_n\right).
$$

除第6.2节的两条直接物理时延基线外，离散搜索候选使用

$$
\tau_n=\frac{j_n}{K\Delta f},
\qquad
j_n\in\{0,\ldots,575\}.
$$

一个 delay index为

$$
\Delta\tau=\frac1{K\Delta f}=57.870370\ {\rm ns}.
$$

第6.2节的基线直接使用连续的物理时延公差，不得先量化到57.87 ns网格；否则1/5/10/30 ns档位会发生无意义的舍入或重合。所有输出保存：

- `delay_ns`；
- `delay_grid_coordinates = K\Delta f\tau_n`；
- 仅对网格对齐候选保存整数 `delay_indices`，直接物理时延基线将该字段记为不适用；
- $K=576$、$\Delta f=30000$ Hz；
- 实际使用的连续相位定义；
- 适用时的 residue、lift和规范代表。

人工 delay是循环移位/频域线性相位，不作为真实传播时延计入 CP预算。

### 3.3 residue、lift、导频密度

频域导频间隔为 $S_f$ 时，

$$
N_p=\frac K{S_f},
\qquad
j_n=r_n+N_pq_n,
\qquad
r_n=j_n\bmod N_p.
$$

默认 comb 24 时 $N_p=24$；comb 12/6 时分别为 $48/96$。lift不改变当前导频行相位，但会改变数据 RE相位、有效相关矩和 $\mathbf R_{DP}$。

上述整数 residue/lift分解只适用于DFT网格对齐候选。`AP_RMS_T1` 和 `AP_TEPS_T1` 使用连续时延，直接由实际 $\mathbf V$ 构造导频矩阵并计算rank、condition和CE，不强行赋予整数 residue/lift。

在独立同分布 Tx分支下，候选按公共循环移位和天线排列等价规范化；若未来引入分支相关或不同 PDP，不得继续使用该等价化。

### 3.4 发射端和接收端知识

| 标签 | 发射端选择器可用信息 | 接收端 |
|---|---|---|
| T1-support | 只接收 RMS delay spread、$\epsilon$、$T_\epsilon$ 和系统几何；禁止读取实际 TDL tap power或完整协方差 | R2 matched，知道实际 PDP和候选 $\mathbf V$ |
| T2-covariance | 接收完整 TDL-A PDP和频域协方差 | R2 matched，知道实际 PDP和候选 $\mathbf V$ |

虽然评价场景固定为 TDL-A，T1选择器仍必须通过接口隔离，不得从 profile名称或标准 tap表恢复完整协方差。

## 4. 分集指标、outage与CE

物理展宽信道下

$$
\rho_{\rm eff}(k,l)
=
R_{{\rm phy},kl}
\frac{(\mathbf V\mathbf V^H)_{kl}}{N_t}.
$$

直接计算

$$
M_2^{\rm eff}
=
\sum_{k\ne l}|\rho_{\rm eff}(k,l)|^2,
\qquad
M_4^{\rm eff}
=
\sum_{k\ne l}|\rho_{\rm eff}(k,l)|^4.
$$

两者使用完整频域矩阵定义，并以按频差汇总的等价实现加速；Phase 0必须验证两种实现逐候选一致。当前有问题的 $J_{\rm CDD}$ 核不参与筛选、冻结或验收。

单位能量 Gray 16QAM BICM目标谱效率为

$$
R=4\times\frac{553}{1024}=2.1602\ {\rm bit/RE}.
$$

ideal-CSI块互信息和outage为

$$
I_{\rm blk}
=
\frac1K\sum_k
I_{\rm QAM}
\left(
\mathrm{snr}\frac{|g_k|^2}{8}
\right),
\qquad
P_{\rm out}=\Pr[I_{\rm blk}<R].
$$

报告 10%和1% outage所需 SNR。outage依赖实际物理信道统计和候选 $\mathbf V$，但在固定数据 RE和目标速率的 ideal-CSI定义下不依赖信道估计器。

matched CE使用

$$
L_{\rm CE}
=
\frac{
\operatorname{tr}\left[
\mathbf R_{DD}
-
\mathbf R_{DP}
(\mathbf R_{PP}+\sigma_{LS}^2\mathbf I)^{-1}
\mathbf R_{PD}
\right]
}{
\operatorname{tr}(\mathbf R_{DD})
}.
$$

本轮不使用80 dB相对门槛。所有候选先检查 pilot rank、条件数和数值有限性；只有正 outage增益候选及其基线在目标 SNR处计算 matched CE NMSE。

系统几何基线增补后，所有冻结基线和候选还需报告复合 PDP 的等效时延扩展。对物理 TDL 抽头 $(\tau_\ell,p_\ell)$ 和8路等功率 CDD时延 $\tau_n^{\rm CDD}$，构造周期 $T_u$ 上的离散复合 PDP：

$$
\tilde\tau_{\ell,n}
=
(\tau_\ell+\tau_n^{\rm CDD})\bmod T_u,
\qquad
\tilde p_{\ell,n}=\frac{p_\ell}{N_t}.
$$

主字段为圆周 RMS delay spread：

$$
\sigma_{\tau,\rm circ}
=
\min_{\mu\in[0,T_u)}
\sqrt{\sum_{\ell,n}\tilde p_{\ell,n}
d_{\rm circ}^2(\tilde\tau_{\ell,n},\mu)}.
$$

同时报告包含至少99%复合功率的最短圆周区间宽度 $T_{99,\rm circ}$。人工CDD是循环移位，因此不得用未经圆周解缠的 `max(delay)-min(delay)` 替代这两个量。输出字段分别为 `equivalent_circular_rms_delay_spread_ns` 和 `equivalent_99pct_circular_support_width_ns`。

## 5. Sidon与厚 Sidon构造

### 5.1 严格 Sidon

对

$$
s_{a,c}=(j_a+j_c)\bmod K,\qquad 0\le a\le c<N_t,
$$

要求全部 $N_t(N_t+1)/2$ 个无序二元和互不相同。构造时固定 $j_0=0$，通过回溯、随机搜索或逐点加入，只有新产生的二元和与已有二元和均不冲突时才接受新 index，直到得到 $N_t$ 个元素；最后消除公共移位和天线排列重复。

本轮固定严格 Sidon参考为

$$
\mathbf j_{\rm S0}=[0,1,3,7,12,20,30,65].
$$

它在 $K=576$ 下的36个无序二元和全部不同，最小 pair-sum圆周间距为1；comb 24下 residue为 `[0,1,3,7,12,20,6,17]`，最小 fold间距也为1。严格 Sidon构造本身依赖 $K,N_t$，不依赖 $T_\epsilon$ 和导频密度；实际 CE还依赖导频折叠。

### 5.2 厚 Sidon

定义 $T_\epsilon$ 为覆盖 $1-\epsilon$ 平均 PDP能量的最短连续区间，本轮固定 $\epsilon=0.01$。一般保护条件为

$$
d_{\rm circ}(s_{a,c},s_{b,d})
>
2T_\epsilon+T_{\rm margin},
$$

$$
d_{\rm fold}(\tau_a,\tau_b)
>
T_\epsilon+T_{\rm sync}.
$$

本轮主值 $T_{\rm margin}=T_{\rm sync}=0$。换算为 index：

$$
g_\Sigma
=
\left\lfloor
\frac{2T_\epsilon}{\Delta\tau}
\right\rfloor+1,
\qquad
g_{\rm fold}
=
\left\lfloor
\frac{T_\epsilon}{\Delta\tau}
\right\rfloor+1.
$$

构造器只接收 $T_\epsilon,K,N_t,N_p,\Delta f$：

1. 在模 $N_p$ 圆周上选择8个 residue，使最小间距不小于 $g_{\rm fold}$；
2. 为 residue分配整数 lift并形成 $j_n=r_n+N_pq_n$；
3. 规范化去重；
4. 要求36个 pair sum模 $K$ 的最小圆周间距不小于 $g_\Sigma$；
5. 同时满足两项者为硬厚 Sidon；硬约束不可行时，按两类最小间距构造软 Pareto候选，但不得标记为硬厚 Sidon。

必要的圆周装填上界为

$$
g_{\Sigma,\max}
\le
\left\lfloor
\frac{K}{N_t(N_t+1)/2}
\right\rfloor,
\qquad
g_{{\rm fold},\max}
\le
\left\lfloor
\frac{N_p}{N_t}
\right\rfloor.
$$

026在 $K=576,N_t=8,N_p=24$ 下的构造结果为：

| 场景 | $T_\epsilon$ | $g_\Sigma/g_{\rm fold}$ | 200,000几何样本中硬候选数 | 已验证有增益的硬候选 |
|---|---:|---:|---:|---|
| TDL-A 10 ns | 47.966 ns | 2 / 1 | 17,065 | `[0,4,69,198,330,489,496,540]` |
| TDL-A 20 ns | 95.932 ns | 4 / 2 | 164 | 未冻结硬候选outage |
| TDL-A 30 ns | 143.898 ns | 5 / 3 | 3 | 未冻结硬候选outage |
| TDL-C 5 ns | 31.533 ns | 2 / 1 | 17,391 | `[0,3,36,101,320,402,526,567]` |
| TDL-C 10 ns | 63.065 ns | 3 / 2 | 1,515 | `[0,3,42,87,105,165,246,300]` |

这些计数来自随机几何样本去重，不是全空间穷举总数。厚 Sidon直接依赖 $T_\epsilon$、导频密度、$N_t$、$K$ 和 $\Delta f$；具体 PDP功率在 T1构造阶段只通过 $T_\epsilon$ 进入，但会在后续 $M_{2/4}^{\rm eff}$ 和outage评价中改变排序。

按标准化 TDL-A PDP随 RMS delay spread线性缩放，027五个主场景的预期几何门槛为：

| TDL-A RMS | 预期 $T_\epsilon$ | $g_\Sigma/g_{\rm fold}$ | comb-24硬可行性预判 |
|---:|---:|---:|---|
| 1 ns | 4.797 ns | 1 / 1 | 可行；S0本身满足两项几何门槛 |
| 5 ns | 23.983 ns | 1 / 1 | 可行；S0本身满足两项几何门槛 |
| 10 ns | 47.966 ns | 2 / 1 | 可行；S0的pair-sum gap不足 |
| 30 ns | 143.898 ns | 5 / 3 | 可行但接近fold装填极限 |
| 100 ns | 479.660 ns | 17 / 9 | 硬不可行：分别超过16/3装填上界 |

该表是基于result-026 TDL-A有效支撑比例的执行前预期。Phase 0必须从本轮实际展开PDP重新计算；若不一致，以实际确定性计算为准并在进入搜索前暂停核对，不允许静默沿用预期值。A100即使硬厚 Sidon不可行，仍保留最大化pair/fold最小距离的软候选并完整报告。

## 6. 候选与基线

本轮代号和比较角色如下：

| 代号 | 发射端知识 | 候选空间/构造 | 实验角色 | 是否为主基线 |
|---|---|---|---|---|
| `B0_QC` | 固定历史设计 | 原QC等差时延 `[0,9,18,27,36,45,54,63]` | 历史固定参考 | 否 |
| `S0_SIDON` | 固定历史设计 | 严格 Sidon `[0,1,3,7,12,20,30,65]` | 历史固定参考 | 否 |
| `AP_RMS_T1` | T1，只知道 RMS delay spread | 直接取等差公差 $T_{\rm RMS}$：$\tau_n=nT_{\rm RMS}$ | 主基线的组成项 | 是 |
| `AP_TEPS_T1` | T1，只知道 $T_\epsilon$ | 直接取等差公差 $T_\epsilon$：$\tau_n=nT_\epsilon$ | 主基线的组成项 | 是 |
| `AP_TU_NT` | 固定系统几何 | 等差公差 $T_u/N_t$，delay index `[0,72,144,216,288,360,432,504]` | 符号周期基线；comb-24 pilot rank 1 | 否，补充基线 |
| `AP_TALIAS_NT` | 固定系统几何 | 等差公差 $\tau_{\rm alias}/N_t$，delay index `[0,3,6,9,12,15,18,21]` | 导频无混叠周期基线；comb-24 pilot rank 8 | 否，补充基线 |
| `AP_T1_BASE` | T1，只知道 RMS delay spread和 $T_\epsilon$ | `AP_RMS_T1` 与 `AP_TEPS_T1` 在各outage目标上的较优包络 | 所有主增益的分母/比较基线 | **是** |
| `AP_T2_CTRL` | T2，知道完整协方差 | 等差 CDD中按有效矩选择 | 评价协方差感知等差设计的总体作用 | 否，是对照 |
| `GEO_T1_CTRL` | T1，只知道 $T_\epsilon$ | 厚 Sidon及pair/fold几何候选 | 与T2非等差候选做知识等级匹配 | 否，是对照 |
| `MEFF_T2_CAND` | T2，知道完整协方差 | 非等差时延中按 $M_2^{\rm eff}/M_4^{\rm eff}$ 搜索 | 本轮主挑战候选 | 否，是待验证候选 |

其中 `AP` 表示 arithmetic progression等差候选，`GEO` 表示只使用支撑和离散几何，`MEFF` 表示使用协方差加权有效矩，`BASE/CTRL/CAND` 分别表示主基线、对照和挑战候选。

### 6.1 固定参考

| 标签 | delay index |
|---|---|
| B0/QC | `[0,9,18,27,36,45,54,63]` |
| S0/严格 Sidon | `[0,1,3,7,12,20,30,65]` |
| AP_TU_NT | `[0,72,144,216,288,360,432,504]` |
| AP_TALIAS_NT | `[0,3,6,9,12,15,18,21]` |

B0、S0和两条系统几何基线在所有场景计算，但不改写 Q1 的 `AP_T1_BASE`。其中

$$
T_u=\frac1{\Delta f}=33.333333\ \mu{\rm s},
\qquad
\tau_{\rm alias}=\frac1{S_f\Delta f}=1.388889\ \mu{\rm s}.
$$

所以 `AP_TU_NT` 公差为 $4.166667\ \mu{\rm s}$、DFT步长72；`AP_TALIAS_NT` 公差为 $173.611111\ {\rm ns}$、DFT步长3。两者均不得按 outage 或 CE 重新选择公差。

### 6.2 T1展宽感知等差主基线

这里的 RMS delay spread和 $T_\epsilon$ 是等差时延本身的两个公差，不是筛选其他等差候选的硬约束或保护门槛。固定 $\tau_0=0$，直接构造

$$
\tau_n^{\rm RMS}=nT_{\rm RMS},
\qquad
\tau_n^\epsilon=nT_\epsilon,
\qquad n=0,\ldots,7.
$$

两者分别记为 `AP_RMS_T1` 和 `AP_TEPS_T1`。构造时直接使用上述 ns值计算

$$
V_{k,n}=\exp(-j2\pi k\Delta f\tau_n),
$$

不得舍入到整数 delay index，也不枚举、筛选或使用 $M_2/M_4$ 选择这两条基线。两条基线均须完整计算和报告。对每个目标 $p\in\{10\%,1\%\}$，主比较基线 `AP_T1_BASE` 定义为二者所需 SNR的较低包络：

$$
\gamma_p({\tt AP\_T1\_BASE})
=
\min\left[
\gamma_p({\tt AP\_RMS\_T1}),
\gamma_p({\tt AP\_TEPS\_T1})
\right].
$$

10%和1%目标允许由不同组成基线胜出，表格必须注明每个目标的胜出者。后续CE比较使用对应outage目标的胜出组成基线，而不是一个不存在的包络时延向量。两条基线的构造均禁止读取完整 TDL-A协方差。

### 6.3 T2等差参考

T2等差对照在整数DFT网格上枚举

$$
j_n=(j_0+na)\bmod K,\qquad a\in\{1,\ldots,K-1\},
$$

消除公共移位、天线排列和相同无序集合重复后，使用实际 $\mathbf R_{\rm phy}$：

1. 选择 $M_2^{\rm eff}$ 最低候选；
2. 在预定 $M_2^{\rm eff}$ 容差带内选择 $M_4^{\rm eff}$ 最低候选；
3. 保留少量 $M_2^{\rm eff}/M_4^{\rm eff}$ Pareto代表。

该组记为 `AP_T2_CTRL`，是协方差感知的等差对照。它相对 `AP_T1_BASE` 的差异同时包含“读取完整协方差”和“在网格等差公差中优化”两部分，因此只能回答协方差感知等差设计的总体价值，不能解释为纯粹的知识变量因果隔离；它不是主增益基线。

### 6.4 T1非等差候选

每个 $T_\epsilon$ 档位生成：

- 硬厚 Sidon候选；
- pair-sum间距最大候选；
- fold间距最大候选；
- 两类间距的软 Pareto候选；
- 严格 Sidon S0。

选择和冻结阶段不读取实际 PDP或协方差。该组记为 `GEO_T1_CTRL`，是支撑知识下的非等差几何对照，不是主增益基线。

### 6.5 T2协方差加权非等差候选

搜索空间为 $K=576$ 上8个唯一 delay index，并至少要求 comb 24导频 residue互异、pilot rank为8。搜索使用实际 TDL-A协方差和以下起点：

- B0、S0、`AP_T2_CTRL`，以及 `AP_RMS_T1`/`AP_TEPS_T1` 的最近网格投影；投影只作为搜索起点并使用独立ID，不替代或改写连续时延主基线；
- 当前场景的硬/软厚 Sidon代表；
- 随机唯一 residue/lift候选；
- 多个确定性均匀 residue和lift起点。

搜索目标以 $M_2^{\rm eff}$ 为主、$M_4^{\rm eff}$ 为辅。严格字典序通常不会让 $M_4^{\rm eff}$ 生效，因此采用预定 $\epsilon$-constraint：

1. 保留全局最低 $M_2^{\rm eff}$ 候选；
2. 在相对最低 $M_2^{\rm eff}$ 不超过 `[0.10,0.25,0.50] dB` 的三个集合中，分别选择最低 $M_4^{\rm eff}$；
3. 保留全局最低 $M_4^{\rm eff}$ 候选作为消融；
4. 最多保留8个 $M_2^{\rm eff}/M_4^{\rm eff}$ Pareto代表；
5. 公共移位或天线排列等价候选只保留一个。

每场景建议搜索不少于200,000个唯一候选状态；固定随机种子和停止预算，不按中间结果临时增加某个场景预算。该组记为 `MEFF_T2_CAND`，是本轮主挑战候选。

`[0.10,0.25,0.50] dB` 容差和200,000状态预算是草案建议值，必须在正式执行前由研究者确认。

## 7. 实验计划

### Phase 0：增量验证与026复用回执

本轮不重复执行026已经通过、且底层代码未改变的完整校准。以下项目直接复用026回执：

- 网格搜索候选的576点相位、ns和 `delay_indices` 换算；
- B0/S0固定delay、严格 Sidon和基础pair/fold统计；
- TDL-A已有档位的频域协方差构造；
- matched CE闭式公式；
- 16QAM互信息表和outage归一化。

复用条件是相关函数内容和固定配置与026一致，并在 `validation/reuse_receipt.json` 中记录源代码hash、026证据路径和逐项复用状态。任一相关函数发生实质修改时，只重跑受影响的原校准项。

027必须新增验证：

1. $M_2^{\rm eff}/M_4^{\rm eff}$ 全矩阵实现与按频差加速实现逐候选一致；
2. 新增 TDL-A 100 ns协方差 Hermitian、半正定、对角归一，且 $T_\epsilon$ 和硬不可行装填判定正确；
3. `AP_RMS_T1` 和 `AP_TEPS_T1` 分别直接以 $T_{\rm RMS}$ 和 $T_\epsilon$ 为等差公差，连续相位实现不做整数index舍入；二者均不能读取 profile、tap power或完整协方差；`GEO_T1_CTRL` 只能读取 $T_\epsilon$ 和系统几何；
4. `AP_T2_CTRL` 和 `MEFF_T2_CAND` 确实读取对应场景的完整协方差；
5. 新搜索器的公共移位/天线排列去重和同seed重放；
6. comb 12/6下的导频数、fold周期、pilot rank、共同数据 RE集合和CE计算；
7. 成对bootstrap的确定性重放和简单合成样本校准；
8. E4候选去重、manifest SHA-256和未确认manifest拒绝运行。

Phase 0只做确定性、小样本或合成验证，不运行性能仿真。任一新增物理定义、知识隔离、指标一致性或门控测试失败时停止，不进入正式搜索。

### E1：协方差加权有效矩搜索与 ideal-CSI outage

对 TDL-A `[1,5,10,30,100] ns` 分别执行：

1. 计算 $T_\epsilon$ 和实际 $\mathbf R_{\rm phy}$；
2. 构造 B0、S0、`AP_RMS_T1`、`AP_TEPS_T1`、`AP_TU_NT`、`AP_TALIAS_NT`、`AP_T2_CTRL`、`GEO_T1_CTRL`；
3. 搜索 `MEFF_T2_CAND`；
4. 计算并保存全部候选的 $M_2^{\rm eff}$、$M_4^{\rm eff}$、rank、condition、等效圆周 RMS delay spread和99%圆周支撑宽度；仅对网格对齐候选另存pair/fold gap；
5. 按第6节规则冻结候选，冻结后不得按outage结果回改；
6. 系统几何基线增补后，每场景进入outage的候选总数不超过18个；
7. 使用每场景200,000个共同信道样本计算 SNR `0:0.5:24 dB` 的outage；
8. 若10%或1%没有形成双侧 bracket，按0.5 dB步长向需要方向扩展，每次最多2 dB，直到闭合或达到预定总扩展4 dB；
9. 使用对数域插值报告10%/1%目标 SNR，并分别取 `AP_RMS_T1` 与 `AP_TEPS_T1` 中所需SNR较低者形成该目标的 `AP_T1_BASE`；
10. 记录10%和1%目标各自的主基线胜出者；
11. 使用候选共享样本的成对 bootstrap报告相对相应胜出基线的目标 SNR差及95%区间，bootstrap重复数建议1,000。

主要报告：

- `MEFF_T2_CAND - AP_T1_BASE`；
- `AP_T2_CTRL - AP_T1_BASE`；
- `MEFF_T2_CAND - GEO_T1_CTRL`；
- `MEFF_T2_CAND - AP_T2_CTRL`；
- B0和S0固定参考。

增益定义为“基线目标 SNR减候选目标 SNR”，正值表示候选更好。点估计为正但区间跨0时标记为方向未确定；95%区间下界大于0时标记为方向确认；增益不低于0.10 dB同时标记为具有可见幅度。该口径是草案建议，正式执行前确认。

### E2：Sidon/厚 Sidon适用条件分析

复用 E1同一批场景、候选和outage样本，形成以下地图：

1. $T_\epsilon/\Delta\tau$；
2. 所需和实际 pair-sum gap；
3. 所需和实际 fold gap；
4. pair packing比例
   $$
   \eta_\Sigma=
   \frac{N_t(N_t+1)g_\Sigma}{2K};
   $$
5. fold packing比例
   $$
   \eta_{\rm fold}=
   \frac{N_tg_{\rm fold}}{N_p};
   $$
6. 硬厚 Sidon可行候选数和软约束最优距离；
7. S0、硬/软厚 Sidon相对 `AP_T1_BASE` 的10%/1% outage增益；
8. $M_2^{\rm eff}/M_4^{\rm eff}$、outage和CE排序是否一致。

严格 Sidon的“适用”只在相应场景下outage或后续BLER相对主基线改善时成立；满足加性构造本身不等于性能成立。厚 Sidon同样区分“几何硬可行”“ideal-CSI outage有利”和“estimated-CSI BLER有利”三个层级。

### E3：目标 SNR处CE与条件触发的导频密度

对 E1中10%或1% outage点估计优于 `AP_T1_BASE` 的候选，以及其基线、B0/S0和两条系统几何基线：

1. 在候选自己的10%/1%目标 SNR计算matched CE NMSE；
2. 在对应基线的10%/1%目标 SNR再次计算双方CE，区分工作点差异和结构差异；
3. 保存 pilot rank、condition number、最小奇异值、$\mathbf R_{PP}$和$\mathbf R_{DP}$诊断；
4. 默认comb 24下，若候选在任一目标的CE相对基线劣化超过 `0.5 dB`，触发comb 12；
5. comb 12后仍劣化超过 `0.5 dB`，触发comb 6；
6. 候选 delay不重新优化，只改变导频集合和matched估计器；
7. 导频集合必须嵌套，并在各密度共同的数据 RE交集上报告CE；
8. 主诊断固定每导频 RE功率，因此更高密度会增加总导频能量和开销；结果只解释为“增加观测后的CE可恢复性”，不冒充同开销系统增益。

`0.5 dB` 触发门槛和每导频 RE功率固定策略是草案建议，执行前确认。若需要比较最终链路净收益，必须另行固定导频总能量、payload和有效谱效率，不在E3后验混入。

### E4：研究者确认后才开启的 estimated-CSI BLER

#### E4.1 候选清单

E1–E3完成后先生成 `e4_candidate_gate.csv` 和 `e4_candidate_gate.md`，字段至少包括：

- 来源实验、场景、原candidate ID和family；
- `delay_ns`、`delay_grid_coordinates`、适用时的 `delay_indices`、residue、lift和规范key；
- 主基线及其delay；
- 10%/1% outage目标和增益；
- 增益95%区间；
- comb 24目标 SNR处CE；
- 可用时的comb 12/6 CE；
- strict/thick Sidon标志；
- $M_2^{\rm eff}$、$M_4^{\rm eff}$；
- 与其他候选的重复关系；
- 建议运行的场景和导频配置。

本轮候选按“任一目标outage点估计为正”全部列出，同时单独标记区间是否确认。result-026从其已冻结outage表重新读取所有正增益非基线候选；预计至少覆盖：

- E1：`E1_0090`、`E1_0006`、`E1_0005`、`E1_0003`、`E1_0079`、`E1_0061`、`E1_0028`；
- E2：`E2-A10_0311`、`E2-A20_0216`、`E2-A30_0053`、`E2-C5_0132`、`E2-C10_0098`。

实际清单必须从
`outputs/experiment026_cdd_design/20260724_main/e1_outage/outage_targets.csv`
和
`outputs/experiment026_cdd_design/20260724_main/e2_outage/outage_targets.csv`
重新确定，不以本段手写列表替代数据读取。

去重规则：

1. 先按公共移位和天线排列规范化后的delay key去重；
2. 若result-026候选与本轮候选相同，本轮条目保留并合并provenance，不重复运行；
3. 同一delay在不同物理场景下不是重复实验，可保留多个场景行；
4. B0、S0和各场景主基线作为比较参考，不计入“新候选”数量。

生成规范化 JSON manifest及 SHA-256。研究者必须确认确切manifest、场景、导频配置和预算；正式BLER入口必须要求显式提供已确认manifest及其hash，不匹配时拒绝运行。

#### E4.2 门控条件

以下全部满足后才允许启动：

1. E1–E3正式输出完成；
2. 两种有效矩实现、outage目标和CE均通过验证；
3. 候选已去重；
4. 研究者书面确认 `e4_candidate_gate` manifest；
5. 研究者确认粗扫/加密trial预算和是否纳入高密度导频变体。

本plan确认不自动等价于E4候选确认；E4需要前序结果之后的第二次确认。

#### E4.3 链路流程

对确认manifest中的每个“场景×候选×基线”：

1. 运行20 paired trials/点的smoke，只验证配置、协方差、共同随机数、CE、LLR、LDPC和输出；
2. 默认comb 24运行 SNR `13:0.5:20 dB`、每点400 paired trials的粗扫；
3. 任一候选未跨越10%时，以0.5 dB向相邻方向扩展，每次最多1 dB，总扩展不超过2 dB；
4. 根据候选和基线各自10%跨越区间的离散并集，以0.25 dB加密，每点3,000 paired trials；
5. 使用预定区间全部点的二项logit拟合10%目标；保存delta-method区间、配对四格计数和McNemar精确检验；
6. 只有预计目标附近每候选累计错误块不少于30时才运行1%加密，否则报告未闭合或先导结果；
7. 正式报告 estimated-CSI BLER、CE、错误块数、trial数、目标SNR和保守95%增益区间。

上述链路网格和trial数沿用result-026量级，仍需在E4第二次确认时结合最终候选数核定。未经确认不得因预算方便擅自删候选。

#### E4.4 2026-07-27 确认后的自适应执行口径

本次链路执行另生成 `e4_link_gate/`，不得覆盖原 `e4_gate/`。执行 manifest 使用原生数值数组保存 `delay_grid_coordinates`，并冻结以下配置：

- DMRS 只使用 comb 24；
- smoke 为首场景 A5 的一次集成 preflight：每候选一个预定 SNR 点、20 trials；后续场景不重复 smoke；
- prescan 每候选每个目标以 `ideal-CSI outage SNR + 5 dB` 为初始中心，按 0.5 dB 网格取中心及相邻点，每点 400 trials；
- 若 10%或1%未形成双侧 bracket，按0.5 dB逐次向所需方向扩展，每个目标相对初始点的总扩展不超过2 dB；
- 10%和1%分别在各候选自己的相邻粗扫 bracket 内以0.25 dB网格加密，每点3,000 trials；正式细扫完成初始三点后，再用正式 BLER 检查是否实际夹住目标，若粗扫小样本导致区间偏移，则按候选、按0.25 dB向所需方向追加，仍受每目标最多2 dB扩展约束，禁止直接外推；不同候选不为形成全局凸包而补算无关 SNR；
- 任何候选在目标附近累计错误块少于30时只报告样本不足，不作确定性排序；目标没有双侧 bracket 时不外推；
- 同一 `scenario_id/SNR/trial` 的候选共享信道、payload、导频噪声和数据噪声；候选分批执行时仍由稳定 trial key 恢复同一随机样本；
- 全部 SNR 点保存逐 trial TB error flag，以便候选与三条基线在共同 SNR 点形成成对四格计数；
- 每个场景完成后输出 BLER 曲线、目标 SNR 表、相对三条基线及 estimated-CSI 三者最优基线的比较表、CE、成对计数、展开配置和场景报告，然后停止等待研究者确认。

原 manifest 包含47个场景×物理候选条目；首轮系统几何基线增补 manifest 为49条，只在 A5 增加 `AP_TU_NT` 和 `AP_TALIAS_NT`。A10/A30增补后，当前执行 manifest 必须包含恰好53个条目，场景条目数为 `A5/A1/A10/A30/A100 = 11/10/12/11/9`，并把此次授权执行顺序冻结为 `A5 → A10 → A30`。同时保存原47条和49条manifest、增补来源、相位分母576、首 active subcarrier 相位参考、trial预算和场景顺序；生成的新 hash及对应回执必须写入输出。

### E5：A30 comb-12/6 高导频密度增补

#### E5.1 研究问题与可证伪假设

本增补回答：A30 下 comb-24 未闭合、但具有正 ideal-CSI outage 增益的结构化 CDD family，增加导频密度后能否形成10%/1% estimated-CSI BLER 双侧 bracket，并在与相同导频密度五基线的公平比较中保留增益。

预定假设为：

1. comb-12 或 comb-6 会降低部分候选的 matched CE NMSE 和高 SNR 误差地板，使原 comb-24 未闭合目标闭合；
2. CE 改善不保证净 BLER 增益，因为同导频密度基线也会改善，且导频开销与总导频能量同时改变；
3. 若候选在 comb-6 下仍未闭合，或闭合后不优于同密度五基线最优者，则该导频密度不能证明其 A30 实用增益。

#### E5.2 固定系统与公平性

除 DMRS 频域间隔外，沿用第3.1节 A30 系统：48 PRB、$K=576$、30 kHz、8 Tx / 1 Rx / 1 layer、TDL-A 30 ns、零速度、16QAM、MCS 8、LDPC 最多8次迭代、known-$\mathbf V$ 和 known-PDP matched 全带频域 LMMSE。相位仍为

$$
V_{k,n}=\exp\left(-j2\pi k\frac{j_n}{576}\right),
$$

相位分母固定576，参考为首个 active subcarrier，连续时延不得量化。两个增补条件为：

| 条件 | DMRS spacing | 每DMRS symbol导频数 | 总pilot RE | data RE | pilot period $N_p=K/S_f$ |
|---|---:|---:|---:|---:|---:|
| comb-12 | 12 | 48 | 96 | 5664 | 48 |
| comb-6 | 6 | 96 | 192 | 5568 | 96 |

主执行保持每个 pilot RE 功率不变，因此 comb-12/6 相对 comb-24 的总导频能量分别为2倍和4倍；本增补的绝对改善只能解释为“更多观测与更多导频能量下的链路表现”。同一密度内所有候选和基线使用相同 DMRS、MCS、数据 RE、trial和共同随机数，允许比较相对 BLER。不得把不同密度间的绝对 BLER差异表述为同开销净收益。

#### E5.3 随导频密度变化的时延规则

每个密度独立生成候选，禁止先看 BLER 再改变时延：

1. `B0_QC=[0,9,18,27,36,45,54,63]`、`S0_SIDON=[0,1,3,7,12,20,30,65]` 和 `AP_TU_NT=[0,72,144,216,288,360,432,504]` 保持固定；
2. `AP_RMS_T1` 仍为连续物理公差30 ns，`AP_TEPS_T1` 仍为连续物理公差本轮 A30 的 $T_\epsilon=143.898$ ns，二者不随导频密度变化；
3. `AP_TALIAS_NT` 使用当前 $\tau_{\rm alias}=1/(S_f\Delta f)$，所以 comb-12 的 delay index 为 `[0,6,12,18,24,30,36,42]`，comb-6 为 `[0,12,24,36,48,60,72,84]`；
4. `AP_T2_CTRL` 在整数等差集合中重新施加“模当前 $N_p$ residue互异、pilot rank 8”的可行性约束，再按第6.3节相同 $M_2^{\rm eff}/M_4^{\rm eff}$ 规则冻结2个代表；
5. `GEO_T1_CTRL` 使用当前 $N_p$ 重算 fold gap、硬/软厚 Sidon 几何和200,000个确定性状态，按第6.4节冻结最多3个代表；
6. `MEFF_T2_CAND` 使用当前 $N_p$ 的唯一 residue约束、相同起点、种子派生规则、200,000个唯一状态和 `[0.10,0.25,0.50] dB` $M_2$ 容差，按第6.5节冻结最多7个代表。

每个密度必须输出所有冻结候选的 `delay_grid_coordinates`、物理 ns、residue、lift、pilot rank、condition number、pair/fold gap和选择规则。result 必须给出 comb-24/12/6 对照表；搜索 family 的编号只在各自密度内有效，不得把相同后缀误写为同一个物理候选。

#### E5.4 ideal-outage、链路 manifest 与正式运行

每个密度分别使用200,000个共同 A30 信道样本计算当前数据 RE集合上的 ideal-CSI 10%/1% outage，SNR网格、扩展、插值和1,000次成对 bootstrap沿用E1。链路 manifest 包含：

- 五基线：`B0_QC`、`AP_RMS_T1`、`AP_TEPS_T1`、`AP_TU_NT`、`AP_TALIAS_NT`；
- 固定参考 `S0_SIDON`；
- `AP_T2_CTRL`、`GEO_T1_CTRL`、`MEFF_T2_CAND` 在10%/1% outage各自胜出代表的去重并集。

comb-12和comb-6分别生成 manifest及 SHA-256；本段研究者确认即为这两个按规则确定的 manifest 执行授权，但实际运行前仍必须由代码核对 schema、密度、A30范围、候选数、时延字段和hash。若同一 family 两个目标由不同候选胜出，两者均运行。

每个密度执行一次每候选20 trials的集成 smoke。正式链路沿用E4.4：

- prescan以各候选当前密度 ideal-outage SNR加5 dB为中心，0.5 dB步长，每点400 paired trials；
- 每目标最多向所需方向扩展2 dB；
- 形成粗扫 bracket 后以0.25 dB、每点3,000 paired trials正式加密；
- 正式点偏移时继续按0.25 dB追加，但仍受2 dB扩展上限；
- 没有正式双侧 bracket时不外推；目标附近累计错误块少于30时标记样本不足。

#### E5.5 输出、图与验收

输出写入 `outputs/experiment027_meff_sidon/<run_id>/e5_dense_dmrs/comb12/` 和 `comb6/`，至少包含 search、outage、manifest、smoke、prescan、refine、final、展开配置、日志、逐trial错误标志和机器可读表。

每个密度必须输出：

1. 全候选 prescan BLER曲线和正式闭合目标 BLER曲线；
2. 10%和1%的“matched CE NMSE at target SNR–target SNR”散点图；
3. 空心点表示当前密度的 ideal-outage目标 SNR及该点matched CE，实心点表示正式闭合的 estimated-CSI BLER目标 SNR及该点matched CE；未闭合目标只画空心点，不外推；
4. 同一候选的空心/实心点用线连接；颜色、线型和 marker 必须调用E4相同 `curve_style` 映射；
5. 图对应 CSV，包含 density、candidate、family、target、point kind、SNR、CE NMSE、BLER区间和闭合状态；
6. comb-24/12/6 精确时延变化表及变化原因。

验收分别报告：

- 各密度10%/1%闭合曲线数；
- 各 family 相对同密度五基线最优者的目标 SNR点估计和95%区间；
- 提高密度是否只实现闭合，还是保留正 BLER增益；
- 导频能量与开销变化造成的解释边界；
- 未闭合、样本不足、rank不足和 CE地板异常。

### E6：A100 comb-6 链路增补

#### E6.1 研究问题、候选与公平性

本增补回答：在TDL-A 100 ns、comb-6下，A100中具有ideal-CSI outage潜力的结构化CDD family能否形成10%/1% estimated-CSI BLER双侧bracket，并相对同一comb-6六基线最优者保留增益。除RMS delay spread固定为100 ns且只执行comb-6外，系统、功率、MCS、接收机、数据RE、共同随机数和公平性约束与E5一致。

候选在看到BLER前冻结：

1. 固定`B0_QC`、`S0_SIDON`和`AP_TU_NT`物理时延；
2. `AP_RMS_T1`使用连续100 ns公差，`AP_TEPS_T1`使用A100实际$T_\epsilon$连续公差；
3. `AP_TALIAS_NT`使用comb-6的当前$\tau_{\rm alias}/N_t$，delay index为`[0,12,24,36,48,60,72,84]`；
4. `AP_T2_CTRL`、`GEO_T1_CTRL`和`MEFF_T2_CAND`按A100完整协方差、comb-6的$N_p=96$和E5相同规则分别重新搜索；
5. 新增`AP_TU_NTM1`，公差为$T_u/(N_t-1)$，连续物理时延为

   $$
   \tau_n=n\frac{T_u}{N_t-1},\qquad n=0,\ldots,N_t-1.
   $$

   在$N_t=8$时，`delay_grid_coordinates=[0,576/7,\ldots,576]`；首末分支在当前循环相位下等价，必须保存连续坐标和实际pilot rank，不得按模576去重；
6. 链路manifest包含六基线、`S0_SIDON`，以及三个搜索family在10%/1% ideal-outage上的胜出代表去重并集。

A100的硬厚Sidon在comb-24下已由装填上界判为不可行；comb-6扩大fold周期但不改变pair-sum装填上界，仍不得把软几何胜出者标记为硬厚Sidon。

#### E6.2 统计预算、输出与验收

搜索使用200,000个唯一T2状态和200,000个确定性几何状态；ideal-CSI outage使用200,000个共同A100信道样本和1,000次成对bootstrap。链路先执行每候选20 trials smoke，再按各候选ideal-outage目标加5 dB初始化，以0.5 dB、每点400 trials粗扫；10%和1%分别以0.25 dB、每点3,000 trials细扫，每目标最多向所需方向扩展2 dB。目标未形成正式双侧bracket时不外推，目标附近累计错误块少于30时标记样本不足。

输出写入`outputs/experiment027_meff_sidon/<run_id>/e6_a100_dense_dmrs/comb6/`，并生成独立manifest、SHA-256、审批回执、展开配置、逐trial错误标志和机器可读表。必须交付：

1. 全候选prescan BLER曲线；
2. 正式闭合目标的10%/1% BLER扫描曲线；
3. 10%和1%各一张“matched CE NMSE at target SNR–target SNR”二维散点图及对应CSV；
4. A100 comb-24/comb-6精确时延变化表；
5. 各family相对同密度六基线最优者的目标SNR点估计、95%区间、闭合数、未闭合与样本不足说明。

二维图继续使用E5语义：空心点为ideal-outage目标及该SNR的matched CE，实心点为正式闭合的estimated-CSI BLER目标及该SNR的matched CE，实心点横向误差棒为目标SNR的95%区间；未闭合目标只有空心点，不外推。颜色、线型和marker必须复用E4/E5的`curve_style`映射。

## 8. 统计、公平性与解释边界

1. 同场景候选共享信道和随机数；
2. 搜索候选在看到outage前冻结，BLER候选在看到BLER前冻结；
3. T1选择器不得读取完整协方差；
4. T2候选按有效矩选择，不使用outage反向优化；
5. 10%和1%目标必须形成双侧bracket；未闭合时不强行外推正式结论；
6. ideal-CSI outage、CE和estimated-CSI BLER分别命名，不互相替代；
7. 多个场景中仅一个正点估计不能自动表述为稳健规律；必须同时报告全部预定场景和区间；
8. 高导频密度CE诊断改变导频能量和开销，除非另行做同有效速率链路，不表述为免费性能改善；
9. TDL-A为零均值NLoS模型；本轮不把结论推广到TDL-D/E、CDL、空间相关或几何LoS；
10. 不使用result-026当前有问题的 $J_{\rm CDD}$ 核。

## 9. 代码修改、测试与输出

### 9.1 预计代码

- `cdd_lls/design/cdd_metrics.py`：补充并验证直接 $M_2^{\rm eff}/M_4^{\rm eff}$ 接口；不得复用未修复JCDD核；
- `cdd_lls/design/cdd_search.py`：T1直接物理公差等差基线、严格/厚 Sidon、T2有效矩搜索和规范去重；
- `tools/run_plan027_meff_design.py`：Phase 0、E1、E2、outage和候选冻结；
- `tools/run_plan027_ce_density.py`：E3目标SNR CE和条件触发导频密度；
- `tools/build_plan027_bler_gate.py`：E4去重清单、manifest和hash；
- `tools/run_plan027_bler.py`：只接受已确认manifest的多候选链路入口；按场景执行一次preflight、自适应prescan、10%/1% refine、断点续跑、逐trial错误标志、目标拟合、比较表和曲线；
- `tools/analyze_plan027.py`：目标SNR、bootstrap、适用条件地图、候选清单和最终汇总；
- `tests/test_plan027_meff.py`、`tests/test_plan027_search.py`、`tests/test_plan027_gate.py`：指标、构造、知识隔离、去重和门控测试。

若现有模块可直接复用，可减少文件；不得把有效矩、规范化或门控逻辑复制到多个脚本。

### 9.2 输出目录

```text
outputs/experiment027_meff_sidon/<run_id>/
```

建议结构：

```text
validation/
e1_search/
e1_outage/
e2_sidon_map/
e3_ce_density/
e4_gate/
e4_link_gate/
e4_link/<scenario>/{smoke,prescan,refine_10pct,refine_1pct,final}/
e4_smoke/
e4_prescan/
e4_refine_10pct/
e4_refine_1pct/
final/
```

必需输出：

- `resolved_experiment.json`、`commands.json`、`environment.json`、日志和代码/工作区标识；
- `e1_all_candidates.csv`、`e1_frozen_candidates.csv`、`e1_meff_pareto.csv`；
- `e1_outage_curves.csv`、`e1_outage_targets.csv`、`e1_paired_bootstrap.csv`；
- `e2_support_geometry.csv`、`e2_sidon_candidates.csv`、`e2_applicability.csv`；
- `e3_ce_at_targets.csv`、`e3_density_trigger.csv`；
- `e1_search/system_baseline_addendum.json`、两条系统几何基线在 `e1_outage_targets.csv` 中的10%/1%目标、各自工作点 CE；
- `final/all_candidate_equivalent_delay_spread.csv`：全部冻结基线/候选的圆周 RMS与99%圆周支撑宽度；
- `e4_candidate_gate.csv`、`e4_candidate_gate.md`、`e4_candidate_manifest.json`和SHA-256；
- 门控开启后的BLER曲线、逐点错误计数、配对四格计数、CE和目标SNR；
- `final/final_summary.json`和图表源数据。

CSV中的ns、SNR、NMSE、增益和概率字段必须带明确单位后缀。

## 10. 执行顺序与停止条件

1. 研究者确认本草案中的主基线、有效矩容差、搜索预算、outage统计和E3触发门槛；
2. 实现并通过Phase 0；
3. 运行E1搜索并冻结候选；
4. 运行E1 10%/1% outage；
5. 完成E2 Sidon适用条件分析；
6. 对正outage候选运行E3目标SNR CE和条件触发导频密度；
7. 生成E4去重候选manifest；
8. 暂停，等待研究者第二次确认；
9. 只有manifest确认后运行E4 smoke、粗扫和加密；
10. 生成 `research/result-027.md` 与 `research/result-027-text.md`；
11. 研究者确认result后再更新 `GOALS.md`、`KNOWLEDGE.md` 和研究索引。

以下情况必须停止对应阶段并保留已有输出：

- 物理时延、相位分母、功率归一化或导频数量不一致；
- T1选择器读取了实际协方差；
- $M_2^{\rm eff}/M_4^{\rm eff}$ 两种实现不一致；
- 搜索、冻结或共同随机数无法重放；
- outage/BLER目标未闭合；
- 需要改变基线、候选集合、阈值或统计预算；
- E4 manifest未确认或hash不匹配。

## 11. result-027必须回答

两版result必须自包含并回答：

1. 每个 TDL-A delay spread的 $T_\epsilon$、协方差和厚 Sidon几何要求；
2. `AP_RMS_T1`、`AP_TEPS_T1`、各outage目标的 `AP_T1_BASE` 胜出者、`AP_T2_CTRL`、`GEO_T1_CTRL`、`MEFF_T2_CAND` 的确切构造和delay；
3. 哪些候选由最低 $M_2^{\rm eff}$、受限最低 $M_4^{\rm eff}$ 或Pareto规则选出；
4. 五个场景完整的10%/1% outage目标、增益和区间；
5. 完整协方差指导设计相对T1等差基线是否有实用价值；
6. 匹配候选空间后能否观察到协方差知识的独立价值；
7. 严格 Sidon和厚 Sidon在哪些几何、outage、CE层级成立或失效；
8. 导频密度提高是否改善CE，改善中包含什么能量和开销变化；
9. E4候选如何从027和026汇总、去重，最终确认manifest是什么；
10. 若E4获准运行，各候选的10%和可用1% estimated-CSI BLER是否保持增益；
11. 全部未闭合目标、负结果、数值异常和适用范围；
12. 哪些结论可以进入 `KNOWLEDGE.md`，哪些仍需后续验证。
13. 两条系统几何基线的10%/1% ideal-CSI outage、各自工作点 CE、pilot rank，以及 A5 estimated-CSI BLER是否闭合；
14. 所有冻结基线和候选在各物理场景下的等效圆周 RMS delay spread与99%圆周支撑宽度。
15. A100 comb-6六基线与重搜候选的精确时延、10%/1%扫描、闭合状态、二维散点、区间，以及`AP_TU_NTM1`的循环端点退化。

## 12. 建议复现入口

实际实现后必须把精确命令和展开参数保存到输出目录。建议接口：

```powershell
& D:\venvs\cdd-s102\Scripts\python.exe -m unittest tests.test_plan027_meff tests.test_plan027_search tests.test_plan027_gate
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan027_meff_design.py --stage validate --run-id <run_id>
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan027_meff_design.py --stage search --run-id <run_id>
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan027_meff_design.py --stage add-system-baselines --run-id <run_id>
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan027_meff_design.py --stage outage --run-id <run_id>
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan027_ce_density.py --stage ce --run-id <run_id>
& D:\venvs\cdd-s102\Scripts\python.exe tools\build_plan027_bler_gate.py --run-id <run_id>
& D:\venvs\cdd-s102\Scripts\python.exe tools\build_plan027_bler_gate.py --run-id <run_id> --scope envelope-link
```

研究者确认E4 manifest后才允许：

```powershell
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan027_bler.py --stage smoke --scenario A5 --run-id <run_id> --approved-manifest <path> --approved-sha256 <sha256> --approval-receipt <approval.json>
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan027_bler.py --stage prescan --scenario A5 --run-id <run_id> --approved-manifest <path> --approved-sha256 <sha256> --approval-receipt <approval.json>
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan027_bler.py --stage refine-10pct --scenario A5 --run-id <run_id> --approved-manifest <path> --approved-sha256 <sha256> --approval-receipt <approval.json>
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan027_bler.py --stage refine-1pct --scenario A5 --run-id <run_id> --approved-manifest <path> --approved-sha256 <sha256> --approval-receipt <approval.json>
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan027_bler.py --stage analyze --scenario A5 --run-id <run_id> --approved-manifest <path> --approved-sha256 <sha256> --approval-receipt <approval.json>
& D:\venvs\cdd-s102\Scripts\python.exe tools\analyze_plan027.py --run-id <run_id>
```

E6 A100 comb-6使用：

```powershell
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan027_dense_dmrs.py --stage search --scenario A100 --spacing 6 --run-id <run_id> --search-states 200000 --geometry-states 200000
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan027_dense_dmrs.py --stage outage --scenario A100 --spacing 6 --run-id <run_id> --outage-samples 200000 --bootstrap-repeats 1000
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan027_dense_dmrs.py --stage manifest --scenario A100 --spacing 6 --run-id <run_id>
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan027_bler.py --stage <smoke|prescan|refine-10pct|refine-1pct|analyze> --scenario A100 --run-id <run_id> --approved-manifest <path> --approved-sha256 <sha256> --approval-receipt <approval.json>
& D:\venvs\cdd-s102\Scripts\python.exe tools\run_plan027_dense_dmrs.py --stage analyze --scenario A100 --spacing 6 --run-id <run_id>
```
