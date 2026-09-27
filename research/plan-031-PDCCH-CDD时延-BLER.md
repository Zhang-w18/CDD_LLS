# plan-031：PDCCH AL2/4/8 CDD 时延 BLER

目标：比较 `B0_QC`、`S0_SIDON`、`AP_RMS_T1`、`AP_TEPS_T1`、`GEO_T1_CTRL`（硬厚 Sidon 不可行时的软几何控制）、`A100_PRG_DFT8_6RB`（precoder cycling）、`A100_SMALL_CDD_QSTEP0P25_MATCHED_CDD` 和 `A100_SMALL_CDD_QSTEP0P25_TRANSPARENT_CDD`。后两个候选使用完全相同的小时延 CDD 发射信号，分别表示 non-transparent 与 transparent 接收；仅输出 estimated-CSI DCI BLER。
注：GEO_T1_CTRL 的全称可以理解为“仅使用 T1 信息的几何控制组”：
- GEO：只优化时延的几何关系。
- T1：发射端只知道 \(T_\epsilon\)，不知道完整 PDP/协方差。
- CTRL：作为对照组，不代表最优方案。
- 它不是硬厚 Sidon，而是硬约束不可行时，从相同准则导出的软折中候选。

配置：固定 A=41、AL=2/4/8、QPSK、8Tx/1Rx、48-RB/1-symbol CORESET、30 kHz、non-interleaved、first CCE 0、L=6、static TDL-A 100 ns。各 AL 实占 12/24/48 RB，$K=144/288/576$，$E=216/432/864$，DMRS 占比 1/4。均单位范数、data/DMRS 同预编码。

`A=41` 沿用 plan-029 的代表性短 DCI operating point，用于固定 payload 做横向比较；它不是 3GPP 规定的统一 DCI 长度，实际 DCI 大小由 DCI format、BWP 和相关配置共同决定（3GPP TS 38.212 §7.3.1）。

计划取值：各 AL 均按当前有效带宽定义

$$
V_{k,n}=\frac1{\sqrt8}\exp(-j2\pi k j_n/K),
\qquad q_K=\frac1{K\Delta f}.
$$

PDCCH DMRS 为 $k=12b+\{1,5,9\}$，等效 $S_f=4$。因此 $N_p=K/4=36/72/144$，三个 AL 的 DMRS 折叠周期均为 $\tau_{\rm alias}=1/(4\Delta f)=8.333333\ \mu\mathrm{s}$；这只是导频上的模周期，不是人工时延的硬上限。允许规则中的 lift $j=r+N_pq$；规范化人工时延只要求 $0\le j<K$，对应全带相位周期 $1/\Delta f=33.333333\ \mu\mathrm{s}$。

| AL | RB | $K$ | $q_K=1/(K\Delta f)$ | $N_p$ |
|---:|---:|---:|---:|---:|
| 2 | 12 | 144 | 231.481481 ns | 36 |
| 4 | 24 | 288 | 115.740741 ns | 72 |
| 8 | 48 | 576 | 57.870370 ns | 144 |

表中 $\mathbf j=\boldsymbol\tau/q_K=K\Delta f\boldsymbol\tau$。规则分别作用于模 $K$ 的 pair-sum 几何和模 $N_p$ 的 DMRS fold 几何；改变 $K$ 后模数也改变，所以某个 $K$ 下满足规则不推出另一个 $K$ 下自动满足，必须重新生成或核验。所有可设计方案要求任意两分支折叠后不同：连续坐标满足 $(j_a-j_b)\bmod N_p\ne0$，整数坐标等价于 $j_a\bmod N_p$ 互异。`AP_RMS_T1` 也做相同检查，但其物理值由 RMS 定义，若失败只能标记失效而不回调时延；本轮该检查通过。

`B0_QC` 不重新搜索，复用 48-RB 历史物理时延；AL2/4 只换算为各自的连续 $j$，从而消除原先直接复用整数 $j$ 导致的 AL2 fold 重合。两个小时延 CDD 候选同样固定历史物理小时延，发射端逐 RE 的 `V` 完全相同；差别只允许来自接收机协方差。`AP_RMS_T1`、`AP_TEPS_T1` 直接固定物理公差；`S0_SIDON` 和 `GEO_T1_CTRL` 则在每个 $K$ 的模环上单独核验或生成。连续方案不得量化到整数栅格。

### AL2：12 RB，$K=144$，$N_p=36$

| 方案 | $\mathbf j=\boldsymbol\tau/q_K$ | $\boldsymbol\tau$ / ns | 最小 fold gap / $q_K$ |
|---|---|---|---:|
| `B0_QC` | `[0,2.25,4.5,6.75,9,11.25,13.5,15.75]` | `[0,520.833,1041.667,1562.5,2083.333,2604.167,3125,3645.833]` | 2.25 |
| `S0_SIDON` | `[0,1,3,7,12,20,30,65]` | `[0,231.481,694.444,1620.370,2777.778,4629.630,6944.444,15046.296]` | 1 |
| `AP_RMS_T1` | `[0,.432,.864,1.296,1.728,2.160,2.592,3.024]` | `[0,100,200,300,400,500,600,700]` | 0.432 |
| `AP_TEPS_T1` | `[0,2.072131,4.144262,6.216394,8.288525,10.360656,12.432787,14.504918]` | `[0,479.660,959.320,1438.980,1918.640,2398.300,2877.960,3357.620]` | 2.072131 |
| `GEO_T1_CTRL` | `[0,3,7,12,22,101,127,133]` | `[0,694.444,1620.370,2777.778,5092.593,23379.630,29398.148,30787.037]` | 3 |
| `A100_SMALL_CDD_QSTEP0P25_MATCHED_CDD` | `[0,.0625,.125,.1875,.25,.3125,.375,.4375]` | `[0,14.468,28.935,43.403,57.870,72.338,86.806,101.273]` | 0.0625 |
| `A100_SMALL_CDD_QSTEP0P25_TRANSPARENT_CDD` | `[0,.0625,.125,.1875,.25,.3125,.375,.4375]` | `[0,14.468,28.935,43.403,57.870,72.338,86.806,101.273]` | 0.0625 |

### AL4：24 RB，$K=288$，$N_p=72$

| 方案 | $\mathbf j=\boldsymbol\tau/q_K$ | $\boldsymbol\tau$ / ns | 最小 fold gap / $q_K$ |
|---|---|---|---:|
| `B0_QC` | `[0,4.5,9,13.5,18,22.5,27,31.5]` | `[0,520.833,1041.667,1562.5,2083.333,2604.167,3125,3645.833]` | 4.5 |
| `S0_SIDON` | `[0,1,3,7,12,20,30,65]` | `[0,115.741,347.222,810.185,1388.889,2314.815,3472.222,7523.148]` | 1 |
| `AP_RMS_T1` | `[0,.864,1.728,2.592,3.456,4.320,5.184,6.048]` | `[0,100,200,300,400,500,600,700]` | 0.864 |
| `AP_TEPS_T1` | `[0,4.144262,8.288525,12.432787,16.577050,20.721312,24.865574,29.009837]` | `[0,479.660,959.320,1438.980,1918.640,2398.300,2877.960,3357.620]` | 4.144262 |
| `GEO_T1_CTRL` | `[0,3,10,40,130,191,205,271]` | `[0,347.222,1157.407,4629.630,15046.296,22106.481,23726.852,31365.741]` | 3 |
| `A100_SMALL_CDD_QSTEP0P25_MATCHED_CDD` | `[0,.125,.25,.375,.5,.625,.75,.875]` | `[0,14.468,28.935,43.403,57.870,72.338,86.806,101.273]` | 0.125 |
| `A100_SMALL_CDD_QSTEP0P25_TRANSPARENT_CDD` | `[0,.125,.25,.375,.5,.625,.75,.875]` | `[0,14.468,28.935,43.403,57.870,72.338,86.806,101.273]` | 0.125 |

### AL8：48 RB，$K=576$，$N_p=144$

| 方案 | $\mathbf j=\boldsymbol\tau/q_K$ | $\boldsymbol\tau$ / ns | 最小 fold gap / $q_K$ |
|---|---|---|---:|
| `B0_QC` | `[0,9,18,27,36,45,54,63]` | `[0,520.833,1041.667,1562.5,2083.333,2604.167,3125,3645.833]` | 9 |
| `S0_SIDON` | `[0,1,3,7,12,20,30,65]` | `[0,57.870,173.611,405.093,694.444,1157.407,1736.111,3761.574]` | 1 |
| `AP_RMS_T1` | `[0,1.728,3.456,5.184,6.912,8.640,10.368,12.096]` | `[0,100,200,300,400,500,600,700]` | 1.728 |
| `AP_TEPS_T1` | `[0,8.288525,16.577050,24.865574,33.154099,41.442624,49.731149,58.019674]` | `[0,479.660,959.320,1438.980,1918.640,2398.300,2877.960,3357.620]` | 8.288525 |
| `GEO_T1_CTRL` | `[0,7,86,101,302,353,511,540]` | `[0,405.093,4976.852,5844.907,17476.852,20428.241,29571.759,31250]` | 7 |
| `A100_SMALL_CDD_QSTEP0P25_MATCHED_CDD` | `[0,.25,.5,.75,1,1.25,1.5,1.75]` | `[0,14.468,28.935,43.403,57.870,72.338,86.806,101.273]` | 0.25 |
| `A100_SMALL_CDD_QSTEP0P25_TRANSPARENT_CDD` | `[0,.25,.5,.75,1,1.25,1.5,1.75]` | `[0,14.468,28.935,43.403,57.870,72.338,86.806,101.273]` | 0.25 |

precoder cycling 不属于 CDD 时延集合：`A100_PRG_DFT8_6RB` 在 AL2/4/8 分别使用 DFT `[0,1]`、`[0,1,2,3]`、`[0,1,2,3,4,5,6,7]`。

`S0_SIDON` 已逐一核验：三个 $K$ 下 36 个模 $K$ 无序二元和均唯一，fold residue 也均互异。所有表列方案的最小 fold gap 均严格大于零；连续方案仍须用真实 DMRS index 构造 $\mathbf V_P$，报告 rank/condition 和 CE floor。

A100 的硬厚 Sidon 要求/装填上界分别为 AL2 `5/4`、AL4 `9/8`、AL8 `17/16`，故均不可行。表中 `GEO_T1_CTRL` 是只读取 $T_{0.01}=479.660$ ns 和系统几何的固定种子软 Pareto 平衡代表，实际 pair/fold gap 分别为 `1/3`、`3/3`、`7/7`；选择分数为 $\min(d_\Sigma/g_\Sigma,d_{\rm fold}/g_{\rm fold})$，冻结后不得按 CE/BLER 回选，也不得称为硬厚 Sidon。

接收机口径冻结如下，不得由 candidate ID 隐式推断：

- `A100_PRG_DFT8_6RB`：透明；UE 不知道各 PRG 的 DFT 向量索引，逐 6-RB PRG 只使用底层物理 PDP 推导的 $\mathbf R_{phy}$ 做 LMMSE，不跨 PRG。
- `A100_SMALL_CDD_QSTEP0P25_TRANSPARENT_CDD`：透明；UE 不知道 CDD delay/`V`，在整个候选占用带宽只使用底层物理 PDP 推导的 $\mathbf R_{phy}$ 做 LMMSE。
- `A100_SMALL_CDD_QSTEP0P25_MATCHED_CDD`：non-transparent；保留原正式数据，UE 使用真实 $\mathbf R_g=\mathbf R_{phy}\odot(\mathbf V\mathbf V^H)$ 做全带 matched LMMSE。
- `B0_QC`、`S0_SIDON`、`AP_RMS_T1`、`AP_TEPS_T1`、`GEO_T1_CTRL`：non-transparent；各自使用真实 $\mathbf R_g$ 做全带 matched LMMSE。

执行：扩展 PDCCH runner；测试资源/预编码、无噪声译码及回归。smoke 后 300 trials/点、1 dB 预扫；正式前冻结 AL 网格，10%/1%附近≤0.25 dB。正式每点≥10,000 trials、200 errors、上限 50,000。原计划要求同 AL/SNR/trial 共享样本；研究者于 2026-09-04 明确允许各时延方案独立并行运行，因此正式执行改用由全局 seed、candidate ID 和固定标签唯一派生的独立随机流。相对增益区间按独立 bootstrap 计算，不声明 paired-sample 方差缩减。报告 BLER、目标 SNR 及 95%区间；未双侧 bracket 不外推。

2026-09-07 补充执行：保留已完成的 matched 小时延旧曲线，只为 transparent 小时延候选独立预扫和正式运行，不重跑其他候选。全局算法字段使用中性的 `channel_estimation: frequency_lmmse`，runner 另增加显式 `receiver_covariance_mode`；是否透明只由后者决定。单元测试核对 transparent filter 的权重与直接由 `R_phy` 构造的全带 LMMSE 完全一致、与 matched-effective filter 不同。正式配置允许该新增候选使用独立 `snr_points_db`，分析器按候选核对完整网格；统计预算、seed 派生、停止条件、目标插值与 bootstrap 口径保持不变。

交付：输出至 `outputs/experiment031_pdcch_cdd/20260903_main/`，生成两版 result。按 AL 报告相对 `B0_QC` 及 cycling 基线的 10%/1% SNR 差；不跨 AL 归因，确认前不更新全局结论。

## 2026-09-07 新增场景：TDL-C 300 ns、4Tx、2-symbol CORESET（待执行）

本节是在原 AL2/4/8、TDL-A 100 ns 正式结果之外新增的独立场景，不覆盖 `20260903_main` 原始数据。研究者明确要求把移动性纳入本轮，因此这是对 `GOALS.md` 当前“暂不纳入移动性”边界的本 plan 授权例外；结果确认前不据此更新全局结论。由于 AL、CORESET duration、TDL profile/DS、速度、载频和天线数同时改变，新旧场景之间只做描述性对照，不把差值归因于其中单一因素。

### 新场景配置与资源计数

除本节明确修改项外沿用本 plan：DCI payload A=41、CRC24C/RNTI `0xFFFF`、QPSK、1Rx、48-RB CORESET、30 kHz、FFT 4096、CP 288、non-interleaved、first CCE 0、L=6、单位总发射功率、data/DMRS 同预编码、仅 estimated-CSI BLER，以及既有 trial、停止、插值和独立 bootstrap 口径。修改为 AL1/2/4、2 CORESET symbols、4Tx、Sionna TDL-C 300 ns、3 km/h、4 GHz、20 sinusoids。速度为 0.833333 m/s，最大 Doppler 为 11.119 Hz。

对任一 AL，PDCCH candidate 固定包含 $6\,\mathrm{AL}$ 个 REG；一个 REG 是 `1 RB × 1 symbol`，即 12 个时频 RE。因此同一 AL 从 1-symbol 改为 2-symbol 时，candidate 的总时频 RE 数不变，只是每 symbol 的频域宽度减半。AL2/4/8 的资源计数对照为：

| AL | 1-symbol：每 symbol 实占 RB | 1-symbol：$K$ | 1-symbol：总时频 RE | 2-symbol：每 symbol 实占 RB | 2-symbol：$K$ | 2-symbol：总时频 RE |
|---:|---:|---:|---:|---:|---:|---:|
| 2 | 12 | 144 | 144 | 6 | 72 | 144 |
| 4 | 24 | 288 | 288 | 12 | 144 | 288 |
| 8 | 48 | 576 | 576 | 24 | 288 | 576 |

这里 $K$ 始终是单个 OFDM symbol 上的唯一子载波数，不把两个 symbol 合并计数；本新增场景实际使用 AL1/2/4，其中 AL1 的 2-symbol 计数为每 symbol 3 RB、$K=36$、两个 symbol 共 72 个时频 RE。

一个 CCE 的 6 REG 在 2-symbol non-interleaved 映射下覆盖 3 个连续 RB；每个 6-REG bundle 覆盖 `3 RB × 2 symbols`。因此

$$
N_{\rm RB,occ}=3\,\mathrm{AL},\qquad
K=36\,\mathrm{AL},\qquad
N_p=K/4,
$$

其中 $N_p$ 是每个 DMRS symbol 的唯一频率导频数；两个 symbol 的实际 DMRS 观测数为 $2N_p=18\,\mathrm{AL}$。data RE 为 $54\,\mathrm{AL}$，故 $E=108\,\mathrm{AL}$。PDCCH DMRS 仍为 comb-4，$\tau_{\rm alias}=8.333333\ \mu\mathrm{s}$。

| AL | 实占 RB | $K$ | $q_K=1/(K\Delta f)$ | 每 symbol $N_p$ | 总 DMRS RE | data RE | $E$ |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 3 | 36 | 925.925926 ns | 9 | 18 | 54 | 108 |
| 2 | 6 | 72 | 462.962963 ns | 18 | 36 | 108 | 216 |
| 4 | 12 | 144 | 231.481481 ns | 36 | 72 | 216 | 432 |

### 时延生成输入与规则

从 Sionna 1.0.2 TDL-C 300 ns 的 24-tap PDP，按最短连续区间覆盖至少 99% 能量计算得

$$
T_{\rm RMS}=300\ \mathrm{ns},\qquad
T_{0.01}=1891.950\ \mathrm{ns},
$$

对应区间 `[0,1891.950] ns`，覆盖功率 `0.99334723`。3 km/h 和 4 GHz 不改变 PDP 时延 candidate，只进入两个 OFDM symbol 的时间协方差。CDD 使用

$$
V_{k,n}=\frac12\exp(-j2\pi k j_n/K),
\qquad \tau_n=j_nq_K.
$$

- `B0_QC`：不重新搜索，取 8Tx 历史物理时延的前四支 `[0,520.833,1041.667,1562.500] ns`，再换算各 AL 的连续 $j$。它保持了原等差公差，但没有保持原 B0_QC 的核心 QC 性质：在旧 comb-24、8Tx 条件下，520.833 ns 等于 \(3\tau_{\rm alias}/8\)，折叠后正好形成均匀的 8 点排列。换成 4Tx、PDCCH comb-4 后，直接截取前四路不再是均匀的 4 点折叠设计，也不能称为 4Tx 优化结果。如果按原 QC 准则为当前 4Tx PDCCH 重新构造，不需要穷搜。PDCCH DMRS 频域间隔为 \(S_f=4\)，因此：
$$
\tau_{\rm alias}=\frac1{4\Delta f}=8.3333\ \mu s,
\qquad
\Delta\tau_{\rm QC}=\frac{\tau_{\rm alias}}{4}=2.0833\ \mu s.
$$
推荐时延为：
$$
[0,2.0833,4.1667,6.2500]\ \mu s.
$$
对应各 AL 的归一化坐标 $j=\tau K\Delta f$：
AL	\(K\)	4Tx QC 坐标
AL1	36	[0, 2.25, 4.5, 6.75]
AL2	72	[0, 4.5, 9, 13.5]
AL4	144	[0, 9, 18, 27]
- `S0_SIDON`：取历史严格 Sidon 的四元素前缀 `[0,1,3,7]`；已分别在 $\mathbb Z_{36}$、$\mathbb Z_{72}$、$\mathbb Z_{144}$ 核验 10 个无序二元和唯一且 fold residue 互异。
- `AP_RMS_T1`：$\tau_n=nT_{\rm RMS}$；`AP_TEPS_T1`：$\tau_n=nT_{0.01}$，均不做整数舍入。
- 两条 `SMALL_CDD_QSTEP0P25`：取历史 48-RB 小时延的前四支 `[0,14.468,28.935,43.403] ns`；matched/transparent 发射 `V` 完全相同。QSTEP0P25 表示“小公差等差 CDD”，其中 0P25 即参考坐标中的 0.25，不是 0.25 ns，也不是 DMRS 折叠周期的 0.25。
原始定义以 48 RB、\(K=576\) 为参考：
$$
q=\frac{1}{K\Delta f}
  =\frac{1}{576\times30\,{\rm kHz}}
  =57.87037\ {\rm ns}.
$$
相邻天线的人工时延公差为：
$$
\Delta\tau=0.25q=14.4676\ {\rm ns}.
$$
因此原来的 8Tx 时延为：
$$
\mathbf j=[0,0.25,0.5,\ldots,1.75],
$$
$$
\boldsymbol\tau
=[0,14.468,28.935,43.403,57.870,72.338,86.806,101.273]\ {\rm ns}.
$$
- `GEO_T1_CTRL`：取 $T_{\rm margin}=T_{\rm sync}=0$。固定 $j_0=0$，穷举其余三个互异整数坐标，要求 fold residue 互异；计算模 $K$ 的 10 个 pair-sum 最小圆周间距 $d_\Sigma$ 和模 $N_p$ 的 $d_{\rm fold}$。分别冻结最大 pair-gap、最大 fold-gap及最大化 $\min(d_\Sigma/g_\Sigma,d_{\rm fold}/g_{\rm fold})$ 的平衡代表；同分依次按另一 gap 较大、字典序较小决胜，等价数组去重。该过程不读取完整 PDP 功率、协方差、CE 或 BLER。

硬厚 Sidon 的门限和装填上界为：

| AL | $g_\Sigma/g_{\rm fold}$ | pair/fold 装填上界 | 判定 |
|---:|---:|---:|---|
| 1 | `5/3` | `3/2` | 两项均不可行 |
| 2 | `9/5` | `7/4` | 两项均不可行 |
| 4 | `17/9` | `14/9` | pair 项不可行 |

因此不存在硬厚 Sidon；以下几何候选全部保持 family 名 `GEO_T1_CTRL` 并明确标为软控制，不能称作硬厚 Sidon。

### AL1 candidate：3 RB，$K=36$，$N_p=9$

| candidate | $\mathbf j=\boldsymbol\tau/q_K$ | $\boldsymbol\tau$ / ns | pair/fold gap | 备注 |
|---|---|---|---:|---|
| `C300_B0_QC` | `[0,.5625,1.125,1.6875]` | `[0,520.833,1041.667,1562.500]` | `—/.5625` | 历史物理参考 |
| `C300_S0_SIDON` | `[0,1,3,7]` | `[0,925.926,2777.778,6481.481]` | `1/1` | 严格 Sidon |
| `C300_AP_RMS_T1` | `[0,.324,.648,.972]` | `[0,300,600,900]` | `—/.324` | 连续 RMS 等差 |
| `C300_AP_TEPS_T1` | `[0,2.043306,4.086612,6.129918]` | `[0,1891.950,3783.900,5675.850]` | `—/2.043306` | 连续 99% 支撑等差 |
| `C300_GEO_T1_CTRL_BAL_PAIR` | `[0,2,6,14]` | `[0,1851.852,5555.556,12962.963]` | `2/1` | 平衡与最大 pair 代表重合 |
| `C300_GEO_T1_CTRL_FOLD` | `[0,2,5,16]` | `[0,1851.852,4629.630,14814.815]` | `1/2` | 最大 fold 代表 |
| `C300_SMALL_CDD_QSTEP0P25_MATCHED_CDD` | `[0,.015625,.03125,.046875]` | `[0,14.468,28.935,43.403]` | `—/.015625` | non-transparent |
| `C300_SMALL_CDD_QSTEP0P25_TRANSPARENT_CDD` | 同上 | 同上 | `—/.015625` | transparent |

### AL2 candidate：6 RB，$K=72$，$N_p=18$

| candidate | $\mathbf j=\boldsymbol\tau/q_K$ | $\boldsymbol\tau$ / ns | pair/fold gap | 备注 |
|---|---|---|---:|---|
| `C300_B0_QC` | `[0,1.125,2.25,3.375]` | `[0,520.833,1041.667,1562.500]` | `—/1.125` | 历史物理参考 |
| `C300_S0_SIDON` | `[0,1,3,7]` | `[0,462.963,1388.889,3240.741]` | `1/1` | 严格 Sidon |
| `C300_AP_RMS_T1` | `[0,.648,1.296,1.944]` | `[0,300,600,900]` | `—/.648` | 连续 RMS 等差 |
| `C300_AP_TEPS_T1` | `[0,4.086612,8.173224,12.259836]` | `[0,1891.950,3783.900,5675.850]` | `—/4.086612` | 连续 99% 支撑等差 |
| `C300_GEO_T1_CTRL_BAL` | `[0,4,12,33]` | `[0,1851.852,5555.556,15277.778]` | `4/3` | 平衡代表 |
| `C300_GEO_T1_CTRL_PAIR` | `[0,5,16,49]` | `[0,2314.815,7407.407,22685.185]` | `5/2` | 最大 pair 代表 |
| `C300_GEO_T1_CTRL_FOLD` | `[0,4,10,32]` | `[0,1851.852,4629.630,14814.815]` | `2/4` | 最大 fold 代表 |
| `C300_SMALL_CDD_QSTEP0P25_MATCHED_CDD` | `[0,.03125,.0625,.09375]` | `[0,14.468,28.935,43.403]` | `—/.03125` | non-transparent |
| `C300_SMALL_CDD_QSTEP0P25_TRANSPARENT_CDD` | 同上 | 同上 | `—/.03125` | transparent |

### AL4 candidate：12 RB，$K=144$，$N_p=36$

| candidate | $\mathbf j=\boldsymbol\tau/q_K$ | $\boldsymbol\tau$ / ns | pair/fold gap | 备注 |
|---|---|---|---:|---|
| `C300_B0_QC` | `[0,2.25,4.5,6.75]` | `[0,520.833,1041.667,1562.500]` | `—/2.25` | 历史物理参考 |
| `C300_S0_SIDON` | `[0,1,3,7]` | `[0,231.481,694.444,1620.370]` | `1/1` | 严格 Sidon |
| `C300_AP_RMS_T1` | `[0,1.296,2.592,3.888]` | `[0,300,600,900]` | `—/1.296` | 连续 RMS 等差 |
| `C300_AP_TEPS_T1` | `[0,8.173224,16.346448,24.519672]` | `[0,1891.950,3783.900,5675.850]` | `—/8.173224` | 连续 99% 支撑等差 |
| `C300_GEO_T1_CTRL_BAL` | `[0,10,30,97]` | `[0,2314.815,6944.444,22453.704]` | `10/5` | 平衡代表 |
| `C300_GEO_T1_CTRL_PAIR` | `[0,11,44,66]` | `[0,2546.296,10185.185,15277.778]` | `11/3` | 最大 pair 代表 |
| `C300_GEO_T1_CTRL_FOLD` | `[0,9,18,27]` | `[0,2083.333,4166.667,6250]` | `0/9` | 最大 fold 代表；非严格 Sidon |
| `C300_SMALL_CDD_QSTEP0P25_MATCHED_CDD` | `[0,.0625,.125,.1875]` | `[0,14.468,28.935,43.403]` | `—/.0625` | non-transparent |
| `C300_SMALL_CDD_QSTEP0P25_TRANSPARENT_CDD` | 同上 | 同上 | `—/.0625` | transparent |

precoder cycling 使用新 ID `C300_PRG_DFT4_6REG`，4 维 DFT 向量归一化为单位范数；每个 6-REG bundle 内 data/DMRS 同向量，AL1/2/4 分别轮询 `[0]`、`[0,1]`、`[0,1,2,3]`。配置中的 bundle size 均为 `L=6 REG`：旧 1-symbol 场景对应 `6 RB × 1 symbol`（频域 72 子载波、共 72 时频 RE），本 2-symbol 场景对应 `3 RB × 2 symbols`（每 symbol 频域 36 子载波、共 72 时频 RE）。2-symbol 的物理 bundle 跨两个 symbol，但不在 symbol 之间切换 precoder；同一 3-RB bundle 的两个 symbol 使用同一个 DFT 向量，向量只在相邻频域 bundle 之间 cycling。

### 接收机、实现与执行

- `C300_PRG_DFT4_6REG` 和 `...TRANSPARENT_CDD` 保持透明；前者逐 6-REG bundle 使用物理 TDL-C 时频协方差，后者全候选带宽使用物理时频协方差。其他 CDD 均为 non-transparent，使用已知 $V$ 构成的真实等效时频协方差。
- 两个 DMRS symbol 不做静态平均；使用真实 `(symbol,subcarrier)` pilot/data 坐标和已知 3 km/h、4 GHz 的时间相关，构造二维 time-frequency LMMSE。CDD 的时间相关沿用物理 TDL-C，频率相关乘 $VV^H$；cycling 不跨 bundle 估计。
- 扩展 `cdd_lls/sim/pdcch_cdd.py` 当前写死的 `8Tx/1-symbol/AL2,4,8/TDL-A100/static` 校验；$K$ 必须由候选实际唯一子载波数得到，2-symbol 下不得继续使用旧式 `72×AL`。复用现有 active-only TDL 和二维 LMMSE 公共实现，不复制链路主循环。
- 测试 AL1/2/4 的 CCE/REG/RB、data/DMRS/E 计数，4Tx 相位与 DFT4 cycling，全部 fold 唯一性、S0 pair-sum、GEO 穷举回执、TDL-C $T_{0.01}$、时频协方差、无噪声 DCI 解码及原 AL2/4/8 回归。smoke 必须实际覆盖两个 symbol 且验证非零 Doppler realization 在 symbol 间变化。
- seed 使用 `20260907`，各 candidate 保持独立派生随机流。先以 300 trials/点、1 dB 步长预扫；正式 trial 1 前按 AL 和 candidate 冻结网格，10%/1% 附近步长不超过 0.25 dB。正式每点至少 10,000 trials、目标 200 errors、上限 50,000；未双侧 bracket 不外推。

新增输出写入 `outputs/experiment031_pdcch_cdd/20260907_c300_4tx_2sym/`，保存展开配置、候选生成与穷举回执、PDP/$T_{0.01}$ 回执、原始 BLER/CE/trial 数据、日志和环境信息。执行后在两版 `result-031` 增加独立补充章节，按 AL 报告相对 `C300_B0_QC` 与 `C300_PRG_DFT4_6REG` 的 10%/1% SNR 差和 95% 区间，并单列 pilot rank/condition、二维 CE NMSE 与零噪声 CE floor。

## 2026-09-10 补充规划：C300 2-symbol 场景 ideal-CSI BLER（已执行，待确认）

### 目的与范围

在已完成的 TDL-C 300 ns、4Tx、2-symbol CORESET 场景中补充 AL1/2/4 的 ideal-CSI DCI BLER 曲线，用于分离“发射波形本身的分集/编码性能”与“真实 DMRS 信道估计损失”。除接收机改为 ideal CSI 外，保持本场景的 4Tx/1Rx、48-RB CORESET、30 kHz、TDL-C 300 ns、3 km/h、4 GHz、A=41、CRC24C、QPSK、资源映射、总发射功率、DMRS 开销、data/DMRS 预编码和各 candidate 时延不变。ideal CSI 直接使用每个 data RE 上生成的真实时变等效信道做 coherent MRC/equalization；DMRS RE 仍保留且不承载数据，但不参与信道估计。

本补充覆盖前述 C300 主实验中的不同发射波形，不纳入后续“严格 Sidon shortlist”新增候选。AL1 使用 `C300_B0_QC`、`C300_S0_SIDON`、`C300_AP_RMS_T1`、`C300_AP_TEPS_T1`、`C300_GEO_T1_CTRL_BAL_PAIR`、`C300_GEO_T1_CTRL_FOLD`、`C300_PRG_DFT4_6REG` 和 small CDD；AL2/4 使用对应的 B0、S0、AP-RMS、AP-$T_{0.01}$、三个 GEO、PRG 和 small CDD。`C300_SMALL_CDD_QSTEP0P25_MATCHED_CDD` 与 `...TRANSPARENT_CDD` 的发射波形完全相同，ideal CSI 下接收机差异消失，因此只运行一条、命名为 `C300_SMALL_CDD_QSTEP0P25_IDEAL_CSI`，不得把同一波形重复计为两条独立结果。其他 family 名和 candidate 名保持既有称呼。

### 实现、统计与执行

- 扩展 `pdcch-cdd-bler-v2` 接受 `receiver.channel_estimation: ideal`。ideal 分支必须令 `data_estimate=data_true`，不构造或调用 LMMSE filter；`receiver_covariance_mode` 对该分支无效。输出仍记录真实 candidate、AL、时频资源与 BLER，CE NMSE 固定为 0，仅作为实现回执，不进入 candidate 排名。
- 测试 AL1/2/4 下 ideal 分支的 CE 误差逐元素为 0、无噪声 DCI 全部通过、MRC 分母与噪声方差正确、2-symbol 非零 Doppler 信道仍随 symbol 变化，以及 small CDD 两个历史接收口径在相同样本下产生完全相同的 ideal-CSI LLR/判决。estimated-CSI 原路径必须回归通过。
- 新建 AL1/2/4 的 prescan/formal YAML，seed 使用 `20260910`，按 AL、candidate 和 stage 稳定派生独立随机流。先在 `-12:2:18 dB` 固定宽网格上每点运行 300 trials；若 10% 或 1% BLER 未被双侧 bracket，则以 2 dB 步长只向缺失方向扩展。随后冻结每个 candidate 的正式网格，两个目标附近步长不超过 0.25 dB。
- 正式每点至少 10,000 trials、累计 200 errors 后可停止、上限 50,000 trials，保存 trial error flags、展开配置、日志和环境回执。目标 SNR 使用与本 plan 既有结果相同的 isotonic 加局部插值；95% 区间继续使用独立 bootstrap。未双侧 bracket 时不外推，只报告未取得目标。

输出写入 `outputs/experiment031_pdcch_cdd/20260910_c300_4tx_2sym_ideal_csi/`。在两版 `result-031` 增加配对章节：每个 AL 分别给出 ideal-CSI BLER 图、10%/1% SNR、相对 `C300_B0_QC` 和 `C300_PRG_DFT4_6REG` 的 SNR 差及 95% 区间；再对具有合格 bracket 的同名 waveform 报告 $\Delta_{\rm CE}=\mathrm{SNR}_{\rm estimated}-\mathrm{SNR}_{\rm ideal}$。该差值用于量化当前接收机造成的总 penalty，不得解释为只由频域、时域或 PRG 边界中的单一因素导致。验收条件是三组 AL 的所有冻结 candidate 完成预算或明确达到上限、无噪声与 ideal-CE 测试通过、两版 result 数值和结论一致；结果确认前不更新 `KNOWLEDGE.md` 或 `GOALS.md`。

## 2026-09-11 新增规划：C300 4Tx/2Rx、2-symbol 的 ideal/estimated-CSI BLER 与 NMSE（部分已执行：AL1/AL2 已完成，AL4 未完成）

截至 2026-09-16，AL1/AL2 的 estimated/ideal-CSI 正式仿真、NMSE 汇总和阶段性分析已经完成，结果见两版 `result-031` 第 10 节；AL4 的 estimated-CSI 目录仍有未完成 candidate，ideal-CSI 正式目录尚不存在。因此本节不能标为“全部已执行”，后续只需完成和分析 AL4，不得重跑或覆盖已经完成的 AL1/AL2 原始数据。

### 目的、唯一物理变量与结论边界

本节在前述 TDL-C 300 ns、4Tx/1Rx、2-symbol CORESET 主实验及其 ideal-CSI 补充的基础上，只把接收天线数从 `n_rx: 1` 改为 `n_rx: 2`，用于观察无接收相关假设下的双分支接收对 ideal-CSI BLER、estimated-CSI BLER 和信道估计 NMSE 的影响。除 `n_rx` 外，发射端、信道每分支的边缘统计、时频资源、编码调制、DMRS、候选和接收算法定义全部保持不变；本节不重新搜索、筛选或修改任何时延候选，也不纳入后续“严格 Sidon shortlist”候选。

物理配置严格沿用本 plan 的“TDL-C 300 ns、4Tx、2-symbol CORESET”新增场景：AL1/2/4、2 个相邻 CORESET symbols、4Tx、48-RB CORESET、30 kHz、FFT 4096、CP 288、non-interleaved、first CCE 0、`L=6 REG`、TDL-C 300 ns、3 km/h、4 GHz、20 sinusoids、A=41、CRC24C/RNTI `0xFFFF`、QPSK、data/DMRS 同预编码，以及单位总发射功率。所有预编码行继续满足 $\sum_{n=1}^{4}|V_{k,n}|^2=1$；不得改成每 Tx 单位功率。SNR 按 `DESIGN.md` 第 3.9 节定义为单位总发射信号功率相对单个 Rx 分支噪声功率的比值，故每个接收分支的复高斯噪声方差仍为 $\sigma^2=10^{-\mathrm{SNR}_{\rm dB}/10}$，不随 `n_rx` 额外除以 2；双 Rx 的收益来自独立分支经 MRC 后的阵列/分集作用，而不是修改输入 SNR 或总发射功率。

本轮固定假设两个 Rx 分支不相关：对每个 Tx-Rx 对独立生成具有相同 TDL-C PDP、Doppler 和平均功率归一化的信道，不引入空间相关矩阵、互耦、分支增益失配或分支噪声相关。因此结果只回答“默认独立双 Rx”模型，不能外推到共址相关天线、极化相关或不等噪声分支。

### 冻结候选与配置继承

estimated-CSI 侧逐字段继承以下 1Rx 正式配置：

- AL1：`configs/pdcch_result031_c300_al1_formal.yaml`；
- AL2：`configs/pdcch_result031_c300_al2_formal.yaml`；
- AL4：`configs/pdcch_result031_c300_al4_formal.yaml`。

ideal-CSI 侧逐字段继承：

- AL1：`configs/pdcch_result031_c300_ideal_al1_formal.yaml`；
- AL2：`configs/pdcch_result031_c300_ideal_al2_formal.yaml`；
- AL4：`configs/pdcch_result031_c300_ideal_al4_formal.yaml`。

候选集合及精确 $\mathbf j$、物理时延、phase denominator 和 cycling order 必须与第 186--227 行的 C300 主候选完全一致。AL1 保留 B0、S0、AP-RMS、AP-$T_{0.01}$、两个 GEO、PRG DFT4、small CDD matched/transparent；AL2/4 保留 B0、S0、AP-RMS、AP-$T_{0.01}$、三个 GEO、PRG DFT4、small CDD matched/transparent。ideal CSI 下 small CDD 的 matched/transparent 仍是同一发射波形，只运行一条 ideal 曲线，不重复计数。

新配置相对对应 1Rx 源配置的允许差异分为两类：

1. 唯一物理差异：`antenna.n_rx: 1 -> 2`。
2. 实验控制差异：新输出目录、seed/随机流 namespace、prescan/formal SNR 网格、trial 预算和 batch size。estimated 配置必须继续使用源配置的 `receiver.channel_estimation: frequency_lmmse`，ideal 配置必须继续使用源配置的 `receiver.channel_estimation: ideal`；这些控制字段不得改变 Tx 数、候选参数、信道、资源、波形或对应接收机口径。

在 smoke 前生成 `config_diff_receipt.json`：对每个 AL 和接收机模式记录源配置、新配置及 SHA-256，并在忽略上述允许字段后递归比较。出现其他字段差异即停止，不进入预扫。特别地，候选 ID、family、scheme、`delay_grid_coordinates`、`q_seconds`、`phase_denominator`、`receiver_covariance_mode` 和 cycling order 必须逐项相同；ideal small CDD 只允许沿用前节已定义的单曲线命名映射。

### 双 Rx 接收、信道估计、MRC 与 NMSE

每个 Rx 分支独立执行与 1Rx 完全相同的信道估计，不在 Rx 维联合估计，也不使用另一分支的 DMRS。estimated-CSI 分支继续使用 `channel_estimation: frequency_lmmse`：matched CDD 使用真实等效时频协方差，transparent full-band/PRG 使用各自既有物理协方差口径，两个 DMRS symbol 使用真实时频坐标且不做静态平均。ideal-CSI 分支直接令每个 Rx、每个 data RE 上的 $\hat g_r=g_r$，不得构造 LMMSE filter。

两种 CSI 模式均在 Rx 维做 coherent MRC：

$$
\hat x=\frac{\sum_{r=1}^{2}\hat g_r^*y_r}
{\max\left(\sum_{r=1}^{2}|\hat g_r|^2,\epsilon\right)},
\qquad
\sigma_{\rm eff}^2=
\frac{\sigma^2}{\max\left(\sum_{r=1}^{2}|\hat g_r|^2,\epsilon\right)}.
$$

demapper、LLR clipping、Polar 解码和 CRC 判错逻辑保持不变。estimated-CSI 的逐 trial NMSE 沿用现有 ratio-of-sums 定义，但在两个 Rx 和全部 data RE 上共同求和：

$$
\mathrm{NMSE}_b=
\frac{\sum_{r=1}^{2}\sum_{q\in\mathcal D}|\hat g_{b,r,q}-g_{b,r,q}|^2}
{\max\left(\sum_{r=1}^{2}\sum_{q\in\mathcal D}|g_{b,r,q}|^2,10^{-30}\right)}.
$$

曲线点报告 $10\log_{10}\!\left(E_b[\mathrm{NMSE}_b]\right)$，不得改成先逐 trial 转 dB 再平均。ideal-CSI NMSE 必须逐元素精确为 0，只作实现回执，不绘为有限 dB 曲线，也不参与候选排序。由于两分支边缘分布相同，增加 Rx 不应系统性改变 NMSE 的总体期望；可能出现的差异首先按有限样本方差和 ratio 统计量解释，不能预设 2Rx 必然改善信道估计器本身。

### 配置、预扫与正式统计

新建 AL1/2/4 × estimated/ideal 的 smoke、prescan 和 formal YAML，统一写入 `outputs/experiment031_pdcch_cdd/20260911_c300_4tx_2rx_2sym/`。全局 seed 使用 `20260911`；estimated 与 ideal 使用不同且固定的 random-stream namespace，各 candidate 继续由全局 seed、AL、candidate ID 和固定标签派生独立随机流。因此候选间、CSI 模式间以及与历史 1Rx 结果之间均按独立样本处理，不声明 paired-sample 方差缩减。`batch_size` 在所有正式配置中冻结为 50；只有 smoke 证明资源充足且在正式 trial 1 前统一修改配置与回执时才可调整，正式开始后不得变更。

1. smoke：每个 AL、每种 CSI 模式、每个 candidate 至少运行 20 trials 和一个有限 SNR 点；另运行 ideal 的无噪声点。smoke 结果与正式结果隔离，不能合并计数。
2. prescan：两种 CSI 模式统一先在 `-12:2:18 dB` 上每点运行 300 trials。若任一 candidate 的 10% 或 1% BLER 未被双侧 bracket，以 2 dB 步长只向缺失方向扩展。NMSE 与 BLER 使用相同 estimated-CSI trial，不另开随机样本。
3. freeze：按 AL、candidate 和 CSI 模式冻结正式 SNR 网格；10%/1% 目标附近的相邻原始点间隔不超过 0.25 dB。正式 trial 1 之前把精确网格、配置 SHA-256 和 prescan 选择回执写入 `prescan_freeze_receipt.json`；冻结前不得开始正式运行。
4. formal：每点至少 10,000 trials，累计 200 个误块后可停止，上限 50,000 trials。已有点追加必须从绝对 trial 区间继续，不能从 trial 1 重跑后相加。达到上限仍不足 200 errors 可以保留为受限点，但目标未被双侧 bracket 时不外推。

各 BLER 点保存 errors/trials、Wilson 95% 区间和逐 trial error flags；各 estimated-CSI 点同时保存逐 trial 线性 NMSE。目标 SNR 继续使用递减 isotonic 加目标邻域局部插值；10%/1% 的 candidate gain、2Rx estimated-ideal penalty 以及 2Rx 相对既有 1Rx 的差值均用独立 bootstrap 给出 95% 区间。历史 1Rx 数据只从第 7、8 节已核验的正式 CSV 读取，不重跑、不与 2Rx trial 合并；若两边任一目标没有合格 bracket，则该差值不报告。

### 实现、测试与 smoke 验收

本节复用 `cdd_lls/sim/pdcch_cdd.py` 的可配置 `n_rx`、逐 Rx LMMSE 和 MRC 实现，不复制链路主循环，也不修改候选生成算法。若现有 runner 或分析器尚不能满足以下检查，只做支撑本节的最小扩展并补充定向测试：

- 配置校验接受正整数 `n_rx=2`，拒绝 0、负数、布尔值和非整数；1Rx 路径数值回归不变。
- active-only TDL 输出和接收张量的 Rx 维为 2；固定 seed 下两个 Rx 分支不逐元素相同，经验交叉相关只作 smoke 回执，不作为总体独立性的证明。
- 固定小型数组手算两分支 MRC 的合并值、分母和 $\sigma_{\rm eff}^2$；单分支置零时退化为另一分支结果。
- 对 estimated CSI，逐 Rx 单独调用同一个 LMMSE 线性算子，测试一个分支的 pilot 扰动不会改变另一个分支的估计；NMSE 与上述双重求和公式一致。
- 对 ideal CSI，两个分支的 `data_estimate` 均与真实等效信道逐元素相等、CE NMSE 精确为 0、无 LMMSE filter，有限噪声和无噪声解调均走同一 MRC 路径。
- 核对任一 data/DMRS RE 的预编码总功率为 1、两个 Rx 分支的噪声方差均为 $10^{-\mathrm{SNR}_{\rm dB}/10}$，且未因 `n_rx=2` 或 `n_tx=4` 再缩放。
- AL1/2/4 的 CCE/REG/RB、data/DMRS/E 计数、两个 symbol 的非零 Doppler 变化、全部候选几何回执和 exact delay 数组与 1Rx 源配置一致。

smoke 验收还要求所有输出无 NaN/Inf，BLER 在低/高 SNR 方向总体合理，ideal CSI 不得系统性劣于同 waveform 的 estimated CSI；若出现反序，先增加 smoke 样本并核查 LLR、MRC 分母、noise scaling 和随机流，不直接据少量样本作物理结论。只有配置差分、单测、无噪声链路和 smoke 全部通过后才允许执行 prescan。

### 产物、分析图与完成判据

输出目录至少包含：展开 YAML、`config_diff_receipt.json`、环境与依赖版本、测试日志、smoke 日志、`prescan_freeze_receipt.json`、每个 AL/candidate/CSI 模式的 `bler_points.csv`、逐 trial error flags、estimated-CSI 逐 trial NMSE、`run_metadata.json`，以及汇总 `formal_points.csv`、`target_snr.csv`、`target_gains.csv`、`delta_ce.csv`、`delta_rx.csv` 和分析元数据。所有 CSV 中必须显式记录 `n_tx=4`、`n_rx=2`、CSI 模式、AL、candidate ID、SNR 定义标签和总发射功率归一化标签。

两版 `result-031` 增加同名补充章节；每个 AL 至少给出三张主图：

1. 2Rx ideal-CSI BLER 曲线；
2. 2Rx estimated-CSI BLER 曲线；
3. 2Rx estimated-CSI NMSE 曲线。

BLER 图必须显示 10%/1% 参考线和原始正式点，不用插值曲线遮蔽非单调或 error floor；NMSE 图以同一 AL 内统一坐标展示所有 estimated-CSI candidate。表格分别报告：2Rx 的 10%/1% 目标 SNR、相对 `C300_B0_QC` 和 `C300_PRG_DFT4_6REG` 的差值、$\Delta_{\rm CE}^{(2Rx)}=\mathrm{SNR}_{\rm estimated,2Rx}-\mathrm{SNR}_{\rm ideal,2Rx}$，以及 $\Delta_{\rm Rx}=\mathrm{SNR}_{\rm 1Rx}-\mathrm{SNR}_{\rm 2Rx}$。$\Delta_{\rm Rx}>0$ 表示 2Rx 所需 SNR 更低；它包含独立分支的 MRC 阵列/分集收益，不能解释为信道估计算法自身获得同等 dB 改善。

完成判据为：6 组 AL×CSI 正式配置的全部冻结 candidate 均完成预算或明确达到 50,000 上限；三个要求的曲线族、原始数据、目标及差值区间齐全；配置差分证明唯一物理变量为 Rx 数；测试和数值审计通过；图文两版 `result-031` 的数字、限制和结论一致。研究者确认结果前，不把本节结果写入 `KNOWLEDGE.md`，也不据此改变 `GOALS.md` 的阶段验收状态。

### 仿真后梳理：信道估计、接收合并与解调实现

本小节记录本轮实际使用的数学模型和代码路径。以下公式中的预编码矩阵均指代码里的实际归一化矩阵 `candidate.precoder.C`，记为 $\mathbf C$；它已经包含 $1/\sqrt{N_t}$，每个 RE 满足 $\sum_n|C_{k,n}|^2=1$，因此公式中不得再除一次 $N_t$。

#### 1. 符号、底层信道与等效信道

- $N_t=4$、$N_r=2$：发射和接收天线数；$b$：Monte Carlo trial 索引；$r\in\{0,\ldots,N_r-1\}$：Rx 索引；$n\in\{0,\ldots,N_t-1\}$：Tx 索引。
- $s$：CORESET 内 OFDM symbol 索引；$k,l$：候选占用带宽内的子载波索引；$q=(s,k)$：一个时频 RE 坐标。
- $H_{b,r,n,s,k}$：第 $b$ 个 trial 中，Tx $n$ 到 Rx $r$ 在 RE $(s,k)$ 上的底层物理频域信道。
- $C_{k,n}$：Tx $n$ 在子载波 $k$ 上的复预编码系数；同一候选的 data 与 DMRS 使用相同 $\mathbf C$。
- $g_{b,r,s,k}$：预编码后、Rx $r$ 实际看到的单层等效信道。
- $x_{b,s,k}$、$y_{b,r,s,k}$ 和 $n_{b,r,s,k}$：发送符号、接收符号和复高斯噪声；上标 $*$ 表示复共轭，$(\cdot)^H$ 表示共轭转置，$j=\sqrt{-1}$。

等效信道和接收模型为

$$
g_{b,r,s,k}
=
\sum_{n=0}^{N_t-1}C_{k,n}H_{b,r,n,s,k},
\qquad
y_{b,r,s,k}=g_{b,r,s,k}x_{b,s,k}+n_{b,r,s,k},
$$

其中 $n_{b,r,s,k}\sim\mathcal{CN}(0,\sigma_n^2)$ 表示零均值、复方差为 $\sigma_n^2=\mathbb E[|n|^2]$ 的圆对称复高斯噪声，且 $\sigma_n^2=10^{-\mathrm{SNR}_{\rm dB}/10}$。当前模型假设不同 Tx、不同 Rx 的底层分支独立且同分布：

$$
\mathbb E\!\left[
H_{r,n,s,k}H_{r',n',s',l}^*
\right]
=
\delta_{r,r'}\delta_{n,n'}
R_t[s,s']R_{\rm phy}[k,l].
$$

$\delta_{a,b}$ 是 Kronecker delta；$R_t[s,s']$ 是两个 OFDM symbol 时刻之间的时间相关；$R_{\rm phy}[k,l]$ 是单个底层 Tx--Rx 分支的物理频率协方差。两个 Rx 没有空间相关，所以不构造跨 Rx 协方差，也不做跨 Rx 联合信道估计。

#### 2. 物理协方差与预编码后的等效协方差

对 TDL 的第 $u$ 条径，记归一化功率为 $p_u$、物理时延为 $\tau_u$，且 $\sum_up_u=1$；子载波 $k$ 的基带频率为 $f_k$。平台直接在候选实际占用的 active subcarriers 上计算

$$
R_{\rm phy}[k,l]
=
\sum_u p_u
\exp\!\left[-j2\pi\tau_u(f_k-f_l)\right].
$$

时间矩阵 $\mathbf R_t$ 由 Sionna 的 TDL time-covariance helper 根据 TDL profile、速度、载频、OFDM symbol duration 和 symbol 数生成：

$$
R_t[s,s']
:=
\mathbb E\!\left[H_{r,n,s,k}H_{r,n,s',k}^*\right].
$$

实现采用 TDL 的时频可分近似，所以底层任意两个 RE $q=(s,k)$、$q'=(s',l)$ 的协方差是 $R_t[s,s']R_{\rm phy}[k,l]$。在 Tx 分支独立的假设下，真实等效信道协方差为

$$
R_g[q,q']
=
R_t[s,s']R_{\rm phy}[k,l]
\sum_{n=0}^{N_t-1}C_{k,n}C_{l,n}^*.
$$

写成频率矩阵即

$$
\mathbf R_{g,f}
=
\mathbf R_{\rm phy}\odot(\mathbf C\mathbf C^H),
$$

其中 $\odot$ 为 Hadamard 逐元素乘法。`pdcch_cdd.py` 先计算 `physical_covariance`，再以 `band_precoder @ band_precoder.conj().T` 构造上述预编码相关因子。需要区分“生成真实信道所用的 $\mathbf C$”和“接收机假设的协方差”：真实接收符号始终由实际 $\mathbf C$ 形成，`receiver_covariance_mode` 只决定 LMMSE 权重使用哪一个协方差。

- `matched_effective`：接收机使用真实 $\mathbf R_g$，知道当前预编码对等效协方差的影响。
- `physical_fullband`：接收机使用 $\mathbf R_{\rm phy}$，忽略人工预编码因子，并允许整个候选带宽内的 pilot 共同估计 data RE。
- `physical_prg`：同样使用物理协方差，但每个 PRG 独立构造局部滤波器；跨 PRG 的权重固定为 0。
- `ideal`：不构造协方差滤波器，直接把生成器给出的真实 data-RE 等效信道作为 $\widehat g$。

#### 3. 从 DMRS 到全部待解调 data RE 的二维 LMMSE

记 $P$ 为两个 CORESET symbols 中全部 DMRS RE 的坐标集合，$D$ 为全部待解调 data RE 的坐标集合。对每个 trial 和每个 Rx 分支，DMRS 观测满足

$$
y^{(P)}_{r,p}=x^{(P)}_p g_{r,p}+n_{r,p},
\qquad p\in P.
$$

先除以已知 DMRS 符号得到 LS 观测

$$
z_{r,p}
:=
\frac{y^{(P)}_{r,p}}{x^{(P)}_p}
=
g_{r,p}+\widetilde n_{r,p}.
$$

本仿真的 DMRS 为单位模长，所以 $\operatorname{var}(\widetilde n_{r,p})=\sigma_n^2$。若一般导频不是单位模长，则相应 LS 噪声方差应为 $\sigma_n^2/|x^{(P)}_p|^2$，不能直接沿用这里的矩阵。

由接收机采用的时频协方差 $R_{\rm ass}(q,q')$，一次性抽取

$$
[\mathbf R_{PP}]_{i,j}
=R_{\rm ass}(p_i,p_j),
\qquad
[\mathbf R_{DP}]_{a,j}
=R_{\rm ass}(d_a,p_j),
$$

其中 $p_i,p_j\in P$，$d_a\in D$；$\mathbf R_{DD}$ 是 data-to-data 协方差，$\mathbf R_{PD}=\mathbf R_{DP}^H$，$\mathbf I$ 是维数为 $|P|$ 的单位矩阵。由于每个坐标同时包含 symbol 和 subcarrier，矩阵元素具体为时间相关与频率相关的乘积；两个 DMRS symbol 没有先做平均，也不是“先频域插值、再时域插值”的两级算法。LMMSE 权重和 data-RE 估计为

$$
\mathbf A_{\rm LMMSE}
=
\mathbf R_{DP}
\left[
\mathbf R_{PP}+(\sigma_n^2+\lambda)\mathbf I
\right]^{-1},
\qquad
\widehat{\mathbf g}_{r,D}
=
\mathbf A_{\rm LMMSE}\mathbf z_{r,P}.
$$

$\lambda\ge0$ 是仅在数值求解需要时增加的 diagonal loading；实现不显式求逆，而是对括号内矩阵做 Cholesky 分解和两次三角方程求解。权重按 candidate 和噪声方差缓存。相同的 $\mathbf A_{\rm LMMSE}$ 分别作用于每个 $b,r$ 的 $\mathbf z_{r,P}$，Rx 之间不交换观测。

这里“全部 RE 的估计”需要精确定义：当前 2-symbol PDCCH 正式路径把目标集合设为全部 data RE，即直接输出 $\widehat{\mathbf g}_{r,D}$；它不额外生成包含 DMRS RE 在内的完整占用网格估计张量。如果以后确实需要任意目标集合 $A$ 上的估计，只需用 $\mathbf R_{AP}$ 替换 $\mathbf R_{DP}$，公式本身不变。

当接收机协方差与真实协方差匹配时，data 集合上的理论误差协方差为

$$
\mathbf R_{e,D}
=
\mathbf R_{DD}
-
\mathbf R_{DP}
(\mathbf R_{PP}+\sigma_n^2\mathbf I)^{-1}
\mathbf R_{PD}.
$$

transparent 模式属于协方差失配；此时上述 matched 简式不能代表真实误差，正式仿真仍用生成的真实 $g$ 与实际 $\widehat g$ 逐 trial 统计 NMSE。

#### 4. Rx 合并、QPSK 软解调与译码

对每个 data RE，平台先在 Rx 维做 coherent MRC。定义

$$
D_{b,q}
=
\max\!\left(
\sum_{r=0}^{N_r-1}|\widehat g_{b,r,q}|^2,
10^{-12}
\right),
$$

则均衡后的符号及传给 QPSK demapper 的有效噪声方差为

$$
\widetilde x_{b,q}
=
\frac{
\sum_{r=0}^{N_r-1}\widehat g_{b,r,q}^*y_{b,r,q}
}{D_{b,q}},
\qquad
\widehat\sigma_{b,q}^2
=
\frac{\sigma_n^2}{D_{b,q}}.
$$

estimated CSI 使用 $\widehat g$；ideal CSI 令 $\widehat g=g$ 后走同一 MRC 路径。$10^{-12}$ 只用于防止数值除零。当前 demapper 把估计信道当作精确信道，$\widehat\sigma^2$ 只包含经 MRC 缩放后的 AWGN，没有把 $g-\widehat g$ 的残余信道估计误差另加到 LLR 方差；因此在 CE floor 明显时可能产生过置信 LLR，这也是解释高 SNR BLER floor 时必须保留的实现边界。

QPSK 映射为

$$
x(b_0,b_1)
=
\frac{(1-2b_0)+j(1-2b_1)}{\sqrt2}.
$$

代码采用 $L_i=\log[P(b_i=1\mid\widetilde x)/P(b_i=0\mid\widetilde x)]$ 的译码器符号约定，因此

$$
L_0
=
-\frac{2\sqrt2\,\Re\{\widetilde x\}}{\widehat\sigma^2},
\qquad
L_1
=
-\frac{2\sqrt2\,\Im\{\widetilde x\}}{\widehat\sigma^2}.
$$

LLR 随后裁剪到配置的 `[-llr_clip,+llr_clip]`，按原 coded-bit 顺序展开并送入 PDCCH 解码链，包括 scrambling 恢复、rate recovery、Polar SCL 解码和 CRC24C/RNTI 检查。一个 trial 只要 CRC 失败或恢复出的 payload 任一 bit 与发送 payload 不同，就计为一个 block error；BLER 为 block errors 除以 trials。

#### 5. NMSE 的逐 trial 与曲线统计

对第 $b$ 个 trial，平台把两个 Rx 和全部 data RE 一起做 ratio-of-sums：

$$
\operatorname{NMSE}_b
=
\frac{
\sum_{r=0}^{N_r-1}\sum_{q\in D}
|\widehat g_{b,r,q}-g_{b,r,q}|^2
}{
\max\!\left(
\sum_{r=0}^{N_r-1}\sum_{q\in D}|g_{b,r,q}|^2,
10^{-30}
\right)
}.
$$

$10^{-30}$ 仅防止极端数值下分母为零。每个 SNR 曲线点先在线性域对全部 trial 求均值，再转成 dB：

$$
\operatorname{NMSE}_{\rm curve,dB}
=
10\log_{10}\!\left[
\frac{1}{B}\sum_{b=1}^{B}\operatorname{NMSE}_b
\right],
$$

其中 $B$ 是该 SNR 点实际完成的 trial 数。不得对每个 $\operatorname{NMSE}_b$ 先取 dB 再平均，因为两种统计量不同。ideal CSI 的 raw $\operatorname{NMSE}_b$ 精确为 0；CSV 写出 dB 字段时的 `max(mean,1e-30)` 保护会对应 -300 dB，但这只是数值占位，不是有限的物理估计误差，也不绘制为 ideal-CSI NMSE 性能曲线。

增加独立 Rx 分支会增加 NMSE 分子、分母中的同分布样本数，通常降低该 ratio 统计量的 Monte Carlo 波动，但不会改变每个分支使用的协方差或 LMMSE 权重。因此 2Rx 的主要端到端收益应从 MRC 后的 BLER 判断，不能把 $\Delta_{\rm Rx}$ 直接解释为信道估计器获得相同 dB 的改善。

## 2026-09-08 新增规划：不同 Tx 数和 AL 下的严格 Sidon 候选（待执行）

### 目的与范围

本节研究同一 Tx 数和 AL 下，不同严格 Sidon 集合的 estimated-CSI BLER 是否不同，并用 NMSE 曲线辅助判断差异是否与信道估计质量有关。场景直接复用本 plan 已定义的两组配置：8Tx 使用 `A100` 的 AL2/4/8，4Tx 使用 `C300` 的 AL1/2/4。每个 AL 按自己的 $K$ 和 $N_p$ 单独生成整数索引，不跨 AL 复用同一物理时延或同一整数集合。

本节只能得到每个既有系统设置下、表列候选中的最优严格 Sidon。4Tx 与 8Tx 两组还同时改变了信道、CORESET duration 和移动性，因此不能把跨组差异单独归因于 Tx 数。

### 已完成的硬约束、生成与精确去重

表列候选已经完成以下离线组合计算；这些计算不是信道仿真，也不使用 CE、NMSE 或 BLER：

1. 对 $N_t$ 个整数索引，检查全部 $N_t(N_t+1)/2$ 个无序二元和在模 $K$ 下两两不同；4Tx 和 8Tx 分别检查 10 和 36 个二元和。
2. 检查 $j_n\bmod N_p$ 两两不同。当前 comb-4 DMRS 下，这等价于各 Tx 占用不同的导频折叠位置；相应 pilot matrix 满列秩且条件数为 1。
3. 对天线排列先排序，再在所有公共模 $K$ 循环平移中取字典序最小代表，删除排列和平移等价重复。
4. 不删除反射 $j\mapsto-j$、模 $K$ 单位乘法或仅有相同 pair/fold gap 的集合，因为这些变换未证明在实际 PDP、有限 PDCCH RE 和 matched LMMSE 下保持相同 NMSE 与 BLER。

4Tx 对固定一个元素为 0 的组合空间做完整枚举；8Tx 完整枚举不可执行，使用 seed `20260908` 的均匀随机整数集合生成固定有效池。硬筛和精确去重后的数量如下：

| 场景 | 候选生成范围 | 满足硬约束、尚未平移去重 | 精确去重后的有效池 | proposal 数 |
|---|---:|---:|---:|---:|
| 4Tx AL1，$K=36,N_p=9$ | 完整枚举 | 1,920 | 480 | 6,545 |
| 4Tx AL2，$K=72,N_p=18$ | 完整枚举 | 34,400 | 8,600 | 57,155 |
| 4Tx AL4，$K=144,N_p=36$ | 完整枚举 | 375,072 | 93,768 | 477,191 |
| 8Tx AL2，$K=144,N_p=36$ | 固定随机池 | — | 20,000 | 1,620,000 |
| 8Tx AL4，$K=288,N_p=72$ | 固定随机池 | — | 20,000 | 140,000 |
| 8Tx AL8，$K=576,N_p=144$ | 固定随机池 | — | 20,000 | 60,000 |

仅满足严格 Sidon 和导频可辨识性仍会留下数百至数万个集合，不能全部运行正式 BLER。为得到固定且可直接仿真的小表，每个场景保留 8 个几何覆盖代表：历史 `S0_SIDON`、最大 pair-gap、最大 fold-gap、最小/最大圆周覆盖跨度、最低/最高低频 lag array-factor 能量，并在角色重合时用规范化几何特征的确定性最远点补足。低频 lag 能量只使用 $1\le d<N_p$ 的 $|A_\tau(d)|^2$，不读取物理 PDP。该步骤是试验设计抽样，不是性能筛选，也不声称被删除候选的 BLER 更差。

生成、校验和去重由 `tools/prepare_plan031_sidon_shortlist.py` 固定重放，完整回执为 `outputs/experiment031_pdcch_cdd/20260908_strict_sidon_search/sidon_shortlist.json`。以下表格是后续仿真的冻结输入；`pair/fold` 均以整数栅格表示，`span` 是规范化集合的最小圆周覆盖跨度。

### 8Tx AL2：`K=144`，`N_p=36`

| candidate | $\mathbf j$ | $\boldsymbol\tau$ / ns | pair/fold gap | span |
|---|---|---|---:|---:|
| `A100_SIDON8_AL2_01`；`S0_SIDON` | `[0,1,3,7,12,20,30,65]` | `[0.000,231.481,694.444,1620.370,2777.778,4629.630,6944.444,15046.296]` | `1/1` | 65 |
| `A100_SIDON8_AL2_02` | `[0,3,7,25,46,57,124,139]` | `[0.000,694.444,1620.370,5787.037,10648.148,13194.444,28703.704,32175.926]` | `1/3` | 77 |
| `A100_SIDON8_AL2_03` | `[0,2,6,14,25,32,47,141]` | `[0.000,462.963,1388.889,3240.741,5787.037,7407.407,10879.630,32638.889]` | `1/1` | 50 |
| `A100_SIDON8_AL2_04` | `[0,7,30,49,71,92,116,136]` | `[0.000,1620.370,6944.444,11342.593,16435.185,21296.296,26851.852,31481.481]` | `1/1` | 120 |
| `A100_SIDON8_AL2_05` | `[0,3,73,84,101,113,121,128]` | `[0.000,694.444,16898.148,19444.444,23379.630,26157.407,28009.259,29629.630]` | `1/1` | 74 |
| `A100_SIDON8_AL2_06` | `[0,1,29,74,90,123,127,132]` | `[0.000,231.481,6712.963,17129.630,20833.333,28472.222,29398.148,30555.556]` | `1/1` | 99 |
| `A100_SIDON8_AL2_07` | `[0,2,49,58,84,92,124,129]` | `[0.000,462.963,11342.593,13425.926,19444.444,21296.296,28703.704,29861.111]` | `1/1` | 97 |
| `A100_SIDON8_AL2_08` | `[0,4,11,31,44,54,68,129]` | `[0.000,925.926,2546.296,7175.926,10185.185,12500.000,15740.741,29861.111]` | `1/1` | 83 |

### 8Tx AL4：`K=288`，`N_p=72`

| candidate | $\mathbf j$ | $\boldsymbol\tau$ / ns | pair/fold gap | span |
|---|---|---|---:|---:|
| `A100_SIDON8_AL4_01`；`S0_SIDON` | `[0,1,3,7,12,20,30,65]` | `[0.000,115.741,347.222,810.185,1388.889,2314.815,3472.222,7523.148]` | `1/1` | 65 |
| `A100_SIDON8_AL4_02` | `[0,4,20,33,96,155,203,253]` | `[0.000,462.963,2314.815,3819.444,11111.111,17939.815,23495.370,29282.407]` | `2/4` | 225 |
| `A100_SIDON8_AL4_03` | `[0,7,54,64,89,181,243,262]` | `[0.000,810.185,6250.000,7407.407,10300.926,20949.074,28125.000,30324.074]` | `1/7` | 196 |
| `A100_SIDON8_AL4_04` | `[0,23,53,99,138,183,215,241]` | `[0.000,2662.037,6134.259,11458.333,15972.222,21180.556,24884.259,27893.519]` | `1/1` | 241 |
| `A100_SIDON8_AL4_05` | `[0,3,11,145,183,252,272,284]` | `[0.000,347.222,1273.148,16782.407,21180.556,29166.667,31481.481,32870.370]` | `1/1` | 154 |
| `A100_SIDON8_AL4_06` | `[0,1,9,14,25,42,90,252]` | `[0.000,115.741,1041.667,1620.370,2893.519,4861.111,10416.667,29166.667]` | `1/1` | 126 |
| `A100_SIDON8_AL4_07` | `[0,1,14,48,107,135,152,219]` | `[0.000,115.741,1620.370,5555.556,12384.259,15625.000,17592.593,25347.222]` | `1/1` | 219 |
| `A100_SIDON8_AL4_08` | `[0,4,48,63,116,152,195,280]` | `[0.000,462.963,5555.556,7291.667,13425.926,17592.593,22569.444,32407.407]` | `1/1` | 203 |

### 8Tx AL8：`K=576`，`N_p=144`

| candidate | $\mathbf j$ | $\boldsymbol\tau$ / ns | pair/fold gap | span |
|---|---|---|---:|---:|
| `A100_SIDON8_AL8_01`；`S0_SIDON` | `[0,1,3,7,12,20,30,65]` | `[0.000,57.870,173.611,405.093,694.444,1157.407,1736.111,3761.574]` | `1/1` | 65 |
| `A100_SIDON8_AL8_02` | `[0,5,124,136,163,403,487,505]` | `[0.000,289.352,7175.926,7870.370,9432.870,23321.759,28182.870,29224.537]` | `5/5` | 336 |
| `A100_SIDON8_AL8_03` | `[0,15,46,100,173,208,516,548]` | `[0.000,868.056,2662.037,5787.037,10011.574,12037.037,29861.111,31712.963]` | `1/14` | 268 |
| `A100_SIDON8_AL8_04` | `[0,44,135,218,292,379,458,528]` | `[0.000,2546.296,7812.500,12615.741,16898.148,21932.870,26504.630,30555.556]` | `1/4` | 485 |
| `A100_SIDON8_AL8_05` | `[0,3,15,23,70,104,139,569]` | `[0.000,173.611,868.056,1331.019,4050.926,6018.519,8043.981,32928.241]` | `1/2` | 146 |
| `A100_SIDON8_AL8_06` | `[0,1,5,160,178,241,524,568]` | `[0.000,57.870,289.352,9259.259,10300.926,13946.759,30324.074,32870.370]` | `1/1` | 293 |
| `A100_SIDON8_AL8_07` | `[0,13,215,259,273,373,406,449]` | `[0.000,752.315,12442.130,14988.426,15798.611,21585.648,23495.370,25983.796]` | `1/3` | 374 |
| `A100_SIDON8_AL8_08` | `[0,1,34,276,297,395,415,476]` | `[0.000,57.870,1967.593,15972.222,17187.500,22858.796,24016.204,27546.296]` | `1/1` | 334 |

### 4Tx AL1：`K=36`，`N_p=9`

| candidate | $\mathbf j$ | $\boldsymbol\tau$ / ns | pair/fold gap | span |
|---|---|---|---:|---:|
| `C300_SIDON4_AL1_01`；`S0_SIDON` | `[0,1,3,7]` | `[0.000,925.926,2777.778,6481.481]` | `1/1` | 7 |
| `C300_SIDON4_AL1_02` | `[0,2,6,14]` | `[0.000,1851.852,5555.556,12962.963]` | `2/1` | 14 |
| `C300_SIDON4_AL1_03` | `[0,2,5,16]` | `[0.000,1851.852,4629.630,14814.815]` | `1/2` | 16 |
| `C300_SIDON4_AL1_04` | `[0,1,4,6]` | `[0.000,925.926,3703.704,5555.556]` | `1/1` | 6 |
| `C300_SIDON4_AL1_05` | `[0,7,15,26]` | `[0.000,6481.481,13888.889,24074.074]` | `1/1` | 25 |
| `C300_SIDON4_AL1_06` | `[0,3,7,15]` | `[0.000,2777.778,6481.481,13888.889]` | `1/1` | 15 |
| `C300_SIDON4_AL1_07` | `[0,1,6,23]` | `[0.000,925.926,5555.556,21296.296]` | `1/1` | 19 |
| `C300_SIDON4_AL1_08` | `[0,6,14,29]` | `[0.000,5555.556,12962.963,26851.852]` | `1/1` | 21 |

### 4Tx AL2：`K=72`，`N_p=18`

| candidate | $\mathbf j$ | $\boldsymbol\tau$ / ns | pair/fold gap | span |
|---|---|---|---:|---:|
| `C300_SIDON4_AL2_01`；`S0_SIDON` | `[0,1,3,7]` | `[0.000,462.963,1388.889,3240.741]` | `1/1` | 7 |
| `C300_SIDON4_AL2_02` | `[0,5,16,49]` | `[0.000,2314.815,7407.407,22685.185]` | `5/2` | 39 |
| `C300_SIDON4_AL2_03` | `[0,4,10,32]` | `[0.000,1851.852,4629.630,14814.815]` | `2/4` | 32 |
| `C300_SIDON4_AL2_04` | `[0,1,4,6]` | `[0.000,462.963,1851.852,2777.778]` | `1/1` | 6 |
| `C300_SIDON4_AL2_05` | `[0,16,33,53]` | `[0.000,7407.407,15277.778,24537.037]` | `1/1` | 52 |
| `C300_SIDON4_AL2_06` | `[0,3,7,64]` | `[0.000,1388.889,3240.741,29629.630]` | `1/3` | 15 |
| `C300_SIDON4_AL2_07` | `[0,1,59,68]` | `[0.000,462.963,27314.815,31481.481]` | `1/1` | 14 |
| `C300_SIDON4_AL2_08` | `[0,11,23,48]` | `[0.000,5092.593,10648.148,22222.222]` | `1/1` | 47 |

### 4Tx AL4：`K=144`，`N_p=36`

| candidate | $\mathbf j$ | $\boldsymbol\tau$ / ns | pair/fold gap | span |
|---|---|---|---:|---:|
| `C300_SIDON4_AL4_01`；`S0_SIDON` | `[0,1,3,7]` | `[0.000,231.481,694.444,1620.370]` | `1/1` | 7 |
| `C300_SIDON4_AL4_02` | `[0,11,44,66]` | `[0.000,2546.296,10185.185,15277.778]` | `11/3` | 66 |
| `C300_SIDON4_AL4_03` | `[0,8,20,64]` | `[0.000,1851.852,4629.630,14814.815]` | `4/8` | 64 |
| `C300_SIDON4_AL4_04` | `[0,1,4,6]` | `[0.000,231.481,925.926,1388.889]` | `1/1` | 6 |
| `C300_SIDON4_AL4_05` | `[0,34,69,107]` | `[0.000,7870.370,15972.222,24768.519]` | `1/1` | 106 |
| `C300_SIDON4_AL4_06` | `[0,3,7,136]` | `[0.000,694.444,1620.370,31481.481]` | `1/3` | 15 |
| `C300_SIDON4_AL4_07` | `[0,1,131,136]` | `[0.000,231.481,30324.074,31481.481]` | `1/1` | 14 |
| `C300_SIDON4_AL4_08` | `[0,1,47,95]` | `[0.000,231.481,10879.630,21990.741]` | `1/1` | 95 |

### 后续只执行表列候选的仿真

- 不再设置“BLER 前的分层筛选”，也不依据解析 NMSE、CE floor 或相关矩继续淘汰候选。上述 48 个候选全部进入 estimated-CSI BLER 预扫和正式仿真。
- 预扫沿用本 plan 的 300 trials/点、1 dB 步长，只用于为每个候选确定正式 SNR 网格，不删除候选。正式目标附近步长不超过 0.25 dB；每点至少 10,000 trials、累计 200 errors、上限 50,000 trials，10%/1% BLER 未双侧 bracket 时不外推。
- 8Tx 使用 `A100` 已定义的频域 matched LMMSE；4Tx 使用 `C300` 已定义的二维时频 matched LMMSE。所有候选均为 non-transparent，接收机使用候选真实 $\mathbf V$ 构造的等效协方差。
- 每个 BLER SNR 点同时保留 data-RE NMSE 的均值、中位数、10%/90%分位数和 95% 区间，并保留零噪声 CE floor、pilot rank/condition。BLER 和 NMSE 使用同一批 trial，结果图按 `(Tx, AL)` 给出颜色一致、横轴对齐的 BLER 与 NMSE 两个面板。
- 每个 `(Tx, AL)` 以 1% estimated-CSI BLER 所需 SNR 最低为主判据、10% 为次判据；报告相对该场景 `S0_SIDON` 的 SNR 差和 95% 区间。无法统计区分时报告并列最优，不强行选择唯一集合。

新增仿真输出继续写入 `outputs/experiment031_pdcch_cdd/20260908_strict_sidon_search/`。执行前只需把表列候选写入正式 YAML/runner 配置，完成 candidate 数量、完整 $\mathbf j$、ns 换算、严格 Sidon、fold residue、pilot rank 和无噪声解码的 validate；验证通过后冻结配置并开始预扫。结果确认前不更新 `KNOWLEDGE.md`。

### 2026-09-10 预扫后正式范围变更：明确 CE/BLER floor 的候选不做细网格

研究者于 2026-09-10 明确将目标收窄为寻找各既有 `(Tx, AL)` 场景下性能最好的表列严格 Sidon 候选，并授权不再为已经显示明确 CE/BLER floor、无法竞争 1% BLER 主判据的候选运行正式细 SNR 网格。本节取代上文“48 个候选全部进入正式仿真”的要求；48 个候选仍全部完成了冻结的 300 trials/点预扫，只有满足以下全部条件的候选才在正式阶段前停止：

1. 在原预扫及冻结的高 SNR 补扫范围内没有 1% BLER 的相邻原始点双侧 bracket；
2. 观测最低 BLER 点的 Wilson 95% 区间下界仍严格高于 1%；
3. 最高 SNR 点的 data-RE CE NMSE 距解析零噪声 CE floor 不超过 3 dB；
4. 最高 5 个 SNR 点的 BLER 中位数高于该候选的观测最低 BLER，表明曲线已进入平台或回升，而不是仍持续向 1% 下降。

CE floor 只作为辅助证据，不能单独推出 BLER floor；筛除必须同时满足上述 BLER 和 CE 条件。判据固定重放入口为 `tools/analyze_plan031_strict_sidon_prescan.py`，逐候选量化记录和带 SHA-256 的回执分别为 `outputs/experiment031_pdcch_cdd/20260908_strict_sidon_search/screening/candidate_screening.csv` 与 `outputs/experiment031_pdcch_cdd/20260908_strict_sidon_search/screening/screening_receipt.json`。

| candidate | 预扫 SNR / dB | 最低 BLER（SNR / dB） | 该点 Wilson 95% 下界 | 高端 5 点 BLER 中位数 | 零噪声 CE floor / dB | 最高 SNR CE 距 floor / dB |
|---|---:|---:|---:|---:|---:|---:|
| `C300_SIDON4_AL1_05` | 4 至 30 | 0.0967（13） | 0.0682 | 0.2633 | -9.245 | 1.069 |
| `C300_SIDON4_AL1_06` | 4 至 30 | 0.0300（14） | 0.0159 | 0.0900 | -11.618 | 1.592 |
| `C300_SIDON4_AL1_08` | 4 至 30 | 0.0267（20） | 0.0136 | 0.0667 | -11.357 | 0.776 |
| `C300_SIDON4_AL2_05` | 0 至 20 | 0.0533（13） | 0.0331 | 0.1033 | -7.719 | 2.191 |
| `C300_SIDON4_AL4_05` | -4 至 20 | 0.0700（11） | 0.0462 | 0.2000 | -4.565 | 1.395 |

以上 5 个候选不进入正式细网格；其预扫与补扫数据继续作为负结果保存，1% 所需 SNR 记为不可用且不外推。其余 43 个候选进入正式阶段，仍使用每点至少 10,000 trials、累计 200 errors、上限 50,000 trials，以及 10%/1% bracket 宽度不超过 0.25 dB 的原统计口径。该筛选只排除明显无法达到主判据的候选，不利用候选之间的细小预扫 BLER 差异选择优胜者。

### 2026-09-10 正式细网格前增加 3,000-trial 粗确认

研究者确认 300 trials/点不足以比较 1% BLER 附近的候选性能，要求先对保留的 43 个候选运行中等预算粗确认，再判断可能最优的候选。300-trial 预扫只定位 SNR 区域，不用于依据候选之间的细小差异筛选优胜者。

对每个保留候选，分别找到 300-trial 原始曲线中 10% 和 1% BLER 的第一个下降向相邻 1 dB bracket；每个 bracket 固定取低端、中点和高端，即 0.5 dB 间隔的 3 点，两项目标的点取并集。冻结结果为 43 个候选、254 个 candidate/SNR 点；每点固定 3,000 trials，共 762,000 candidate-trials。所有点保存 error flags 和逐 trial data-RE CE NMSE，candidate 仍使用由全局 seed、candidate ID 和固定标签派生的独立随机流。

配置生成入口为 `tools/prepare_plan031_strict_sidon_coarse_configs.py`，配置冻结回执为 `outputs/experiment031_pdcch_cdd/20260908_strict_sidon_search/coarse_confirmation_config_receipt.json`，六份配置为 `configs/pdcch_result031_strict_sidon_{a100_al2,a100_al4,a100_al8,c300_al1,c300_al2,c300_al4}_coarse.yaml`。仿真输出写入 `outputs/experiment031_pdcch_cdd/20260908_strict_sidon_search/coarse_confirmation/`。

粗确认只用于排除明显较差候选和形成“仍可能最优”的候选集合，不作为最终最优声明。每个 `(Tx, AL)` 仍以 1% BLER 所需 SNR 为主、10% 为辅；候选间比较必须报告二项统计不确定性。只有粗确认后仍可能最优的候选才进入不超过 0.25 dB、至少 10,000 trials/点、目标 200 errors 的最终正式仿真。

### 2026-09-10 粗确认完成记录与初步筛选

六个场景的粗确认已完成，共 43 个候选、254 个 candidate/SNR 点、762,000 candidate-trials，每点均为 3,000 trials。逐点 CSV、error flags、逐 trial data-RE CE NMSE、展开配置、验证回执、日志和环境信息保存在 `outputs/experiment031_pdcch_cdd/20260908_strict_sidon_search/coarse_confirmation/`。分析程序逐点复算 error flags 的错误数及 CE NMSE 的线性均值，两者与 CSV 的最大差分别为 0 和 0 dB；日志未发现运行错误。

固定分析入口为 `tools/analyze_plan031_strict_sidon_coarse.py`。对存在原始双侧 bracket 的目标，先用 Jeffreys 平滑后的 BLER 做 trial 数加权递减 isotonic fit，再在局部 `log10(BLER)` 上插值得到门限；二项不确定性由每点独立 Bernoulli bootstrap 4,000 次给出。86 个候选/目标组合中有 15 个未形成双侧 bracket，因此只记录单侧界限且禁止外推；其中 11 个是 1% 目标，另外 4 个是 10% 目标。

“仍可能最优”的初步规则只使用 1% 主目标：候选的 bootstrap `P(best)` 不低于 5%，或其 95% 区间与点估计领先者的 95% 区间重叠；若门限被网格截断，但其单侧界限仍可能优于领先者区间上界，则保守保留。该规则得到以下 18 个候选：

- `a100_al2`：`A100_SIDON8_AL2_02`；其 1% 门限点估计为 4.794 dB，95% 区间为 [4.516, 4.973] dB。
- `a100_al4`：`A100_SIDON8_AL4_01`、`A100_SIDON8_AL4_02`、`A100_SIDON8_AL4_03`、`A100_SIDON8_AL4_06`；领先点估计为候选 01 的 0.702 dB，95% 区间为 [0.435, 0.935] dB。
- `a100_al8`：`A100_SIDON8_AL8_01`、`A100_SIDON8_AL8_02`、`A100_SIDON8_AL8_05`、`A100_SIDON8_AL8_06`；领先点估计为候选 01 的 -2.824 dB，95% 区间为 [-3.009, -2.674] dB。
- `c300_al1`：`C300_SIDON4_AL1_04`；其 1% 门限点估计为 11.824 dB，95% 区间为 [11.461, 11.975] dB。
- `c300_al2`：`C300_SIDON4_AL2_01`、`C300_SIDON4_AL2_03`、`C300_SIDON4_AL2_04`、`C300_SIDON4_AL2_06`；候选 04 的领先点估计为 5.585 dB，95% 区间为 [5.309, 5.786] dB，候选 01 和 06 的 1% 门限均被截断为 `>5 dB`，因该界限仍低于领先者区间上界而保守保留。
- `c300_al4`：`C300_SIDON4_AL4_01`、`C300_SIDON4_AL4_02`、`C300_SIDON4_AL4_03`、`C300_SIDON4_AL4_06`；领先点估计为候选 06 的 0.683 dB，95% 区间为 [0.453, 0.880] dB。

逐候选门限、排序、诊断、样式和六场景 BLER/CE NMSE 图分别保存在 `outputs/experiment031_pdcch_cdd/20260908_strict_sidon_search/coarse_analysis/target_snr_coarse.csv`、`candidate_ranking_coarse.csv`、`diagnostics.csv`、`curve_styles.json` 和 `*_coarse_bler_nmse.png`，带参数和候选集合的回执为 `analysis_receipt.json`。以上只是决定最终细网格范围的探索性筛选，不构成最终最优声明，也不据此更新 `KNOWLEDGE.md`。

### 2026-09-10 最终细扫描配置冻结

上述 18 个候选进入最终细扫描。对于粗确认中已有双侧 bracket 的目标，填满该原始 bracket 内所有 0.25 dB 点；对于被网格截断的目标，在低 SNR 侧向下扩展 0.5 dB，或在高 SNR 侧向上扩展 1.0 dB，并保持 0.25 dB 间隔。两个目标的点取并集。冻结结果为 149 个 candidate/SNR 点；每点至少 10,000 trials，累计达到 200 errors 后停止，最多 50,000 trials，因此总预算范围为 1,490,000 至 7,450,000 candidate-trials。

配置生成入口为 `tools/prepare_plan031_strict_sidon_fine_configs.py`，配置冻结回执为 `outputs/experiment031_pdcch_cdd/20260908_strict_sidon_search/fine_config_receipt.json`，六份配置为 `configs/pdcch_result031_strict_sidon_{a100_al2,a100_al4,a100_al8,c300_al1,c300_al2,c300_al4}_fine.yaml`。细扫描使用独立的 `plan031_strict_sidon_final_fine_v1` 随机流命名空间，不重复使用粗确认样本。输出写入 `outputs/experiment031_pdcch_cdd/20260908_strict_sidon_search/fine/`，编排状态、验证回执和日志写入相邻的 `fine_orchestration/`。六个场景已通过严格 Sidon、fold residue、pilot rank/condition、CE floor、完整几何和无噪声编解码验证；正式仿真仍由研究者手动运行。

### 2026-09-10 最终分析前的边界补点

149 点细扫描完成后，原始数组与 CSV 一致性审计通过，但 36 个候选/目标组合中只有 25 个达到“原始双侧 bracket 且至少 95% bootstrap 重采样仍有双侧 bracket”的正式区间要求。另有 8 个目标的交点紧邻网格边缘、3 个目标没有双侧 bracket；`A100_SIDON8_AL8_01` 的 10% 原始 bracket 还跨越未采样区间 `[-4.5,-3.0]`，宽 1.5 dB。不得用条件化在仍有 bracket 的 bootstrap 子样本所得区间替代正式 95% 区间，也不得跨该 1.5 dB 空隙插值形成正式结论。

因此只为受影响的 10 个候选追加 23 个新 SNR 点，不重跑原 149 点。边缘目标向缺失侧增加 0.25/0.5 dB 两点；1.5 dB 空隙在预期交点相邻处增加 `-4.25 dB`。每点仍为至少 10,000 trials、累计 200 errors、上限 50,000 trials，并沿用 `plan031_strict_sidon_final_fine_v1` 随机流命名空间。配置生成入口为 `tools/prepare_plan031_strict_sidon_fine_supplement.py`，冻结回执为 `outputs/experiment031_pdcch_cdd/20260908_strict_sidon_search/fine_supplement_config_receipt.json`，四份配置为 `configs/pdcch_result031_strict_sidon_{a100_al4,a100_al8,c300_al2,c300_al4}_fine_supplement.yaml`。补点输出独立写入 `fine_supplement/`，分析时与原细扫描按候选和 SNR 合并。配置、验证和补点仿真均已完成：23/23 个 candidate/SNR 点累计 627,700 trials。

### 2026-09-10 最终统计分析完成

最终分析合并原细扫描 149 点与补点 23 点，共 172 个 candidate/SNR 点、2,644,000 trials、92,766 errors。逐 trial error flags 与 CSV 的长度和错误数完全一致，CE NMSE 线性均值复算最大误差为 0 dB，所有 36 个 10%/1% 目标均由宽度不超过 0.25 dB 的原始双侧点夹住，并且至少 95% bootstrap 重采样保持 bracket。统计入口、目标插值、10,000 次 bootstrap、候选间配对比较和最终图表详见 `tools/analyze_plan031_strict_sidon_fine.py` 及 `research/result-031-PDCCH-CDD时延-BLER.md`/`research/result-031-PDCCH-CDD时延-BLER-text.md` 第 9 节；本轮结果仍待研究者确认，因此未更新 `KNOWLEDGE.md`/`GOALS.md`。

## 2026-09-14 补充规划：2-symbol、AL1 的 CDD911 estimated/ideal-CSI BLER（已执行，待确认）

### 目的、范围与既有计划的优先级

本补充只处理 TDL-C 300 ns、4Tx、2-symbol CORESET、AL1 场景，新增 4Tx 人工时延 `[0,0,911,911] ns`，candidate ID 和图例名称均固定为 `CDD911`。交付目标是在既有 1Rx、AL1 的 estimated-CSI 和 ideal-CSI BLER 图上分别补入 CDD911，并额外给出 2Rx 的 CDD911 曲线。每种 CSI 模式各生成一张对比图，每张图严格只含以下四条曲线：

| 曲线 | Rx 数 | 数据动作 |
|---|---:|---|
| `CDD911` | 1 | 新运行 |
| `CDD911` | 2 | 新运行 |
| `C300_S0_SIDON` | 1 | 复用既有正式结果，只重绘 |
| `C300_PRG_DFT4_6REG` | 1 | 复用既有正式结果，只重绘 |

estimated-CSI 与 ideal-CSI 均采用这四条曲线的相同集合。不得重跑 1Rx Sidon 或 1Rx precoder cycling，不运行也不绘制 2Rx Sidon 或 2Rx precoder cycling；B0、AP、GEO、small-CDD 及其他严格 Sidon 候选均不进入本补充的图、表或新仿真。这里的“同一张图”是指同一 CSI 模式内四条曲线共轴；estimated-CSI 与 ideal-CSI 分成两张图，避免把两种接收机定义混在同一坐标图中。

本节是前述“C300 4Tx/2Rx、2-symbol 的 ideal/estimated-CSI BLER 与 NMSE”规划对本次 AL1 交付范围的后续收窄。执行本节时不得调用其中的 AL1 全候选 2Rx 配置，也不得把已经存在的其他 2Rx candidate 数据加入本节结果。该旧节的 AL2/4 状态不由本节改变；本节不授权运行 AL2/4。

### 冻结系统配置与 CDD911 定义

四个新运行组合为 `n_rx in {1,2}` × `channel_estimation in {frequency_lmmse,ideal}`，每个组合只含 `CDD911`。除 Rx 数和 CSI 模式外，物理链路逐字段冻结为：4Tx、单层、48-RB CORESET、2 个相邻 CORESET symbols、AL1、每 symbol 实占 3 RB、$K=36$、每 symbol 9 个 comb-4 DMRS、总 DMRS RE 18、data RE 54、$E=108$、30 kHz SCS、FFT 4096、CP 288、non-interleaved、first CCE 0、`L=6 REG`、A=41、CRC24C/RNTI `0xFFFF`、QPSK、Polar list size 8、TDL-C 300 ns、3 km/h、4 GHz、20 sinusoids、归一化信道，以及 data/DMRS 使用同一预编码。

CDD911 的物理时延和 AL1 相位坐标固定为

$$
\boldsymbol\tau=[0,0,911,911]\ \mathrm{ns},
\qquad
q_K=\frac{1}{36\times30\ \mathrm{kHz}}
=925.925926\ \mathrm{ns},
$$

$$
\mathbf j=\frac{\boldsymbol\tau}{q_K}
=[0,0,0.98388,0.98388].
$$

发射预编码为

$$
W_{k,n}=\frac12\exp\!\left(-j2\pi k j_n/36\right),
$$

其中 $k$ 是 candidate 内从 0 开始的本地唯一子载波索引，phase denominator 固定为 36。配置必须使用 `delay_grid_coordinates: [0,0,0.98388,0.98388]`；生成回执另存精确的 `[0,0,911,911] ns`，并逐项验证由坐标反算的时延。两个 0 和两个 911 ns 是有意的重复时延，不得去重、排序成其他集合或解释为四个互异 CDD 分支。每个 RE 满足 $\sum_n|W_{k,n}|^2=1$，SNR 定义为单位总发射功率相对单个 Rx 分支噪声功率，故每个 Rx 的复高斯噪声方差均为 $10^{-\mathrm{SNR}_{\rm dB}/10}$，不随 `n_rx` 缩放。

CDD911 配置另须显式设置 `allow_duplicate_delays: true`；该开关只放宽默认的时延互异校验，不改变发射、信道或接收机计算，其他候选仍默认拒绝重复时延。

两个 Rx 分支继续使用 `DESIGN.md` 第 3.9 节的独立同分布假设：各 Tx--Rx 信道具有相同 TDL-C 边缘统计，Rx 间不相关且噪声独立，不引入空间相关、分支失配或跨 Rx 联合估计。结果不得外推到相关 Rx 场景。

### estimated-CSI、ideal-CSI 与合并口径

estimated-CSI 的 CDD911 固定为 non-transparent：`receiver.channel_estimation: frequency_lmmse` 且 `receiver_covariance_mode: matched_effective`。接收机知道 CDD911 的真实 $\mathbf V$，以真实 `(symbol,subcarrier)` pilot/data 坐标、TDL-C 时间相关和 CDD911 等效频率协方差构造二维时频 LMMSE；两个 DMRS symbol 不做静态平均。2Rx 时对每个 Rx 独立应用同一个 LMMSE 线性算子，不跨 Rx 联合估计。

ideal-CSI 直接令每个 Rx、每个 data RE 的 $\hat g_r=g_r$，不构造 LMMSE filter；DMRS 保留但不参与估计，CE NMSE 必须逐元素精确为 0。ideal-CSI 下不存在透明或非透明信道估计的差别，但发射波形仍严格使用 CDD911。

两种 CSI 模式都按第 3.9 节公式在 Rx 维 coherent MRC。estimated-CSI 同时保存逐 trial 的线性 data-RE NMSE，按全部 Rx 与全部 data RE 的误差能量之和除以真实信道能量之和，再先在线性域跨 trial 求均值、最后转 dB；ideal-CSI 的零 NMSE 只作回执，不绘有限 dB 曲线。

### 既有曲线的数据复用与可比性审计

estimated-CSI 的历史 1Rx 曲线只从 `outputs/experiment031_pdcch_cdd/20260907_c300_4tx_2sym/analysis/formal_points.csv` 读取 `C300_S0_SIDON` 和 `C300_PRG_DFT4_6REG` 的 AL1 行；ideal-CSI 只从 `outputs/experiment031_pdcch_cdd/20260910_c300_4tx_2sym_ideal_csi/analysis/formal_points.csv` 读取同名两个 candidate 的 AL1 行。分析器必须保留 `source_csv` 指针，并联合对应的 `resolved_config.yaml` 和 `run_metadata.json` 核对 candidate、AL、`n_tx=4`、`n_rx=1`、CSI 模式、资源计数、信道配置、SNR 定义和功率归一化；任一关键字段不一致即停止合图。历史点及其 trial flags 不复制成新的仿真结果，也不与新 trial 合并计数。

历史 1Rx 曲线与新 CDD911 使用独立随机样本；四个 CDD911 新运行组合也使用互异的固定随机流 namespace。因此所有跨曲线目标 SNR 差均按独立样本 bootstrap，不声明 paired-sample 方差缩减。历史数据只重绘，不因图形整合而改变其既有 BLER、errors、trials 或 Wilson 区间。

### 固定 SNR 网格、正式统计与停止条件

本补充不做 prescan，也不根据 CDD911 的初步结果选择或追加 SNR 点。每个 Rx/CSI 组合直接复用对应既有 AL1 正式配置中 Sidon 与 precoder cycling 两条曲线的 SNR 点并集；这样既不重新设计网格，也覆盖两条保留基线原先使用的工作区间。记 `{a:0.25:b}` 为从 $a$ 到 $b$、包含两端的 0.25 dB 等间隔点，四组冻结网格为：

| CSI | Rx | CDD911 固定 SNR / dB | 来源配置 |
|---|---:|---|---|
| estimated | 1 | `{8:0.25:9} ∪ {10:0.25:11} ∪ {13:0.25:14} ∪ {16.5:0.25:17.5}` | `configs/pdcch_result031_c300_al1_formal.yaml` |
| estimated | 2 | `{2.75:0.25:6.5} ∪ {8.25:0.25:10}` | `configs/pdcch_result031_c300_2rx_estimated_al1_formal.yaml` |
| ideal | 1 | `{5.75:0.25:7.5} ∪ {8.75:0.25:10.75} ∪ {14.5:0.25:16.25}` | `configs/pdcch_result031_c300_ideal_al1_formal.yaml` |
| ideal | 2 | `{0:0.25:5} ∪ {7.75:0.25:9.5}` | `configs/pdcch_result031_c300_2rx_ideal_al1_formal.yaml` |

新增四份正式配置，每份只允许包含 `CDD911`。统一 seed 为 `20260914`，并按 Rx 数、CSI 模式和 candidate ID 派生稳定且互异的随机流；输出写入 `outputs/experiment031_pdcch_cdd/20260914_c300_al1_cdd911/`。四份配置分别复制上表对应来源配置的系统字段、trial 预算和 batch size，只修改 candidate、输出目录、seed/namespace，并按表中并集设置 SNR。配置完成后调用 runner 的 `--stage validate`，核对 CDD911 的 `[0,0,911,911] ns` 回算、phase denominator、重复时延未被去重、资源计数、接收机模式和 1Rx/2Rx MRC 路径；validate 通过后直接正式运行，不另建 smoke 或 prescan 配置。

正式统计完全沿用既有 AL1 曲线：每点至少 10,000 trials，累计达到 200 errors 后可停止，最多 50,000 trials。保存 errors/trials、Wilson 95% 区间、逐 trial error flags、展开配置、日志和环境信息；estimated-CSI 还保存逐 trial 线性 NMSE。中断恢复必须从绝对 trial 区间继续，禁止从 trial 1 重跑后相加。

目标 SNR 使用预先规定的递减 isotonic fit 加目标邻域 `log10(BLER)` 局部插值。只有固定网格中的原始正式点对 10%/1% 形成双侧 bracket 时才报告目标值和 bootstrap 95% 区间；若 CDD911 的目标落在上述既有网格空隙或范围之外，则明确报告“未取得”或单侧界限，不外推，也不自动补点。是否再增加 SNR 点由研究者看过首轮结果后另行决定，不属于本补充的默认执行范围。

### 实现范围、测试与复现入口

复用 `cdd_lls/sim/pdcch_cdd.py` 和 `tools/run_plan031_candidates.py` 的多 Rx、二维 LMMSE、ideal CSI、MRC 与可恢复运行实现，不复制链路主循环。最低限度只新增：

- `configs/pdcch_result031_cdd911_{1rx,2rx}_{estimated,ideal}_formal.yaml` 四份单 candidate 配置；
- `tools/analyze_plan031_cdd911.py`，用于读取四条 CDD911 新曲线和按 CSI 模式复用的历史 1Rx Sidon/cycling 曲线，完成来源核对、目标统计和两张图。

若现有公共实现已经满足要求，不修改 `cdd_lls/`，也不新增专用 prepare 脚本或专用测试文件。配置的精确时延换算、重复时延、matched-effective 接收机和四曲线白名单由 `--stage validate` 与分析器的启动检查完成；仅当公共实现确有缺口时才做最小修改并补相应定向测试。由于本场景是动态 TDL-C、4Tx、QPSK PDCCH，超出 `tools/run_bler_curves.py` 文档规定的 static TDL-A、8Tx、16QAM 固定范围，本补充不得直接套用该入口。

建议复现顺序为：

```powershell
python tools/run_plan031_candidates.py --config <四份 formal YAML> --stage validate
python tools/run_plan031_candidates.py --config <四份 formal YAML> --stage run --max-workers 4
python tools/analyze_plan031_cdd911.py --bootstrap-repeats 4000
python -m pytest tests/test_pdcch.py tests/test_rmmse_time_frequency.py tests/test_plan031_analysis.py -q
```

尖括号表示上述四份明确命名的 YAML，应在实际命令中逐一展开，不是可原样传给 shell 的字面参数。正式执行前只要求配置审阅、定向测试和 validate 通过。

### 输出、图表、结果问题与验收

输出目录至少保存展开 YAML、环境信息、测试与运行日志、四个组合的 `bler_points.csv`、逐 trial error flags、estimated-CSI NMSE 数组、`formal_points.csv`、`target_snr.csv`、`target_differences.csv`、`diagnostics.csv`、`analysis_metadata.json` 和历史数据来源回执。所有新 CSV 必须显式记录 `candidate_id=CDD911`、`delay_ns=[0,0,911,911]`、`delay_grid_coordinates=[0,0,0.98388,0.98388]`、`n_tx=4`、`n_rx`、CSI 模式、AL1、`matched_effective` 或 ideal 标签、SNR 定义和单位总发射功率标签。

最终只生成以下两张主图并写入 `docs/figures/result-031/`：

1. AL1 estimated-CSI BLER：1Rx CDD911、2Rx CDD911、1Rx Sidon、1Rx precoder cycling；
2. AL1 ideal-CSI BLER：1Rx CDD911、2Rx CDD911、1Rx Sidon、1Rx precoder cycling。

两图均使用对数 BLER 纵轴，显示原始正式点、Wilson 95% 误差棒及 10%/1% 参考线；颜色区分 waveform，线型或 marker 区分 Rx 数，图例显式写 `1Rx`/`2Rx` 和 CSI 模式。不得用插值曲线遮蔽原始非单调点。横轴范围取四条正式曲线 SNR 并集，不能因裁剪隐藏任一已有 bracket。

在 `research/result-031-PDCCH-CDD时延-BLER.md` 与 `research/result-031-PDCCH-CDD时延-BLER-text.md` 增加同名配对章节，必须回答：CDD911 在 1Rx 下相对历史 Sidon 与 precoder cycling 的 10%/1% SNR 差及独立 bootstrap 95% 区间；CDD911 从 1Rx 到 2Rx 的 $\Delta_{\rm Rx}=\mathrm{SNR}_{1Rx}-\mathrm{SNR}_{2Rx}$ 及区间；estimated 与 ideal 的 $\Delta_{\rm CE}=\mathrm{SNR}_{estimated}-\mathrm{SNR}_{ideal}$ 及区间；estimated-CSI 的 CE NMSE、解析零噪声 CE floor、pilot rank/condition 与重复时延造成的可辨识性限制。$\Delta_{\rm Rx}>0$ 只表示当前独立双 Rx MRC 模型下降低了所需输入 SNR，不得解释为信道估计器本身获得同等 dB 改善；ideal/estimated 差值也不得归因于单一的时间或频率插值因素。

验收要求为：四个新组合在上述固定网格上均完成预算或明确达到上限；两张图的曲线白名单和历史数据来源审计通过；CDD911 的 10%/1% 目标按固定网格如实报告，未闭合时不把补点作为本轮完成前提；trial flags、errors、trials、Wilson 区间和 estimated-CSI NMSE 线性均值复算一致；图文两版 result 的数字、限制和结论一致。结果仍属于 plan-031 的待确认补充，研究者确认前不更新 `KNOWLEDGE.md` 或 `GOALS.md`，也不创建 Git checkpoint。

### 2026-09-14 首轮分析后定向补点（已执行）

首轮固定网格已全部完成，但 CDD911 1Rx 的 10%/1% 和两种 Rx 的 1% 交点落在既有网格空隙中。研究者确认不做新预扫，直接为这些目标补充以下 20 个 0.25 dB 局部点：

| CSI | Rx | 新增 SNR / dB | 目标 |
|---|---:|---|---|
| estimated | 1 | `9.25, 9.50, 9.75`；`14.25, 14.50, 14.75` | 10%；1% |
| estimated | 2 | `6.75, 7.00, 7.25` | 1% |
| ideal | 1 | `7.75, 8.00, 8.25, 8.50`；`12.50, 12.75, 13.00, 13.25` | 10%；1% |
| ideal | 2 | `5.25, 5.50, 5.75` | 1% |

四份既有 formal YAML 直接加入表列点，保持原 output directory、seed、random-stream namespace、candidate 和 trial 预算不变。`resume: true` 必须跳过已完成点，只运行新增 SNR；中断后重复同一命令继续。新增点仍为至少 10,000 trials、200 errors 后停止、上限 50,000，并保存逐 trial flags 和 estimated-CSI NMSE。补点完成后重新运行 `tools/analyze_plan031_cdd911.py`，正式目标仍要求不宽于 0.25 dB 的原始双侧 bracket 和至少 95% 的 4,000 次 bootstrap replicate 有效。

为展示历史 Sidon 与 precoder cycling 的完整已有曲线，最终图允许在高预算 formal 点之外叠加其原 prescan 点：estimated 数据来自 `outputs/experiment031_pdcch_cdd/20260907_c300_4tx_2sym/prescan/al1/`，ideal 数据来自 `outputs/experiment031_pdcch_cdd/20260910_c300_4tx_2sym_ideal_csi/prescan/al1/`。prescan 每点只有 300 trials，必须使用较淡线条并在图注中明确标识；它们只用于显示宽 SNR 趋势，不进入目标 SNR、差值、置信区间或正式样本审计。

## 2026-09-14 补充规划：CDD911/CDD130 transparent estimated-CSI BLER（已执行，待确认）

### 目的与范围

本补充继续使用第 11 节的 TDL-C 300 ns、4Tx、2-symbol CORESET、AL1 场景，只新增 transparent estimated-CSI 接收结果，不新增 ideal-CSI 仿真。需要运行以下四条新曲线：`CDD911_transparent` 的 1Rx/2Rx，以及 `CDD130_transparent` 的 1Rx/2Rx。`CDD911_transparent` 与既有 `CDD911` 使用完全相同的 `[0,0,911,911] ns` 发射波形，唯一接收机差异是由 `matched_effective` 改为不知道 CDD 的 `physical_fullband` 协方差；`CDD130_transparent` 使用 `[0,0,130,130] ns` 发射时延并采用同一 transparent 接收口径。

本轮要回答两个问题：同一 CDD911 波形在接收端知道/不知道 CDD 时，1Rx 和 2Rx 的 estimated-CSI BLER 与 CE NMSE 如何变化；在相同 transparent 接收机下，把重复时延对从 911 ns 缩短为 130 ns 后，1Rx 和 2Rx 是否都能在指定网格内观察到不高于 $10^{-2}$ 的原始 BLER 点。本补充不重跑第 11 节既有的 CDD911 matched、CDD911 ideal、Sidon 或 precoder-cycling 数据，也不改变第 11 节首轮补点的执行状态。

### 冻结物理配置、candidate 与 transparent 接收机

除本节明确列出的 candidate、接收机协方差口径、SNR 网格、输出目录和随机流外，系统配置逐字段继承第 11 节：4Tx、单层、AL1、48-RB CORESET、2 个相邻 CORESET symbols、每 symbol 实占 3 RB、$K=36$、18 个 DMRS RE、54 个 data RE、$E=108$、30 kHz SCS、FFT 4096、CP 288、non-interleaved、first CCE 0、`L=6 REG`、A=41、CRC24C/RNTI `0xFFFF`、QPSK、Polar list size 8、TDL-C 300 ns、3 km/h、4 GHz、20 sinusoids、归一化信道，以及 data/DMRS 同预编码。预编码继续满足每个 RE 的总发射功率为 1；每个 Rx 分支的噪声方差为 $10^{-\mathrm{SNR}_{\rm dB}/10}$，不随 Rx 数缩放。2Rx 分支独立同分布，分别估计后按既有公式做 coherent MRC，不引入 Rx 相关或跨 Rx 联合估计。

两条 candidate 冻结为：

| candidate ID/图例 | 物理人工时延 / ns | `delay_grid_coordinates` | `phase_denominator` | `receiver_covariance_mode` |
|---|---|---|---:|---|
| `CDD911_transparent` | `[0,0,911,911]` | `[0,0,0.98388,0.98388]` | 36 | `physical_fullband` |
| `CDD130_transparent` | `[0,0,130,130]` | `[0,0,0.1404,0.1404]` | 36 | `physical_fullband` |

两者均使用

$$
W_{k,n}=\frac12\exp\!\left(-j2\pi k j_n/36\right),
$$

其中 $k$ 为 candidate 内从 0 开始的本地唯一子载波索引，$q_K=925.925926$ ns，故 130 ns 对应精确坐标 $130/q_K=0.1404$。两条配置都必须显式设置 `allow_duplicate_delays: true`；重复的 0 ns 和重复的非零时延是预定波形的一部分，不得去重或改写。

全局 `receiver.channel_estimation` 保持 `frequency_lmmse`。`physical_fullband` 表示 UE 不知道 CDD delay 或 $\mathbf V$，仅使用底层 TDL-C 物理时频协方差在整个 AL1 candidate 占用带宽构造二维 LMMSE；两个 DMRS symbol 使用真实时间坐标，不做静态平均。该口径与 `C300_SMALL_CDD_QSTEP0P25_TRANSPARENT_CDD` 的 full-band transparent 定义一致。它是相对于真实等效协方差的失配估计，不能在分析时改用 candidate 的人工时延修正协方差。

### 固定 SNR 网格、统计预算与 $10^{-2}$ 判据

两条 candidate 使用相同的逐 Rx 固定网格，不做 prescan，也不以首批结果删点：

| Rx | 固定 SNR / dB | 每条 candidate 点数 |
|---:|---|---:|
| 1 | `10, 11, 12, 14, 16, 16.5, 17, 17.5, 18, 18.5, 19, 19.5` | 12 |
| 2 | `4, 5, 6, 7, 7.5, 8, 8.5, 9, 9.5, 10` | 10 |

合计新增 44 个 candidate/SNR 点。沿用既有 CDD911 estimated 配置的资源预算，1Rx/2Rx 的 `batch_size` 分别固定为 100/50；每点至少运行 10,000 trials，累计达到 200 errors 后可停止，最多 50,000 trials；保存逐 trial error flags、errors/trials、Wilson 95% 区间和逐 trial 线性 data-RE CE NMSE。中断恢复必须沿用相同 seed 和 random-stream namespace，从绝对 trial 区间继续，禁止从 trial 1 重跑后与旧数据相加。

“至少看到 $10^{-2}$ BLER”在本节中定义为：四条新曲线中的每一条都至少有一个完成预算的原始正式点满足 `errors / trials <= 0.01`；插值值、isotonic 拟合值或 Wilson 区间越过 0.01 均不能替代这个原始点判据。所有点仍须显示 Wilson 区间。若任一曲线在表列最高 SNR 完成 50,000 trials 后仍没有原始 BLER 不高于 0.01，则本补充的该项验收记为未满足，保留全部负结果且不得外推；由于研究者只冻结了上述 SNR 点，本 plan 不自动增加表外点，后续扩展须先补充明确的 SNR 与上限。

10%/1% 目标 SNR 仍使用第 11 节的递减 isotonic fit、局部 `log10(BLER)` 插值和 4,000 次独立 Bernoulli bootstrap，但只有相邻原始正式点形成不宽于 0.25 dB 的双侧 bracket，且至少 95% bootstrap replicate 保持 bracket 时才作为合格目标报告。由于本节网格按曲线观察而不是精确门限定位设计，宽于 0.25 dB 的跨点区间只能在 `target_snr.csv` 中记录 raw bracket/不合格状态，不得给出正式插值门限。matched 与 transparent、CDD911 与 CDD130、1Rx 与 2Rx 使用独立随机流，所有差值区间按独立样本 bootstrap，不声明 paired-sample 方差缩减。

### 实现、验证与复现入口

新增两份 formal YAML：

- `configs/pdcch_result031_cdd_transparent_1rx_formal.yaml`：只含 `CDD911_transparent`、`CDD130_transparent`，使用上述 1Rx 网格；
- `configs/pdcch_result031_cdd_transparent_2rx_formal.yaml`：只含相同两条 candidate，使用上述 2Rx 网格。

两份配置的全局 seed 固定为 `20260914`，分别使用新的 `plan031-cdd-transparent-1rx-formal-v1` 和 `plan031-cdd-transparent-2rx-formal-v1` namespace；runner 继续按 candidate ID 派生互异随机流。输出写入 `outputs/experiment031_pdcch_cdd/20260914_c300_al1_cdd911/transparent_supplement/{1rx,2rx}/`，不得写入或覆盖既有 `formal/{estimated,ideal}/` shard。

复用 `cdd_lls/sim/pdcch_cdd.py` 与 `tools/run_plan031_candidates.py` 的 duplicate-delay、`physical_fullband`、二维 LMMSE、多 Rx MRC 和可恢复运行实现；不复制主循环。扩展 `tools/analyze_plan031_cdd911.py`，使其在保留既有第 11.2 节 matched/ideal 数据和来源审计的同时读取四条 transparent 新曲线。分析器必须核对每个新 shard 的 expanded/resolved config、candidate ID、精确 delay、`allow_duplicate_delays=true`、`receiver_covariance_mode=physical_fullband`、Rx 数、AL、资源计数、信道、SNR 定义、功率归一化和完整固定网格，任一字段不一致即停止合图。

validate 和定向测试至少覆盖：两个 candidate 的物理时延回算与重复值未被去重；transparent LMMSE 权重与直接由底层物理时频协方差构造的 full-band 权重一致，并与同波形 `matched_effective` 权重不同；一个 Rx 的 pilot 扰动不改变另一 Rx 的估计；每个 RE 总预编码功率为 1；每分支噪声方差不随 `n_rx` 缩放；无噪声链路可解码；输出无 NaN/Inf。若公共实现无需修改，可把配置审计补入现有分析测试，不新增专用算法测试文件。

建议复现顺序为：

```powershell
python tools/run_plan031_candidates.py --config configs/pdcch_result031_cdd_transparent_1rx_formal.yaml --config configs/pdcch_result031_cdd_transparent_2rx_formal.yaml --stage validate
python -m pytest tests/test_pdcch.py tests/test_rmmse_time_frequency.py tests/test_plan031_analysis.py -q
python tools/run_plan031_candidates.py --config configs/pdcch_result031_cdd_transparent_1rx_formal.yaml --config configs/pdcch_result031_cdd_transparent_2rx_formal.yaml --stage run --max-workers 4
python tools/analyze_plan031_cdd911.py --bootstrap-repeats 4000
```

### 图表、result 与完成条件

更新第 11.2 节现有 AL1 estimated-CSI BLER 图，在原四条曲线基础上加入 `CDD911_transparent` 1Rx/2Rx 和 `CDD130_transparent` 1Rx/2Rx，共八条曲线：既有 `CDD911` matched 1Rx/2Rx、历史 Sidon 1Rx、历史 precoder cycling 1Rx，以及本节四条 transparent 曲线。颜色区分 waveform/时延，线型或 marker 同时区分 `matched`/`transparent` 与 1Rx/2Rx；图例必须显式写明 CSI 为 estimated、Rx 数和 transparent 口径。保留原始点、Wilson 95% 误差棒、10%/1% 参考线，并把横轴扩展到全部绘制点的 SNR 并集。第 11.2 节 ideal-CSI 图不增加本节曲线、不改变原数据。

每个新 shard 还必须保存实际展开配置、`resolved_config.yaml`、`run_metadata.json`、环境信息和运行日志。同步更新 `research/result-031-PDCCH-CDD时延-BLER.md` 与 `research/result-031-PDCCH-CDD时延-BLER-text.md` 的第 11 节，并在现有 analysis 目录追加或更新 `formal_points.csv`、`target_snr.csv`、`target_differences.csv`、`same_snr_comparisons.csv`、`diagnostics.csv`、`source_receipt.json`、`curve_styles.json` 和 `analysis_metadata.json`。新行必须显式记录 candidate ID、物理 delay 数组与单位、grid coordinates、`n_tx=4`、`n_rx`、`receiver=estimated`、`receiver_covariance_mode=physical_fullband`、AL1、SNR 定义和单位总发射功率标签。result 必须分别报告四条曲线是否达到原始 $10^{-2}$ 点、其最低原始 BLER/Wilson 区间/错误数/试验数、CE NMSE 与 zero-noise CE floor，并在具备合格 bracket 时报告 CDD911 matched-transparent、CDD911-CDD130 transparent 和 1Rx-2Rx 的目标 SNR 差及 95% 区间；不合格目标只报告 raw bracket 或单侧界限。

完成条件为：44 个新增点全部完成预算或明确达到 50,000 上限；四条新曲线分别通过原始 BLER 不高于 0.01 的判据，否则明确记录未满足；逐 trial flags、errors、trials、Wilson 区间和 CE NMSE 线性均值复算一致；第 11.2 节 estimated 图正好包含上述八条白名单曲线且 ideal 图不变；来源与配置审计、定向测试和 UTF-8/diff 检查通过；两版 result 的数字、限制和结论一致。结果经研究者确认前不更新 `KNOWLEDGE.md` 或 `GOALS.md`，也不创建 Git checkpoint。

执行记录：44/44 个新增点已完成，共 1,224,350 trials、14,938 errors；四条曲线均取得原始 BLER 不高于 0.01 的点。分析器已完成逐 trial、Wilson 区间和 CE NMSE 复算，并生成第 11.2 节八曲线 estimated-CSI 图、保持 ideal-CSI 图为原四条。精确数据、目标状态、来源和脚本哈希见 `outputs/experiment031_pdcch_cdd/20260914_c300_al1_cdd911/analysis/`；结论见两版 `research/result-031` 第 11 节，仍待研究者确认。

## 2026-09-15 补充规划：Sidon transparent、CDD130 ideal 与两类 codebook（部分已执行，待确认）

### 目的、候选与冻结配置

本补充继续使用第 11 节的 TDL-C 300 ns、4Tx、2-symbol CORESET、AL1 场景，共纳入七条曲线：`C300_S0_SIDON_TRANSPARENT` 1Rx estimated-CSI，`CDD130` 1Rx/2Rx ideal-CSI，原 precoder cycling 的 2Rx transparent estimated-CSI 与 ideal-CSI，以及新增固定 LTE 预编码的 2Rx transparent estimated-CSI 与 ideal-CSI。原 precoder cycling 使用 4 维 DFT 矩阵，candidate ID 和图例名称在本补充中固定为 `DFTcodebook`；新增方案命名为 `LTEcodebook`，其 codebook 为 `[1,-1,-1,1]`、`[1,-1,1,-1]`、`[1,1,-1,-1]`、`[1,1,1,1]`，本轮只选择单位范数向量 `[1,-1,1,-1]/2`，不在四个向量之间轮换。Sidon transparent 与既有 `C300_S0_SIDON` 使用完全相同的 `[0,1,3,7]` grid coordinates 和发射波形，只把接收机由 `matched_effective` 改为不知道 CDD 的 `physical_fullband` 二维时频 LMMSE。CDD130 ideal 使用与昨日 `CDD130_transparent` 完全相同的 `[0,0,130,130] ns` 发射时延和 `[0,0,0.1404,0.1404]` grid coordinates，但直接使用每个 data RE 的真实等效信道，不构造 LMMSE filter。

系统参数冻结为：4Tx、单层、AL1、48-RB CORESET、2 个相邻 symbols、每 symbol 3 RB、$K=36$、18 个 DMRS RE、54 个 data RE、$E=108$、30 kHz SCS、FFT 4096、CP 288、non-interleaved、first CCE 0、`L=6 REG`、A=41、CRC24C/RNTI `0xFFFF`、QPSK、Polar list size 8、TDL-C 300 ns、3 km/h、4 GHz、20 sinusoids、归一化信道和单位总发射功率。2Rx 仍为独立同分布分支并做 coherent MRC；单个 Rx 分支噪声方差不随 Rx 数缩放。CDD130 必须设置 `allow_duplicate_delays: true`；Sidon 不设置该放宽开关。两类 codebook 均保持每个 6-REG bundle 内 data/DMRS 同预编码；estimated-CSI 使用 precoder-cycling 既有的 transparent `physical_prg` 口径，ideal-CSI 直接使用每个 data RE 的真实等效信道。AL1 只有一个 6-REG bundle：`DFTcodebook` 固定使用 DFT 索引 `[0]`，即 `[1,1,1,1]/2`；`LTEcodebook` 固定使用 `[1,-1,1,-1]/2`。两者是不同的固定发射波形，本轮不发生 bundle 间 cycling。

### 固定 SNR、预算与输出

| 曲线 | CSI/接收口径 | 固定 SNR / dB | 点数 |
|---|---|---|---:|
| `C300_S0_SIDON_TRANSPARENT`, 1Rx | estimated；`frequency_lmmse` + `physical_fullband` | `10, 11, 12, 15, 16, 17, 18` | 7 |
| `CDD130`, 1Rx | ideal | `10, 11, 12, 15, 16, 17, 18` | 7 |
| `CDD130`, 2Rx | ideal | `5, 6, 7, 7.5, 8, 8.5, 9, 10` | 8 |
| `DFTcodebook`, 2Rx | estimated；transparent `frequency_lmmse` + `physical_prg` | `5, 6, 7, 7.5, 8, 8.5, 9, 10` | 8 |
| `DFTcodebook`, 2Rx | ideal | `5, 6, 7, 7.5, 8, 8.5, 9, 10` | 8 |
| `LTEcodebook`, 2Rx | estimated；transparent `frequency_lmmse` + `physical_prg` | `5, 6, 7, 7.5, 8, 8.5, 9, 10` | 8 |
| `LTEcodebook`, 2Rx | ideal | `5, 6, 7, 7.5, 8, 8.5, 9, 10` | 8 |

冻结网格合计 54 个 candidate/SNR 点。可复用原三条曲线已经完成且落入新网格的点，网格外历史点保留但不计入本补充的冻结点数；缺失点从新的绝对 trial 区间运行，不得重跑后相加。1Rx/2Rx 的 `batch_size` 分别为 100/50；每点至少 10,000 trials，达到 200 errors 后可停止，上限 50,000 trials。保存 errors/trials、Wilson 95% 区间、逐 trial error flags、展开配置、日志和环境信息；所有 estimated-CSI 曲线另存逐 trial 线性 CE NMSE，ideal-CSI 的 CE NMSE 必须逐元素为 0。七条曲线均须至少取得一个完成预算且 `errors/trials <= 0.01` 的原始正式点；若最高 SNR 仍未达到，则保留负结果且不自动增加表外 SNR。

配置沿用 `configs/pdcch_result031_sidon_transparent_1rx_formal.yaml`、`configs/pdcch_result031_cdd130_1rx_ideal_formal.yaml` 和 `configs/pdcch_result031_cdd130_2rx_ideal_formal.yaml`；两类 codebook 的 2Rx estimated/ideal 分别写入 `configs/pdcch_result031_codebooks_2rx_estimated_formal.yaml` 和 `configs/pdcch_result031_codebooks_2rx_ideal_formal.yaml`。两份 codebook 配置只列尚未完成的 SNR：`DFTcodebook` estimated 为 `6,7,7.5,8`，ideal 为 `5,6,7,7.5,10`；其目标网格中已由 `C300_PRG_DFT4_6REG` 完成的 estimated `5,8.5,9,10` 和 ideal `8,8.5,9` 直接复用，禁止重跑。`LTEcodebook` 两种 CSI 均运行完整 Rx2 网格。所有配置的 seed 固定为 `20260915`，各 shard 使用互异 namespace。输出写入 `outputs/experiment031_pdcch_cdd/20260915_c300_al1_sidon_transparent_cdd130_ideal/` 下按 candidate/CSI/Rx 分开，不得覆盖昨日结果。中断后沿相同绝对 trial 区间恢复。

本补充不修改链路主算法；原三份配置已直接复用此前通过配置校验和正式运行的 `physical_fullband`、ideal-CSI、duplicate-delay 与多 Rx 路径。新增 codebook shard 在正式运行前必须 validate：核对上述四个 LTE 向量及选择索引、所选 LTE 向量精确等于 `[1,-1,1,-1]/2`、两类向量均为单位范数、AL1 只有一个 bundle、estimated/ideal 接收路径和多 Rx 功率口径；若公共 codebook 参数化存在缺口，只做最小扩展并补定向测试。正式 runner 仍须在 trial 开始前完成内置配置校验并写出 resolved config；校验失败即不得产生或合并 trial。

数据完成后扩展 `tools/analyze_plan031_cdd911.py`：把 Sidon transparent 1Rx 及 `DFTcodebook`/`LTEcodebook` 2Rx 加入 estimated-CSI 对比，把 CDD130 ideal 1Rx/2Rx 及两类 codebook 2Rx 加入 ideal-CSI 对比；保留既有曲线和统一样式，图例明确 CSI、Rx 数及 codebook 名称。目标 SNR 仍要求相邻 raw bracket 不宽于 0.25 dB 且至少 95% 的 4,000 次 bootstrap replicate 保持 bracket；本轮给定粗网格若不满足，只报告 raw bracket/单侧界限，不强行插值。最终同步更新两版 `result-031`、分析表、来源回执和约 13 cm 预览，并报告两种固定预编码波形的 BLER 与 estimated-CSI CE NMSE 差异；研究者确认前不更新 `KNOWLEDGE.md`/`GOALS.md`，不创建 Git checkpoint。

执行记录：原三条曲线的旧网格 29/29 点已完成，共 770,800 trials、61,660 errors，逐 trial、Wilson 区间和 CE NMSE 审计无差异；CDD130 ideal 的 1Rx/2Rx 均取得原始 BLER 不高于 0.01 的点，Sidon transparent 1Rx 在原最高 13 dB 的最低 BLER 为 0.7514。两份 codebook 配置的 25/25 个缺失点也已完成，共 504,200 trials、6,648 errors；结合历史 `C300_PRG_DFT4_6REG` 的 7 个复用点，两类 codebook 的四条 2Rx 曲线均形成完整 8 点网格并已更新第 11.2 节图表。按本次统一 54 点网格现已有 49 点，尚需补 Sidon transparent 1Rx 的 `15,16,17,18` 和 CDD130 ideal 1Rx 的 `15`；完成这 5 点并重新审计前，本补充仍标记为部分已执行、待确认。

## 2026-09-20 补充规划：CDD911 non-transparent 与 CDD130 transparent 的 4Tx/2Rx AL2/AL4 estimated-CSI BLER（已执行，数据待分析）

现有数据只覆盖这两种口径的 AL1，尚无 4Tx/2Rx、AL2/AL4 的 estimated-CSI 曲线。本补充一次性增加 `CDD911`（`[0,0,911,911] ns`，`matched_effective`，即 non-transparent）和 `CDD130_transparent`（`[0,0,130,130] ns`，`physical_fullband`）在 AL2、AL4 下的四条曲线，不重复规划相同场景，也不新增 ideal-CSI 数据。

除 AL 及其派生的 candidate 带宽、$K$、coded bits 和 phase coordinates 外，统一沿用本 plan 的 C300 4Tx/2Rx、2-symbol estimated-CSI 正式场景：48-RB CORESET、30 kHz、FFT 4096、CP 288、TDL-C 300 ns、3 km/h、4 GHz、20 sinusoids、A=41、CRC24C/RNTI `0xFFFF`、QPSK、Polar list size 8、单位总发射功率及独立 Rx 分支 MRC。两条 candidate 均设置 `allow_duplicate_delays: true`；AL2 的 $K=72$，CDD911/CDD130 coordinates 分别为 `[0,0,1.96776,1.96776]`、`[0,0,0.2808,0.2808]`；AL4 的 $K=144$，分别为 `[0,0,3.93552,3.93552]`、`[0,0,0.5616,0.5616]`。

每个 AL 合并为一份 formal 配置。SNR 网格、`batch_size=50`、每点至少 10,000 trials、达到 200 errors 后可停止、最多 50,000 trials，以及逐 trial error flags 和线性 CE NMSE 保存要求，分别直接继承对应 AL 的既有 2Rx estimated-CSI formal 配置中相同接收口径的 small-CDD candidate；不做新的后验选点。配置为 `configs/pdcch_result031_cdd911_cdd130_2rx_estimated_al2_formal.yaml` 和 `configs/pdcch_result031_cdd911_cdd130_2rx_estimated_al4_formal.yaml`，输出统一写入 `outputs/experiment031_pdcch_cdd/20260920_c300_4tx_2rx_al2_al4_cdd911_cdd130/estimated/`，seed 为 `20260920`，两个 AL 使用互异 namespace。

依次执行 runner 的 `--stage validate` 和 `--stage run`；validate 必须核对 AL、Rx 数、coded bits、物理时延回算、重复时延保留、接收协方差口径、单位功率和输出路径。运行可中断恢复，必须从相同 namespace 的绝对 trial 区间继续。完成条件是 63 个 candidate/SNR 点均达到停止条件或 50,000 上限，并保存 resolved config、run metadata、`bler_points.csv`、逐 trial flags 和 estimated-CSI NMSE；本补充当前只保存原始数据，不生成图、不分析目标 SNR，也不修改两版 result、`KNOWLEDGE.md` 或 `GOALS.md`。

复现命令：

```powershell
python tools/run_plan031_candidates.py --config configs/pdcch_result031_cdd911_cdd130_2rx_estimated_al2_formal.yaml --config configs/pdcch_result031_cdd911_cdd130_2rx_estimated_al4_formal.yaml --stage validate
python tools/run_plan031_candidates.py --config configs/pdcch_result031_cdd911_cdd130_2rx_estimated_al2_formal.yaml --config configs/pdcch_result031_cdd911_cdd130_2rx_estimated_al4_formal.yaml --stage run --max-workers 4
```

执行记录：四个 candidate/AL shard 的 validate 均通过，63/63 个固定 candidate/SNR 点全部完成停止条件，共 1,370,950 trials、29,167 errors；四个 shard 的 resolved config、run metadata、运行日志、逐点 error flags 和 estimated-CSI NMSE 文件齐全。原始数据已保存到上述输出目录。研究者随后要求把 AL2/AL4 的 CDD911、CDD130、B0 QC、S0 Sidon 和 PRG DFT4 estimated-CSI BLER 画为一张双子图；已由 `tools/plot_plan031_cdd911_cdd130_al2_al4.py` 读取正式 CSV 生成，图中不显示误差棒，竖向主/次网格间隔分别为 0.5/0.25 dB。尚未分析目标 SNR 或修改两版 result。
