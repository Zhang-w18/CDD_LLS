# plan-033：4Rx 下 4/8Tx、3/60 km/h 的 CDD、透明 PRG 与过时 CSI MRT 对比

> 状态：执行中，正式扫描已获准开始。研究者于 2026-09-15 确认执行主机没有 GPU；033 专用 runner、配置和定向测试已实现，定向测试已通过。两个 smoke 场景均通过审计，四个初始 prescan 均完成且 48 个“场景内方案×目标 BLER”都有真实采样双侧 bracket，无需范围扩展。研究者已于 2026-09-15 确认第 10 节冻结项。`A100_NT8_NR4_V3`、`A100_NT8_NR4_V60` 与 `A100_NT4_NR4_V60` 已完成每点 10,000 trials 的初始正式预算并生成部分 result；`A100_NT8_NR4_V60` 的 5 个待追加 SNR 已完成第一轮 10,001--11,000 trials，8 个冻结端点仍待继续。其余 1%端点追加以及 `A100_NT4_NR4_V3` 正式场景仍未完成。

## 1. 目的与研究问题

本轮把 result-032 的移动 TDL-A、二维时频 LMMSE 和 5 ms aged-CSI MRT，与 result-028 已实现的 4Rx MRC 及小时延透明/非透明 CDD 放入同一链路，回答：

1. 在 8Tx/4Rx 与 4Tx/4Rx 下，`B0_QC`、严格 `SIDON`、transparent PRG6、aged-CSI MRT PRG6、small-delay transparent CDD 和 small-delay non-transparent CDD 的 10%/1% estimated-CSI BLER 如何；
2. 同一 Tx/Rx 条件从 3 km/h 增至 60 km/h 后，各方案的门限变化如何；
3. 同一速度从 8Tx/4Rx 改为 4Tx/4Rx 后，CDD 几何、透明基线和 aged MRT 的相对排序如何；
4. 小时延 CDD 中，UE 知道与不知道 CDD delay 所产生的 matched/transparent 接收机差异是否仍存在。

本轮主体只做单层 estimated-CSI 链路和 data-RE CE NMSE；第 11 节例外增加一个不进入正式门限验收的 8Tx/4Rx ideal-CSI 轻量诊断图。本轮不搜索新的 outage 最优候选，也不把结果推广到相关天线、CDL、ICI、实际码本 PMI 或多层检测。

本轮对应 `GOALS.md` 的多 Rx、移动性和发射天线数适用范围检查；候选定义遵循 `DESIGN.md` §3.4、§3.7--3.9 与 `docs/design/CDD_DELAY_SELECTION_RULES.md`。最终判据是实际 LDPC BLER，delay 几何、pilot rank、condition number 和 CE NMSE只作候选合法性与机理诊断。

## 2. 冻结系统条件

| 参数 | 冻结取值 |
|---|---|
| 波形 | 48 PRB，$K=576$ 个连续有效子载波，30 kHz SCS，FFT 4096，CP 288 samples，10 个 PDSCH symbols |
| 天线与层数 | `8Tx/4Rx/1 layer` 与 `4Tx/4Rx/1 layer` |
| 信道 | Sionna 1.0.2 TDL-A，RMS delay spread 100 ns，3.5 GHz，20 sinusoids；Tx--Rx 分支独立、同 PDP、无空间相关 |
| 速度 | 3 km/h 与 60 km/h；`min_speed=max_speed`，每个场景速度固定 |
| DMRS | symbols `[2,7]`，comb-6；每 symbol 96 pilot RE，共 192 个独立 LS 观测；5568 data RE |
| 接收机 | 两个 DMRS 不平均；用场景真实速度的时间协方差和指定频域协方差做二维时频 LMMSE；各 Rx 独立估计后在 data RE 做 MRC |
| 调制编码 | 16QAM；NR 256QAM MCS table 的 MCS 8；目标码率 553/1024；Sionna LDPC 最多 8 次迭代 |
| CSI age | 10 slots，名义 5 ms；旧 CSI 为第 0 个 symbol sample，当前 PDSCH 为第 140--149 个 sample，实际 age 4.994791667 ms |
| 随机性 | base seed `20260727`；absolute-trial 可重放；同一场景、SNR、trial 的六条曲线共享当前 TDL、payload、编码比特、data AWGN 和 DMRS AWGN |

这里的 TDL-A 是 3GPP TDL-A 功率时延轮廓；3/60 km/h 两个场景都具有随 OFDM symbol 演化的时变衰落，不能称为零速度 block-static 信道。只采用 OFDM-symbol 采样的时变频域乘法模型，假设定时和载波同步，不模拟 CFO、ICI 或 ISI；多普勒只表现为 symbol 间信道变化和 CSI 老化。

`num_sinusoids=20` 是 Sionna 的 sum-of-sinusoids 时变衰落近似阶数，用于合成每个 TDL diffuse tap 的 Doppler 时间过程；它不是 TDL tap 数、传播路径数、天线数或 trial 数。该值保持 result-032 的取值，不在本轮调参。

## 3. 归一化、SNR 与多 Rx 定义

遵循 `DESIGN.md` 的新仿真入口规范，发射端显式使用单位总功率预编码。若 $\mathbf V$ 的元素模长为 1，则实际预编码向量为 $\mathbf W=\mathbf V/\sqrt{N_t}$，每个频点满足 $\lVert\mathbf W_k\rVert_2^2=1$，并令

$$
\sigma_n^2=\frac{1}{10^{\mathrm{SNR}_{\rm dB}/10}}.
$$

SNR 定义为每根 Rx 分支的平均接收 SNR；不得因 $N_r=4$ 再乘或除 4。因此 4Rx MRC 包含约 6.02 dB 的平均合并功率增益和接收分集，但横轴不是四分支合计 SNR。历史入口中“$\lVert\mathbf V_k\rVert_2^2=N_t$ 且噪声方差为 $N_t/\mathrm{SNR}$”的内部缩放在完整链路上一致应用时与本口径数学等价，但 033 不采用该内部表示，功率、协方差和噪声审计均直接记录单位总发射功率口径。

033 的 CDD、DFT-PRG 和 aged-MRT 三条发射路径都必须显式调用单位范数模式，不允许先生成平方范数为 $N_t$ 的向量后只在噪声或接收机一侧补偿。由此冻结以下缩放：

- matched CDD 使用 $\mathbf R_g=\mathbf R_{\rm phy}\odot(\mathbf W\mathbf W^H)$，其对角功率为 1；
- transparent PRG、transparent CDD 和 aged MRT 使用 $\mathbf R_{\rm phy}$，不再使用 $N_t\mathbf R_{\rm phy}$；
- data AWGN 和每个独立 DMRS LS 观测的噪声方差均为 $1/\mathrm{SNR}_{\rm linear}$；
- aged-MRT、DFT 和 CDD 的每个频点/PRG 权值均满足平方范数 1。

归一化不改变 delay residue、Sidon 二元和、pilot rank 或 condition number。若未归一化的 pilot 矩阵满足 $\mathbf V_P^H\mathbf V_P=N_p\mathbf I$，则本轮实际矩阵满足 $\mathbf W_P^H\mathbf W_P=(N_p/N_t)\mathbf I$。它会改变协方差和绝对幅度的书写方式，但在信号、协方差和噪声三处一致缩放时不应单独造成 BLER 位移。033 禁止混用两种表示，也不得直接合并 result-028/032 的历史 trial；历史门限只用于预扫描范围参考。

estimated-CSI MRC 对每个 data RE 使用

$$
\widehat x=
\frac{\sum_{r=0}^{3}\widehat g_r^*y_r}
{\sum_{r=0}^{3}|\widehat g_r|^2},
\qquad
\widehat\sigma_{\rm eff}^2=
\frac{\sigma_n^2}{\sum_{r=0}^{3}|\widehat g_r|^2}.
$$

LLR 噪声方差不额外加入 CE-error-aware 项，保持现平台口径。CE NMSE 对每个 trial 先在四根 Rx 和全部 5568 个 data RE 上求误差能量/真实能量比，再在线性域跨 trial 平均并转 dB；保存和、平方和及 95% Monte Carlo 区间。

## 4. 4Tx/8Tx 时延推导与完整冻结表

当前参数给出

$$
N_p=K/6=96,
\qquad
\frac{1}{K\Delta f}=57.870370\ {\rm ns},
\qquad
\tau_{\rm alias}=\frac{1}{6\Delta f}=5.555556\ \mu{\rm s}.
$$

表中所有 CDD 物理时延统一按 $\tau_n=j_n/(K\Delta f)$ 换算；这里 $K=576$、$\Delta f=30$ kHz，因此一个 delay-grid 单位等于 57.870370 ns。

4Tx 候选不能直接携带或截取 8 维 delay。时延只按解析几何和基本导频可辨识性确定，不使用 CE、outage 或 BLER 仿真筛选。`B0_QC` 取 comb-6 折叠圆上等间隔且最大分离的四个 residue，即公差

$$
\Delta j_{\rm QC}=\frac{N_p}{N_t}=\frac{96}{4}=24.
$$

严格 Sidon 则使用 8Tx 历史集合的四项嵌套构造。下列两张表是本轮六方案的完整发射端冻结表；`不适用（无人工 CDD delay）` 是明确配置值，不是待补字段。

| 8Tx 方案/curve ID | delay grid $\mathbf j$ | 人工 delay (ns) | 其余发射端冻结值 |
|---|---|---|---|
| `A100_NT8_B0_QC` | `[0,9,18,27,36,45,54,63]` | `[0,520.833333,1041.666667,1562.500000,2083.333333,2604.166667,3125.000000,3645.833333]` | matched CDD；UE 知道 delay |
| `A100_NT8_S0_SIDON` | `[0,1,3,7,12,20,30,65]` | `[0,57.870370,173.611111,405.092593,694.444444,1157.407407,1736.111111,3761.574074]` | matched CDD；UE 知道 delay |
| `A100_NT8_TRANSPARENT_PRG6` | 不适用（无人工 CDD delay） | 不适用（无人工 CDD delay） | 8 个单位范数 DFT8 向量；PRG index `[0,1,2,3,4,5,6,7]`；UE transparent |
| `A100_NT8_AGED_MRT_PRG6` | 不适用（无人工 CDD delay） | 不适用（无人工 CDD delay） | 每 trial、每 PRG 由旧信道计算单位范数主特征向量；CSI age 140 symbols = 4.994791667 ms；UE transparent |
| `A100_NT8_SMALL_CDD_QSTEP0P25_TRANSPARENT` | `[0,0.25,0.5,0.75,1,1.25,1.5,1.75]` | `[0,14.467593,28.935185,43.402778,57.870370,72.337963,86.805556,101.273148]` | small-delay CDD；UE 不知道 delay |
| `A100_NT8_SMALL_CDD_QSTEP0P25_MATCHED` | `[0,0.25,0.5,0.75,1,1.25,1.5,1.75]` | `[0,14.467593,28.935185,43.402778,57.870370,72.337963,86.805556,101.273148]` | 与上一行发射向量逐 RE 完全相同；仅 UE 改为知道 delay |

| 4Tx 方案/curve ID | delay grid $\mathbf j$ | 人工 delay (ns) | 其余发射端冻结值/comb-6 检查 |
|---|---|---|---|
| `A100_NT4_B0_QC` | `[0,24,48,72]` | `[0,1388.888889,2777.777778,4166.666667]` | matched CDD；residue `[0,24,48,72]` 等间隔；fold gap 24；pilot rank 4；condition number 1 |
| `A100_NT4_S0_SIDON` | `[0,1,3,7]` | `[0,57.870370,173.611111,405.092593]` | matched CDD；residue 互异；pilot rank 4；condition number 1；10 个模 576 无序二元和全部唯一 |
| `A100_NT4_TRANSPARENT_PRG6` | 不适用（无人工 CDD delay） | 不适用（无人工 CDD delay） | 4 个单位范数 DFT4 向量；PRG index `[0,1,2,3,0,1,2,3]`；UE transparent |
| `A100_NT4_AGED_MRT_PRG6` | 不适用（无人工 CDD delay） | 不适用（无人工 CDD delay） | 每 trial、每 PRG 由旧信道计算单位范数主特征向量；CSI age 140 symbols = 4.994791667 ms；UE transparent |
| `A100_NT4_SMALL_CDD_QSTEP0P25_TRANSPARENT` | `[0,0.25,0.5,0.75]` | `[0,14.467593,28.935185,43.402778]` | small-delay CDD；UE 不知道 delay；连续坐标不得舍入 |
| `A100_NT4_SMALL_CDD_QSTEP0P25_MATCHED` | `[0,0.25,0.5,0.75]` | `[0,14.467593,28.935185,43.402778]` | 与上一行发射向量逐 RE 完全相同；仅 UE 改为知道 delay；连续坐标不得舍入 |

4Tx `B0_QC` 是按当前 $N_p$ 和 $N_t$ 重新生成的等差 QC 参考，不是 8Tx `[0,9,\ldots,63]` 的截断。它与 `docs/design/CDD_DELAY_SELECTION_RULES.md` 的 $\tau_{\rm alias}/N_t$ 等差构造相同；选择依据只是最大化四个导频折叠中心的等间隔，不声称已经由链路性能筛选。等差集合不要求满足 Sidon 二元和唯一性；其 10 个无序二元和只有 7 个不同值，必须如实写入 delay audit。`S0_SIDON` 使用 8Tx 历史集合的前四项，也是按升序逐项接受的严格 Sidon 构造；其二元和为 `[0,1,2,3,4,6,7,8,10,14]`，满足 `docs/design/CDD_DELAY_SELECTION_RULES.md` §3 的严格 Sidon 条件。

TDL-A 100 ns 的 $T_{0.01}=479.66$ ns。在 $T_{\rm margin}=T_{\rm sync}=0$ 的诊断口径下，$g_\Sigma=17$、$g_{\rm fold}=9$。上述 `S0_SIDON` 的实际 pair/fold gap 均为 1，因此它是严格 Sidon，但不是硬厚 Sidon；本轮不得将其标记为 thick Sidon。4Tx 的硬厚装填必要上界为 pair/fold `57/24`，只说明硬厚候选没有被必要条件排除，不构成本轮新增候选的授权。

所有 CDD 均按

$$
V_{k,n}=\exp\left(-j2\pi k j_n/576\right)
$$

构造，以首个 active subcarrier 为相位参考。人工 delay 是数字相位/循环移位，不计入真实传播时延或 CP。manifest 必须保存 grid 坐标、ns、FFT-sample 等效值、residue、无序二元和、pair/fold gap、pilot rank、condition number 和选择规则。

## 5. 每个场景的六条曲线

四个场景各运行以下六条 estimated-CSI 曲线；同一场景使用相同公共 SNR 网格和 paired trials。

| 曲线 | 发射端 | 接收端频域协方差/知识 |
|---|---|---|
| `B0_QC` | 第 4 节对应 $N_t$ 的 B0 delay | UE 知道 delay；全带 matched $\mathbf R_g=\mathbf R_{\rm phy}\odot(\mathbf W\mathbf W^H)$ |
| `S0_SIDON` | 第 4 节对应 $N_t$ 的严格 Sidon delay | UE 知道 delay；全带 matched $\mathbf R_g$ |
| `transparent PRG6` | 48 PRB 分为 8 个 6-RB PRG | UE 不知道 DFT index；每 PRG 使用 $\mathbf R_{\rm phy}$，不跨 PRG 插值 |
| `aged-CSI MRT PRG6` | 5 ms 旧 CSI 生成的未量化 PRG 主特征向量 | UE transparent；每 PRG 使用 $\mathbf R_{\rm phy}$，不跨 PRG 插值 |
| `small-delay CDD, transparent` | 第 4 节对应 $N_t$ 的 q-step 0.25 CDD | UE 不知道 delay；全带使用 $\mathbf R_{\rm phy}$ |
| `small-delay CDD, non-transparent` | 与上一行完全相同 | UE 知道 delay；全带使用 matched $\mathbf R_g$ |

PRG DFT 映射固定为：

- 8Tx：8 个 8 维 DFT 向量按 PRG `[0,1,2,3,4,5,6,7]` 各使用一次；
- 4Tx：4 个 4 维 DFT 向量按 PRG `[0,1,2,3,0,1,2,3]` 使用两轮。

每个 DFT/MRT 向量平方范数为 1。4Rx aged MRT 对每个 6-RB PRG 构造

$$
\mathbf G_b=
\sum_{k\in b}\mathbf H_{{\rm old},k}^{H}\mathbf H_{{\rm old},k},
\qquad
\mathbf w_b=\mathbf u_{\max}(\mathbf G_b),
$$

其中 $\mathbf H_{{\rm old},k}\in\mathbb C^{4\times N_t}$，并在当前 10 个 PDSCH symbols 内保持不变。这是过时、未量化 subband CSI 的 MRT 上界型基线，不是 3GPP Type-I 码本 PMI；result 不得混称。

## 6. 四个场景与比较公平性

场景 ID 冻结为：

1. `A100_NT8_NR4_V3`；
2. `A100_NT8_NR4_V60`；
3. `A100_NT4_NR4_V3`；
4. `A100_NT4_NR4_V60`。

同一场景内六条曲线严格 paired，可用 paired bootstrap 比较目标 SNR。不同速度的信道时间协方差不同，不宣称跨速度 paired 方差缩减；不同 $N_t$ 的信道维度、总 delay span 和 DFT codebook 也不同，跨 $N_t$ 差值使用独立场景不确定性并明确不能归因于单一因素。

3 km/h 时仍使用 5 ms aged CSI，不把 aged MRT 替换为当前 CSI。两档速度均从同一 Sionna SoS 定义生成连续 realization；old/current 样本不得来自两个独立 realization。

## 7. 实现范围、测试与 smoke

现有 `tools/run_plan032_tdl_mobility.py` 把 `n_tx=8`、`n_rx=1`、速度 60 km/h、功率 8 和 DFT8 映射写死；`tools/run_bler_curves.py` 的 4Rx 支持则主要覆盖 zero-speed/static 路径。033 需要把两条已验证能力组合为一个可配置但受本 plan schema 严格限制的入口，不能只修改 YAML 假定现有 runner 已支持。

预计工作：

1. 把移动 TDL 历史/当前抽样、二维时频 LMMSE、PRG、aged MRT 和 Rx-MRC 的可复用部分放入 `cdd_lls/`，或在不复制链路主循环的前提下抽取公共 runner；保留 result-032 旧入口行为；
2. 新增 `tools/run_plan033_tdl_mobility_mimo.py`，schema 显式包含 `n_tx`、`n_rx=4`、`speed_kmh`、candidate manifest、SNR/trial 区间和输出目录；CDD、DFT-PRG、aged-MRT 均强制 `normalize=true`，展开配置不得接受隐式默认值；
3. 新增 4Tx/8Tx delay manifest 生成与审计，禁止把 8 维 delay 截断发生在仿真内核的隐式路径；
4. 新增 prescan、formal、merge/analyze 配置和可恢复调度入口；已有点追加只能使用连续且不重叠的 absolute-trial 区间；
5. 新增 `tools/analyze_result033_tdl_mobility_mimo.py`，只读原始 CSV/NPY 汇总生成数值审计、paired/independent bootstrap 和本地图。

测试遵循“只覆盖 033 增量，不复制既有测试”的原则。result-028 已覆盖 4Rx CE/MRC、transparent/matched CDD、PRG 边界和 interval 合并；result-032 已覆盖移动 TDL 历史抽样、二维时频 LMMSE、时间相关、aged MRT、seed 重放和恢复编排。本轮不为这些行为另写同义测试，只运行相关既有测试作为回归。

新增测试仅覆盖：

1. 033 配置能精确展开四个冻结的 `n_tx/n_rx/speed` 组合，单位总功率、噪声和协方差缩放在 4Tx/8Tx 下没有遗留硬编码常数；
2. 4Tx B0/Sidon/small-delay 的冻结数组及其解析审计，以及 DFT4 的 PRG 映射 `[0,1,2,3,0,1,2,3]`；
3. 4Rx aged MRT 的 $\mathbf G_b$ 跨四根 Rx 和 PRG 子载波求和，与小型手算矩阵一致，且 $\lVert\mathbf w_b\rVert_2^2=1$；
4. 一个轻量 runner 集成测试确认同一场景六条曲线共享 trial key，并把 `n_tx/n_rx/speed/candidate/receiver` 写入互不混合的输出记录。

回归只定向运行已有的 4Rx MRC、history sampling、PRG 二维滤波和 aged-MRT 测试，再对一份 result-028 4Rx 配置与一份 result-032 mobility 配置执行 `--stage validate`；不运行整个历史测试集。

smoke 只跑两个覆盖组合边界的场景：`A100_NT4_NR4_V3` 与 `A100_NT8_NR4_V60`。每场景六条曲线共享 20 trials，并取一个低 SNR 和一个高 SNR；其余两个场景由随后的正式 prescan 首点承担集成检查，不再重复 smoke。执行主机为 CPU-only；smoke 检查六曲线落盘、数组有限、功率/维度、old/current 同 realization、CPU placement、逐 batch 采样的进程 RSS 内存峰值和墙钟耗时，不重复做已由单元测试覆盖的统计精度或中断恢复压力测试。batch size 先从 result-028 4Rx 已验证的 25 起步；若移动历史张量导致 OOM，只允许降低 batch，不改变 trial、seed 或物理定义，并把选择写入展开配置。先完成 `A100_NT4_NR4_V3` smoke 并依据实际耗时和 RSS 决定是否继续第二个 smoke 与 prescan，不根据历史结果猜测 CPU 预算。

实际 smoke 审计：两个场景均得到 6 条曲线、2 个 SNR、24 个有限逐 trial 数组和2份功率/MRC 诊断。`A100_NT8_NR4_V60` 的两个点墙钟耗时分别为 1.999 s 和 1.843 s，逐 batch 采样的进程 RSS 峰值为 4,136,407,040 bytes，执行设备为 CPU；MRT 权值平方范数范围为 `[0.9999999999999984,1.0000000000000020]`。`A100_NT4_NR4_V3` 的两个点墙钟耗时分别为 2.423 s 和 1.869 s；其 trial 完成于 CPU runtime 诊断字段加入前，因此峰值 RSS 不可追溯，兼容审计明确标记该缺失，不重跑已完成 trial。证据分别位于 `outputs/experiment033_tdl_mobility_mimo/20260915_main/smoke/A100_NT4_NR4_V3/` 和 `outputs/experiment033_tdl_mobility_mimo/20260915_main/smoke/A100_NT8_NR4_V60/`。

## 8. 预扫描、正式 SNR 冻结与统计预算

预扫描不并入正式结果。每点 400 个共同 paired trials，初始网格为：

- 3 km/h 两个场景：`[-2,0,2,4,6,8,10,12] dB`；
- 60 km/h 两个场景：`[2,4,6,8,10,12,14,16] dB`。

选择依据是 result-028 的 4Rx 接收增益和 result-032 的移动门限；这些链路预扫描只用于定位正式 SNR 网格，不参与 4Tx delay 选择，也不是正式结果。若任一曲线尚未同时观察到 10%上下侧或 1%上下侧点，只向缺失方向按 2 dB 扩展；每个场景最多扩展到 `[-6,20] dB`。达到边界仍未 bracket 时暂停并报告，不凭外推冻结正式 crossing。

`A100_NT8_NR4_V3` 初始 prescan 在旧 curve ID 下完成后，使用 `tools/migrate_plan033_curve_ids.py` 将 aged MRT 和 small-delay matched 的 ID 迁移到第 4 节冻结名称。迁移覆盖本轮两个 smoke 和该 prescan：48 个 NPY 文件只重命名且迁移前后 SHA-256 不变，3 个 interval 文件的 144 条数组引用均存在，72 个 interval record 的 error count、CE sum 和 CE sumsq 与逐 trial 数组一致，candidate manifest 与 resolved-run hash 已重算。完整收据为 `outputs/experiment033_tdl_mobility_mimo/20260915_main/curve_id_migration_receipt.json`；该迁移不改变任何信道、payload、噪声、error flag 或 CE 数值。

预扫描后、formal trial 1 前，为每个场景冻结一套六曲线共同正式网格：10%和 1% crossing 邻域最大间隔 0.25 dB，过渡区最大 0.5 dB；只保留形成双侧 bracket 所需范围及两侧最多一个保护点。不得为描绘低于 1% 的尾部继续加高 SNR；若公共网格因最困难曲线延伸，使其他曲线远低于 1%，这些共享点只按公共最低预算执行，不为其单独追误块。

正式预算：

1. 每个冻结点先运行至少 10,000 个共同 trials，单次可恢复 interval 不超过 1,000；
2. 每条曲线实际 1% bracket 的两个端点若 BLER 位于 0.5%--2%，追加至至少 200 个错误块或 50,000 trials 上限；低于 0.5%的保护/共享点不追 200 errors；
3. 10% bracket 端点在 10,000 trials 下自然具有充足错误数；仍报告 Wilson 95%区间；
4. 达到 50,000 trials 仍不足时保留原始点并标注样本上限，不平滑、不替换为伪计数；
5. 每条曲线一旦具有真实采样点形成的 10%和 1%双侧 bracket，即完成该曲线目标；全场景六条均闭合后停止向高 SNR 扩展。

端点追加由 `tools/run_plan033_append.py` 按 `append_requirements_initial.csv` 中冻结的 1% bracket 端点调度。每次调用对仍未达到 200 errors 且未到 50,000 trials 的 SNR 只推进一个 1,000-trial paired interval；首次追加前将原始配置、展开配置、hash、manifest、审计和代码状态复制到场景目录的 `append_scheduler/baseline/`，每轮请求保存在 `append_scheduler/requests/`，从而避免追加配置覆盖初始 10,000-trial 证据。重复同一命令根据最新合并 CSV 重新判停，不重跑已有 absolute trial。

逐点 BLER 报告 Wilson 95%区间。10%/1%目标 SNR 仅在相邻真实采样点的双侧 bracket 内按 log-BLER 线性插值；同一场景的方案差值用 absolute-trial paired bootstrap。图上的 marker 和线段使用原始 Monte Carlo 点，不用 PAVA、PCHIP、平滑、单调修正或外推。

## 9. 输出、分析与 result 要求

正式产物写入

`outputs/experiment033_tdl_mobility_mimo/<run_id>/<stage>/<scenario_id>/`，

至少保存：原 YAML、展开配置及 SHA-256、代码版本/工作区变更标识、candidate/receiver manifest 及 hash、delay audit、逐点 BLER/CE CSV、error flags、CE trial arrays、absolute-trial intervals、seed/pairing audit、滤波器诊断、MRC/功率诊断、old/current 时间相关、运行日志、耗时、CPU placement 和逐 batch 采样的进程 RSS 内存峰值。

分析必须逐场景报告六条曲线的：

- 10%/1% crossing、bracket、错误数、trials、Wilson 区间和 crossing bootstrap 区间；
- 相对 `transparent PRG6` 的 SNR 差；
- `S0_SIDON - B0_QC`、small-delay transparent/non-transparent、aged MRT/transparent PRG6 的成对差值；
- data-RE CE NMSE 及 95%区间；
- 3→60 km/h 和 8→4Tx 的目标变化，并明确跨场景不是 paired 比较。

本地脚本生成四张 estimated-CSI BLER 主图（每场景一张）和必要的 CE 图，纵轴重点覆盖 $10^{-1}$ 与 $10^{-2}$，不为显示更低 BLER 扩展仿真。另生成速度对比和 Tx 数对比的目标 SNR 汇总表；只有表格过密时才生成汇总图。所有图从保存的本地 CSV 生成，精确图数据另存 CSV，样式、尺寸和 13 cm 预览遵循 `docs/agent/RESULT_SPEC.md`；Agent 不加载结果图片，研究者人工核对预览。

正式完成后成对生成：

- `research/result-033-PDSCH-4Rx-MIMO移动性.md`；
- `research/result-033-PDSCH-4Rx-MIMO移动性-text.md`。

两版必须给出相同的配置、关键原始数值、异常、判定和证据路径；无图版不得嵌图。结果未经研究者确认前不更新 `KNOWLEDGE.md`/`GOALS.md`，不创建 Git checkpoint。

## 10. 正式 SNR 网格冻结区

四个场景冻结如下；数组单位均为 dB，同一场景六条曲线共享全部点和 absolute-trial key：

| 场景 | 正式公共 SNR 数组 | 点数 |
|---|---|---:|
| `A100_NT8_NR4_V3` | `[-2,-1.75,-1.5,-1.25,-1,-0.75,-0.5,-0.25,0,0.5,1,1.5,2,2.5,3,3.5,4,4.25,4.5,4.75,5,5.25,5.5,5.75,6,6.25,6.5,6.75,7,7.25,7.5,7.75,8]` | 33 |
| `A100_NT8_NR4_V60` | `[4,4.25,4.5,4.75,5,5.25,5.5,5.75,6,6.25,6.5,6.75,7,7.25,7.5,7.75,8]` | 17 |
| `A100_NT4_NR4_V3` | `[0,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.5,3,3.5,4,4.25,4.5,4.75,5,5.25,5.5,5.75,6,6.25,6.5,6.75,7,7.25,7.5,7.75,8]` | 29 |
| `A100_NT4_NR4_V60` | `[4,4.25,4.5,4.75,5,5.25,5.5,5.75,6,6.25,6.5,6.75,7,7.25,7.5,7.75,8]` | 17 |

取点依据为 `outputs/experiment033_tdl_mobility_mimo/20260915_main/prescan/prescan_bracket_audit.json`（SHA-256 `5ef1c3dba65a33e9ae1f5ef1d45e558bba553b3bfb555776c7b12fa3564651a7`）。四份 prescan 原始汇总依次为：

- `outputs/experiment033_tdl_mobility_mimo/20260915_main/prescan/A100_NT8_NR4_V3/final/estimated_csi_bler_points.csv`；
- `outputs/experiment033_tdl_mobility_mimo/20260915_main/prescan/A100_NT8_NR4_V60/final/estimated_csi_bler_points.csv`；
- `outputs/experiment033_tdl_mobility_mimo/20260915_main/prescan/A100_NT4_NR4_V3/final/estimated_csi_bler_points.csv`；
- `outputs/experiment033_tdl_mobility_mimo/20260915_main/prescan/A100_NT4_NR4_V60/final/estimated_csi_bler_points.csv`。

每份 CSV 均为 6 条曲线 × 8 个 SNR × 400 paired trials；共 384 份 error/CE trial 数组已复算，error count、CE sum 和 CE sumsq 全部与 interval CSV 一致。四场景的 10%/1% 共 48 个目标均为 `bracketed`，没有执行范围扩展。正式网格在每个已观察 crossing bracket 内取 0.25 dB 间隔；`A100_NT8_NR4_V3` 和 `A100_NT4_NR4_V3` 的两个分离 crossing 区之间取 0.5 dB 过渡点，不增加低于 1% 的高 SNR 尾部或额外保护点。

formal YAML 与文件 SHA-256 冻结为：

| YAML | SHA-256 |
|---|---|
| `configs/plan033_nt8_nr4_v3_formal.yaml` | `88ed2830ec2d9b3976f3183d7b34400922614c14f818332351ba9c9b5cde635f` |
| `configs/plan033_nt8_nr4_v60_formal.yaml` | `a59b244f08b2b1e8ec6b2beb1aadb8ea56a0d2fc013f6dfe0dd8ec12827b37a4` |
| `configs/plan033_nt4_nr4_v3_formal.yaml` | `eb595f4033b51a73e4b895e0cfaad1c0cfaf44a8be6a0faa2c4415193bd00983` |
| `configs/plan033_nt4_nr4_v60_formal.yaml` | `e4098cec779232c8639b0e62d2ab046439832e20de5bba50b26a599d5593bf00` |

`batch_size=25`，每点初始目标为 10,000 paired trials。formal runner 每次调用对每个尚未完成的 SNR 只推进一个不超过 1,000 trials 的连续 absolute-trial interval并立即落盘；因此初始阶段重复调用 10 次即可完成，意外中断后重发同一命令会从已落盘边界恢复。初始公共预算为 96 点 × 10,000 = 960,000 paired trials（5,760,000 curve-trials）。追加只针对实际 1% bracket 中 BLER 位于 0.5%--2% 的端点，逐 SNR 把六曲线共同目标提高至获得至少 200 errors 或 50,000 trials；按每场景六曲线最多 12 个互异端点计算，严格总上界为 2,880,000 paired trials（17,280,000 curve-trials）。

研究者确认日期：**2026-09-15**。研究者已明确确认本节冻结的正式公共 SNR 网格、配置 hash 和预算，允许开始 formal trial 1；后续不得在不留下修订记录的情况下改变这些冻结项。

## 11. 8Tx/4Rx ideal-CSI Sidon 与 PRG6 轻量诊断图

本节只回答一个诊断问题：在 `A100_NT8_NR4_V3` 的单层 8Tx/4Rx 链路上，`S0_SIDON` 与 `transparent PRG6` 的 ideal-CSI BLER 原始曲线如何。选择 3 km/h 是为了在一个 PDSCH 内尽量减少时间分集对观察的干扰；本节不比较速度，也不把这张图解释为 60 km/h 结论。

发射定义、单位总发射功率、每根 Rx 分支 SNR、TDL-A 100 ns、DMRS/data RE、MCS 和 LDPC 与本 plan 前文相同。`S0_SIDON` 使用 `j=[0,1,3,7,12,20,30,65]`；`transparent PRG6` 的 8 个 6-RB PRG 依次使用 DFT8 向量 `[0,1,2,3,4,5,6,7]`。ideal-CSI MRC 在每个 data RE 直接使用真实等效信道；不运行 CE，不使用 DMRS 估计值，但仍实际生成 payload、AWGN、软解调和 LDPC 译码。

固定 seed `20260727`、batch size `25`、SNR 网格 `[3.5,3.75,4,4.25,4.5,4.75,5,5.25,5.5,5.75,6] dB`，每点 `1000` 个 paired trials。同一 `SNR + absolute trial` 下两条曲线共享底层时变信道、payload 和原始 data noise。这个预算只用于快速观察曲线间隔：CSV 保留 trials、error count、BLER 和 Wilson 95% 区间，零误块点保留真实 `BLER=0`，仅画图时使用 `0.5/trials` 下界。不拟合 10%/1% crossing，不进行 bootstrap，不追加到正式错误数，不得将点估计的微小差异写成已验证全局结论。

新增 `tools/plot_plan033_ideal_sidon_prg6.py`，单次命令完成链路运行、可恢复 interval 落盘、原始 CSV 合并和本地 PNG 生成。产物写入 `outputs/experiment033_tdl_mobility_mimo/20260916_ideal_sidon_prg6/A100_NT8_NR4_V3/`，至少包含展开配置、interval CSV、逐 trial error flags、`ideal_csi_bler_points.csv`、`plot_data.csv`、样式 JSON、常规 PNG 和 13 cm 预览 PNG。本节不要求 smoke 或新增测试；执行前只需做 Python 语法检查，实际图由研究者运行后人工核对。
