# plan-023：Track B 先导——滑窗透明接收机下「局部聚簇 + 全局轮转」V 的纯数值 Pareto 扫描

> 本计划自包含：实现 agent 只读本文件即可开工，所有符号在本文件内定义。
> 本轮**不跑 estimated-CSI BLER 链路仿真**（只有 §6 校准需要 2 条 ideal-CSI BLER 曲线），主体是确定性闭式计算 + Monte Carlo outage，与 plan-022 互相独立，可并行执行。

---

## 1. 仿真目的

研究目标（可追溯到 `GOALS.md` 北极星问题与 `DESIGN.md` §5.1 研究者判据）：寻找可参数化的频域相位预编码矩阵 `V`，满足 ① 同 DMRS 开销下明显优于**透明基线**（透明小 delay CDD、PRG precoder cycling）；② 性能超过或接近 **QC 非透明方案**（显式 CDD + V-aware matched LMMSE）。

本轮验证三个假设（对应 `DESIGN.md` §7 的可证伪预测 P5 / P4 / P6+P7）：

- **H1（主假设，P5，Track B 存亡判定）**：在固定的 V-agnostic **4RB 滑窗 MMSE 参考接收机**下，「局部聚簇 + 全局轮转」类透明设计（§4 的 CC / CN 家族）在（闭式失配 NMSE，MC outage SNR）二维平面上**严格支配**透明 CDD 与 PRG cycling 构成的基线 Pareto 前沿。
- **H2（P4）**：在**宽带鲁棒 Wiener** 参考接收机下，透明设计空间塌缩——CE 代价只由等效信道复合时延跨度决定，与相位结构无关（同跨度候选 NMSE 打平）。
- **H3（附件，P6+P7）**：平坦模型 MC outage 下，栅格等差 CDD 家族各步长的 outage 曲线彼此重合（同轨道论证）；Sidon 型 delay 集的 1% 尾部不劣于任何等差集。

同时给出判据②的数值级预判：最优透明候选与 QC 非透明基准点在两轴上的差距。

---

## 2. 系统模型与符号（自包含）

- 发射分支数 `N_t = 8`，单接收天线，单层。
- 子载波间隔 `Δf = 30 kHz`；有效子载波数 `K`：48 PRB → `K=576`，24 PRB → `K=288`（本轮只跑这两个带宽）。
- 频域相位矩阵 `V ∈ C^{K×N_t}`，恒模：`V[k,n] = exp(j·φ[k,n])`，`k=0..K-1` 子载波索引，`n=0..N_t-1` 分支索引。
- 底层分支信道 `h ~ CN(0, I_8)`（每分支一个复标量；物理 delay spread 5 ns ≪ 1/(K·Δf)，带内近平坦，故用平坦模型）。
- 等效频域信道 `g = V h ∈ C^K`；其协方差 `R_g = V V^H`（`R_h = I` 时）。
- **DMRS**：comb 间隔 `S_f = 24` 子载波，导频位置 `k = 24·p`，`p = 0..N_p-1`，`N_p = K/24`（48PRB→24，24PRB→12）；2 个 DMRS 符号时域平均，静态信道下等效为单组频域观测、导频 LS 噪声方差 `σ_LS² = σ²/2`。所有候选、所有基线共用同一 DMRS 配置（同开销比较的公平性由此保证）。
- **SNR 归一化**：定义平均每 RE 接收 SNR `snr̄`，即 `σ² = N_t / snr̄`（因 `E|g_k|² = N_t`）。全文 NMSE 和 outage 都以 `snr̄`（dB）为横轴。注意此归一化与链路仿真的 Es/N0 定义可能差一个常数，只做候选间相对比较。
- DFT 栅格 delay：`τ = j/(K·Δf)`，`j` 整数（栅格步长 48PRB 57.87 ns，24PRB 115.74 ns）；栅格 delay `j` 的相位为 `φ[k,n] = -2π·k·j_n/K`。
- 导频无混叠时延周期 `τ_alias = 1/(S_f·Δf) = 1.389 µs`。
- 谱效率 `R = 4 × 553/1024 = 2.1602 bit/RE`（16QAM，code rate 553/1024）。

---

## 3. 两轴评价指标 + 基准点（全部确定性或低成本 MC，无 BLER 链路）

### 3.1 分集轴：MC outage（每候选亚秒~秒级）

1. **预计算 `I_QAM` 查找表（只算一次）**：单位能量 16QAM 在 AWGN 下的 BICM 互信息 `I_QAM(ρ)`（bit/RE），`ρ` 为线性接收 SNR。在 `ρ_dB ∈ [-20, 35]`、步长 0.1 dB 的网格上用 MC 积分（每点 ≥ 4×10^5 噪声样本）或 Gauss–Hermite 求积计算；查表时对 dB 轴线性插值，越界处 clamp 到 [0, 4]。
2. **抽样**：`h ~ CN(0, I_8)` 共 `N_mc = 2×10^5` 个样本，随机种子固定 `20260709`；**同一带宽内所有候选共用同一组 h 样本**（common random numbers，候选间成对比较）。
3. 对每个样本算块平均互信息 `I(h) = (1/K)·Σ_k I_QAM(snr̄·|g_k|²/N_t)`。
4. `P_out(snr̄) = Pr[I(h) < R]`，`snr̄ ∈ [0, 24]` dB 步长 0.5 dB；对数域插值报告 **10% 与 1% outage SNR**。

### 3.2 CE 轴：固定参考接收机的闭式失配 NMSE（确定性，零 trial）

参考接收机是**问题定义的一部分**，估计器矩阵由"假设先验"定死、不依赖真实 `V`；真实性能用精确闭式失配 MSE 计算。

**通用公式**：给定窗口内导频集合 `P`、目标数据 RE 集合 `D`（数据 RE = 全部非导频子载波），接收机用**假设协方差** `R̃` 构造估计器
```
A = R̃_DP · (R̃_PP + σ_LS²·I)^{-1}
```
真实 NMSE 用**真实协方差** `R_g = V V^H` 的子块精确算出：
```
MSE(D) = tr( A(R_PP + σ_LS²I)A^H − A·R_PD − R_DP·A^H + R_DD )
NMSE   = Σ_windows MSE(D) / Σ_windows tr(R_DD)
```
其中 `R_PP, R_PD, R_DP, R_DD` 是 `R_g` 在导频/数据行列上的子块。NMSE 在 `snr̄ ∈ {8, 16, 24}` dB 三点计算。

**假设协方差（鲁棒 Wiener 先验）**：均匀 PDP `[0, τ_w]`，
```
R̃[k,l] = N_t · exp(−jπ(k−l)Δf·τ_w) · sinc((k−l)Δf·τ_w),   sinc(x)=sin(πx)/(πx)
```
（接收机被授予真实平均功率 `N_t` 与真实 `σ_LS²`，这是对透明接收机的合理让步。）

**三个接收机配置**：

| 编号 | 接收机 | 窗口规则 | 假设先验 | 用途 |
|---|---|---|---|---|
| **R1（主）** | 4RB 滑窗 MMSE | 把带宽按 1RB（12 sc）分块；块 `b` 的数据 RE 用覆盖块 `b`、尽量居中的连续 4RB 窗口估计（带边截断平移，始终保持 4RB 宽）；窗口内导频 = 落在窗内的 comb 导频（每 DMRS 符号 2 个） | 均匀 PDP，`τ_w = 300 ns`（主值）；敏感性附件再跑 `τ_w ∈ {100, 694}` ns | H1 主验收指标 |
| **R2（次）** | 宽带鲁棒 Wiener | 单一全带窗口，全部 `N_p` 导频 | 均匀 PDP，`τ_w ∈ {200, 694, 1389}` ns 三档 | H2（P4 塌缩检验） |
| **R3（基准）** | V-aware matched LMMSE（= Algorithm 1） | 全带 | `R̃ = R_g`（真实协方差，UE 知道 `V`） | 只用于 QC 非透明基准点定位 |

**特殊规则**：PRG cycling 候选（§4 家族 B2）的接收机窗口**额外被 PRG 边界截断**（真实 UE 知道 PRG 大小、不跨界联合处理）；相位连续的候选（B1/CC/CN）窗口自由滑动。这条不对称正是"相位连续性"价值的来源，必须严格实现。

### 3.3 每候选附带输出的结构诊断量

- **复合时延跨度 `τ_span`**：对每列 `v_n` 做 K 点 IDFT 得时延分布 `|IDFT(v_n)|²`，叠加所有列后取包含 95% 能量的最小区间宽度（秒）。用于 H2 塌缩图横轴。
- R1 逐窗口 NMSE 的最大值（定位"最坏窗口"）。

---

## 4. 候选族（全部参数化构造，恒模）

每个带宽独立生成。分支索引 `n = 0..7`。

### B1. 透明小 delay CDD（基线 1）

`N_eff` 个不同 delay（多分支共享同一 delay 合并为一支）：分支 `n` 取 delay `τ = (n mod N_eff)·δ`，栅格取整 `j_n = round((n mod N_eff)·δ·K·Δf)`，相位 `φ[k,n] = −2πk·j_n/K`。

- 扫描：`N_eff ∈ {2, 4, 8}` × 复合跨度 `(N_eff−1)·δ ∈ {50, 100, 200, 400, 800, 1389} ns` → 18 个候选/带宽。

### B2. PRG precoder cycling（基线 2）

PRG 大小 `P_RB ∈ {2, 4, 8}` RB。第 `s` 个 PRG 内相位常数：`φ[k,n] = 2π·n·(s mod 8)/8`（逐 PRG 轮转 8 点 DFT 向量，PRG 边界处相位跳变、不连续）。→ 3 个候选/带宽。接收机窗口按 §3.2 特殊规则截断。

### B3. CC：相位连续 precoder cycling（挑战者主家族，「局部聚簇 + 全局轮转」）

参数 `(N_seg, T, schedule)`：

- 段长 `L_seg = K/N_seg`，段边界 `k = s·L_seg`。
- **段内**（`k ∈ [s·L_seg + T/2, (s+1)·L_seg − T/2)`）：常相位 `φ[k,n] = c[s,n] = 2π·n·r_s/8`（局部 group delay = 0，所有分支完美聚簇；段内等效信道为常数）。
- **过渡带**（内部边界 ±T/2）：线性斜坡 `φ[k,n] = c[s,n] + wrap(c[s+1,n] − c[s,n]) · (k − (s+1)L_seg + T/2)/T`，`wrap(·)` 取到 `(−π, π]`（最短路径，最小化过渡带局部斜率）。
- **轮转日程 `r_s`** 两种：sequential `r_s = s mod 8`；bit-reversal `r_s = [0,4,2,6,1,5,3,7][s mod 8]`。
- 扫描：48PRB `N_seg ∈ {8, 12, 16, 24}`，24PRB `N_seg ∈ {4, 6, 8, 12}`；`T ∈ {6, 12, 24}` 子载波（约束 `T ≤ L_seg/3`，不满足的组合跳过）× 2 种 schedule → 每带宽约 20 个候选。

### B4. CN：封顶斜率 N-series（挑战者次家族）

分段线性连续相位：8 段，每段每分支的栅格斜率 `j[s,n]` 独立均匀抽自 `{0, 1, ..., a}`（栅格单位，保证段内分支 group delay 聚在 `a` 个栅格步内），相位递推 `φ[k+1,n] = φ[k,n] − 2π·j[s(k),n]/K`，`φ[0,n]=0`。

- 扫描：封顶 `a ∈ {1, 2}` × 每档 50 个随机候选（种子 `23`，候选编号与种子偏移一一对应，可复现）→ 100 个候选/带宽。

### B5. QC 非透明基准点（folded-grid CDD，`DESIGN.md` §4.2 构造）

- 48PRB：`j_n = 9n`（即历史 C7，step 520.8 ns；折叠后混叠圆均匀 8 点，导频域精确正交）。
- 24PRB：`j_n = n`（step 115.7 ns；导频相位 `exp(−j2πpn/12)` 为 12 点 DFT 不同列，精确正交）。
- 用 **R3** 计算 matched NMSE、同 §3.1 计算 outage，作为基准点画在 Pareto 图上；同时也用 R1/R2 各算一遍（展示 QC 设计在透明接收机下的表现，供对照）。

### B6. 附件：P6/P7 检验集（只算 outage，不算 NMSE；仅 48PRB）

- 栅格等差 CDD：`j_n = σ·n`，`σ ∈ {1, 2, 4, 9}`。
- Sidon 型栅格集：`j_n ∈ {0, 1, 3, 7, 12, 20, 30, 65}`。
- 全部用 §3.1 的同一组 h 样本（成对比较，差值方差极小）。

候选总量约 145/带宽，加附件 5 个。

---

## 5. 仿真条件汇总

| 项目 | 取值 |
|---|---|
| 带宽 | 48 PRB（K=576）与 24 PRB（K=288），独立跑 |
| N_t / Rx | 8 Tx / 1 Rx，单层 |
| Δf | 30 kHz |
| 信道模型 | 平坦模型 `g = Vh`，`h ~ CN(0, I_8)`（理由见 §2；物理 5 ns 展宽的鲁棒性不在本轮） |
| DMRS | comb `S_f=24`，2 符号平均，`σ_LS² = σ²/2` |
| R（谱效率门限） | 2.1602 bit/RE |
| MC outage | `N_mc = 2×10^5`，种子 20260709，common random numbers |
| NMSE SNR 点 | `snr̄ ∈ {8, 16, 24}` dB |
| outage SNR 网格 | 0~24 dB，步长 0.5 dB |
| 数值防护 | 闭式公式中矩阵求逆若条件数 > 1e10，加 `1e-12·tr/dim` 对角加载并在异常节记录该候选 |

---

## 6. 校准实验（Phase 0：先跑，两个门槛都过了才跑主扫描）

- **0a（闭式 NMSE 实现正确性）**：取 2 个候选（B5 的 48PRB QC 点 + B3 的一个 CC 候选 `N_seg=12, T=12, sequential`），用已有 MC 估计器工具（`tools/search_precoder_design_alg1.py::make_alg1_full_cov_estimator`，300 trials）在 `snr̄ ∈ {8,16,24}` dB 复算 matched NMSE，与 R3 闭式值对比。**门槛：三点均一致到 ±0.3 dB。**
- **0b（outage 与 ideal-CSI BLER 的平行性，即 P3）**：对同样 2 个候选，用现有链路平台（`run.py` + YAML，参照 `outputs/` 下既有 ideal-CSI 实验配置）跑 ideal-CSI BLER（SNR 网格覆盖 10%~1% 区间，400 trials/点），与 §3.1 outage 曲线对比水平间距。**门槛：10% 与 1% 两个工作点的水平间距（coding gap）之差 < 0.4 dB。**
- 任一门槛不过：**停止主扫描**，把偏差写进 result-023 异常节并结束本轮（这说明 `DESIGN.md` §2/§3 的指标框架本身要修，比扫描更优先）。

---

## 7. 输出指标与格式

### 7.1 数据文件（`outputs/track_b_pilot_scan/<run_id>/`）

- `candidates_<K>prb.csv`：每行一个候选，列 = `id, family, 参数(N_eff/span/P_RB/N_seg/T/schedule/a/seed_offset/j_n列表), tau_span_ns, outage10_dB, outage1_dB, nmse_R1_{8,16,24}dB, nmse_R1_worst_window_16dB, nmse_R2_tw{200,694,1389}_16dB, [nmse_R3_{8,16,24}dB 仅B5], cond_flag`。
- `calibration/`：Phase 0 的对照数字与曲线。
- `figures/`：见 7.2。

### 7.2 图表（result-023 中相对路径引用并各配一句话说明）

1. **主图（每带宽一张）**：Pareto 散点，x = R1 闭式 NMSE @16 dB（dB），y = 10% outage SNR（dB）；家族用颜色区分，B1∪B2 的 Pareto 前沿画成阶梯线，B5 QC 基准点用星标；另出 1% outage 版本。
2. **塌缩图（H2）**：x = `τ_span`，y = R2 NMSE @16 dB（`τ_w=1389 ns` 档），按家族着色。
3. **附件图（H3）**：B6 各候选 outage 曲线叠画（含与等差集的成对差值曲线）。
4. Phase 0 校准对比图 2 张。

### 7.3 `research/result-023.md` 结构（固定）

TL;DR（一句话结论 + 主图指针）→ 仿真目的与验收标准（引用本 plan §8）→ 配置回执 → 关键指标表（表1：H1 支配性判定表——最优 CC/CN 候选 vs 基线前沿 vs QC 点的两轴数字；表2：H2 同跨度分组的 NMSE 极差；表3：H3 成对 outage 差值；表4：Phase 0 校准回执）→ 与预期的偏差（对照 §8 逐条）→ 异常现象（条件数、插值失稳、MC 尾部样本不足等）→ 原始数据路径。

---

## 8. 预期结论与量化验收标准

**H1（主判定，Track B 存亡）**：

- **成立**：存在 ≥1 个 B3/B4 候选，相对 B1∪B2 基线 Pareto 前沿满足其一——(i) 在 10% outage SNR 不劣于前沿点（容差 +0.05 dB）的前提下 R1 NMSE@16dB 好 **≥ 1.0 dB**；或 (ii) 在 R1 NMSE 不劣（容差 +0.2 dB）的前提下 10% outage SNR 好 **≥ 0.3 dB**。→ 下一轮（plan-024）对 top≤5 候选 + 两基线 + QC 点上 estimated-CSI BLER 链路终审。
- **不成立**：无候选达标 → Track B 在其自己的主场（滑窗透明接收机）被判死，`KNOWLEDGE.md` 固化"透明方向无结构性空间"，下一轮转向 §4.3 缺陷 regime 1（导频稀缺：`S_f ≥ 48` 或 8~12 PRB 或 `N_t=16`）或 DMRS pattern 联合设计。
- **预期**：CC 家族段内 CE 几乎无代价（局部单抽头），NMSE 应显著优于同 outage 的透明 CDD；outage 侧 `N_seg=8` 的 8 点离散轮转会比 CDD 连续扫相差一些，`N_seg` 增大应收敛——预计 (i) 路径达标。

**判据②数值预判（报告项，非本轮验收）**：报告最优透明候选与 QC 点的差距 `(ΔNMSE, Δoutage10)`。暂定"接近"参考线：ΔNMSE ≤ 2 dB 且 Δoutage10 ≤ 0.5 dB（阈值待研究者确认；最终判定以 plan-024 链路 BLER 为准）。

**H2（P4）**：把全部候选按 `τ_span` 分桶（±10%），R2（`τ_w=1389 ns` 档）NMSE 桶内极差 < 0.5 dB → 塌缩确认；若某结构在同跨度下稳定好 > 1 dB → P4 证伪，回改 `DESIGN.md` §3.4/§7。

**H3（P6+P7）**：等差集间成对 outage SNR 差 ≤ 0.05 dB → 同轨道预测成立，result-021 的 C7 尾部优势归因改判为平坦模型外效应；若 `σ=9` 复现 > 0.1 dB 优势 → P7 证伪，回改 `DESIGN.md` §2.4.4。Sidon 集 1% outage 劣于最优等差集 > 0.05 dB → P6 第二层机理存疑，记录待查。

---

## 9. 备注

- **复用**：候选构造可参考 `tools/run_v_design_piecewise_tradeoff.py` 的分段线性递推；R3/校准用 `tools/search_precoder_design_alg1.py::make_alg1_full_cov_estimator`；链路平台 `run.py` + YAML。闭式失配 MSE 与 MC outage 为**新实现**，务必先过 Phase 0 门槛。
- **计算量**：闭式 NMSE 每候选毫秒级；MC outage 每候选 `2×10^5 × K×8` 复乘，向量化后全部候选预计 CPU 数十分钟~数小时，无需 GPU。
- **公平性红线**：所有候选同一 DMRS、同一 `σ_LS²`、同一 h 样本；PRG 截断规则只对 B2 生效。R1 的 `τ_w=300 ns` 对所有候选统一（不许对某家族调优接收机）。
- **待研究者确认后再交实现 agent**：① R1 参考接收机规格（4RB 窗 + `τ_w=300 ns` 主值）；② "接近 QC"参考线（2 dB / 0.5 dB）；③ 是否需要补 36 PRB（本轮为省算力砍掉，若 H1 成立在 plan-024 补）。
