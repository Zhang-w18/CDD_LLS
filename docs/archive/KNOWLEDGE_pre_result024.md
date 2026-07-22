# KNOWLEDGE —— 智能体记忆体（导航式）

> **用法**：本文件是给分析 agent 每轮读的精简记忆。按下面「目录」定位需要的小节即可，不必逐字通读。
> 每条结论只写「一句话主张 + 一句话原因 + 指针」。需要**数字证据/图表**看给人的详尽文档 `docs/FINDINGS.md`（不必每轮读）；需要**原始推导**看 `docs/archive/`。
> `[E##]` 是历史实验编号，仅供在 FINDINGS/archive 中检索定位，本文件不展开。

## 当前状态（一句话）

CDD（含大 delay 变体）在当前 8Tx/1Rx + 稀疏 comb DMRS 体制下是 estimated-CSI BLER 的实践最优；分段线性相位 N-series 只在矩阵层面和 ideal-CSI 层面赢，一进真实稀疏导频估计就倒扣。plan-022 正在验证「导频感知联合目标」能否翻盘。

## 目录

- **A. 平台与算法命名** — 跑什么、算法怎么叫（改动实现前先看）
- **B. 已确立结论** — B1~B8，每条一句话主张
- **C. 已排除的死路** — 不要重复尝试的方向
- **D. 当前开放问题** — plan-022 及之后要回答的
- **E. 历史文档索引** — 需要溯源时才打开

---

## A. 平台与算法命名

- **平台**：Sionna LDPC + 自研 link-level 链路，YAML 配置驱动，`run.py` 入口，输出到 `outputs/<name>/sim_*/`。
- **已实现**：static TDL、单层 PDSCH、2/4/8 Tx → 1/4 Rx、CDD/PRG/一般 `V` 预编码、IDEAL / RMMSE(4RB & wideband) / reconstruction 系列 CE、Sionna LDPC bit-level BLER。
- **未实现（暂不在范围）**：UE 移动性、Doppler、time-varying TDL、宽带 massive-MIMO 预编码降维。
- **算法统一命名**（后续一律沿用）：
  - **Algorithm 1** — Direct equivalent-channel RMMSE。CDD 用 shifted-PDP Toeplitz 协方差；一般 `V` 用完整非平稳协方差 `R_g = R_phy ⊙ (CC^H)`。
  - **Algorithm 2B** — multi-pilot delay-domain basis LMMSE：先估底层物理分支信道，再重构等效信道。
  - **Algorithm 2C** — local-window deterministic：在相干带宽窗口内把底层信道当常数。
  - **Algorithm 3** — non-CDD per-port DMRS：导频不加 CDD，正交估各端口物理信道，再用已知 delay 合成等效信道。

---

## B. 已确立结论

> 每条：主张 → 原因 → 证据指针。数字与图在 `docs/FINDINGS.md` 同名小节。

### B1. Algorithm 1 是最强通用 baseline
只要所有估计器共享同一组 CDD-combined DMRS 观测，matched direct RMMSE 就最强。它直接估计检测器需要的等效信道、用上了已知 CDD 统计、不做多余的分支分离。证据：FINDINGS §B1。

### B2. Algorithm 2B/2C 在统一网格上没有稳健增益
最佳观测增益仅 +0.014 dB，处于有限 trial 噪声内。根因是**观测受限**：从单层 CDD-combined 导频反推多个物理分支 tap，敏感矩阵 `Phi` 条件数 1e12~1e13、nullity 可达未知量 70%+，重构不增加独立观测。证据：FINDINGS §B2。

### B3. Algorithm 3 在「稀疏导频 + 大 CDD delay」能拿大增益
最大 +18.237 dB（delay=512 samples、DMRS spacing=24 时 15 个组合全正）。原因是它**改了参考信号设计本身**（正交估各端口），不是换估计器。代价：额外 DMRS 开销 / 降每端口密度，比较必须核算开销公平性。证据：FINDINGS §B3。

### B4. N-series 分段线性相位在矩阵层面确实改善 Pareto 前沿
相位连续的局部 slope 切换（尤其带符号非均匀字母表 + 拉丁方布局）在分集/相干带宽平面上优于常规 CDD；ideal-CSI BLER 也观测到增益（48PRB 下 10%/1% 增益达 +0.36~+1.55 dB）。证据：FINDINGS §B4。

### B5.【核心痛点】ideal-CSI 分集增益换不到 estimated-CSI BLER 增益
稀疏 DMRS + Algorithm 1 下，N-series 的 matched RMMSE NMSE 在 20dB 处倒扣 1.81~3.26 dB（常规候选），24PRB N8 倒扣 8.04 dB 并出现 ~11-15% BLER floor。原因：N-series 制造的**非平稳等效协方差**在稀疏导频下比 CDD 的平稳 shifted-PDP 更难插值。**这就是"结果不符合预期"的确切位置——是折中的真实体现,不是 bug。** 证据：FINDINGS §B5。

### B6. 矩阵指标只能筛选,不能当最终目标
`log det(V^H V)`、`B_0.5`（相干带宽）都无法决定完整非平稳协方差在稀疏导频下的可估计性,也决定不了编码 BLER。任何候选 `V` 必须过两关：matched full-cov RMMSE NMSE + estimated-CSI 编码 BLER。证据：FINDINGS §B6。

### B7. CDD delay 增大对 BLER 非单调
等效信道相干带宽随 delay 增大而压缩,是否更优取决于与 DMRS spacing 的交互,存在非平凡最优点。不要假设"delay 越大越好"。证据：FINDINGS §B7。

### B8. 接收机 LLR 未注入 CE 误差自噪声（潜在修正点）
当前解调 LLR 的噪声方差没有加入信道估计误差项。这可能是大 NMSE 候选（如 24PRB N8）出现 BLER floor 而非平滑退化的部分原因。尚未验证,列为候选修正。证据：FINDINGS §B8。

---

## C. 已排除的死路（不要重复尝试）

1. **单纯调参救 Algorithm 2B/2C 的数值病态** —— 已验证无效的有：缩短 delay support（90/95/99% 能量）、加大 diagonal loading（1e-6~1e-2）、pairwise 解耦（条件数同样爆炸）。除非**改变导频观测方式本身**（如走 Algorithm 3 思路），否则别再动这些旋钮。见 B2。
2. **用矩阵指标本身作为 `V` 设计的最终验收** —— 必须过 RMMSE NMSE + BLER 关。见 B5、B6。
3. **假设"CDD delay 越大分集越好"** —— 必须结合 DMRS spacing 联合看。见 B7。

---

## D. 当前开放问题

- **[plan-022 主问题]** 联合目标 `J(V) = J_div(V) - mu * L_CE(V; 实际DMRS pattern)` 能否筛出「保留 ideal-CSI 分集增益、且 estimated-CSI BLER 不劣于 CDD」的候选？验收标准见 `research/plan-022.md` §7 与 `GOALS.md`。
- **[plan-023 主问题，Track B 先导]** 在固定 4RB 滑窗透明接收机下，「局部聚簇 + 全局轮转」（相位连续 cycling / 封顶 N-series）能否在（闭式失配 NMSE，MC outage）平面严格支配透明 CDD 与 PRG cycling（DESIGN §7 P5，Track B 存亡判定）？附带检验 P4（宽带窗塌缩）、P6/P7（等差同轨道 / Sidon 尾部）。纯数值扫描、与 plan-022 并行，见 `research/plan-023.md`。**待研究者确认参考接收机规格与"接近 QC"阈值后再交实现 agent。**
- **[失败分支]** 若联合目标筛完仍全劣于 CDD → 转向"联合优化 DMRS pattern 本身"（而非只调 `V`），或转 Algorithm 3 类"改参考信号"方向（见 B3）。
- **[候选修正,不阻塞]** 接收机 LLR 是否应加 CE-error-aware 噪声项?会不会改变 N8 的 BLER floor 结论?见 B8,暂未排期。

---

## E. 历史文档索引（按需溯源,默认不读）

| 想找什么 | 去哪 |
|---|---|
| **顶层理论框架 / `V` 设计空间 / 度量指标 / 两条赛道**（做战略规划时对照） | **`DESIGN.md`**（v2,顶层,非归档） |
| v1 旧设计文档（Gram/平滑性旧框架,已被 v2 取代;候选结构原始细节仍在此） | `docs/archive/DESIGN_v1.md` |
| 结论 B1~B8 的数字证据 + 图表（给人看的详尽版） | `docs/FINDINGS.md` |
| QC 系统模型、Algorithm 1/2B/2C/3 完整数学推导 | `docs/archive/QC_explicit_CDD_link_level_simulation_plan.md` |
| 实验 1–19 完整过程记录（含病态诊断细节） | `docs/archive/experiment_record_20260608.md` |
| N-series 矩阵构造、Pareto 扫描细节 | `docs/archive/v_design_piecewise_tradeoff_experiment.md` |
| 结论 B1~B8 的英文综合出处 | `docs/archive/cdd_nontransparent_channel_estimation_and_frequency_domain_precoder_report.md` |
