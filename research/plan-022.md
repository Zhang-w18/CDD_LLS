# plan-022：导频感知联合目标筛选 V，验证能否同时保留分集增益与 estimated-CSI BLER

## 背景（自包含，无需查历史文档）

系统研究比较不同频域相位预编码矩阵 `V`（CDD 是其特例）在分集增益与信道估计代价之间的折中。上一轮实验（`research/result-021.md`）用矩阵指标（log-det、相干带宽 `B_0.5`）筛出的"N-series"分段线性相位候选，在 ideal-CSI BLER 上确实优于 CDD，但换算到 estimated-CSI（用同一组稀疏 DMRS 做 Algorithm 1 matched full-covariance RMMSE）后，NMSE 全部倒扣（20dB 处 +1.8~+8.0 dB），BLER 没有一个候选稳定超过 CDD。

本轮改变筛选目标本身：不再只用矩阵代理指标筛选，而是把**实际 DMRS pattern 下的 matched RMMSE NMSE** 直接纳入筛选目标，验证是否存在候选在这个更真实的目标下仍能超过 CDD。

## 系统参数（与 result-021 保持一致，便于直接对比）

| 参数 | 取值 |
|---|---|
| 物理信道 | Static TDL，指数 PDP，delay spread = 5 ns |
| 天线 | 8 Tx / 1 Rx，单层 |
| SCS / FFT size | 30 kHz / 4096（采样率 122.88 MHz） |
| PDSCH 时域 | 10 OFDM symbols |
| 带宽 | 48 PRB (576 sc) / 36 PRB (432 sc) / 24 PRB (288 sc)，三个带宽独立跑 |
| Segment 数 | 8（每段长度 = 带宽/8，即 72/54/36 subcarriers） |
| DMRS（主实验） | symbol [2, 7]，频域间隔 `S_f = 24` subcarriers，两个 DMRS symbol 合并平均，`sigma_LS^2 = sigma_w^2 / 2` |
| DMRS（鲁棒性子实验，见第 5 节） | 同上，仅 `S_f = 12` |
| MCS | NR 256QAM table MCS 8，实际调制 16QAM，code rate 553/1024 |
| LDPC | 8 iterations |

## 1. 候选池生成

复用 `tools/run_v_design_piecewise_tradeoff.py` 中已有的连续分段线性相位构造（相位连续、按段整数 slope 递推，见 `docs/archive/v_design_piecewise_tradeoff_experiment.md` §3A.4 "PW cont random slope"）。每个带宽独立生成候选池：

- 候选数：每个带宽至少 200 个随机 slope 组合候选（8 段 × 8 发射分支的 slope 矩阵随机采样），复用已有的 slope alphabet 采样方式，不需要重新设计生成算法。
- 额外固定加入的参考点（每个带宽都要有）：
  - `result-021.md` 中该带宽表现最好的 N-series 候选（例如 48PRB 的 N6）；
  - `result-021.md` 中该带宽匹配 `B_0.5` 的 CDD 参考点（例如 48PRB 的 C0/C1/C2）；
  - `result-021.md` 异常现象里提到的大 delay CDD 参考 `C7`（`step=64 samples`），因为它在 ideal-CSI 1% BLER 尾部表现最好，之前被矩阵匹配方法漏掉，本轮必须显式纳入对比，不能只用 `B_0.5` 最接近原则重新筛一遍。

## 2. 联合目标与筛选

对候选池中每个 `V`（含参考点），计算两个量：

1. **矩阵分集指标** `J_div(V) = log det((1/K) V^H V + eps*I)`，`eps = 1e-9`，`K` 为该带宽 active subcarrier 数。（沿用已有实现。）
2. **导频感知 CE 代价** `L_CE(V)`：用**实际 DMRS pattern**（主实验 `S_f=24`，symbol [2,7]）算 matched full-covariance RMMSE NMSE，SNR 取 `{8, 16, 24}` dB 三点各 300 trials（复用 `tools/search_precoder_design_alg1.py::make_alg1_full_cov_estimator`），取三点 NMSE(dB) 的均值作为 `L_CE(V)`。这是本轮和上一轮"只看矩阵指标"的关键区别——筛选阶段就用真实 RMMSE NMSE，不是等到最后验证阶段才发现代价。

联合得分：

```
J(V) = J_div(V) - mu * L_CE_dB(V)
```

`mu` 取 `{0, 0.1, 0.3, 1.0, 3.0}` 分别筛一遍（`mu=0` 退化为纯分集筛选，作为对照）。每个 `(带宽, mu)` 组合选出 `J(V)` 最高的 5 个候选。

**输出要求**：每个带宽输出一张候选池全量的 `(J_div, L_CE_dB, J at each mu)` CSV，一张 `L_CE_dB vs J_div` 散点图（标出被各 `mu` 选中的点）。

## 3. 候选去重与验证集

把 5 个 `mu` × 3 个带宽选出的候选取并集去重，加上第 1 节的固定参考点，形成最终验证集（预计每带宽 10-20 个候选，含 CDD 参考）。

## 4. 链路级验证（主实验）

对验证集中每个候选、每个带宽，运行：

- **Ideal-CSI BLER**：SNR grid `[-6, -4, -2, 0, 2, 4, 6, 8, 10, 12]` dB，400 trials/点，或达到 100 个 TB error 提前停止（取先满足者）。
- **Estimated-CSI BLER（Algorithm 1 matched full-covariance RMMSE）**：同一 SNR grid，同样的 400-trial / 100-error 停止规则，DMRS pattern 用主实验 `S_f=24`。
- 两者都要输出 10% BLER 和 1% BLER 插值 SNR。

## 5. 鲁棒性子实验（次要，规模更小）

取第 3 节验证集中综合表现最好的 3 个非 CDD 候选 + 1 个最强 CDD 参考，在 `S_f=12`（更密导频）下重复第 4 节的 estimated-CSI BLER（不需要重跑 ideal-CSI，因为 ideal-CSI 不依赖 DMRS pattern）。目的：判断"联合目标筛选是否只在稀疏导频下失效，加密导频后是否能挽回"。

## 6. 输出格式要求（写入 `research/result-022.md`，严格遵循此结构）

本结果既供分析 agent 迭代、也供研究者本人审阅,因此要求：**TL;DR 结论前置**;关键指标用表格**对照计划里的预期**;异常单列;所有图表在 `research/result-022.md` 里用相对路径引用并配一句话说明（不要只丢路径）。图表本身放 `outputs/<name>/.../figures/` 或 `docs/figures/`,BLER 曲线务必标注 10%/1% 目标线、CDD 参考曲线用同色系区分候选。

```markdown
# result-022
## TL;DR（一句话结论 + 最关键的 1~2 张图指针）

## 配置回执（实际使用的参数）
（列出实际跑的 mu 值、候选数、SNR grid、trial 数等，与本计划确认一致或说明偏差）

## 关键指标（表格，对照计划里的预期值）
- 表1：各带宽下，验证集候选的 estimated-CSI 10%/1% BLER 相对最强 CDD 参考的增益（dB，正值=候选更好）
- 表2：各带宽下，验证集候选的 matched RMMSE NMSE（20dB 处）相对 CDD 的倒扣（dB）
- 表3：S_f=12 鲁棒性子实验的 estimated-CSI 10% BLER 增益

## 与预期的偏差
（对照第 7 节验收标准逐条判断）

## 异常现象
（数值病态、条件数异常、trial 数不足导致的不可靠点等）

## 原始数据路径
（raw/ 或 outputs/ 下的实际路径，含候选池 CSV、散点图、BLER 曲线图）
```

## 7. 验收标准（量化，判断路径要短）

- **假设成立**：存在至少一个非 CDD 候选，在某个带宽 + 主 DMRS pattern（`S_f=24`）下，estimated-CSI 10% BLER 所需 SNR ≤ 该带宽最强 CDD 参考（含 C7 大 delay 点）+ 0.2 dB。→ 记录该候选的 `(mu, J_div, L_CE)`，作为下一轮设计的种子。
- **假设部分成立**：仅在 `S_f=12`（更密导频）子实验中出现上述优势，`S_f=24` 下没有。→ 结论固化为"联合目标只在导频密度足够时有效，需要联合优化 DMRS 密度而非只调 V"，下一轮转向该方向。
- **假设不成立**：所有 `mu`、所有带宽、两种导频密度下都没有候选达标。→ 结论固化为"在当前 8Tx/1Rx + 稀疏 comb DMRS 体制下 CDD（含大 delay 变体）是实践最优", 下一轮把研究重心转向 DMRS pattern 本身的联合设计，或转向 Algorithm 3 类"改变参考信号"的方向（见 `KNOWLEDGE.md` 已确立结论第 3 条）。

## 备注

- 沿用 common random numbers：同一 SNR 点下所有候选共享相同 channel/payload/noise realization，减少候选间比较方差。
- 若某候选联合协方差矩阵条件数超过 `1e10`，需要在 result-022 的异常现象里单独标出，不要让它悄悄拉低平均 NMSE 掩盖其他候选的真实表现。
