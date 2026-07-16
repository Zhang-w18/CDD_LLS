# 研究目标

## 北极星问题

能否设计一个频域相位预编码矩阵 `V`（CDD 是其中一种特例）及配套信道估计算法，在**相同 DMRS 开销**下，同时获得比传统全带 CDD + direct RMMSE 更好的分集增益和更低的 estimated-CSI BLER？还是说 CDD 已经是这个折中问题的实践最优解？

> 顶层理论框架、`V` 的设计空间与候选结构（分段线性/chirp/平滑随机相位/frame 设计）、理论命题见 **`DESIGN.md`**。本文件只写目标与阶段验收；做战略规划时对照 `DESIGN.md`。

## 背景（一句话版本，细节见 `docs/archive/`）

- QC 提案：显式 CDD 相位连续，允许 UE 做 wideband channel estimation，相对 PRG-level precoder cycling 有分集增益；本项目已复现该趋势并建好 link-level 仿真平台（Sionna LDPC，2Tx/4Tx/8Tx→4Rx/1Rx，TDL 信道，可配置 CE 算法）。
- 在此基础上把问题推广到任意频域相位矩阵 `V`：`V` 的列正交性（分集）和频域平滑性（信道估计处理增益）天然冲突。CDD 是"每个分支只有一个全带固定 delay"约束下的自然最优解，但不是一般相位矩阵空间里的全局最优。

## 当前阶段（从 plan-022 起）

上一阶段（实验 1–21，详见 `KNOWLEDGE.md`）证明了：分段线性相位 N-series 设计确实能在**矩阵层面**和 **ideal-CSI BLER** 上超过 CDD，但这个优势在**稀疏导频 + Algorithm 1 RMMSE** 的 estimated-CSI 链路上消失甚至倒扣。

当前阶段目标：不再只用矩阵指标（log-det、B_0.5）筛选 `V`，而是把导频可估计性直接纳入设计目标，验证联合目标

```
J(V) = J_div(V) - mu * L_CE(V; 实际DMRS pattern)
```

是否能筛出「ideal-CSI 分集增益保留、estimated-CSI BLER 不倒扣」的候选。

## 验收标准（决定本阶段成败）

- **假设成立**：存在候选 `V`，在与 CDD 相同 DMRS 开销/密度下，estimated-CSI 10% BLER 所需 SNR 不高于最优 CDD 参考点（容差 ±0.2 dB，即在 400-trial 分辨率内不可分辨也算持平）。
- **假设不成立**：扫描过 `mu` 的合理范围（见 plan-022）后，所有候选的 estimated-CSI BLER 仍全部劣于 CDD → 结论固化为"CDD 在当前稀疏导频体制下是实践最优"，转向探索联合优化 DMRS pattern 本身（而不只是 `V`）。

## 不在当前阶段范围内

- UE 移动性 / Doppler / time-varying TDL（仍是 Phase B，未启动）。
- 宽带 massive-MIMO 预编码降维到 effective ports 后再叠加 CDD（后续扩展方向，未启动）。
