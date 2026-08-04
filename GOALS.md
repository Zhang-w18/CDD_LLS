# 研究目标

## 1. 总体研究问题

能否设计频域相位预编码矩阵 `\mathbf V` 及配套信道估计方法，在相同 DMRS 开销下，相对等差 CDD 获得更低的 estimated-CSI BLER，并明确结论成立的物理信道、带宽、接收机知识和失配范围？

`\mathbf V` 是有效子载波到发射分支的恒模频域相位矩阵，CDD 是其线性相位特例。系统模型和符号见 `DESIGN.md`。

## 2. 背景

项目已建立 Sionna LDPC 链路级仿真平台，支持 static TDL、2/4/8 Tx、1/4 Rx、CDD/PRG/一般 `\mathbf V`、多种信道估计器以及 ideal-CSI 和 estimated-CSI BLER。

既有实验说明：矩阵指标和 ideal-CSI 增益不能替代真实 DMRS 条件下的估计与 BLER；分段线性 N-series 在已测稀疏导频 estimated-CSI 条件下整体劣化；改变参考信号结构时必须核算 DMRS 开销。详见 `KNOWLEDGE.md`。

## 3. result-024 后的进展

固定 4RB V-agnostic 接收机下，CC/CN 分段候选没有支配透明 CDD/PRG 基线。接收机获得分段边界后，48 PRB 仍无候选通过，24 PRB 仅有一个半透明诊断点通过。该方向当前不作为主要后续方向。

在 48 PRB、8 Tx / 1 Rx、平坦分支信道、DMRS comb 24、相同开销和双方 matched LMMSE 条件下，Sidon delay 索引 `[0,1,3,7,12,20,30,65]` 相对等差 QC CDD：

- 10% BLER 所需 SNR 改善 0.33 dB，保守 95% 区间 [0.20, 0.45] dB；
- 1% BLER 所需 SNR 改善 0.85 dB，保守 95% 区间 [0.63, 1.07] dB；
- matched CE NMSE 最大差 0.06 dB。

因此，“存在一般 `\mathbf V` 在同 DMRS 开销下优于等差 CDD”的假设已在上述限定条件内成立。该结论不能直接推广到频率选择性 TDL、定时误差、协方差失配、其他带宽或移动信道。

result-025 的 delay-matched `/576` 补做只把底层平坦分支信道替换为 Sionna TDL-A 5 ns，并保持双方使用各自真实复合协方差的 matched LMMSE。该条件下 Sidon 相对 QC 的 10% BLER 改善为 0.156 dB，保守 95% 区间 [0.045, 0.266] dB；1% 改善点估计为 0.090 dB，保守 95% 区间 [-0.210, 0.391] dB。10% 小幅优势得到确认，1% 排序不能确定。TDL 加权后二阶相关代理变为 QC 略低，Sidon 四阶代理仍明显较低；Sidon 的 data-RE CE NMSE 较 QC 差 0.754–1.271 dB。该结果只适用于 TDL-A 5 ns、零速度、当前 DMRS/MCS 和 known-delay matched 接收机。

## 4. 当前阶段目标

在 result-025 已确认的基础上，寻找物理展宽信道下更有潜力的结构化 CDD 时延设计，并继续确定原 Sidon 优势的适用范围：

1. 在完整物理协方差已知时联合设计 CDD residue/lift，并检查分集代理、matched CE 和 outage 的 Pareto 前沿；
2. 在发射端只知道有效时延支撑时，探索厚 Sidon 的可行条件和潜力区域；
3. 在近似平坦 TDL 与 LoS TDL 条件下复核原 Sidon 的 estimated-CSI BLER；
4. 后续检查协方差失配、定时误差和 24 / 36 / 48 PRB 一致性；
5. 汇总保持优势、持平和反向的物理条件。

新 plan 的阈值、PDP 集合、trial 数和停止条件必须在实施前与研究者确认。

## 5. 阶段验收要求

1. 与等差 QC CDD 使用相同 DMRS、开销、MCS、trial 和随机样本；
2. 每个条件报告 10% BLER，样本允许时报告 1% BLER；
3. 报告错误块计数、目标 SNR 方法和置信区间；
4. 同时报告导频矩阵条件数或最小奇异值、matched CE NMSE 和 BLER；
5. 区分 matched 与 mismatched 接收机；
6. 报告全部预先规定条件；
7. 完成稳健性检查前，结论必须包含当前限定条件。

判定分为：

- **稳健成立**：主要预定条件内，10% BLER 优势统计区间保持为正且无未解释 CE 失稳；
- **条件成立**：只在部分条件保持优势，明确条件集合和失效条件；
- **不成立**：主要实际条件下优势稳定反向，更新 `KNOWLEDGE.md` 并停止把该 Sidon 集作为主要候选。

具体数值门槛由下一份经研究者确认的 plan 固定，不得在看到正式结果后修改。

## 6. 暂不纳入本阶段

- UE 移动性、Doppler、time-varying TDL；
- massive-MIMO 预编码降维后再叠加 CDD；
- 标准化、信令开销和硬件实现评估；
- 与 CDD residue/lift 无关的一般高维 `\mathbf V` 大规模搜索。

## 7. 最新证据

`research/plan-025.md`、`research/result-025.md`、`research/result-025-text.md`、`outputs/experiment025_sionna_tdl_rmmse/20260723_delay_matched/`。
