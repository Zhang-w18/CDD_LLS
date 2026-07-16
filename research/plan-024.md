# plan-024：分段边界已知接收机诊断，以及 Sidon / QC 的 V-aware 链路终审

## 1. 目的与结论边界

本轮回答两个相互独立的问题。

1. **E1（分段估计诊断）**：plan-023 的 B3/B4 候选在 4RB 滑窗跨越分段边界时产生较大失配。若接收机知道候选的频域分段边界，并禁止跨段联合估计，闭式 NMSE 是否显著下降？
2. **E2（Sidon 链路终审）**：H3 在平坦模型的块互信息 outage 上发现 Sidon delay 集优于 QC 等差 delay CDD。两者均采用 V-aware matched LMMSE 时，该优势是否传递到实际有限码长 LDPC estimated-CSI BLER？

E1 的接收机使用候选的 `N_seg` 和边界，因此不再是 plan-023 H1 的完全 V-agnostic R1，而是**边界已知、相位值未知的半透明诊断接收机**。E1 若改善，只能说明 H1 的失败部分来自窗口与边界不对齐；不能直接恢复“零信令透明方案”的结论。

## 2. E1：分段对齐 MMSE（R4）

### 2.1 系统与候选

- 与 plan-023 主扫描一致：`N_t=8`、`Δf=30 kHz`、平坦模型 `g=Vh`、DMRS comb `S_f=24`、两个 DMRS 符号平均、`snr̄=16 dB`。
- 带宽：48 PRB（`K=576`）与 24 PRB（`K=288`）。
- 候选：plan-023 的全部 B3（CC）和 B4（CN）。
- B3 使用候选自己的 `N_seg`；B4 固定 `N_seg=8`。

### 2.2 R4定义

把第 `s` 段定义为

```text
W_s = {sL, ..., (s+1)L-1},  L=K/N_seg.
```

目标数据为该段全部非导频子载波 `D_s`，观测只使用该段内的导频 `P_s`。估计器仍不知道真实 `V`，使用与 R1 相同的均匀 PDP 先验 `τ_w=300 ns`：

```text
A_s = Rtilde_DsPs (Rtilde_PsPs + sigma_LS^2 I)^-1,
g_hat_Ds = A_s y_Ps.
```

真实 MSE 用 `R_g=VV^H` 的子块闭式计算。R4不跨段，但B3的线性过渡带位于边界两侧，因此过渡带数据仍计入相邻段，不能通过删除困难RE虚增性能。

### 2.3 输出与判据

- 每候选输出 `R1_NMSE`、`R4_NMSE`、`gain_R4_vs_R1_dB = NMSE_R1(dB)-NMSE_R4(dB)`、每段导频数最小/最大值。
- 报告各家族改善候选比例、中位改善和最大改善。
- 使用 plan-023 已保存的 outage SNR，重新检查挑战者是否能在“R4诊断坐标”下越过原 B1∪B2 R1 Pareto 前沿。此比较仅为诊断性上界，不作为公平产品判定。
- **有继续价值**：至少一个 B3/B4 候选改善 ≥1 dB，且在诊断坐标下满足 plan-023 H1 的任一支配门槛。

## 3. E2：Sidon与QC的V-aware estimated-CSI BLER

### 3.1 候选与公平性

仅跑48 PRB：

- `QC_arith_s9`：`j_n=[0,9,18,27,36,45,54,63]`；
- `Sidon`：`j_n=[0,1,3,7,12,20,30,65]`。

两者均为恒模DFT-grid CDD，满足全带列正交，且 `j_n mod 24` 互异，因此导频域列正交。两者使用完全相同的：

- DMRS pattern与开销；
- 平坦 `h~CN(0,I_8)` 样本、payload、导频噪声和数据噪声（common random numbers）；
- 真实 `R_g=VV^H` 构造的全带 matched LMMSE；
- 16QAM、MCS 8（码率553/1024）、8次LDPC迭代。

### 3.2 链路流程

每个trial：

1. 抽取一组平坦8分支信道 `h`；
2. 对每个候选计算 `g=Vh`；
3. 在导频位置加入方差 `sigma_LS^2=(N_t/snr̄)/2` 的LS噪声；
4. 用候选自己的matched LMMSE矩阵估计全带 `g_hat`；
5. 在数据RE发送同一编码TB，加入相同数据噪声；
6. 使用 `g_hat` 做单层均衡、16QAM软解调和LDPC译码。

SNR网格初值 `13.0:0.5:18.0 dB`。快速主扫描每点400 trials；若运行时间允许，对10%/1%交叉附近加密。所有结论必须同时报告错误块数和二项统计不确定性，样本不足时标为“先导、不可终审”。

### 3.3 判据

- 主判据：Sidon的10% BLER SNR相对QC改善 ≥0.15 dB，且交叉区置信区间不支持反向排序。
- 次判据：报告1% BLER SNR差；若目标附近任一候选累计错误块 <30，不做确定性1%结论。
- 同时报告matched CE NMSE差，验证性能差异不是由信道估计条件不公平造成。

## 4. 输出

固定目录：`outputs/experiment024_segment_sidon_qc/20260716_main/`。

- `segment_nmse.csv`
- `segment_summary.json`
- `sidon_qc_bler.csv`
- `sidon_qc_summary.json`
- `figures/segment_nmse_gain.png`
- `figures/sidon_qc_bler.png`
- `research/result-024.md`

复现入口：`tools/run_experiment024_segment_sidon_qc.py`。
