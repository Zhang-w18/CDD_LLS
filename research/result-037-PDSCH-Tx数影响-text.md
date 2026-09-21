# result-037：PDSCH Tx数影响（4T1R 部分结果，无图版）

> 状态：**4T1R 正式结果完成；16T/32T 未纳入本版结论**。研究者在实现阶段把接收天线数从
> plan-037 草案中的 2Rx 改为 1Rx；本结果以实际展开配置 `A100_NT4_NR1_V60` 为准。
> 未经研究者确认，本结果不更新 `KNOWLEDGE.md` 或 `GOALS.md`。

## 1. 结论

4T1R formal 自适应状态为 `complete`，无 capped endpoint。estimated CSI 下：

- B0 相对 PRG6 的 10%/1% BLER 目标 SNR 增益为 0.452 dB
  （95% paired-bootstrap CI `[0.365,0.582]` dB）和 1.102 dB（`[0.940,1.282]` dB）；
- Sidon 相对 PRG6 的对应增益为 0.595 dB（`[0.484,0.705]` dB）和
  1.216 dB（`[1.054,1.391]` dB）。

四个主比较区间均完全大于 0，因此支持当前冻结 4Tx B0、Sidon 相对 PRG6 的正增益。
该结论只适用于 4Tx/1Rx、TDL-A 100 ns、60 km/h 和当前接收机，不支持跨 Tx 趋势或纯 $N_t$
因果解释。

## 2. 配置、数据与复现

实际场景为 4Tx/1Rx、1 layer、TDL-A 100 ns、3.5 GHz、60 km/h、20 sinusoids；48 PRB、
576 active subcarriers、30 kHz SCS；DMRS symbols `[2,7]`、comb-6；16QAM、NR 256QAM table
MCS 8、目标码率 553/1024。estimated 接收机使用两个 DMRS 的二维时频 RMMSE，ideal 接收机使用
真实 data-RE 等效信道。SNR 网格为 13--20 dB、0.25 dB 间隔，base seed 为 `20260727`，
batch size 为 25。各点实际使用 1,000--30,000 paired trials。

冻结构造：B0 `[0,18,36,54]`；严格 Sidon `[0,24,72,240]`；PRG6 使用 DFT4 index
`[0,1,2,3,0,1,2,3]`。逐 RE 总发射功率为 1。B0/Sidon pilot rank 均为 4，condition number
约为 1；Sidon 的 10 个模 576 无序二元和全异。

正式配置：`configs/plan037_nt4_nr1_v60_formal.yaml`，SHA-256：
`929e716e043f7b46faa53de46e0e5734671d77e10387ba3b3b87eb38606bc384`。

原始与展开数据根目录：
`outputs/experiment037_pdsch_tx_scaling/20260920/formal/A100_NT4_NR1_V60/`。

复现分析：

```powershell
python tools/analyze_result037_pdsch_tx_scaling.py --scene outputs/experiment037_pdsch_tx_scaling/20260920/formal/A100_NT4_NR1_V60 --bootstrap-replicates 1000
```

## 3. Crossing、增益与接收机 penalty

| 方案 | CSI | BLER | crossing SNR (dB) | crossing 95% CI (dB) | 相对 PRG6 增益 (dB) | 增益 95% CI (dB) |
|---|---|---:|---:|---:|---:|---:|
| B0 | estimated | 10% | 15.098 | [14.972,15.165] | 0.452 | [0.365,0.582] |
| Sidon | estimated | 10% | 14.956 | [14.854,15.054] | 0.595 | [0.484,0.705] |
| PRG6 | estimated | 10% | 15.550 | [15.495,15.587] | 0 | reference |
| B0 | estimated | 1% | 17.008 | [16.949,17.071] | 1.102 | [0.940,1.282] |
| Sidon | estimated | 1% | 16.895 | [16.841,16.954] | 1.216 | [1.054,1.391] |
| PRG6 | estimated | 1% | 18.111 | [17.966,18.278] | 0 | reference |
| B0 | ideal | 10% | 13.881 | [13.790,13.987] | 0.801 | [0.661,0.925] |
| Sidon | ideal | 10% | 13.765 | [13.693,13.841] | 0.916 | [0.804,1.039] |
| PRG6 | ideal | 10% | 14.682 | [14.606,14.784] | 0 | reference |
| B0 | ideal | 1% | 15.666 | [15.616,15.716] | 1.596 | [1.485,1.672] |
| Sidon | ideal | 1% | 15.529 | [15.460,15.591] | 1.733 | [1.616,1.825] |
| PRG6 | ideal | 1% | 17.262 | [17.159,17.317] | 0 | reference |

estimated-minus-ideal penalty（10%/1%）分别为：B0 1.217/1.342 dB，Sidon
1.190/1.366 dB，PRG6 0.869/0.848 dB。对应 95% CI 见
`analysis/estimated_minus_ideal_penalties_bootstrap.csv`。

完整 crossing bracket、端点 errors/trials、bootstrap 有效数：
`analysis/target_crossings_bootstrap.csv`。增益数据：`analysis/target_gains_vs_prg6_bootstrap.csv`。
三类曲线的逐点文字数据为 `analysis/estimated_csi_bler.csv`、`analysis/ideal_csi_bler.csv` 和
`analysis/estimated_csi_ce_nmse.csv`。上述路径均相对于本节给出的正式输出根目录。

## 4. 解释、完整性与未完成项

**事实：** Sidon 的 estimated/ideal crossing 点估计均最低；本轮预声明主比较为各方案相对 PRG6，
因此不声明 Sidon 相对 B0 的统计显著性。

**推断：** CDD 的 ideal 增益较大，但 estimated-minus-ideal penalty 也比 PRG6 大，说明额外 CE
损失抵消了部分分集收益；在 4T1R 下仍保留正的 estimated-CSI 净增益。NMSE 只作机理诊断，
不单独决定胜负。

数据审计核验了 87 个 estimated 点、87 个 ideal 点、1,446 条 interval 和 2,169 个数组；没有
发现 shape、错误计数、非有限值或 absolute-trial 间隙/重叠。运行使用 CPU，峰值 RSS
2,756,280,320 bytes。12 个 crossing 及全部 gain/penalty 均有 1,000/1,000 个有效 paired-bootstrap
重复。

未完成：16T/32T formal、跨 Tx independent bootstrap、跨 Tx 汇总图，以及 plan 中最终完整结论。
plan 的 2Rx→1Rx 文字修订也尚待同步。

