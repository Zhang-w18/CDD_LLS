# result-035：PDSCH 2Rx、60 km/h 下 4/8Tx 的 estimated/ideal CSI 与 CE NMSE

> 状态：正式实验已执行，结果待研究者确认。完整配置、原始数字、bracket、errors/trials、统计区间、限制和复现命令见同内容无图版 `research/result-035-PDSCH 2Rx仿真-text.md`。

## 1. 数据审计

本轮按 `research/plan-035-PDSCH 2Rx仿真.md` 完成 4Tx/2Rx 与8Tx/2Rx、60 km/h 的
estimated/ideal CSI paired formal。4Tx 共321,000 paired trials，8Tx 共333,000；两场景
自适应状态均为 `complete`，无 capped endpoint。逐 interval、error flag、NMSE 数组、trial
连续性、配对身份和 CPU placement 的全量审计均通过。具体哈希和审计数量见无图版第1节。

## 2. 主图

### 4Tx/2Rx

![4Tx estimated-CSI BLER](../outputs/experiment035_pdsch_2rx/20260917_fix1/formal/A100_NT4_NR2_V60/analysis/estimated_csi_bler.png)

![4Tx ideal-CSI BLER](../outputs/experiment035_pdsch_2rx/20260917_fix1/formal/A100_NT4_NR2_V60/analysis/ideal_csi_bler.png)

![4Tx estimated-CSI CE NMSE](../outputs/experiment035_pdsch_2rx/20260917_fix1/formal/A100_NT4_NR2_V60/analysis/estimated_csi_ce_nmse.png)

### 8Tx/2Rx

![8Tx estimated-CSI BLER](../outputs/experiment035_pdsch_2rx/20260917_fix1/formal/A100_NT8_NR2_V60/analysis/estimated_csi_bler.png)

![8Tx ideal-CSI BLER](../outputs/experiment035_pdsch_2rx/20260917_fix1/formal/A100_NT8_NR2_V60/analysis/ideal_csi_bler.png)

![8Tx estimated-CSI CE NMSE](../outputs/experiment035_pdsch_2rx/20260917_fix1/formal/A100_NT8_NR2_V60/analysis/estimated_csi_ce_nmse.png)

## 3. 定量结果与结论

全部12条曲线的10%/1% crossing、500次 paired bootstrap 区间、真实 bracket 和端点
errors/trials 见 `research/result-035-PDSCH 2Rx仿真-text.md` 第2节。关键 estimated-CSI 差值如下，负值表示前者更优：

| 场景 | 比较 | 10%差值 dB [95%] | 1%差值 dB [95%] |
|---|---|---:|---:|
| 4Tx | S0_SIDON − B0_QC | -0.226 [-0.291,-0.151] | +0.083 [+0.024,+0.127] |
| 8Tx | S0_SIDON − B0_QC | -0.337 [-0.403,-0.282] | -0.252 [-0.302,-0.203] |
| 4Tx | small matched − transparent | -0.149 [-0.190,-0.105] | -0.197 [-0.307,-0.112] |
| 8Tx | small matched − transparent | -0.471 [-0.546,-0.423] | -0.705 [-0.773,-0.624] |
| 4Tx | aged MRT − transparent PRG6 | +0.665 [+0.604,+0.722] | +1.284 [+1.200,+1.380] |
| 8Tx | aged MRT − transparent PRG6 | +0.451 [+0.336,+0.517] | +1.333 [+1.219,+1.436] |

结论：8Tx estimated 下 SIDON 在两个目标均优于 B0；4Tx 只在10%目标优于 B0。
matched small-delay 的优势只出现在 estimated CSI，ideal 下与 transparent 完全相同，说明优势来自
接收端协方差知识。约5 ms aged MRT 在两个场景和两种接收机下均差于 transparent PRG6。
所有方案均存在正的 estimated-vs-ideal penalty，8Tx B0_QC 最大，约1.49/1.57 dB。

NMSE、适用边界、2Rx/4Rx 尚未完成的比较、作废 prescan 记录及完整复现路径见无图版。
结果待研究者确认，尚未更新 `KNOWLEDGE.md`/`GOALS.md`。
