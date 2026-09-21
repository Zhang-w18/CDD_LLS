# plan-030：PDCCH AL8 Polar repetition 与 BLER 验证

## 1. 目标

在 result-029 的 Sionna PDCCH 分支上增加 aggregation level 8 端到端能力。固定 DCI payload A=41、CRC24C/RNTI `0xFFFF`，验证 `K=65`、Polar mother `N=512`、rate-matched `E=864` 的标准 repetition、LLR 合并和 estimated-CSI BLER 曲线。本轮不比较 CDD 时延或不同 REG bundle 方案。

## 2. 固定物理口径

- 48 PRB、2 CORESET symbols、30 kHz、FFT 4096、CP 288；8 Tx/1 Rx。
- non-interleaved CCE–REG mapping，first CCE 0，AL8，REG bundle L=6；候选为 8 CCE、48 REG、8 bundles。
- 每 REG 为 9 data+3 DMRS RE，因此 data/DMRS 为 432/144 RE，QPSK `E=864`。
- A=41，CRC24C 24 bits，RNTI `0xFFFF`，`K=65`，`N=512`，hybrid CA-SCL list 8。
- `E>N` 时对 sub-block interleaving 后的 N bits 按 `k mod N` repetition；接收端对重复位置 LLR 求和后交给 Sionna rate recovery/decoder。
- 单位范数 DFT8 precoder 在 8 个物理 bundle 上依次轮询 0–7；同一 bundle 的 data/DMRS 使用同一向量。
- static Sionna 1.0.2 TDL-A，DS 100 ns、3.5 GHz、0 km/h；bundle 内联合时域 LS+频域线性插值、MRC。
- SNR 为单位总发射功率、单位能量 active PDCCH RE 的 Es/N0；seed `20260903`。

## 3. 执行与判据

1. 扩展 `SionnaPDCCHPolarCodec`，不修改安装目录中的 Sionna；配置入口开放 AL8，同时继续拒绝未验证的 AL16。
2. 测试 AL8 RE/E/N 计数、标准 repetition 序列、重复 LLR 合并、无噪声 DCI CRC/Polar 解码和高 SNR TDL 小批量链路；回归 AL4。
3. smoke/预扫定位 1% 区域，随后冻结正式 SNR 网格。正式每点至少 2,000 trials、目标 100 errors、最多 20,000 trials。
4. 正式输出写入 `outputs/experiment030_pdcch_al8/20260903_main/`，包含展开配置、CSV、逐 trial flags/CE、元数据、日志和曲线。
5. 生成 `research/result-030-PDCCH-AL8-Polar.md` 与 `research/result-030-PDCCH-AL8-Polar-text.md`；更新 `docs/PDCCH_SIMULATION.md`，保持简短。

验收通过条件：测试全部通过；AL8 无噪声可完整解码；正式 BLER 随 SNR 总体单调，并至少得到一个 BLER 不高于 0.01 且具有非零错误计数的实测点。结果未经研究者确认前不更新 `KNOWLEDGE.md` 或 `GOALS.md`。
