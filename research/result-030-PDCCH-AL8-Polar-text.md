# result-030-text：PDCCH AL8 Polar repetition 与 BLER（无图版）

> 对应 `research/plan-030-PDCCH-AL8-Polar.md`。实现和正式运行已完成，结论待研究者确认。

## 1. 结论

事实：PDCCH 分支现支持 A=41、AL8、E=864。实现不修改 Sionna 安装文件：先用 Sionna 生成 N=512 的 DCI Polar 序列，再按 38.212 circular-buffer 顺序 repetition 至 E=864；接收端先对重复位置 LLR 求和，再进入 Sionna rate recovery 与 hybrid CA-SCL-8。

事实：正式 BLER 从 -4 dB 的 0.2445 单调降至 -0.5 dB 的 0.00465。-1 dB 为 101/8800=0.011477；-0.5 dB 为 93/20000=0.00465，Wilson 95% 区间 `[0.003798,0.005693]` 完全低于 0.01，验收通过。

## 2. 配置与实现口径

- 48 PRB、2 CORESET symbols、30 kHz、FFT 4096、CP 288；8 Tx/1 Rx；static Sionna 1.0.2 TDL-A，DS 100 ns、3.5 GHz、0 km/h。
- non-interleaved mapping、first CCE 0、AL8、REG bundle L=6：8 CCE=48 REG=8 bundles，432 data RE、144 DMRS RE，QPSK 得 E=864。
- A=41、CRC24C、RNTI `0xFFFF`，K=65、N=512。8 个 bundle 依次使用单位范数 DFT8 向量 0–7；data/DMRS 共用 bundle 向量。
- estimated CSI 为 bundle 内跨两个 symbol 联合 LS 后频域线性插值；MRC。SNR 为单位总发射功率、单位能量 active PDCCH RE 的 Es/N0。
- seed `20260903`；每点至少 2000 trials、目标 100 errors、上限 20000。

展开配置：`outputs/experiment030_pdcch_al8/20260903_main/resolved_config.yaml`。正式配置 SHA-256：`68B55862D7FD02B36D3FAF115FC8FD153402AE720F9155B84D797AC0ACB4291E`；基准提交 `83b25bc2ea77324ae98f7dc67136fc1e93b1f32c`，代码为未提交的 result-029/030 工作区增量。

## 3. 原始结果

| SNR dB | errors/trials | BLER | Wilson 95% | CE NMSE dB |
|---:|---:|---:|---:|---:|
| -4.0 | 489/2000 | 0.244500 | [0.226165, 0.263814] | -0.459 |
| -3.0 | 213/2000 | 0.106500 | [0.093727, 0.120782] | -1.416 |
| -2.0 | 100/2900 | 0.034483 | [0.028434, 0.041764] | -2.505 |
| -1.5 | 102/4400 | 0.023182 | [0.019134, 0.028062] | -3.044 |
| -1.0 | 101/8800 | 0.011477 | [0.009455, 0.013926] | -3.493 |
| -0.5 | 93/20000 | 0.004650 | [0.003798, 0.005693] | -3.995 |

共 40100 trials、1098 errors，用时 930.8 s。精确 CSV、元数据、日志和逐 trial flags/CE 位于 `outputs/experiment030_pdcch_al8/20260903_main/`；flags 长度和错误和均与 CSV 一致，本无图版不依赖图片即可核验。

复现：

```powershell
python tools/run_pdcch_bler_curves.py --config configs/pdcch_result030_al8_formal.yaml --stage validate
python tools/run_pdcch_bler_curves.py --config configs/pdcch_result030_al8_formal.yaml --stage run
```

## 4. 验收与限制

测试覆盖 AL8 RE/E/N、完整 DFT8 bundle 轮询、标准 repetition、LLR 合并、无噪声 DCI 解码及高 SNR TDL 链路。正式 BLER 与 CE NMSE 均无相邻反向，1% 具有上下实测点。

推断：结果支持 AL8 rate-matching 与现有链路接口正确工作，不支持将 AL4/AL8 差异全部归因于编码增益，因为两者占用的 REG 数、DMRS 数和 DFT 向量覆盖也不同。AL16、候选散列、CFO/ICI 和 PDCCH 专用 CDD 优化仍未支持。未经确认，不更新 `KNOWLEDGE.md`/`GOALS.md`。
