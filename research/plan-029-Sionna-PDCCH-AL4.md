# plan-029：Sionna PDCCH 链路分支开发与 AL4 BLER 验收（已确认）

## 1. 交付目标

在现有 PDSCH 链路之外新增可复用 PDCCH 分支，支持标准 CCE/REG/DMRS 资源映射、DCI CRC24C 与 RNTI 掩码、Sionna Polar 编译码、QPSK、REG-bundle 常量 DFT precoder cycling、Sionna TDL 信道、ideal/estimated CSI 接收和 bit-level DCI BLER。当前轮只验证平台能力，不搜索或比较 PDCCH 专用 CDD 时延。

## 2. 系统模型与公平性

- PDCCH 单层；1 CCE 含 6 REG，每个 REG 含 9 data RE 和 3 DMRS RE，故 aggregation level 4 时 data/DMRS RE 为 216/72，rate-matched bits 为 $E=432$。
- DCI payload 为 41 bits，CRC24C 为 24 bits，CRC mask RNTI 为 `0xFFFF`；调制固定 QPSK；Polar decoder 为 CRC-aided hybrid SCL，list size 8（SC 预解码，CRC 失败时回退 SCL）。
- 标准 CCE-to-REG mapping 支持 non-interleaved 与 interleaved；本轮验收使用 non-interleaved、REG bundle size 6。预编码在完整 REG bundle 内恒定，相邻 bundle 按配置顺序轮询归一化空间 DFT 向量；DMRS 与 data 使用同一向量。
- SNR 定义为单位总发射功率、单位 QPSK RE 能量下的 $E_s/N_0$；所有 DFT 向量范数为 1。噪声在 data 和 DMRS RE 上使用相同复高斯方差。
- 已知 candidate、已知 RNTI；DCI block error 定义为 payload 不一致或 CRC failure。当前轮不模拟 blind decoding false alarm。

## 3. 验收配置与预算

| 参数 | smoke | 正式验收 |
|---|---|---|
| CORESET | 48 PRB，2 symbols，30 kHz，4096 FFT，288 CP | 同 smoke |
| 天线/预编码 | 8 Tx / 1 Rx / 1 layer；6-REG bundle；8 维 DFT `[0..7]` 轮询 | 同 smoke |
| 信道 | Sionna 1.0.2 static TDL-A，DS 100 ns，3.5 GHz，0 km/h | 同 smoke |
| DCI | A=41，CRC24C，RNTI `0xFFFF`，AL4，E=432，QPSK，hybrid SCL-8 | 同 smoke |
| 接收机 | noiseless/ideal 与 bundle 内 DMRS LS+插值 | estimated CSI；另保留 ideal CSI 功能测试 |
| SNR/trials | 低预算预扫确定闭合区间 | 冻结 SNR 网格；每点至少 2,000 trials，目标 100 errors，最多 20,000 trials |
| 种子 | 20260902 | 20260902；absolute-trial 派生 |

正式网格在 smoke 后、正式 trial 1 前冻结并写入展开配置。必须有相邻实测点分别高于和不高于 1% BLER；1%附近点至少 30 errors，否则追加至上限并标记样本边界。BLER 报告 Wilson 95% 区间，不平滑原始点。

## 4. 实现与执行步骤

1. 实现 PDCCH config 校验、CCE/REG bundle mapping、candidate data/DMRS 坐标和计数。
2. 实现基于 Sionna `Polar5GEncoder/Decoder` 的 DCI codec；用确定性 affine offset 补齐标准 DCI CRC all-one 初始化与 RNTI mask；实现 PDCCH scrambling/descrambling。
3. 实现 REG-bundle 常量 DFT cycling、Sionna TDL 频域链路、MRC、ideal CSI 和 bundle-aware LS/interpolation estimated CSI。
4. 增加 YAML 入口、CSV/JSON/PNG 输出、resolved config、日志与可恢复的 deterministic trial 口径。
5. 单元测试计数、mapping、DFT 恒定性、noiseless Polar/CRC、配置拒绝和小批量链路；随后 smoke、冻结正式 SNR 网格并运行至 1%双侧闭合。
6. 生成 `research/result-029-Sionna-PDCCH-AL4.md`、`research/result-029-Sionna-PDCCH-AL4-text.md` 和不超过 1000 字的 `docs/PDCCH_SIMULATION.md`。

## 5. 修改范围、输出与复现

预计新增 `cdd_lls/phy/pdcch.py`、`cdd_lls/phy/pdcch_codec.py`、`cdd_lls/sim/pdcch.py`、`tools/run_pdcch_bler_curves.py`、PDCCH YAML 与测试；只在确有必要时小范围修改公共接口，PDSCH 默认行为必须回归通过。

正式输出固定为 `outputs/experiment029_pdcch/20260902_main/`，至少包含 `resolved_config.yaml`、`bler_points.csv`、`trial_error_flags/`、`bler_curve.png`、`run_metadata.json` 和日志。result 必须核对 A/K/E/N、REG/data/DMRS 数、bundle 映射、归一化、SNR 定义、每点 errors/trials、Wilson 区间、1%闭合、复现命令和限制。

验收通过条件：全部测试通过；指定配置无噪声可完整解码；正式 estimated-CSI 曲线无无法解释的明显反向，并以实测点双侧闭合 1% BLER。若无法闭合，则保存负结果并定位是 Polar/LLR、CE、资源映射还是预算问题，不得后验改变物理口径。
