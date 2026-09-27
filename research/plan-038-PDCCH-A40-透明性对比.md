# plan-038：PDCCH A=40 CDD 透明性与 precoder cycling 对比

状态：计划、配置与执行脚本已准备；正式仿真尚未运行。

## 1. 目的与研究问题

本轮在同一 4Tx/2Rx、2-symbol PDCCH 链路中，以 payload size $A=40$ 比较 AL1/2/4 下四种 CDD 发射波形的 transparent/non-transparent estimated-CSI BLER 与信道估计 NMSE，并保留 transparent precoder cycling 基线。需要回答：

1. 对同一发射波形，接收机使用真实等效协方差（non-transparent）与只使用物理 TDL-C 协方差（transparent）时，BLER 和 data-RE CE NMSE 相差多少；
2. Sidon、B0 QC、CDD911、CDD130 和 precoder cycling 的排序是否随 AL 改变；
3. 平均 CE NMSE 的排序是否与 BLER 排序一致。

本轮是独立的 A=40 正式实验，不复用或混合 plan-031 的 A=41 trial。结论只适用于本文冻结的信道、资源、功率和接收机口径；不得直接推广为其他 PDP、速度、payload 或 AL 下的普适结论。

## 2. 冻结系统配置

- 单层 4Tx/2Rx；两个 Rx 分支独立同分布，检测时 coherent MRC。
- 48-RB CORESET、两个相邻 OFDM symbols、30 kHz SCS、FFT 4096、CP 288。
- non-interleaved CCE-to-REG mapping、first CCE 0、6-REG bundle。
- AL1/2/4 分别占用每 symbol 3/6/12 RB，唯一频率子载波数 $K=36/72/144$；DMRS RE 数为 18/36/72，data RE 数为 54/108/216，coded bits $E=108/216/432$。
- payload $A=40$、CRC24C、RNTI `0xFFFF`、QPSK、Polar `hybSCL`、list size 8。
- Sionna TDL-C 300 ns、3 km/h、4 GHz、20 sinusoids、归一化信道。
- 每个 active RE 的总发射功率为 1；四个 Tx 上的预编码向量单位范数。每个 Rx 分支的噪声方差为 $10^{-\mathrm{SNR}_{\rm dB}/10}$，不因 `n_rx=2` 额外缩放。
- estimated-CSI 使用二维时频 LMMSE。NMSE 按每个 trial 的全部 Rx 分支与全部 data RE 聚合，并先在线性域求均值、最后转 dB。

## 3. 候选和接收机知识

每个 AL 固定 9 条曲线。CDD 的 data/DMRS 使用同一频域预编码。

| 波形 | 物理人工时延 | non-transparent | transparent |
|---|---|---|---|
| Sidon | 各 AL 的 grid coordinates `[0,1,3,7]` | `matched_effective` | `physical_fullband` |
| B0 QC | `[0,520.833,1041.667,1562.500] ns` | `matched_effective` | `physical_fullband` |
| CDD911 | `[0,0,911,911] ns` | `matched_effective` | `physical_fullband` |
| CDD130 | `[0,0,130,130] ns` | `matched_effective` | `physical_fullband` |
| PRG DFT4 | 每个 6-REG bundle 内固定 DFT4 向量 | 不运行 | `physical_prg` |

`matched_effective` 表示接收机使用由真实 $\mathbf V$ 与物理 PDP 构成的等效时频协方差；`physical_fullband` 表示接收机不知道 CDD delay/$\mathbf V$，在整个 candidate 带宽上只使用物理 TDL-C 协方差；`physical_prg` 表示接收机知道 6-REG bundle 边界，但不知道各 bundle 的 DFT4 索引，每个 bundle 独立使用物理协方差。

CDD911/CDD130 的重复时延是波形定义的一部分，必须设置 `allow_duplicate_delays: true`。precoder cycling 的 AL1/2/4 顺序分别为 `[0]`、`[0,1]`、`[0,1,2,3]`；AL1 只有一个 bundle，因此实际是固定 DFT index 0，不发生 bundle 间切换。

## 4. 公平性与随机流

- 所有候选具有完全相同的 payload、编码、DMRS、时频资源、总发射功率、信道模型、trial 预算和译码器。
- 全局 seed 为 `20260922`。AL、阶段和 candidate ID 进入固定 random-stream namespace；候选间、AL 间和 prescan/formal 间使用独立随机流。
- 因此候选差值按独立样本处理，不声明 paired-sample 方差缩减。
- A=40 与 plan-031 的 A=41 数据不得合并、续跑或用于正式统计；旧结果只用于确定足够宽的预扫边界。

## 5. 预扫、正式网格与预算

### 5.1 预扫

预扫只运行四条 non-transparent CDD 和 `PRG_DFT4_TRANSPARENT`，不运行四条 transparent CDD；后者直接使用配对 non-transparent 曲线冻结出的正式 SNR 点。每个预扫点固定 300 trials，不提前按错误数停止：

| AL | SNR 网格 / dB |
|---:|---|
| 1 | `0:1:10` |
| 2 | `-2:1:5` |
| 4 | `-6:1:2` |

预扫只用于冻结正式网格，不作为正式门限结果。本轮不扫描 10% BLER 位置，只要求看到 1% BLER crossing。Sidon、B0 QC、CDD911 和 CDD130 的 SNR 范围一律以各自 non-transparent（`matched_effective`）曲线为准：若其 1% 存在相邻预扫点双侧 bracket，正式网格只覆盖该 1 dB bracket，并以 0.25 dB 步长加密，不向两端额外扩展；若没有 bracket，则围绕 log-BLER 距离 1% 最近的预扫点取 $\pm0.5$ dB、步长 0.25 dB 的最小诊断窗口。对应 transparent 曲线直接复用同一波形 non-transparent 曲线的完整正式 SNR 网格，不运行 transparent 预扫，也不因其目标 crossing 状态扩展范围。precoder cycling 没有 non-transparent 配对，按自身预扫采用同一 1% bracket/最近点规则。自动选择及完整网格写入 `formal_grid_receipt.json`。

如果 transparent 曲线在由 non-transparent 冻结的范围内没有穿越 1% BLER，不追加 SNR 点，也不要求取得 1% bracket；只报告 raw BLER、单侧界限或 error floor。non-transparent 曲线或 precoder cycling 在预扫全范围内仍未接近目标时，其正式诊断窗口只用于量化负结果，不得外推目标门限。若研究者随后要求获得更低目标的 bracket，必须先补充 plan、向缺失方向扩展预扫，再重新冻结正式网格。

### 5.2 正式运行

每个 candidate/SNR 点至少 10,000 trials；达到 200 block errors 后允许停止；最多 50,000 trials。一个 trial 只要 CRC 失败或 payload 任一 bit 错误即计为 block error。所有点保存逐 trial error flags 和逐 trial CE NMSE。

正式 1% 目标只在相邻原始正式点形成不宽于 0.25 dB 的双侧 bracket 时报告；目标插值使用按 trials 加权的单调递减 isotonic BLER 和局部 $\log_{10}(\mathrm{BLER})$ 插值，95% 区间使用至少 4,000 次独立逐点 Bernoulli bootstrap。无合格 bracket 的目标只报告 raw bracket、单侧界限或未达到状态，不外推；transparent 曲线未穿越 1% 不构成实验未完成。10% BLER 不安排专用正式点、不要求 bracket，也不作为验收指标；若现有点偶然形成合格 bracket，只能作为附带结果。

## 6. 实现与执行步骤

1. 使用 `configs/pdcch_result038_al{1,2,4}_prescan.yaml` 定义三组配置；运行前仍对全部 9 个 candidate 执行公共 PDCCH runner 的完整 validate。
2. 只对四条 non-transparent CDD 和 precoder cycling 运行预扫，并保存展开配置、日志、BLER CSV、error flags 和 CE NMSE。
3. `tools/prepare_plan038_formal.py` 只读取保存的预扫 CSV，按第 5.1 节规则生成三份 formal YAML 和选择回执。
4. 再次 validate 自动生成的 formal YAML，随后通过 `tools/run_plan038.py` 最多并行运行三个互不重叠 candidate shard；相同命令可安全续跑。
5. `tools/analyze_plan038.py` 审计 27 个 candidate/AL shard，合并 prescan 与 formal 点（相同 SNR 优先 formal），生成统一绘图输入和来源哈希。
6. 正式完成后按 `docs/agent/RESULT_SPEC.md` 生成 `research/result-038-PDCCH-A40-透明性对比.md` 与对应 `-text.md`。

本轮不修改 `cdd_lls/` 链路算法；若 validate 暴露公共实现缺口，应暂停正式运行，先做最小修复和定向回归测试。

## 7. 图、数据与输出

根目录为 `outputs/experiment038_pdcch_a40/20260922_4tx_2rx_2sym/`：

- `prescan/al{1,2,4}/<candidate>/`：预扫原始数据；
- `formal/al{1,2,4}/<candidate>/`：正式原始数据；
- `formal_grid_receipt.json`：冻结网格依据；
- `analysis/plot_points.csv`、`source_receipt.json`：合图输入和来源审计；
- `analysis/plan038_all_al_bler_nmse.png`：唯一总图。

总图固定为一张 2×3 figure：三列依次为 AL1/2/4，上排为 logarithmic BLER，下排为 data-RE CE NMSE；每个 panel 含该 AL 的全部 9 条曲线。颜色区分发射波形，实线/点线区分 non-transparent/transparent，precoder cycling 单独使用虚线。图中保留原始点和 1% BLER 参考线，不用平滑曲线遮蔽非单调或 error floor。

## 8. 验收条件

- 27 个 candidate/AL 组合均通过配置、资源、A=40 Polar 无噪声编解码、2Rx、2-symbol 非零 Doppler、功率和候选几何校验。
- Sidon/B0/CDD911/CDD130 的 matched 与 transparent 配置除 `candidate_id`、label 和 `receiver_covariance_mode` 外具有相同发射波形；对应 delay seconds 在浮点容差内一致。
- 所有正式点满足 trial/error 停止规则；CSV errors/trials 与逐 trial flags 完全一致，CE NMSE 在线性域复算一致且无 NaN/Inf。
- 总图恰有 3 个 BLER panel、3 个 NMSE panel，每个 panel 9 条白名单曲线；来源回执指向本轮 A=40 数据，不读取 plan-031 A=41 CSV。
- result 分 AL 报告各曲线的 1% 目标及 95% 区间或不合格原因；同时报告 matched-transparent 差值、代表 SNR 的 NMSE、zero-noise CE floor，以及 BLER 与 NMSE 排序不一致的事实。10% 只在现有正式点自然形成合格 bracket 时作为附带结果，不为其补点。
- 结果经研究者确认前不更新 `KNOWLEDGE.md`/`GOALS.md`，不创建 Git checkpoint。

## 9. 直接运行命令

建议分阶段执行，以便预扫后审阅冻结网格：

```powershell
python tools/run_plan038.py --stage validate --max-workers 3
python tools/run_plan038.py --stage prescan --max-workers 3
python tools/run_plan038.py --stage prepare
python tools/run_plan038.py --stage validate-formal --max-workers 3
python tools/run_plan038.py --stage formal --max-workers 3
python tools/run_plan038.py --stage analyze
```

确认无需在预扫与正式运行之间人工审阅时，可一次执行：

```powershell
python tools/run_plan038.py --stage all --max-workers 3
```

`--stage all` 可能运行很久；中断后重复相同命令会跳过已经完成的点并续跑。正式运行前不得改变已经生成的 formal YAML 或 `formal_grid_receipt.json`；需要改变网格时应删除尚未产生正式 trial 的 formal 配置、重新执行 prepare，并记录原因。

## 10. 2026-09-23 均匀采样补点与最终绘图窗口

研究者检查首轮正式图后，要求在以下窗口内按 0.25 dB 等间隔补齐正式点，并使最终 BLER/NMSE 图只显示这些点。补点复用原 formal 的 seed、random-stream namespace、candidate ID、停止规则和输出目录；已有 SNR 点由 resume 跳过，不重跑、不覆盖。AL1 Sidon transparent 不新增 trial，保留已有 4--5 dB 的均匀正式点，以维持每个 AL 九条方案均在总图中。

| AL | candidate | 最终绘图/补点范围 / dB |
|---:|---|---:|
| 1 | Sidon non-transparent | 4--6 |
| 1 | Sidon transparent | 4--5；不补点 |
| 1 | B0 QC non-transparent / transparent | 5--8 |
| 1 | CDD911 non-transparent / transparent | 6--8 |
| 1 | CDD130 non-transparent / transparent | 7--9 |
| 1 | precoder cycling transparent | 7--10 |
| 2 | Sidon non-transparent | 0--1 |
| 2 | Sidon transparent | 4--5 |
| 2 | B0 QC non-transparent | 0--1.5 |
| 2 | B0 QC transparent | 1.5--2.5 |
| 2 | CDD911 non-transparent | 0.5--2.5 |
| 2 | CDD911 transparent | 2--3 |
| 2 | CDD130 non-transparent / transparent | 2.5--4 |
| 2 | precoder cycling transparent | 1.5--3 |
| 4 | Sidon non-transparent / B0 QC non-transparent | -4-- -2.5 |
| 4 | Sidon transparent | -2-- -1 |
| 4 | B0 QC transparent | -1--0 |
| 4 | CDD911 non-transparent | -3.5-- -2 |
| 4 | CDD911 transparent | -0.5--0.5 |
| 4 | CDD130 non-transparent / transparent | -2.5-- -1 |
| 4 | precoder cycling transparent | -3.5-- -2 |

`tools/run_plan038_refinement.py` 从三份原 formal YAML 生成 refinement 配置，严格校验窗口端点位于 0.25 dB 栅格，并把新增点追加到原 formal candidate 目录。`tools/analyze_plan038.py` 不再读取 prescan 点；它只读取 formal CSV，按上表逐曲线裁剪，并要求窗口内点集精确等于完整 0.25 dB 栅格，缺点或多点均停止绘图。总图仍为一张 2×3 figure，上排 BLER、下排 NMSE。
