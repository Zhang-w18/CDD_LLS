# docs 目录导航

## 当前应优先使用的文档

根目录文档是项目当前规范来源：

- `../AGENTS.md`：Agent 工作流、文档职责和实验边界；
- `../GOALS.md`：研究问题、阶段进展和验收状态；
- `../DESIGN.md`：当前理论框架、公式、物理定义和规范符号；
- `../KNOWLEDGE.md`：已有正式结果支持的全局结论、负结果和开放问题。

`docs/` 下当前可直接使用的专题材料：

- `design/CDD_DELAY_SELECTION_RULES.md`：严格 Sidon、RMS 等差、$99\%$ 有效支撑等差和厚 Sidon 的参数化生成规则；
- `reports/`：面向研究者的阶段或专题报告；使用其结论时核对对应 `research/result-NNN-text.md`；
- `figures/`：被 Markdown 结果和报告引用的图，不作为脱离原始数据的独立证据。

## 背景和历史追溯材料

- `FINDINGS.md`：实验 1--21 的详细证据；用于追溯早期结论，不是当前设计规范；
- `V_design_progress_report_H1_H3.md`：result-023/024 阶段报告；用于理解 Track B 和严格 Sidon 起点；
- `design/DESIGN_ANNOTATED.md`：历史注释版设计文档；仅按 `DESIGN.md` 或其他当前文档的明确指针追溯；
- `archive/`：已被当前文档替代的历史材料，不应作为当前事实或配置来源。

判断某个实验结论是否有效时，优先读取相应的 `../research/result-NNN-text.md`、机器可读 CSV/JSON 和日志；不要仅依据阶段报告、图片或 archive 中的旧描述。
