# magnus-skills

AI 自主蒸馏科学论文的复现代码仓库。

每个分支对应一篇论文（以 `<首作者姓氏小写><年份>` 命名），包含从论文中提取的数值计算复现脚本。这些代码由 Magnus 平台的知识蒸馏蓝图（`distill-knowledge`）驱动 Claude Code 自主生成，仅依赖 `numpy`、`scipy`、`matplotlib`。

配套的蓝图和技能注册在 Magnus 平台上，可通过 `magnus run <blueprint-id>` 直接运行。
