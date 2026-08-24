# A. 总体评价

本次审阅以以下证据优先级进行：最新 LaTeX 源文件与实际 CSV 输出；最新代码与原始目录；导师批注与学校要求；交接说明；旧稿及早期 Mc=1.7 文件。GitHub 远端 `main`、本地 `HEAD` 与 `IC_Research_Project-main.zip` 均指向提交 `4c0741c5621a0968f31ffda2229396dc1e569a28`，关键代码和 Mc=2.0 输出的哈希一致。`IC_Research_Project.zip` 同时含早期 Mc=1.7 与后期 Mc=2.0 文件，因此不能仅凭文件名或压缩包内顺序判断最终分析。

当前研究设计是清楚且有价值的：在同一 (M_c=2.0) catalogue、同一 penalised non-stationary temporal ETAS likelihood 和同一 triggering starts 下，只改变三种 declustering 所构造的背景率初值，从而隔离 initialisation sensitivity。实际结果支持一条连贯主线：三种方法产生明显的 event-classification 和 starting-background 差异，但这些差异在 baseline ETAS 优化后大幅衰减，最终 triggering parameters 与 fitted background functions 几乎一致。

不过，上传的论文还不是可提交稿。最新 `main.tex` 已改为 inference-sensitivity 题目并含较完整的 Methods/Results，但 Abstract、Introduction、Discussion、Endmatter、作者信息仍是模板或空白，正文与补充材料仍保留大量学校模板提示。压缩包内 `main.pdf` 是 2025 年的通用模板 PDF，与 2026-08-23 的 `main.tex` 不一致，不能作为论文现状依据。当前 `main.tex` 还引用了两个不存在的结果表标签，使用了 `\toprule` 等命令却未加载 `booktabs`，而 supplementary 又引用了压缩包内不存在的图像。

本次已经完成并交付：

- 可直接替换的 English LaTeX Method、Results、Discussion；
- 逐项问题分级、评分和修改记录；
- 对 raw catalogue、Mc、declustering、ETAS、propagation 与 robustness 输出的重新核对；
- 同一日网格上的 initial-to-final propagation 修正分析；
- NN mixture threshold 的可复现持久化输出。

修订后的核心结论应表述为：在本 catalogue、三种已实现的 declustering-derived starting functions 与 baseline penalised non-stationary temporal ETAS specification 下，fitted inference 对 initialisation 高度不敏感；这不等于 declustering classification 不重要，也不能外推为 forecasting equivalence、global optimiser robustness 或所有 ETAS 使用方式的一般结论。

# B. 审阅标准及修改前后评分

以下量表根据应用统计硕士论文要求为本项目定制。它是编辑性评估，不是学校正式成绩。“修改后”评分假设本次修订章节被完整合入论文；尚未完成的 Abstract、Introduction、Endmatter、全稿编译和代码复现仍会限制最终成绩。

| 维度 | 权重 | 修改前 | 修改后 | 评价依据 |
|---|---:|---:|---:|---|
| 研究问题、研究设计与贡献一致性 | 15 | 10 | 13 | 最新题目和主实验已一致；修订后明确 controlled initialisation design、适用范围与贡献 |
| 统计方法正确性、完整性与可复现性 | 20 | 10 | 14 | 修正数据范围、Mc 解释、NN threshold、ETAS bounds、likelihood/integration；代码 clean-run 问题仍待修 |
| Results 与 Method/实际输出对应 | 20 | 11 | 15 | 补齐表格、2019 class imbalance、ETAS/robustness 数值；传播指标改为同网格 |
| 统计解释、效应量与不确定性 | 15 | 7 | 11 | 移除过度“保守/完整/等价”语言；明确描述性指标；仍无 CI、bootstrap 或 equivalence margin |
| Discussion、文献联系与结论克制性 | 15 | 0 | 10 | 原稿无 Discussion；现已提供完整 Discussion，并补入 ETAS 优化与初值敏感性的直接文献 |
| 学术写作、结构、图表、引用与呈现 | 15 | 6 | 9 | 修订章节结构和术语已统一；全稿仍含模板、缺 Introduction/Abstract、存在构建问题 |
| **总分** | **100** | **44** | **72** | 修改后达到“有条件的 Distinction 边缘”核心章节水平；整篇论文尚未达到 submission-ready |

若未修复 Critical 项，当前完整论文仍可能因不完整、不可复现和构建不一致被压到 Pass/Fail 边界；若本次章节被正确整合并补齐 Introduction、Abstract、Endmatter、复现说明和必要诊断，核心研究内容具备 70+ 的潜力。更高分数取决于对 original contribution、现有工作联系、统计诊断和独立研究判断的进一步展示。

# C. 按 Critical、Major、Moderate、Minor 分类的问题清单

## Critical

1. **论文源文件与 PDF 不一致。** 最新 `main.tex` 的题目为 inference sensitivity，但 zip 内 `main.pdf` 是 2025 年通用模板，含占位姓名、旧日期和空章节。后续必须由最新源文件重新编译，不能提交现有 PDF。
2. **完整论文尚未完成。** Abstract、Introduction、Discussion、Endmatter、作者/CID/导师/声明仍为模板；supplementary 也保留教学示例。当前状态不可提交。
3. **原始数据空间范围写错。** 原稿写 (30^circ)–(39^circ) N、(124^circ)–(111^circ) W；实际 `SearchResults.txt` 的观测范围为 (34.500^circ)–(36.500^circ) N、(118.99983^circ)–(116.00400^circ) W。查询边界没有保存在文件中，不能把观测极值当作查询条件。
4. **(M_c=2.0) 的解释过强。** 2019 的 MBS-WW 为 3.3，因此 (M_c=2.0) 不能写成“保证全时期完整”或“conservative threshold”。应写为统一 catalogue 的 operational compromise，并在 limitations 中讨论 2019 residual incompleteness。
5. **原 propagation ratios 比较了不同时间网格。** 原 CSV 用月度 raw rates 作分母、日网格 final spline curves 作分子。修订后在共同日网格比较实际 spline starts，MAE ratios 为 0.003835、0.008216、0.011744，对应 98.83%–99.62% attenuation。核心结论不变，但原“所有均小于 1%”不再完全成立。
6. **最新 `project_code.R` 不能从干净 R session 顺序运行。** rolling-Mc 段无条件调用未定义的 `estimate_maxc()`；该函数只存在于旧 `projectcode.R`。因此现有 CSV 可核对，但代码尚不满足学校的端到端 reproducibility 要求。

## Major

1. 主代码是 monolithic script，执行时安装包、重新运行昂贵模型并覆盖输出；没有 package lockfile、session information、stage cache 或清晰 reproduction command。
2. optional declustering sensitivity 中调用未定义的 `reasenberg_style_declustering()`，而实际函数名为 `reasenberg_style_declustering_fast()`；虽因 `RUN_SENSITIVITY <- FALSE` 未触发，仍属于潜在失败路径。
3. final ETAS 只有三种 scientifically structured starts。`RUN_EXTRA_STARTS <- FALSE`，且没有 `ETAS_extra_starting_value_robustness_Mc_2_0.csv`；因此不能声称已证明 global optimum 或 generic optimiser robustness。
4. 没有 time-rescaling/residual diagnostics、influence diagnostics 或其他 absolute model-fit 检查。AIC/BIC 只能支持相对模型选择，不能证明 final ETAS adequacy。
5. 没有 parameter standard errors、confidence intervals、bootstrap 或 formal equivalence margins。结果可写 numerical similarity/low sensitivity，不能写 statistically insignificant 或 statistically equivalent。
6. robustness 只固定 NN start，且只改变 spline df 和 cutoff；没有 penalty sensitivity，也没有三种 initialisation × alternative specifications 的完整组合。
7. 原 Results 引用了 `tab:initial_background_metrics` 与 `tab:etas_parameters`，但源文件中没有对应表；robustness 也缺精确参数变化和 background-curve metrics。
8. GFT 代码正式保存的是 (M_c=1.7) 诊断，而 (M_c=2.0) 是综合时序异质性后的人为共同阈值。Methods/Results 必须区分 estimator、diagnostic 与 operational choice。
9. NN mixture threshold 只打印到 console，未保存；本次重建并保存为 (-2.739873)。正式复现流程应将 threshold、component means、seed 与 model form 一并输出。
10. current supplementary 引用多张 zip 内不存在的图，并保留学校示例章节；必须重建或删除，否则论文不能完整编译。

## Moderate

1. 原稿把 167,584 个 pre-(M_c) records 与 15,207 个 main-analysis events 的描述混在一起；2019 的 50,276 raw events（30.0%）与 6,865 filtered events（45.1%）也需区分。
2. raw monthly background rates 与实际 spline starting functions 不是同一对象；修订后分别报告，并以实际 spline starts 进行传播比较。
3. Reasenberg 实现必须始终称为 Reasenberg-style；NN 也不是 Zaliapin--Ben-Zion 2020 stochastic thinning algorithm 的精确复现。
4. 2019 raw agreement 高主要由 clustered 类占比极高造成；不能写“2019 overall agreement lower”。
5. background-rate MAE/RMSE 单位必须统一为 events/day；correlation、Jaccard、(kappa) 和 ratios 无单位。
6. diagnostic stationary/piecewise ETAS 与 final spline ETAS 的 parameter transform、history treatment 和 model purpose 不完全相同；只能用于 screening，不应与 final fit 作逐参数比较。
7. main ETAS 不使用空间位置，且没有 catalogue-start 前的 triggering history；这两个假设需在 limitations 中明确。
8. 现有 bibliographic keys 均存在，抽查的 DOI 与标题一致；但关于 ETAS likelihood starting values/local optima 的直接支持仍缺。

## Minor

1. 全文需统一 British spelling：catalogue、modelling、initialisation；正式题目可保留已确定的 “Initialization”。
2. 删除 “this study first…”, “the method of…”, 重复的 NN–Reasenberg 句子和模板化表达。
3. `\toprule`、`\midrule`、`\bottomrule` 需要在 preamble 加 `\usepackage{booktabs}`。
4. supplementary 中出现两次 “Supplementary Materials” 并保留统计模板示例，应清理。
5. figure captions 需说明 EDA 图是 pre-(M_c) catalogue，并确认所有图在正文首次引用顺序正确。

# D. 修改记录表

| 位置 | 修改前 | 修改后 | 修改原因 |
|---|---|---|---|
| 论文定位 | inference + forecasting 的历史残留 | 仅研究 fitted non-stationary ETAS inference 对 declustering-derived initialisation 的敏感性 | 与最新题目、实际分析和无 forecast outputs 的事实一致 |
| Study design | 三种 declustering 方法介绍为三个分析 catalogue | 明确同一完整 (M_c=2.0) catalogue，declustering 只提供 background starts | 隔离 initialisation effect，避免混淆数据过滤效应 |
| Data range | 30–39 N、124–111 W | 报告实际观测范围 34.500–36.500 N、118.99983–116.00400 W，并把 query bounds 标为待确认 | 原稿与 `SearchResults.txt` 冲突 |
| Missing/duplicates | 泛称删除 missing，未说明实际数量 | 明确 required-field filtering 后仍为 167,584；event IDs 唯一；同 timestamp 不自动删除 | 提高可复现性，避免虚构去重步骤 |
| EDA | 未区分 pre-(M_c) 与 final catalogue | 明确 EDA 在 threshold selection 前；分别报告 raw 与 filtered 2019 counts | 修复研究流程和样本量混淆 |
| (M_c) 方法 | MBS-WW、GFT、阈值选择关系模糊 | MBS-WW 为 estimator，GFT 为 diagnostic，(M_c=2.0) 为统一 operational choice | 与代码和输出一致 |
| (M_c) 解释 | “conservative” 并暗示全期 complete | 明确低于 2019 estimate 3.3，并保留 residual incompleteness limitation | 防止过度陈述 |
| NN 方法 | threshold 与 catalogue (b)-value 未报告 | 增加 (hat b=0.8152)、(d_f=1.6)、V-mixture means 与 threshold (-2.740) | 使方法可复现 |
| Declustering terminology | “Reasenberg”/“nearest-neighbour declustering” 容易暗示 canonical implementation | 改为 Reasenberg-style 和 deterministic NN cluster classification，并说明与原算法差异 | 避免错误归因 |
| ETAS likelihood | 只给 intensity，缺 integration、boundary、bounds | 增加 penalised likelihood、daily trapezoid background integral、analytic truncated trigger integral、prehistory limitation、参数边界 | 满足统计硕士 Methods 可复现性 |
| Initial backgrounds | raw monthly series 与 spline starts 混称 | 分开报告 raw monthly rates 与实际 smooth spline starts | 保证 Methods–Results 对应 |
| Propagation | 月度分母/日度分子 ratios 0.00372–0.00980 | 同日网格 ratios 0.00384–0.01174，attenuation 98.83%–99.62% | 修正比较口径 |
| 2019 agreement | 仅说 raw agreement 很高或“overall agreement lower” | 同时报 background fractions、raw agreement、Jaccard 和 (kappa) | 处理严重类别不平衡 |
| ETAS Results | 只用文字说参数相似，引用缺失表 | 增加 log-likelihood、(K,alpha,c,p)、convergence 和 penalised-objective 范围 | 提供可核对数值 |
| Robustness | 仅报 log-likelihood，称 broadly robust | 增加 parameter percentage changes、background MAE/RMSE/correlation，并限定固定 NN start | 防止超出分析范围 |
| Interpretation | “no effect”, “equivalent”, “significant” 风险 | 改为 low numerical sensitivity、almost completely attenuated、within the starts/specification examined | 无 formal inference，不得使用显著性/等价性措辞 |
| Discussion | 空白模板 | 增加主要发现、解释、文献比较、优势、局限、实践意义和克制结论 | 完成用户要求和学校结构要求 |

# E. 修订后的 Method

完整、可直接替换的 English LaTeX Method 位于：

[`revised_method_results_discussion.tex`](./revised_method_results_discussion.tex)

该文件从 `\section{Methods}` 开始，包含 study design、data/preprocessing、EDA、magnitude completeness、model screening、three declustering methods、output harmonisation、construction of starts、penalised ETAS likelihood、optimisation、same-grid propagation metrics 和 limited robustness。加入主论文前，应在 preamble 添加 `\usepackage{booktabs}`。

# F. 修订后的 Results

同一文件中的 `\section{Results}` 是修订后的完整 Results：

[`revised_method_results_discussion.tex`](./revised_method_results_discussion.tex)

Results 已加入所有关键 headline numbers 和可直接使用的 LaTeX tables，包括 model comparison、declustering fractions、pairwise agreement、actual spline starting functions、three ETAS fits、same-grid propagation 和 specification robustness。新计算的可复现 CSV 为：

- [`verified_daily_grid_initial_background_metrics_Mc_2_0.csv`](./verified_daily_grid_initial_background_metrics_Mc_2_0.csv)
- [`verified_daily_grid_propagation_Mc_2_0.csv`](./verified_daily_grid_propagation_Mc_2_0.csv)
- [`verified_nn_mixture_diagnostic_Mc_2_0.csv`](./verified_nn_mixture_diagnostic_Mc_2_0.csv)

对应的轻量复现脚本为 [`compute_verified_initialisation_metrics.R`](./compute_verified_initialisation_metrics.R)。该脚本只读取已有 catalogue/background/final-fit CSV，不重新拟合昂贵的 ETAS 模型。

# G. 完整 Discussion

同一文件中的 `\section{Discussion}` 是完整、无小标题的 English LaTeX Discussion：

[`revised_method_results_discussion.tex`](./revised_method_results_discussion.tex)

Discussion 只解释 Results 中已出现的发现；使用现有且已核对的 `perry2024comparative`、`mizrahi2021effect`、`kattamanchi2017nonstationary` 和 `zaliapin2020declustering` 等 citation keys；没有新增虚构文献或 forecast 结果。关于 ETAS likelihood 的数值不稳定、平坦或多峰目标函数及起始值依赖，已补入 `veen2008estimation` 与 `lombardi2015simulated`，BibTeX 条目位于 [`additional_references.bib`](./additional_references.bib)。

# H. 仍需我确认或补充的信息

1. `[需要作者确认：SCEDC 查询时设置的 latitude、longitude、depth、magnitude 上下限，以及实际下载日期。]` 当前只能可靠报告文件内实际观测范围。
2. `[需要作者确认：学校最终 page limit。]` 2026-08 marking criteria 与 writing-guide 第 14 页写 35 pages，但当前 LaTeX template 与 writing-guide 第 16 页写 30 pages。按“最新原始文件优先”应采用 35，但提交前应向课程方确认；若无法确认，按更严格的 30 pages 排版最安全。
3. `[需要作者确认：作者姓名、CID、supervisor、submission date、declaration date、acknowledgements。]`
4. `[需要作者补充：Abstract、Introduction 与 Endmatter。]` 本次任务按要求重点修改 Method/Results 并撰写 Discussion；完整提交还需要 research gap、contribution、data/code availability 和 AI-use statement。
5. `[需要作者确认：是否在最终提交前执行 neutral/jittered multiple starts。]` 目前代码开关为 FALSE，不能把局部最优稳健性写成已证实结果。
6. `[需要作者确认：是否增加 time-rescaling residual diagnostics。]` 这是当前最重要的未完成 model-assumption check。
7. `[需要作者确认：是否增加 parameter uncertainty 或预先定义的 equivalence margins。]` 若时间不足，至少保留本次 Discussion 中“描述性、非正式等价检验”的限制。
8. `[需要作者处理：将 additional_references.bib 中两条已核验文献合并进主论文 refs.bib。]`
9. `[需要作者处理：将 analysis repository 中使用的 figures 复制到 LaTeX project，并清理 supplementary 模板和失效图片引用。]`
10. `[需要作者处理：修复 project_code.R 中未定义函数、移除运行时安装包、记录 R/package versions，并给出一次干净环境的 reproduction command。]`

# I. 最终一致性检查结果

| 检查项 | 结果 | 说明 |
|---|---|---|
| 研究问题与设计 | 通过 | 统一为 initialisation sensitivity，不再包含 forecasting |
| 证据版本 | 通过 | GitHub remote/local/最新 main zip 均为 commit `4c0741c`; LaTeX source 比 zip 内 stale PDF 新 |
| Raw sample | 通过 | 167,584 events；2019 raw (n=50,276) |
| Main sample | 通过 | (M_c=2.0), (n=15,207)；2019 (n=6,865) |
| 数据范围 | 已修正 | 采用 raw file 实际观测范围，不再使用错误的 30–39/124–111 |
| Mc 结果 | 通过但有限制 | overall 1.7、2019 3.3、non-2019 1.9；2.0 明确为折中阈值 |
| Declustering counts | 通过 | GK 2,599；NN 4,573；R 3,847 background-like |
| Agreement metrics | 通过 | overall 和 2019/non-2019 均逐项对应 CSV；类别不平衡解释已修正 |
| Initial background metrics | 通过 | raw monthly metrics 与 actual spline-start metrics 已分开 |
| ETAS parameters/logLik | 通过 | 三 fit convergence code 0；数值逐项对应 CSV |
| Propagation ratios | 已修正并通过 | initial/final 现均在共同 daily grid；核心结论不变 |
| Robustness | 通过但范围有限 | 精确数值对应三份 CSV；明确只固定 NN start |
| Causal/显著性语言 | 通过 | 未使用因果、不当 significance 或未经检验的 equivalence |
| 引用 | 通过 | 原有 citation keys 均存在；抽查 DOI/题名一致；ETAS 优化与初值敏感性文献已补入 |
| LaTeX fragment syntax | 通过 | 用 minimal article + amsmath/natbib/booktabs/graphicx 完成 draft-mode syntax check |
| 完整论文编译 | 未通过 | 当前原始 main.tex/supplementary 有模板、missing assets 和缺失 `booktabs`；需整合后再编译 |
| 原始文件保护 | 通过 | 未修改或删除原论文 zip、PDF、raw data、已有代码/output；保留用户已有 PNG 工作树改动 |

最终判断：修订后的 Method、Results 和 Discussion 在研究逻辑、数值、术语和适用范围上已一致；尚未解决的问题已集中在 H，不再隐含猜测。当前最优先的下一步不是继续润色，而是把修订章节合入最新 `main.tex`、修复 clean-run reproducibility、补齐 Introduction/Abstract/Endmatter，并完成一次全稿编译与视觉检查。
