# Method、Results 与 Discussion 全面审阅及修改说明

## 1. 审阅对象、证据优先级与版本核对

本次审阅以 `Statistics_Research_Project_Title__1_ (2).zip` 中 2026-08-24 保存的 `main.tex` 为最新论文版本。该压缩包中的 `main.pdf` 仍是 2025 年生成的学校模板 PDF，与 740 行的最新 `main.tex` 不一致，因此不能用于判断正文是否已更新。这里的“最新源文件”专指该 `main.tex`，并不指压缩包中时间较旧的 PDF。

证据优先级为：

1. 原始 catalogue、最新 Mc=2.0 CSV、ETAS 输出和可重算结果；
2. `(2).zip` 中的最新 `main.tex` 与 `refs.bib`；
3. 学校模板、写作说明和 marking criteria；
4. 导师第一版批注；
5. `论文交接说明.md` 和较早压缩包。

已确认的版本差异与处理如下：

| 项目 | 文件中的冲突 | 本次采用的处理 |
|---|---|---|
| SCEDC 空间、震级和深度范围 | `main.tex` 把 30--39°N、124--111°W、−1--9、−5--30 km 写成 catalogue 的实际范围；原始文件中的观测极值明显更窄 | 作者已确认前一组是 SCEDC 页面检索条件。Method 应同时区分“检索条件”和“下载文件中的实际观测范围” |
| 主论文与 PDF | `main.tex` 已有完整 Method/Results 和 Discussion；`main.pdf` 仍为旧模板 | 只以最新 `main.tex` 为论文内容依据；合并修改后重新编译 PDF |
| propagation ratios | 交接说明及旧结果把月度起始率与日度最终曲线混作分母、分子 | 采用共同日网格上的实际平滑起始函数和最终曲线：MAE ratios 0.003835、0.008216、0.011744 |
| page limit | 模板及写作指南部分页面写 30 页；较新的 marking criteria/说明出现 35 页 | 未确认前按更严格的 30 页控制正文最安全 |

## 2. 针对本论文的审阅标准

本量表不是学校正式评分，而是依据应用统计硕士论文要求建立的编辑性检查标准。

| 维度 | 权重 | 当前 `(2).zip` | 完成本文件所列修改后 | 主要判断依据 |
|---|---:|---:|---:|---|
| 研究问题、设计与贡献的一致性 | 15 | 11 | 14 | 受控 initialisation experiment 已成形；需在 Introduction 和 Discussion 中更直接限定贡献 |
| Method 的统计正确性、清晰度与可复现性 | 20 | 13 | 16 | 流程较完整；ETAS 解释、history cutoff 公式、query/observed range 和数值归属仍需修正 |
| Results 与 Method、代码输出的一一对应 | 20 | 11 | 17 | 主数值大多正确；raw monthly 与 smooth starts 混用是当前最严重的不一致 |
| 效应量、不确定性与克制解释 | 15 | 8 | 12 | 已使用 MAE/RMSE/Jaccard/κ；仍无 CI、bootstrap、equivalence margin 或 residual diagnostics |
| Introduction、Discussion 与文献联系 | 15 | 7 | 12 | Introduction 仍是模板；Discussion 缺 GFT、研究空白和部分限制；本次已提供补充材料 |
| 结构、语言、图表、引用和提交完整性 | 15 | 7 | 10 | 表格顺序/标签、模板文字、图像与 stale PDF 尚需整合处理 |
| **总分** | **100** | **57** | **81** | 修改后核心统计叙事可达到较强 Distinction 水平；可复现性和正式模型诊断仍限制更高评价 |

## 3. 问题分级

### Critical

1. Abstract 和 Introduction 仍是模板，当前 PDF 也不是由最新源文件生成；现状不可提交。
2. Results 第一个 initial-background 段落报告的是 raw monthly metrics，但引用的表格报告 smooth spline starts，文字和表格不是同一对象。
3. `M_c=2.0` 被称为 “conservative” 并被写成可“ensure”全时期 completeness；这与 2019 MBS-WW estimate 3.3 冲突。
4. 主代码在全新 R session 中仍有未定义函数路径，尚不能声称端到端可复现；这不会推翻现有 CSV，但必须在提交前修复或明确给出可运行入口。

### Major

1. ETAS 公式未把 $H$ 写入求和条件，但正文称使用 3,650-day triggering-history cutoff。
2. `Non-stationary ETAS model` 对只学过统计学核心课程的读者仍解释不足：conditional intensity、self-excitation、background/triggering decomposition、log link、spline penalty 和 point-process likelihood 需要直观解释。
3. 2019 GFT 已在 Results 报告，却没有进入 Discussion；因此 completeness limitation 与主要实验之间的联系不完整。
4. Method 得到的 common-catalogue $b=0.8152$、Gaussian-mixture component means −5.120/−1.432 和 threshold −2.740 没有在 Results 中报告，影响复现。
5. “Clustered” 同时包含 foreshocks 和 aftershocks，对统计背景读者含义不清，也容易被误解为已证明的物理关系。
6. 只有三种有科学依据的起始值，没有执行 neutral/jittered multi-start；不能声称 global optimum 或对任意初值都稳健。
7. 无 time-rescaling/residual goodness-of-fit、参数置信区间、bootstrap 或 formal equivalence analysis；不能使用 “statistically equivalent/insignificant” 等措辞。

### Moderate

1. diagnostic piecewise-background ETAS 与最终 spline ETAS 容易被当成同一个 non-stationary model，应明确前者只用于模型选择动机。
2. 多个表格在 `\caption` 前放置 `\label`，与学校模板要求相反；declustering summary 还有重复 caption。
3. “convergence code 0 indicates convergence” 过强；它只表示 optimiser 报告 successful termination，不证明全局最优。
4. Discussion 重复 Method 的方法细节和研究设计，但没有始终说明这些重复内容用于支持哪一项结论。
5. robustness 只固定 NN start，因此只能回答 specification sensitivity，不能证明 alternative specifications 下的 initialisation robustness。
6. temporal ETAS likelihood 没有使用位置和深度，且 catalogue 开始前的 triggering history 不可见；应进入限制。

### Minor

1. 统一 British spelling：catalogue、modelling、initialisation、optimisation；题目若已正式确定可保留 `Initialization`。
2. 修复 Discussion 的单复数、逗号拼接句和 `significantly reduced` 等无检验依据的措辞。
3. 删除 main text、Endmatter 和 Supplementary 中的学校提示文字和示例。
4. 所有图表 caption 应说明数据阶段、单位、比较对象和一句主要信息。

## 4. Method：在哪里改、改什么、为什么

| 位置 | 当前问题 | 应进行的修改 | 原因 |
|---|---|---|---|
| `Experiment Overview` | 标题较口语化；第一段虽正确但对“为何只改变初值”强调仍可更早 | 改名 `Study design`；第一句直接定义 controlled initialisation experiment，并明确不是比较三个删减后的 catalogue | 使研究问题、设计和后续结果从开头一致 |
| `Earthquake catalogue and preprocessing` 首段 | 把检索条件写成实际 catalogue 范围 | 写成：SCEDC query limits were 30--39°N, 124--111°W, magnitude −1--9 and depth −5--30 km；随后另报 observed ranges 34.500--36.500°N、−118.99983--−116.00400°、−0.98--7.10、−2.3--28.9 km | 检索条件定义研究范围；观测极值描述返回数据，二者不能混用 |
| 同上 | 只引用搜索网页 | 同时引用 `scedc2013` 数据 DOI 和 `scedc_catalog_search2026` 搜索界面 | 数据来源和访问方式均可追溯 |
| `Magnitude completeness assessment` | GFT 与 MBS-WW 的角色基本正确，但原始方法文献不足 | 为 GFT/Mc 补 `wiemer2000minimum`、`woessner2005assessing`；保留 95% primary、90% secondary diagnostic 的定义 | 对课程外方法给出原始来源 |
| 同上 | 数值应放在哪里尚不一致 | Method 只说明预先采用统一 operational threshold $M_c=2.0$ 作为分析规则；catalogue-wide/2019/non-2019 estimates 和 GFT percentages 全部放 Results | 方法部分报告规则和步骤；由数据得到的估计属于结果 |
| `Assessment of triggering and background non-stationarity` 开头 | “diagnostic ETAS” 对非地震统计读者不直观 | 在模型比较前用 3--4 句解释 homogeneous Poisson、ETAS self-excitation、stationary background 和 time-varying background；明确 “diagnostic” 指 deliberately simplified comparison model，只回答是否允许背景率随时间变化会改善相对拟合 | 满足学校要求：课程外方法必须在正文解释 |
| 同上结尾 | “not a residual goodness-of-fit test” 技术性过强 | 改成普通语言：AIC/BIC comparison can tell which candidate model fits relatively better, but cannot show that the final spline model reproduces the observed event-time pattern adequately; that requires separate residual checks | 保留统计含义，同时让本科统计读者可理解 |
| `Declustering design` 和表格 | `Clustered` 含义模糊 | 二元类别统一为 `Background/mainshock-like` 与 `Foreshock/aftershock-like sequence member`；首次出现时说明是 algorithmic labels，不是物理真值或因果判定 | 避免把所有非 background 事件笼统称为 clustered |
| 三种 declustering 方法 | Discussion 与 Method 重复解释 | 完整算法、参数和步骤只保留在 Method；Discussion 只需一句概括三者使用不同 cluster definitions | 避免重复且保留解释力 |
| NN method 最后一段 | 目前已正确说不是 exact stochastic algorithm，应保留 | 保持 “deterministic nearest-neighbour cluster classification”；把 $\hat b=0.8152$、component means 和 threshold 从方法描述中的最终数值移到 Results，Method 仅说明这些量如何估计 | 既防止错误归因，也遵守 Methods/Results 分工 |
| `Comparison ... and construction ...` 第一段 | background count 被重复描述两次 | 删除重复的 “For each method, the number and proportion ...” 句，只保留一次；随后依次定义 raw agreement、Jaccard、κ 和 three-method disagreement | 减少冗余 |
| 同一 subsection 的 log-rate 段落 | 最新版解释已基本充分 | 保留四个关键解释：log 与 final log-link 同尺度并确保正值；0.1 只为避免 log(0)；$E_m$ 是 catalogue 第一/最后保留事件之间在该月实际观测的天数；spline regression 把每月一个率转换为 optimiser 需要的 8 个系数 | 这些内容是复现 initialization 的必要信息，不应删 |
| 同一 subsection 最后一段 | 用户曾询问是否有必要保留 | **保留，但缩短。** 必须明确 declustering 只产生 alternative starts，foreshocks/aftershocks 没有从 ETAS catalogue 删除，三次拟合使用同一 15,207 events 和同一 likelihood | 这是整个 controlled design 能够隔离 initialisation effect 的关键条件 |
| `Non-stationary ETAS model` 开头 | 直接进入公式；background 的统计含义和 non-stationary 区别不够清楚 | 公式前增加直观说明：$\lambda(t)$ 是给定过去事件后单位时间内的即时预期事件率；$\mu(t)$ 是模型未归因于可见早期事件的部分；求和项是既往事件的 self-exciting contribution；non-stationary 仅指 $\mu(t)$ 允许随时间变化 | 让统计读者先理解模型再读符号 |
| ETAS intensity 公式 | 求和写 `\sum_{j<i}`，未体现 $H$ | 改为 `\sum_{j<i:\,0<t_i-t_j\le H}` | 与 3,650-day implementation 一致 |
| background spline 与 penalty | 只解释 positive rate，没有解释 penalty 的读者意义 | 补一句：8-df spline provides a smooth low-dimensional curve; the second-difference penalty discourages neighbouring coefficients from changing too abruptly and protects against an unnecessarily wiggly background | 解释未知方法而不重复线性模型课程内容 |
| likelihood 段 | 数学正确，但可读性不足 | 在公式后解释第一项奖励 observed event times with high predicted intensity，integral term penalises predicting too many events over the whole interval；说明 daily trapezoidal integration 和 analytic triggering integral | 提高可复现性和统计直觉 |
| optimisation | code 0 被后文当作证明收敛 | Method 说明记录 `optimiser status, gradient/diagnostic if available, objective and parameter bounds`；Results 使用 “reported successful termination” | 防止将软件状态码等同于全局收敛 |

### 对用户之前重点句子的明确处理

- `The fitted unequal-variance model had component means ...`、`threshold=-2.740`、`b=0.8152`：它们是数据/拟合产生的结果，应放 Results；Method 只写估计规则、fixed seed 和模型形式。
- `M_c=2.0`：作为后续分析的 inclusion rule 可以在 Method 中出现；“overall estimate=1.7、2019=3.3、non-2019=1.9、GFT=92.6%/90.8%”应放 Results。
- `Importantly, these declustering outputs were used only ...`：有必要保留，因为它定义了研究的对照条件；但可压缩为 2--3 句，避免重复。

## 5. Results：在哪里改、改什么、正确数值是什么

| 位置 | 当前问题 | 应进行的修改 | 正确证据/原因 |
|---|---|---|---|
| `Exploratory catalogue characteristics...` | 2019 被突出但未命名具体序列 | 将 2019 明确联系到 Ridgecrest sequence，并引用 `ross2019hierarchical`；说明它占 raw catalogue 50,276/167,584=30.0%，且包含 catalogue 最大 M7.1 事件 | 解释为什么不是任意挑选一年 |
| 同一 subsection | `single conservative threshold`、`ensuring ... completeness` 过强 | 改成 `common operational threshold`/`design compromise`；明确 $M_c=2.0$ 低于 2019 estimate 3.3，因此不保证 2019 complete | 与实际诊断一致 |
| 同一 subsection | GFT 报了数值但解释不足 | 保留：at $M_c=1.7$, GFT=92.6% in 2019 and 90.8% outside 2019；明确两者达到 secondary 90% criterion、均未达到 95% criterion | 这既支持 catalogue-wide threshold 附近有合理拟合，也显示 residual uncertainty |
| `Evidence for triggering...` | 把 piecewise diagnostic model 直接称作 `the non-stationary ETAS model`，容易与最终 spline fit 混淆 | 全段及表格统一为 `diagnostic piecewise-background ETAS`；最后一句只说它 motivated a time-varying background in the main model | screening logLik 27,673.39 与最终 spline logLik 约 27,767.5 不是同一个模型 |
| model-comparison table | `\label` 在 caption 前 | 每个 table 改成 `\caption{...}` 后立即 `\label{...}` | 学校模板明确要求 |
| declustering summary table | 重复 caption；`Clustered` 含义模糊 | 删除第二个 caption；列名改 `Foreshock/aftershock-like` | 避免 LaTeX 编号/目录异常和术语歧义 |
| classification Results | “almost all events were labelled clustered” | 改成 “almost all events were classified as foreshock/aftershock-like sequence members” | 清楚说明 dominant class 是什么 |
| classification Results 之后 | NN 的关键可复现结果缺失 | 新增一小段或 supplementary table：common-catalogue $\hat b=0.8152$；mixture means −5.120 and −1.432；posterior-equality threshold −2.740 | 数值来自 verified NN mixture output |
| `Declustering-derived starting backgrounds` 第一段 | 文字给 raw monthly metrics，却引用 smooth-start table；NN--Reasenberg 句子重复 | 拆为两个段落。第一段明确是 raw monthly rates，不引用 smooth table；第二段引出实际输入 optimiser 的 smooth starts 并引用表格 | 当前文字与表格的对象和数值不同 |
| raw monthly 段 | 需保留正确 headline values | GK--NN: MAE 0.209, RMSE 0.324, $r=0.490$；GK--R: 0.133, 0.180, 0.687；NN--R: 0.091, 0.185, 0.841。2019 GK--NN: 0.502, 0.632, −0.167 | 对应 `diagnostic_initialisation_difference_metrics*.csv` |
| smooth-start table | 表格数值正确，但正文使用了 raw values | 正文改为 GK--NN 0.203052/0.294838/$r=0.302010$；GK--R 0.129409/0.155881/0.503749；NN--R 0.076162/0.158183/0.930644 | 这些才是实际送入 optimiser 的函数差异 |
| `Sensitivity of fitted ETAS inference...` | “code 0, indicating all have converged” 过强 | 改成 “The optimiser returned status code 0 for all three fits, indicating successful termination under its stopping rule.” 后面继续用 close objectives/parameters 支持 low numerical sensitivity | 不把状态码解释成 global convergence |
| same-grid propagation | 文字和表格目前正确 | 保留 MAE ratios 0.003835、0.008216、0.011744，对应 99.62%、99.18%、98.83% attenuation；不要再写 “all less than 1% remained” | NN--R 保留 1.1744%，不能四舍五入为 <1% |
| robustness | 宽表格过小 | 将表拆成 (a) logLik + parameters 与 (b) background comparison，或使用 `tabularx/siunitx`；caption 明确 NN start fixed | 提高可读性并限制结论范围 |

### Results 中建议增加的核心图

新图 `initial_vs_final_background_revised.pdf` 直接展示三条 smooth starting functions 明显不同，而三条 final fitted functions 几乎重合。它比单独列 propagation table 更直观地支撑主结论，且完全由已有 CSV 重建，不属于新增分析。建议放在 `Sensitivity of fitted ETAS inference to initialisation` 中，在传播表之前首次引用。

建议 caption：

```latex
\caption{Declustering-derived smooth background-rate starting functions
(top) and the corresponding final fitted non-stationary ETAS background
functions (bottom), evaluated on the same daily grid. The shaded band marks
2019. The starting curves differ in level and temporal pattern, whereas the
three final curves are almost coincident.}
```

## 6. Discussion：在哪里改、改什么、为什么

| 位置 | 当前问题 | 应进行的修改 | 原因 |
|---|---|---|---|
| 开头 2--5 句 | 语法错误；scope 与 general implication 混在一个长段 | 第一段只报告主发现：labels 和 starts materially different；final parameters/background nearly identical；0.38%--1.17% initial MAE remained；结论限定在 this catalogue, three starts and baseline model | 符合学校要求并避免外推 |
| 主结论句 | “no effect” 容易被理解为 declustering 普遍无关 | 用 `highly robust to the three declustering-derived initialisations examined`；紧接一句 `This is not evidence that declustering is generally unimportant.` | 最终结果确实是对最终拟合影响很小，但范围只覆盖本实验 |
| 方法差异解释段 | 若重复三种算法完整定义会与 Method 重复 | Discussion 只用一句比较性概括：`The disagreement is consistent with the methods' different definitions of a sequence: fixed magnitude-dependent windows, data-derived nearest-neighbour links, and expanding interaction groups.` | Discussion 的作用是解释为什么会不同，不是重新教学算法步骤 |
| 2019 agreement 段 | `class-aware agreement measures` 不直观 | 直接写：overall agreement was high because the methods nearly always agreed on the very common foreshock/aftershock-like class; Jaccard was needed to ask whether they selected the same much smaller background-like set | 不用术语也能解释类别不平衡 |
| attenuation 段 | `Strong attenuation` 未先定义；`significantly reduced` 暗示检验 | 先定义 attenuation 是 initial pairwise curve difference 在 fitting 后的减少；用 `substantially/almost completely reduced`；给 98.83%--99.62% 数值 | 无显著性检验，不能使用 `significantly` |
| 同一段 | 现有解释正确但可更直接 | 明确 declustering 只给 starting coefficients；fitting 后背景和 triggering parameters 都按同一 penalised likelihood 重新估计；起始曲线并未作为约束保留 | 这就是差异衰减的直接设计原因 |
| optimisation interpretation | 可能让读者误以为证明 unique optimum | 保留 Veen & Schoenberg (2008)、Lombardi (2015)；明确三种 starts 到达 essentially same numerical solution，但 neutral/jittered starts 未运行，不能证明 global optimum | 与实际代码开关一致 |
| model-screening/GFT 段 | Discussion 完全未讨论 GFT | 新增一段：GFT at 1.7 reached 90% but not 95%，且 2019 MBS-WW=3.3；因此 $M_c=2.0$ 保持三个实验使用同一事件，但不消除 high-rate period 的 completeness uncertainty | 回应用户和学校对主要限制的要求 |
| 文献比较/研究空白 | 目前只说已有 declustering 比较，不够突出本文位置 | 衔接三组文献：Perry/Mizrahi 说明算法改变 labels/catalogue summaries；Ogata/Zhuang 说明 declustering 可用于 ETAS preliminary background；Kattamanchi/Kumazawa 说明 spline non-stationary ETAS；再说明本文隔离的是 starts 是否传到同一 final likelihood | 明确 original contribution，不夸称首次 |
| strengths 段 | 与 Method 的 propagation chain 有重复 | 保留，但每一点都明确“为什么增强结论”：same events isolate initialisation；Jaccard prevents dominant class masking；same daily grid makes attenuation denominator/numerator comparable；specification checks delimit robustness | 说明存在意义，而不是罗列做过的步骤 |
| limitations | 当前已有一部分，但不完整 | 加入：one region；2019 residual incompleteness；one parameter set/method；custom NN 和 Reasenberg-style；temporal not spatial ETAS；no pre-catalogue history；fixed penalty/bounds；limited specification grid；no arbitrary multistart；no CI/equivalence/residual diagnostics/forecast assessment | 诚实回答方法和发现对哪些限制稳健 |
| 最后一段 | 目前仍有 ETAS/inference/initialisation 等新闻读者术语 | 改为 plain-language press-office paragraph，见下方建议 | 满足学校明确要求 |

### 建议在 Discussion 中加入的 GFT 段落内容

该段不引入新结果，只解释 Results 中已有数值：

> The completeness diagnostics temper the interpretation of the sensitivity
> experiment. At the catalogue-wide MBS-WW estimate of $M_c=1.7$, the GFT was
> 92.6\% for 2019 and 90.8\% for the remaining years: both met the secondary
> 90\% criterion but neither met 95\%. Moreover, the separate 2019 MBS-WW
> estimate was 3.3. The common $M_c=2.0$ threshold therefore ensured that all
> three initialisation experiments used the same events, but it did not remove
> uncertainty about small-event detection during the busiest sequence.

### 建议的最后一段（Press Office/journalist audience）

> Three approaches to sorting Southern California earthquakes disagreed about
> which events represented the underlying background activity, especially
> during the intense 2019 sequence. Yet when those alternatives were used only
> as starting suggestions for the same statistical model, which then
> re-estimated the rates from all earthquakes, the final answers were almost
> identical. For this dataset, the initial disagreement therefore had little
> influence on the fitted model. Analysts should still test other regions and
> modelling choices before assuming the same reassuring result will always
> hold.

## 7. 回答此前逐句修改问题

1. **“最终结论是不是没什么影响？”** 是，但要精确写成：对本 catalogue、baseline model 和这三种 declustering-derived starts，最终 fitted ETAS inference 的数值影响很小；不能写成 declustering 对所有研究都没有影响。
2. **是否只能得到这个地区的结果？** 直接证据只来自 Southern California catalogue。研究可以提出一个更广的 methodological hypothesis，但不能声称已在其他地区得到验证。
3. **Discussion 是否要再次完整解释三种方法？** 不需要。Method 给完整定义；Discussion 用一句对比其不同假设来解释为何标签会不同。
4. **2019 段是否需要 `class-aware`？** 不需要。直接解释 dominant class 让 raw agreement 看起来很高，而 Jaccard 关注较小的 background-like set，更容易理解。
5. **`attenuation` 指什么？** 指三个 starting background functions 的 pairwise differences 在模型拟合后缩小；本研究 MAE difference 减少 98.83%--99.62%。
6. **strengths 段是否重复？** 有保留价值，但必须从“列步骤”改成“解释为何 controlled design、minority-class agreement 和 same-grid metrics 使主结论可识别、可比较”。
7. **initialisation construction 的解释是否过多？** log、pseudo-count、partial-month exposure 和 spline projection 都是复现所必需，不建议删除；只需避免和 Results 重复数值。
8. **comparison subsection 末尾的 controlled-design 段是否需要？** 需要，它是论文最关键的因果隔离条件之一：只比较 initialisation，不比较不同训练 catalogue。

## 8. 图表修改与可直接替换文件

已生成 PDF（首选，矢量）和 PNG（备用）版本：

- `figures/annual_event_counts_revised.pdf`：以色盲友好颜色突出 2019，并标注 50,276；caption 需说明 2000/2026 为 partial years。
- `figures/rolling_mc_revised.pdf`：移除无信息量 legend，直接标出 operational $M_c=2.0$、catalogue-wide estimate 1.7 和 2019；说明每个点是 1,000 successive events，间隔 250 events。
- `figures/initial_vs_final_background_revised.pdf`：同一日网格上对比 smooth starts 与 final fits；用颜色和线型同时区分方法，并标注最终曲线几乎重合。

对应可复现脚本为 `make_revised_figures.R`。脚本只读取现有 CSV，不重新运行 declustering 或 ETAS。

## 9. 建议的整合顺序

1. 用 `abstract_and_introduction.tex` 替换模板 Abstract 和 Introduction，并删除 `General tips`、`Tips for the Introduction`。
2. 把 `references_to_add.bib` 中八个新条目合并进 `refs.bib`；已有 citation keys 不要重复加入。
3. 按第 4--6 节逐项修改 Method、Results、Discussion，不直接覆盖原始 ZIP。
4. 替换两张 EDA 图，加入 initial-vs-final 主结论图。
5. 将所有 table 调整为 caption 后 label，删除 duplicate captions；确认 `\usepackage{booktabs}` 保留。
6. 清理 Endmatter/Supplementary 模板；补 data availability、reproducibility、GitHub URL 和 AI-use statement。
7. 修复 clean-run 入口后，从最新 `main.tex` 重新编译 PDF，并逐页检查交叉引用、字号、表格宽度和 page limit。

## 10. 仍未由现有分析解决的事项

以下项目不应在论文中写成已经完成：

- arbitrary neutral/jittered multi-start optimisation；
- formal global-optimum analysis；
- time-rescaling or other point-process residual diagnostics；
- parameter standard errors、confidence intervals、bootstrap 或 equivalence margins；
- penalty sensitivity 或完整 `initialisation × specification` grid；
- spatio-temporal ETAS fit；
- forecasting or out-of-sample predictive comparison。

如果提交前不新增这些分析，当前设计仍可成立，但 Discussion 必须把它们保留为限制，而不能暗示已经验证。
