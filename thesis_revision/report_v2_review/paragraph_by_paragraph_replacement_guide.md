# 逐段替换指南：最新版 `main.tex`

本文件不是一般性修改建议。每一项均说明在 `main.tex` 中搜索什么文字、删除到哪里，以及应直接粘贴的完整 LaTeX。未列出的段落可暂时保留。

## 使用规则

1. 先备份当前 `main.tex`。
2. 按本文件编号从前到后替换，避免交叉引用标签重复。
3. 所有 `\caption{}` 后立即放 `\label{}`。
4. 百分号必须写成 `\%`；数值范围使用 `--`；数学量放在 `$...$` 中。
5. preamble 保留 `\usepackage{booktabs}` 和 `\usepackage{graphicx}`。

---

## A. Abstract 与 Introduction

### A1. Abstract

定位：从 `\begin{abstract}` 删除到对应的 `\end{abstract}`。

替换为工作区文件 `abstract_and_introduction.tex` 中的完整 abstract。该版本已经包含：

- 对 triggering parameters 与 estimated background function 的区分；
- 对三种已检验 initialisations 的实际意义；
- 对其他 algorithms、specifications 和 arbitrary starts 的范围限制。

### A2. 删除模板提示并替换 Introduction

定位：删除 `\paragraph{General tips.}` 开始的全部模板文字，以及当前 `\section{Introduction}` 下的 `Tips for the Introduction`。

从 `\section{Introduction}` 到 `\section{Methods}` 之前的内容，替换为 `abstract_and_introduction.tex` 中的完整 Introduction。不要保留模板段落，也不要重复加入第二个 `\section{Introduction}`。

---

## B. Method 的逐段替换

### M1. `Experiment Overview` 整个 subsection

定位：从

```latex
\subsection{Experiment Overview}
```

删除到 `\subsection{Earthquake catalogue and preprocessing}` 之前。

替换为：

```latex
\subsection{Study design}
\label{sec:study_design}

This study assessed the sensitivity of fitted non-stationary temporal ETAS
inference to declustering-derived initial values for the background-rate
function. It was designed as a controlled initialisation experiment rather
than as a comparison of three reduced earthquake catalogues. A single
magnitude-thresholded catalogue, likelihood, background spline, penalty,
triggering-history rule, parameter bounds, optimiser, and set of triggering
parameter starting values were used throughout the main experiment. The only
quantity varied between the three main fits was the initial background-rate
function constructed from Gardner--Knopoff, nearest-neighbour, or
Reasenberg-style declustering. Events classified as foreshocks or aftershocks
were not removed from the ETAS input catalogue.

The analysis followed four linked stages. First, exploratory analysis and
magnitude-completeness diagnostics were used to define a common analytical
catalogue. Second, three structurally different declustering schemes were
applied to that catalogue and their event-level classifications were compared.
Third, the resulting background-like events were converted into three smooth
background-rate starting functions. Finally, the same non-stationary temporal
ETAS model was fitted from each starting function, and the initial and final
between-method differences were compared. A limited model-specification
analysis examined spline flexibility and the length of the triggering history
while holding the nearest-neighbour initialisation fixed.
```

### M2. `Earthquake catalogue and preprocessing` 的第一个数据段落

定位：搜索

```text
In this report, an earthquake catalogue refers to
```

删除该段直至 `All records were labelled as local earthquakes...` 所在段落结束。

替换为：

```latex
An earthquake catalogue is a structured record in which each row represents a
detected event and contains its origin time, location, depth, magnitude, and
associated identifiers. The catalogue used in this study was obtained through
the Southern California Earthquake Data Center (SCEDC) SCSN catalogue-search
interface \citep{scedc2013,scedc_catalog_search2026}. The query covered 20 May
2000 to 20 May 2026, with requested limits of $30.0^\circ$--$39.0^\circ$ N,
$124.0^\circ$--$111.0^\circ$ W, magnitude $-1.0$--$9.0$, and focal depth
$-5.0$--$30.0$ km. These values were search limits rather than the observed
extrema of the returned events.

The downloaded file contained 167,584 records. Its observed ranges were
$34.500^\circ$--$36.500^\circ$ N for latitude,
$-118.99983^\circ$--$-116.00400^\circ$ for longitude, $-0.98$--$7.10$ for
magnitude, and $-2.3$--$28.9$ km for depth. All records were labelled as local
earthquakes, and event identifiers were unique.
```

### M3. 同一 subsection 的 EDA 段落

定位：搜索

```text
Exploratory data analysis was performed before selecting
```

整段替换为：

```latex
Exploratory data analysis was performed before selecting the completeness
threshold. Annual and monthly event counts and the cumulative event count were
used to identify temporal changes in recorded activity. Daily counts in a
30-day window around the catalogue's largest event were used to describe
short-term concentration, and epicentral maps were used to examine geographic
coverage. The exceptionally large number of events recorded in 2019, which
includes the Ridgecrest sequence, motivated separate 2019 and non-2019
completeness and declustering summaries \citep{ross2019hierarchical}. These
analyses were descriptive and were used to identify catalogue features that
required explicit treatment in the subsequent sensitivity experiment.
```

### M4. `Magnitude completeness assessment` 整个 subsection

定位：从

```latex
\subsection{Magnitude completeness assessment}
```

删除到 `\subsection{Assessment of triggering and background non-stationarity}` 之前。

替换为：

```latex
\subsection{Magnitude completeness assessment}
\label{method23}

The magnitude of completeness, $M_c$, is the lowest magnitude above which an
earthquake catalogue is treated as sufficiently complete for a specified
region and observation period. It is an observational threshold rather than a
physical lower limit on earthquake occurrence. Events below $M_c$ may be
preferentially missed when monitoring performance changes or when an intense
earthquake sequence temporarily reduces the detectability of small events.
Such omissions can distort estimated event rates, declustering relationships,
and ETAS parameters. Catalogue completeness was therefore assessed before
constructing the common analytical catalogue.

The primary estimator was the Woessner--Wiemer $b$-value stability method
(MBS-WW) \citep{woessner2005assessing}. This method is based on the
Gutenberg--Richter frequency--magnitude relationship
\citep{gutenberg1944frequency},
\[
\log_{10}N(M)=a-bM,
\]
where $N(M)$ is the cumulative number of earthquakes with magnitude at least
$M$, and $b$ describes how rapidly earthquake frequency decreases as magnitude
increases. Above an appropriate completeness threshold, the estimated
$b$-value should remain reasonably stable as the lower magnitude threshold is
increased. Comparative work has shown that the performance of catalogue-based
$M_c$ estimators can vary under heterogeneous completeness, while MBS-WW can
provide comparatively reliable estimates when sufficiently dense data are
available \citep{wang2025magnitude}.

Candidate thresholds were considered in increments of $\Delta M=0.1$. For a
candidate threshold $m$, the $b$-value was estimated from earthquakes with
$M\geq m$ using
\[
\hat b(m)
=
\frac{\log_{10}(e)}
{\bar M_m-\left(m-\Delta M/2\right)},
\]
where $\bar M_m$ is the mean magnitude of the $n_m$ events at or above the
candidate threshold. Candidates with fewer than 30 events were excluded from
the stability assessment.

The uncertainty of the estimated $b$-value was calculated using the
Shi--Bolt standard-error formula \citep{shi1982standard},
\[
\delta b(m)
=
\ln(10)\,\hat b(m)^2
\sqrt{
\frac{
\sum_{i=1}^{n_m}\left(M_i-\bar M_m\right)^2
}{
n_m(n_m-1)
}
}.
\]
The factor $\ln(10)\approx2.3026$ arises because the Gutenberg--Richter
relationship is expressed using base-10 logarithms.

For each candidate threshold $m$, the estimated $\hat b(m)$ was compared with
the mean of the estimated $b$-values between $m$ and $m+0.5$. Let
\[
\bar b_{0.5}(m)
=
\frac{1}{|\mathcal{M}_m|}
\sum_{m'\in\mathcal{M}_m}\hat b(m'),
\qquad
\mathcal{M}_m
=
\{m':m\leq m'\leq m+0.5\}.
\]
A candidate was classified as stable when
\[
\left|\hat b(m)-\bar b_{0.5}(m)\right|
\leq
\delta b(m).
\]
MBS-WW selected the lowest candidate threshold satisfying this condition. In
practical terms, this identifies the point at which further increases in the
magnitude threshold no longer produce changes in the estimated $b$-value that
are large relative to its estimated uncertainty.

The MBS-WW candidate was then assessed using the Gutenberg--Richter
goodness-of-fit test (GFT) of \citet{wiemer2000minimum}. The GFT compares the
observed cumulative frequency--magnitude distribution with the distribution
predicted by the fitted Gutenberg--Richter relationship. For magnitude bins
$M_k\geq m$, the fitted cumulative count was
\[
N_{\mathrm{GR}}(M_k)
=
N_{\mathrm{obs}}(m)
10^{-\hat b(m)(M_k-m)},
\]
where $N_{\mathrm{obs}}(M_k)$ is the observed cumulative number of events with
magnitude at least $M_k$. Agreement was summarised by
\[
R(m)
=
100\left[
1-
\frac{
\sum_k
\left|
N_{\mathrm{obs}}(M_k)-N_{\mathrm{GR}}(M_k)
\right|
}{
\sum_k N_{\mathrm{obs}}(M_k)
}
\right].
\]
A value close to 100\% indicates close agreement between the observed and
fitted cumulative distributions. The GFT was not calculated when fewer than
50 events were available above a candidate threshold.

Following the conventional interpretation of the GFT, 95\% was treated as the
primary criterion for strong agreement with the Gutenberg--Richter
distribution, while 90\% was used as a less stringent secondary diagnostic
\citep{wiemer2000minimum,woessner2005assessing}. Passing the 90\% criterion was
not treated as equivalent to passing the 95\% criterion. The GFT was used to
evaluate the MBS-WW candidate rather than as an independent estimator that
automatically replaced it.

The exploratory analysis identified 2019 as an unusually active period.
MBS-WW and the GFT were therefore evaluated separately for earthquakes
recorded during 2019 and for those recorded in all other years. This comparison
assessed whether the catalogue-wide estimate concealed a period-specific
change in small-event detection.

Shorter-term variation was examined using chronologically ordered windows of
1,000 successive earthquakes, with each window advanced by 250 events. Within
each window, the MBS-WW estimate was recalculated, together with the
corresponding $b$-value and its uncertainty. Because these were event-based
rather than fixed-duration windows, they covered shorter calendar periods when
earthquake activity was high and longer periods when activity was lower.

A common operational threshold of $M_c=2.0$ was used for all subsequent
declustering comparisons and ETAS fits so that every method was applied to the
same set of events. This threshold was a design compromise between retaining
earthquakes and reducing temporal heterogeneity in completeness. It was not
interpreted as evidence that detection was strictly complete above
$M_c=2.0$ during every part of the observation period, particularly during
2019.
```

### M5. `Assessment of triggering and background non-stationarity` 整个 subsection

定位：从该 subsection 标题删除到 `\subsection{Earthquake declustering}` 之前。

替换为：

```latex
\subsection{Assessment of triggering and background non-stationarity}
\label{sec:nonstationarity_screening}

A homogeneous Poisson process assumes that events occur independently at one
constant rate. In contrast, a temporal ETAS model is a self-exciting point
process: its instantaneous event rate contains a background component together
with additional contributions from earlier earthquakes. A stationary ETAS
model keeps the background component constant, whereas a non-stationary ETAS
model allows it to change over time \citep{ogata1988statistical}.

Two preliminary model comparisons were used to motivate the final model.
First, a homogeneous Poisson process with one rate parameter was compared with
a stationary temporal ETAS model containing one constant background rate and
four triggering parameters. Second, the stationary ETAS model was compared
with a deliberately simplified diagnostic model in which the background rate
was piecewise constant over approximately five-year intervals, while the four
triggering parameters were shared across intervals. The 26-year observation
period produced six background intervals and ten fitted parameters in total.
Here, ``diagnostic'' means that this piecewise model was used only to determine
whether allowing temporal variation in the background improved relative fit;
its interval-specific background estimates were not used as the final
background function.

All comparison models used the same $M\geq2.0$ catalogue and the available
within-catalogue triggering history. Models were compared using maximised
log-likelihood, Akaike information criterion (AIC), and Bayesian information
criterion (BIC). Lower AIC or BIC indicates better relative support after
accounting for the number of fitted parameters. These comparisons can show
whether one candidate specification fits relatively better than another, but
they cannot establish that the final spline model reproduces the observed
event-time pattern adequately. That question would require separate
point-process residual diagnostics.
```

### M6. `Declustering design` 中二元标签定义段

定位：从

```text
For subsequent comparison, the output of each method was reduced
```

删除到 `not observed ground truth.` 所在段结束。

替换为：

```latex
For a common comparison, the method-specific outputs were reduced to the
binary indicator
\[
B_i^{(d)}=
\begin{cases}
1, & \text{event $i$ is background or the mainshock of a sequence under method $d$},\\
0, & \text{event $i$ is classified as a foreshock or aftershock under method $d$},
\end{cases}
\]
where $d\in\{\mathrm{GK},\mathrm{NN},\mathrm{R}\}$. Events with
$B_i^{(d)}=1$ are described as background/mainshock-like, and events with
$B_i^{(d)}=0$ as foreshock/aftershock-like sequence members. These are
algorithmic classifications used to harmonise the three outputs; they are not
observed causal labels or physical ground truth.
```

### M7. Gardner--Knopoff 公式后的解释句

定位并删除：

```latex
The M represent the magnitude of earthquake.
```

替换为：

```latex
Here, $M$ denotes the magnitude of the potential mainshock, $D(M)$ is the
spatial window in kilometres, and $T(M)$ is the temporal window in days.
```

### M8. Nearest-neighbour mixture 段落

定位：搜索

```text
A two-component Gaussian mixture was fitted
```

将该段替换为：

```latex
A two-component Gaussian mixture was fitted to the finite
$\log_{10}\eta_i$ values using a fixed random seed of 123. The two components
represented relatively close and relatively distant nearest-neighbour links.
The value at which the fitted components had equal posterior probability was
used as the classification threshold. Links below the threshold were retained
as sequence links. The fitted component means and resulting threshold were
treated as analysis results and are reported in
Section~\ref{results_declustering}.
```

### M9. `Comparison of declustering outputs and construction...` 整个 subsection

定位：从该 subsubsection 标题删除到 `\subsection{Non-stationary ETAS model}` 之前。

替换为：

```latex
\subsubsection{Comparison of declustering outputs and construction of
non-stationary ETAS initialisations}
\label{method_declustering_comparison}

For each method, the number and proportion of background/mainshock-like events
were recorded. Pairwise raw agreement measured the proportion of events given
the same binary label. Because a high raw agreement can arise when nearly all
events belong to one category, two additional summaries were calculated.
Jaccard similarity measured overlap between the two sets of
background/mainshock-like events, and Cohen's $\kappa$ measured agreement
beyond that expected by chance. The proportion of events without a unanimous
three-method label was also calculated. These comparisons were repeated for
2019 and non-2019 events because the exploratory analysis identified unusually
high activity during 2019.

For method $d$, background/mainshock-like events were aggregated by calendar
month to give a raw monthly rate in events per day. Pairwise mean absolute
error (MAE), root mean squared error (RMSE), and Pearson correlation were
calculated for the raw monthly series over the full period and within 2019.
These were descriptive summaries of difference in level and temporal pattern;
no null hypothesis or $P$-value was attached to them.

For declustering method $d$ and month $m$, let $n_{dm}$ be the number of
background/mainshock-like events and let $E_m$ be the number of observed days
in that month. The ETAS observation interval ran from the first to the last
retained $M\geq M_c$ event. Thus, $E_m$ was the full number of days for an
interior month but only the observed part of the first and last months. This
exposure adjustment prevented partial boundary months from being treated as
complete calendar months.

The monthly value used to construct the starting background was
\[
y_{dm}
=
\log\left(\frac{n_{dm}+0.1}{E_m}\right).
\]
The logarithm placed the initial rates on the same scale as the final model,
which represents $\log\mu(t)$ as a spline and therefore guarantees
$\mu(t)>0$. Months with no background/mainshock-like events had an unadjusted
rate of zero, whose logarithm is undefined. Adding 0.1 made every initial rate
finite; it was a numerical pseudo-count used only to construct initial values
and was not treated as an observed earthquake.

The ETAS optimiser required spline coefficients rather than a separate rate
for every month. The monthly log rates were therefore projected onto the same
natural-cubic-spline basis used in the final model. With $t_m$ denoting the
midpoint of the observed portion of month $m$, the starting coefficients were
\[
\boldsymbol\beta_d^{(0)}
=
\arg\min_{\boldsymbol\beta}
\left\{
\sum_m
\left[y_{dm}-B(t_m)\boldsymbol\beta\right]^2
+10^{-4}\lVert\boldsymbol\beta\rVert_2^2
\right\}.
\]
This regression did not estimate a separate scientific relationship. It
converted the monthly summaries into the eight coefficients required to form
the smooth starting function
$\mu_d^{(0)}(t)=\exp\{B(t)\boldsymbol\beta_d^{(0)}\}$. The small ridge term
stabilised the matrix calculation if spline-basis columns were nearly
dependent. Both the pseudo-count and ridge term affected starting values only;
they did not alter the catalogue or ETAS likelihood.

For a like-for-like propagation comparison, the three smooth starting
functions were evaluated on the same daily grid as the final fitted background
functions. MAE, RMSE, and correlation were calculated on this common grid.
This replaced an earlier diagnostic that used monthly rates in the denominator
and daily fitted rates in the numerator.

Declustering was used only to construct the alternative starting background
functions. No event classified as a foreshock or aftershock was removed. All
three ETAS fits used the same complete $M\geq2.0$ catalogue and the same
likelihood; only the starting background function differed. This controlled
design isolated whether declustering-induced differences persisted after ETAS
estimation.
```

### M10. `Non-stationary ETAS model` 整个 subsection

定位：从该 subsection 标题删除到 `\subsection{Model-specification robustness}` 之前。

替换为：

```latex
\subsection{Non-stationary temporal ETAS model}
\label{method_etas}

The temporal epidemic-type aftershock sequence (ETAS) model is a self-exciting
point-process model. Its conditional intensity $\lambda(t)$ is the expected
instantaneous event rate, in events per day, given the earthquakes observed
before time $t$ \citep{ogata1988statistical}. The intensity is separated into
a background component and contributions from earlier events. In this model,
``background'' means the part of the fitted rate not attributed to observed
earlier catalogue events; it is not a directly observed physical category.

For event times $t_i$ and magnitudes $M_i\geq M_c$, the fitted intensity was
\begin{equation}
\lambda(t_i)
=
\mu(t_i)
+
\sum_{j<i:\,0<t_i-t_j\leq H}
K\exp\{\alpha(M_j-M_c)\}
(t_i-t_j+c)^{-p},
\label{eq:etas_intensity}
\end{equation}
where $\mu(t_i)$ is the background rate, $K$ is triggering productivity,
$\alpha$ controls how productivity changes with the earlier event's magnitude,
$c$ regulates the model at very short time lags, $p$ is the temporal decay
exponent, and $H$ is the maximum triggering-history length. The sum represents
the combined contribution of all eligible earlier events. Locations and depths
were used by the declustering methods but did not enter this temporal ETAS
likelihood.

The main analysis used $H=3,650$ days. Events before the start of the supplied
catalogue were unavailable, so only earlier events observed within the
catalogue contributed to triggering. The cutoff was a modelling and
computational approximation, not a claim that older earthquakes have no
physical effect.

The model was non-stationary because the background rate was allowed to vary
with time. Its logarithm was represented using a natural cubic spline,
\begin{equation}
\log\mu(t)
=
B(t)\boldsymbol\beta,
\label{eq:etas_background}
\end{equation}
where $B(t)$ was an eight-degree-of-freedom spline basis, including an
intercept and boundary knots at the first and last event times. The log link
ensured that $\mu(t)$ remained positive. The spline provided a smooth,
low-dimensional representation of temporal variation rather than estimating
an unrelated rate at every time point
\citep{kumazawa2014nonstationary,kattamanchi2017nonstationary}.

To discourage an unnecessarily irregular background curve, the objective
included the fixed roughness penalty
\[
P(\boldsymbol\beta)
=
\lambda_{\mathrm{pen}}
\sum_k(\Delta^2\beta_k)^2,
\qquad
\lambda_{\mathrm{pen}}=1.
\]
The second differences $\Delta^2\beta_k$ measure changes between neighbouring
spline coefficients, so the penalty discourages abrupt changes in the fitted
background while still allowing gradual temporal variation.

Let $T$ be the endpoint of the observation period. The unpenalised
point-process log-likelihood was
\[
\ell(\theta)
=
\sum_i\log\lambda(t_i)
-
\int_0^T\lambda(s)\,\mathrm{d}s.
\]
The first term rewards models assigning high intensity to observed event
times, while the integral term penalises models that predict too many events
over the observation period. The background integral was evaluated on a daily
grid using trapezoidal integration. For each parent event, the triggering
integral was evaluated analytically up to the earlier of $T$ and $H$ days after
that event. The optimised objective was
$-\ell(\theta)+P(\boldsymbol\beta)$; the unpenalised log-likelihood and penalty
were retained separately for reporting.

The parameters $K$, $\alpha$, and $c$ were estimated on logarithmic scales,
and $p=1+\exp(q)$ enforced $p>1$. Bounds on the natural scales were
$10^{-6}\leq K\leq10$, $0.05\leq\alpha\leq5$,
$10^{-4}\leq c\leq10$ days, and $1.001\leq p\leq6$; each spline coefficient
was bounded between $-20$ and 5. Optimisation used L-BFGS-B with at most 300
iterations. The common triggering starting values were $K=0.02$,
$\alpha=1.0$, $c=0.01$ days, and $p=1.1$. Only the spline coefficient vector
$\boldsymbol\beta_d^{(0)}$ differed between the three main fits.

The optimiser status, unpenalised log-likelihood, penalty, penalised objective,
and final triggering estimates were compared. Final background functions were
evaluated on the common daily grid. For each pair of declustering methods, the
propagation ratios were
\[
R_{\mathrm{MAE}}
=
\frac{
\mathrm{MAE}\{\hat\mu_{d_1},\hat\mu_{d_2}\}
}{
\mathrm{MAE}\{\mu^{(0)}_{d_1},\mu^{(0)}_{d_2}\}
},
\qquad
R_{\mathrm{RMSE}}
=
\frac{
\mathrm{RMSE}\{\hat\mu_{d_1},\hat\mu_{d_2}\}
}{
\mathrm{RMSE}\{\mu^{(0)}_{d_1},\mu^{(0)}_{d_2}\}
}.
\]
Values near zero indicate that little of the initial between-method difference
remained after fitting. These ratios measure numerical sensitivity for the
observed catalogue; they are not formal equivalence tests and do not represent
sampling uncertainty.
```

### M11. `Model-specification robustness` 的两个正文段落

定位：保留 subsection 标题，但将其下两个现有正文段落全部替换为：

```latex
A limited robustness analysis assessed sensitivity to two modelling choices:
the flexibility of the spline background and the length of the triggering
history. The nearest-neighbour-derived initialisation and
$\lambda_{\mathrm{pen}}=1$ were held fixed. The baseline used eight spline
degrees of freedom and $H=3,650$ days. Three alternatives used six or ten
spline degrees of freedom with the baseline cutoff, or eight degrees of
freedom with $H=7,300$ days.

For each specification, optimiser status, unpenalised log-likelihood,
penalised objective, and triggering parameters $K$, $\alpha$, $c$, and $p$
were compared with the baseline. Because spline coefficient vectors have
different dimensions when the degrees of freedom change, they were not
compared coefficient by coefficient. Instead, the estimated background
functions were evaluated on the same daily grid and compared using MAE, RMSE,
and Pearson correlation. As all alternative specifications were fitted only
from the nearest-neighbour start, this analysis assessed specification
sensitivity and did not repeat the full initialisation experiment under every
specification.
```

---

## C. Results 的逐段替换

### R1. `Exploratory catalogue characteristics and magnitude completeness` 整个 subsection

定位：从该 subsection 标题删除到 `\subsection{Evidence for triggering and a time-varying background}` 之前。

替换为：

```latex
\subsection{Exploratory catalogue characteristics and magnitude completeness}
\label{results_eda_mc}

The pre-threshold catalogue contained 167,584 earthquakes. The year 2019
contained 50,276 events, corresponding to 30.0\% of the raw catalogue. The
largest recorded event had magnitude 7.1 and occurred on 6 July 2019 during
the Ridgecrest sequence \citep{ross2019hierarchical}. Annual counts show that
2019 was exceptional relative to the surrounding years
(Figure~\ref{fig:eda_temporal}), motivating period-specific completeness and
declustering summaries.

\begin{figure}[htbp]
  \centering
  \includegraphics[width=0.88\textwidth]{annual_event_counts_revised.pdf}
  \caption{Annual earthquake counts in the pre-threshold catalogue. The year
  2019 contained 50,276 recorded events and is highlighted. The years 2000 and
  2026 are partial observation years.}
  \label{fig:eda_temporal}
\end{figure}

The catalogue-wide MBS-WW estimate was approximately $M_c=1.7$. The
period-specific estimates were $M_c=3.3$ for 2019 and $M_c=1.9$ for all other
years combined. At $M_c=1.7$, the GFT values were 92.6\% for 2019 and 90.8\%
for non-2019 events. Both met the secondary 90\% criterion but neither met the
primary 95\% criterion. The rolling analysis showed that the highest MBS-WW
estimates were concentrated around 2019 rather than remaining uniformly high
throughout the catalogue (Figure~\ref{fig:rolling_mc}).

\begin{figure}[htbp]
  \centering
  \includegraphics[width=0.88\textwidth]{rolling_mc_revised.pdf}
  \caption{MBS-WW estimates of magnitude completeness in windows of 1,000
  successive earthquakes, advanced by 250 events. The shaded band marks 2019;
  the horizontal lines show the catalogue-wide estimate of 1.7 and the common
  operational threshold of 2.0. Closely spaced window midpoints during 2019
  reflect the unusually high event rate.}
  \label{fig:rolling_mc}
\end{figure}

A common operational threshold of $M_c=2.0$ was adopted to preserve one
analytical catalogue for the controlled comparison. This was above the
catalogue-wide estimate of 1.7 but below the separate 2019 estimate of 3.3. It
therefore represented a compromise between retaining events and reducing
temporal heterogeneity; it did not establish complete detection throughout
2019. Applying the threshold retained 15,207 earthquakes, of which 6,865
(45.1\%) occurred during 2019.
```

### R2. `Evidence for triggering and a time-varying background` 整个 subsection

定位：从该 subsection 标题删除到 `\subsection{Declustering-induced classification differences}` 之前。

替换为：

```latex
\subsection{Evidence for triggering and a time-varying background}
\label{results_nonstationarity}

The preliminary comparisons supported both earthquake-to-earthquake
triggering and temporal variation in the background rate. The homogeneous
Poisson model had log-likelihood $-8{,}046.72$, whereas the stationary ETAS
model had log-likelihood $27{,}577.37$ (Table~\ref{tab:model_comparison}). The
corresponding AIC and BIC differences were 71,240.18 and 71,209.66,
respectively, in favour of stationary ETAS. A constant-rate process without
triggering was therefore a poor relative description of the event times.

Allowing the diagnostic ETAS background to vary across approximately
five-year intervals increased the log-likelihood to $27{,}673.39$. Relative to
stationary ETAS, AIC decreased by 182.03 and BIC by 143.89 despite the five
additional parameters. This supported using a time-varying background in the
main spline model. It did not, by itself, establish absolute goodness of fit
for that final model.

\begin{table}[htbp]
  \centering
  \caption{Preliminary model comparisons used to motivate earthquake
  triggering and temporal variation in the background rate. Lower AIC and BIC
  indicate better relative fit. The diagnostic piecewise-background ETAS was
  not the final spline model.}
  \label{tab:model_comparison}
  \begin{tabular}{lrrrr}
    \toprule
    Model & Log-likelihood & $k$ & AIC & BIC \\
    \midrule
    Homogeneous Poisson
      & $-8046.72$ & 1 & $16095.43$ & $16103.06$ \\
    Stationary ETAS
      & $27577.37$ & 5 & $-55144.75$ & $-55106.60$ \\
    Diagnostic piecewise-background ETAS
      & $27673.39$ & 10 & $-55326.78$ & $-55250.49$ \\
    \bottomrule
  \end{tabular}
\end{table}
```

### R3. `Declustering-induced classification differences` 整个 subsection

定位：从该 subsection 标题删除到 `\subsection{Declustering-derived starting backgrounds}` 之前。

替换为：

```latex
\subsection{Declustering-induced classification differences}
\label{results_declustering}

The three methods assigned different background/mainshock-like fractions to
the same 15,207 earthquakes (Table~\ref{tab:declustering_summary}).
Gardner--Knopoff classified 2,599 events (17.1\%) as
background/mainshock-like, compared with 4,573 (30.1\%) for nearest-neighbour
and 3,847 (25.3\%) for Reasenberg-style. The largest difference, between
Gardner--Knopoff and nearest-neighbour, was 13.0 percentage points. Across all
three methods, 2,514 events (16.5\%) did not receive a unanimous label.

\begin{table}[htbp]
  \centering
  \caption{Binary classifications from the three declustering methods using
  the common $M\geq2.0$ catalogue.}
  \label{tab:declustering_summary}
  \begin{tabular}{lrrr}
    \toprule
    Method
      & Background/mainshock-like
      & Foreshock/aftershock-like
      & Background/mainshock-like (\%) \\
    \midrule
    Gardner--Knopoff & 2599 & 12608 & 17.1 \\
    Nearest-neighbour & 4573 & 10634 & 30.1 \\
    Reasenberg-style & 3847 & 11360 & 25.3 \\
    \bottomrule
  \end{tabular}
\end{table}

Pairwise raw agreement ranged from 0.852 to 0.915
(Table~\ref{tab:declustering_agreement}). Agreement on the less common
background/mainshock-like category was weaker: Jaccard similarity ranged from
0.523 to 0.734. Nearest-neighbour and Reasenberg-style were most similar
($\kappa=0.789$, Jaccard $=0.734$), whereas Gardner--Knopoff and
nearest-neighbour were furthest apart ($\kappa=0.600$, Jaccard $=0.523$).

\begin{table}[htbp]
  \centering
  \caption{Pairwise event-level agreement. Jaccard similarity is calculated
  for the background/mainshock-like sets.}
  \label{tab:declustering_agreement}
  \begin{tabular}{lrrr}
    \toprule
    Comparison & Raw agreement & Jaccard & Cohen's $\kappa$ \\
    \midrule
    GK vs NN & 0.852 & 0.523 & 0.600 \\
    GK vs Reasenberg & 0.902 & 0.624 & 0.709 \\
    NN vs Reasenberg & 0.915 & 0.734 & 0.789 \\
    \bottomrule
  \end{tabular}
\end{table}

Class imbalance was particularly strong in 2019. The
background/mainshock-like fractions were 0.9\% for Gardner--Knopoff, 3.6\%
for nearest-neighbour, and 2.0\% for Reasenberg-style. Pairwise raw agreement
therefore exceeded 0.972 because the methods almost always agreed on the much
larger foreshock/aftershock-like category. In contrast, the 2019 Jaccard
similarities were 0.235 for GK--NN, 0.420 for GK--Reasenberg, and 0.512 for
NN--Reasenberg, with corresponding $\kappa$ values of 0.371, 0.586, and
0.668. High raw agreement therefore did not imply that the methods selected
the same small set of background/mainshock-like events.

For the nearest-neighbour construction, the Gutenberg--Richter estimate from
the common catalogue was $\hat b=0.8152$. The fitted unequal-variance
two-component Gaussian mixture had component means $-5.120$ and $-1.432$ on
the $\log_{10}\eta$ scale. The equal-posterior-probability rule gave the
classification threshold $-2.740$. These values were used to create the
nearest-neighbour links reported above.
```

### R4. `Declustering-derived starting backgrounds` 整个 subsection

定位：从该 subsection 标题删除到 `\subsection{Sensitivity of fitted ETAS inference to initialisation}` 之前。

替换为：

```latex
\subsection{Declustering-derived starting backgrounds}
\label{results_initialisation}

The classification differences propagated into the raw monthly background
rates. Over the full period, the GK--NN comparison had MAE 0.209 events/day,
RMSE 0.324 events/day, and correlation 0.490. GK--Reasenberg had MAE 0.133,
RMSE 0.180, and correlation 0.687, while NN--Reasenberg had MAE 0.091, RMSE
0.185, and correlation 0.841. During 2019, the GK--NN differences increased to
MAE 0.502 and RMSE 0.632 events/day, with correlation $-0.167$.

Table~\ref{tab:initial_background_metrics} reports the smooth functions that
were actually supplied to the ETAS optimiser, evaluated on the common daily
grid after the pseudo-count, exposure adjustment, spline projection, and ridge
stabilisation. The largest difference remained GK--NN, with MAE 0.203 and RMSE
0.295 events/day. Its correlation was 0.302, showing that the two starting
functions differed in temporal shape as well as level. NN--Reasenberg had the
smallest MAE (0.076) and the highest correlation (0.931). The main experiment
therefore did not begin from nearly identical background functions.

\begin{table}[htbp]
  \centering
  \caption{Pairwise differences between the smooth background functions used
  as ETAS starting values, evaluated on the common daily grid. MAE and RMSE
  are in events/day; correlation is dimensionless.}
  \label{tab:initial_background_metrics}
  \begin{tabular}{lrrr}
    \toprule
    Comparison & MAE & RMSE & Correlation \\
    \midrule
    GK vs NN & 0.203052 & 0.294838 & 0.302010 \\
    GK vs Reasenberg & 0.129409 & 0.155881 & 0.503749 \\
    NN vs Reasenberg & 0.076162 & 0.158183 & 0.930644 \\
    \bottomrule
  \end{tabular}
\end{table}
```

### R5. `Sensitivity of fitted ETAS inference to initialisation` 整个 subsection

定位：从该 subsection 标题删除到 `\subsection{Robustness to ETAS model specification}` 之前。

替换为：

```latex
\subsection{Sensitivity of fitted ETAS inference to initialisation}
\label{results_etas}

The optimiser returned status code 0 for all three main fits, indicating
successful termination under its stopping rule. The triggering estimates and
unpenalised log-likelihoods were very similar
(Table~\ref{tab:etas_parameters}). The largest pairwise ranges were
$2.09\times10^{-5}$ for $K$, $5.53\times10^{-4}$ for $\alpha$,
$2.50\times10^{-5}$ days for $c$, and $3.98\times10^{-4}$ for $p$. The
unpenalised log-likelihoods differed by 0.121. The penalised objectives were
$-27765.576$, $-27765.584$, and $-27765.582$ for the GK, NN, and
Reasenberg-style starts, respectively, a total range below 0.009. These
differences do not constitute a formal equivalence test, but they show low
numerical sensitivity across the three starts examined.

\begin{table}[htbp]
  \centering
  \caption{Non-stationary temporal ETAS fits from the three
  declustering-derived background starts. Log-likelihood is unpenalised.}
  \label{tab:etas_parameters}
  \begin{tabular}{lrrrrr}
    \toprule
    Initialisation & Log-likelihood & $K$ & $\alpha$ & $c$ & $p$ \\
    \midrule
    Gardner--Knopoff
      & 27767.432 & 0.020154 & 1.588437 & 0.009616 & 1.122273 \\
    Nearest-neighbour
      & 27767.466 & 0.020157 & 1.588518 & 0.009610 & 1.122048 \\
    Reasenberg-style
      & 27767.553 & 0.020175 & 1.587965 & 0.009591 & 1.121875 \\
    \bottomrule
  \end{tabular}
\end{table}

The final fitted background functions were also nearly coincident. Pairwise
correlations exceeded 0.999991 and final MAEs ranged from 0.000779 to 0.001063
events/day. Comparing the starting and final functions on the same daily grid
gave MAE ratios of 0.00384, 0.00822, and 0.01174
(Table~\ref{tab:propagation}). Thus, 0.38\%--1.17\% of the initial pairwise
MAE remained, corresponding to attenuation of 98.83\%--99.62\%. RMSE
attenuation ranged from 99.10\% to 99.68\%.

\begin{figure}[htbp]
  \centering
  \includegraphics[width=0.90\textwidth]{initial_vs_final_background_revised.pdf}
  \caption{Declustering-derived smooth background-rate starting functions
  (top) and the corresponding final fitted non-stationary ETAS background
  functions (bottom), evaluated on the same daily grid. The shaded band marks
  2019. The starting functions differ in level and temporal pattern, whereas
  the three final functions are almost coincident.}
  \label{fig:initial_final_background}
\end{figure}

\begin{table}[htbp]
  \centering
  \caption{Propagation of differences between the starting and final
  background-rate functions on a common daily grid. MAE and RMSE are in
  events/day; ratios are final difference divided by starting difference.}
  \label{tab:propagation}
  \resizebox{\textwidth}{!}{%
  \begin{tabular}{lrrrrrrr}
    \toprule
    Comparison
      & Initial MAE & Final MAE & MAE ratio
      & Initial RMSE & Final RMSE & RMSE ratio & Final $r$ \\
    \midrule
    GK vs NN
      & 0.203052 & 0.000779 & 0.003835
      & 0.294838 & 0.000951 & 0.003224 & 0.999991 \\
    GK vs Reasenberg
      & 0.129409 & 0.001063 & 0.008216
      & 0.155881 & 0.001396 & 0.008957 & 0.999993 \\
    NN vs Reasenberg
      & 0.076162 & 0.000894 & 0.011744
      & 0.158183 & 0.001177 & 0.007442 & 0.999993 \\
    \bottomrule
  \end{tabular}}
\end{table}

Taken together, the main experiment showed substantial differences in event
classifications and starting background functions, followed by only very small
differences in the fitted triggering parameters and background functions under
the baseline model.
```

### R6. `Robustness to ETAS model specification` 整个 subsection

定位：从该 subsection 标题删除到 `\section{Discussion}` 之前。

替换为：

```latex
\subsection{Robustness to ETAS model specification}
\label{results_robustness}

The optimiser returned status code 0 for all three alternative specifications
(Table~\ref{tab:robustness_parameters}). Relative to the eight-degree-of-freedom
baseline, six degrees of freedom reduced the log-likelihood by 15.20, ten
degrees of freedom increased it by 2.98, and extending the triggering-history
cutoff to 7,300 days changed it by $-0.44$. Across the alternatives, the
largest percentage changes were 0.41\% for $K$, 0.12\% for $\alpha$, 3.38\%
for $c$, and 0.42\% for $p$. The triggering-parameter estimates therefore
remained similar under the tested changes.

\begin{table}[htbp]
  \centering
  \caption{Triggering-parameter estimates in the limited model-specification
  analysis. All models used the nearest-neighbour-derived starting function.}
  \label{tab:robustness_parameters}
  \begin{tabular}{lrrrrr}
    \toprule
    Specification & Log-likelihood & $K$ & $\alpha$ & $c$ & $p$ \\
    \midrule
    Baseline (df 8; 3650 d)
      & 27767.466 & 0.020157 & 1.588518 & 0.009610 & 1.122048 \\
    df 6; 3650 d
      & 27752.268 & 0.020202 & 1.587419 & 0.009935 & 1.126762 \\
    df 10; 3650 d
      & 27770.450 & 0.020239 & 1.586574 & 0.009495 & 1.120410 \\
    df 8; 7300 d
      & 27767.029 & 0.020155 & 1.588306 & 0.009673 & 1.123139 \\
    \bottomrule
  \end{tabular}
\end{table}

The estimated background functions also remained strongly correlated with the
baseline function, although their shapes were more affected by spline
flexibility than by the triggering-history length
(Table~\ref{tab:robustness_background}). The six-degree-of-freedom function
had MAE 0.0280 events/day and correlation 0.9825 relative to baseline; the
ten-degree-of-freedom function had MAE 0.0128 and correlation 0.9966. Extending
the cutoff produced the smallest difference, with MAE 0.00536 and correlation
0.99979. Because only the nearest-neighbour start was used, these results do
not establish initialisation robustness under the alternative specifications.

\begin{table}[htbp]
  \centering
  \caption{Differences between each alternative fitted background function and
  the baseline background function on the common daily grid. MAE and RMSE are
  in events/day.}
  \label{tab:robustness_background}
  \begin{tabular}{lrrr}
    \toprule
    Specification & MAE & RMSE & Correlation \\
    \midrule
    df 6; 3650 d & 0.027989 & 0.040999 & 0.982543 \\
    df 10; 3650 d & 0.012818 & 0.017872 & 0.996613 \\
    df 8; 7300 d & 0.005357 & 0.006712 & 0.999792 \\
    \bottomrule
  \end{tabular}
\end{table}
```

---

## D. Discussion：逐段替换

当前 Discussion 的主要段落均需调整。最安全的操作是保留
`\section{Discussion}` 和 `\label{sec:discussion}`，然后删除其后至
`\section{Endmatter}` 之前的全部正文，按以下 D1--D11 顺序粘贴。每个编号对应一个独立段落。

### D1. 主要发现

```latex
The three declustering schemes produced materially different event
classifications and background-rate starting functions, but these differences
were almost completely attenuated when the same baseline non-stationary
temporal ETAS model was fitted to the same catalogue. The three fits produced
very similar triggering parameters, unpenalised log-likelihoods, penalised
objectives, and fitted background functions. On the common daily grid, only
0.38\%--1.17\% of the initial pairwise MAE remained. For this catalogue and
baseline specification, fitted ETAS inference was therefore highly robust to
the three declustering-derived initialisations examined.
```

### D2. 三种分类为何不同

```latex
The disagreement is consistent with the methods' different definitions of an
earthquake sequence: Gardner--Knopoff uses fixed magnitude-dependent
space--time windows, the nearest-neighbour construction uses data-derived
space--time--magnitude links, and the Reasenberg-style procedure builds
interaction groups that can expand or merge. Comparative studies have likewise
shown that declustering algorithms and their parameter settings can change the
number and characteristics of retained events
\citep{perry2024comparative,mizrahi2021effect}. In the present analysis, the
background/mainshock-like fraction ranged from 17.1\% to 30.1\%, and 16.5\%
of events lacked a unanimous label. These are meaningful classification
differences even though the final ETAS solution was stable.
```

### D3. 2019 与类别不平衡

```latex
The 2019 results show why the overall percentage of matching labels was not
sufficient to describe agreement. The methods agreed on more than 97\% of all
2019 events because they almost always agreed on the much larger
foreshock/aftershock-like category. They agreed much less on which events
belonged to the smaller background/mainshock-like category. For example, the
GK--NN Jaccard index was only 0.235, meaning that only 23.5\% of the events
selected as background/mainshock-like by either method were selected by both.
The largest raw monthly starting-rate differences also occurred during 2019.
Thus, the exceptional activity identified in the exploratory analysis had a
clear consequence for the construction of the alternative initial values.
```

### D4. Attenuation 的含义与原因

```latex
Here, attenuation refers to the reduction in the differences between the
starting background functions after the ETAS model was fitted. Pairwise MAEs
between the final functions were only 0.38\%--1.17\% of those between the
corresponding starting functions, equivalent to reductions of
98.83\%--99.62\%. Declustering supplied only the initial spline coefficients;
it did not constrain the fitted background to remain close to its starting
curve. After initialisation, every fit re-estimated the background and
triggering components by optimising the same penalised likelihood with the
same events, spline basis, penalty, parameter bounds, and triggering-parameter
starts. The results are therefore consistent with the three
declustering-derived starts leading to essentially the same fitted numerical
solution.
```

### D5. 不能证明 global optimum

```latex
This result does not establish that the penalised likelihood has a unique
global optimum. ETAS likelihood estimation can be computationally difficult
when background and triggering parameters are strongly interdependent, and
performance can depend on parameterisation and starting values
\citep{veen2008estimation,lombardi2015simulated}. Neutral and randomly
perturbed multi-start fits were available as optional code but were not
executed. The experiment therefore demonstrates robustness to the three
scientifically motivated starts examined, rather than independence from all
possible starting values or proof of a global optimum.
```

### D6. GFT 与 completeness limitation

```latex
The completeness diagnostics temper the interpretation of the sensitivity
experiment. At the catalogue-wide MBS-WW estimate of $M_c=1.7$, the GFT was
92.6\% for 2019 and 90.8\% for the remaining years: both met the secondary
90\% criterion but neither met the primary 95\% criterion. Moreover, the
separate 2019 MBS-WW estimate was 3.3. The common $M_c=2.0$ threshold ensured
that all three initialisation experiments used the same events, but it did not
remove uncertainty about small-event detection during the busiest sequence.
Residual incompleteness could affect nearest-neighbour links, window
membership, and the division of ETAS intensity into background and triggering
components.
```

### D7. 与现有研究和 research gap 的联系

```latex
Previous work establishes the two components of the present problem
separately. Comparative declustering studies show that alternative algorithms
can produce different catalogues and derived seismicity summaries
\citep{perry2024comparative,mizrahi2021effect}. Conventional declustering has
also been used to obtain a preliminary background estimate before fitting an
ETAS model to the full catalogue \citep{ogata1998space,zhuang2002stochastic},
while non-stationary ETAS studies have represented the background rate using
penalised splines \citep{kumazawa2014nonstationary,kattamanchi2017nonstationary}.
The present study addresses the narrower question connecting these areas:
whether differences from alternative declustering procedures persist when
they are used only as starting functions for the same penalised
non-stationary temporal ETAS likelihood.
```

### D8. 研究优势

```latex
A principal strength is the controlled comparison. Because the same 15,207
events and the same likelihood were used in all main fits, differences in the
final solutions can be attributed to the starting background functions rather
than to different training catalogues. Jaccard similarity and $\kappa$
supplemented raw agreement so that the dominant foreshock/aftershock-like
category did not conceal disagreement about the smaller
background/mainshock-like set. Comparing the actual smooth starting functions
and final functions on the same daily grid also ensured that the attenuation
ratios had comparable numerators and denominators. These design features make
the low observed sensitivity interpretable rather than merely showing that
three reported parameter tables look similar.
```

### D9. 方法和外推限制

```latex
Several limitations restrict the scope of the findings. The study used one
regional catalogue, so generalisability across tectonic settings, network
histories, and magnitude ranges is unknown. Each declustering method used one
baseline parameterisation, so the experiment did not cover uncertainty across
all plausible parameter values. The nearest-neighbour classification was a
deterministic construction rather than the complete stochastic, spatially
adaptive procedure of \citet{zaliapin2020declustering}, and the
Reasenberg-style implementation was not an exact reproduction of
\texttt{CLUSTER2000}. The fitted ETAS model was temporal rather than
spatio-temporal: locations informed declustering but did not enter the ETAS
background or triggering kernel.
```

### D10. Specification、uncertainty 与 model fit 限制

```latex
The spline degrees of freedom, fixed penalty, parameter bounds, and
3,650-day triggering cutoff were modelling choices. The limited specification
analysis changed only the spline degrees of freedom and cutoff and used only
the nearest-neighbour start. Triggering-parameter estimates changed little,
and the alternative background functions remained strongly correlated with
the baseline, although the background was more responsive to spline
flexibility. This does not establish robustness to the penalty or to every
initialisation--specification combination. In addition, pre-catalogue
triggering history was unavailable. No standard errors, confidence intervals,
bootstrap analysis, or formal equivalence margins were calculated, and no
time-rescaling or other point-process residual diagnostic was reported.
Consequently, close parameter estimates and improved AIC/BIC should not be
interpreted as formal equivalence or proof of absolute model adequacy. The
study also did not evaluate forecasts, so fitted similarity cannot be extended
to predictive equivalence.
```

### D11. 实际意义与最后的 Press Office 段落

```latex
Three approaches to sorting Southern California earthquakes disagreed about
which events represented the underlying background activity, especially
during the intense 2019 sequence. Yet when those alternatives were used only
as starting suggestions for the same statistical model, which then
re-estimated the rates from all earthquakes, the final answers were almost
identical. For this dataset, analysts therefore have practical flexibility in
selecting among the three tested procedures when their outputs are used only
to start this type of model. Other regions and modelling choices should still
be tested before assuming that the same reassuring result will always hold.
```

---

## E. References 的精确新增内容

将 `references_to_add.bib` 中的新条目复制进主 `refs.bib`。本轮 Method 特别需要确认存在以下 keys：

```text
gutenberg1944frequency
shi1982standard
wiemer2000minimum
woessner2005assessing
```

不要在 `refs.bib` 中重复添加已经存在的 citation key。

---

## F. 图文件位置

将下列 PDF 复制到 LaTeX project 与 `main.tex` 相同的图像搜索目录，或在
`\includegraphics{}` 中使用正确相对路径：

```text
thesis_revision/report_v2_review/figures/annual_event_counts_revised.pdf
thesis_revision/report_v2_review/figures/rolling_mc_revised.pdf
thesis_revision/report_v2_review/figures/initial_vs_final_background_revised.pdf
```

如果使用 PNG 版本，必须同步把上述 `\includegraphics` 扩展名改为 `.png`。
