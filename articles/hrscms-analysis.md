# HRS-CMS Contingency Table Analysis

[![Open In
Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/postlink-group/postlink/blob/main/notebooks/hrscms_contingency.ipynb)

## **1. Linked Data Set: HRS-CMS Example**

The HRS-CMS (Health and Retirement Study linked with Centers for
Medicare & Medicaid Services) dataset combines survey and administrative
records to study health and socioeconomic outcomes among older adults.
This linkage enables researchers to analyze the relationships between
health indicators and conditions recorded in different data sources
while accounting for potential errors in record linkage.

A subset of records is considered high-quality based on stricter linkage
criteria, whereas the remaining records are obtained through automated
matching procedures that may introduce linkage error. Consequently,
discrepancies between variables recorded in the two sources can occur,
motivating the use of statistical methods that adjust for
misclassification.

For demonstration purposes, we examine a small aggregated subset of the
HRS-CMS dataset. Each row represents a cell count for one combination of
HRS self-report, CMS administrative record, and linkage mismatch status.

``` r

agg <- data.frame(
  NRSHOM     = as.integer(c(0, 1, 0, 1, 0, 0)),
  mds_nrshom = as.integer(c(0, 0, 1, 1, 0, 1)),
  mismatch   = as.integer(c(0, 0, 0, 0, 1, 1)),
  count      = as.integer(c(271, 9, 5, 15, 56, 3))
)
agg # Data input (aggregated counts only)
```

| NRSHOM | mds_nrshom | mismatch | count |
|-------:|-----------:|---------:|------:|
|      0 |          0 |        0 |   271 |
|      1 |          0 |        0 |     9 |
|      0 |          1 |        0 |     5 |
|      1 |          1 |        0 |    15 |
|      0 |          0 |        1 |    56 |
|      0 |          1 |        1 |     3 |

Aggregated HRS-CMS data {.table .table
style="width: max-content !important; max-width: none !important; margin-bottom: 20px;"}

The dataset includes the following variables:

- `NRSHOM`: indicator of nursing home residence based on HRS
  self-reports (1: yes, 0: no).  
- `mds_nrshom`: indicator of nursing home residence based on CMS
  administrative records (1: yes, 0: no).  
- `mismatch`: indicator of linkage mismatch status (1: mismatched link,
  0: correctly linked record).  
- `count`: number of records in each aggregated cell.

The overall linkage mismatch rate is computed as a weighted proportion
using the aggregated counts:

``` r

mm_rate <- sum(agg$mismatch * agg$count) / sum(agg$count)
data.frame(`Mismatch rate` = round(mm_rate, 4))
```

| Mismatch.rate |
|--------------:|
|        0.1643 |

Overall linkage mismatch rate {.table .table
style="width: max-content !important; max-width: none !important; margin-bottom: 20px;"}

The observed mismatch rate motivates the use of statistical methods that
explicitly adjust for linkage error. Naive analyses that ignore these
mismatches may underestimate or distort associations between variables.
In later sections, we demonstrate the workflow for contingency table
analysis under linkage error using the `postlink` package.

## **2. Naive Approach**

We first use the observed aggregated data without accounting for
potential linkage mismatches. This approach treats the linked data as if
all records were correctly matched.

``` r

# Total number of records represented in the aggregated data
n <- sum(agg$count)

# Naive contingency table in counts
ctab <- xtabs(count ~ NRSHOM + mds_nrshom, data = agg)
as.data.frame.matrix(ctab)
```

|     |   0 |   1 |
|:----|----:|----:|
| 0   | 327 |   8 |
| 1   |   9 |  15 |

Naive contingency table in counts {.table .table
style="width: max-content !important; max-width: none !important; margin-bottom: 20px;"}

``` r

# Naive contingency table in proportions
ctab_prop <- ctab / sum(ctab)
round(as.data.frame.matrix(ctab_prop), 4)
```

|     |      0 |      1 |
|:----|-------:|-------:|
| 0   | 0.9109 | 0.0223 |
| 1   | 0.0251 | 0.0418 |

Naive contingency table in proportions {.table .table
style="width: max-content !important; max-width: none !important; margin-bottom: 20px;"}

The naive table summarizes the observed joint distribution of HRS
self-reported nursing home residence (`NRSHOM`) and CMS administrative
nursing home residence (`mds_nrshom`) without adjusting for linkage
error.

## **3. Exact Approach**

We next restrict the analysis to aggregated cells corresponding to
correctly linked records, where `mismatch = 0`. This provides a
comparison table based only on records without linkage mismatches.

``` r

# Keep only correctly linked records
agg_exact <- subset(agg, mismatch == 0)

# Exact-match contingency table in counts
ctab_exact <- xtabs(count ~ NRSHOM + mds_nrshom, 
                    data = agg_exact)
as.data.frame.matrix(ctab_exact)
```

|     |   0 |   1 |
|:----|----:|----:|
| 0   | 271 |   5 |
| 1   |   9 |  15 |

Exact-match contingency table in counts {.table .table
style="width: max-content !important; max-width: none !important; margin-bottom: 20px;"}

``` r

# Exact-match contingency table in proportions
ctab_exact_prop <- ctab_exact / sum(ctab_exact)
round(as.data.frame.matrix(ctab_exact_prop), 4)
```

|     |      0 |      1 |
|:----|-------:|-------:|
| 0   | 0.9033 | 0.0167 |
| 1   | 0.0300 | 0.0500 |

Exact-match contingency table in proportions {.table .table
style="width: max-content !important; max-width: none !important; margin-bottom: 20px;"}

This table shows the joint distribution of HRS self-reports and CMS
administrative records among correctly linked records only. Comparing
this table with the naive table illustrates how linkage mismatches may
affect the observed relationship between the two variables.

## **4. Adjustment Approach**

We adjust for potential mismatches between HRS self-reports (`NRSHOM`)
and CMS records (`mds_nrshom`) using the observed mismatch rate.

``` r

library(postlink)

# Expand the aggregated counts only for compatibility with postlink
hrscms_clean <- agg[
  rep(seq_len(nrow(agg)), agg$count),
  c("NRSHOM", "mds_nrshom")
]

rownames(hrscms_clean) <- NULL

# Mismatch rate from aggregated counts
mm_rate <- sum(agg$mismatch * agg$count) / sum(agg$count)

# Create adjustment object
adj <- adjMixture(
  linked.data = hrscms_clean,
  m.rate = mm_rate
)

# Adjusted contingency table
adjusted_table <- plctable(
  formula = ~ NRSHOM + mds_nrshom,
  adjustment = adj
)

# Adjusted probability table
round(as.data.frame.matrix(adjusted_table$phat), 4)
```

|     |      0 |      1 |
|:----|-------:|-------:|
| 0   | 0.9139 | 0.0170 |
| 1   | 0.0197 | 0.0494 |

Adjusted contingency table in proportions {.table .table
style="width: max-content !important; max-width: none !important; margin-bottom: 20px;"}

We also demonstrate an alternative lower mismatch rate (not all
non-exact matches are mismatches):

``` r

# Alternative lower mismatch rate
lower_mm_rate <- 0.5 * mm_rate

adj_lower <- adjMixture(
  linked.data = hrscms_clean,
  m.rate = lower_mm_rate
)

adjusted_table_lower <- plctable(
  formula = ~ NRSHOM + mds_nrshom,
  adjustment = adj_lower
)

round(as.data.frame.matrix(adjusted_table_lower$phat), 4)
```

|     |      0 |      1 |
|:----|-------:|-------:|
| 0   | 0.9132 | 0.0194 |
| 1   | 0.0222 | 0.0452 |

Adjusted contingency table using lower mismatch rate {.table .table
style="width: max-content !important; max-width: none !important; margin-bottom: 20px;"}

These adjusted tables estimate the joint distribution of the two nursing
home residence indicators after accounting for possible linkage errors.

## **5. Evaluation**

We evaluate the adjusted contingency table for the HRSCMS dataset,
quantifying how well it recovers the joint distribution of nursing home
residence indicators across HRS and CMS records while accounting for
mismatches.

``` r

# Use the adjusted table from the main adjustment approach
phat_adj <- as.vector(t(adjusted_table_lower$phat)) 

# Use the exact-match table from the aggregated data
exactprobs <- t(ctab_exact_prop)

# MRAE
mrae <- mean(abs((phat_adj - exactprobs) / exactprobs))

# KLD
kld <- sum(exactprobs * log(exactprobs / phat_adj))

# GOF
gof <- n * sum((exactprobs - phat_adj)^2 / exactprobs)

# Association
adj_counts <- matrix(
  n * phat_adj,
  byrow = TRUE,
  nrow = 2,
  ncol = 2
)

c2 <- suppressWarnings(chisq.test(adj_counts))

# Kappa
p_o <- sum(diag(adjusted_table_lower$phat))
p_e <- sum(rowSums(adjusted_table_lower$phat) * colSums(adjusted_table_lower$phat))
k <- (p_o - p_e) / (1 - p_e)

# Combine results into table
data.frame(
  Metric = c("MRAE", "KLD", "GOF", "Chi-square", "Kappa"),
  Value  = round(c(mrae, kld, gof,
                   as.numeric(c2$statistic), 
                   k), 4)
)
```

| Metric     |    Value |
|:-----------|---------:|
| MRAE       |   0.1334 |
| KLD        |   0.0017 |
| GOF        |   1.0995 |
| Chi-square | 147.0581 |
| Kappa      |   0.6625 |

Evaluation metrics for the adjusted contingency table {.table .table
style="width: max-content !important; max-width: none !important; margin-bottom: 20px;"}

These metrics summarize the behavior of the `postlink` adjustment
procedure relative to the exact-match contingency table. MRAE quantifies
the average relative difference in joint cell probabilities, KLD
measures distributional divergence, GOF summarizes the overall
discrepancy between adjusted and exact-match probabilities, the
chi-square statistic measures association between the HRS and CMS
indicators, and Kappa evaluates agreement beyond chance.

## **References**

Centers for Medicare & Medicaid Services (CMS). 2026. [CMS
Data](https://www.cms.gov/), U.S. Department of Health & Human Services,
Washington, DC, US.

Health and Retirement Study (HRS). 2026. [HRS
Dataset](https://hrs.isr.umich.edu/), University of Michigan, Ann Arbor,
MI, US.

Slawski, Martin, Brady T. West, Priyanjali Bukke, Zhenbang Wang, Guoqing
Diao, and Emanuel Ben-David. 2024. “A General Framework for Regression
with Mismatched Data Based on Mixture Modelling.” *Journal of the Royal
Statistical Society Series A: Statistics in Society*, qnae083.
