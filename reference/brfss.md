# Behavioral Risk Factor Surveillance System (BRFSS) 2013 Subsample

The `brfss` data set contains a randomly selected subsample of 2,000
respondents from the 2013 Behavioral Risk Factor Surveillance System
(BRFSS) survey. The BRFSS is an ongoing, state-based telephone survey
conducted by the U.S. Centers for Disease Control and Prevention (CDC)
that collects information on health-related risk behaviors, chronic
health conditions, and the use of preventive health services from
non-institutionalized adults.

## Usage

``` r
data(brfss)
```

## Format

A data frame with 2,000 rows and 8 variables.

## Details

In the context of post-linkage data analysis (PLDA), this data set is
partitioned into two separate files to simulate a secondary analysis
scenario where disparate data sources must be linked using
quasi-identifiers (e.g., state, month of interview, sex, and marital
status) prior to performing regression analyses (e.g., predicting weight
based on height, physical health, and mental health).

Pre-processing steps have been applied to this subsample to remove
missing values ("Refused" or `NA`) and to filter out extreme outliers in
physical measurements.

- `X`: State of residence (categorical).

- `imonth`: Month of interview (categorical).

- `Weight`: Weight of the respondent in pounds. Records with a weight
  less than 50 lbs or greater than 650 lbs were excluded.

- `Height`: Height of the respondent, formulated as (feet \* 100 +
  inches). Records with a height less than 54 inches or greater than 84
  inches were excluded.

- `Physhlth`: Number of days with poor physical health. Records marked
  "Refused" or NA were excluded.

- `Menthlth`: Number of days with poor mental health. Records marked
  "Refused" or NA were excluded.

- `Exerany`: Indicator of whether the respondent engaged in physical
  activity. Records marked "Refused" or NA were excluded.

- `m.rate`: Block-wise (by interview month) mismatch rate.

## References

Centers for Disease Control and Prevention (CDC) (2013). "Behavioral
Risk Factor Surveillance System Survey Data." U.S. Department of Health
and Human Services, Centers for Disease Control and Prevention.
Available at <https://www.cdc.gov/brfss/annual_data/annual_2013.html>.
