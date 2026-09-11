version 17.0
clear all
set more off

/*
Purpose
-------
Reproduce Table 8 from analysis_sample.xlsx:
  - Pearson correlation
  - ICC(A,1): two-way mixed-effects, absolute-agreement, single-measure ICC
  - Lin's concordance correlation coefficient (CCC)
  - Bland-Altman mean bias and 95% limits of agreement
  - OLS calibration regression: Human = intercept + slope x GPT_H

Run from the folder containing analysis_sample.xlsx:
    do table8_agreement_statistics.do

Optional arguments:
    do table8_agreement_statistics.do "path/to/analysis_sample.xlsx" ///
        "path/to/Table8_agreement_statistics.rtf"
*/

args datafile outfile

if `"`datafile'"' == "" {
    local datafile "analysis_sample.xlsx"
}
if `"`outfile'"' == "" {
    local outfile "Table8_agreement_statistics.rtf"
}

capture confirm file "`datafile'"
if _rc {
    display as error "Input file not found: `datafile'"
    exit 601
}

/* Import the three data columns without relying on Stata's header conversion. */
import excel using "`datafile'", sheet("Final") cellrange(A2) clear
rename (A B C) (firm_year gpt_h human)
keep firm_year gpt_h human
drop if missing(firm_year) & missing(gpt_h) & missing(human)

replace firm_year = upper(strtrim(firm_year))

foreach v in gpt_h human {
    capture confirm numeric variable `v'
    if _rc {
        destring `v', replace ignore(",")
    }
}

/* Sector classification based on the firm codes in the supplied workbook. */
generate str12 sector = ""
replace sector = "Airlines"    if regexm(firm_year, "^(DELTA|LH|SA)")
replace sector = "Oil and gas" if regexm(firm_year, "^(BP|CHE|ENBRIDGE|SHELL)")

assert !missing(firm_year, gpt_h, human)
assert inlist(sector, "Airlines", "Oil and gas")
isid firm_year

/* Integrity checks for the supplied 2014-2023 analysis sample. */
quietly count
assert r(N) == 70
quietly count if sector == "Airlines"
assert r(N) == 30
quietly count if sector == "Oil and gas"
assert r(N) == 40

capture program drop agreement_stats
program define agreement_stats, rclass
    version 17.0
    syntax [if]
    marksample touse

    preserve
        keep if `touse'
        drop if missing(gpt_h, human)

        quietly count
        local n = r(N)
        if `n' < 3 {
            display as error "At least three complete observations are required."
            restore
            exit 2001
        }

        tempname R Cov Nval Rval ICCval CCCval Biasval Lowval Highval
        tempname Aval Bval MeanX MeanY VarX VarY Grand SSR SSC SSE MSR MSC MSE

        scalar `Nval' = `n'

        /* Pearson correlation. */
        quietly correlate gpt_h human
        matrix `R' = r(C)
        scalar `Rval' = `R'[1,2]

        /* Means and sample variances. */
        quietly summarize gpt_h, meanonly
        scalar `MeanX' = r(mean)
        scalar `VarX'  = r(Var)

        quietly summarize human, meanonly
        scalar `MeanY' = r(mean)
        scalar `VarY'  = r(Var)

        /*
        Lin's CCC. The sample covariance and sample variances reproduce the
        estimator and rounding used in the manuscript table.
        */
        quietly correlate gpt_h human, covariance
        matrix `Cov' = r(C)
        scalar `CCCval' = (2 * `Cov'[1,2]) / ///
            (`VarX' + `VarY' + (`MeanX' - `MeanY')^2)

        /*
        ICC(A,1), following the two-way absolute-agreement ANOVA formula:
        (MSR-MSE) / [MSR + (k-1)MSE + k(MSC-MSE)/n], with k=2 raters.
        */
        tempvar rowmean ssrow sse
        scalar `Grand' = (`MeanX' + `MeanY') / 2
        generate double `rowmean' = (gpt_h + human) / 2
        generate double `ssrow' = (`rowmean' - `Grand')^2
        quietly summarize `ssrow', meanonly
        scalar `SSR' = 2 * r(sum)

        scalar `SSC' = `Nval' * ///
            ((`MeanX' - `Grand')^2 + (`MeanY' - `Grand')^2)

        generate double `sse' = ///
            (gpt_h - `rowmean' - `MeanX' + `Grand')^2 + ///
            (human - `rowmean' - `MeanY' + `Grand')^2
        quietly summarize `sse', meanonly
        scalar `SSE' = r(sum)

        scalar `MSR' = `SSR' / (`Nval' - 1)
        scalar `MSC' = `SSC' / (2 - 1)
        scalar `MSE' = `SSE' / ((`Nval' - 1) * (2 - 1))
        scalar `ICCval' = (`MSR' - `MSE') / ///
            (`MSR' + `MSE' + (2 / `Nval') * (`MSC' - `MSE'))

        /* Bland-Altman bias and 95% limits of agreement. */
        tempvar difference
        generate double `difference' = gpt_h - human
        quietly summarize `difference'
        scalar `Biasval' = r(mean)
        scalar `Lowval'  = r(mean) - invnormal(.975) * r(sd)
        scalar `Highval' = r(mean) + invnormal(.975) * r(sd)

        /* Calibration regression. */
        quietly regress human gpt_h
        scalar `Aval' = _b[_cons]
        scalar `Bval' = _b[gpt_h]

    restore

    return scalar N         = `Nval'
    return scalar pearson   = `Rval'
    return scalar icca1     = `ICCval'
    return scalar ccc       = `CCCval'
    return scalar bias      = `Biasval'
    return scalar loa_low   = `Lowval'
    return scalar loa_high  = `Highval'
    return scalar intercept = `Aval'
    return scalar slope     = `Bval'
end

/* Calculate and store pooled and sector-specific results. */
tempname posth
tempfile table8data
postfile `posth' str12 sample int N double pearson icca1 ccc bias ///
    loa_low loa_high intercept slope using "`table8data'", replace

quietly agreement_stats
post `posth' ("Pooled") (r(N)) (r(pearson)) (r(icca1)) (r(ccc)) ///
    (r(bias)) (r(loa_low)) (r(loa_high)) (r(intercept)) (r(slope))

quietly agreement_stats if sector == "Airlines"
post `posth' ("Airlines") (r(N)) (r(pearson)) (r(icca1)) (r(ccc)) ///
    (r(bias)) (r(loa_low)) (r(loa_high)) (r(intercept)) (r(slope))

quietly agreement_stats if sector == "Oil and gas"
post `posth' ("Oil and gas") (r(N)) (r(pearson)) (r(icca1)) (r(ccc)) ///
    (r(bias)) (r(loa_low)) (r(loa_high)) (r(intercept)) (r(slope))

postclose `posth'
use "`table8data'", clear

format N %4.0f
format pearson icca1 ccc %5.3f
format bias loa_low loa_high intercept slope %7.2f

display as text _newline "Table 8 calculation results"
list sample N pearson icca1 ccc bias loa_low loa_high intercept slope, ///
    noobs abbreviate(18)

/*
Write a Word-compatible RTF table. RTF Unicode controls are used for the
curly apostrophe, en dash and multiplication sign.
*/
file open rtf using "`outfile'", write text replace
file write rtf "{\rtf1\ansi\ansicpg1252\deff0" _n
file write rtf "{\fonttbl{\f0 Times New Roman;}}" _n
file write rtf ///
    "\landscape\paperw15840\paperh12240\margl540\margr540" ///
    "\margt720\margb720\fs20" _n
file write rtf ///
    "\pard\sa40\b Table 8: Pooled and sector-specific agreement, bias, " ///
    "and calibration statistics for GPT_H and human scores\b0\par" _n

/* Header row: top and bottom rules, no vertical gridlines. */
file write rtf "\trowd\trgaph80\trleft0" ///
    "\clbrdrt\brdrs\brdrw10\clbrdrb\brdrs\brdrw10\cellx1300" ///
    "\clbrdrt\brdrs\brdrw10\clbrdrb\brdrs\brdrw10\cellx2000" ///
    "\clbrdrt\brdrs\brdrw10\clbrdrb\brdrs\brdrw10\cellx3300" ///
    "\clbrdrt\brdrs\brdrw10\clbrdrb\brdrs\brdrw10\cellx4500" ///
    "\clbrdrt\brdrs\brdrw10\clbrdrb\brdrs\brdrw10\cellx5700" ///
    "\clbrdrt\brdrs\brdrw10\clbrdrb\brdrs\brdrw10\cellx8600" ///
    "\clbrdrt\brdrs\brdrw10\clbrdrb\brdrs\brdrw10\cellx11000" ///
    "\clbrdrt\brdrs\brdrw10\clbrdrb\brdrs\brdrw10\cellx14760" _n

file write rtf ///
    "\pard\intbl\ql\b Sample\b0\cell" ///
    "\pard\intbl\qc\b N\b0\cell" ///
    "\pard\intbl\qc\b Pearson r\b0\cell" ///
    "\pard\intbl\qc\b ICC(A,1)\b0\cell" ///
    "\pard\intbl\qc\b Lin\u8217?s CCC\b0\cell" ///
    "\pard\intbl\qc\b Mean bias (GPT \u8211? Human)\b0\cell" ///
    "\pard\intbl\qc\b 95% limits of agreement\b0\cell" ///
    "\pard\intbl\ql\b Calibration regression\b0\cell\row" _n

forvalues i = 1/`=_N' {
    local s    = sample[`i']
    local ns   = strtrim(string(N[`i'], "%9.0f"))
    local ps   = strtrim(string(pearson[`i'], "%9.3f"))
    local is   = strtrim(string(icca1[`i'], "%9.3f"))
    local cs   = strtrim(string(ccc[`i'], "%9.3f"))
    local bs   = strtrim(string(bias[`i'], "%9.2f"))
    local los  = "[" + strtrim(string(loa_low[`i'], "%9.2f")) + ///
        ", " + strtrim(string(loa_high[`i'], "%9.2f")) + "]"
    local cal  = "Human = " + strtrim(string(intercept[`i'], "%9.2f")) + ///
        " + " + strtrim(string(slope[`i'], "%9.2f")) + " \u215? GPT"

    local bottom ""
    if `i' == _N {
        local bottom "\clbrdrb\brdrs\brdrw10"
    }

    file write rtf "\trowd\trgaph80\trleft0" ///
        "`bottom'\cellx1300"  "`bottom'\cellx2000" ///
        "`bottom'\cellx3300"  "`bottom'\cellx4500" ///
        "`bottom'\cellx5700"  "`bottom'\cellx8600" ///
        "`bottom'\cellx11000" "`bottom'\cellx14760" _n

    file write rtf ///
        "\pard\intbl\ql `s'\cell" ///
        "\pard\intbl\qc `ns'\cell" ///
        "\pard\intbl\qc `ps'\cell" ///
        "\pard\intbl\qc `is'\cell" ///
        "\pard\intbl\qc `cs'\cell" ///
        "\pard\intbl\qc `bs'\cell" ///
        "\pard\intbl\qc `los'\cell" ///
        "\pard\intbl\ql `cal'\cell\row" _n
}

file write rtf "}" _n
file close rtf

display as result _newline "Created: `outfile'"

