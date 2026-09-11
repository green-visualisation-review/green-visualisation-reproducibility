version 17.0
clear all
set more off

/*
Purpose
-------
Reproduce Tables 9 and 10 from the supplied analysis sample.

Table 9 reports pooled and sector-specific paired-sample comparisons.
Table 10 reports the main agreement statistics, exclusion of BP,
leave-one-firm-out and leave-one-year-out ranges, and a 5% winsorised
sensitivity analysis of the paired differences.

The code uses only official Stata commands and functions. No community-
contributed packages are required.

Default run, with this do-file and the workbook in the same folder:
    do tables9_10_paired_robustness.do

Optional arguments:
    do tables9_10_paired_robustness.do "path/to/analysis_sample(1).xlsx" ///
        "path/to/Tables9_10_paired_robustness.rtf"
*/

args datafile outfile

if `"`datafile'"' == "" {
    local datafile "analysis_sample(1).xlsx"
}
if `"`outfile'"' == "" {
    local outfile "Tables9_10_paired_robustness.rtf"
}

capture confirm file "`datafile'"
if _rc {
    display as error "Input file not found: `datafile'"
    exit 601
}

/* Import the three columns without relying on Stata's header conversion. */
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

/* Recover firm, year, and sector from the firm-year identifier. */
generate int year = real(regexs(1)) ///
    if regexm(firm_year, "([0-9][0-9][0-9][0-9])$")
generate str12 firm = regexr(firm_year, ///
    "[ ]*[0-9][0-9][0-9][0-9]$", "")
replace firm = strtrim(firm)

generate str12 sector = "Oil and gas"
replace sector = "Airlines" if inlist(firm, "DELTA", "LH", "SA")

assert firm_year != "" & firm != "" & sector != ""
assert !missing(year, gpt_h, human)
assert inlist(firm, "BP", "CHE", "DELTA", "ENBRIDGE", "LH", "SA", "SHELL")
isid firm year

/* Integrity checks for the supplied 2014-2023 sample. */
quietly count
assert r(N) == 70
quietly count if sector == "Airlines"
assert r(N) == 30
quietly count if sector == "Oil and gas"
assert r(N) == 40

/* ---------------------------------------------------------------------- */
/* Paired-sample statistics, including the noncentral-t CI for Cohen's dz. */
/* ---------------------------------------------------------------------- */
capture program drop paired_stats
program define paired_stats, rclass
    version 17.0
    syntax [if]
    marksample touse

    preserve
        keep if `touse'
        drop if missing(gpt_h, human)

        tempvar difference
        generate double `difference' = gpt_h - human

        quietly summarize `difference'

        tempname Nval Meanval SDval DFval Tval Pval Critval
        tempname Lowval Highval DZval DZLowval DZHighval

        scalar `Nval' = r(N)
        scalar `Meanval' = r(mean)
        scalar `SDval' = r(sd)
        scalar `DFval' = `Nval' - 1
        scalar `Tval' = `Meanval' / (`SDval' / sqrt(`Nval'))
        scalar `Pval' = 2 * ttail(`DFval', abs(`Tval'))
        scalar `Critval' = invttail(`DFval', .025)
        scalar `Lowval' = `Meanval' - ///
            `Critval' * `SDval' / sqrt(`Nval')
        scalar `Highval' = `Meanval' + ///
            `Critval' * `SDval' / sqrt(`Nval')

        /* Cohen's dz and its 95% CI from noncentral-t inversion. */
        scalar `DZval' = `Meanval' / `SDval'
        scalar `DZLowval' = npnt(`DFval', `Tval', .975) / sqrt(`Nval')
        scalar `DZHighval' = npnt(`DFval', `Tval', .025) / sqrt(`Nval')

    restore

    return scalar N       = `Nval'
    return scalar mean    = `Meanval'
    return scalar sd      = `SDval'
    return scalar df      = `DFval'
    return scalar t       = `Tval'
    return scalar p       = `Pval'
    return scalar ci_low  = `Lowval'
    return scalar ci_high = `Highval'
    return scalar dz      = `DZval'
    return scalar dz_low  = `DZLowval'
    return scalar dz_high = `DZHighval'
end

/* ---------------------------------------------------------------------- */
/* Pearson r, ICC(A,1), and Lin's CCC for Table 10.                       */
/* ---------------------------------------------------------------------- */
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

        tempname R Cov Nval Rval ICCval CCCval
        tempname MeanX MeanY VarX VarY Grand SSR SSC SSE MSR MSC MSE

        scalar `Nval' = `n'

        quietly correlate gpt_h human
        matrix `R' = r(C)
        scalar `Rval' = `R'[1,2]

        quietly summarize gpt_h, meanonly
        scalar `MeanX' = r(mean)
        scalar `VarX' = r(Var)

        quietly summarize human, meanonly
        scalar `MeanY' = r(mean)
        scalar `VarY' = r(Var)

        /* Lin's CCC, using sample covariance and sample variances. */
        quietly correlate gpt_h human, covariance
        matrix `Cov' = r(C)
        scalar `CCCval' = (2 * `Cov'[1,2]) / ///
            (`VarX' + `VarY' + (`MeanX' - `MeanY')^2)

        /* ICC(A,1): two-way absolute-agreement, single-measure ICC. */
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

    restore

    return scalar N       = `Nval'
    return scalar pearson = `Rval'
    return scalar icca1   = `ICCval'
    return scalar ccc     = `CCCval'
end

/* ---------------------------------------------------------------------- */
/* Table 9 calculations.                                                  */
/* ---------------------------------------------------------------------- */
tempname p9
tempfile table9data
postfile `p9' str36 comparison int N double mean_diff sd_diff t df p ///
    ci_low ci_high dz dz_low dz_high using "`table9data'", replace

quietly paired_stats
post `p9' ("GPT_H - Human (pooled)") (r(N)) (r(mean)) (r(sd)) ///
    (r(t)) (r(df)) (r(p)) (r(ci_low)) (r(ci_high)) ///
    (r(dz)) (r(dz_low)) (r(dz_high))

quietly paired_stats if sector == "Airlines"
post `p9' ("GPT_H - Human (airlines)") (r(N)) (r(mean)) (r(sd)) ///
    (r(t)) (r(df)) (r(p)) (r(ci_low)) (r(ci_high)) ///
    (r(dz)) (r(dz_low)) (r(dz_high))

quietly paired_stats if sector == "Oil and gas"
post `p9' ("GPT_H - Human (oil and gas)") (r(N)) (r(mean)) (r(sd)) ///
    (r(t)) (r(df)) (r(p)) (r(ci_low)) (r(ci_high)) ///
    (r(dz)) (r(dz_low)) (r(dz_high))

postclose `p9'

/* ---------------------------------------------------------------------- */
/* Table 10 calculations.                                                 */
/* ---------------------------------------------------------------------- */

/* Main pooled result. */
quietly agreement_stats
local main_r   = r(pearson)
local main_icc = r(icca1)
local main_ccc = r(ccc)

quietly paired_stats
local main_bias = r(mean)
local main_t    = r(t)
local main_df   = r(df)
local main_lo   = r(ci_low)
local main_hi   = r(ci_high)

/* Excluding BP. */
quietly agreement_stats if firm != "BP"
local bp_r   = r(pearson)
local bp_icc = r(icca1)
local bp_ccc = r(ccc)

quietly paired_stats if firm != "BP"
local bp_bias = r(mean)
local bp_t    = r(t)
local bp_df   = r(df)
local bp_lo   = r(ci_low)
local bp_hi   = r(ci_high)

/* Leave-one-firm-out agreement ranges. */
tempname plofo
tempfile lofodata
postfile `plofo' str12 omitted_firm double pearson icca1 ccc ///
    using "`lofodata'", replace

quietly levelsof firm, local(firms)
foreach f of local firms {
    quietly agreement_stats if firm != "`f'"
    post `plofo' ("`f'") (r(pearson)) (r(icca1)) (r(ccc))
}
postclose `plofo'

preserve
    use "`lofodata'", clear
    quietly summarize pearson, meanonly
    local lofo_r_lo = r(min)
    local lofo_r_hi = r(max)
    quietly summarize icca1, meanonly
    local lofo_icc_lo = r(min)
    local lofo_icc_hi = r(max)
    quietly summarize ccc, meanonly
    local lofo_ccc_lo = r(min)
    local lofo_ccc_hi = r(max)
restore

/* Leave-one-year-out agreement ranges. */
tempname ployo
tempfile loyodata
postfile `ployo' int omitted_year double pearson icca1 ccc ///
    using "`loyodata'", replace

quietly levelsof year, local(years)
foreach y of local years {
    quietly agreement_stats if year != `y'
    post `ployo' (`y') (r(pearson)) (r(icca1)) (r(ccc))
}
postclose `ployo'

preserve
    use "`loyodata'", clear
    quietly summarize pearson, meanonly
    local loyo_r_lo = r(min)
    local loyo_r_hi = r(max)
    quietly summarize icca1, meanonly
    local loyo_icc_lo = r(min)
    local loyo_icc_hi = r(max)
    quietly summarize ccc, meanonly
    local loyo_ccc_lo = r(min)
    local loyo_ccc_hi = r(max)
restore

/* 5% winsorisation of the paired differences, not of each score variable. */
preserve
    generate double difference = gpt_h - human
    quietly _pctile difference, p(5 95)
    local win_lo_cut = r(r1)
    local win_hi_cut = r(r2)

    generate double difference_w = ///
        min(max(difference, `win_lo_cut'), `win_hi_cut')
    quietly summarize difference_w

    local win_n    = r(N)
    local win_bias = r(mean)
    local win_sd   = r(sd)
    local win_df   = `win_n' - 1
    local win_t    = `win_bias' / (`win_sd' / sqrt(`win_n'))
    local win_crit = invttail(`win_df', .025)
    local win_ci_lo = `win_bias' - ///
        `win_crit' * `win_sd' / sqrt(`win_n')
    local win_ci_hi = `win_bias' + ///
        `win_crit' * `win_sd' / sqrt(`win_n')
restore

/* Display compact verification results in Stata's Results window. */
preserve
    use "`table9data'", clear
    format N df %4.0f
    format mean_diff sd_diff t ci_low ci_high dz dz_low dz_high %7.3f
    format p %9.6f
    display as text _newline "Table 9 calculation results"
    list, noobs abbreviate(18)
restore

display as text _newline "Table 10 range checks"
display as text "LOFO Pearson r: " %6.3f `lofo_r_lo' " to " %6.3f `lofo_r_hi'
display as text "LOFO ICC(A,1): " %6.3f `lofo_icc_lo' " to " %6.3f `lofo_icc_hi'
display as text "LOFO Lin's CCC: " %6.3f `lofo_ccc_lo' " to " %6.3f `lofo_ccc_hi'
display as text "LOYO Pearson r: " %6.3f `loyo_r_lo' " to " %6.3f `loyo_r_hi'
display as text "LOYO ICC(A,1): " %6.3f `loyo_icc_lo' " to " %6.3f `loyo_icc_hi'
display as text "LOYO Lin's CCC: " %6.3f `loyo_ccc_lo' " to " %6.3f `loyo_ccc_hi'
display as text "Winsor cutoffs: " %6.2f `win_lo_cut' " and " %6.2f `win_hi_cut'

/* ---------------------------------------------------------------------- */
/* Export both tables to a Word-compatible RTF file.                      */
/* ---------------------------------------------------------------------- */
file open rtf using "`outfile'", write text replace
file write rtf "{\rtf1\ansi\ansicpg1252\deff0" _n
file write rtf "{\fonttbl{\f0 Times New Roman;}}" _n
file write rtf ///
    "\landscape\paperw15840\paperh12240\margl540\margr540" ///
    "\margt540\margb540\fs20" _n

/* ------------------------------- Table 9 ------------------------------ */
file write rtf ///
    "\pard\sa40\b Table 9: Pooled and sector-specific paired-sample " ///
    "comparisons of GPT_H and human scores\b0\par" _n

file write rtf "\trowd\trgaph70\trleft0" ///
    "\clbrdrt\brdrs\brdrw10\clbrdrb\brdrs\brdrw10\cellx3000" ///
    "\clbrdrt\brdrs\brdrw10\clbrdrb\brdrs\brdrw10\cellx3600" ///
    "\clbrdrt\brdrs\brdrw10\clbrdrb\brdrs\brdrw10\cellx5100" ///
    "\clbrdrt\brdrs\brdrw10\clbrdrb\brdrs\brdrw10\cellx6600" ///
    "\clbrdrt\brdrs\brdrw10\clbrdrb\brdrs\brdrw10\cellx7400" ///
    "\clbrdrt\brdrs\brdrw10\clbrdrb\brdrs\brdrw10\cellx8000" ///
    "\clbrdrt\brdrs\brdrw10\clbrdrb\brdrs\brdrw10\cellx8750" ///
    "\clbrdrt\brdrs\brdrw10\clbrdrb\brdrs\brdrw10\cellx11000" ///
    "\clbrdrt\brdrs\brdrw10\clbrdrb\brdrs\brdrw10\cellx12150" ///
    "\clbrdrt\brdrs\brdrw10\clbrdrb\brdrs\brdrw10\cellx14760" _n

file write rtf ///
    "\pard\intbl\ql\b Comparison\b0\cell" ///
    "\pard\intbl\qc\b Obs\b0\cell" ///
    "\pard\intbl\qc\b Mean difference\b0\cell" ///
    "\pard\intbl\qc\b SD of difference\b0\cell" ///
    "\pard\intbl\qc\b t\b0\cell" ///
    "\pard\intbl\qc\b df\b0\cell" ///
    "\pard\intbl\qc\b p\b0\cell" ///
    "\pard\intbl\qc\b Two-sided 95% CI for\line mean difference\b0\cell" ///
    "\pard\intbl\qc\b Cohen\u8217?s dz\b0\cell" ///
    "\pard\intbl\qc\b 95% CI for dz\b0\cell\row" _n

preserve
    use "`table9data'", clear

    forvalues i = 1/`=_N' {
        local comp = comparison[`i']
        local comp = subinstr("`comp'", " - ", " \u8211? ", .)
        local ns   = strtrim(string(N[`i'], "%9.0f"))
        local ms   = strtrim(string(mean_diff[`i'], "%9.2f"))
        local sds  = strtrim(string(sd_diff[`i'], "%9.2f"))
        local ts   = strtrim(string(t[`i'], "%9.2f"))
        local dfs  = strtrim(string(df[`i'], "%9.0f"))
        if p[`i'] < .001 {
            local ps "<0.001"
        }
        else {
            local ps = strtrim(string(p[`i'], "%9.3f"))
        }
        local cis = "[" + strtrim(string(ci_low[`i'], "%9.2f")) + ///
            ", " + strtrim(string(ci_high[`i'], "%9.2f")) + "]"
        local dzs = strtrim(string(dz[`i'], "%9.2f"))
        local dzcis = "[" + strtrim(string(dz_low[`i'], "%9.2f")) + ///
            ", " + strtrim(string(dz_high[`i'], "%9.2f")) + "]"

        local bottom ""
        if `i' == _N {
            local bottom "\clbrdrb\brdrs\brdrw10"
        }

        file write rtf "\trowd\trgaph70\trleft0" ///
            "`bottom'\cellx3000"  "`bottom'\cellx3600" ///
            "`bottom'\cellx5100"  "`bottom'\cellx6600" ///
            "`bottom'\cellx7400"  "`bottom'\cellx8000" ///
            "`bottom'\cellx8750"  "`bottom'\cellx11000" ///
            "`bottom'\cellx12150" "`bottom'\cellx14760" _n

        file write rtf ///
            "\pard\intbl\ql `comp'\cell" ///
            "\pard\intbl\qc `ns'\cell" ///
            "\pard\intbl\qc `ms'\cell" ///
            "\pard\intbl\qc `sds'\cell" ///
            "\pard\intbl\qc `ts'\cell" ///
            "\pard\intbl\qc `dfs'\cell" ///
            "\pard\intbl\qc `ps'\cell" ///
            "\pard\intbl\qc `cis'\cell" ///
            "\pard\intbl\qc `dzs'\cell" ///
            "\pard\intbl\qc `dzcis'\cell\row" _n
    }
restore

/* Prepare the formatted Table 10 strings. */
local main_r_s   = strtrim(string(`main_r', "%9.3f"))
local main_icc_s = strtrim(string(`main_icc', "%9.3f"))
local main_ccc_s = strtrim(string(`main_ccc', "%9.3f"))
local main_bias_s = strtrim(string(`main_bias', "%9.2f"))
local main_tdf_s = strtrim(string(`main_t', "%9.2f")) + ///
    " (" + strtrim(string(`main_df', "%9.0f")) + ")"
local main_ci_s = "[" + strtrim(string(`main_lo', "%9.2f")) + ///
    ", " + strtrim(string(`main_hi', "%9.2f")) + "]"

local bp_r_s   = strtrim(string(`bp_r', "%9.3f"))
local bp_icc_s = strtrim(string(`bp_icc', "%9.3f"))
local bp_ccc_s = strtrim(string(`bp_ccc', "%9.3f"))
local bp_bias_s = strtrim(string(`bp_bias', "%9.2f"))
local bp_tdf_s = strtrim(string(`bp_t', "%9.2f")) + ///
    " (" + strtrim(string(`bp_df', "%9.0f")) + ")"
local bp_ci_s = "[" + strtrim(string(`bp_lo', "%9.2f")) + ///
    ", " + strtrim(string(`bp_hi', "%9.2f")) + "]"

local lofo_r_s = strtrim(string(`lofo_r_lo', "%9.3f")) + ///
    "\u8211?" + strtrim(string(`lofo_r_hi', "%9.3f"))
local lofo_icc_s = strtrim(string(`lofo_icc_lo', "%9.3f")) + ///
    "\u8211?" + strtrim(string(`lofo_icc_hi', "%9.3f"))
local lofo_ccc_s = strtrim(string(`lofo_ccc_lo', "%9.3f")) + ///
    "\u8211?" + strtrim(string(`lofo_ccc_hi', "%9.3f"))

local loyo_r_s = strtrim(string(`loyo_r_lo', "%9.3f")) + ///
    "\u8211?" + strtrim(string(`loyo_r_hi', "%9.3f"))
local loyo_icc_s = strtrim(string(`loyo_icc_lo', "%9.3f")) + ///
    "\u8211?" + strtrim(string(`loyo_icc_hi', "%9.3f"))
local loyo_ccc_s = strtrim(string(`loyo_ccc_lo', "%9.3f")) + ///
    "\u8211?" + strtrim(string(`loyo_ccc_hi', "%9.3f"))

local win_bias_s = strtrim(string(`win_bias', "%9.2f"))
local win_tdf_s = strtrim(string(`win_t', "%9.2f")) + ///
    " (" + strtrim(string(`win_df', "%9.0f")) + ")"
local win_ci_s = "[" + strtrim(string(`win_ci_lo', "%9.2f")) + ///
    ", " + strtrim(string(`win_ci_hi', "%9.2f")) + "]"
local dash "\u8212?"

/* ------------------------------ Table 10 ------------------------------ */
file write rtf "\pard\sa180\par" _n
file write rtf ///
    "\pard\sa40\b Table 10: Robustness checks for agreement between " ///
    "GPT_H and human scores\b0\par" _n

file write rtf "\trowd\trgaph70\trleft0" ///
    "\clbrdrt\brdrs\brdrw10\clbrdrb\brdrs\brdrw10\cellx3300" ///
    "\clbrdrt\brdrs\brdrw10\clbrdrb\brdrs\brdrw10\cellx4800" ///
    "\clbrdrt\brdrs\brdrw10\clbrdrb\brdrs\brdrw10\cellx6200" ///
    "\clbrdrt\brdrs\brdrw10\clbrdrb\brdrs\brdrw10\cellx7600" ///
    "\clbrdrt\brdrs\brdrw10\clbrdrb\brdrs\brdrw10\cellx10500" ///
    "\clbrdrt\brdrs\brdrw10\clbrdrb\brdrs\brdrw10\cellx12300" ///
    "\clbrdrt\brdrs\brdrw10\clbrdrb\brdrs\brdrw10\cellx14760" _n

file write rtf ///
    "\pard\intbl\ql\b Robustness specification\b0\cell" ///
    "\pard\intbl\qc\b Pearson r\b0\cell" ///
    "\pard\intbl\qc\b ICC(A,1)\b0\cell" ///
    "\pard\intbl\qc\b Lin\u8217?s CCC\b0\cell" ///
    "\pard\intbl\qc\b Mean bias (GPT_H \u8211? Human)\b0\cell" ///
    "\pard\intbl\qc\b t(df)\b0\cell" ///
    "\pard\intbl\qc\b 95% CI\b0\cell\row" _n

/* Main pooled result. */
file write rtf "\trowd\trgaph70\trleft0\cellx3300\cellx4800" ///
    "\cellx6200\cellx7600\cellx10500\cellx12300\cellx14760" _n
file write rtf ///
    "\pard\intbl\ql Main pooled result\cell" ///
    "\pard\intbl\qc `main_r_s'\cell" ///
    "\pard\intbl\qc `main_icc_s'\cell" ///
    "\pard\intbl\qc `main_ccc_s'\cell" ///
    "\pard\intbl\qc `main_bias_s'\cell" ///
    "\pard\intbl\qc `main_tdf_s'\cell" ///
    "\pard\intbl\qc `main_ci_s'\cell\row" _n

/* Excluding BP. */
file write rtf "\trowd\trgaph70\trleft0\cellx3300\cellx4800" ///
    "\cellx6200\cellx7600\cellx10500\cellx12300\cellx14760" _n
file write rtf ///
    "\pard\intbl\ql Excluding BP\cell" ///
    "\pard\intbl\qc `bp_r_s'\cell" ///
    "\pard\intbl\qc `bp_icc_s'\cell" ///
    "\pard\intbl\qc `bp_ccc_s'\cell" ///
    "\pard\intbl\qc `bp_bias_s'\cell" ///
    "\pard\intbl\qc `bp_tdf_s'\cell" ///
    "\pard\intbl\qc `bp_ci_s'\cell\row" _n

/* Leave-one-firm-out range. */
file write rtf "\trowd\trgaph70\trleft0\cellx3300\cellx4800" ///
    "\cellx6200\cellx7600\cellx10500\cellx12300\cellx14760" _n
file write rtf ///
    "\pard\intbl\ql Leave-one-firm-out range\cell" ///
    "\pard\intbl\qc `lofo_r_s'\cell" ///
    "\pard\intbl\qc `lofo_icc_s'\cell" ///
    "\pard\intbl\qc `lofo_ccc_s'\cell" ///
    "\pard\intbl\qc `dash'\cell" ///
    "\pard\intbl\qc `dash'\cell" ///
    "\pard\intbl\qc `dash'\cell\row" _n

/* Leave-one-year-out range. */
file write rtf "\trowd\trgaph70\trleft0\cellx3300\cellx4800" ///
    "\cellx6200\cellx7600\cellx10500\cellx12300\cellx14760" _n
file write rtf ///
    "\pard\intbl\ql Leave-one-year-out range\cell" ///
    "\pard\intbl\qc `loyo_r_s'\cell" ///
    "\pard\intbl\qc `loyo_icc_s'\cell" ///
    "\pard\intbl\qc `loyo_ccc_s'\cell" ///
    "\pard\intbl\qc `dash'\cell" ///
    "\pard\intbl\qc `dash'\cell" ///
    "\pard\intbl\qc `dash'\cell\row" _n

/* 5% winsorised sensitivity, with bottom rule. */
file write rtf "\trowd\trgaph70\trleft0" ///
    "\clbrdrb\brdrs\brdrw10\cellx3300" ///
    "\clbrdrb\brdrs\brdrw10\cellx4800" ///
    "\clbrdrb\brdrs\brdrw10\cellx6200" ///
    "\clbrdrb\brdrs\brdrw10\cellx7600" ///
    "\clbrdrb\brdrs\brdrw10\cellx10500" ///
    "\clbrdrb\brdrs\brdrw10\cellx12300" ///
    "\clbrdrb\brdrs\brdrw10\cellx14760" _n
file write rtf ///
    "\pard\intbl\ql 5% winsorized sensitivity\cell" ///
    "\pard\intbl\qc `dash'\cell" ///
    "\pard\intbl\qc `dash'\cell" ///
    "\pard\intbl\qc `dash'\cell" ///
    "\pard\intbl\qc `win_bias_s'\cell" ///
    "\pard\intbl\qc `win_tdf_s'\cell" ///
    "\pard\intbl\qc `win_ci_s'\cell\row" _n

file write rtf ///
    "\pard\sa0\i Note: detailed results can be found in Appendix 4.\i0\par" _n
file write rtf "}" _n
file close rtf

display as result _newline "Created: `outfile'"
