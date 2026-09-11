version 17.0
clear all
set more off

/*
Purpose
-------
Reproduce manuscript Table 11:
  "Nested fixed-effects regressions of blinded human rubric scores."

Model 1:
  Human = DeltaGreen + firm fixed effects + year fixed effects

Model 2:
  Human = DeltaGreen + GPT_H + firm fixed effects + year fixed effects

The displayed coefficient standard errors are conventional OLS standard
errors. Inference for GPT_H is supplemented with a null-imposed, firm-level
wild-cluster bootstrap using Webb six-point weights, symmetric two-sided
p-values, 99,999 replications, and seed 20260729.

REQUIRED INPUT
--------------
This script requires the canonical analysis workbook prepared from:
  1. scoresforanalysis_sample.xlsx (GPT_H and Human); and
  2. Appendix 3 pixel results (sustainability, baseline, and adjusted
     green_coverage_all_mean values).

The workbook contains 70 observations on sheet "Final" and uses the
machine-readable variables firm_year, gpt_h, human,
green_sustainability_share, green_baseline_share, delta_green_share,
delta_green_pp, firm, and year.

Repository setup
----------------
Place the supplied canonical workbook at:
    data/analysis/analysis_sample.xlsx

Install the required community-contributed bootstrap command once:
    ssc install boottest

The calculation log records the installed boottest build reported by
"which boottest". Preserve that version information with the repository.

Default run from the repository root:
    do "code/stata/06_table11_nested_fixed_effects.do"

Optional arguments:
    do "code/stata/06_table11_nested_fixed_effects.do" ///
        "path/to/analysis_sample.xlsx" ///
        "path/to/Table11_nested_fixed_effects.rtf" ///
        "path/to/Table11_nested_fixed_effects.log"

The script computes all table entries from the data. It does not hard-code
the reported coefficients or bootstrap results.
*/

args datafile outfile logfile

if `"`datafile'"' == "" {
    local datafile "data/analysis/analysis_sample.xlsx"
}
if `"`outfile'"' == "" {
    local outfile "outputs/tables/Table11_nested_fixed_effects.rtf"
}
if `"`logfile'"' == "" {
    local logfile "outputs/logs/Table11_nested_fixed_effects.log"
}

capture mkdir "outputs"
capture mkdir "outputs/tables"
capture mkdir "outputs/logs"

capture confirm file "`datafile'"
if _rc {
    display as error "Required canonical analysis file not found: `datafile'"
    display as error "Place analysis_sample.xlsx in data/analysis/."
    exit 601
}

capture log close table11log
log using "`logfile'", text replace name(table11log)

display as text "Input:  `datafile'"
display as text "Output: `outfile'"
display as text "Stata:  " c(stata_version)

/* ---------------------------------------------------------------------- */
/* 1. Import and validate the 70-observation held-out panel.              */
/* ---------------------------------------------------------------------- */

capture noisily import excel using "`datafile'", ///
    sheet("Final") firstrow case(lower) clear
if _rc {
    display as error "Could not import sheet 'Final' from the canonical workbook."
    log close table11log
    exit 198
}

local required firm_year gpt_h human green_sustainability_share ///
    green_baseline_share delta_green_share delta_green_pp firm year
foreach v of local required {
    capture confirm variable `v'
    if _rc {
        display as error "Required variable missing from analysis_sample.xlsx: `v'"
        log close table11log
        exit 111
    }
}

keep `required'

foreach v in gpt_h human green_sustainability_share ///
    green_baseline_share delta_green_share delta_green_pp year {
    capture confirm numeric variable `v'
    if _rc {
        capture noisily destring `v', replace ignore(",")
        if _rc {
            display as error "Variable `v' could not be converted to numeric form."
            log close table11log
            exit 109
        }
    }
}

replace firm_year = upper(strtrim(itrim(firm_year)))
replace firm = subinstr(strtrim(firm), " ", "", .)
replace firm = upper(firm)

/* Drop only fully blank formatted rows; reject every partial record below. */
drop if firm_year == "" & firm == "" & missing(year) & missing(gpt_h) & ///
    missing(human) & missing(green_sustainability_share) & ///
    missing(green_baseline_share) & missing(delta_green_share) & ///
    missing(delta_green_pp)

assert firm_year != ""
assert firm != ""
assert !missing(year, gpt_h, human, green_sustainability_share, ///
    green_baseline_share, delta_green_share, delta_green_pp)
assert inlist(firm, "BP", "CHE", "DELTA", "ENBRIDGE", "LH", "SA", "SHELL")
assert year == floor(year) & inrange(year, 2014, 2023)
assert inrange(gpt_h, 0, 100) & inrange(human, 0, 100)

/* Verify the firm-year text against the separately supplied key columns. */
generate str20 expected_key = firm + string(year, "%4.0f")
generate str20 observed_key = subinstr(firm_year, " ", "", .)
assert observed_key == expected_key
drop expected_key observed_key

/* Verify both the DeltaGreen construction and its percentage-point scale. */
assert abs(delta_green_share - ///
    (green_sustainability_share - green_baseline_share)) < 1e-10
assert abs(delta_green_pp - 100 * delta_green_share) < 1e-9

drop delta_green_pp
generate double delta_green_pp = 100 * delta_green_share

egen firm_id = group(firm), label

label variable delta_green_pp "DeltaGreen (percentage points)"
label variable gpt_h "GPT_H rubric score"
label variable human "Blinded human rubric score"

isid firm_id year

quietly count
assert r(N) == 70

quietly levelsof firm_id, local(firm_levels)
local n_firms : word count `firm_levels'
assert `n_firms' == 7

quietly summarize year, meanonly
assert r(min) == 2014
assert r(max) == 2023

bysort firm_id: assert _N == 10
bysort year: assert _N == 7

/* Fingerprint the exact DeltaGreen series supporting the manuscript. */
quietly summarize delta_green_pp
if abs(r(mean) - 1.129758729) > 1e-6 | ///
   abs(r(sd)   - 1.984169712) > 1e-6 | ///
   abs(r(min)  + 2.447213458) > 1e-6 | ///
   abs(r(max)  - 6.856712486) > 1e-6 {
    display as error "DeltaGreen does not match the frozen Table 11 analysis sample."
    log close table11log
    exit 459
}

/* ---------------------------------------------------------------------- */
/* 2. Model 1: DeltaGreen plus firm and year fixed effects.               */
/* ---------------------------------------------------------------------- */

quietly regress human delta_green_pp i.firm_id i.year
estimates store Table11_Model1

scalar m1_b_delta  = _b[delta_green_pp]
scalar m1_se_delta = _se[delta_green_pp]
scalar m1_t_delta  = m1_b_delta / m1_se_delta
scalar m1_p_delta  = 2 * ttail(e(df_r), abs(m1_t_delta))
scalar m1_n        = e(N)
scalar m1_r2       = e(r2)
scalar m1_r2a      = e(r2_a)
scalar m1_dfr      = e(df_r)

/* ---------------------------------------------------------------------- */
/* 3. Model 2: add GPT_H.                                                 */
/* ---------------------------------------------------------------------- */

quietly regress human delta_green_pp gpt_h i.firm_id i.year
estimates store Table11_Model2

scalar m2_b_delta  = _b[delta_green_pp]
scalar m2_se_delta = _se[delta_green_pp]
scalar m2_t_delta  = m2_b_delta / m2_se_delta
scalar m2_p_delta  = 2 * ttail(e(df_r), abs(m2_t_delta))

scalar m2_b_gpt    = _b[gpt_h]
scalar m2_se_gpt   = _se[gpt_h]
scalar m2_t_gpt    = m2_b_gpt / m2_se_gpt
scalar m2_p_gpt    = 2 * ttail(e(df_r), abs(m2_t_gpt))

scalar m2_n        = e(N)
scalar m2_r2       = e(r2)
scalar m2_r2a      = e(r2_a)
scalar m2_dfr      = e(df_r)
scalar delta_r2a   = m2_r2a - m1_r2a

quietly test gpt_h
scalar partial_f     = r(F)
scalar partial_p     = r(p)
scalar partial_dfnum = r(df)
scalar partial_dfden = r(df_r)

/* VIFs for the two continuous regressors, conditional on both FE sets. */
quietly regress delta_green_pp gpt_h i.firm_id i.year
scalar vif_delta = 1 / (1 - e(r2))

quietly regress gpt_h delta_green_pp i.firm_id i.year
scalar vif_gpt = 1 / (1 - e(r2))

/* Confirm the conventional results before running the simulation. */
if abs(m1_b_delta  - .970527109) > 1e-6 | ///
   abs(m1_se_delta - .375571575) > 1e-6 | ///
   abs(m2_b_delta  - .367635488) > 1e-6 | ///
   abs(m2_se_delta - .181409797) > 1e-6 | ///
   abs(m2_b_gpt    - .764434374) > 1e-6 | ///
   abs(m2_se_gpt   - .055568589) > 1e-6 | ///
   abs(m1_r2       - .818441995) > 1e-6 | ///
   abs(m2_r2       - .960865249) > 1e-6 {
    display as error "The conventional estimates do not match the frozen Table 11 sample."
    log close table11log
    exit 459
}

/* ---------------------------------------------------------------------- */
/* 4. Null-imposed firm-level wild-cluster bootstrap for GPT_H.           */
/* ---------------------------------------------------------------------- */

capture which boottest
if _rc {
    display as error "The required boottest package is not installed."
    display as error "Install it once with: ssc install boottest"
    log close table11log
    exit 499
}
which boottest

/* Refit Model 2 because the VIF auxiliary regressions changed e(). */
quietly regress human delta_green_pp gpt_h i.firm_id i.year

capture noisily boottest gpt_h, cluster(firm_id) bootcluster(firm_id) ///
    weight(webb) statistic(t) ptype(symmetric) ///
    reps(99999) seed(20260729) level(95) ///
    ptolerance(1e-6) format(%9.6f) nograph
local boot_rc = _rc
if `boot_rc' {
    display as error "boottest failed; inspect `logfile' for the diagnostic output."
    log close table11log
    exit `boot_rc'
}

scalar wcb_p      = r(p)
scalar wcb_reps   = r(reps)
scalar wcb_null   = r(null)
matrix WCB_CI     = r(CI)

if rowsof(WCB_CI) != 1 | colsof(WCB_CI) < 2 {
    display as error "boottest returned a disjoint or malformed confidence set."
    display as error "Inspect the log and do not report a single interval."
    log close table11log
    exit 498
}

scalar wcb_lo = WCB_CI[1,1]
scalar wcb_hi = WCB_CI[1,2]

if wcb_null != 1 {
    display as error "Warning: boottest did not report null-imposed inference."
}
if round(wcb_p, .001) != .002 | ///
   round(wcb_lo, .001) != .613 | ///
   round(wcb_hi, .001) != .841 {
    display as error "Warning: bootstrap results differ from the manuscript after rounding."
    display as error "Use the values generated below and update Table 11 if necessary."
}

/* ---------------------------------------------------------------------- */
/* 5. Display an auditable calculation summary.                           */
/* ---------------------------------------------------------------------- */

display as text _newline "Table 11 calculation results"
display as text "Model 1 DeltaGreen: b=" %9.6f m1_b_delta ///
    "  SE=" %9.6f m1_se_delta "  t=" %7.3f m1_t_delta ///
    "  p=" %8.6f m1_p_delta
display as text "Model 2 DeltaGreen: b=" %9.6f m2_b_delta ///
    "  SE=" %9.6f m2_se_delta "  t=" %7.3f m2_t_delta ///
    "  p=" %8.6f m2_p_delta
display as text "Model 2 GPT_H:      b=" %9.6f m2_b_gpt ///
    "  SE=" %9.6f m2_se_gpt "  t=" %7.3f m2_t_gpt ///
    "  conventional p=" %8.6f m2_p_gpt
display as text "GPT_H firm-level WCB: p=" %8.6f wcb_p ///
    "  95% CI=[" %9.6f wcb_lo ", " %9.6f wcb_hi "]"
display as text "Model 1 R2/adjusted R2: " %8.6f m1_r2 ///
    " / " %8.6f m1_r2a
display as text "Model 2 R2/adjusted R2: " %8.6f m2_r2 ///
    " / " %8.6f m2_r2a
display as text "Change in adjusted R2:  " %8.6f delta_r2a
display as text "Partial F(" %2.0f partial_dfnum "," ///
    %2.0f partial_dfden ")=" %9.4f partial_f ///
    "  p=" %8.6f partial_p
display as text "VIF DeltaGreen/GPT_H:   " %6.3f vif_delta ///
    " / " %6.3f vif_gpt
display as text "Bootstrap replications: " %9.0f wcb_reps

/* ---------------------------------------------------------------------- */
/* 6. Format calculated values for the Word-compatible RTF table.        */
/* ---------------------------------------------------------------------- */

local emdash "\u8212?"

local m1_delta = strtrim(string(m1_b_delta, "%9.3f"))
local m2_delta = strtrim(string(m2_b_delta, "%9.3f"))
if m1_p_delta < .01 {
    local m1_delta "`m1_delta'***"
}
else if m1_p_delta < .05 {
    local m1_delta "`m1_delta'**"
}
if m2_p_delta < .01 {
    local m2_delta "`m2_delta'***"
}
else if m2_p_delta < .05 {
    local m2_delta "`m2_delta'**"
}

local m1_delta_se = "(" + strtrim(string(m1_se_delta, "%9.3f")) + ")"
local m2_delta_se = "(" + strtrim(string(m2_se_delta, "%9.3f")) + ")"
local m1_delta_t  = strtrim(string(m1_t_delta, "%9.2f"))
local m2_delta_t  = strtrim(string(m2_t_delta, "%9.2f"))
local m1_delta_p  = cond(m1_p_delta < .001, "<0.001", ///
    strtrim(string(m1_p_delta, "%9.3f")))
local m2_delta_p  = cond(m2_p_delta < .001, "<0.001", ///
    strtrim(string(m2_p_delta, "%9.3f")))

local m2_gpt = strtrim(string(m2_b_gpt, "%9.3f"))
if m2_p_gpt < .01 {
    local m2_gpt "`m2_gpt'***"
}
else if m2_p_gpt < .05 {
    local m2_gpt "`m2_gpt'**"
}

local m2_gpt_se = "(" + strtrim(string(m2_se_gpt, "%9.3f")) + ")"
local m2_gpt_t  = strtrim(string(m2_t_gpt, "%9.2f"))
local m2_gpt_p  = cond(m2_p_gpt < .001, "<0.001", ///
    strtrim(string(m2_p_gpt, "%9.3f")))
local wcb_p_s   = cond(wcb_p < .001, "<0.001", ///
    strtrim(string(wcb_p, "%9.3f")))
local wcb_ci_s = "[" + strtrim(string(wcb_lo, "%9.3f")) + ///
    ", " + strtrim(string(wcb_hi, "%9.3f")) + "]"

local n1_s  = strtrim(string(m1_n, "%9.0f"))
local n2_s  = strtrim(string(m2_n, "%9.0f"))
local r21_s = strtrim(string(m1_r2, "%9.4f"))
local r22_s = strtrim(string(m2_r2, "%9.4f"))
local ra1_s = strtrim(string(m1_r2a, "%9.4f"))
local ra2_s = strtrim(string(m2_r2a, "%9.4f"))
local dra_s = strtrim(string(delta_r2a, "%9.4f"))
local r21_s = subinstr("`r21_s'", "0.", ".", 1)
local r22_s = subinstr("`r22_s'", "0.", ".", 1)
local ra1_s = subinstr("`ra1_s'", "0.", ".", 1)
local ra2_s = subinstr("`ra2_s'", "0.", ".", 1)
local dra_s = subinstr("`dra_s'", "0.", ".", 1)
local pf_s = "F(" + strtrim(string(partial_dfnum, "%9.0f")) + ///
    "," + strtrim(string(partial_dfden, "%9.0f")) + ")=" + ///
    strtrim(string(partial_f, "%9.2f"))
local pp_s  = cond(partial_p < .001, "<0.001", ///
    strtrim(string(partial_p, "%9.3f")))
local vif_delta_s = strtrim(string(vif_delta, "%9.2f"))
local vif_gpt_s   = strtrim(string(vif_gpt, "%9.2f"))

/* Helper for a three-column RTF row. */
capture program drop table11_rtf_row
program define table11_rtf_row
    version 17.0
    syntax, Handle(name) Label(string asis) ///
        Modelone(string asis) Modeltwo(string asis) [Bottom]

    local border ""
    if `"`bottom'"' != "" {
        local border "\clbrdrb\brdrs\brdrw10"
    }

    file write `handle' "\trowd\trgaph80\trleft0" ///
        "`border'\cellx5200" ///
        "`border'\cellx7600" ///
        "`border'\cellx10080" _n
    file write `handle' ///
        "\pard\intbl\ql `label'\cell" ///
        "\pard\intbl\qc `modelone'\cell" ///
        "\pard\intbl\qc `modeltwo'\cell\row" _n
end

file open rtf using "`outfile'", write text replace
file write rtf "{\rtf1\ansi\ansicpg1252\deff0" _n
file write rtf "{\fonttbl{\f0 Times New Roman;}}" _n
file write rtf ///
    "\paperw12240\paperh15840\margl1080\margr1080" ///
    "\margt900\margb900\fs20" _n
file write rtf ///
    "\pard\sa80\b Table 11: Nested fixed-effects regressions of " ///
    "blinded human rubric scores\b0\par" _n

/* Header row with top and bottom rules. */
file write rtf "\trowd\trgaph80\trleft0" ///
    "\clbrdrt\brdrs\brdrw10\clbrdrb\brdrs\brdrw10\cellx5200" ///
    "\clbrdrt\brdrs\brdrw10\clbrdrb\brdrs\brdrw10\cellx7600" ///
    "\clbrdrt\brdrs\brdrw10\clbrdrb\brdrs\brdrw10\cellx10080" _n
file write rtf ///
    "\pard\intbl\ql\b Variable\b0\cell" ///
    "\pard\intbl\qc\b Model 1\b0\cell" ///
    "\pard\intbl\qc\b Model 2\b0\cell\row" _n

table11_rtf_row, handle(rtf) label("\u916?Green") ///
    modelone("`m1_delta'") modeltwo("`m2_delta'")
table11_rtf_row, handle(rtf) label("SE") ///
    modelone("`m1_delta_se'") modeltwo("`m2_delta_se'")
table11_rtf_row, handle(rtf) label("t") ///
    modelone("`m1_delta_t'") modeltwo("`m2_delta_t'")
table11_rtf_row, handle(rtf) label("p") ///
    modelone("`m1_delta_p'") modeltwo("`m2_delta_p'")

table11_rtf_row, handle(rtf) label("GPT_H") ///
    modelone("`emdash'") modeltwo("`m2_gpt'")
table11_rtf_row, handle(rtf) label("SE") ///
    modelone("`emdash'") modeltwo("`m2_gpt_se'")
table11_rtf_row, handle(rtf) label("t") ///
    modelone("`emdash'") modeltwo("`m2_gpt_t'")
table11_rtf_row, handle(rtf) label("Conventional p") ///
    modelone("`emdash'") modeltwo("`m2_gpt_p'")
table11_rtf_row, handle(rtf) label("WCB p") ///
    modelone("`emdash'") modeltwo("`wcb_p_s'")
table11_rtf_row, handle(rtf) label("WCB 95% CI") ///
    modelone("`emdash'") modeltwo("`wcb_ci_s'")

table11_rtf_row, handle(rtf) label("Firm FE") ///
    modelone("Yes") modeltwo("Yes")
table11_rtf_row, handle(rtf) label("Year FE") ///
    modelone("Yes") modeltwo("Yes")
table11_rtf_row, handle(rtf) label("N") ///
    modelone("`n1_s'") modeltwo("`n2_s'")
table11_rtf_row, handle(rtf) label("R{\super 2}\nosupersub") ///
    modelone("`r21_s'") modeltwo("`r22_s'")
table11_rtf_row, handle(rtf) label("Adjusted R{\super 2}\nosupersub") ///
    modelone("`ra1_s'") modeltwo("`ra2_s'")
table11_rtf_row, handle(rtf) ///
    label("Change in adjusted R{\super 2}\nosupersub") ///
    modelone("`emdash'") modeltwo("`dra_s'")
table11_rtf_row, handle(rtf) label("Partial F-test for GPT_H") ///
    modelone("`emdash'") modeltwo("`pf_s'")
table11_rtf_row, handle(rtf) label("Partial F-test p") ///
    modelone("`emdash'") modeltwo("`pp_s'") bottom

file write rtf ///
    "\pard\sa40\fs18\i Notes:\i0 The dependent variable is the " ///
    "blinded human rubric score. \u916?Green is baseline-adjusted green " ///
    "coverage expressed in percentage points. Model 1 includes \u916?Green; " ///
    "Model 2 adds GPT_H. Conventional OLS standard errors are reported in " ///
    "parentheses. Both models include firm and year fixed effects. " ///
    "GPT_H inference is supplemented with a null-imposed firm-level wild-" ///
    "cluster bootstrap (WCB) using Webb six-point weights, symmetric " ///
    "two-sided p-values, 99,999 replications, and seed 20260729. " ///
    "Untabulated VIFs are `vif_delta_s' for \u916?Green and `vif_gpt_s' " ///
    "for GPT_H. ** p < 0.05; *** p < 0.01.\par" _n
file write rtf "}" _n
file close rtf

display as result _newline "Created: `outfile'"
log close table11log
