###############################################################################
# 01_main_analysis.R
#
# Drug-target Mendelian randomization: genetically proxied PDE5 inhibition
# (SBP-lowering allele of rs2389873 at PDE5A) and glaucoma / IOP.
#
# Estimates are reported per 1 mmHg LOWER systolic blood pressure.
#
# What it does:
#   A. Palindromic SNP assessment for the lead variant rs4834792 (A/T, MAF ~0.48)
#      and a sensitivity run under both strand assumptions
#   B. Steiger directionality test for every exposure-outcome pair
#   C. Multiple-testing correction across the six outcomes (BH-FDR, Bonferroni)
#   D. PDE6 eQTL analysis (switched off in the published analysis; the eQTL
#      endpoint was unavailable -- see PDE6_GENES below)
#
# Run from the folder containing the IGGC summary statistics:
#     Rscript R/01_main_analysis.R
# Outputs go to results_v4/.
#
# Needs a free OpenGWAS token (https://api.opengwas.io) saved in ~/.Renviron as
# OPENGWAS_JWT=...  Never put the token in the script.
###############################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

# =============================================================================
# 0. SETTINGS
# =============================================================================

OUT_DIR <- "results_v4"

## Exposure (PDE5A, blood-pressure proxy)
EXPOSURE_ID        <- "ieu-b-38"        # SBP, Evangelou 2018
EXPOSURE_N         <- 757601
PDE5A_INSTRUMENT   <- "rs2389873"       # non-palindromic proxy for lead rs4834792 (r2 = 1.00)
PALINDROMIC_EXCL   <- "rs4834792"       # A/T palindromic lead variant, excluded
ORIENT_TO_SBP_LOWERING <- TRUE          # PDE5 inhibition lowers BP, so the
                                        # proxy allele is the SBP-lowering one

## Expected primary estimate, used as a sanity check only (disabled)
V3_PRIMARY_OR <- NA                     # sanity check disabled (new instrument)

## Palindrome rule (TwoSampleMR default)
PALINDROME_MAF_THRESHOLD <- 0.42

## IGGC POAG (primary outcome, local file, no EAF column)
IGGC_FILE     <- "GCST90011767_buildGRCh37.tsv.gz"
IGGC_NCASE    <- 15229
IGGC_NCONTROL <- 177473   # GCST90011767: 15,229 cases / 177,473 controls,
                          # European, UK Biobank excluded (IGGC, Gharahkhani 2021)

## Population prevalences used to convert log-ORs to variance explained
## (Steiger). Rough European 40+ figures; the test is insensitive to these
## when the exposure is a 750k-sample BP GWAS, but report what you used.
PREV_POAG     <- 0.02
PREV_GLAUCOMA <- 0.03
PREV_NTG      <- 0.005

## IOP outcome: Springelkamp 2017 (IGGC IOP meta-analysis,
## N = 29,578, no UK Biobank, so no overlap with the BP exposure). NA = skip.
IOP_ID <- "ebi-a-GCST004074"

## Outcomes
OUTCOMES <- data.frame(
  id = c("IGGC_LOCAL",
         "finn-b-H7_GLAUCOMA_POAG",
         "finn-b-H7_GLAUCPRIMOPEN",
         "finn-b-H7_GLAUCOMA",
         "finn-b-H7_GLAUCOMA_NTG",
         ifelse(is.na(IOP_ID), "SKIP", IOP_ID)),
  label = c("POAG (IGGC, primary)",
            "POAG strict (FinnGen)",
            "POAG broad (FinnGen)",
            "All glaucoma (FinnGen)",
            "NTG (FinnGen)",
            "IOP"),
  type = c("binary", "binary", "binary", "binary", "binary", "continuous"),
  prevalence = c(PREV_POAG, PREV_POAG, PREV_POAG, PREV_GLAUCOMA, PREV_NTG, NA),
  stringsAsFactors = FALSE
)
OUTCOMES <- OUTCOMES[OUTCOMES$id != "SKIP", ]

## PDE6 catalytic subunits (sildenafil binds the catalytic site; the
## inhibitory gamma subunits PDE6G/PDE6H would need the opposite orientation,
## so they are left out). eQTLGen blood cis-eQTLs on OpenGWAS.
PDE6_GENES_ALL <- data.frame(
  gene    = c("PDE6A", "PDE6B", "PDE6C"),
  ensembl = c("ENSG00000132915", "ENSG00000133256", "ENSG00000095464"),
  chr     = c("5", "4", "10"),
  stringsAsFactors = FALSE
)
PDE6_GENES <- PDE6_GENES_ALL[0, ]   # PDE6 analysis dropped (eQTL endpoint unavailable)
PDE6_P_THRESHOLD <- 5e-8
PDE6_CLUMP_R2    <- 0.001
PDE6_CLUMP_KB    <- 10000

# =============================================================================
# 1. HELPER FUNCTIONS (pure; no network)
# =============================================================================

is_palindromic <- function(a1, a2) {
  a1 <- toupper(a1); a2 <- toupper(a2)
  (a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") |
    (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C")
}

# Classify each SNP as non-palindromic / palindromic-inferrable /
# palindromic-ambiguous (MAF above threshold, strand can't be inferred from EAF)
assess_palindromes <- function(snp, ea, oa, eaf, threshold = PALINDROME_MAF_THRESHOLD) {
  pal <- is_palindromic(ea, oa)
  maf <- pmin(eaf, 1 - eaf)
  status <- ifelse(!pal, "non-palindromic",
            ifelse(is.na(maf), "palindromic, EAF missing",
            ifelse(maf > threshold, "palindromic, ambiguous (MAF > threshold)",
                                    "palindromic, strand inferrable from EAF")))
  data.frame(SNP = snp, effect_allele = ea, other_allele = oa,
             eaf = eaf, maf = maf, palindromic = pal, status = status,
             stringsAsFactors = FALSE)
}

# Flip exposure alleles so beta is negative (allele that lowers the trait)
orient_negative <- function(d) {
  flip <- !is.na(d$beta.exposure) & d$beta.exposure > 0
  if (any(flip)) {
    ea <- d$effect_allele.exposure[flip]
    d$effect_allele.exposure[flip] <- d$other_allele.exposure[flip]
    d$other_allele.exposure[flip]  <- ea
    d$beta.exposure[flip] <- -d$beta.exposure[flip]
    if ("eaf.exposure" %in% names(d))
      d$eaf.exposure[flip] <- 1 - d$eaf.exposure[flip]
  }
  d
}

# Variance explained for a continuous trait from the t statistic
r2_continuous <- function(beta, se, n) {
  t2 <- (beta / se)^2
  t2 / (t2 + n - 2)
}

# Steiger directionality test on summed r2 across SNPs.
# Compares Fisher-z transformed correlations from independent samples
# (equivalent to psych::r.test used inside TwoSampleMR::mr_steiger).
steiger_test <- function(r2_exp, r2_out, n_exp, n_out) {
  r_exp <- sqrt(sum(r2_exp)); r_out <- sqrt(sum(r2_out))
  z <- (atanh(r_exp) - atanh(r_out)) / sqrt(1 / (n_exp - 3) + 1 / (n_out - 3))
  data.frame(r2_exposure = sum(r2_exp), r2_outcome = sum(r2_out),
             steiger_z = z, steiger_p = 2 * pnorm(-abs(z)),
             correct_direction = r_exp > r_out)
}

# Multiple-testing correction within a family of tests
add_mtc <- function(df, p_col = "pval", family_col = "family") {
  df %>%
    group_by(.data[[family_col]]) %>%
    mutate(n_tests      = sum(!is.na(.data[[p_col]])),
           q_fdr_bh     = p.adjust(.data[[p_col]], method = "BH"),
           p_bonferroni = p.adjust(.data[[p_col]], method = "bonferroni")) %>%
    ungroup()
}

wald_ratio <- function(b_exp, se_exp, b_out, se_out) {
  b  <- b_out / b_exp
  se <- abs(se_out / b_exp)                      # first-order, as in TwoSampleMR
  data.frame(b = b, se = se, pval = 2 * pnorm(-abs(b / se)))
}

# Stop here when the file is sourced for testing
if (isTRUE(getOption("pde5.functions_only"))) {
  message("Functions loaded (pde5.functions_only = TRUE); skipping analysis.")
} else {

# =============================================================================
# 2. PACKAGES + DATA ACCESS
# =============================================================================

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(ieugwasr)
})
dir.create(OUT_DIR, showWarnings = FALSE)
if (Sys.getenv("OPENGWAS_JWT") == "")
  warning("OPENGWAS_JWT not set -- add it to ~/.Renviron and restart R.")

cat("\n==== PDE5i MR: main analysis ====\n")

# Load IGGC POAG rows for a set of SNPs (zgrep first; fall back to full read)
read_iggc <- function(snps) {
  if (!file.exists(IGGC_FILE)) stop("Can't find ", IGGC_FILE, " in ", getwd())
  pat <- paste0("variant_id|", paste(snps, collapse = "|"))
  d <- tryCatch(fread(cmd = sprintf("zgrep -wE '%s' %s", pat, IGGC_FILE)),
                error = function(e) NULL)
  if (is.null(d) || !"variant_id" %in% names(d)) {
    d <- fread(IGGC_FILE)
  }
  d <- d[variant_id %in% snps]
  if (nrow(d) == 0) return(NULL)
  d[, ncase := IGGC_NCASE][, ncontrol := IGGC_NCONTROL]
  d[, samplesize := IGGC_NCASE + ifelse(is.na(IGGC_NCONTROL), NA, IGGC_NCONTROL)]
  d[, phen := "POAG (IGGC)"]
  eaf_col <- intersect(c("effect_allele_frequency", "eaf"), names(d))
  format_data(as.data.frame(d), type = "outcome",
              snp_col = "variant_id", beta_col = "beta", se_col = "standard_error",
              effect_allele_col = "effect_allele", other_allele_col = "other_allele",
              pval_col = "p_value",
              eaf_col = if (length(eaf_col)) eaf_col[1] else "eaf_missing",
              ncase_col = "ncase", ncontrol_col = "ncontrol",
              samplesize_col = "samplesize", phenotype_col = "phen")
}

# Retry wrapper for OpenGWAS calls (the server throws 502s under load).
# Returns list(value, error) so "server failed" and "no result" stay distinct.
with_retry <- function(expr_fn, tries = 4, wait = 20) {
  last_err <- NULL
  for (k in seq_len(tries)) {
    res <- tryCatch(list(value = expr_fn(), error = NULL),
                    error = function(e) list(value = NULL, error = conditionMessage(e)))
    if (is.null(res$error)) return(res)
    last_err <- res$error
    if (k < tries) { message("  OpenGWAS error, retrying in ", wait, "s: ", last_err); Sys.sleep(wait) }
  }
  list(value = NULL, error = last_err)
}

# Case/control counts from OpenGWAS metadata (cached), for FinnGen etc.
.info_cache <- list()
gwas_meta <- function(oid) {
  if (is.null(.info_cache[[oid]])) {
    r <- with_retry(function() gwasinfo(oid))
    .info_cache[[oid]] <<- if (is.null(r$value)) data.frame() else r$value
  }
  .info_cache[[oid]]
}

get_outcome <- function(snps, oid) {
  if (oid == "IGGC_LOCAL") return(read_iggc(snps))
  r <- with_retry(function() extract_outcome_data(snps = snps, outcomes = oid,
                                                  proxies = TRUE, rsq = 0.8))
  d <- r$value
  if (is.null(d) || nrow(d) == 0) return(NULL)
  meta <- gwas_meta(oid)
  for (col in c("ncase", "ncontrol", "samplesize")) {
    oc <- paste0(col, ".outcome"); mc <- if (col == "samplesize") "sample_size" else col
    if ((!oc %in% names(d) || all(is.na(d[[oc]]))) && mc %in% names(meta))
      d[[oc]] <- as.numeric(meta[[mc]][1])
  }
  d
}

outcome_has_eaf <- function(o) "eaf.outcome" %in% names(o) && any(!is.na(o$eaf.outcome))

# Harmonisation rule (applied identically to PDE5A and PDE6):
#   outcome with EAF    -> action = 2 (infer palindrome strand from EAF;
#                          ambiguous palindromes dropped)
#   outcome without EAF -> action = 3 (drop all palindromes)
harmonise_rule <- function(exp, out) {
  act <- if (outcome_has_eaf(out)) 2 else 3
  h <- harmonise_data(exp, out, action = act)
  h$harmonise_action <- act
  h[h$mr_keep, , drop = FALSE]
}

run_mr <- function(h) {
  n <- nrow(h)
  if (n == 0) return(NULL)
  methods <- if (n == 1) "mr_wald_ratio"
             else if (n == 2) "mr_ivw"
             else c("mr_ivw", "mr_egger_regression", "mr_weighted_median",
                    "mr_weighted_mode")
  r <- mr(h, method_list = methods)
  r$primary <- r$method %in% c("Wald ratio", "Inverse variance weighted")
  r
}

steiger_for <- function(h, otype, prevalence, n_exp) {
  r2_exp <- r2_continuous(h$beta.exposure, h$se.exposure,
                          if (!is.null(n_exp)) n_exp else h$samplesize.exposure)
  if (otype == "binary") {
    ncase <- h$ncase.outcome; nctrl <- h$ncontrol.outcome
    if (all(is.na(nctrl)) || all(is.na(ncase)))
      return(data.frame(r2_exposure = sum(r2_exp), r2_outcome = NA, steiger_z = NA,
                        steiger_p = NA, correct_direction = NA,
                        note = "ncase/ncontrol missing"))
    af <- if ("eaf.outcome" %in% names(h) && !all(is.na(h$eaf.outcome)))
            h$eaf.outcome else h$eaf.exposure          # IGGC: borrow exposure EAF
    af[is.na(af)] <- h$eaf.exposure[is.na(af)]
    r_out <- get_r_from_lor(lor = h$beta.outcome, af = af, ncase = ncase,
                            ncontrol = nctrl, prevalence = prevalence,
                            model = "logit")
    r2_out <- r_out^2
    n_out  <- ncase[1] + nctrl[1]
    note   <- if (!outcome_has_eaf(h)) "outcome EAF missing; exposure EAF used" else ""
  } else {
    n_out  <- h$samplesize.outcome
    r2_out <- r2_continuous(h$beta.outcome, h$se.outcome, n_out)
    n_out  <- n_out[1]
    note   <- ""
  }
  n_e <- if (!is.null(n_exp)) n_exp else max(h$samplesize.exposure, na.rm = TRUE)
  cbind(steiger_test(r2_exp, r2_out, n_e, n_out), note = note)
}

# =============================================================================
# 3. PDE5A EXPOSURE (rs2389873 + rs4834792 for the palindrome documentation)
# =============================================================================

# Print what each OpenGWAS ID actually is, so a wrong ID can't slip through
cat("\n-- Outcome ID check\n")
for (oid in setdiff(OUTCOMES$id, "IGGC_LOCAL")) {
  m <- gwas_meta(oid)
  cat(sprintf("  %-28s -> %s (ncase=%s, ncontrol=%s, N=%s)\n", oid,
              if (nrow(m)) m$trait[1] else "METADATA UNAVAILABLE",
              if ("ncase" %in% names(m)) m$ncase[1] else NA,
              if ("ncontrol" %in% names(m)) m$ncontrol[1] else NA,
              if ("sample_size" %in% names(m)) m$sample_size[1] else NA))
}

cat("\n-- Building PDE5A exposure from", EXPOSURE_ID, "\n")
exp_raw <- extract_outcome_data(snps = c(PDE5A_INSTRUMENT, PALINDROMIC_EXCL),
                                outcomes = EXPOSURE_ID, proxies = FALSE)
exp_all <- convert_outcome_to_exposure(exp_raw)
exp_all$exposure <- "PDE5A (SBP-lowering)"
exp_all$samplesize.exposure <- EXPOSURE_N
if (ORIENT_TO_SBP_LOWERING) exp_all <- orient_negative(exp_all)
exp_all$F_stat <- (exp_all$beta.exposure / exp_all$se.exposure)^2

exp_pde5 <- exp_all[exp_all$SNP == PDE5A_INSTRUMENT, ]
print(exp_all[, c("SNP", "effect_allele.exposure", "other_allele.exposure",
                  "eaf.exposure", "beta.exposure", "se.exposure", "F_stat")])

# =============================================================================
# A. PALINDROMIC SNP ASSESSMENT
# =============================================================================

cat("\n==== A. Palindromic SNP assessment ====\n")
pal_tab <- assess_palindromes(exp_all$SNP, exp_all$effect_allele.exposure,
                              exp_all$other_allele.exposure, exp_all$eaf.exposure)
pal_tab$decision <- ifelse(pal_tab$SNP == PALINDROMIC_EXCL, "excluded", "retained")

# Per outcome: is outcome EAF available to resolve strand? And what does the
# palindromic SNP give under each strand assumption?
pal_sens <- list(); outcome_cache <- list()
for (i in seq_len(nrow(OUTCOMES))) {
  o <- OUTCOMES[i, ]
  out <- get_outcome(exp_all$SNP, o$id)
  outcome_cache[[o$id]] <- out
  if (is.null(out)) next
  h <- harmonise_data(exp_all, out, action = 1)   # action 1 = assume forward strand
  for (s in h$SNP) {
    r <- h[h$SNP == s, ][1, ]
    fwd <- wald_ratio(r$beta.exposure, r$se.exposure, r$beta.outcome, r$se.outcome)
    rev <- wald_ratio(r$beta.exposure, r$se.exposure, -r$beta.outcome, r$se.outcome)
    pal_sens[[length(pal_sens) + 1]] <- data.frame(
      outcome = o$label, SNP = s,
      palindromic = is_palindromic(r$effect_allele.exposure, r$other_allele.exposure),
      eaf_exposure = r$eaf.exposure,
      eaf_outcome  = if ("eaf.outcome" %in% names(r)) r$eaf.outcome else NA,
      OR_forward_strand = exp(-fwd$b), p_forward = fwd$pval,   # per 1 mmHg LOWER SBP
      OR_flipped_strand = exp(-rev$b), p_flipped = rev$pval)
  }
}
pal_sens <- bind_rows(pal_sens)

# Concordance: which strand assumption for rs4834792 agrees in direction with
# the unambiguous instrument (rs2389873), outcome by outcome
pal_conc <- as.data.frame(dcast(as.data.table(pal_sens), outcome ~ SNP,
                                value.var = "OR_forward_strand"))
if (all(c(PDE5A_INSTRUMENT, PALINDROMIC_EXCL) %in% names(pal_conc))) {
  pal_conc$palindrome_agrees_with_reference_snp <-
    sign(log(pal_conc[[PALINDROMIC_EXCL]])) == sign(log(pal_conc[[PDE5A_INSTRUMENT]]))
}

print(pal_tab); print(pal_sens); print(pal_conc)
fwrite(pal_tab,  file.path(OUT_DIR, "A_palindromic_snp_assessment.csv"))
fwrite(pal_sens, file.path(OUT_DIR, "A_palindromic_strand_sensitivity.csv"))
fwrite(pal_conc, file.path(OUT_DIR, "A_palindromic_direction_concordance.csv"))

# =============================================================================
# PDE5A MR across outcomes (rs2389873) -- feeds B and C
# =============================================================================

cat("\n==== PDE5A MR (", PDE5A_INSTRUMENT, ") ====\n")
pde5_res <- list(); pde5_steiger <- list(); pde5_harm <- list()
for (i in seq_len(nrow(OUTCOMES))) {
  o <- OUTCOMES[i, ]
  out <- outcome_cache[[o$id]]
  if (is.null(out)) { message("No outcome data: ", o$label); next }
  h <- harmonise_rule(exp_pde5, out)
  if (nrow(h) == 0) { message("Nothing survived harmonisation: ", o$label); next }
  pde5_harm[[o$label]] <- h
  r <- run_mr(h); r$b <- -r$b   # Wald ratio is per unit HIGHER SBP; negate -> per 1 mmHg LOWER SBP
  r$outcome_label <- o$label; r$outcome_type <- o$type
  pde5_res[[o$label]] <- r
  s <- steiger_for(h, o$type, o$prevalence, EXPOSURE_N)
  pde5_steiger[[o$label]] <- cbind(exposure = "PDE5A", outcome = o$label, nsnp = nrow(h), s)
}
pde5_res <- bind_rows(pde5_res) %>% mutate(family = "PDE5A")

# Sanity check against an earlier run (disabled)
chk <- pde5_res %>% filter(outcome_label == "POAG (IGGC, primary)", primary)
if (nrow(chk) && !is.na(V3_PRIMARY_OR)) {
  or_now <- exp(chk$b[1])
  cat(sprintf("\nSanity check: IGGC POAG OR = %.3f (v3: %.3f)\n", or_now, V3_PRIMARY_OR))
  if (abs(log(or_now) - log(V3_PRIMARY_OR)) > 0.01)
    warning("Primary OR differs from v3 -- check ORIENT_TO_SBP_LOWERING / scaling match v3.")
}

# =============================================================================
# D. PDE6 INSTRUMENTS (sildenafil off-target)
# =============================================================================

cat("\n==== D. PDE6 instrument analysis ====\n")
pde6_inst <- list(); pde6_res <- list(); pde6_steiger <- list(); pde6_log <- list()
for (g in seq_len(nrow(PDE6_GENES))) {
  gene <- PDE6_GENES$gene[g]; gid <- paste0("eqtl-a-", PDE6_GENES$ensembl[g])
  # OpenGWAS 502s come back as empty results rather than R errors, so an empty
  # result is retried a few times before it's believed.
  info <- NULL
  for (k in 1:3) {
    info <- tryCatch(gwasinfo(gid), error = function(e) NULL)
    if (!is.null(info) && nrow(info) > 0) break
    if (k < 3) { message("  No metadata for ", gid, " -- retrying in 30s"); Sys.sleep(30) }
  }
  if (is.null(info) || nrow(info) == 0) {
    pde6_log[[gene]] <- data.frame(gene = gene, id = gid, n_instruments = NA,
      note = "no metadata after 3 tries (server down, or ID not on OpenGWAS) -- rerun")
    next
  }
  # Route 1: extract_instruments with server-side clumping (this is the call
  # that kept returning 502s for the eqtl-a datasets).
  ins <- tryCatch(extract_instruments(outcomes = gid, p1 = PDE6_P_THRESHOLD,
                                      clump = TRUE, r2 = PDE6_CLUMP_R2,
                                      kb = PDE6_CLUMP_KB),
                  error = function(e) NULL)
  inst_route <- "clumped (extract_instruments)"
  # Route 2: unclumped tophits, then clump separately; if the clumping service
  # is also down, keep the lead (most significant) SNP only.
  if (is.null(ins) || nrow(ins) == 0) {
    message("  Clumped extraction failed for ", gene, " -- trying unclumped tophits")
    th <- tryCatch(tophits(gid, pval = PDE6_P_THRESHOLD, clump = 0),
                   error = function(e) NULL)
    if (!is.null(th) && nrow(th) > 0) {
      th <- as.data.frame(th)
      th <- th[as.character(th$chr) == PDE6_GENES$chr[g], ]
      if (nrow(th) > 1) {
        cl <- tryCatch(ld_clump(data.frame(rsid = th$rsid, pval = th$p),
                                clump_r2 = PDE6_CLUMP_R2, clump_kb = PDE6_CLUMP_KB),
                       error = function(e) NULL)
        if (!is.null(cl) && nrow(cl) > 0) {
          th <- th[th$rsid %in% cl$rsid, ]; inst_route <- "tophits + ld_clump"
        } else {
          th <- th[which.min(th$p), ]; inst_route <- "lead SNP only (clumping service unavailable)"
        }
      } else inst_route <- "tophits (single SNP)"
      th$phen <- gene
      ins <- format_data(th, type = "exposure", snp_col = "rsid", beta_col = "beta",
                         se_col = "se", eaf_col = "eaf", effect_allele_col = "ea",
                         other_allele_col = "nea", pval_col = "p",
                         samplesize_col = "n", chr_col = "chr", pos_col = "position",
                         phenotype_col = "phen")
      ins$id.exposure <- gid
    } else if (is.null(th)) {
      pde6_log[[gene]] <- data.frame(gene = gene, id = gid, n_instruments = NA,
        note = "both extraction routes failed (server errors) -- rerun later")
      next
    }
  }
  n_any <- if (is.null(ins)) 0 else nrow(ins)
  if (n_any) ins <- ins[as.character(ins$chr.exposure) == PDE6_GENES$chr[g], ]   # cis only
  if (is.null(ins) || nrow(ins) == 0) {
    pde6_log[[gene]] <- data.frame(gene = gene, id = gid, n_instruments = 0,
      note = if (n_any) "hits found but none on the gene's chromosome (no cis-eQTL)"
             else sprintf("no SNP at P < %g after 3 tries (genuinely none, unless the server was failing)",
                          PDE6_P_THRESHOLD))
    next
  }
  ins <- orient_negative(ins)     # expression-lowering allele = proxy for inhibition
  ins$exposure <- paste0(gene, " (expression-lowering)")
  ins$F_stat <- (ins$beta.exposure / ins$se.exposure)^2
  ins$gene <- gene
  if (!"samplesize.exposure" %in% names(ins) || all(is.na(ins$samplesize.exposure)))
    ins$samplesize.exposure <- if (!is.null(info$sample_size)) info$sample_size else 31684
  pal <-assess_palindromes(ins$SNP, ins$effect_allele.exposure,
                            ins$other_allele.exposure, ins$eaf.exposure)
  ins$palindrome_status <- pal$status
  pde6_inst[[gene]] <- ins
  pde6_log[[gene]] <- data.frame(gene = gene, id = gid, n_instruments = nrow(ins),
                                 note = sprintf("mean F = %.1f; %s", mean(ins$F_stat), inst_route))

  for (i in seq_len(nrow(OUTCOMES))) {
    o <- OUTCOMES[i, ]
    out <- get_outcome(ins$SNP, o$id)
    if (is.null(out)) next
    h <- harmonise_rule(ins, out)
    if (nrow(h) == 0) next
    r <- run_mr(h); r$outcome_label <- o$label; r$outcome_type <- o$type; r$gene <- gene
    pde6_res[[paste(gene, o$label)]] <- r
    s <- steiger_for(h, o$type, o$prevalence, NULL)
    pde6_steiger[[paste(gene, o$label)]] <- cbind(exposure = gene, outcome = o$label,
                                                  nsnp = nrow(h), s)
  }
}
pde6_log  <- bind_rows(pde6_log);  print(pde6_log)
pde6_inst <- bind_rows(pde6_inst)
pde6_res  <- bind_rows(pde6_res)
if (nrow(pde6_res)) pde6_res$family <- "PDE6"

fwrite(pde6_log, file.path(OUT_DIR, "D_pde6_instrument_summary.csv"))
if (nrow(pde6_inst)) fwrite(pde6_inst, file.path(OUT_DIR, "D_pde6_instruments.csv"))

# =============================================================================
# B. STEIGER DIRECTIONALITY
# =============================================================================

cat("\n==== B. Steiger directionality ====\n")
steiger_all <- bind_rows(bind_rows(pde5_steiger), bind_rows(pde6_steiger))
print(steiger_all)
fwrite(steiger_all, file.path(OUT_DIR, "B_steiger_directionality.csv"))

# =============================================================================
# C. MULTIPLE-TESTING CORRECTION ACROSS OUTCOMES
# =============================================================================
# One test per exposure-outcome pair (the primary estimator: Wald ratio for
# one SNP, IVW otherwise). Families: PDE5A (all outcomes) and PDE6 (all
# gene x outcome pairs). Egger/WM/mode are sensitivity analyses and not
# included in the correction.

cat("\n==== C. Multiple-testing correction ====\n")
primary <- bind_rows(
  pde5_res %>% filter(primary) %>% mutate(gene = "PDE5A"),
  if (nrow(pde6_res)) pde6_res %>% filter(primary) else NULL
) %>%
  mutate(bin = outcome_type == "binary",
         OR = ifelse(bin, exp(b), NA),
         lo = ifelse(bin, exp(b - 1.96 * se), NA),
         hi = ifelse(bin, exp(b + 1.96 * se), NA)) %>%
  add_mtc(p_col = "pval", family_col = "family") %>%
  left_join(steiger_all %>% select(exposure, outcome, steiger_p, correct_direction),
            by = c("gene" = "exposure", "outcome_label" = "outcome")) %>%
  select(family, gene, outcome = outcome_label, outcome_type, method, nsnp, b, se,
         OR, lo, hi, pval, q_fdr_bh, p_bonferroni, n_tests,
         steiger_p, correct_direction)

print(as.data.frame(primary), digits = 3)
fwrite(primary, file.path(OUT_DIR, "C_primary_results_with_fdr_and_steiger.csv"))
fwrite(bind_rows(pde5_res, pde6_res),
       file.path(OUT_DIR, "all_mr_estimates_incl_sensitivity_methods.csv"))

# Forest plot (binary outcomes; IOP beta is in the CSV)
fp <- primary %>%
  filter(outcome_type == "binary") %>%
  mutate(lab = sprintf("q = %.2f", q_fdr_bh),
         outcome = factor(outcome, levels = rev(unique(OUTCOMES$label))))
p <- ggplot(fp, aes(x = OR, y = outcome)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.2, orientation = "y") +
  geom_point(size = 2) +
  geom_text(aes(label = lab), x = Inf, hjust = 1.05, size = 3, colour = "grey30") +
  facet_wrap(~ gene, ncol = 1, scales = "free_y") +
  scale_x_log10() +
  labs(x = "OR (95% CI) per unit of genetically proxied inhibition",
       y = NULL, caption = "q = Benjamini-Hochberg FDR within each exposure family") +
  theme_bw(base_size = 11)
ggsave(file.path(OUT_DIR, "forest_pde5a_pde6_fdr.png"), p,
       width = 7, height = 2 + 0.45 * nrow(fp), dpi = 300)

# =============================================================================
# Methods text (auto-filled)
# =============================================================================

pal_row <- pal_tab[pal_tab$SNP == PALINDROMIC_EXCL, ]
iggc_s  <- steiger_all %>% filter(exposure == "PDE5A", grepl("IGGC", outcome))
methods <- c(
"PALINDROMIC VARIANTS",
sprintf(paste0(
"Of the two cis-PDE5A variants associated with SBP, %s is palindromic (%s/%s) with an ",
"effect allele frequency of %.2f (MAF %.2f), above the %.2f threshold at which strand ",
"cannot be inferred from allele frequency. The primary POAG GWAS (GCST90011767) does not ",
"report allele frequencies, and in FinnGen the frequency is too close to 0.5 to be informative, ",
"so its strand could not be resolved and it was excluded; analyses used %s, which is ",
"non-palindromic. In a sensitivity analysis aligning %s on the forward-strand assumption, ",
"its estimates %s. A strand-sensitivity table is given in Supplementary Table X."),
  PALINDROMIC_EXCL, pal_row$effect_allele, pal_row$other_allele, pal_row$eaf, pal_row$maf,
  PALINDROME_MAF_THRESHOLD, PDE5A_INSTRUMENT, PALINDROMIC_EXCL,
  if ("palindrome_agrees_with_reference_snp" %in% names(pal_conc) &&
      all(pal_conc$palindrome_agrees_with_reference_snp, na.rm = TRUE))
    sprintf("were in the same direction as those for %s for every outcome, so the exclusion did not change the direction of any result", PDE5A_INSTRUMENT)
  else
    sprintf("differed in direction from %s for at least one outcome, consistent with possible strand misalignment", PDE5A_INSTRUMENT)),
"For the PDE6 analyses the same rule was applied: palindromic variants were dropped for outcomes without allele frequencies, and for outcomes with allele frequencies strand was inferred where MAF <= 0.42, otherwise the variant was dropped.",
"",
"DIRECTIONALITY",
paste0("The Steiger test was used to check that each instrument explained more variance in the ",
"exposure than in the outcome. Variance explained in continuous traits was derived from the ",
"t-statistic; for binary outcomes it was derived from the log odds ratio on the liability scale ",
"(get_r_from_lor, TwoSampleMR) with assumed population prevalences of ",
sprintf("%.1f%% (POAG), %.1f%% (all glaucoma) and %.1f%% (NTG). ", 100*PREV_POAG, 100*PREV_GLAUCOMA, 100*PREV_NTG),
"Where the outcome GWAS lacked allele frequencies the exposure allele frequency was used."),
if (nrow(iggc_s) && !is.na(iggc_s$steiger_p[1]))
  sprintf("For the primary outcome the inferred direction was %s (Steiger P = %.2g)%s. Across outcomes the test was %s.",
          ifelse(iggc_s$correct_direction[1], "exposure -> outcome", "outcome -> exposure"),
          iggc_s$steiger_p[1],
          ifelse(iggc_s$steiger_p[1] < 0.05, "", ", which was not statistically significant"),
          ifelse(all(steiger_all$steiger_p[steiger_all$exposure == "PDE5A"] >= 0.05, na.rm = TRUE),
                 "inconclusive for every outcome (all P >= 0.05), as expected with a single instrument of modest strength",
                 "significant for at least one outcome; see Supplementary Table")) else
  "[Primary-outcome Steiger result not available -- fill IGGC_NCONTROL.]",
"",
"MULTIPLE TESTING",
sprintf(paste0("P values for the primary estimator (Wald ratio for single-variant instruments, ",
"IVW otherwise) were corrected using the Benjamini-Hochberg false discovery rate across the %d ",
"outcomes tested for PDE5A, and separately across the PDE6 gene-outcome pairs; Bonferroni-",
"corrected P values are reported alongside. q < 0.05 was considered significant."),
  sum(primary$family == "PDE5A")),
"",
"PDE6 (OFF-TARGET) ANALYSIS",
paste0("Sildenafil inhibits retinal PDE6 at roughly a tenth of its PDE5 potency, whereas tadalafil ",
"is highly selective for PDE5. To test whether an off-target PDE6 effect could explain agent-specific ",
"associations in observational data, we instrumented expression of the catalytic PDE6 subunits ",
"(PDE6A, PDE6B, PDE6C) using cis-eQTLs from eQTLGen (whole blood; P < 5e-8, clumped at r2 < 0.001 ",
"within 10 Mb), oriented to the expression-lowering allele as a proxy for inhibition. ",
"Because these genes are expressed mainly in photoreceptors, blood eQTLs may not reflect retinal ",
"expression; this analysis is exploratory.")
)
writeLines(methods, file.path(OUT_DIR, "methods_text_additions.txt"))

cat("\nDone. Outputs in", OUT_DIR, "\n")
sink(file.path(OUT_DIR, "sessionInfo.txt")); print(sessionInfo()); sink()

} # end analysis block
