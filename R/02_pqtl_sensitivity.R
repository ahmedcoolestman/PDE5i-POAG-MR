###############################################################################
# 02_pqtl_sensitivity.R
#
# Sensitivity analysis: PDE5A plasma protein (Suhre et al. 2017, KORA F4,
# OpenGWAS prot-c-5256_86_3) as the exposure, cis-window PDE5A +/- 150 kb,
# against the same six outcomes as the main analysis. Estimates are per 1 SD
# LOWER protein. Also runs the PheWAS lookup for rs2389873 (Supp. Table S5)
# and prints the strand-sensitivity and Steiger outputs from script 01
# (Supp. Tables S2, S3).
#
# Run after 01_main_analysis.R, from the same folder:
#     Rscript R/02_pqtl_sensitivity.R
# Outputs go to results_pqtl/.
###############################################################################

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(ieugwasr)
  library(data.table)
  library(dplyr)
})

OUT_DIR <- "results_pqtl"
dir.create(OUT_DIR, showWarnings = FALSE)

# ---------------------------------------------------------------- settings ---
PQTL_ID     <- "prot-c-5256_86_3"
PQTL_N      <- 1000                      # KORA F4 discovery (Suhre 2017)
REG_START   <- 120265550                 # PDE5A +/- 150 kb, GRCh37 (as main analysis)
REG_END     <- 120700146
REGION      <- paste0("4:", REG_START, "-", REG_END)
P_MAIN      <- 5e-8                      # preferred threshold
P_RELAXED   <- 5e-6                      # used only if nothing passes P_MAIN
CLUMP_R2    <- 0.001
CLUMP_KB    <- 250
MAIN_SNP    <- "rs2389873"               # main-analysis instrument (v4)

IGGC_FILE     <- "GCST90011767_buildGRCh37.tsv.gz"
IGGC_NCASE    <- 15229
IGGC_NCONTROL <- 177473

OUTCOMES <- data.frame(
  id    = c("IGGC_LOCAL", "finn-b-H7_GLAUCOMA_POAG", "finn-b-H7_GLAUCPRIMOPEN",
            "finn-b-H7_GLAUCOMA", "finn-b-H7_GLAUCOMA_NTG", "ebi-a-GCST004074"),
  label = c("POAG (IGGC, primary)", "POAG strict (FinnGen)", "POAG broad (FinnGen)",
            "All glaucoma (FinnGen)", "NTG (FinnGen)", "IOP"),
  type  = c(rep("binary", 5), "continuous"),
  stringsAsFactors = FALSE
)

is_palindromic <- function(a1, a2) {
  a1 <- toupper(a1); a2 <- toupper(a2)
  (a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") |
    (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C")
}

cat("\n==== PDE5A pQTL sensitivity analysis ====\n")

# ----------------------------------------------- 1. cis-pQTL instruments ---
cat("\n-- Extracting PDE5A region from", PQTL_ID, "\n")
reg <- tryCatch(associations(variants = REGION, id = PQTL_ID, proxies = 0),
                error = function(e) { message("OpenGWAS error: ", conditionMessage(e)); NULL })
if (is.null(reg) || nrow(reg) == 0) {
  cat("Region query returned nothing; falling back to tophits on chr4\n")
  reg <- tryCatch(tophits(PQTL_ID, pval = P_RELAXED, clump = 0), error = function(e) NULL)
  if (!is.null(reg)) reg <- reg[as.character(reg$chr) == "4" &
                                  reg$position >= REG_START & reg$position <= REG_END, ]
}
if (is.null(reg) || nrow(reg) == 0) stop("No PDE5A-region data returned for ", PQTL_ID)
reg <- as.data.frame(reg)
cat("SNPs in region:", nrow(reg), "| min P:", signif(min(reg$p), 3), "\n")

sig <- reg[reg$p < P_MAIN, ]
P_USED <- P_MAIN
if (nrow(sig) == 0) {
  sig <- reg[reg$p < P_RELAXED, ]
  P_USED <- P_RELAXED
  cat("No cis-pQTL at P <", P_MAIN, "- using P <", P_RELAXED, "\n")
}
if (nrow(sig) == 0) stop("No cis-pQTL at P < ", P_RELAXED)
sig <- sig[order(sig$p), ]
cat("Candidate cis-pQTLs at P <", P_USED, ":", nrow(sig), "\n")

# Clump, then double-check with pairwise LD (the clumping server kept two
# SNPs with r2 = 0.80 in the main analysis, so don't rely on it alone)
cl <- tryCatch(ld_clump(data.frame(rsid = sig$rsid, pval = sig$p),
                        clump_kb = CLUMP_KB, clump_r2 = CLUMP_R2, pop = "EUR"),
               error = function(e) NULL)
keep <- if (!is.null(cl) && nrow(cl) > 0) cl$rsid else sig$rsid[1]
keep <- sig$rsid[sig$rsid %in% keep]
if (length(keep) > 1) {
  ldm <- tryCatch(ld_matrix(keep, pop = "EUR"), error = function(e) NULL)
  if (!is.null(ldm)) {
    nm <- sub("_.*", "", rownames(ldm)); r2 <- ldm^2; dimnames(r2) <- list(nm, nm)
    final <- character(0)
    for (s in keep) {
      if (!s %in% nm) next
      if (length(final) == 0 || all(r2[s, final] < 0.01)) final <- c(final, s)
    }
    keep <- final
  } else keep <- keep[1]
}
ins_raw <- sig[sig$rsid %in% keep, ]
ins_raw$pal <- is_palindromic(ins_raw$ea, ins_raw$nea)
cat("Independent cis-pQTLs:", nrow(ins_raw), "\n")

if (!"n" %in% names(ins_raw) || all(is.na(ins_raw$n))) ins_raw$n <- PQTL_N
ins_raw$n[is.na(ins_raw$n)] <- PQTL_N
ins_raw$phen <- "PDE5A protein (lower = inhibition)"
exp_dat <- format_data(ins_raw, type = "exposure", snp_col = "rsid",
                       beta_col = "beta", se_col = "se", eaf_col = "eaf",
                       effect_allele_col = "ea", other_allele_col = "nea",
                       pval_col = "p", samplesize_col = "n", phenotype_col = "phen")
exp_dat$id.exposure <- PQTL_ID

# Orient to the protein-LOWERING allele (proxy for PDE5 inhibition)
flip <- exp_dat$beta.exposure > 0
if (any(flip)) {
  ea <- exp_dat$effect_allele.exposure[flip]
  exp_dat$effect_allele.exposure[flip] <- exp_dat$other_allele.exposure[flip]
  exp_dat$other_allele.exposure[flip]  <- ea
  exp_dat$beta.exposure[flip] <- -exp_dat$beta.exposure[flip]
  exp_dat$eaf.exposure[flip]  <- 1 - exp_dat$eaf.exposure[flip]
}
exp_dat$F_stat <- (exp_dat$beta.exposure / exp_dat$se.exposure)^2
inst_out <- exp_dat[, c("SNP", "effect_allele.exposure", "other_allele.exposure",
                        "eaf.exposure", "beta.exposure", "se.exposure",
                        "pval.exposure", "F_stat")]
inst_out$is_palindromic <- is_palindromic(inst_out$effect_allele.exposure,
                                          inst_out$other_allele.exposure)
cat("\nInstrument(s), oriented to protein-lowering allele (beta in SD units):\n")
print(inst_out, row.names = FALSE)
fwrite(inst_out, file.path(OUT_DIR, "S4a_pqtl_instruments.csv"))

# Is the pQTL the same signal as the BP instrument?
ldchk <- tryCatch(ld_matrix(unique(c(exp_dat$SNP, MAIN_SNP)), pop = "EUR"),
                  error = function(e) NULL)
if (!is.null(ldchk)) {
  cat("\nLD (r2) between pQTL instrument(s) and", MAIN_SNP, ":\n")
  print(round(ldchk^2, 3))
}

# ------------------------------------------------------------ 2. outcomes ---
read_iggc <- function(snps) {
  pat <- paste0("variant_id|", paste(snps, collapse = "|"))
  d <- tryCatch(fread(cmd = sprintf("zgrep -wE '%s' %s", pat, IGGC_FILE)),
                error = function(e) NULL)
  if (is.null(d) || !"variant_id" %in% names(d)) d <- fread(IGGC_FILE)
  d <- d[variant_id %in% snps]
  if (nrow(d) == 0) return(NULL)
  d[, `:=`(ncase = IGGC_NCASE, ncontrol = IGGC_NCONTROL,
           samplesize = IGGC_NCASE + IGGC_NCONTROL, phen = "POAG (IGGC)")]
  format_data(as.data.frame(d), type = "outcome", snp_col = "variant_id",
              beta_col = "beta", se_col = "standard_error",
              effect_allele_col = "effect_allele", other_allele_col = "other_allele",
              pval_col = "p_value", ncase_col = "ncase", ncontrol_col = "ncontrol",
              samplesize_col = "samplesize", phenotype_col = "phen")
}
get_outcome <- function(snps, oid) {
  if (oid == "IGGC_LOCAL") return(read_iggc(snps))
  d <- tryCatch(extract_outcome_data(snps = snps, outcomes = oid, proxies = TRUE, rsq = 0.8),
                error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(NULL)
  d
}
has_eaf <- function(o) "eaf.outcome" %in% names(o) && any(!is.na(o$eaf.outcome))

res <- list(); used <- list()
for (i in seq_len(nrow(OUTCOMES))) {
  o <- OUTCOMES[i, ]
  out <- get_outcome(exp_dat$SNP, o$id)
  if (is.null(out)) { message("No outcome data: ", o$label); next }
  act <- if (has_eaf(out)) 2 else 3
  h <- harmonise_data(exp_dat, out, action = act)
  h <- h[h$mr_keep, , drop = FALSE]
  if (nrow(h) == 0) { message("Nothing survived harmonisation: ", o$label); next }
  m <- if (nrow(h) == 1) "mr_wald_ratio" else "mr_ivw"
  r <- mr(h, method_list = m)
  proxy <- if ("proxy.outcome" %in% names(h)) any(h$proxy.outcome %in% TRUE) else FALSE
  # Wald ratio = per unit HIGHER exposure; negate so estimates are per SD LOWER protein
  res[[o$label]] <- data.frame(outcome = o$label, type = o$type, method = r$method,
                               nsnp = r$nsnp, b = -r$b, se = r$se, pval = r$pval,
                               proxy_used = proxy)
}
res <- bind_rows(res)
res$q_fdr_bh     <- p.adjust(res$pval, method = "BH")
res$p_bonferroni <- p.adjust(res$pval, method = "bonferroni")
res$estimate <- ifelse(res$type == "binary", exp(res$b), res$b)
res$lo       <- ifelse(res$type == "binary", exp(res$b - 1.96 * res$se), res$b - 1.96 * res$se)
res$hi       <- ifelse(res$type == "binary", exp(res$b + 1.96 * res$se), res$b + 1.96 * res$se)
res$scale    <- ifelse(res$type == "binary", "OR per SD lower PDE5A protein",
                       "mmHg IOP per SD lower PDE5A protein")

cat("\n==== pQTL results (Supp. Table S4) ====\n")
cat("Threshold used for instrument selection: P <", P_USED, "\n")
print(res[, c("outcome", "method", "nsnp", "estimate", "lo", "hi", "pval",
              "q_fdr_bh", "p_bonferroni", "proxy_used", "scale")],
      row.names = FALSE, digits = 3)
fwrite(res, file.path(OUT_DIR, "S4b_pqtl_results.csv"))

# ------------------------------------------------- 3. PheWAS (Supp. S5) ---
cat("\n==== PheWAS for", MAIN_SNP, "(Supp. Table S5) ====\n")
ph <- tryCatch(as.data.frame(phewas(MAIN_SNP, pval = 5e-8)), error = function(e) NULL)
if (!is.null(ph) && nrow(ph)) {
  genes <- c(ENSG00000164096 = "ARLN", ENSG00000145390 = "USP53",
             ENSG00000138735 = "PDE5A", ENSG00000150961 = "SEC24D",
             ENSG00000245958 = "SEPTIN7P14", ENSG00000260091 = "lncRNA (unnamed)",
             ENSG00000248280 = "retired Ensembl ID")
  ph$gene <- ifelse(ph$trait %in% names(genes), genes[ph$trait], "")
  keepc <- intersect(c("trait", "gene", "id", "ea", "nea", "beta", "se", "p", "n"), names(ph))
  ph <- ph[order(ph$p), keepc]
  print(ph, row.names = FALSE)
  fwrite(ph, file.path(OUT_DIR, "S5_phewas_rs2389873.csv"))
}

# ----------------------------------- 4. existing S2 / S3 (for the tables) ---
for (f in c("results_v4/A_palindromic_strand_sensitivity.csv",
            "results_v4/B_steiger_directionality.csv")) {
  cat("\n====", f, "====\n")
  if (file.exists(f)) print(fread(f), digits = 3) else cat("not found\n")
}

cat("\nDone. Outputs in", OUT_DIR, "\n")
