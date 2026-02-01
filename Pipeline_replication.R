library(data.table)
f1 <- "data_raw/34017140-GCST90013862-EFO_0005842.h.tsv.gz"
f2 <- "data_raw/GCST90475566.h.tsv.gz"
f3 <- "data_raw/GCST90255675.h.tsv.gz"

dir.create("data_clean", showWarnings = FALSE)

# Quick sanity check inputs
peek_cols <- function(f) names(fread(f, nrows = 1))
cat("\nschema:\n")
cat("f1 cols:", length(peek_cols(f1)), "\n")
cat("f2 cols:", length(peek_cols(f2)), "\n")
cat("f3 cols:", length(peek_cols(f3)), "\n")

# Chromosome and allele harmonization
clean_chr <- function(x){
  x <- gsub("^chr", "", as.character(x), ignore.case = TRUE)
  x <- toupper(x)
  x[x == "X"] <- "23"
  x[x == "Y"] <- "24"
  x[x %in% c("MT","M")] <- "25"
  suppressWarnings(as.integer(x))
}

is_snp_allele <- function(a) a %chin% c("A","C","G","T")

is_pal <- function(a,b){
  (a=="A" & b=="T") | (a=="T" & b=="A") | (a=="C" & b=="G") | (a=="G" & b=="C")
}

#Harmonize GWAS summary statistics
read_one <- function(path, prefer_hm = TRUE) {
  dt <- fread(path)
  
  if (prefer_hm && all(c("hm_chrom","hm_pos","hm_effect_allele","hm_other_allele") %in% names(dt))) {
    chr_col <- "hm_chrom"; pos_col <- "hm_pos"
    ea_col  <- "hm_effect_allele"; oa_col <- "hm_other_allele"
    src <- "hm_*"
  } else if (all(c("chromosome","base_pair_location","effect_allele","other_allele") %in% names(dt))) {
    chr_col <- "chromosome"; pos_col <- "base_pair_location"
    ea_col  <- "effect_allele"; oa_col <- "other_allele"
    src <- "plain"
  } else {
    stop("No chr/pos/alleles columns in: ", path)
  }
  
  setnames(dt, old=c(chr_col,pos_col,ea_col,oa_col), new=c("chr","pos","EA","OA"))
  dt[, chr := clean_chr(chr)]
  dt[, pos := suppressWarnings(as.integer(pos))]
  dt <- dt[!is.na(chr) & !is.na(pos)]
  
  dt[, EA := toupper(as.character(EA))]
  dt[, OA := toupper(as.character(OA))]
  dt <- dt[is_snp_allele(EA) & is_snp_allele(OA) & EA != OA]
  
  if (all(c("beta","standard_error") %in% names(dt))) {
    dt[, beta := suppressWarnings(as.numeric(beta))]
    dt[, se   := suppressWarnings(as.numeric(standard_error))]
    dt <- dt[!is.na(beta) & !is.na(se) & se > 0]
    dt[, z := beta / se]
    zsrc <- "beta/se"
    
  } else if (all(c("odds_ratio","ci_upper","ci_lower") %in% names(dt))) {
    dt[, OR   := suppressWarnings(as.numeric(odds_ratio))]
    dt[, ci_u := suppressWarnings(as.numeric(ci_upper))]
    dt[, ci_l := suppressWarnings(as.numeric(ci_lower))]
    dt <- dt[!is.na(OR) & OR > 0 & !is.na(ci_u) & ci_u > 0 & !is.na(ci_l) & ci_l > 0]
    
    dt[, se := (log(ci_u) - log(ci_l)) / (2 * 1.96)]
    dt <- dt[!is.na(se) & se > 0]
    dt[, z := log(OR) / se]
    zsrc <- "logOR/SE(from CI)"
    
  } else if (all(c("odds_ratio","standard_error") %in% names(dt))) {
    dt[, OR := suppressWarnings(as.numeric(odds_ratio))]
    dt[, se := suppressWarnings(as.numeric(standard_error))]
    dt <- dt[!is.na(OR) & OR > 0 & !is.na(se) & se > 0]
    dt[, z := log(OR) / se]
    zsrc <- "logOR/se(standard_error)"
    
  } else {
    stop("No usable effect+uncertainty columns in: ", path)
  }
  
  if ("variant_id" %in% names(dt)) dt[, variant_id := as.character(variant_id)] else dt[, variant_id := NA_character_]
  if ("rsid" %in% names(dt)) dt[, rsid := as.character(rsid)] else dt[, rsid := NA_character_]
  
  dt[, key_cp := paste0(chr, ":", pos)]
  
  out <- dt[, .(key_cp, chr, pos, EA, OA, z, variant_id, rsid)]
  attr(out, "pos_source") <- src
  attr(out, "z_source") <- zsrc
  out
}


d1 <- read_one(f1, prefer_hm = TRUE)    
d2 <- read_one(f2, prefer_hm = FALSE)   
d3 <- read_one(f3, prefer_hm = FALSE)   

# Overlap diagnostics across studies
cat("\nread sizes:\n")
cat("d1:", nrow(d1), "pos:", attr(d1,"pos_source"), "z:", attr(d1,"z_source"), "\n")
cat("d2:", nrow(d2), "pos:", attr(d2,"pos_source"), "z:", attr(d2,"z_source"), "\n")
cat("d3:", nrow(d3), "pos:", attr(d3,"pos_source"), "z:", attr(d3,"z_source"), "\n")

k1 <- unique(d1$key_cp); k2 <- unique(d2$key_cp); k3 <- unique(d3$key_cp)
cat("overlap d1&d2:", length(intersect(k1,k2)), "\n")
cat("overlap d1&d3:", length(intersect(k1,k3)), "\n")
cat("overlap d2&d3:", length(intersect(k2,k3)), "\n")
cat("overlap all3 :", length(Reduce(intersect, list(k1,k2,k3))), "\n")

saveRDS(d1, "data_clean/h1_clean.rds")
saveRDS(d2, "data_clean/h2_clean.rds")
saveRDS(d3, "data_clean/h3_clean.rds")

# Tag columns by study index so merges keep study-specific fields distinct
tag_dt <- function(dt, tag){
  setnames(dt,
           old=c("chr","pos","EA","OA","z","variant_id","rsid"),
           new=paste0(c("chr","pos","EA","OA","z","variant_id","rsid"), tag))
  dt
}

d1t <- tag_dt(copy(d1), "1")
d2t <- tag_dt(copy(d2), "2")
d3t <- tag_dt(copy(d3), "3")

# Merge studies on chr:pos key
m <- merge(d1t, d2t, by="key_cp")
cat("\nmerge12:", nrow(m), "\n")
m <- merge(m, d3t, by="key_cp")
cat("merge123:", nrow(m), "\n")

# Harmonize alleles and build the final 3-study Z-score matrix
cat("\nstart:", nrow(m), "\n")

m <- m[!(is_pal(EA1,OA1) | is_pal(EA2,OA2) | is_pal(EA3,OA3))]
cat("after pal filter:", nrow(m), "\n")

ok2_same <- (m$EA2 == m$EA1 & m$OA2 == m$OA1)
ok2_swap <- (m$EA2 == m$OA1 & m$OA2 == m$EA1)
m[ok2_swap, z2 := -z2]
m <- m[ok2_same | ok2_swap]
cat("after align2:", nrow(m), "\n")

ok3_same <- (m$EA3 == m$EA1 & m$OA3 == m$OA1)
ok3_swap <- (m$EA3 == m$OA1 & m$OA3 == m$EA1)
m[ok3_swap, z3 := -z3]
m <- m[ok3_same | ok3_swap]
cat("after align3:", nrow(m), "\n")

Zmat <- as.matrix(m[, .(z1, z2, z3)])
saveRDS(Zmat, "data_clean/Zmat_3study_harmonised.rds")
cat("\nFinal merged SNPs:", nrow(m), "\n")

# Basic diagnostics on Zmat 
dim(Zmat)
sum(!is.finite(Zmat))
summary(as.vector(Zmat))

cor(Zmat[,1], Zmat[,2], use="pairwise.complete.obs")
cor(Zmat[,1], Zmat[,3], use="pairwise.complete.obs")
cor(Zmat[,2], Zmat[,3], use="pairwise.complete.obs")

colMeans(Zmat)
cor(abs(Zmat[,1]), abs(Zmat[,3]))

saveRDS(m, "data_clean/m_3study_harmonised_full.rds")
saveRDS(Zmat, "data_clean/Zmat_3study_harmonised.rds")


library(csmGmm)
ls("package:csmGmm")
help(package = "csmGmm")

# Sanity check
m <- readRDS("data_clean/m_3study_harmonised_full.rds")
cat("m dimensions:", dim(m), "\n")
cat("m column names:\n")
print(names(m))
cat("First few rows:\n")
print(head(m, 3))

# Create standard output directories for results/logs
create_project_directories <- function() {
  dirs <- c("results", "results/figures", "results/tables", "logs")
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
      cat("Created directory:", dir, "\n")
    }
  }
}

# Check basic consistency between Zmat and merged table m
check_data_consistency <- function(Zmat, m) {
  cat("\n=== DATA CONSISTENCY CHECK ===\n")
  cat("Zmat dimensions:", dim(Zmat), "\n")
  cat("m dimensions:", dim(m), "\n")
  
  if (nrow(Zmat) != nrow(m)) {
    stop("ERROR: Zmat and m have different number of rows!")
  }
  
  cat("Missing values in Zmat:", sum(is.na(Zmat)), "\n")
  cat("Missing values in m:", sum(is.na(m)), "\n")
  
  cat("\nZ-score summary:\n")
  print(summary(as.vector(Zmat)))
  
  same_sign <- apply(Zmat, 1, function(x) all(x > 0) | all(x < 0))
  cat("\nSNPs with consistent direction in all 3 studies:", 
      round(mean(same_sign) * 100, 2), "%\n")
  cat("All positive:", round(mean(apply(Zmat, 1, function(x) all(x > 0))) * 100, 2), "%\n")
  cat("All negative:", round(mean(apply(Zmat, 1, function(x) all(x < 0))) * 100, 2), "%\n")
  
  return(TRUE)
}

create_directories <- function() {
  dirs <- c("results", "results/figures", "results/tables", "logs")
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
      cat("Created directory:", dir, "\n")
    }
  }
}

dir.create("results", showWarnings=FALSE)
dir.create("results/figures", recursive=TRUE, showWarnings=FALSE)
dir.create("results/tables", recursive=TRUE, showWarnings=FALSE)

Zmat <- readRDS("data_clean/Zmat_3study_harmonised.rds")
m    <- readRDS("data_clean/m_3study_harmonised_full.rds")

stopifnot(nrow(Zmat) == nrow(m), ncol(Zmat) == 3)

# Create a smaller subset for a quick EM test run
test_size <- min(10000, nrow(Zmat))
Ztest <- Zmat[1:test_size, ]

# Initialize symmetric mixture component means for 3-study setting
initMuList <- list(
  matrix(0, 3, 1),  # b=0: (0,0,0)
  
  cbind(c(0,0,2.0),  c(0,0,-2.0)),      # b=1: (0,0,±)
  cbind(c(0,2.0,0),  c(0,-2.0,0)),      # b=2: (0,±,0)
  cbind(c(0,2.5,2.5),c(0,-2.5,-2.5)),   # b=3: (0,±,±)
  
  cbind(c(2.0,0,0),  c(-2.0,0,0)),      # b=4: (±,0,0)
  cbind(c(2.5,0,2.5),c(-2.5,0,-2.5)),   # b=5: (±,0,±)
  cbind(c(2.5,2.5,0),c(-2.5,-2.5,0)),   # b=6: (±,±,0)
  
  cbind(c(3.0,3.0,3.0), c(-3.0,-3.0,-3.0)) # b=7: (±,±,±) 
)

# Initialize mixture weights per component set
initPiList <- list(
  c(0.85),
  c(0.0075, 0.0075),
  c(0.0075, 0.0075),
  c(0.005,  0.005),
  c(0.0075, 0.0075),
  c(0.005,  0.005),
  c(0.005,  0.005),
  c(0.01,   0.01)
)

for (i in 1:8) stopifnot(ncol(initMuList[[i]]) == length(initPiList[[i]]))

# Fit csmGmm symmetric independent EM on the test subset
cat("Test run on", test_size, "SNPs\n")
fit_test <- symm_fit_ind_EM(
  testStats = Ztest,
  initMuList = initMuList,
  initPiList = initPiList,
  sameDirAlt = TRUE,
  eps = 1e-4,
  checkpoint = TRUE
)
saveRDS(fit_test, "results/fit_test_10k.rds")
cat("Test done. iter=", fit_test$iter, " lfdr range=", range(fit_test$lfdrResults), "\n")
