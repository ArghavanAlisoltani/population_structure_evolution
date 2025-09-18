## ==============================================================
## Most-distinct k selection + certification (PC-space distances)
## ==============================================================
setwd("~/Desktop/OSU_projects/conifers/LP/sample_selection/ARGweaver_sample_selection/final_script3")
## ---------------- USER SETTINGS ----------------
k_target  <- 100          # how many to select (e.g., 100, 200)
pc_rank   <- 30           # number of PCs to use for distance
min_sep_q <- 0.05         # initial min-sep quantile (e.g., 0.05 = 5th pct of all pairwise distances)
relax     <- 0.85         # relax factor for adaptive min-sep
seed      <- 1            # RNG seed (only used for optional pool step)

## ---------------- 0) READ GENOTYPES + PCs ----------------
# Expect a plain numeric table. If your first column is IDs, read with row.names=1.
snp<-read.table("../../snp_1490x25099_CharlesSEP2025.txt")
snp[1:5,1:5]
snpm= as.matrix(snp)
snpm[1:5,1:5]
var_snp <- apply(snpm, 2, var, na.rm = TRUE)
X_filtered <- snpm[, var_snp > 0]
X <- as.matrix(X_filtered)

# PCs
Xc <- scale(X, center = TRUE, scale = TRUE)
pc <- prcomp(Xc, center = FALSE, scale. = FALSE, rank. = pc_rank)
S  <- pc$x[, 1:pc_rank, drop = FALSE]

# Distance (Euclidean in PC space)
D  <- as.matrix(dist(S))
diag(D) <- 0
ids <- rownames(D)   # assign IDs if none
rownames(D) <- colnames(D) <- ids

## ---------------- 1) GREEDY MAX–MIN (ADAPTIVE MIN-SEP) ----------------
greedy_maxmin_adaptive <- function(D, k, min_sep_init = 0, relax = 0.8, min_sep_floor = 0) {
  ids <- rownames(D)
  diag(D) <- NA
  # Seed = most isolated (largest mean distance)
  seed <- ids[ which.max(rowMeans(D, na.rm = TRUE)) ]
  sel  <- seed
  nn   <- D[, seed]; names(nn) <- ids
  
  min_sep <- min_sep_init
  while (length(sel) < k) {
    cand <- setdiff(ids, sel)
    if (!length(cand)) break
    
    cand_ok <- cand[ nn[cand] >= (min_sep - 1e-12) ]
    if (!length(cand_ok)) {
      if (min_sep <= min_sep_floor + .Machine$double.eps) {
        j <- cand[ which.max(nn[cand]) ]
        sel <- c(sel, j)
        nn  <- pmin(nn, D[, j])
      } else {
        min_sep <- max(min_sep * relax, min_sep_floor)
      }
      next
    }
    
    j <- cand_ok[ which.max(nn[cand_ok]) ]
    sel <- c(sel, j)
    nn  <- pmin(nn, D[, j])
  }
  sel
}

## ---------------- 2) 1-SWAP POLISH (MAX–MIN) ----------------
polish_1swap_maxmin <- function(D, sel, max_iter = 50) {
  ids <- rownames(D)
  for (it in 1:max_iter) {
    ii <- match(sel, ids)
    M  <- D[ii, ii]; diag(M) <- NA
    cur <- suppressWarnings(min(M, na.rm = TRUE))
    improved <- FALSE
    uns <- setdiff(ids, sel)
    for (s in sel) {
      sm <- setdiff(sel, s)
      for (t in uns) {
        S2 <- c(sm, t)
        jj <- match(S2, ids)
        Mc <- D[jj, jj]; diag(Mc) <- NA
        cand <- suppressWarnings(min(Mc, na.rm = TRUE))
        if (is.finite(cand) && cand > cur + 1e-12) {
          sel <- S2; improved <- TRUE; break
        }
      }
      if (improved) break
    }
    if (!improved) break
  }
  sel
}

## ---------------- 3) RUN SELECTION ----------------
all_d   <- D[upper.tri(D)]
eps0    <- as.numeric(quantile(all_d, min_sep_q, na.rm = TRUE))
sel0    <- greedy_maxmin_adaptive(D, k = k_target, min_sep_init = eps0, relax = relax, min_sep_floor = 0)
sel_fin <- polish_1swap_maxmin(D, sel0, max_iter = 100)

cat("Selected:", length(sel_fin), "IDs\n")
M <- D[sel_fin, sel_fin, drop = FALSE]; diag(M) <- NA
r_set <- suppressWarnings(min(M, na.rm = TRUE))
cat("Greedy + polish radius (min pairwise distance):", r_set, "\n")

library(mixOmics)
meta <- read.csv("../../metadata.matchesCappa_GRM.csv", stringsAsFactors = FALSE)
snp<-read.table("../../snp_1490x25099_CharlesSEP2025.txt")
snp[1:5,1:5]
snpm= as.matrix(snp)
snpm[1:5,1:5]
var_snp <- apply(snpm, 2, var, na.rm = TRUE)
X_filtered <- snpm[, var_snp > 0]
X_scaled <- scale(X_filtered, center = TRUE, scale = TRUE)
pca.res<- pca(X_scaled,center = F, ncomp =3, scale = F)
df<-data.frame(IDs=row.names(snp), ordersnp=c(1:1490))
merged<-merge(df, meta, by.x="IDs",by.y= "name", all.x=T)
merged<-merged[order(merged$ordersnp),]

merged$selected<-ifelse(merged$IDs %in% sel_fin,
                          "selected", "not_selected" )

plotIndiv(
  pca.res,
  group = merged$selected,
  ind.names = rownames(X_scaled),   # or rownames(snpm)
  label = TRUE,
  col = c("red","gray"),
  legend = TRUE, legend.title = "selected"
)
merged$selected<-ifelse(merged$IDs %in% sel_fin,
                        "selected", "not_selected" )
tbl<-data.frame(table(merged$mum,merged$selected))
tbl1<-data.frame(table(merged$proc,merged$selected))
tbl2<-data.frame(table(merged$site,merged$selected))

plotIndiv(
  pca.res,
  group = merged$selected,
  ind.names = FALSE,
  style = "ggplot2",
  col = c("gray","red"),
  pch = 16,    # still works
  cex = 1.4,
  legend = TRUE,
  legend.title = " "
)
merged_filt<-merged[merged$selected=="selected",]
library(data.table)
fwrite(merged_filt,paste("merged",k_target,".tsv",sep = ""),sep="\t")

## ---------------- 4) MILP FEASIBILITY (LAZY CUTS, GLPK) ----------------
suppressPackageStartupMessages({
  library(ompr); library(ompr.roi); library(ROI); library(ROI.plugin.glpk)
})

edge_count <- function(D, r, tol = 1e-9) sum(D[upper.tri(D)] < (r - tol))

# exact feasibility at radius r (time-limited)
feasible_kcenter_radius_exact <- function(D, k, r, tol = 1e-9, tm_ms = 120000, verbose = FALSE) {
  n <- nrow(D)
  A <- which(D < (r - tol) & upper.tri(D), arr.ind = TRUE)
  model <- MIPModel()
  model <- add_variable(model, x[i], i = 1:n, type = "binary")
  model <- add_constraint(model, sum_expr(x[i], i = 1:n) == k)
  if (nrow(A)) model <- add_constraint(model, x[A[c,1]] + x[A[c,2]] <= 1, c = 1:nrow(A))
  model <- set_objective(model, 0)
  
  res <- solve_model(model, with_ROI(solver = "glpk",
                                     control = list(tm_limit = tm_ms, presolve = TRUE),
                                     verbose = verbose))
  list(ok = solver_status(res) %in% c("optimal","feasible"),
       status = solver_status(res))
}

# Pick an "easy" infeasible upper bound: step r up until conflict edges exceed a target
pick_easy_hi <- function(D, r_start, step = 1.0, target_edges = 2e5, max_steps = 50) {
  r <- r_start
  for (s in 1:max_steps) {
    r <- r + step
    if (edge_count(D, r) >= target_edges) return(r)
  }
  r
}

# Binary search with time-limited exact feasibility
kcenter_certify_exact <- function(D, k, lo, hi, tol = 1e-3, tm_ms = 120000) {
  best <- lo
  for (it in 1:35) {
    mid <- (lo + hi) / 2
    ans <- feasible_kcenter_radius_exact(D, k, mid, tm_ms = tm_ms)
    if (ans$ok) { best <- mid; lo <- mid } else { hi <- mid }
    if ((hi - lo) < tol) break
  }
  best
}

# ---- Use with your selection ----
sel_idx <- match(sel_fin, rownames(D))
M <- D[sel_fin, sel_fin, drop = FALSE]; diag(M) <- NA
r_set <- min(M, na.rm = TRUE)

# Fast sanity: only the FEASIBLE side
below <- feasible_kcenter_radius_exact(D, length(sel_fin), r_set - 1e-6, tm_ms = 60000)
cat("Sanity below r_set: ", below$ok, " (status=", below$status, ")\n", sep = "")

# Find an easy infeasible upper bound (far above boundary)
hi_easy <- pick_easy_hi(D, r_start = r_set, step = 1.0, target_edges = 2e5)
cat("Chosen easy upper bound (hi): ", sprintf("%.4f", hi_easy),
    " with edges=", edge_count(D, hi_easy), "\n", sep = "")

# Now certify r* by binary search between r_set*0.8 and hi_easy (time-limited)
r_star <- kcenter_certify_exact(D, k = length(sel_fin),
                                lo = max(0, r_set * 0.8),
                                hi = hi_easy,
                                tol = 1e-3,
                                tm_ms = 60000)
cat("Certified r* (time-limited exact) ≈ ", sprintf("%.4f", r_star),
    "   (your set radius = ", sprintf("%.4f", r_set), ")\n", sep = "")
#########################################################
# r_set: your achieved min pairwise distance
#########################################################
M <- D[sel_fin, sel_fin]; diag(M) <- NA
r_set <- min(M, na.rm = TRUE)

set.seed(1)
B <- 300
r_rand <- replicate(B, {
  idx <- sample(rownames(D), length(sel_fin))
  M <- D[idx, idx]; diag(M) <- NA
  suppressWarnings(min(M, na.rm = TRUE))
})

cat("Your radius:", r_set, "\n",
    "Random mean radius:", mean(r_rand), "\n",
    "Random 95th percentile:", quantile(r_rand, 0.95), "\n",
    "Fraction random < yours:", mean(r_rand < r_set), "\n")

# Coverage check: how well selected points "cover" everyone
nn_selected <- apply(D[, sel_fin, drop = FALSE], 1, min)  # distance of each sample to nearest selected
cover_radius <- max(nn_selected)                           # smaller = better coverage
cat("Coverage radius (max nearest-to-selected):", cover_radius, "\n")

