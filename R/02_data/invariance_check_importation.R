# =====================================================================
# invariance_check_importation.R   (base R, dependency-free)
#
# Purpose: prove that switching the importation formula from
#   OLD (K1):  V_in_day = V_TR_d × N_total × S_m / 365
#   NEW (M1):  V_in_day = A_d × S_m / 365
# is a PURE COMMON RESCALING (a single global scalar c) and therefore
# leaves the district RANKING, threshold-crossing years and season
# lengths UNCHANGED; only the absolute P_ufuk scale shifts.
#
# Algebraic core:
#   V_in_new / V_in_old = A_d / (V_TR_d × N_total)
#                       = A_d / ((A_d / Σ_j A_j) × N_total)
#                       = (Σ_j A_j) / N_total  =  c     (independent of d)
#
# Since c > 0 is identical for every district/month/scenario, and
#   P_ufuk = 1 - exp(-Λ·P_est)  is strictly increasing in Λ,
# multiplying all Λ_d by the same c preserves the ordering of {Λ_d·P_est_d}.
# Downstream, using published P_ufuk_old:
#   P_ufuk_new = 1 - (1 - P_ufuk_old)^c        (common exponent c ∈ (0,1))
# and x ↦ 1 - (1-x)^c is strictly increasing on (0,1) ⟹ rank preserved.
#
# P_est, threshold year (R0=1) and season length contain NO importation
# term and are therefore literally identical under both formulas.
# =====================================================================

options(stringsAsFactors = FALSE)

# ---------------------------------------------------------------------
# CONFIG — point these at your project, or leave defaults to run standalone
# ---------------------------------------------------------------------
TW_RDS    <- Sys.getenv("TW_RDS",  unset = "data_processed/travel_weights_static.rds")
N_TOTAL   <- as.numeric(Sys.getenv("N_TOTAL", unset = "40900000"))  # Σ mean_arrivals (~40.9M)
TOL       <- 1e-9

# ---------------------------------------------------------------------
# PART A — load A_d (raw district foreign arrivals) and prove ratio = c
# ---------------------------------------------------------------------
if (file.exists(TW_RDS)) {
  tw <- readRDS(TW_RDS)
  stopifnot("gelis_yabanci" %in% names(tw), "travel_weight" %in% names(tw))
  dsts <- data.frame(
    district = if ("district_name" %in% names(tw)) tw$district_name else tw$district_id,
    A_d      = as.numeric(tw$gelis_yabanci),
    V_TR     = as.numeric(tw$travel_weight)
  )
  cat("A_d kaynagi: ", TW_RDS, "\n", sep = "")
} else {
  # Standalone fallback (illustrative TÜİK-2022-style values; shares ~ audit read).
  cat("[uyari] ", TW_RDS, " bulunamadi — ORNEK A_d ile calisiyor.\n", sep = "")
  dsts <- data.frame(
    district = c("Fethiye","Kartal","Hopa","Egirdir","Zonguldak"),
    A_d      = c(1000000,  197000, 76000, 1150,     760),
    V_TR     = NA_real_
  )
  dsts$V_TR <- dsts$A_d / sum(dsts$A_d)
}

sumA <- sum(dsts$A_d)
c_scale <- sumA / N_TOTAL

dsts$V_in_old <- dsts$V_TR * N_TOTAL      # eski (K1): 5-ilce payi × ulusal toplam
dsts$V_in_new <- dsts$A_d                 # yeni (M1): dogrudan ham ilce gelisi
dsts$ratio    <- dsts$V_in_new / dsts$V_in_old

cat("\n================ PART A: ortak olcek =================\n")
cat(sprintf("Sigma A_d (5 ilce)   = %s\n", format(sumA, big.mark = ",")))
cat(sprintf("N_total (ulusal)     = %s\n", format(N_TOTAL, big.mark = ",")))
cat(sprintf("c = SigmaA / N_total = %.6f   (Lambda ~%.1f kat kuculur)\n\n", c_scale, 1/c_scale))
print(dsts[, c("district","A_d","V_TR","ratio")], row.names = FALSE)

ratio_constant <- max(abs(dsts$ratio - c_scale)) < TOL
cat(sprintf("\n[TEST A] V_in(NEW)/V_in(OLD) tum ilcelerde sabit (= c)?  %s\n",
            ifelse(ratio_constant, "GECTI", "KALDI")))
if (!ratio_constant) stop("Oran sabit degil — saf yeniden-olcekleme bozulmus.", call. = FALSE)

# ---------------------------------------------------------------------
# PART B — apply the common rescaling to published P_ufuk and check rank
#   Replace `pufuk` below with your real results table if desired:
#     pufuk <- read.csv("results/pufuk_50yr.csv")  # cols: district, ssp, P_ufuk_old
# Seeded here with Tablo 6.5.1.1 (tezim.docx) values so it runs standalone.
# ---------------------------------------------------------------------
pufuk <- data.frame(
  district = rep(c("Egirdir","Fethiye","Hopa","Kartal","Zonguldak"), each = 3),
  ssp      = rep(c("ssp126","ssp245","ssp585"), times = 5),
  P_ufuk_old = c(
    5.58e-9, 8.40e-8, 1.97e-8,   # Egirdir
    1.42e-4, 3.55e-3, 3.96e-4,   # Fethiye
    1.00e-2, 7.11e-2, 1.50e-1,   # Hopa
    2.67e-1, 3.29e-1, 4.46e-1,   # Kartal
    9.14e-5, 4.50e-4, 1.04e-3    # Zonguldak
  )
)

# Common-scalar transform (guaranteed rank-preserving for c in (0,1)):
pufuk$P_ufuk_new <- 1 - (1 - pufuk$P_ufuk_old)^c_scale

cat("\n================ PART B: siralama =================\n")
rank_ok_all <- TRUE
for (s in unique(pufuk$ssp)) {
  sub <- pufuk[pufuk$ssp == s, ]
  r_old <- sub$district[order(-sub$P_ufuk_old)]
  r_new <- sub$district[order(-sub$P_ufuk_new)]
  ok <- identical(r_old, r_new)
  rank_ok_all <- rank_ok_all && ok
  cat(sprintf("\n[%s]  siralama ozdes mi? %s\n", s, ifelse(ok, "EVET", "HAYIR")))
  print(data.frame(
    district   = sub$district,
    P_old      = formatC(sub$P_ufuk_old, format = "e", digits = 2),
    P_new_M1   = formatC(sub$P_ufuk_new, format = "e", digits = 2)
  )[order(-sub$P_ufuk_old), ], row.names = FALSE)
}
cat(sprintf("\n[TEST B] Tum SSP'lerde ilce siralamasi degismedi mi?  %s\n",
            ifelse(rank_ok_all, "GECTI", "KALDI")))

# ---------------------------------------------------------------------
# PART C — importation-free metrics are unchanged (no computation needed)
# ---------------------------------------------------------------------
cat("\n================ PART C: degismeyenler =================\n")
cat("P_est,tau  : ithalat terimi icermez  -> birebir AYNI\n")
cat("Esik yili  : R0 = lambda_local/gamma -> ithalattan bagimsiz, AYNI\n")
cat("Sezon boyu : aktif ay (P_est>0)      -> ithalattan bagimsiz, AYNI\n")

cat("\n================ SONUC =================\n")
cat(sprintf("Yontem 1 (dogrudan A_d), eski formulun yalnizca sabit bir\n"))
cat(sprintf("carpanla (c=%.4f) yeniden-olceklenmis halidir. Ilce siralamasi,\n", c_scale))
cat(sprintf("esik-gecis yillari ve sezon uzunlugu DEGISMEZ; yalnizca mutlak\n"))
cat(sprintf("P_ufuk olcegi makullesir (ornek SSP2-4.5 Kartal: %.3f -> %.4f).\n",
            pufuk$P_ufuk_old[pufuk$district=="Kartal" & pufuk$ssp=="ssp245"],
            pufuk$P_ufuk_new[pufuk$district=="Kartal" & pufuk$ssp=="ssp245"]))

invisible(list(c_scale = c_scale, test_A = ratio_constant, test_B = rank_ok_all))
