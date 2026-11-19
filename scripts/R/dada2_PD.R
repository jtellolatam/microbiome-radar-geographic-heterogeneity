#!/usr/bin/env Rscript
# =============================================================================
# Microbiome Radar - Pipeline DADA2 unificado
# Brazo Alzheimer, bloque A (China): Li 2019, Ling 2021, Zhuang 2019
#
# Region 16S: V3-V4 en las tres.
# Los primers vienen recortados en Ling y presentes con barcode variable
# en Li y Zhuang, por eso el recorte es POR SECUENCIA, nunca por posicion fija.
#
# Uso:
#   Rscript dada2_AD_china.R Li
#   Rscript dada2_AD_china.R Ling
#   Rscript dada2_AD_china.R Zhuang
#
# Ejecutar una cohorte a la vez. Cada una tarda 1-3 h en un M4.
# =============================================================================

suppressPackageStartupMessages({
  library(dada2)
  library(Biostrings)
  library(ShortRead)
})

set.seed(42)   # invariante MR-30. Unico punto estocastico del pipeline.

# =============================================================================
# CONFIG
# =============================================================================
BASE <- path.expand("~/Documents/MicrobiomeRadar_Research")
AD   <- file.path(BASE, "02_Parkinson_Neurodegenerative")
SILVA <- file.path(BASE, "_reference", "silva_nr99_v138.2_toSpecies_trainset.fa.gz")

# Primers V3-V4. Li y Zhuang usan 338F/806R; Ling ya viene recortada.
FWD_338 <- "ACTCCTACGGGAGGCAGCAG"   # 338F (Li, Zhuang)
REV_806 <- "GGACTACHVGGGTWTCTAAT"   # 806R
FWD_341 <- "CCTACGGGNGGCWGCAG"      # 341F (Ling)
REV_785 <- "ACTACHVGGGTATCTAATCC"   # 785R

COHORTES <- list(
  Finlandia = list(
    dir       = file.path(AD, "02e_Finlandia"),
    ciudad    = "Finlandia",
    primers   = TRUE,
    fwd       = FWD_341, rev = REV_785,
    trunc     = c(280, 250),
    solo_heces = FALSE
  ),
  Rumania = list(
    dir       = file.path(AD, "02f_Rumania"),
    ciudad    = "Rumania",
    primers   = TRUE,
    patF      = "_R1_001\\.fastq\\.gz$",
    patR      = "_R2_001\\.fastq\\.gz$",
    fwd       = FWD_341, rev = REV_785,
    trunc     = c(275, 230),
    solo_heces = FALSE
  )
)

# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1 || !(args[1] %in% names(COHORTES))) {
  stop("Uso: Rscript dada2_PD.R [Finlandia|Rumania]")
}
NOM <- args[1]
CF  <- COHORTES[[NOM]]

RAW <- file.path(CF$dir, "01_raw")
OUT <- file.path(CF$dir, "02_processed")
TMP <- file.path(CF$dir, "tmp_dada2")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(TMP, showWarnings = FALSE, recursive = TRUE)

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
  flush.console()
}

log_msg("=== COHORTE: ", NOM, " (", CF$ciudad, ") ===")

# =============================================================================
# 1. Localizar ficheros
# =============================================================================
pat_F <- if (!is.null(CF$patF)) CF$patF else "_1\\.fastq\\.gz$"
pat_R <- if (!is.null(CF$patR)) CF$patR else "_2\\.fastq\\.gz$"
fnF <- sort(list.files(RAW, pattern = pat_F, full.names = TRUE))
fnR <- sort(list.files(RAW, pattern = pat_R, full.names = TRUE))
stopifnot(length(fnF) == length(fnR), length(fnF) > 0)

muestras <- sub(pat_F, "", basename(fnF))
log_msg("Ficheros encontrados: ", length(fnF), " pares")

# --- Li: quedarse solo con heces (prefijo F_ en library_name) ---
if (CF$solo_heces) {
  mapa <- file.path(CF$dir, "00_screening", "mapa_heces.tsv")
  if (!file.exists(mapa)) stop("Falta ", mapa, " para filtrar heces")
  heces <- read.delim(mapa, stringsAsFactors = FALSE)
  keep  <- muestras %in% heces$run_accession
  log_msg("Filtro heces: ", length(fnF), " -> ", sum(keep), " pares")
  fnF <- fnF[keep]; fnR <- fnR[keep]; muestras <- muestras[keep]
}

# =============================================================================
# 2. Recorte de primers POR SECUENCIA (no por posicion)
# =============================================================================
recortar_primers <- function(fF, fR, oF, oR, pF, pR, mm = 2, chunk = 1e5) {
  sF <- FastqStreamer(fF, chunk); sR <- FastqStreamer(fR, chunk)
  on.exit({ close(sF); close(sR) })
  if (file.exists(oF)) file.remove(oF)
  if (file.exists(oR)) file.remove(oR)
  tot <- 0L; kept <- 0L
  repeat {
    rF <- yield(sF); rR <- yield(sR)
    if (length(rF) == 0) break
    tot <- tot + length(rF)
    mF <- vmatchPattern(pF, sread(rF), max.mismatch = mm, fixed = FALSE)
    mR <- vmatchPattern(pR, sread(rR), max.mismatch = mm, fixed = FALSE)
    hitF <- elementNROWS(mF) > 0
    hitR <- elementNROWS(mR) > 0
    ok <- hitF & hitR
    if (!any(ok)) next
    endF <- rep(NA_integer_, length(rF)); endR <- endF
    endF[hitF] <- sapply(end(mF[hitF]), `[`, 1)
    endR[hitR] <- sapply(end(mR[hitR]), `[`, 1)
    rF2 <- narrow(rF[ok], start = endF[ok] + 1L)
    rR2 <- narrow(rR[ok], start = endR[ok] + 1L)
    largo <- width(rF2) >= 50 & width(rR2) >= 50
    rF2 <- rF2[largo]; rR2 <- rR2[largo]
    if (length(rF2) == 0) next
    writeFastq(rF2, oF, mode = "a", compress = TRUE)
    writeFastq(rR2, oR, mode = "a", compress = TRUE)
    kept <- kept + length(rF2)
  }
  c(entrada = tot, salida = kept)
}

if (CF$primers) {
  log_msg("Recortando primers: ", CF$fwd, " / ", CF$rev)
  cutF <- file.path(TMP, paste0(muestras, "_cut_1.fastq.gz"))
  cutR <- file.path(TMP, paste0(muestras, "_cut_2.fastq.gz"))
  stats <- matrix(0L, nrow = length(muestras), ncol = 2,
                  dimnames = list(muestras, c("entrada", "salida")))
  for (i in seq_along(muestras)) {
    stats[i, ] <- recortar_primers(fnF[i], fnR[i], cutF[i], cutR[i], CF$fwd, CF$rev)
    if (i %% 20 == 0) log_msg("  ", i, "/", length(muestras), " muestras")
  }
  pct <- 100 * sum(stats[, "salida"]) / sum(stats[, "entrada"])
  log_msg("  lecturas retenidas: ", round(pct, 1), "%")
  if (pct < 50) warning("Se perdio >50% en el recorte. Revisar primers.")
  ok <- stats[, "salida"] > 0
  log_msg("  muestras con datos: ", sum(ok), "/", length(ok))
  cutF <- cutF[ok]; cutR <- cutR[ok]; muestras <- muestras[ok]
  write.csv(stats, file.path(OUT, paste0(NOM, "_recorte_primers.csv")))
} else {
  log_msg("Primers ya recortados en origen: se omite el recorte")
  cutF <- fnF; cutR <- fnR
}

# =============================================================================
# 3. Perfiles de calidad (para auditar los parametros de truncado)
# =============================================================================
pdf(file.path(OUT, paste0(NOM, "_calidad.pdf")), width = 10, height = 5)
print(plotQualityProfile(cutF[1:min(4, length(cutF))]))
print(plotQualityProfile(cutR[1:min(4, length(cutR))]))
dev.off()
log_msg("Perfiles de calidad -> ", NOM, "_calidad.pdf")

# =============================================================================
# 4. Filtrado y truncado
# =============================================================================
filtF <- file.path(TMP, "filtered", paste0(muestras, "_F.fastq.gz"))
filtR <- file.path(TMP, "filtered", paste0(muestras, "_R.fastq.gz"))
names(filtF) <- muestras; names(filtR) <- muestras

log_msg("Filtrando (truncLen = ", CF$trunc[1], "/", CF$trunc[2], ")...")
out <- filterAndTrim(cutF, filtF, cutR, filtR,
                     truncLen = CF$trunc,
                     maxN = 0, maxEE = c(2, 5), truncQ = 2,
                     rm.phix = TRUE, compress = TRUE,
                     multithread = TRUE, verbose = FALSE)
log_msg("  retenidas: ", round(100 * sum(out[,2]) / sum(out[,1]), 1), "%")

# quitar muestras sin lecturas tras el filtro
existe <- file.exists(filtF) & file.exists(filtR)
filtF <- filtF[existe]; filtR <- filtR[existe]
muestras <- muestras[existe]
log_msg("  muestras con datos: ", length(muestras))

# =============================================================================
# 5. Modelo de error e inferencia de ASV
# =============================================================================
log_msg("Aprendiendo modelo de error (forward)...")
errF <- learnErrors(filtF, multithread = TRUE, randomize = TRUE, verbose = 0)
log_msg("Aprendiendo modelo de error (reverse)...")
errR <- learnErrors(filtR, multithread = TRUE, randomize = TRUE, verbose = 0)

pdf(file.path(OUT, paste0(NOM, "_errores.pdf")), width = 8, height = 6)
print(plotErrors(errF, nominalQ = TRUE))
print(plotErrors(errR, nominalQ = TRUE))
dev.off()

log_msg("Inferencia de ASV...")
ddF <- dada(filtF, err = errF, multithread = TRUE, verbose = 0)
ddR <- dada(filtR, err = errR, multithread = TRUE, verbose = 0)

log_msg("Fusionando pares...")
mrg <- mergePairs(ddF, filtF, ddR, filtR, minOverlap = 12, verbose = FALSE)

seqtab <- makeSequenceTable(mrg)
log_msg("  ASVs antes de quimeras: ", ncol(seqtab))

# =============================================================================
# 6. Quimeras
# =============================================================================
log_msg("Eliminando quimeras...")
seqtab.nc <- removeBimeraDenovo(seqtab, method = "consensus",
                                multithread = TRUE, verbose = FALSE)
log_msg("  ASVs finales: ", ncol(seqtab.nc),
        "  (", round(100 * sum(seqtab.nc) / sum(seqtab), 1), "% de lecturas)")

# =============================================================================
# 7. Taxonomia con SILVA
# =============================================================================
log_msg("Asignando taxonomia (SILVA)... esto tarda")
taxa <- assignTaxonomy(seqtab.nc, SILVA, multithread = TRUE,
                       tryRC = TRUE, verbose = FALSE)

# =============================================================================
# 8. Colapso a genero
# =============================================================================
gen <- taxa[, "Genus"]
val <- !is.na(gen) & nzchar(gen)
log_msg("ASVs con genero asignado: ", sum(val), " / ", length(gen))

mat_gen <- t(rowsum(t(seqtab.nc[, val, drop = FALSE]), group = gen[val]))
log_msg("Generos distintos: ", ncol(mat_gen))

# =============================================================================
# 9. Tabla de seguimiento (para Methods)
# =============================================================================
getN <- function(x) sum(getUniques(x))
track <- data.frame(
  muestra   = muestras,
  filtradas = out[existe, 2],
  denoisedF = sapply(ddF, getN),
  denoisedR = sapply(ddR, getN),
  fusionadas = sapply(mrg, getN),
  sin_quimeras = rowSums(seqtab.nc),
  a_genero  = rowSums(mat_gen)
)

# =============================================================================
# 10. Salidas
# =============================================================================
write.csv(mat_gen,  file.path(OUT, paste0(NOM, "_genus_counts.csv")))
write.csv(track,    file.path(OUT, paste0(NOM, "_track.csv")), row.names = FALSE)
saveRDS(seqtab.nc,  file.path(OUT, paste0(NOM, "_seqtab.rds")))
saveRDS(taxa,       file.path(OUT, paste0(NOM, "_taxa.rds")))
writeLines(capture.output(sessionInfo()),
           file.path(OUT, paste0(NOM, "_sessionInfo.txt")))

log_msg("=== LISTO ===")
log_msg("Muestras: ", nrow(mat_gen), "  Generos: ", ncol(mat_gen))
log_msg("Profundidad mediana: ", median(rowSums(mat_gen)))
log_msg("Salidas en: ", OUT)
log_msg("Revisar ", NOM, "_track.csv antes de modelar: si se pierde ",
        ">50% en algun paso, ajustar parametros.")









q()
