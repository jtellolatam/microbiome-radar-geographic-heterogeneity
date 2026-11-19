#!/usr/bin/env Rscript
# Microbiome Radar - Pipeline DADA2 brazo TEA, region V4
# Las tres cohortes en modo SINGLE-END forward, truncLen=145,
# para garantizar tratamiento identico entre cohortes.
suppressPackageStartupMessages({
  library(dada2); library(Biostrings); library(vegan)
})
set.seed(42)
BASE  <- path.expand("~/Documents/MicrobiomeRadar_Research")
TEA   <- file.path(BASE, "01_TEA_Neurodivergent")
SILVA <- file.path(BASE, "_reference", "silva_nr99_v138.2_toSpecies_trainset.fa.gz")
TRUNC <- 145; MAXEE <- 2
COHORTES <- list(
  Dan = list(dir=file.path(TEA,"01f_Dan_Nanjing"), ciudad="Nanjing", pais="China",
             patron="\\.fastq\\.gz$", sufijo="\\.fastq\\.gz$"),
  Xiamen = list(dir=file.path(TEA,"01g_Xiamen"), ciudad="Xiamen", pais="China",
             patron="_1\\.fastq\\.gz$", sufijo="_1\\.fastq\\.gz$"),
  Caltech = list(dir=file.path(TEA,"01h_Caltech_USA"), ciudad="Caltech", pais="EEUU",
             patron="_1\\.fastq\\.gz$", sufijo="_1\\.fastq\\.gz$")
)
args <- commandArgs(trailingOnly=TRUE)
if (length(args)<1 || !(args[1] %in% names(COHORTES)))
  stop("Uso: Rscript dada2_TEA_v4.R [Dan|Xiamen|Caltech]")
NOM <- args[1]; CF <- COHORTES[[NOM]]
RAW <- file.path(CF$dir,"01_raw"); OUT <- file.path(CF$dir,"02_processed")
TMP <- file.path(CF$dir,"tmp_dada2")
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
dir.create(TMP, showWarnings=FALSE, recursive=TRUE)
log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(),"%H:%M:%S"), paste0(...))); flush.console()
}
log_msg("=== COHORTE: ",NOM," (",CF$ciudad,", ",CF$pais,") ===")
log_msg("Modo single-end forward, truncLen = ",TRUNC)
fn <- sort(list.files(RAW, pattern=CF$patron, full.names=TRUE))
muestras <- sub(CF$sufijo,"",basename(fn))
log_msg("Ficheros encontrados: ",length(fn))
if (length(fn)==0) stop("Sin ficheros que coincidan")
pdf(file.path(OUT,paste0(NOM,"_calidad.pdf")), width=10, height=5)
print(plotQualityProfile(fn[1:min(4,length(fn))])); dev.off()
log_msg("Perfil de calidad guardado")
filt <- file.path(TMP,"filtered",paste0(muestras,"_F.fastq.gz"))
names(filt) <- muestras
log_msg("Filtrando...")
out <- filterAndTrim(fn, filt, truncLen=TRUNC, maxN=0, maxEE=MAXEE, truncQ=2,
                     rm.phix=TRUE, compress=TRUE, multithread=TRUE, verbose=FALSE)
log_msg("  retenidas: ", round(100*sum(out[,2])/sum(out[,1]),1), "%")
existe <- file.exists(filt); filt <- filt[existe]; muestras <- muestras[existe]
log_msg("  muestras con datos: ",length(muestras),"/",nrow(out))
log_msg("Aprendiendo modelo de error...")
err <- learnErrors(filt, multithread=TRUE, randomize=TRUE, verbose=0)
pdf(file.path(OUT,paste0(NOM,"_errores.pdf")), width=8, height=6)
print(plotErrors(err, nominalQ=TRUE)); dev.off()
log_msg("Inferencia de ASV...")
dd <- dada(filt, err=err, multithread=TRUE, verbose=0)
seqtab <- makeSequenceTable(dd)
log_msg("  ASVs antes de quimeras: ",ncol(seqtab))
log_msg("Eliminando quimeras...")
seqtab.nc <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=FALSE)
log_msg("  ASVs finales: ",ncol(seqtab.nc),"  (",
        round(100*sum(seqtab.nc)/sum(seqtab),1),"% de lecturas)")
log_msg("Asignando taxonomia (SILVA)...")
taxa <- assignTaxonomy(seqtab.nc, SILVA, multithread=TRUE, tryRC=TRUE, verbose=FALSE)
gen <- taxa[,"Genus"]; val <- !is.na(gen) & nzchar(gen)
log_msg("ASVs con genero asignado: ",sum(val)," / ",length(gen))
mat_gen <- t(rowsum(t(seqtab.nc[,val,drop=FALSE]), group=gen[val]))
log_msg("Generos distintos: ",ncol(mat_gen))
getN <- function(x) sum(getUniques(x))
track <- data.frame(muestra=muestras, entrada=out[existe,1], filtradas=out[existe,2],
  denoised=sapply(dd,getN), sin_quimeras=rowSums(seqtab.nc), a_genero=rowSums(mat_gen))
write.csv(mat_gen, file.path(OUT,paste0(NOM,"_genus_counts.csv")))
write.csv(track, file.path(OUT,paste0(NOM,"_track.csv")), row.names=FALSE)
saveRDS(seqtab.nc, file.path(OUT,paste0(NOM,"_seqtab.rds")))
saveRDS(taxa, file.path(OUT,paste0(NOM,"_taxa.rds")))
writeLines(capture.output(sessionInfo()), file.path(OUT,paste0(NOM,"_sessionInfo.txt")))
log_msg("=== LISTO ===")
log_msg("Muestras: ",nrow(mat_gen),"  Generos: ",ncol(mat_gen))
log_msg("Profundidad mediana: ",median(rowSums(mat_gen)))
log_msg("Retencion global: ",round(100*sum(track$a_genero)/sum(track$entrada),1),"%")
log_msg("Salidas en: ",OUT)
