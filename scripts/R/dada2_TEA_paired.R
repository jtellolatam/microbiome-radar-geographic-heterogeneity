#!/usr/bin/env Rscript
# M3 — reprocesamiento paired-end de Xiamen y Caltech
# Compara con el single-end usado en el manuscrito.
# Unica variable que cambia: modo de procesamiento.
suppressPackageStartupMessages({ library(dada2) })
set.seed(42)
BASE <- path.expand("~/Documents/MicrobiomeRadar_Research")
TEA  <- file.path(BASE,"01_TEA_Neurodivergent")
SILVA<- file.path(BASE,"_reference","silva_nr99_v138.2_toSpecies_trainset.fa.gz")
COH <- list(
  Xiamen  = list(dir=file.path(TEA,"01g_Xiamen"),
                 patF="_1\\.fastq\\.gz$", patR="_2\\.fastq\\.gz$"),
  Caltech = list(dir=file.path(TEA,"01h_Caltech_USA"),
                 patF="_1\\.fastq\\.gz$", patR="_2\\.fastq\\.gz$")
)
args <- commandArgs(trailingOnly=TRUE)
if (length(args)<1 || !(args[1] %in% names(COH)))
  stop("Uso: Rscript dada2_TEA_paired.R [Xiamen|Caltech]")
NOM <- args[1]; CF <- COH[[NOM]]
RAW <- file.path(CF$dir,"01_raw")
OUT <- file.path(CF$dir,"02_processed_paired")
TMP <- file.path(CF$dir,"tmp_paired")
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
dir.create(TMP, showWarnings=FALSE, recursive=TRUE)
log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(),"%H:%M:%S"), paste0(...))); flush.console()
}
log_msg("=== ",NOM," · modo PAIRED-END ===")
fnF <- sort(list.files(RAW, pattern=CF$patF, full.names=TRUE))
fnR <- sort(list.files(RAW, pattern=CF$patR, full.names=TRUE))
muestras <- sub(CF$patF,"",basename(fnF))
log_msg("Pares encontrados: ",length(fnF))
if (length(fnF)==0 || length(fnF)!=length(fnR)) stop("pares no coinciden")
filtF <- file.path(TMP,"filt",paste0(muestras,"_F.fastq.gz"))
filtR <- file.path(TMP,"filt",paste0(muestras,"_R.fastq.gz"))
names(filtF) <- muestras; names(filtR) <- muestras
log_msg("Filtrando (truncLen = 145/140)...")
out <- filterAndTrim(fnF, filtF, fnR, filtR, truncLen=c(145,140),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, verbose=FALSE)
log_msg("  retenidas: ", round(100*sum(out[,2])/sum(out[,1]),1), "%")
ok <- file.exists(filtF) & file.exists(filtR)
filtF <- filtF[ok]; filtR <- filtR[ok]; muestras <- muestras[ok]
log_msg("  muestras con datos: ",length(muestras),"/",nrow(out))
log_msg("Modelo de error...")
errF <- learnErrors(filtF, multithread=TRUE, randomize=TRUE, verbose=0)
errR <- learnErrors(filtR, multithread=TRUE, randomize=TRUE, verbose=0)
log_msg("Inferencia de ASV...")
ddF <- dada(filtF, err=errF, multithread=TRUE, verbose=0)
ddR <- dada(filtR, err=errR, multithread=TRUE, verbose=0)
log_msg("Fusionando pares...")
mer <- mergePairs(ddF, filtF, ddR, filtR, minOverlap=12, verbose=FALSE)
seqtab <- makeSequenceTable(mer)
log_msg("  ASVs antes de quimeras: ",ncol(seqtab))
log_msg("  fusion exitosa: ", round(100*sum(seqtab)/sum(sapply(ddF,function(x)sum(getUniques(x)))),1), "%")
log_msg("Quimeras...")
seqtab.nc <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=FALSE)
log_msg("  ASVs finales: ",ncol(seqtab.nc)," (",
        round(100*sum(seqtab.nc)/sum(seqtab),1),"% de lecturas)")
log_msg("Taxonomia SILVA...")
taxa <- assignTaxonomy(seqtab.nc, SILVA, multithread=TRUE, tryRC=TRUE, verbose=FALSE)
gen <- taxa[,"Genus"]; val <- !is.na(gen) & nzchar(gen)
mat <- t(rowsum(t(seqtab.nc[,val,drop=FALSE]), group=gen[val]))
log_msg("Generos: ",ncol(mat))
write.csv(mat, file.path(OUT,paste0(NOM,"_genus_counts.csv")))
saveRDS(seqtab.nc, file.path(OUT,paste0(NOM,"_seqtab.rds")))
getN <- function(x) sum(getUniques(x))
track <- data.frame(muestra=muestras, entrada=out[ok,1], filtradas=out[ok,2],
  denoisedF=sapply(ddF,getN), fusionadas=sapply(mer,function(x)sum(x$abundance)),
  sin_quimeras=rowSums(seqtab.nc), a_genero=rowSums(mat))
write.csv(track, file.path(OUT,paste0(NOM,"_track.csv")), row.names=FALSE)
log_msg("=== LISTO ===")
log_msg("Muestras: ",nrow(mat),"  Generos: ",ncol(mat))
log_msg("Profundidad mediana: ",median(rowSums(mat)))
log_msg("Retencion global: ",round(100*sum(track$a_genero)/sum(track$entrada),1),"%")
