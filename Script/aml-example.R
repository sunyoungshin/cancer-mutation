if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("atSNP")
library(atSNP)
BiocManager::install("GenomicRanges")
library(GenomicRanges)
BiocManager::install("MotIV")
library(MotIV)
install.packages("data.table")
library(data.table)
install.packages("ggplot2")
library(ggplot2)

system("wget http://159.226.67.237/sun/oncobase/assets/data/download/ICGC+TCGA+COSMIC+ClinVar.txt.gz")
system("gunzip ICGC+TCGA+COSMIC+ClinVar.txt.gz")
system("awk '/AML/' ICGC+TCGA+COSMIC+ClinVar.txt > somatic-mutations.txt")

somatic.mutations <- fread("somatic-mutations.txt")
setnames(somatic.mutations, names(somatic.mutations), c("chr", "start", "end", "ref", "alt", "gene_region", "gene_symbol", "effect", "mutation_type", "SIFT_score", "polyphen2_HDIV_score", "gnomAD_exome_ALL", "gnomAD_genome_ALL", "avsnp147", "ExAC_ALL", "1000g2015aug_all", "esp6500siv2_all", "CG69", "cosmic70", "TCGA_Occurrence", "ICGC_Occurrence", "CLNDBN", "GWAS", "MGI_Phenotype", "HGMD_GeneBased", "Four_DB_Record"))
table(somatic.mutations$mutation_type)

snv.dt0 <- somatic.mutations[mutation_type %in% "SNV"]
snv.dt0[chr %in% "chrMT"]$chr <- "chrM"
full.snv.dt0 <- snv.dt0[!(chr %in% "chrUnknown")]
novel.snv.dt0 <- full.snv.dt0[avsnp147 %in% "-"]
novel.snv.dt0$snpid <- mapply(function(x, y) paste0(x, ":", y),  novel.snv.dt0$chr, novel.snv.dt0$start)
novel.snvids <- unique(novel.snv.dt0$snpid)
novel.snv.dt <- novel.snv.dt0[snpid %in% novel.snvids]

str(novel.snv.dt)

no.novel.snvs <- length(novel.snvids)
no.novel.snvs
full.snvids <- unique(full.snv.dt0$avsnp147)
no.full.snvs <- length(full.snvids) + length(novel.snvids) - 1
no.full.snvs
no.novel.snvs/no.full.snvs

novel.snv.input.dt <- novel.snv.dt[, c("chr", "start", "snpid", "ref", "alt")]
setnames(novel.snv.input.dt, names(novel.snv.input.dt), c("chr", "snp", "snpid", "a1", "a2"))
novel.snv.input.dt
write.table(novel.snv.input.dt, file = "novel-snv_input.txt", row.names = FALSE, quote = FALSE)

novel.snv.info <- LoadSNPData("novel-snv_input.txt", genome.lib = "BSgenome.Hsapiens.UCSC.hg19", half.window.size = 30, default.par = FALSE, mutation = FALSE)

str(novel.snv.info)
length(unique(novel.snv.info$snpid))

system('wget "http://jaspar.genereg.net/download/CORE/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.txt"')
jaspar.vert.motifs <- LoadMotifLibrary(filename="JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.txt")
atsnp.scores <- ComputeMotifScore(jaspar.vert.motifs, novel.snv.info, ncores = 30)

str(atsnp.scores$snp.tbl)
str(atsnp.scores$motif.scores)

atsnp.result <- ComputePValues(motif.lib = jaspar.vert.motifs, snp.info = novel.snv.info, motif.scores = atsnp.scores$motif.scores, ncores = 30)
atsnp.result <- as.data.table(atsnp.result)
atsnp.result[, pval_rank_bh := p.adjust(pval_rank, method = "BH"), by = motif]
atsnp.result[, pval_snp_bh := p.adjust(pval_snp, method = "BH"), by = motif]
atsnp.result[, pval_ref_bh := p.adjust(pval_ref, method = "BH"), by = motif]
str(atsnp.result)

candidate.result <- rbind(atsnp.result[pval_rank_bh<0.05 & pval_snp_bh<0.05 & pval_ref_bh>0.05], atsnp.result[pval_rank_bh<0.05 & pval_snp_bh>0.05 & pval_ref_bh<0.05])

system('grep -i "MOTIF" JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.txt > jaspar-vert_tfs.txt')
jaspar.motif.tf <- fread("jaspar-vert_tfs.txt", header = FALSE)[,2:3]
setnames(jaspar.motif.tf, 1:2, c("motif", "tf"))
jaspar.motif.tf
candidate.result.tf0 <- merge(candidate.result, jaspar.motif.tf, by="motif")

jaspar.pwm<-mapply(function(x) makePWM(t(x)), jaspar.vert.motifs)
jaspar.ic<-mapply(function(x) x@ic, jaspar.pwm, SIMPLIFY=FALSE)
jaspar.ic.median<-mapply(median, jaspar.ic)
candidate.result.tf <- candidate.result.tf0[motif %in% jaspar.motif.tf$motif[which(jaspar.ic.median >=1)]]

novel.nonexonic.snv.dt <- novel.snv.dt[!(gene_region %in% "exonic")]
setnames(candidate.result.tf, "snpbase", "alt")
setkey(novel.nonexonic.snv.dt, "snpid", "alt")
setkey(candidate.result.tf, "snpid", "alt")
candidate.snv.dt <- merge(novel.nonexonic.snv.dt, candidate.result.tf, by = c("snpid", "alt"))

candidate.snvids <- unique(candidate.snv.dt$snpid)
candidate.snv.chr.start <- unique(candidate.snv.dt[, c("chr", "start")])
candidate.snv.gr <- GRanges(seqnames=candidate.snv.chr.start$chr, IRanges(candidate.snv.chr.start$start, width=1))
candidate.snv.gr

novel.snv.chr.start <- unique(novel.snv.input.dt[, c("chr","snp")])
novel.snv.gr <- GRanges(seqnames=novel.snv.chr.start$chr, IRanges(novel.snv.chr.start$snp, width=1))
novel.snv.gr

system("wget https://www.encodeproject.org/files/ENCFF925ZDU/@@download/ENCFF925ZDU.bed.gz")
system("wget https://www.encodeproject.org/files/ENCFF041CSB/@@download/ENCFF041CSB.bed.gz")
system("wget https://www.encodeproject.org/files/ENCFF185LPE/@@download/ENCFF185LPE.bed.gz")
system("wget https://www.encodeproject.org/files/ENCFF505KHS/@@download/ENCFF505KHS.bed.gz")
system("gunzip *.gz")

da.peak.files <- c("ENCFF925ZDU.bed", "ENCFF041CSB.bed", "ENCFF185LPE.bed", "ENCFF505KHS.bed")
da.peaks.lt <- lapply(da.peak.files, function(peak.file) {peak.dt <- fread(peak.file); setnames(peak.dt, names(peak.dt), c("chrID", "peakStart", "peakStop", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"))})
da.peaks.coord.lt <- lapply(da.peaks.lt, function(peaks) peaks[,c(1:3)])
da.peaks.grange.lt <- lapply(da.peaks.coord.lt, function(peak.coord)  {setnames(peak.coord, names(peak.coord), c("chrom", "start", "end")); as(peak.coord, "GRanges")})
names(da.peaks.grange.lt) <- da.peak.files
da.peaks.grange.union <- Reduce(union, da.peaks.grange.lt)
da.peaks.grange.union

o.candidate.snv.peaks <- findOverlaps(candidate.snv.gr, da.peaks.grange.union)
o.novel.snv.peaks <- findOverlaps(novel.snv.gr, da.peaks.grange.union)
candidate.regulatory.snv.gr <- candidate.snv.gr[unique(queryHits(o.candidate.snv.peaks)),]
novel.regulatory.snv.gr <- novel.snv.gr[unique(queryHits(o.novel.snv.peaks)),]
candidate.regulatory.snv.len <- length(candidate.regulatory.snv.gr)
candidate.snv.len <- length(candidate.snv.gr)
novel.regulatory.snv.len <- length(novel.regulatory.snv.gr)
novel.snv.len <- length(novel.snv.gr)
regulatory.snv.tab <- matrix(c(candidate.regulatory.snv.len, candidate.snv.len, novel.regulatory.snv.len, novel.snv.len), 2, 2)
regulatory.snv.tab

prop.test(c(candidate.regulatory.snv.len, novel.regulatory.snv.len), c(candidate.snv.len, novel.snv.len))
phyper(candidate.regulatory.snv.len, novel.regulatory.snv.len, novel.snv.len-novel.regulatory.snv.len, candidate.snv.len, lower.tail=FALSE)

candidate.regulatory.snvids0 <- candidate.snvids[unique(queryHits(o.candidate.snv.peaks))]
candidate.regulatory.snv.dt0 <- candidate.snv.dt[snpid %in% candidate.regulatory.snvids0]

candidate.regulatory.tf <- candidate.regulatory.snv.dt0[, .N, by=tf]
candidate.regulatory.tf$tf <- toupper(candidate.regulatory.tf$tf)
de.tf <- fread("table_degenes.txt")
setnames(de.tf, 1:6, c("tf", "gene_id", "median_tumor", "median_normal", "log2FC", "adj_p"))
candidate.regulatory.de.tf0 <- merge(candidate.regulatory.tf, de.tf, by = "tf")
candidate.regulatory.de.tf <- candidate.regulatory.de.tf0[log2FC > 0 & adj_p < 0.1]
candidate.regulatory.de.tf[order(-N), ]

pdf("bargraph_tf.pdf", width = 8, height = 5, useDingbats = FALSE)
ggplot(candidate.regulatory.de.tf, aes(reorder(tf, -N), N, fill = -log10(adj_p))) +
  geom_bar(stat = "identity") +
  scale_fill_continuous(low="blue", high="red") +
  labs(x= "TF", y="SNV count")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

candidate.regulatory.snv.dt <- candidate.regulatory.snv.dt0[ tf %in% candidate.regulatory.de.tf$tf]
candidate.regulatory.snvids <- unique(candidate.regulatory.snv.dt$snpid)
candidate.regulatory.snv.dt[, list(snpid, motif, tf, ref, alt, gene_region, 
                                          gene_symbol, Four_DB_Record, pval_rank_bh, 
                                          pval_snp_bh, pval_ref_bh)]

sessionInfo()