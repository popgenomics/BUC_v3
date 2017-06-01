#!/usr/bin/Rscript
for(i in commandArgs()){
	tmp = strsplit(i, "=")
	if(tmp[[1]][1] == "species"){dataset = tmp[[1]][2]}
}

x = read.table("output_ENCprime_JoNo.txt", h=T)
y = read.table("coverage.txt", h=F)

#traits = read.csv("../../scripts/life_history_traits.csv", h=T)

res = merge(y, x, by.x=1, by.y=1)

cov = res[,2]
GCutr = res$f_GC_tot 
modL = loess(cov~GCutr)

GCcorrected_expression = residuals(modL)

res = cbind(res, GCcorrected_expression)
colnames(res)[1] = "contigs"
colnames(res)[2] = "coverage"

# statistics for all contigs
write.table(res, col.names=T, row.names=F, quote=F, file="output_contigs_v1.txt", sep="\t", na = "nan")

# statistics for the whole species
nLoci = nrow(res)

Lcds_avg = mean(res$n_codons*3)
Lcds_std = sd(res$n_codons*3)

Lutr_avg = mean(res$L_utr)
Lutr_std = sd(res$L_utr)

GCutr_avg = mean(res$f_GC_utr)
GCutr_std = sd(res$f_GC_utr)

GCcds_avg = mean(res$f_GC_tot)
GCcds_std = sd(res$f_GC_tot)

GC12_avg = mean((res$f_GC_1 + res$f_GC_2)/2)
GC12_std = sd((res$f_GC_1 + res$f_GC_2)/2)

GC3_avg = mean(res$f_GC_3)
GC3_std = sd(res$f_GC_3)

expression_ENc = cor.test(res$GCcorrected_expression, res$Nc)
r_expression_ENc = expression_ENc$estimate
pval_expression_ENc = expression_ENc$p.value

expression_ENcP = cor.test(res$GCcorrected_expression, res$Ncp)
r_expression_ENcP = expression_ENcP$estimate
pval_expression_ENcP = expression_ENcP$p.value

nReads_ENc = cor.test(res$coverage, res$Nc)
r_nReads_ENc = nReads_ENc$estimate
pval_nReads_ENc = nReads_ENc$p.value

nReads_ENcP = cor.test(res$coverage, res$Ncp)
r_nReads_ENcP = nReads_ENcP$estimate
pval_nReads_ENcP = nReads_ENcP$p.value

GC3_GCutr = cor.test(res$f_GC_3, res$f_GC_utr)
r_GC3_GCutr = GC3_GCutr$estimate
pval_GC3_GCutr = GC3_GCutr$p.value

GC12_GCutr = cor.test(res$f_GC_1 + res$f_GC_2, res$f_GC_utr)
r_GC12_GCutr = GC12_GCutr$estimate
pval_GC12_GCutr = GC12_GCutr$p.value

tmp = round(c(nLoci, Lcds_avg, Lcds_std, GCcds_avg, GCcds_std, GC12_avg, GC12_std, GC3_avg, GC3_std, Lutr_avg, Lutr_std, GCutr_avg, GCutr_std, r_expression_ENc, r_expression_ENcP, r_nReads_ENc, r_nReads_ENcP, r_GC3_GCutr, r_GC12_GCutr), 5)
res2 = c(dataset, tmp)
names(res2) = c("dataset", "nLoci", "length_cds_avg", "length_cds_std", "GCcds_avg", "GCcds_std", "GC12_avg", "GC12_std", "GC3_avg", "GC3_std", "length_utr_avg", "length_utr_std", "GCutr_avg", "GCutr_std", "r_expression_ENc", "r_expression_ENcP", "r_nReads_ENc", "r_nReads_ENcP", "r_GC3_GCutr", "r_GC12_GCutr")

write.table(t(res2), file = "output_summarized_v1.txt", quote = F, sep = "\t", row.names = F, col.names = T)

