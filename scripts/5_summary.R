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
GC12 = res$f_GC_1 + res$f_GC_2 
modL = loess(cov~GC12)

GCcorrected_expression = residuals(modL)

res = cbind(res, GCcorrected_expression)
colnames(res)[1] = "contigs"
colnames(res)[2] = "coverage"

# statistics for all contigs
write.table(res, col.names=T, row.names=F, quote=F, file="output_contigs_v1.txt", sep="\t", na = "nan")

# statistics for the whole species
nLoci = nrow(res)

L_avg = mean(res$n_codons*3)
L_std = sd(res$n_codons*3)

GCtot_avg = mean(res$f_GC_tot)
GCtot_std = sd(res$f_GC_tot)

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

tmp = round(c(nLoci, L_avg, L_std, GCtot_avg, GCtot_std, GC12_avg, GC12_std, GC3_avg, GC3_std, r_expression_ENc, r_expression_ENcP, r_nReads_ENc, r_nReads_ENcP), 5)
res2 = c(dataset, tmp)
names(res2) = c("dataset", "nLoci", "length_avg", "length_std", "GCtot_avg", "GCtot_std", "GC12_avg", "GC12_std", "GC3_avg", "GC3_std", "r_expression_ENc", "r_expression_ENcP", "r_nReads_ENc", "r_nReads_ENcP")

write.table(t(res2), file = "output_summarized_v1.txt", quote = F, sep = "\t", row.names = F, col.names = T)

