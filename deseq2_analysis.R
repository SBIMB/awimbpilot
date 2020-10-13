#####This was the code used for the differential abundance  analyses section (using DESeq2) of the paper by: 
##Oduaran, OH., Tamburini, FB., Sahibdeen, V., Brewster, R., Gómez-Olivé, FX., Kahn, K., Norris, SA., Tollman, SM., Twine, R., Wade, AN., Wagner, RG., Lombard, Z., Bhatt, AS., Hazelhurst, S.
##Gut Microbiome Profiling of a Rural and Urban South African Cohort Reveals Biomarkers of a Population in Lifestyle Transition

##Adapted from: https://bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-mixture-models.html
##http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
##and https://www.huber.embl.de/users/klaus/Teaching/DESeq2Predoc2014.html 


dds_all<-phyloseq_to_deseq2(filt_agt_swt_sub, ~ batch + bmi_info + site)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds_all), 1, gm_mean)
dds_all = estimateSizeFactors(dds_all, geoMeans = geoMeans)
dds_all = estimateDispersions(dds_all, fitType="local")
as_all<-nbinomWaldTest(dds_all)
res_all_test<-results(as_all, pAdjustMethod = "BH")
table(res_all_test$padj < 0.1)
res_all_test<-res_all_test[ !is.na(res_all_test$padj), ]
res_all_test<-res_all_test[ !is.na(res_all_test$pvalue), ]
res_all_test <- res_all_test[, -which(names(res_all_test) == "padj")]
FDR.res_all <- fdrtool(res_all_test$stat, statistic= "normal", plot = T)
FDR.res_all$param[1, "sd"]
res_all_test[,"padj"]  <- p.adjust(FDR.res_all$pval, method = "BH")
res_all_p<-as.data.frame(res_all_test)
sigtab<-res_all_p
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(filt_agt_swt_sub)[rownames(sigtab),], "matrix"))
alpha=0.1
sigtab = sigtab[(sigtab$padj < alpha), ]
sigtabgen<-sigtab
theme_set(theme_classic())
v_plot_agt_swt_all<-ggplot(sigtabgen) + geom_point(aes(y=-log10(padj), x=log2FoldChange, color=Phylum)) + geom_hline(yintercept=-log10(alpha), col="brown", linetype ="dashed") + geom_vline(xintercept= 0, col="brown", linetype ="dashed") + labs(x = "Log2 Fold Ratio", title = "Site - All" ) + scale_color_manual(values=c(phy_col1, phy_col2, phy_col3), na.value= "#DEAA03") + theme(legend.text=element_text(face = "italic"))
sig_gen_as<- sigtab[(sigtab$padj < alpha), ]
write.csv(sig_gen_as, "sig_gen_as_A.csv")

dds_a_lnob = phyloseq_to_deseq2(filt_agt_lnob, ~ batch + bmi_info) #Agt lean vs obese comparisons
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds_a_lnob), 1, gm_mean)
dds_a_lnob = estimateSizeFactors(dds_a_lnob, geoMeans = geoMeans)
dds_a_lnob = estimateDispersions(dds_a_lnob, fitType="local")
as_all<-nbinomWaldTest(dds_a_lnob)
res_a_lnob_test<-results(as_all, pAdjustMethod = "BH")
table(res_a_lnob_test$padj < 0.1)
res_a_lnob_test<-res_a_lnob_test[ !is.na(res_a_lnob_test$padj), ]
res_a_lnob_test<-res_a_lnob_test[ !is.na(res_a_lnob_test$pvalue), ]
res_a_lnob_test <- res_a_lnob_test[, -which(names(res_a_lnob_test) == "padj")]
FDR.res_a_lnob <- fdrtool(res_a_lnob_test$stat, statistic= "normal", plot = T)
FDR.res_a_lnob$param[1, "sd"]
res_a_lnob_test[,"padj"]  <- p.adjust(FDR.res_a_lnob$pval, method = "BH")
res_a_lnob_p<-as.data.frame(res_a_lnob_test)
sigtab<-res_a_lnob_p
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(filt_agt_lnob)[rownames(sigtab),], "matrix"))
alpha=0.1
sigtab = sigtab[(sigtab$padj < alpha), ]
sigtabgen<-sigtab
theme_set(theme_classic())
v_plot_agt_lnob<-ggplot(sigtabgen) + geom_point(aes(y=-log10(padj), x=log2FoldChange, color=Phylum)) + geom_hline(yintercept=-log10(alpha), col="brown", linetype ="dashed") + geom_vline(xintercept= 0, col="brown", linetype ="dashed") + labs(x = "Log2 Fold Ratio", title = "BMI - BBR") + scale_color_manual(values= c(phy_col1, phy_col2, phy_col3), na.value= "#DEAA03") + theme(legend.text=element_text(face = "italic")) ##BBR (Obese vs Lean)")
sig_gen_alnob<- sigtab[(sigtab$padj < alpha), ]
write.csv(sig_gen_alnob, "sig_gen_agt_lnob_E.csv")

dds_s_lnob = phyloseq_to_deseq2(filt_swt_lnob, ~ bmi_info)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds_s_lnob), 1, gm_mean)
dds_s_lnob = estimateSizeFactors(dds_s_lnob, geoMeans = geoMeans)
dds_s_lnob = estimateDispersions(dds_s_lnob, fitType="local")
as_all<-nbinomWaldTest(dds_s_lnob)
res_s_lnob_test<-results(as_all, pAdjustMethod = "BH")
table(res_s_lnob_test$padj < 0.1) 
res_s_lnob_test<-res_s_lnob_test[ !is.na(res_s_lnob_test$padj), ]
res_s_lnob_test<-res_s_lnob_test[ !is.na(res_s_lnob_test$pvalue), ]
res_s_lnob_test <- res_s_lnob_test[, -which(names(res_s_lnob_test) == "padj")]
FDR.res_s_lnob <- fdrtool(res_s_lnob_test$stat, statistic= "normal", plot = T)
FDR.res_s_lnob$param[1, "sd"]

res_s_lnob_test[,"padj"]  <- p.adjust(FDR.res_s_lnob$pval, method = "BH")
res_s_lnob_p<-as.data.frame(res_s_lnob_test)
sigtab<-res_s_lnob_p
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(filt_swt_lnob)[rownames(sigtab),], "matrix"))
alpha=0.1
sigtab = sigtab[(sigtab$padj < alpha), ]
sigtabgen<-sigtab
theme_set(theme_classic())
v_plot_swt_lnob<-ggplot(sigtabgen) + geom_point(aes(y=-log10(padj), x=log2FoldChange, color=Phylum)) + geom_hline(yintercept=-log10(alpha), col="brown", linetype ="dashed") + geom_vline(xintercept= 0, col="brown", linetype ="dashed") + labs(x = "Log2 Fold Ratio", title = "BMI - SWT") + scale_color_manual(values= c(phy_col1, phy_col2, phy_col3), na.value= "#DEAA03") + theme(legend.text=element_text(face = "italic"))## SWT (Obese vs Lean)
sig_gen_slnob<- sigtab[(sigtab$padj < alpha), ]
write.csv(sig_gen_slnob, "sig_gen_swt_lnob_F.csv")

dds_as_lnob = phyloseq_to_deseq2(filt_agt_swt_lnob, ~ batch + bmi_info) #Lean vs obese comparisons of pilot cohort excluding overweight samples whilst controlling for batch effects
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds_as_lnob), 1, gm_mean)
dds_as_lnob = estimateSizeFactors(dds_as_lnob, geoMeans = geoMeans)
dds_as_lnob = estimateDispersions(dds_as_lnob, fitType="local")
as_all<-nbinomWaldTest(dds_as_lnob)
res_as_lnob_test<-results(as_all, pAdjustMethod = "BH")
table(res_as_lnob_test$padj < 0.1)
res_as_lnob_test<-res_as_lnob_test[ !is.na(res_as_lnob_test$padj), ]
res_as_lnob_test<-res_as_lnob_test[ !is.na(res_as_lnob_test$pvalue), ]
res_as_lnob_test <- res_as_lnob_test[, -which(names(res_as_lnob_test) == "padj")]
FDR.res_as_lnob <- fdrtool(res_as_lnob_test$stat, statistic= "normal", plot = T)
FDR.res_as_lnob$param[1, "sd"]
res_as_lnob_test[,"padj"]  <- p.adjust(FDR.res_as_lnob$pval, method = "BH")
res_as_lnob_p<-as.data.frame(res_as_lnob_test)
sigtab<-res_as_lnob_p
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(filt_agt_swt_lnob)[rownames(sigtab),], "matrix"))
alpha=0.1
sigtab = sigtab[(sigtab$padj < alpha), ]
sigtabgen<-sigtab
theme_set(theme_classic())
v_plot_agt_swt_lnob<-ggplot(sigtabgen) + geom_point(aes(y=-log10(padj), x=log2FoldChange, color=Phylum)) + geom_hline(yintercept=-log10(alpha), col="brown", linetype ="dashed") + geom_vline(xintercept= 0, col="brown", linetype ="dashed") + labs(x = "Log2 Fold Ratio",  title = "BMI - All") + scale_color_manual(values= c(phy_col1, phy_col2, phy_col3), na.value= "#DEAA03") + theme(legend.text=element_text(face = "italic")) ##ALL – Obese vs Lean
sig_gen_aslnob<- sigtab[(sigtab$padj < alpha), ]
write.csv(sig_gen_aslnob, "sig_gen_agt_swt_lnob_D.csv")

dds_as_ob = phyloseq_to_deseq2(filt_agt_swt_ob, ~ site) #Comparing compositional differencc
s between sites in only obese samples
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds_as_ob), 1, gm_mean)
dds_as_ob = estimateSizeFactors(dds_as_ob, geoMeans = geoMeans)
dds_as_ob = estimateDispersions(dds_as_ob, fitType="local")
as_ob<-nbinomWaldTest(dds_as_ob)
res_as_ob_test<-results(as_ob, pAdjustMethod = "BH")
table(res_as_ob_test$padj < 0.1)
res_as_ob_test<-res_as_ob_test[ !is.na(res_as_ob_test$padj), ]
res_as_ob_test<-res_as_ob_test[ !is.na(res_as_ob_test$pvalue), ]
res_as_ob_test <- res_as_ob_test[, -which(names(res_as_ob_test) == "padj")]
FDR.res_as_ob <- fdrtool(res_as_ob_test$stat, statistic= "normal", plot = T)
FDR.res_as_ob$param[1, "sd"]
res_as_ob_test[,"padj"]  <- p.adjust(FDR.res_as_ob$pval, method = "BH")
res_as_ob_p<-as.data.frame(res_as_ob_test)
sigtab<-res_as_ob_p
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(filt_agt_swt_ob)[rownames(sigtab),], "matrix"))
alpha=0.1
sigtab = sigtab[(sigtab$padj < alpha), ]
sigtabgen<-sigtab
theme_set(theme_classic())
v_plot_agt_swt_ob<-ggplot(sigtabgen) + geom_point(aes(y=-log10(padj), x=log2FoldChange, color=Phylum)) + geom_hline(yintercept=-log10(alpha), col="brown", linetype ="dashed") + geom_vline(xintercept= 0, col="brown", linetype ="dashed") + labs(x = "Log2 Fold Ratio", title = "Site - Obese") + scale_color_manual(values= c(phy_col1, phy_col2, phy_col3), na.value= "#DEAA03")  + theme(legend.text=element_text(face = "italic"))## Obese (SWT vs BBR)
sig_gen_as_ob<- sigtab[(sigtab$padj < alpha), ]
write.csv(sig_gen_as_ob, "sig_gen_agt_swt_ob_C.csv")

dds_as_ln = phyloseq_to_deseq2(filt_agt_swt_ln, ~ batch + site) #Comparing compositional differences between sites in only lean samples whilst controlling for potential batch effects
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(dds_as_ln), 1, gm_mean)
dds_as_ln = estimateSizeFactors(dds_as_ln, geoMeans = geoMeans)
dds_as_ln = estimateDispersions(dds_as_ln, fitType="local")
as_ln<-nbinomWaldTest(dds_as_ln)
res_as_ln_test<-results(as_ln, pAdjustMethod = "BH")
table(res_as_ln_test$padj < 0.1)
res_as_ln_test<-res_as_ln_test[ !is.na(res_as_ln_test$padj), ]
res_as_ln_test<-res_as_ln_test[ !is.na(res_as_ln_test$pvalue), ]
res_as_ln_test <- res_as_ln_test[, -which(names(res_as_ln_test) == "padj")]
FDR.res_as_ln <- fdrtool(res_as_ln_test$stat, statistic= "normal", plot = T)
FDR.res_as_ln$param[1, "sd"]

res_as_ln_test[,"padj"]  <- p.adjust(FDR.res_as_ln$pval, method = "BH")
res_as_ln_p<-as.data.frame(res_as_ln_test)
sigtab<-res_as_ln_p
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(filt_agt_swt_ln)[rownames(sigtab),], "matrix"))
alpha=0.1
sigtab = sigtab[(sigtab$padj < alpha), ]
sigtabgen<-sigtab
theme_set(theme_classic())
v_plot_agt_swt_ln<-ggplot(sigtabgen) + geom_point(aes(y=-log10(padj), x=log2FoldChange, color=Phylum)) + geom_hline(yintercept=-log10(alpha ), col="brown", linetype ="dashed") + geom_vline(xintercept= 0, col="brown", linetype ="dashed") + labs(x = "Log2 Fold Ratio", title = "Site - Lean") + scale_color_manual(values= c(phy_col1, phy_col2, phy_col3), na.value= "#DEAA03") + theme(legend.text=element_text(face = "italic")) ##Lean (SWT vs BBR)
sig_gen_as_ln<- sigtab[(sigtab$padj < alpha), ]
write.csv(sig_gen_as_ln, "sig_gen_agt_swt_ln_B.csv")

png("Fig5_deseq2_volcano_plots2.png", units="mm", width=250, height=200, res=300)
plot_grid(v_plot_agt_swt_all, v_plot_agt_swt_ln, v_plot_agt_swt_ob, v_plot_agt_swt_lnob,  v_plot_agt_lnob, v_plot_swt_lnob, labels = c('A', 'B', 'C', 'D', 'E', 'F'), align = 'vh')
dev.off()                            

fig5_grid<-plot_grid(v_plot_agt_swt_all + theme(legend.position = "none"), v_plot_agt_swt_lnob + theme(legend.position = "none"),
v_plot_agt_lnob + theme(legend.position = "none"),
v_plot_swt_lnob  + theme(legend.position = "none"),
v_plot_agt_swt_ln + theme(legend.position = "none"),
v_plot_agt_swt_ob  + theme(legend.position = "none"),
labels = c('A', 'B', 'C', 'D', 'E', 'F'), align = 'vh')

#created mock data set with phyla represented in analysis and plot data 
#extract the legend from mock plot containing all phyla - see below.
# create some space to the left of the legend

legend<-get_legend(test + theme(legend.box.margin = margin(0, 0, 0, 12)))
png("fig5_ed_volcano.png", units="mm", width=250, height=200, res=300) 
plot_grid(fig5_grid, legend, rel_widths = c(3, .9))
dev.off()

##Mock Plot
test<-ggplot(sigtabgen_test) + geom_point(aes(y=-log10(padj), x=log2FoldChange, color=Phylum)) + geom_hline(yintercept=-log10(alpha ), col="brown", linetype ="dashed") + labs(x = "Log2 Fold Ratio", title = "BMI: Lean ( SWT vs BBR)") + scale_color_manual(values= c(phy_col1, phy_col2, phy_col3), na.value= "#DEAA03") + theme(legend.text=element_text(face = "italic"))



