###This was the code used for the marker taxa analyses section of the paper by:
##Oduaran, OH., Tamburini, FB., Sahibdeen, V., Brewster, R., Gómez-Olivé, FX., Kahn, K., Norris, SA., Tollman, SM., Twine, R., Wade, AN., Wagner, RG., Lombard, Z., Bhatt, AS., Hazelhurst, S.
##Gut Microbiome Profiling of a Rural and Urban South African Cohort Reveals Biomarkers of a Population in Lifestyle Transition 

##Adapted from:
##https://f1000research.com/articles/5-1492/v2#f19
##https://github.com/LangilleLab/microbiome_helper/wiki/Random-Forest-in-R-with-Large-Sample-Sizes
##https://www.listendata.com/2014/11/random-forest-with-r.html
##https://github.com/jsescobar/westernization/blob/master/curated_metagenomes.R


library(fdrtool)
library(curatedMetagenomicData)
library("phyloseq")
library(vegan)
library(microbiome)
library("ggplot2")
library(plyr)

#Download count data
Asnicar_est = curatedMetagenomicData("AsnicarF_2017.metaphlan_bugs_list.stool",counts = TRUE, dryrun = FALSE)[[1]]
Brito_est = curatedMetagenomicData("BritoIL_2016.metaphlan_bugs_list.stool",counts = TRUE, dryrun = FALSE)[[1]]
Feng_est = curatedMetagenomicData("FengQ_2015.metaphlan_bugs_list.stool",counts = TRUE, dryrun = FALSE)[[1]]
Heitz_Buschart_est = curatedMetagenomicData("Heitz-BuschartA_2016.metaphlan_bugs_list.stool", counts = TRUE, dryrun = FALSE)[[1]]
HMP_est = curatedMetagenomicData("HMP_2012.metaphlan_bugs_list.stool",counts = TRUE, dryrun = FALSE)[[1]]
Karlsson_est = curatedMetagenomicData("KarlssonFH_2013.metaphlan_bugs_list.stool",counts = TRUE, dryrun = FALSE)[[1]]
LeChatelier_est = curatedMetagenomicData("LeChatelierE_2013.metaphlan_bugs_list.stool",counts = TRUE, dryrun = FALSE)[[1]]
Liu_est = curatedMetagenomicData("LiuW_2016.metaphlan_bugs_list.stool",counts = TRUE, dryrun = FALSE)[[1]]
Loman_est = curatedMetagenomicData("LomanNJ_2013.metaphlan_bugs_list.stool",counts = TRUE, dryrun = FALSE)[[1]]
Nielsen_est = curatedMetagenomicData("NielsenHB_2014.metaphlan_bugs_list.stool",counts = TRUE, dryrun = FALSE)[[1]]
Obregon_Tito_est = curatedMetagenomicData("Obregon-TitoAJ_2015.metaphlan_bugs_list.stool",counts = TRUE, dryrun = FALSE)[[1]]
Pasolli_est = curatedMetagenomicData("PasolliE_2018.metaphlan_bugs_list.stool",counts = TRUE, dryrun = FALSE)[[1]]
QinJ_est = curatedMetagenomicData("QinJ_2012.metaphlan_bugs_list.stool",counts = TRUE, dryrun = FALSE)[[1]]
QinN_est = curatedMetagenomicData("QinN_2014.metaphlan_bugs_list.stool",counts = TRUE, dryrun = FALSE)[[1]]
Rampelli_est = curatedMetagenomicData("RampelliS_2015.metaphlan_bugs_list.stool",counts = TRUE, dryrun = FALSE)[[1]]
Raymond_est = curatedMetagenomicData("RaymondF_2016.metaphlan_bugs_list.stool",counts = TRUE, dryrun = FALSE)[[1]]
Schirmer_est = curatedMetagenomicData("SchirmerM_2016.metaphlan_bugs_list.stool",counts = TRUE, dryrun = FALSE)[[1]]
Vatanen_est = curatedMetagenomicData("VatanenT_2016.metaphlan_bugs_list.stool",counts = TRUE, dryrun = FALSE)[[1]]
Vincent_est = curatedMetagenomicData("VincentC_2016.metaphlan_bugs_list.stool",counts = TRUE, dryrun = FALSE)[[1]]
Vogtmann_est = curatedMetagenomicData("VogtmannE_2016.metaphlan_bugs_list.stool",counts = TRUE, dryrun = FALSE)[[1]]
Xie_est = curatedMetagenomicData("XieH_2016.metaphlan_bugs_list.stool",counts = TRUE, dryrun = FALSE)[[1]]
Yu_est = curatedMetagenomicData("YuJ_2015.metaphlan_bugs_list.stool",counts = TRUE, dryrun = FALSE)[[1]]
Zeller_est = curatedMetagenomicData("ZellerG_2014.metaphlan_bugs_list.stool",counts = TRUE, dryrun = FALSE)[[1]]

#Convert data to Phyloseq format
Asnicar = ExpressionSet2phyloseq(Asnicar_est)
Brito = ExpressionSet2phyloseq(Brito_est)
Feng = ExpressionSet2phyloseq(Feng_est)
Heitz_Buschart = ExpressionSet2phyloseq(Heitz_Buschart_est)
HMP = ExpressionSet2phyloseq(HMP_est)
Karlsson = ExpressionSet2phyloseq(Karlsson_est)
LeChatelier = ExpressionSet2phyloseq(LeChatelier_est)
Liu = ExpressionSet2phyloseq(Liu_est)
Loman = ExpressionSet2phyloseq(Loman_est)
Nielsen = ExpressionSet2phyloseq(Nielsen_est)
Obregon_Tito = ExpressionSet2phyloseq(Obregon_Tito_est)
Pasolli = ExpressionSet2phyloseq(Pasolli_est)
QinJ = ExpressionSet2phyloseq(QinJ_est)
QinN = ExpressionSet2phyloseq(QinN_est)
Rampelli = ExpressionSet2phyloseq(Rampelli_est)
Raymond = ExpressionSet2phyloseq(Raymond_est)
Schirmer = ExpressionSet2phyloseq(Schirmer_est)
Vatanen = ExpressionSet2phyloseq(Vatanen_est)
Vincent = ExpressionSet2phyloseq(Vincent_est)
Vogtmann = ExpressionSet2phyloseq(Vogtmann_est)
Xie = ExpressionSet2phyloseq(Xie_est)
Yu = ExpressionSet2phyloseq(Yu_est)
Zeller = ExpressionSet2phyloseq(Zeller_est)


Asnicar_clean  = subset_samples(Asnicar,
                                disease == 'healthy' &
                                  is.na(antibiotics_current_use) &
                                  age_category=='adult')

Brito_clean  = subset_samples(Brito,
                              disease == 'healthy' &
                                is.na(antibiotics_current_use) &
                                age_category=='adult')

Feng_clean  = subset_samples(Feng,
                             is.na(disease) &
                               is.na(antibiotics_current_use) &
                               age_category=='adult')

Heitz_Buschart_clean  = subset_samples(Heitz_Buschart,
                                       is.na(disease) &
                                         antibiotics_current_use=='no' &
                                         age_category=='adult')

HMP_clean  = subset_samples(HMP,
                            disease=='healthy' &
                              is.na(antibiotics_current_use) &
                              age_category=='adult')

LeChatelier_clean  = subset_samples(LeChatelier,
                                    study_condition ='control' &
                                      antibiotics_current_use=='no' &
                                      age_category=='adult')
Liu_nW_clean  = subset_samples(Liu,
                               non_westernized=='yes' &
                                 study_condition =='control' &
                                 is.na(antibiotics_current_use) &
                                 age_category=='adult')

Liu_W_clean  = subset_samples(Liu,
                              non_westernized=='no' &
                                study_condition =='control' &
                                is.na(antibiotics_current_use) &
                                age_category=='adult')

Nielsen_clean  = subset_samples(Nielsen,
                                disease=='healthy' &
                                  is.na(antibiotics_current_use) &
                                  age_category=='adult')

Obregon_Tito_W_clean  = subset_samples(Obregon_Tito,
                                       non_westernized=='no' &
                                         disease=='healthy' &
                                         is.na(antibiotics_current_use) &
                                         age_category=='adult')

Obregon_Tito_nW_clean  = subset_samples(Obregon_Tito,
                                        non_westernized=='yes' &
                                          disease=='healthy' &
                                          is.na(antibiotics_current_use) &
                                          age_category=='adult')

Pasolli_clean  = subset_samples(Pasolli,
disease=='healthy' &
  antibiotics_current_use=='no' &
  age_category=='adult')

QinJ_clean  = subset_samples(QinJ,
                             is.na(disease) &
                               is.na(antibiotics_current_use) &
                               age_category=='adult')

QinN_clean  = subset_samples(QinN,
                             is.na(disease) &
                               antibiotics_current_use=='no' &
                               age_category=='adult')

Rampelli_clean  = subset_samples(Rampelli,
                                 study_condition ='control' &
                                   antibiotics_current_use=='no' &
                                   age_category=='adult')

Raymond_clean  = subset_samples(Raymond,
                                disease=='healthy' &
                                  antibiotics_current_use=='no' &
                                  age_category=='adult')

Schirmer_clean  = subset_samples(Schirmer,
                                 disease=='healthy' &
                                   antibiotics_current_use=='no' &
                                   age_category=='adult')

Vogtmann_clean  = subset_samples(Vogtmann,
                                 is.na(disease) &
                                   is.na(antibiotics_current_use) &
                                   age_category=='adult')

Xie_clean  = subset_samples(Xie,
                            disease=='healthy' &
                              is.na(antibiotics_current_use) &
                              age_category=='adult')

Zeller_clean  = subset_samples(Zeller,
                               disease=='healthy' &
                                 is.na(antibiotics_current_use) &
                                 age_category=='adult')

# Merge all clean datasets

merged_datasets = merge_phyloseq(Asnicar_clean, Brito_clean, Feng_clean,
                                 Heitz_Buschart_clean, HMP_clean,
                                 Liu_nW_clean, Liu_W_clean, Nielsen_clean,
                                 Obregon_Tito_W_clean, Obregon_Tito_nW_clean,
                                 QinJ_clean, QinN_clean, Rampelli_clean,
                                 Raymond_clean, Schirmer_clean, Vogtmann_clean,
                                 Xie_clean, Zeller_clean, Pasolli_clean) ## 3,514 taxa and 1,763 samples


###Subset cMD data to 8 genera
pv_phy<- subset_taxa(merged_datasets, Genus == "Prevotella")
tr_phy<- subset_taxa(merged_datasets, Genus == "Treponema")
bf_phy<- subset_taxa(merged_datasets, Genus == "Bifidobacterium")
bn_phy<- subset_taxa(merged_datasets, Genus == "Barnesiella")
ak_phy<- subset_taxa(merged_datasets, Genus == "Akkermansia")
al_phy<- subset_taxa(merged_datasets, Genus == "Alistipes")
ba_phy<- subset_taxa(merged_datasets, Genus == "Bacteroides")
pa_phy<- subset_taxa(merged_datasets, Genus == "Parabacteroides")


cmd_physeq<-merge_phyloseq(pv_phy,tr_phy,bf_phy,bn_phy,ak_phy,al_phy,ba_phy,pa_phy) ## 289 taxa and 1,763 samples
cmd_physeq_w<-subset_samples(cmd_physeq, non_westernized == 'no') ##1,433 samples
cmd_physeq_nw<-subset_samples(cmd_physeq, non_westernized == 'yes') ##330 samples

##Load agt_swt_sub (South African dataset) unfiltered data and subset 8 genera
pv_as_phy<- subset_taxa(agt_swt_sub, Genus == "Prevotella")
tr_as_phy<- subset_taxa(agt_swt_sub, Genus == "Treponema")
bf_as_phy<- subset_taxa(agt_swt_sub, Genus == "Bifidobacterium")
bn_as_phy<- subset_taxa(agt_swt_sub, Genus == "Barnesiella")
ak_as_phy<- subset_taxa(agt_swt_sub, Genus == "Akkermansia")
al_as_phy<- subset_taxa(agt_swt_sub, Genus == "Alistipes")
ba_as_phy<- subset_taxa(agt_swt_sub, Genus == "Bacteroides")
pa_as_phy<- subset_taxa(agt_swt_sub, Genus == "Parabacteroides")

as_physeq<-merge_phyloseq(pv_as_phy,tr_as_phy,bf_as_phy,bn_as_phy,ak_as_phy,al_as_phy,ba_as_phy,pa_as_phy) ##2,164 taxa in 167 samples

as_physeq_tax_glom<-tax_glom(as_physeq, "Genus")
#filt_as_physeq <-filter_taxa(as_physeq, function(x) mean(x) > 0, TRUE)#Only keep ASVs/taxa with mean greater than zero.
##save phyloseq objects (write_phyloseq --microbiome package) and re-import to phyloseq object

=====================================================================================================================================================
##Edit cMD data to remove extra rows and columns with strain and species data - repetitive abundance values
cmd_rf<- read_csv2phyloseq(
  otu.file = "/Users/voke/Desktop/16S_Analysis/Pilot_Data/june_2020_revision/revfinal/orig_pool/final_files/cmd/cmd_otu_table_ed.csv",
  taxonomy.file = "/Users/voke/Desktop/16S_Analysis/Pilot_Data/june_2020_revision/revfinal/orig_pool/final_files/cmd/cmd_taxonomy_table_ed.csv",
  metadata.file = "/Users/voke/Desktop/16S_Analysis/Pilot_Data/june_2020_revision/revfinal/orig_pool/final_files/cmd/cmd_metadata_table.csv",
  sep = ","
)##8 taxa and 1763 samples

as_rf<- read_csv2phyloseq(
  otu.file = "/Users/voke/Desktop/16S_Analysis/Pilot_Data/june_2020_revision/revfinal/orig_pool/final_files/cmd/as_glom_otu_table_ed.csv",
  taxonomy.file = "/Users/voke/Desktop/16S_Analysis/Pilot_Data/june_2020_revision/revfinal/orig_pool/final_files/cmd/as_glom_taxonomy_table_ed.csv",
  metadata.file = "/Users/voke/Desktop/16S_Analysis/Pilot_Data/june_2020_revision/revfinal/orig_pool/final_files/cmd/as_glom_metadata_table.csv",
  sep = ","
)##8 taxa and 167 samples

#merge re-imported cmd data(cmd_rf) with as_physeq data
as_cmd_physeq<-merge_phyloseq(cmd_rf,as_rf)

##Export as_cmd_physeq to edit metatable

cmd_as_rf<- read_csv2phyloseq(
  otu.file = "/Users/voke/Desktop/16S_Analysis/Pilot_Data/june_2020_revision/revfinal/orig_pool/final_files/cmd/cmd_as_otu_table.csv",
  taxonomy.file = "/Users/voke/Desktop/16S_Analysis/Pilot_Data/june_2020_revision/revfinal/orig_pool/final_files/cmd/cmd_as_taxonomy_table.csv",
  metadata.file = "/Users/voke/Desktop/16S_Analysis/Pilot_Data/june_2020_revision/revfinal/orig_pool/final_files/cmd/cmd_as_metadata_table_ed.csv",
  sep = ","
)

=====================================================================================================================================================

##The enclosed section could be done with a simple function:

#extract_taxa = function(dataset){
#  # Filter at a given taxonomic rank
#  genus_table = subset_taxa(dataset, !is.na(Genus) & is.na(Species))
#  # Extract the desired phylotypes
#  # Non-western
#  Prevotella = otu_table(subset_taxa(genus_table, Genus=="Prevotella"))
#  Treponema = otu_table(subset_taxa(genus_table, Genus=="Treponema"))
#  # Western
#  Bacteroides = otu_table(subset_taxa(genus_table, Genus=="Bacteroides"))
#  Bifidobacterium = otu_table(subset_taxa(genus_table, Genus=="Bifidobacterium"))
#  Barnesiella = otu_table(subset_taxa(genus_table, Genus=="Barnesiella"))
#  Alistipes = otu_table(subset_taxa(genus_table, Genus=="Alistipes"))
#  Akkermansia = otu_table(subset_taxa(genus_table, Genus=="Akkermansia"))
#  Parabacteroides = otu_table(subset_taxa(genus_table, Genus=="Parabacteroides"))


filt_cmd_as_rf <-filter_taxa(cmd_as_rf, function(x) mean(x) > 0, TRUE) #Only keep ASVs/taxa with mean greater than zero ==>  8 taxa and 1930 samples
ps_w_as <-subset_samples(filt_cmd_as_rf, non_westernized != 'Non-western')
ps_nw_as <-subset_samples(filt_cmd_as_rf, non_westernized != 'Western')

##Log-transform data

pslog_as_cmd<-transform_sample_counts(filt_cmd_as_rf, function(x) log(1 + x))
##Subset western samples to 167 to merge with RSA data:
pslog_w_cmd <- subset_samples(pslog_as_cmd, non_westernized == 'Western')
randomSample167_w <- sample(sample_names(pslog_w_cmd), size= 167, replace=FALSE)
pslog_w_167<-prune_samples(randomSample167_w, pslog_w_cmd)
pslog_rsa <- subset_samples(pslog_as_cmd, non_westernized == 'RSA')
pslog_w_as_167 <- merge_phyloseq(pslog_rsa,pslog_w_167)
##Get thesame random set from the untransformed data
ps_w_167 <-prune_samples(randomSample167_w, ps_w_as)
ps_rsa <- subset_samples(ps_w_as, non_westernized == 'RSA')
ps_w_as_167 <-merge_phyloseq(ps_rsa,ps_w_167)

##Non-western Samples
pslog_nw_cmd <- subset_samples(pslog_as_cmd, non_westernized == 'Non-western')
randomSample167_nw <- sample(sample_names(pslog_nw_cmd), size= 167, replace=FALSE)
pslog_nw_167<-prune_samples(randomSample167_nw, pslog_nw_cmd)
#pslog_rsa <- subset_samples(pslog_as_cmd, non_westernized == 'RSA')
pslog_nw_as_167 <- merge_phyloseq(pslog_rsa,pslog_nw_167)
##Get thesame random set from the untransformed data
ps_nw_167 <-prune_samples(randomSample167_nw, ps_nw_as)
ps_rsa <- subset_samples(ps_nw_as, non_westernized == 'RSA')
ps_nw_as_167 <-merge_phyloseq(ps_rsa,ps_nw_167)


dataMatrix_w <- data.frame(westernization = sample_data(pslog_w_as_167)$non_westernized, otu_table(t(pslog_w_as_167)))
trainingSample_w <- sample(rownames(sample_data(pslog_w_as_167)), size = 234, replace = FALSE)# take 234  (70 percent)individuals at random to be the training set, and the remaining 99 the test set
inTrain_w <- which(rownames(sample_data(pslog_w_as_167))%in% trainingSample_w)
training_w <- dataMatrix_w[inTrain_w,]
testing_w <- dataMatrix_w[-inTrain_w,]

##Random Forest 
rfFit_w <- train(westernization ~ ., data = training_w,
          method = "rf", preProc = "center")
rfClasses_w <- predict(rfFit_w, newdata = testing_w)
table(rfClasses_w, testing_w$westernization)

##Edited Random Forest plot function to allow for labeling x-axis
varImpPlot1<-function (x, sort = TRUE, n.var = min(30, nrow(x$importance)),
    type = NULL, class = NULL, scale = TRUE, main = deparse(substitute(x)),
    xlab=NULL,...)
{
    if (!inherits(x, "randomForest"))
        stop("This function only works for objects of class `randomForest'")
    imp <- importance(x, class = class, scale = scale, type = type,
        ...)
    if (ncol(imp) > 2)
        imp <- imp[, -(1:(ncol(imp) - 2))]
    nmeas <- ncol(imp)
    if (nmeas > 1) {
        op <- par(mfrow = c(1, 2), mar = c(4, 5, 4, 1), mgp = c(2,
            0.8, 0), oma = c(0, 0, 2, 0), no.readonly = TRUE)
        on.exit(par(op))
    }
    for (i in 1:nmeas) {
        ord <- if (sort)
            rev(order(imp[, i], decreasing = TRUE)[1:n.var])
        else 1:n.var
        xmin <- if (colnames(imp)[i] %in% c("IncNodePurity",
            "MeanDecreaseGini"))
            0
        else min(imp[ord, i])
        if(is.null(xlab)){
        xlab = colnames(imp)[i]
}
        dotchart(imp[ord, i], xlab = xlab, ylab = "",
            main = if (nmeas == 1)
                main
            else NULL, xlim = c(xmin, max(imp[, i])), ...)
    }
    if (nmeas > 1)
        mtext(outer = TRUE, side = 3, text = main, cex = 1.2)
    invisible(imp)
}


##Variable importance plot --use edited function (varImpPlot1) to allow for labeling xlab
western_vip<-varImpPlot1(rfFit_w$finalModel, main = "Western-RSA Comparison", xlab = "Variable Importance Scores")

##For non-western:
dataMatrix_nw <- data.frame(westernization = sample_data(pslog_nw_as_167)$non_westernized, otu_table(t(pslog_nw_as_167)))
trainingSample_nw <- sample(rownames(sample_data(pslog_nw_as_167)), size = 234, replace = FALSE)# take 234  (70 percent)individuals at random to be the training set, and the remaining 99 the test set
inTrain_nw <- which(rownames(sample_data(pslog_nw_as_167))%in% trainingSample_nw)
training_nw <- dataMatrix_nw[inTrain_nw,]
testing_nw <- dataMatrix_nw[-inTrain_nw,]

##Random Forest
rfFit_nw <- train(westernization ~ ., data = training_nw,
          method = "rf", preProc = "center")
rfClasses_nw <- predict(rfFit_nw, newdata = testing_nw)
table(rfClasses_nw, testing_nw$westernization)

##Variable importance plot --use edited function (varImpPlot1) to allow for labeling xlab
non_western_vip<-varImpPlot1(rfFit_nw$finalModel, main = "Non-western-RSA Comparison", xlab = "Variable Importance Scores")

##For Western and non-Western Comparisons:
ps_w_nw_167<-merge_phyloseq(ps_w_167, ps_nw_167)
pslog_w_nw_167<-merge_phyloseq(pslog_w_167, pslog_nw_167)
dataMatrix_w_nw <- data.frame(westernization = sample_data(pslog_w_nw_167)$non_westernized, otu_table(t(pslog_w_nw_167)))
trainingSample_w_nw <- sample(rownames(sample_data(pslog_w_nw_167)), size = 234, replace = FALSE)# take 234  (70 percent)individuals at random to be the training set, and the remaining 99 the test set
inTrain_w_nw <- which(rownames(sample_data(pslog_w_nw_167))%in% trainingSample_w_nw)
training_w_nw <- dataMatrix_nw[inTrain_w_nw,]
testing_w_nw <- dataMatrix_w_nw[-inTrain_w_nw,]
##RF      
rfFit_w_nw <- train(westernization ~ ., data = training_w_nw,
          method = "rf", preProc = "center")
rfClasses_w_nw <- predict(rfFit_w_nw, newdata = testing_w_nw)
table(rfClasses_w_nw, testing_w_nw$westernization)

##Variable importance plot --use edited function (varImpPlot1) to allow for labeling xlab
west_non_western_vip<-varImpPlot1(rfFit_w_nw$finalModel, main = "Western-non-Western Comparison", xlab = "Variable Importance Scores")

##Save plot as grid
png("fig6_randomforest_VIP2.png", units="mm", width=250, height=200, res=300)
par(mfrow=c(1,3))
 varImpPlot1(rfFit_w$finalModel, main = "Western-RSA Comparison", xlab = "Variable Importance Scores")
 varImpPlot1(rfFit_nw$finalModel, main = "non-Western-RSA Comparison", xlab ="Variable Importance Scores")
 varImpPlot1(rfFit_w_nw$finalModel, main = "Western-non-Western Comparison", xlab = "Variable Importance Scores")
 dev.off()

abundance_ps_w_as_167<- ps_w_as_167 %>% transform_sample_counts(function(x) {x/sum(x)} ) %>% psmelt() %>% arrange(Genus)

filt_ps_w_as_167<-abundance_ps_w_as_167%>% filter(Abundance > 0)
all_filt_ps_w_as_167<- filt_ps_w_as_167 %>%
  select(Genus, non_westernized, country, Abundance) %>%
  group_by(Genus, non_westernized)
head(all_filt_ps_w_as_167)

tr_as_w<- all_filt_ps_w_as_167 %>% filter(Genus == "Treponema")
bf_as_w<- all_filt_ps_w_as_167 %>% filter(Genus == "Bifidobacterium")
bn_as_w<- all_filt_ps_w_as_167 %>% filter(Genus == "Barnesiella")
ak_as_w<- all_filt_ps_w_as_167 %>% filter(Genus == "Akkermansia")
al_as_w<- all_filt_ps_w_as_167 %>% filter(Genus == "Alistipes")
ba_as_w<- all_filt_ps_w_as_167 %>% filter(Genus == "Bacteroides")
pa_as_w<- all_filt_ps_w_as_167 %>% filter(Genus == "Parabacteroides")

##Kruskal-Wallis test
kw_pv_as_w<-  kruskal.test(Abundance ~ non_westernized, data =pv_as_w)  #p-value < 2.2e-16
kw_tr_as_w<-  kruskal.test(Abundance ~ non_westernized, data =tr_as_w)  #p-value = 0.6299
kw_bf_as_w<-  kruskal.test(Abundance ~ non_westernized, data =bf_as_w)  #p-value = 3.087e-16
kw_bn_as_w<-  kruskal.test(Abundance ~ non_westernized, data =bn_as_w)  #p-value = 0.1389
kw_ak_as_w<-  kruskal.test(Abundance ~ non_westernized, data =ak_as_w)  #p-value = 0.9392
kw_al_as_w<-  kruskal.test(Abundance ~ non_westernized, data =al_as_w)  #p-value = 2.2e-16
kw_ba_as_w<-  kruskal.test(Abundance ~ non_westernized, data =ba_as_w)  #p-value = 1.076e-08
kw_pa_as_w<-  kruskal.test(Abundance ~ non_westernized, data =pa_as_w)  #p-value = 7.207e-05

##Summary statistics for fold change calculation
pv_as_w_summary<-pv_as_w_summary<- by(pv_as_w$Abundance, pv_as_w$non_westernized, summary)
tr_as_w_summary<- by(tr_as_w$Abundance, tr_as_w$non_westernized, summary)
bf_as_w_summary<- by(bf_as_w$Abundance, bf_as_w$non_westernized, summary)
bn_as_w_summary<- by(bn_as_w$Abundance, bn_as_w$non_westernized, summary)
ak_as_w_summary<- by(ak_as_w$Abundance, ak_as_w$non_westernized, summary)
al_as_w_summary<- by(al_as_w$Abundance, al_as_w$non_westernized, summary)
ba_as_w_summary<- by(ba_as_w$Abundance, ba_as_w$non_westernized, summary)
pa_as_w_summary<- by(pa_as_w$Abundance, pa_as_w$non_westernized, summary)

##Non-western:
abundance_ps_nw_as_167<- ps_nw_as_167 %>% transform_sample_counts(function(x) {x/sum(x)} ) %>% psmelt() %>% arrange(Genus)

filt_ps_nw_as_167<-abundance_ps_nw_as_167%>% filter(Abundance > 0)
all_filt_ps_nw_as_167<- filt_ps_nw_as_167 %>%
  select(Genus, non_westernized, country, Abundance) %>%
  group_by(Genus, non_westernized)
head(all_filt_ps_nw_as_167)

pv_as_nw<- all_filt_ps_nw_as_167 %>% filter(Genus == "Prevotella")
tr_as_nw<- all_filt_ps_nw_as_167 %>% filter(Genus == "Treponema")
bf_as_nw<- all_filt_ps_nw_as_167 %>% filter(Genus == "Bifidobacterium")
bn_as_nw<- all_filt_ps_nw_as_167 %>% filter(Genus == "Barnesiella")
ak_as_nw<- all_filt_ps_nw_as_167 %>% filter(Genus == "Akkermansia")
al_as_nw<- all_filt_ps_nw_as_167 %>% filter(Genus == "Alistipes")
ba_as_nw<- all_filt_ps_nw_as_167 %>% filter(Genus == "Bacteroides")
pa_as_nw<- all_filt_ps_nw_as_167 %>% filter(Genus == "Parabacteroides")


##Kruskal-Wallis test
kw_pv_as_nw<-  kruskal.test(Abundance ~ non_westernized, data =pv_as_nw)  #p-value = 0.00102
kw_tr_as_nw<-  kruskal.test(Abundance ~ non_westernized, data =tr_as_nw)  #p-value = 0.5182
kw_bf_as_nw<-  kruskal.test(Abundance ~ non_westernized, data =bf_as_nw)  #p-value = 2.2e-16
kw_bn_as_nw<-  kruskal.test(Abundance ~ non_westernized, data =bn_as_nw)  #p-value = 6.279e-14
kw_ak_as_nw<-  kruskal.test(Abundance ~ non_westernized, data =ak_as_nw)  #p-value = 0.009197
kw_al_as_nw<-  kruskal.test(Abundance ~ non_westernized, data =al_as_nw)  #p-value = 0.00359
kw_ba_as_nw<-  kruskal.test(Abundance ~ non_westernized, data =ba_as_nw)  #p-value = 2.2e-16
kw_pa_as_nw<-  kruskal.test(Abundance ~ non_westernized, data =pa_as_nw)  #p-value = 2.2e-16

##Summary statistics for fold change calculation
pv_as_nw_summary<-pv_as_nw_summary<- by(pv_as_nw$Abundance, pv_as_nw$non_westernized, summary)
tr_as_nw_summary<- by(tr_as_nw$Abundance, tr_as_nw$non_westernized, summary)
bf_as_nw_summary<- by(bf_as_nw$Abundance, bf_as_nw$non_westernized, summary)
bn_as_nw_summary<- by(bn_as_nw$Abundance, bn_as_nw$non_westernized, summary)
ak_as_nw_summary<- by(ak_as_nw$Abundance, ak_as_nw$non_westernized, summary)
al_as_nw_summary<- by(al_as_nw$Abundance, al_as_nw$non_westernized, summary)
ba_as_nw_summary<- by(ba_as_nw$Abundance, ba_as_nw$non_westernized, summary)
pa_as_nw_summary<- by(pa_as_nw$Abundance, pa_as_nw$non_westernized, summary)

