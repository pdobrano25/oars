### 2025_11_04  ML_Redux-Deux

# Goal: make the most meaningful predictor

# i.e.g AUCs for predicting delta-pH
# and specifically +/- 1.27 (or, from -1.18 to -1.52)

# 11_07 version uses greengenes2
# 11_05 version uses greengenes13.8

# save.image(file =  "./2025_11_07_oars_ml_env.Rdata")

load(file =  "./ml_git_data/ml_loocv_data/2025_11_07_oars_ml_env.Rdata")


library("ggplot2"); library("ggridges"); library("dplyr"); library("tidyverse")

setwd("./ml_git_data/ml_loocv_data")

# set vectors for RS names +/- PBS
rs.names.pbs <- c("PBS", "Authentic", "BobsRedMill", "MSPrebiotic", "LetsDoOrganic", "HiMaize260", "Novelose330", "ActistarRT", "FibersymRW", "Versafibe1490")
rs.names <- c("Authentic", "BobsRedMill", "MSPrebiotic", "LetsDoOrganic", "HiMaize260", "Novelose330", "ActistarRT", "FibersymRW", "Versafibe1490")

# give 9 RS consistent colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
labelcolors <- data.frame(cols = c(gg_color_hue(9), "#000000"),
                          label = c(seq(1:9), "P"))
cols = c(gg_color_hue(5), "#000000")
labelcolors = data.frame(cols = c(
  cols[1], # authentic (red)
  cols[1], # bobsredmill (red)
  cols[1], # msprebiotic (red)
  cols[2], # letsdoorgnaic (gold)
  cols[3], # himaize (green)
  cols[3], # novelose (green)
  cols[4], # actistar (blue)
  cols[4], # fibersym (green)
  cols[5], # versafibe (purple)
  "#000000" #pbs
))

labelcolors.rs =labelcolors$cols[c(1:9)]
names(labelcolors.rs) = rs.names


# load rapidaim data ------------------------------------------------------

# load rapidaim asv data (median)
ml.rapidaim.median.data <- readRDS("../../ml_git_data/2025_11_07_n117_s25_asv_final_median_glom_50k.rds")
# subset to RapidAIM samples for analysis
ml.rapidaim.median.data <- ml.rapidaim.median.data[!grepl("Stool", rownames(ml.rapidaim.median.data)),]
# remove FOS
ml.rapidaim.median.data <- ml.rapidaim.median.data[!grepl("FOS", rownames(ml.rapidaim.median.data)),]
# good

# load rapidaim asv data (all replicates)
ml.rapidaim.data <- readRDS("../../ml_git_data/2025_11_07_n117_s25_asv_final_glom_50k.rds")
# subset to RapidAIM samples for analysis
ml.rapidaim.data <- phyloseq::subset_samples(ml.rapidaim.data, RS_Name %in% rs.names.pbs)
# melt
ml.rapidaim.data <- speedyseq::psmelt(ml.rapidaim.data)
ml.rapidaim.data$code <- paste(ml.rapidaim.data$HM, ml.rapidaim.data$RS_Name, ml.rapidaim.data$Replicate, sep="_")
# cast
ml.rapidaim.data <- reshape2::acast(ml.rapidaim.data, code ~ Glom, value.var="Abundance")
dim(ml.rapidaim.data)
# 3430 samples x 1656 taxa



# load mapping file -------------------------------------------------------

# load mapping file // defunct, but mapping file is fine
# (local, with covariates:)
# ml.rapidaim.mapping <- read.csv(paste0(home.dir,"ml_git_data/2024_04_15_full_mapping_rapidaim_covariates.csv"), sep=",")
# ml.rapidaim.phyloseq <- readRDS(paste0(home.dir,"ml_git_data/2024_09_13_n117_asv_final_50k.rds", sep=""))
# ml.rapidaim.mapping <- subset(ml.rapidaim.mapping, Barcode %in% phyloseq::sample_data(ml.rapidaim.phyloseq)$Barcode)
# # make rownames HM_RS_Replicate
# ml.rapidaim.mapping$code <- paste(ml.rapidaim.mapping$HM, ml.rapidaim.mapping$RS_Name, ml.rapidaim.mapping$Replicate, sep="_")
# # first, remove duplicated samples
# ml.rapidaim.mapping <- subset(ml.rapidaim.mapping, !code %in% c("HM0912.00_RawStool_A","HM0831.00_WaterB_C"))
# # apply rownames
# rownames(ml.rapidaim.mapping) <- ml.rapidaim.mapping$code
# unique(phyloseq::sample_data(ml.rapidaim.phyloseq)$HM_noR) %in%  unique(ml.rapidaim.mapping$HM_noR)
# saveRDS(ml.rapidaim.mapping, paste0(home.dir,"ml_git_data/2024_09_13_n117_rapidaim_mapping_covariates.rds", sep=""))

ml.rapidaim.mapping <- readRDS("../../ml_git_data/2024_09_13_n117_rapidaim_mapping_covariates.rds")


# :: load raw pH ----------------------------------------------------------

# raw pH values (not delta)

ph.data <- readRDS("../2024_03_13_n117_asv_final_mapping.rds")

# subset to important RS
ph.data <- subset(ph.data, RS_Name %in% c(rs.names, "PBS"))

# calculate median value
ph.data <- ph.data %>%
  group_by(HM, RS_Name) %>%
  mutate(med.ph = median(na.omit(pH))) %>% 
  dplyr::select(HM, RS_Name, med.ph) %>% distinct() %>% data.frame()

saveRDS(ph.data, "../2025_11_07_ml_ph_target.Rds")

# load ml data ----------------------------------------------------------

# load ML ASV data:
ml.input.data <- readRDS("../../ml_git_data/2025_11_07_n117_s25_asv_final_median_glom_50k.rds")

# subset to (slurry) stool samples for ml input (and remove raw stool)
ml.input.data <- ml.input.data[grepl("Stool", rownames(ml.input.data)),]

ml.input.data <- ml.input.data[!grepl("RawStool", rownames(ml.input.data)),]
print(paste0("pre-filtering: ", dim(ml.input.data)[2], " features"))
# good

# remove "_Stool" from rownames
rownames(ml.input.data) <- gsub("\\_Stool", "", rownames(ml.input.data))

# confirm ml.input.data and ml.target.data have the same HMs
sum(rownames(ml.input.data) %in% unique(ml.target.data$HM)) / nrow(ml.input.data) * 100
# good

## check sparsity
ml.input.data.sparsity = ml.input.data
ml.input.data.sparsity[ml.input.data.sparsity > 0] <- 1
# check column sparsity
# colSums(ml.input.data.sparsity)
# remove 0 columns
ml.input.data <- ml.input.data[,colSums(ml.input.data.sparsity)!=0]
dim(ml.input.data)
print(paste0("sparsity: ", round(sum(ml.input.data) / length(ml.input.data[ml.input.data==0]), digits=2), "%"))
# ~46% zeroes (vs ~90% without removing 0 columns)
# reduced from 1656 features to 905

saveRDS(ml.input.data, "../2025_11_07_n117_s25_asv_input_data.rds")

print(paste0("post-filtering: ", dim(ml.input.data)[2], " features"))


# 0. DESCRIPTIVE ----------------------------------------------------------------

# plot pH values

# stats; technical replicates --> use lme; alternative would be to create 1 model, but this already yields padj = 0 for all RS
ml.rapidaim.ph.stats = do.call(rbind, lapply(rs.names, function(x){
  #x="Authentic"
  data.ph.all.subset = subset(ml.rapidaim.mapping, RS_Name %in% c(x, "PBS"))
  data.ph.all.subset.stats = lmerTest::lmer(pH ~ RS_Name + age + gender + diagnosis + STUDY + (1|HM), data.ph.all.subset)
  data.ph.all.subset.stats = data.frame(RS_Name = x,
                                        total.rsq = rsq::rsq.lmm(data.ph.all.subset.stats)$model,
                                        fixed.rsq = rsq::rsq.lmm(data.ph.all.subset.stats)$fixed,
                                        pval = coef(summary(data.ph.all.subset.stats))[2,5])
  data.ph.all.subset.stats
}))
# Collate, apply FDR, and assign significance levels
ml.rapidaim.ph.stats = ml.rapidaim.ph.stats %>%
  mutate(padj = p.adjust(pval, method="bonferroni")) %>%
  mutate(sig = ifelse(padj < 0.001, "***", ifelse(padj < 0.01, "**", ifelse(padj < 0.05, "*", "")))) %>%
  mutate(RS_Name = factor(RS_Name, levels=rs.names))

# :: :: PLOT: Δ pH ----------------------------------------------

ml.rapidaim.ph.plot <-  ggplot(ml.target.data %>% subset(variable == "delta_pH") %>%
                                 mutate(target = "Δ pH"),
                               aes(x=(RS_Name), y=(value), fill=(RS_Name)))+
  ggbeeswarm::geom_beeswarm(shape=21, cex=0.75,
                            aes(fill=RS_Name), size=2)+
  geom_boxplot(notch = TRUE, outlier.shape=NA, aes(fill=RS_Name), alpha=0.3)+
  scale_fill_manual(values=labelcolors$cols[c(1:9)])+
  scale_color_manual(values=labelcolors$cols[c(1:9)])+
  scale_alpha_manual(values=c(0.3,1))+
  scale_y_continuous(limits=c(-2.9, 0.6))+
  geom_text(data=ml.rapidaim.ph.stats, aes(x=RS_Name, 
                                           label=ifelse(sig == "***", "*", "")), y=Inf, vjust=1.4, size=6)+
  theme_classic()+theme(legend.position="none",
                        axis.text.x=element_text(angle=45,vjust=1, hjust=1),
                        legend.title = element_text(hjust=0.5),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  facet_wrap(~"Measured Δ pH")+
  labs(x=NULL, y="Measured Δ pH versus PBS")
ml.rapidaim.ph.plot
# good


# :: :: dbRDA ----------------------------------------------------------------
# These axes are, “in successive order, a series of linear combination 
# of the explanatory variables that best explain the variation of 
# the response matrix”

# plot dbRDA
ml.dbrda.data.ph = reshape2::acast(subset(ml.target.data, variable =="delta_pH"), 
                                   HM ~ RS_Name, value.var="value")
# perform Bray-Curtis on stool microbiome
ml.dbrda.data.bray = vegan::vegdist(ml.input.data, method="bray")

# proper order
ml.dbrda.data.ph = ml.dbrda.data.ph[rownames(ml.input.data),]

nrow(ml.dbrda.data.ph)
nrow(ml.input.data)

# distance-based-RDA
dbrda_model_ph <- vegan::dbrda(ml.dbrda.data.bray ~ ., data = as.data.frame(ml.dbrda.data.ph))

# sig
set.seed(25)
dbrda_model_ph_permanova = anova(dbrda_model_ph, by="margin")
dbrda_model_ph_permanova
# very sig

# ggplot
dbrda_model_ph_scores <- vegan::scores(dbrda_model_ph, scaling = 2)

sites_ph <- as.data.frame(dbrda_model_ph_scores$sites)
sites_ph$Sample <- rownames(sites_ph)

# Variable scores (arrows for environmental variables)
vars_ph <- as.data.frame(dbrda_model_ph_scores$biplot)
vars_ph$RS_Name <- rownames(vars_ph)

dbrda_model_ph_permanova.df = data.frame(dbrda_model_ph_permanova) %>% 
  mutate(RS_Name = rownames(.))%>%
  subset(RS_Name != "Residual") %>%
  mutate(sig = ifelse(`Pr..F.` < 0.05, "*", "")) %>%
  mutate(boldness = ifelse(sig == "*", "bold", "plain"))

vars_ph = merge(vars_ph,
                dbrda_model_ph_permanova.df, by="RS_Name")
vars_ph$rs_sig = ifelse(vars_ph$sig == "*", 
                        paste(vars_ph$RS_Name, vars_ph$sig, sep=""),
                        vars_ph$RS_Name)

# :: :: PLOT: Δ pH dbRDA ----------------------------------------------

dbrda_model_ph_plot = ggplot() +
  # Plot samples as points
  geom_point(data = sites_ph, aes(x = dbRDA1, y = dbRDA2), shape=21, fill = "white", size = 2.5, alpha = 1) +
  # Add arrows for variables
  geom_segment(data = vars_ph %>% mutate(RS_Name = factor(RS_Name, levels=rs.names)), 
               aes(x = 0, y = 0, xend = dbRDA1*4, yend = dbRDA2*4,
                   color=RS_Name),
               arrow = arrow(length = unit(0.2, "cm"), type="closed"), linewidth = 0.5) +
  ggnetwork::geom_nodetext_repel(data=vars_ph, aes(x=dbRDA1*4, y=dbRDA2*4, 
                                                   color=RS_Name, label=rs_sig,
                                                   fontface = boldness), size=3)+
  scale_color_manual(values=labelcolors$cols)+
  # Customize axes
  labs(x = paste0("Axis 1 (", round(dbrda_model_ph$CCA$eig[1]/sum(dbrda_model_ph$CCA$eig)*100, 1), "%)"),
       y = paste0("Axis 2 (", round(dbrda_model_ph$CCA$eig[2]/sum(dbrda_model_ph$CCA$eig)*100, 1), "%)"))+
  theme_classic() +
  theme(legend.position = "none",
        strip.text = element_text(size=10))+
  facet_wrap(~"Bray-Curtis distance-based RDA")
dbrda_model_ph_plot



# :: :: Classic stats -----------------------------------------------------

# which taxa are sig associated with response using classical stats
# (Maaslin3)

ml.input.data

ml.target.data

ml.input.data.for.wilcox = ml.input.data %>% as.data.frame() %>%
  mutate(HM = rownames(.)) %>%
  merge(subset(ml.target.data, variable=="delta_pH"), by="HM") %>%
  mutate(response = ifelse(value < -1.27, "strong", "weak")) 

ml.input.data.for.wilcox.output = do.call(rbind, lapply(colnames(ml.input.data), function(taxa){
  # present in N% of samples

  do.call(rbind, lapply(rs.names, function(rs){
    presence.absence = ml.input.data[,colnames(ml.input.data)==taxa]
    presence.absence[presence.absence!=0] = 1
    if(sum(presence.absence) / length(presence.absence) > 0.20){
      
    data.subset = ml.input.data.for.wilcox[,colnames(ml.input.data.for.wilcox) %in%
                                             c(taxa, "RS_Name", "value", "response")]
    colnames(data.subset)[1] = "taxa"
    data.subset = subset(data.subset, RS_Name == rs)
    test.output = wilcox.test(subset(data.subset, response == "strong")$taxa,
                              subset(data.subset, response == "weak")$taxa)
    # calculate LogFC
    pseudo = min(ml.input.data[ml.input.data!=0])/2
    lfc = log2(mean(subset(data.subset, response == "strong")$taxa+pseudo)/
                 mean(subset(data.subset, response == "weak")$taxa+pseudo))
    data.frame(RS_Name = rs,
               taxa = taxa,
               lfc = lfc,
               pval = test.output$p.value)
    }else{
      data.frame(RS_Name = rs,
                 taxa = taxa,
                 lfc = NA,
                 pval = NA)
    }
  }))
  }))

ml.input.data.for.wilcox.output$padj = p.adjust(ml.input.data.for.wilcox.output$pval, method="BH")
ml.input.data.for.wilcox.output$sig = ifelse(ml.input.data.for.wilcox.output$padj < 0.20, "*", "")
subset(ml.input.data.for.wilcox.output, sig == "*")

ml.wilcox.taxa.order = subset(ml.input.data.for.wilcox.output, sig == "*") %>% 
  group_by(taxa) %>% 
  mutate(sum.lfc = sum(lfc)) %>% dplyr::select(taxa, sum.lfc) %>% distinct() %>% arrange(sum.lfc)

# plot
ml.input.data.for.wilcox.output %>%
  subset(taxa %in% subset(ml.input.data.for.wilcox.output, sig == "*")$taxa)%>%
  mutate(RS_Name = factor(RS_Name, levels=rs.names)) %>%
  mutate(taxa = factor(taxa, levels=ml.wilcox.taxa.order$taxa))%>%
  ggplot(aes(x=RS_Name, y=taxa))+
  geom_tile(aes(fill=lfc), color="white")+
  scale_fill_gradient2(low="blue", mid="white", high="red")+
  theme_classic()

ml.wilcox.heatmap.data = ml.input.data.for.wilcox.output %>%
  subset(taxa %in% subset(ml.input.data.for.wilcox.output, sig == "*")$taxa)%>%
  mutate(RS_Name = factor(RS_Name, levels=rs.names)) %>%
  mutate(taxa = factor(taxa, levels=ml.wilcox.taxa.order$taxa))

reshape2::acast(ml.wilcox.heatmap.data,
                RS_Name ~ taxa, value.var="lfc") %>%
  t() %>%
  pheatmap::pheatmap(color = colorRampPalette(c("blue", "white", "red"))(100),
                     annotation_col = data.frame(RS = rs.names.pbs)%>% `rownames<-`(rs.names.pbs),
                     annotation_colors = list(RS = c(labelcolors.rs, "PBS" = "black")),
                     annotation_legend = F,
                     #clustering_distance_cols="correlation",
                     #clustering_distance_rows="correlation",
                     display_numbers = reshape2::acast(ml.wilcox.heatmap.data, taxa ~ RS_Name, value.var="sig")[,rs.names],
                     fontsize_number = 15,number_color="black",
                     breaks=c(seq(min(na.omit(ml.wilcox.heatmap.data$lfc)), 0, length.out=ceiling(100/2) + 1), 
                              seq(max(na.omit(ml.wilcox.heatmap.data$lfc))/100, max(na.omit(ml.wilcox.heatmap.data$lfc)), length.out=floor(100/2))),
                     border_color = "white") 
# color = LogFC of Strong versus Weak fermenters


 # 1. META-ANALYSIS --------------------------------------------------------


# :: :: objective ------------------------------------------------------------


# Chen et al conducted a [meta-analysis](https://pubs.rsc.org/en/content/articlelanding/2023/fo/d3fo00845b) examining the impact of resistant starch on the microbiome across 6 studies.
# I've downloaded SRA files and processed the fastq files through QIIME2 and DADA2.
# In "sra_16s_process.Rmd", I assembled phyloseq objects (with metadata), assigned taxonomy, and previewed their quality.
# Now, we will merge datasets, apply a prevalance/abundance filter to taxa, and conduct some preliminary analyses: Beta-diversity and Maaslin2.
# After, we will apply machine learning (leave-one-dataset-out cross-validation) to predict responses to RS using baseline stool compositions. Responses will be based on Butyrogens. A "responder" will be defined based on median or comparison to PBS. **Importantly**, since this is a validation project, I can apply feature selection using the taxa selected as important in my ML project.


# :: :: processing ---------------------------------------------------------

# save phyloseq + clean

venk.phyloseq.df = readRDS("../../chen_validation/processed_files/venk.phyloseq.gg2.Rds")
venk.meta.clean = readRDS("../../chen_validation/processed_files/venk.meta.clean.Rds")

dee.phyloseq.df = readRDS("../../chen_validation/processed_files/dee.phyloseq.gg2.Rds")
dee.meta.clean = readRDS("../../chen_validation/processed_files/dee.meta.clean.Rds")
dee.meta.clean = subset(dee.meta.clean, RS_Name != "Control")
dee.meta.clean$RS_Name = paste("Cross-linked ", dee.meta.clean$RS_Name, sep="")
# remove control

dem.phyloseq.df = readRDS("../../chen_validation/processed_files/dem.phyloseq.gg2.Rds")
dem.meta.clean = readRDS("../../chen_validation/processed_files/dem.meta.clean.Rds")

upa.phyloseq.df = readRDS("../../chen_validation/processed_files/upa.phyloseq.gg2.Rds")
upa.meta.clean = readRDS("../../chen_validation/processed_files/upa.meta.clean.Rds")

mai.phyloseq.df.mat = readRDS("../../chen_validation/processed_files/mai.phyloseq.gg2.Rds")
mai.meta.clean.2 = readRDS("../../chen_validation/processed_files/mai.meta.clean.Rds")

hug.phyloseq.df = readRDS("../../chen_validation/processed_files/hug.phyloseq.gg2.Rds")
hug.meta.clean = readRDS("../../chen_validation/processed_files/hug.meta.clean.Rds")

han.phyloseq.df = readRDS("../../chen_validation/processed_files/han.phyloseq.gg2.Rds")
han.meta.clean = readRDS("../../chen_validation/processed_files/han.meta.clean.Rds")

flo.phyloseq.df.mat = readRDS("../../chen_validation/processed_files/flo.phyloseq.gg2.Rds")
flo.meta.clean.2 = readRDS("../../chen_validation/processed_files/flo.meta.clean.Rds")

# Merge data into a single ASV table

unified.df = rbind(reshape2::melt(venk.phyloseq.df),
                   reshape2::melt(dee.phyloseq.df),
                   reshape2::melt(dem.phyloseq.df),
                   reshape2::melt(upa.phyloseq.df),
                   reshape2::melt(mai.phyloseq.df.mat),
                   reshape2::melt(hug.phyloseq.df),
                   reshape2::melt(han.phyloseq.df),
                   reshape2::melt(flo.phyloseq.df.mat)) %>% data.frame()
# cast
unified.df = reshape2::acast(unified.df,
                             Var1 ~ Var2, value.var="value") 
# replace NA with 0
unified.df[is.na(unified.df)] = 0

# also, merge metadata
unified.meta = rbind((venk.meta.clean),
                     (dee.meta.clean),
                     (dem.meta.clean),
                     (upa.meta.clean),
                     (mai.meta.clean.2),
                     (hug.meta.clean),
                     (han.meta.clean),
                     (flo.meta.clean.2)) %>% data.frame()
# make rownames = sample
rownames(unified.meta) = unified.meta$sample
nrow(unified.meta)
# 457 samples

# keep subjects with matching ASV data
unified.meta = subset(unified.meta, sample %in% rownames(unified.df))
# keep subjects with both timepoints
unified.meta = subset(unified.meta, 
                      subject %in% subset(data.frame(table(unified.meta$subject)), Freq >1)$Var1)

# check
# unified.meta$subject %>% table() %>% table()

unified.df = unified.df[unified.meta$sample,]
print(paste0(nrow(unified.df), " samples remain"))

rs.studies = sort(c("venkataraman", "deehan","demartino","hanes","hughes","flowers","maier","upadhyaya"))

unified.meta$group = paste(unified.meta$study,
                           unified.meta$RS_Name, sep="_")
rs.groups = sort(unique(unified.meta$group))


# assign RS_type to group

rs.types = data.frame(group = rs.groups,
                      RS_type = c("Green", "Purple", "Blue", "Red", "Red", "Gold", "Gold", "Gold", "Red", "Blue", "Green", "Green", "Blue", "Red"))
rs.types$RS_type = factor(rs.types$RS_type, levels=c("Red", "Gold", "Green", "Blue", "Purple"))
# assign colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
labelcolors.meta <- data.frame(cols = c(gg_color_hue(5)),
                               RS_type = c("Red", "Gold", "Green", "Blue", "Purple"))
rs.types = merge(rs.types,
                 labelcolors.meta, by="RS_type")
rs.types.vector = as.vector(rs.types$cols) 
names(rs.types.vector) = rs.types$group

rs.types$RS_type = factor(rs.types$RS_type, levels=c("Red", "Gold", "Green", "Blue", "Purple"))

# :: :: filter ---------------------------------------------------------------

unified.df.prev = unified.df

unified.df.prev = 
  # select prevalence thresholds
  do.call(rbind, lapply(seq(1:50), function(y){
    # select minimum abundance threshold
    do.call(rbind, lapply(c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1), function(x){
      unified.df.prev[unified.df.prev >= x] <- 1
      unified.df.prev[unified.df.prev < x] <- 0
      prev = colSums(unified.df.prev) 
      # tabulate number of taxa remaining after abundance + prevalence thresholds
      data.frame(prev = sum(prev >= round(y*nrow(unified.df.prev)/100, digits=0)),
                 abun_filt = x,
                 prev_filt = y)
    }))}))

ggplot(unified.df.prev,
       aes(x=(abun_filt*100), y=(prev)))+
  geom_point(aes(fill=prev_filt), shape=21, size=3)+
  scale_fill_gradient2(low=("blue"),
                       mid="white",midpoint=25,
                       high=("red"))+
  scale_x_log10()+
  scale_y_log10()+
  theme_minimal()+
  labs(x="Abundance % Threshold",
       y="Taxa Remaining",
       fill="Prevalence\nThreshold\n(% samples)")

# We can preserve  around 100 taxa if we apply a filter that keeps taxa "detected at 0.1% abundance in 10%+ samples"

unified.df.prev = unified.df

unified.df.prev[unified.df.prev >= 0.01] <- 1
unified.df.prev[unified.df.prev < 0.01] <- 0
prev = colSums(unified.df.prev) 

unified.df.filt = unified.df[,prev >= round(nrow(unified.df)*10/100, digits=0)]

print(paste0(ncol(unified.df.filt), " taxa remain"))

# This includes g__Ruminococcus_E_s__bromii_B
# B. adolescentis is no longer present with gg2
# But, B. faecale is indistinguishable from B. adolescentis using 16S


# :: MMUPHIn -------------------------------------------------------------

# meta
unified.meta
# asv
dim(unified.df.filt)

# Apply batch correction using MMUPHIn

# needs % data
unified.df.filt.1 = sweep(unified.df.filt, 1, rowSums(unified.df.filt), FUN = "/") * 1


unified.df.filt.badj = MMUPHin::adjust_batch(feature_abd=t(unified.df.filt.1),
                                             batch = "study",
                                             data=unified.meta)

# :: LASSO ----------------------------------------------------------------

# goal, use LASSO to calculate a response signature

# pseudo
pseudo.lasso = min(unified.df.filt.badj$feature_abd_adj[unified.df.filt.badj$feature_abd_adj!=0])/2

unified.df.lasso.data = t(unified.df.filt.badj$feature_abd_adj) %>% as.data.frame()
unified.df.lasso.data$sample = rownames(unified.df.lasso.data)
unified.df.lasso.data = unified.df.lasso.data %>%
  merge(unified.meta[,c("sample", "timepoint")], by="sample")%>%
  mutate(timepoint = ifelse(timepoint == "baseline", 0, 1)) 
rownames(unified.df.lasso.data) = unified.df.lasso.data$sample
unified.df.lasso.data$sample = NULL

set.seed(25)
cv_lasso = glmnet::cv.glmnet(
  y = unified.df.lasso.data$timepoint,
  x = scale(as.matrix(unified.df.lasso.data[,!colnames(unified.df.lasso.data) %in% c("timepoint")])),
  alpha = 1,          # 1 = LASSO, 0 = Ridge
  nfolds = 10,        # 10-fold CV
  type.measure = "mse"  # Mean squared error
)
plot(cv_lasso)

# Best lambda
lambda_min = cv_lasso$lambda.min
lasso_coefs = coef(cv_lasso, s = "lambda.min") %>% as.matrix() %>% as.data.frame() %>%
  mutate(feature = rownames(.)) %>%
  subset(feature != "(Intercept)") %>%
  subset(lambda.min !=0) %>%
  mutate(coef = lambda.min) %>%
  dplyr::select(-lambda.min) %>%
  arrange(coef)
lasso_coefs_sparse = coef(cv_lasso, s = "lambda.1se") %>% as.matrix() %>% as.data.frame() %>%
  mutate(feature = rownames(.)) %>%
  subset(feature != "(Intercept)") %>%
  subset(lambda.1se !=0) %>%
  mutate(coef = lambda.1se) %>%
  dplyr::select(-lambda.1se) %>%
  arrange(coef)

# plot features
lasso_coefs_plot = ggplot(lasso_coefs %>%
         mutate(direction = gsub(1, "Positive", gsub(-1, "Negative",  sign(coef)))) %>%
         mutate(direction = factor(direction, levels=rev(c("Negative", "Positive")))),
       aes(y=reorder(feature,coef), x=abs(coef)))+
  geom_vline(xintercept=0, linetype=2, linewidth=0.2)+
  geom_bar(stat="identity",  aes(fill= (coef)), 
           color="black", linewidth=0.4)+
  scale_fill_gradient2(low="blue", mid="white", high="red")+
  theme_classic(base_size = 16)+theme(legend.position="none")+
  facet_grid(direction~., scales="free_y", space="free")+
  labs(x="Abs. LASSO Coefficient", y=NULL)
lasso_coefs_plot

# calculate score by log-ratio
lasso_ratio = do.call(rbind, lapply(1:nrow(unified.df.lasso.data), function(s){
  print(s)
  
  pred = predict(cv_lasso, as.matrix(unified.df.lasso.data[s,!colnames(unified.df.lasso.data) %in% c("timepoint")]))[,1]
  
  ratio = log2(rowSums(unified.df.lasso.data[s,subset(lasso_coefs, coef > 0)$feature]+pseudo.lasso)/
    rowSums(unified.df.lasso.data[s,subset(lasso_coefs, coef < 0)$feature]+pseudo.lasso))
  
  ratio.sparse = log2(rowSums(unified.df.lasso.data[s,subset(lasso_coefs_sparse, coef > 0)$feature]+pseudo.lasso)/
                 rowSums(unified.df.lasso.data[s,subset(lasso_coefs_sparse, coef < 0)$feature]+pseudo.lasso))
  
  sample = rownames(unified.df.lasso.data[s,])
  
  data.frame(sample = sample,
             pred = pred,
             ratio = ratio,
             ratio.sparse = ratio.sparse) %>%
    merge(unified.meta, by="sample")

  }))

ggplot(lasso_ratio, 
       aes(x=timepoint, y=(ratio)))+
  geom_line(aes(group=subject), alpha=0.2)+
  geom_boxplot()+
  geom_point()+
  theme_classic()

wilcox.test(subset(lasso_ratio, timepoint == "baseline")$pred,
            subset(lasso_ratio, timepoint == "RS")$pred, paired=T) # p 2.9e-15
wilcox.test(subset(lasso_ratio, timepoint == "baseline")$ratio,
            subset(lasso_ratio, timepoint == "RS")$ratio, paired=T) # p 3.5e-16
wilcox.test(subset(lasso_ratio, timepoint == "baseline")$ratio.sparse,
            subset(lasso_ratio, timepoint == "RS")$ratio.sparse, paired=T) # p 2e-9
# ratio is optimal

# save these features
lasso_coefs$feature
saveRDS(lasso_coefs$feature, "../2025_11_07_meta_rs_features.Rds")

# :: Score change ---------------------------------------------------------

# plot change in Ratio per individual, plot arrows
# order based on starting value

lasso_scores_plot = lasso_ratio %>% dplyr::select(subject, RS_Name, study, group, timepoint, ratio) %>% pivot_wider(
  names_from = timepoint,
  values_from = ratio,
  names_prefix = "timepoint_") %>%
  mutate(delta = timepoint_RS - timepoint_baseline) %>%
  mutate(d.order = abs(delta)*sign(delta)) %>%
  ggplot(aes(x=reorder(subject, d.order), y=subject))+
  geom_hline(yintercept=0, linetype=2, linewidth=0.2)+
  geom_segment(aes(x=reorder(subject, d.order),
                   xend=reorder(subject, d.order),
                   y=timepoint_baseline, yend=timepoint_RS,
                   color = sign(delta)),
               linewidth=0.3,
               arrow = arrow(length=unit(0.10,"cm"), ends="last", type = "closed"))+
  scale_color_gradient2(low="blue", mid="grey", high="red")+
  coord_flip()+
  theme_classic(base_size = 16)+theme(axis.text.y = element_blank(),
                        axis.ticks.y = element_blank(),
                        legend.position="none")+
  labs(x="Subject", y=expression(Log[2]*FC~RS~Response~Ratio))
lasso_scores_plot


lasso_ratio_wilcox = 
  wilcox.test(subset(lasso_ratio, timepoint == "baseline")$ratio,
              subset(lasso_ratio, timepoint == "RS")$ratio, paired=T)
  
lasso_scores_wilcox_plot = lasso_ratio %>% 
  group_by(subject) %>%
  mutate(direction = sign(ratio - ratio[timepoint == "baseline"])) %>%
  group_by(subject) %>%
  mutate(direction = as.factor(ifelse(direction == 0, direction[timepoint == "RS"], direction)))%>%
  mutate(timepoint = gsub("baseline", "Pre-RS", gsub("RS", "Post-RS", timepoint)))%>%
  mutate(timepoint = factor(timepoint, levels=c("Pre-RS", "Post-RS")))%>%
  arrange(rev(direction))%>%
  ggplot(aes(x=timepoint, y=ratio))+
 # geom_boxplot(width=0.2, outlier.shape=NA)+
  geom_line(aes(color = direction, group=subject),
               linewidth=0.3,
               arrow = arrow(length=unit(0.10,"cm"), ends="last", type = "closed"))+
  geom_segment(x=1, xend=2, y = 10, yend=10, linewidth=0.3)+
  scale_y_continuous(limits=c(-12,11))+
  annotate(geom="text", x=1.5, y=10.5, label="*", size=6)+
  scale_color_manual(values=c("red","blue"))+
  theme_classic(base_size = 12)+theme(#axis.text.y = element_blank(),
                                      axis.ticks.y = element_blank(),
                                      legend.position="none")+
  labs(x=NULL, y=NULL) #expression(Log[2]*FC~RS~Response~Ratio))
lasso_scores_wilcox_plot

lasso_coefs_plot+lasso_scores_plot


lasso_ratio_direction = lasso_ratio %>% 
  group_by(subject) %>%
  mutate(direction = sign(ratio - ratio[timepoint == "baseline"])) %>%
  group_by(subject) %>%
  mutate(direction = as.factor(ifelse(direction == 0, direction[timepoint == "RS"], direction)))%>%
  mutate(timepoint = gsub("baseline", "Pre-RS", gsub("RS", "Post-RS", timepoint)))%>%
  mutate(timepoint = factor(timepoint, levels=c("Pre-RS", "Post-RS")))%>%
  arrange(rev(direction)) %>%
  subset(timepoint == "Post-RS") %>%
  dplyr::select(direction)
table(lasso_ratio_direction$direction)[1] / 208


# 2. RF Δ pH -----------------------------------------------------------

# predict delta pH (no scaling)

# :: a) Default Approach ------------------------------------------------------------------

t1 <- Sys.time()
ml.rf.regression.oob.oars.ph <- 
  do.call(rbind, lapply(rs.names, function(x){
    # y = unique(ml.target.data$variable)[1]
    # x = rs.names[6]
    y = "delta_pH"
    ml.target.data.subset = subset(ml.target.data, variable == y & RS_Name == x)
    rownames(ml.target.data.subset) = ml.target.data.subset$HM
    # do not scale delta, we need absolute measurement
    ml.target.data.subset$value = (ml.target.data.subset$value)
    # repeat 15 iterations
    do.call(rbind, lapply(1:15, function(z){
      # z = 1
      
      full.data = merge(data.frame(ml.input.data), ml.target.data.subset[,c("HM", "value")], by="row.names")
      rownames(full.data) <- full.data$HM
      full.data$Row.names <- NULL
      full.data$HM <- NULL
      full.data
      
      ## Main RF
      set.seed(z)
      rf.results = ranger::ranger(value ~., full.data)
      # output predictions
      rf.predictions = as.vector(rf.results$predictions)
      
      ## Shadow RF (because assuming null data would yield 50% accuracy is not empirical)
      shadow.data = full.data
      set.seed(z)
      shadow.data$value <- sample(shadow.data$value, size = nrow(shadow.data))
      set.seed(z)
      rf.shadow.results = ranger::ranger(value ~., shadow.data)
      # output predictions
      rf.shadow.predictions = as.vector(rf.shadow.results$predictions)
      
      ## collate data
      print(paste(x, " ", y, " ", z))
      rf.data = data.frame(
        HM = row.names(full.data),
        true = full.data$value,
        pred = (rf.predictions),
        shadow = (rf.shadow.predictions)) %>%
        mutate(RS_Name = x,
               iter = z,
               target = y)
      rf.data
    }))}))
t2 <- Sys.time()
t2 - t1 # ~1 mins; 15x faster with ranger

# Overall correlation: 0.556
cor.test(ml.rf.regression.oob.oars.ph$pred,
         ml.rf.regression.oob.oars.ph$true, method="spearman")
ggplot(ml.rf.regression.oob.oars.ph %>%
         mutate(RS_Name = factor(RS_Name, levels=rs.names)),
       aes(x=true, y=pred))+
  geom_point(shape=21, aes(fill=RS_Name))+
  geom_smooth(method="lm", color="black")+
  geom_hline(yintercept=-1.27, linetype=2)+
  geom_vline(xintercept=-1.27, linetype=2)+
  ggpubr::stat_cor(method="spearman")+
  scale_fill_manual(values=labelcolors$cols)+
  theme_classic()+
  labs(x="Measured Δ pH",
       y="Predicted Δ pH")

# RS-specific correlations
do.call(rbind, lapply(rs.names, function(x){data.frame(RS_Name = x, cor = cor.test(subset(ml.rf.regression.oob.oars.ph, RS_Name == x)$pred,
                                                            subset(ml.rf.regression.oob.oars.ph, RS_Name == x)$true, method="spearman")$estimate)}))
ggplot(ml.rf.regression.oob.oars.ph %>%
         #subset(RS_Name %in% c("LetsDoOrganic", "Versafibe1490")) %>%
         mutate(RS_Name = factor(RS_Name, levels=rs.names)),
       aes(x=true, y=pred))+
  geom_point(shape=21, aes(fill=RS_Name))+
  geom_smooth(method="lm", color="black")+
  geom_hline(yintercept=-1.27, linetype=2)+
  geom_vline(xintercept=-1.27, linetype=2)+
  ggpubr::stat_cor(method="spearman")+
  scale_fill_manual(values=labelcolors$cols)+
  theme_classic()+
  facet_wrap(~RS_Name)+
  labs(x="Measured Δ pH",
       y="Predicted Δ pH")


# calculate AUC
ml.rf.regression.oob.oars.ph.auc = ml.rf.regression.oob.oars.ph %>%
  #subset(RS_Name %in% c("LetsDoOrganic", "Versafibe1490")) %>%
  mutate(true.response = ifelse(true < -1.27, "strong", "weak")) %>%
  group_by(RS_Name, iter) %>%
  mutate(cor = cor.test(true, pred)$estimate,
         cor.shadow = cor.test(true, shadow)$estimate) %>%
  mutate(auc = pROC::auc(as.factor(true.response), pred, 
                         direction = ">", levels=c("weak", "strong")) %>% mean()) %>%
  mutate(auc.shadow = pROC::auc(as.factor(true.response), shadow, 
                         direction = ">", levels=c("weak", "strong")) %>% mean()) %>%
  dplyr::select(RS_Name, iter, auc, auc.shadow, cor, cor.shadow) %>% distinct() %>%
  mutate(mean.auc = mean(auc),
         mean.auc.shadow = mean(auc.shadow),
         mean.cor = mean(cor),
         mean.cor.shadow = mean(cor.shadow)) %>%
  dplyr::select(RS_Name, mean.auc, mean.auc.shadow, mean.cor, mean.cor.shadow) %>% distinct()  %>%
  mutate(RS_Name = factor(RS_Name, levels=rs.names))
  
ml.rf.regression.oob.oars.ph.auc.p = do.call(rbind, lapply(rs.names, function(rs){
  data.subset = subset(ml.rf.regression.oob.oars.ph.auc, RS_Name == rs)
  data.frame(RS_Name = rs,
             pval = wilcox.test(data.subset$mean.auc,
                                data.subset$mean.auc.shadow)$p.value)
})) %>%
  mutate(padj = p.adjust(pval, method = "bonferroni")) %>%
  mutate(sig = ifelse(padj < 0.05, "*", "")) %>%
  mutate(RS_Name = factor(RS_Name, levels=rs.names))

ml.rf.regression.oob.oars.ph.auc.mean = ml.rf.regression.oob.oars.ph.auc %>%
  mutate(RS_Name = factor(RS_Name, levels=rs.names))%>%
  group_by(RS_Name) %>%
  mutate(mean.auc = mean(mean.auc),
         mean.shadow = mean(mean.auc.shadow)) %>%
  dplyr::select(RS_Name, mean.auc, mean.shadow) %>% distinct() 

ml.rf.regression.oob.oars.ph.auc.mean = merge(ml.rf.regression.oob.oars.ph.auc.mean,
                                              ml.rf.regression.oob.oars.ph.auc.p, by="RS_Name") %>%
  mutate(RS_Name = factor(RS_Name, levels=rs.names))


ml.rf.regression.oob.oars.ph.auc.plot <- ggplot()+
  # plot shadow
  geom_violin(data=ml.rf.regression.oob.oars.ph.auc, aes(x=RS_Name, y=mean.auc.shadow, fill=RS_Name),fill="grey", alpha=0.4, color=NA)+
  geom_point(data=ml.rf.regression.oob.oars.ph.auc.mean, aes(x=RS_Name, y=mean.shadow), color="grey", size=3, alpha=0.8)+
  # plot real
  geom_violin(data=ml.rf.regression.oob.oars.ph.auc, aes(x=RS_Name, y=mean.auc, fill=RS_Name), color="black")+
  geom_point(data=ml.rf.regression.oob.oars.ph.auc.mean, aes(x=RS_Name, y=mean.auc), shape=21, fill="white", size=2.5)+
  # plot median value
  geom_text(data=ml.rf.regression.oob.oars.ph.auc.mean, aes(x=RS_Name, y=(mean.auc+0.05), label=round(mean.auc, digits=2)), nudge_y=0.04, size=3)+
  # add sig vs shadow
  geom_text(data=ml.rf.regression.oob.oars.ph.auc.mean, aes(x=RS_Name, y=Inf, label=sig), vjust=1.5, size=6)+
  geom_hline(yintercept=0.5, linetype=2, alpha=0.5)+
  scale_y_continuous(breaks=seq(from=0.1, to=0.9, by=0.1), limits=c(0.1,0.9))+
  scale_fill_manual(values=labelcolors$cols[c(1:9)])+
  scale_color_manual(values=labelcolors$cols[c(1:9)])+
  theme_classic()+theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1),
                        legend.position="none",
                        panel.grid.major.y = element_line(color="grey", linewidth=0.2),
                        legend.title = element_text(hjust=0.5),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  facet_wrap(~"Random Forest Out-of-bag AUC")+
  labs(x=NULL, y="AUC", title = "ML")
ml.rf.regression.oob.oars.ph.auc.plot




# :: b) Feature Selection ----------------------------------------------

# Apply Feature Selection from Meta-Analysis (includes B. adol and R. bro)

# n features
n.features = 15

t1 <- Sys.time()
ml.rf.regression.oob.oars.ph.fs <- 
  do.call(rbind, lapply(rs.names, function(x){
    # y = unique(ml.target.data$variable)[1]
    # x = rs.names[6]
    y = "delta_pH"
    ml.target.data.subset = subset(ml.target.data, variable == y & RS_Name == x)
    rownames(ml.target.data.subset) = ml.target.data.subset$HM
    # do not scale delta, we need absolute measurement
    ml.target.data.subset$value = (ml.target.data.subset$value)
    # repeat 15 iterations
    do.call(rbind, lapply(1:15, function(z){
      # z = 1
      
      ml.input.data = ml.input.data[,colnames(ml.input.data) %in% rev(ml_rs_features)[c(1:n.features)]]
      
      full.data = merge(data.frame(ml.input.data), ml.target.data.subset[,c("HM", "value")], by="row.names")
      rownames(full.data) <- full.data$HM
      full.data$Row.names <- NULL
      full.data$HM <- NULL
      full.data
      
      ## Main RF
      set.seed(z)
      rf.results = ranger::ranger(value ~., full.data)
      # output predictions
      rf.predictions = as.vector(rf.results$predictions)
      
      ## Shadow RF (because assuming null data would yield 50% accuracy is not empirical)
      shadow.data = full.data
      set.seed(z)
      shadow.data$value <- sample(shadow.data$value, size = nrow(shadow.data))
      set.seed(z)
      rf.shadow.results = ranger::ranger(value ~., shadow.data)
      # output predictions
      rf.shadow.predictions = as.vector(rf.shadow.results$predictions)
      
      ## collate data
      print(paste(x, " ", y, " ", z))
      rf.data = data.frame(
        HM = row.names(full.data),
        true = full.data$value,
        pred = (rf.predictions),
        shadow = (rf.shadow.predictions)) %>%
        mutate(RS_Name = x,
               iter = z,
               target = y)
      rf.data
    }))}))
t2 <- Sys.time()
t2 - t1 # ~1 mins; 15x faster with ranger

# Overall correlation: 0.56
cor.test(ml.rf.regression.oob.oars.ph.fs$pred,
         ml.rf.regression.oob.oars.ph.fs$true, method="spearman")
ggplot(ml.rf.regression.oob.oars.ph.fs %>%
         mutate(RS_Name = factor(RS_Name, levels=rs.names)),
       aes(x=true, y=pred))+
  geom_point(shape=21, aes(fill=RS_Name))+
  geom_smooth(method="lm", color="black")+
  geom_hline(yintercept=-1.27, linetype=2)+
  geom_vline(xintercept=-1.27, linetype=2)+
  ggpubr::stat_cor(method="spearman")+
  scale_fill_manual(values=labelcolors$cols)+
  theme_classic()+
  labs(x="Measured Δ pH",
       y="Predicted Δ pH")

# RS-specific correlations
do.call(rbind, lapply(rs.names, function(x){data.frame(RS_Name = x, cor = cor.test(subset(ml.rf.regression.oob.oars.ph.fs, RS_Name == x)$pred,
                                                                                   subset(ml.rf.regression.oob.oars.ph.fs, RS_Name == x)$true, method="spearman")$estimate)}))
ggplot(ml.rf.regression.oob.oars.ph.fs %>%
         #subset(RS_Name %in% c("LetsDoOrganic", "Versafibe1490")) %>%
         mutate(RS_Name = factor(RS_Name, levels=rs.names)),
       aes(x=true, y=pred))+
  geom_point(shape=21, aes(fill=RS_Name))+
  geom_smooth(method="lm", color="black")+
  geom_hline(yintercept=-1.27, linetype=2)+
  geom_vline(xintercept=-1.27, linetype=2)+
  ggpubr::stat_cor(method="spearman")+
  scale_fill_manual(values=labelcolors$cols)+
  theme_classic()+
  facet_wrap(~RS_Name)+
  labs(x="Measured Δ pH",
       y="Predicted Δ pH")


# calculate AUC
ml.rf.regression.oob.oars.ph.fs.auc = ml.rf.regression.oob.oars.ph.fs %>%
  #subset(RS_Name %in% c("LetsDoOrganic", "Versafibe1490")) %>%
  mutate(true.response = ifelse(true < -1.27, "strong", "weak")) %>%
  group_by(RS_Name, iter) %>%
  mutate(auc = pROC::auc(as.factor(true.response), pred, 
                         direction = ">", levels=c("weak", "strong")) %>% mean()) %>%
  mutate(auc.shadow = pROC::auc(as.factor(true.response), shadow, 
                                direction = ">", levels=c("weak", "strong")) %>% mean()) %>%
  dplyr::select(RS_Name, iter, auc, auc.shadow) %>% distinct() %>%
  mutate(mean.auc = mean(auc),
         mean.auc.shadow = mean(auc.shadow)) %>%
  dplyr::select(RS_Name, mean.auc, mean.auc.shadow) %>% distinct()  %>%
  mutate(RS_Name = factor(RS_Name, levels=rs.names))

ml.rf.regression.oob.oars.ph.fs.auc.p = do.call(rbind, lapply(rs.names, function(rs){
  data.subset = subset(ml.rf.regression.oob.oars.ph.fs.auc, RS_Name == rs)
  data.frame(RS_Name = rs,
             pval = wilcox.test(data.subset$mean.auc,
                                data.subset$mean.auc.shadow)$p.value)
})) %>%
  mutate(padj = p.adjust(pval, method = "bonferroni")) %>%
  mutate(sig = ifelse(padj < 0.05, "*", ""))

ml.rf.regression.oob.oars.ph.fs.auc.mean = ml.rf.regression.oob.oars.ph.fs.auc %>%
  group_by(RS_Name) %>%
  mutate(mean.auc = mean(mean.auc),
         mean.shadow = mean(mean.auc.shadow)) %>%
  dplyr::select(RS_Name, mean.auc, mean.shadow) %>% distinct()

ml.rf.regression.oob.oars.ph.fs.auc.mean = merge(ml.rf.regression.oob.oars.ph.fs.auc.mean,
                                              ml.rf.regression.oob.oars.ph.fs.auc.p, by="RS_Name") %>%
  mutate(RS_Name = factor(RS_Name, levels=rs.names))


ml.rf.regression.oob.oars.ph.fs.auc.plot <- ggplot()+
  # plot shadow
  geom_violin(data=ml.rf.regression.oob.oars.ph.fs.auc, aes(x=RS_Name, y=mean.auc.shadow, fill=RS_Name),fill="grey", alpha=0.4, color=NA)+
  geom_point(data=ml.rf.regression.oob.oars.ph.fs.auc.mean, aes(x=RS_Name, y=mean.shadow), color="grey", size=3, alpha=0.8)+
  # plot real
  geom_violin(data=ml.rf.regression.oob.oars.ph.fs.auc, aes(x=RS_Name, y=mean.auc, fill=RS_Name), color="black")+
  geom_point(data=ml.rf.regression.oob.oars.ph.fs.auc.mean, aes(x=RS_Name, y=mean.auc), shape=21, fill="white", size=2.5)+
  # plot median value
  geom_text(data=ml.rf.regression.oob.oars.ph.fs.auc.mean, aes(x=RS_Name, y=(mean.auc+0.05), label=round(mean.auc, digits=2)), nudge_y=0.04, size=3)+
  # add sig vs shadow
  geom_text(data=ml.rf.regression.oob.oars.ph.fs.auc.mean, aes(x=RS_Name, y=Inf, label=sig), vjust=1.5, size=6)+
  geom_hline(yintercept=0.5, linetype=2, alpha=0.5)+
  scale_y_continuous(breaks=seq(from=0.1, to=0.9, by=0.1), limits=c(0.1,0.9))+
  scale_fill_manual(values=labelcolors$cols[c(1:9)])+
  scale_color_manual(values=labelcolors$cols[c(1:9)])+
  scale_alpha_manual(values=c(0.3, 0.8))+
  theme_classic()+theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1),
                        legend.position="none",
                        panel.grid.major.y = element_line(color="grey", linewidth=0.2),
                        legend.title = element_text(hjust=0.5),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  facet_wrap(~"Random Forest Out-of-bag AUC")+
  labs(x=NULL, y="AUC", title = "ML FS")
ml.rf.regression.oob.oars.ph.fs.auc.plot


ml.rf.regression.oob.oars.ph.auc.plot + 
  ml.rf.regression.oob.oars.ph.fs.auc.plot


# :: c) PBS ** ------------------------------------------------------------------

# No feature selection 

t1 <- Sys.time()
ml.rf.regression.oob.oars.ph.pbs <- 
  do.call(rbind, lapply(rs.names.pbs, function(x){
    # y = unique(ml.target.data$variable)[1]
    # x = rs.names[6]
    y = "pH"
    ml.target.data.subset = subset(ph.data, RS_Name == x)
    rownames(ml.target.data.subset) = ml.target.data.subset$HM
    # do not scale delta, we need absolute measurement
    ml.target.data.subset$value = (ml.target.data.subset$med.ph)
    # repeat 15 iterations
    do.call(rbind, lapply(1:15, function(z){
      # z = 1
      
      full.data = merge(data.frame(ml.input.data), ml.target.data.subset[,c("HM", "value")], by="row.names")
      rownames(full.data) <- full.data$HM
      full.data$Row.names <- NULL
      full.data$HM <- NULL
      full.data
      
      ## Main RF
      set.seed(z)
      rf.results = ranger::ranger(value ~., full.data)
      # output predictions
      rf.predictions = as.vector(rf.results$predictions)
      
      ## Shadow RF (because assuming null data would yield 50% accuracy is not empirical)
      shadow.data = full.data
      set.seed(z)
      shadow.data$value <- sample(shadow.data$value, size = nrow(shadow.data))
      set.seed(z)
      rf.shadow.results = ranger::ranger(value ~., shadow.data)
      # output predictions
      rf.shadow.predictions = as.vector(rf.shadow.results$predictions)
      
      ## collate data
      print(paste(x, " ", y, " ", z))
      rf.data = data.frame(
        HM = row.names(full.data),
        true = full.data$value,
        pred = (rf.predictions),
        shadow = (rf.shadow.predictions)) %>%
        mutate(RS_Name = x,
               iter = z,
               target = y)
      rf.data
    }))}))
t2 <- Sys.time()
t2 - t1 # ~1 mins; 15x faster with ranger

# first, check correlation of raw pH predictions (not delta)
ml.rf.regression.oob.oars.ph.w.pbs = ml.rf.regression.oob.oars.ph.pbs
ggplot(ml.rf.regression.oob.oars.ph.w.pbs %>%
         mutate(RS_Name = factor(RS_Name, levels=rs.names.pbs)),
       aes(x=true, y=pred))+
  geom_point(shape=21, aes(fill=RS_Name))+
  geom_smooth(method="lm", color="black")+
  #geom_hline(yintercept=-1.27, linetype=2)+
  #geom_vline(xintercept=-1.27, linetype=2)+
  ylim(5.2,7.8)+
  ggpubr::stat_cor(method="spearman")+
  scale_fill_manual(values=c("grey", labelcolors$cols[c(1:9)]))+
  theme_classic()+
  facet_wrap(~RS_Name, scales="free")+
  labs(title = "PBS",
       x="Measured pH",
       y="Predicted pH")

# now predict change in pH using predictions
ml.rf.regression.oob.oars.ph.pbs = ml.rf.regression.oob.oars.ph.pbs %>%
  group_by(HM, iter) %>%
  mutate(true.delta = true - true[RS_Name == "PBS"],
         pred.delta = pred - pred[RS_Name == "PBS"],
         shadow.delta = shadow - shadow[RS_Name == "PBS"]) %>%
  subset(RS_Name != "PBS")

# Overall correlation: 0.589
cor.test(ml.rf.regression.oob.oars.ph.pbs$pred.delta,
         ml.rf.regression.oob.oars.ph.pbs$true.delta, method="spearman")
ggplot(ml.rf.regression.oob.oars.ph.pbs %>%
         mutate(RS_Name = factor(RS_Name, levels=rs.names)),
       aes(x=true.delta, y=pred.delta))+
  geom_point(shape=21, aes(fill=RS_Name))+
  geom_smooth(method="lm", color="black")+
  geom_hline(yintercept=-1.27, linetype=2)+
  geom_vline(xintercept=-1.27, linetype=2)+
  ggpubr::stat_cor(method="spearman")+
  scale_fill_manual(values=labelcolors$cols)+
  theme_classic()+
  labs(title = "PBS",
       x="Measured Δ pH",
       y="Predicted Δ pH")

# RS-specific correlations
do.call(rbind, lapply(rs.names, function(x){data.frame(RS_Name = x, cor = cor.test(subset(ml.rf.regression.oob.oars.ph.pbs, RS_Name == x)$pred.delta,
                                                                                   subset(ml.rf.regression.oob.oars.ph.pbs, RS_Name == x)$true.delta, method="spearman")$estimate)}))

# make neater
ml.rf.regression.oob.oars.ph.pbs.mean = ml.rf.regression.oob.oars.ph.pbs %>%
  group_by(HM, RS_Name) %>%
  mutate(mean.pred = mean(pred.delta),
         max.pred = max(pred.delta),
         min.pred = min(pred.delta)) %>%
  dplyr::select(HM, RS_Name, mean.pred, max.pred, min.pred, true.delta) %>% distinct()
   

# :: :: PLOT: Scatterplots ------------------------------------------------

ml.rf.regression.oob.oars.ph.pbs.plot.neat = ggplot(ml.rf.regression.oob.oars.ph.pbs %>%
         subset(RS_Name != "PBS") %>%
         mutate(RS_Name = factor(RS_Name, levels=rs.names)),
       aes(x=true.delta, y=pred.delta))+
  geom_segment(data= ml.rf.regression.oob.oars.ph.pbs.mean%>%
                 subset(RS_Name != "PBS") %>%
                 mutate(RS_Name = factor(RS_Name, levels=rs.names)),
               aes(x = true.delta, xend = true.delta, y=min.pred, yend=max.pred, color=RS_Name))+
  geom_point(data= ml.rf.regression.oob.oars.ph.pbs.mean%>%
                 subset(RS_Name != "PBS") %>%
                 mutate(RS_Name = factor(RS_Name, levels=rs.names)),
               shape=21, aes(x = true.delta, y=mean.pred, fill=RS_Name))+
  geom_smooth(method="lm", color="black")+
  geom_hline(yintercept=-1.27, linetype=2, linewidth=0.2)+
  geom_vline(xintercept=-1.27, linetype=2, linewidth=0.2)+
  scale_y_continuous(limits=c(-2, 0.25))+
  ggpubr::stat_cor(method="spearman",aes(label = ..r.label..), size=3)+
  scale_fill_manual(values=labelcolors$cols)+
  scale_color_manual(values=labelcolors$cols)+
  theme_classic()+theme(legend.position="none")+
  facet_wrap(~RS_Name)+
  labs(x="Measured Δ pH",
       y="Predicted Δ pH")
ml.rf.regression.oob.oars.ph.pbs.plot.neat

ggplot(ml.rf.regression.oob.oars.ph.pbs %>%
         group_by(HM, RS_Name) %>%
         mutate(mean.pred = mean(pred.delta))%>%
         group_by(HM) %>%
         slice_min(mean.pred) %>%
         dplyr::select(HM, mean.pred, RS_Name, true.delta) %>% distinct(),
       aes(x=true.delta, y=mean.pred))+
  geom_point(shape=21)+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method="spearman", color="black")+
  theme_classic()

# calculate AUC
ml.rf.regression.oob.oars.ph.pbs.auc = ml.rf.regression.oob.oars.ph.pbs %>%
  subset(RS_Name != "PBS") %>%
  #subset(RS_Name %in% c("LetsDoOrganic", "Versafibe1490")) %>%
  mutate(true.response = ifelse(true.delta < -1.27, "strong", "weak")) %>%
  group_by(RS_Name, iter) %>%
  mutate(cor = cor.test(true.delta, pred.delta)$estimate,
         cor.shadow = cor.test(true.delta, shadow.delta)$estimate) %>%
  mutate(auc = pROC::auc(as.factor(true.response), pred.delta, 
                         direction = ">", levels=c("weak", "strong")) %>% mean()) %>%
  mutate(auc.shadow = pROC::auc(as.factor(true.response), shadow.delta, 
                                direction = ">", levels=c("weak", "strong")) %>% mean()) %>%
  dplyr::select(RS_Name, iter, auc, auc.shadow, cor, cor.shadow) %>% distinct() %>%
  mutate(mean.auc = mean(auc),
         mean.auc.shadow = mean(auc.shadow),
         mean.cor = mean(cor),
         mean.cor.shadow = mean(cor.shadow)) %>%
  dplyr::select(RS_Name, mean.auc, mean.auc.shadow, mean.cor, mean.cor.shadow) %>% distinct()  %>%
  mutate(RS_Name = factor(RS_Name, levels=rs.names))


ml.rf.regression.oob.oars.ph.pbs.auc.p = do.call(rbind, lapply(rs.names, function(rs){
  data.subset = subset(ml.rf.regression.oob.oars.ph.pbs.auc, RS_Name == rs)
  data.frame(RS_Name = rs,
             pval = wilcox.test(data.subset$mean.auc,
                                data.subset$mean.auc.shadow)$p.value)
})) %>%
  mutate(padj = p.adjust(pval, method = "bonferroni")) %>%
  mutate(sig = ifelse(padj < 0.05, "*", "")) %>%
  mutate(RS_Name = factor(RS_Name, levels=rs.names))

ml.rf.regression.oob.oars.ph.pbs.auc.mean = ml.rf.regression.oob.oars.ph.pbs.auc %>%
  mutate(RS_Name = factor(RS_Name, levels=rs.names))%>%
  group_by(RS_Name) %>%
  mutate(mean.auc = mean(mean.auc),
         mean.shadow = mean(mean.auc.shadow)) %>%
  dplyr::select(RS_Name, mean.auc, mean.shadow) %>% distinct() 

ml.rf.regression.oob.oars.ph.pbs.auc.mean = merge(ml.rf.regression.oob.oars.ph.pbs.auc.mean,
                                                  ml.rf.regression.oob.oars.ph.pbs.auc.p, by="RS_Name") %>%
  mutate(RS_Name = factor(RS_Name, levels=rs.names))


ml.rf.regression.oob.oars.ph.pbs.auc.plot <- ggplot()+
  # plot shadow
  geom_violin(data=ml.rf.regression.oob.oars.ph.pbs.auc, aes(x=RS_Name, y=mean.auc.shadow, fill=RS_Name),fill="grey", alpha=0.4, color=NA)+
  geom_point(data=ml.rf.regression.oob.oars.ph.pbs.auc.mean, aes(x=RS_Name, y=mean.shadow), color="grey", size=3, alpha=0.8)+
  # plot real
  geom_violin(data=ml.rf.regression.oob.oars.ph.pbs.auc, aes(x=RS_Name, y=mean.auc, fill=RS_Name), color="black")+
  geom_point(data=ml.rf.regression.oob.oars.ph.pbs.auc.mean, aes(x=RS_Name, y=mean.auc), shape=21, fill="white", size=2.5)+
  # plot median value
  geom_text(data=ml.rf.regression.oob.oars.ph.pbs.auc.mean, aes(x=RS_Name, y=(mean.auc+0.05), label=round(mean.auc, digits=2)), nudge_y=0.04, size=3)+
  # add sig vs shadow
  geom_text(data=ml.rf.regression.oob.oars.ph.pbs.auc.mean, aes(x=RS_Name, y=Inf, label=sig), vjust=1.5, size=6)+
  geom_hline(yintercept=0.5, linetype=2, alpha=0.5)+
  scale_y_continuous(breaks=seq(from=0.1, to=0.9, by=0.1), limits=c(0.1,0.9))+
  scale_fill_manual(values=labelcolors$cols[c(1:9)])+
  scale_color_manual(values=labelcolors$cols[c(1:9)])+
  theme_classic()+theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1),
                        legend.position="none",
                        panel.grid.major.y = element_line(color="grey", linewidth=0.2),
                        legend.title = element_text(hjust=0.5),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  facet_wrap(~"Random Forest Out-of-bag AUC")+
  labs(x=NULL, y="AUC", title = "PBS")
ml.rf.regression.oob.oars.ph.pbs.auc.plot

# :: d) FS + PBS ** -----------------------------------------------------

# Apply feature selection to all except PBS

t1 <- Sys.time()
ml.rf.regression.oob.oars.ph.pbs.fs <- 
  do.call(rbind, lapply(rs.names.pbs, function(x){
    # y = unique(ml.target.data$variable)[1]
    # x = rs.names[6]
    y = "pH"
    ml.target.data.subset = subset(ph.data, RS_Name == x)
    rownames(ml.target.data.subset) = ml.target.data.subset$HM
    # do not scale delta, we need absolute measurement
    ml.target.data.subset$value = (ml.target.data.subset$med.ph)
    # repeat 15 iterations
    do.call(rbind, lapply(1:15, function(z){
      # z = 1
      
      if(x != "PBS"){
        ml.input.data = ml.input.data[,colnames(ml.input.data) %in% lasso_coefs$feature]
        #ml.input.data = ml.input.data[,colnames(ml.input.data) %in% ml_rs_features_top]
      }
      full.data = merge(data.frame(ml.input.data), ml.target.data.subset[,c("HM", "value")], by="row.names")
      rownames(full.data) <- full.data$HM
      full.data$Row.names <- NULL
      full.data$HM <- NULL
      full.data
      
      ## Main RF
      set.seed(z)
      rf.results = ranger::ranger(value ~., full.data)
      # output predictions
      rf.predictions = as.vector(rf.results$predictions)
      
      ## Shadow RF (because assuming null data would yield 50% accuracy is not empirical)
      shadow.data = full.data
      set.seed(z)
      shadow.data$value <- sample(shadow.data$value, size = nrow(shadow.data))
      set.seed(z)
      rf.shadow.results = ranger::ranger(value ~., shadow.data)
      # output predictions
      rf.shadow.predictions = as.vector(rf.shadow.results$predictions)
      
      ## collate data
      print(paste(x, " ", y, " ", z))
      rf.data = data.frame(
        HM = row.names(full.data),
        true = full.data$value,
        pred = (rf.predictions),
        shadow = (rf.shadow.predictions)) %>%
        mutate(RS_Name = x,
               iter = z,
               target = y)
      rf.data
    }))}))
t2 <- Sys.time()
t2 - t1 # ~1 mins; 15x faster with ranger


# now predict change in pH using predictions
ml.rf.regression.oob.oars.ph.pbs.fs = ml.rf.regression.oob.oars.ph.pbs.fs %>%
  group_by(HM, iter) %>%
  mutate(true.delta = true - true[RS_Name == "PBS"],
         pred.delta = pred - pred[RS_Name == "PBS"],
         shadow.delta = shadow - shadow[RS_Name == "PBS"]) %>%
  subset(RS_Name != "PBS")


# Overall correlation: 0.52
cor.test(ml.rf.regression.oob.oars.ph.pbs.fs$pred.delta,
         ml.rf.regression.oob.oars.ph.pbs.fs$true.delta, method="spearman")
ggplot(ml.rf.regression.oob.oars.ph.pbs.fs %>%
         mutate(RS_Name = factor(RS_Name, levels=rs.names)),
       aes(x=true.delta, y=pred.delta))+
  geom_point(shape=21, aes(fill=RS_Name))+
  geom_smooth(method="lm", color="black")+
  geom_hline(yintercept=-1.27, linetype=2)+
  geom_vline(xintercept=-1.27, linetype=2)+
  ggpubr::stat_cor(method="spearman")+
  scale_fill_manual(values=labelcolors$cols)+
  theme_classic()+
  labs(title = "FS + PBS",
       x="Measured Δ pH",
       y="Predicted Δ pH")

# RS-specific correlations
do.call(rbind, lapply(rs.names, function(x){data.frame(RS_Name = x, cor = cor.test(subset(ml.rf.regression.oob.oars.ph.pbs.fs, RS_Name == x)$pred.delta,
                                                                                   subset(ml.rf.regression.oob.oars.ph.pbs.fs, RS_Name == x)$true.delta, method="spearman")$estimate)}))

# make neater
ml.rf.regression.oob.oars.ph.pbs.fs.mean = ml.rf.regression.oob.oars.ph.pbs.fs %>%
  group_by(HM, RS_Name) %>%
  mutate(mean.pred = mean(pred.delta),
         max.pred = max(pred.delta),
         min.pred = min(pred.delta)) %>%
  dplyr::select(HM, RS_Name, mean.pred, max.pred, min.pred, true.delta) %>% distinct()


ml.rf.regression.oob.oars.ph.pbs.fs.cor.plot = ggplot(ml.rf.regression.oob.oars.ph.pbs.fs %>%
         subset(RS_Name != "PBS") %>%
         mutate(RS_Name = factor(RS_Name, levels=rs.names)),
       aes(x=true.delta, y=pred.delta))+
  geom_segment(data= ml.rf.regression.oob.oars.ph.pbs.fs.mean%>%
                 subset(RS_Name != "PBS") %>%
                 mutate(RS_Name = factor(RS_Name, levels=rs.names)),
               aes(x = true.delta, xend = true.delta, y=min.pred, yend=max.pred, color=RS_Name))+
  geom_point(data= ml.rf.regression.oob.oars.ph.pbs.fs.mean%>%
               subset(RS_Name != "PBS") %>%
               mutate(RS_Name = factor(RS_Name, levels=rs.names)),
             shape=21, aes(x = true.delta, y=mean.pred, fill=RS_Name))+
  geom_smooth(method="lm", color="black")+
  scale_y_continuous(limits=c(-2, 0.25))+
  geom_hline(yintercept=-1.27, linetype=2, linewidth=0.2)+
  geom_vline(xintercept=-1.27, linetype=2, linewidth=0.2)+
  ggpubr::stat_cor(method="spearman",aes(label = ..r.label..))+
  scale_fill_manual(values=labelcolors$cols)+
  theme_classic()+theme(legend.position="none")+
  facet_wrap(~RS_Name)+
  labs(x="Measured Δ pH",
       y="Predicted Δ pH")
ml.rf.regression.oob.oars.ph.pbs.fs.cor.plot

# calculate AUC
ml.rf.regression.oob.oars.ph.pbs.fs.auc = ml.rf.regression.oob.oars.ph.pbs.fs %>%
  subset(RS_Name != "PBS") %>%
  #subset(RS_Name %in% c("LetsDoOrganic", "Versafibe1490")) %>%
  mutate(true.response = ifelse(true.delta < -1.27, "strong", "weak")) %>%
  group_by(RS_Name, iter) %>%
  mutate(cor = cor.test(true.delta, pred.delta)$estimate,
         cor.shadow = cor.test(true.delta, shadow.delta)$estimate) %>%
  mutate(auc = pROC::auc(as.factor(true.response), pred.delta, 
                         direction = ">", levels=c("weak", "strong")) %>% mean()) %>%
  mutate(auc.shadow = pROC::auc(as.factor(true.response), shadow, 
                                direction = ">", levels=c("weak", "strong")) %>% mean()) %>%
  dplyr::select(RS_Name, iter, auc, auc.shadow, cor, cor.shadow) %>% distinct() %>%
  mutate(mean.auc = (auc),
         mean.auc.shadow = (auc.shadow),
         mean.cor = (cor),
         mean.cor.shadow = (cor.shadow)) %>%
  dplyr::select(RS_Name, mean.auc, mean.auc.shadow, mean.cor, mean.cor.shadow) %>% distinct()  %>%
  mutate(RS_Name = factor(RS_Name, levels=rs.names))

ml.rf.regression.oob.oars.ph.pbs.fs.auc.p = do.call(rbind, lapply(rs.names, function(rs){
  data.subset = subset(ml.rf.regression.oob.oars.ph.pbs.fs.auc, RS_Name == rs)
  data.frame(RS_Name = rs,
             cor.pval = wilcox.test(data.subset$mean.cor,
                                    data.subset$mean.cor.shadow)$p.value,
             auc.pval = wilcox.test(data.subset$mean.auc,
                                data.subset$mean.auc.shadow)$p.value)
})) %>%
  mutate(auc.padj = p.adjust(auc.pval, method = "bonferroni")) %>%
  mutate(cor.padj = p.adjust(cor.pval, method = "bonferroni")) %>%
  mutate(auc.sig = ifelse(auc.padj < 0.05, "*", "")) %>%
  mutate(cor.sig = ifelse(cor.padj < 0.05, "*", "")) %>%
  
  mutate(RS_Name = factor(RS_Name, levels=rs.names))

ml.rf.regression.oob.oars.ph.pbs.fs.auc.mean = ml.rf.regression.oob.oars.ph.pbs.fs.auc %>%
  mutate(RS_Name = factor(RS_Name, levels=rs.names))%>%
  group_by(RS_Name) %>%
  mutate(mean.auc = mean(mean.auc),
         mean.auc.shadow = mean(mean.auc.shadow)) %>%
  mutate(mean.cor = mean(mean.cor),
         mean.cor.shadow = mean(mean.cor.shadow)) %>%
  dplyr::select(RS_Name, mean.auc, mean.cor, mean.cor.shadow, mean.auc.shadow) %>% distinct() 

ml.rf.regression.oob.oars.ph.pbs.fs.auc.mean = merge(ml.rf.regression.oob.oars.ph.pbs.fs.auc.mean,
                                                  ml.rf.regression.oob.oars.ph.pbs.fs.auc.p, by="RS_Name") %>%
  mutate(RS_Name = factor(RS_Name, levels=rs.names))

# :: :: PLOT: RF FS PBS Performance ------------------------------------------------------

ml.rf.regression.oob.oars.ph.pbs.fs.cor.plot <- ggplot()+
  # plot shadow
  geom_violin(data=ml.rf.regression.oob.oars.ph.pbs.fs.auc, aes(x=RS_Name, y=mean.cor.shadow, fill=RS_Name),fill="grey", alpha=0.4, color=NA)+
  geom_point(data=ml.rf.regression.oob.oars.ph.pbs.fs.auc.mean, aes(x=RS_Name, y=mean.cor.shadow), color="grey", size=3, alpha=0.8)+
  # plot real
  geom_violin(data=ml.rf.regression.oob.oars.ph.pbs.fs.auc, aes(x=RS_Name, y=mean.cor, fill=RS_Name), color="black")+
  geom_point(data=ml.rf.regression.oob.oars.ph.pbs.fs.auc.mean, aes(x=RS_Name, y=mean.cor, fill=RS_Name), shape=21,  size=2.5)+
  # plot median value
  geom_text(data=ml.rf.regression.oob.oars.ph.pbs.fs.auc.mean, aes(x=RS_Name, y=mean.cor+0.05, label=round(mean.cor, digits=2)), nudge_y=0.04, size=4)+
  # add sig vs shadow
  geom_text(data=ml.rf.regression.oob.oars.ph.pbs.fs.auc.mean, aes(x=RS_Name, y=Inf, label=cor.sig), vjust=1.5, size=7)+
  geom_hline(yintercept=0, linetype=2, alpha=0.5)+
  scale_y_continuous(breaks=seq(from=-0.15, to=0.65, by=0.1), limits=c(-0.15,0.65))+
  scale_fill_manual(values=labelcolors$cols[c(1:9)])+
  scale_color_manual(values=labelcolors$cols[c(1:9)])+
  theme_classic(base_size = 16)+theme(#axis.text.x=element_text(angle=45,vjust=1, hjust=1),
                        legend.position="none",
                        panel.grid.major.y = element_line(color="grey", linewidth=0.2),
                        legend.title = element_text(hjust=0.5),
                        strip.text = element_text(size=16),
                        axis.text.x = element_blank(),
                        strip.background = element_rect(
                          color="black"))+
  facet_wrap(~"RF OOB with Feature Selection")+
  labs(x=NULL, y="Spearman ρ")
ml.rf.regression.oob.oars.ph.pbs.fs.cor.plot

ml.rf.regression.oob.oars.ph.pbs.fs.auc.plot <- ggplot()+
  # plot shadow
  geom_violin(data=ml.rf.regression.oob.oars.ph.pbs.fs.auc, aes(x=RS_Name, y=mean.auc.shadow, fill=RS_Name),fill="grey", alpha=0.4, color=NA)+
  geom_point(data=ml.rf.regression.oob.oars.ph.pbs.fs.auc.mean, aes(x=RS_Name, y=mean.auc.shadow), color="grey", size=3, alpha=0.8)+
  # plot real
  geom_violin(data=ml.rf.regression.oob.oars.ph.pbs.fs.auc, aes(x=RS_Name, y=mean.auc, fill=RS_Name), color="black")+
  geom_point(data=ml.rf.regression.oob.oars.ph.pbs.fs.auc.mean, aes(x=RS_Name, y=mean.auc, fill=RS_Name), shape=21, size=2.5)+
  # plot median value
  geom_text(data=ml.rf.regression.oob.oars.ph.pbs.fs.auc.mean, aes(x=RS_Name, y=ifelse(RS_Name != "Versafibe1490", mean.auc+0.05, mean.auc + 0.15), label=round(mean.auc, digits=2)), nudge_y=0.04, size=4)+
  # add sig vs shadow
  geom_text(data=ml.rf.regression.oob.oars.ph.pbs.fs.auc.mean, aes(x=RS_Name, y=Inf, label=auc.sig), vjust=2, size=7)+
  geom_hline(yintercept=0.5, linetype=2, alpha=0.5)+
  scale_y_continuous(breaks=seq(from=0.1, to=1., by=0.1), limits=c(0.1, 1.))+
  scale_fill_manual(values=labelcolors$cols[c(1:9)])+
  scale_color_manual(values=labelcolors$cols[c(1:9)])+
  theme_classic(base_size = 16)+theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1),
                        legend.position="none",
                        panel.grid.major.y = element_line(color="grey", linewidth=0.2),
                        legend.title = element_text(hjust=0.5),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  #facet_wrap(~"RF OOB with Feature Selection")+
  labs(x=NULL, y="AUC")
ml.rf.regression.oob.oars.ph.pbs.fs.auc.plot

ml.rf.regression.oob.oars.ph.pbs.auc.plot+
  ml.rf.regression.oob.oars.ph.pbs.fs.auc.plot

# :: FS + PBS ROC ---------------------------------------------------------

ml.rf.oob.auc.ph.pbs.fs.roc = do.call(rbind, lapply(rs.names, function(rs){
  # take mean prediction
  data.subset = subset(ml.rf.regression.oob.oars.ph.pbs.fs, RS_Name == rs) %>%
    group_by(HM) %>%
    mutate(mean.pred = mean(pred.delta)) %>%
    dplyr::select(HM, RS_Name, true.delta, mean.pred) %>% distinct() %>%
    mutate(true.response = factor(ifelse(true.delta < -1.27, "strong", "weak"), levels=c("strong", "weak")))
  roc.curve = pROC::roc(data.subset$true.response, data.subset$mean.pred, direction = "<") %>% suppressMessages()
  # extract coordinates
  roc.curve <- data.frame(
    TPR=rev(roc.curve$sensitivities), 
    FPR=rev(1-roc.curve$specificities)) %>%
    mutate(RS_Name = rs)
}))
  

# :: :: PLOT: ROC Curve ---------------------------------------------------

ml.rf.oob.auc.ph.pbs.fs.roc.plot = ggplot(ml.rf.oob.auc.ph.pbs.fs.roc %>%
         mutate(RS_Name = factor(RS_Name, levels=rs.names)),
       aes(x=FPR, y=TPR))+
  geom_line(aes(color=RS_Name, linetype=RS_Name), linewidth = 0.75)+
  scale_color_manual(values=labelcolors.rs)+
  scale_linetype_manual(values=c(4,8,7,6,5,1,3,2,9))+
  geom_abline(slope=1, linetype=2, linewidth=0.2)+
  theme_classic(base_size = 16)+
                        theme(strip.text = element_text(size=16),
                        legend.position = c(0.85, 0.35),
                        legend.text = element_text(size=7),
                        legend.title=element_blank(),
                        legend.key.spacing.y = unit(-0.2, "cm"))+
  facet_wrap(~"ROC Curves")+
  labs(x="False Positive Rate", y="True Positive Rate")
ml.rf.oob.auc.ph.pbs.fs.roc.plot

# 3. IMPORTANCE --------------------------------------------------


# :: PBS ------------------------------------------------------------------


# use PBS method without feature selection first

t1 <- Sys.time()
ml.rf.regression.oob.oars.ph.pbs.importances <- 
  do.call(rbind, lapply(rs.names.pbs, function(x){
    # y = unique(ml.target.data$variable)[1]
    # x = rs.names[6]
    y = "pH"
    ml.target.data.subset = subset(ph.data, RS_Name == x)
    rownames(ml.target.data.subset) = ml.target.data.subset$HM
    # do not scale delta, we need absolute measurement
    ml.target.data.subset$value = (ml.target.data.subset$med.ph)
    # repeat 15 iterations
    do.call(rbind, lapply(1:15, function(z){
      # z = 1
      
      full.data = merge(data.frame(ml.input.data), ml.target.data.subset[,c("HM", "value")], by="row.names")
      rownames(full.data) <- full.data$HM
      full.data$Row.names <- NULL
      full.data$HM <- NULL
      full.data
      
      ## Main RF
      set.seed(z)
      rf.results = ranger::ranger(value ~., full.data, importance="permutation")
      # output importances
      rf.predictions = data.frame(imp = rf.results$variable.importance)
      rf.predictions$feature = rownames(rf.predictions)
      # calculate correlation coefficient (and sign)
      rf.feature.cor = reshape2::melt(cor(as.matrix(full.data), method="spearman")) %>%
        subset(Var1 == "value") %>%
        subset(Var2 != "value") %>%
        dplyr::select(Var2, value)
      colnames(rf.feature.cor) = c("feature", "cor")
      rf.predictions = merge(rf.predictions, rf.feature.cor, type="feature")
      
      ## Shadow RF (because assuming null data would yield 50% accuracy is not empirical)
      shadow.data = full.data
      set.seed(z)
      shadow.data$value <- sample(shadow.data$value, size = nrow(shadow.data))
      set.seed(z)
      rf.shadow.results = ranger::ranger(value ~., shadow.data, importance="permutation")
      # output importances
      rf.shadow.predictions = data.frame(imp = rf.shadow.results$variable.importance)
      rf.shadow.predictions$feature = rownames(rf.shadow.predictions)
      # calculate correlation coefficient (and sign)
      rf.shadow.feature.cor = reshape2::melt(cor(as.matrix(shadow.data), method="spearman")) %>%
        subset(Var1 == "value") %>%
        subset(Var2 != "value") %>%
        dplyr::select(Var2, value)
      colnames(rf.shadow.feature.cor) = c("feature", "cor")
      rf.shadow.predictions = merge(rf.shadow.predictions, rf.shadow.feature.cor, type="feature")
      
      # merge
      predictions = rbind(rf.predictions %>% mutate(type = "real"),
                          rf.shadow.predictions%>% mutate(type = "shadow")) %>% as.data.frame()
      
      ## collate data
      print(paste(x, " ", y, " ", z))
      predictions = predictions %>%
        mutate(RS_Name = x,
               iter = z,
               target = y)
      predictions
    }))}))
t2 <- Sys.time()
t2 - t1 # ~1 mins; 15x faster with ranger


# calculate mean per RS (real and shadow)
ml.rf.regression.oob.oars.ph.pbs.importances.mean = ml.rf.regression.oob.oars.ph.pbs.importances %>%
  group_by(type, feature, RS_Name, iter) %>%
  group_by(type, RS_Name, feature) %>%
  mutate(mean.cor = mean(cor)) %>%
  mutate(mean.imp = mean(imp)) %>%
  mutate(imp.low = mean(imp) - (sd(imp)/sqrt(n()) * 1.96),
         imp.high = mean(imp) + (sd(imp)/sqrt(n()) * 1.96))%>%
  dplyr::select(type, RS_Name, feature, iter, mean.cor, mean.imp, imp, imp.low, imp.high) %>% distinct() %>% data.frame()

# run stats
ml.rf.regression.oob.oars.ph.pbs.importances.pval = 
  do.call(rbind, lapply(rs.names.pbs, function(rs){
  do.call(rbind, lapply(unique(ml.rf.regression.oob.oars.ph.pbs.importances.mean$feature), function(z){
      data.subset = subset(ml.rf.regression.oob.oars.ph.pbs.importances.mean, feature == z & RS_Name == rs)
      ttest.results = wilcox.test(subset(data.subset, type=="real")$imp,
                                  subset(data.subset, type=="shadow")$imp)$p.value
      ttest.results = data.frame(pval = ttest.results,
                                 RS_Name = rs,
                                 feature = z)
      ttest.results
    }))}))

ml.rf.regression.oob.oars.ph.pbs.importances.pval
ml.rf.regression.oob.oars.ph.pbs.importances.pval$padj = p.adjust(ml.rf.regression.oob.oars.ph.pbs.importances.pval$pval, method="BH")
ml.rf.regression.oob.oars.ph.pbs.importances.pval$sig = ifelse(ml.rf.regression.oob.oars.ph.pbs.importances.pval$padj < 0.05, "*", "")

# merge
ml.rf.regression.oob.oars.ph.pbs.importances.mean = merge(ml.rf.regression.oob.oars.ph.pbs.importances.mean,
                                           ml.rf.regression.oob.oars.ph.pbs.importances.pval, by=c("RS_Name", "feature"))
# delete if mean real < mean null
ml.rf.regression.oob.oars.ph.pbs.importances.mean = ml.rf.regression.oob.oars.ph.pbs.importances.mean %>%
  group_by(RS_Name, feature) %>%
  mutate(sig = ifelse(mean.imp > mean.imp[type == "shadow"], sig, ""))

# subset to sig 
ml.rf.regression.oob.oars.ph.pbs.importances.sig = subset(ml.rf.regression.oob.oars.ph.pbs.importances,
                                                          feature %in% subset(ml.rf.regression.oob.oars.ph.pbs.importances.mean, sig == "*" & mean.imp > 0.005)$feature)
ml.rf.regression.oob.oars.ph.pbs.importances.mean.sig = subset(ml.rf.regression.oob.oars.ph.pbs.importances.mean,
                                                               feature %in% subset(ml.rf.regression.oob.oars.ph.pbs.importances.mean, sig == "*" & mean.imp > 0.005)$feature) %>% dplyr::select(-iter, -imp) %>% distinct()

unique(ml.rf.regression.oob.oars.ph.pbs.importances.sig$feature) # 91

ml.rf.regression.oob.oars.ph.pbs.importances.mean.sig$cor.sign = sign(ml.rf.regression.oob.oars.ph.pbs.importances.mean.sig$mean.cor)

# adjust importance for sign
ml.rf.regression.oob.oars.ph.pbs.importances.mean.sig$adj.imp = ml.rf.regression.oob.oars.ph.pbs.importances.mean.sig$mean.imp * ml.rf.regression.oob.oars.ph.pbs.importances.mean.sig$cor.sign*100

# :: :: PLOT: feature importance ------------------------------------------

reshape2::acast(subset(ml.rf.regression.oob.oars.ph.pbs.importances.mean.sig, type == "real"),
                RS_Name ~ feature, value.var="adj.imp") %>%
  t() %>%
  pheatmap::pheatmap(color = colorRampPalette(c("red", "white", "blue"))(100),
                     annotation_col = data.frame(RS = rs.names.pbs)%>% `rownames<-`(rs.names.pbs),
                     annotation_colors = list(RS = c(labelcolors.rs, "PBS" = "black")),
                     annotation_legend = F,
                     # display_numbers = ml.rapidaim.data.ph.lfc.sig.mat.stars[,rs.names],
                     fontsize_number = 10,number_color="white",
                     breaks=c(seq(min(na.omit(ml.rf.regression.oob.oars.ph.pbs.importances.mean.sig$adj.imp)), 0, length.out=ceiling(100/2) + 1), 
                              seq(max(na.omit(ml.rf.regression.oob.oars.ph.pbs.importances.mean.sig$adj.imp))/100, max(na.omit(ml.rf.regression.oob.oars.ph.pbs.importances.mean.sig$adj.imp)), length.out=floor(100/2))),
                     border_color = "white") 
# value = importance (% decrease in accuracy, e.g. 1 = 1%) * correlation coefficient sign

# Blue = negatively correlated with fermentation (higher baseline = higher pH)
# Red = negatively correlated with fermentation (higher baseline = lower pH)


# :: FS + PBS ------------------------------------------------------------------


# use PBS method WITH feature selection

t1 <- Sys.time()
ml.rf.regression.oob.oars.ph.pbs.fs.importances <- 
  do.call(rbind, lapply(rs.names.pbs, function(x){
    # y = unique(ml.target.data$variable)[1]
    # x = rs.names[6]
    y = "pH"
    ml.target.data.subset = subset(ph.data, RS_Name == x)
    rownames(ml.target.data.subset) = ml.target.data.subset$HM
    # do not scale delta, we need absolute measurement
    ml.target.data.subset$value = (ml.target.data.subset$med.ph)
    # repeat 15 iterations
    do.call(rbind, lapply(1:15, function(z){
      # z = 1
      
      if(x != "PBS"){
        ml.input.data = ml.input.data[,colnames(ml.input.data) %in% lasso_coefs$feature]
        #ml.input.data = ml.input.data[,colnames(ml.input.data) %in% ml_rs_features_top]
      }
      
      full.data = merge(data.frame(ml.input.data), ml.target.data.subset[,c("HM", "value")], by="row.names")
      rownames(full.data) <- full.data$HM
      full.data$Row.names <- NULL
      full.data$HM <- NULL
      full.data
      
      ## Main RF
      set.seed(z)
      rf.results = ranger::ranger(value ~., full.data, importance="permutation")
      # output importances
      rf.predictions = data.frame(imp = rf.results$variable.importance)
      rf.predictions$feature = rownames(rf.predictions)
      # calculate correlation coefficient (and sign)
      rf.feature.cor = reshape2::melt(cor(as.matrix(full.data), method="spearman")) %>%
        subset(Var1 == "value") %>%
        subset(Var2 != "value") %>%
        dplyr::select(Var2, value)
      colnames(rf.feature.cor) = c("feature", "cor")
      rf.predictions = merge(rf.predictions, rf.feature.cor, type="feature")
      
      ## Shadow RF (because assuming null data would yield 50% accuracy is not empirical)
      shadow.data = full.data
      set.seed(z)
      shadow.data$value <- sample(shadow.data$value, size = nrow(shadow.data))
      set.seed(z)
      rf.shadow.results = ranger::ranger(value ~., shadow.data, importance="permutation")
      # output importances
      rf.shadow.predictions = data.frame(imp = rf.shadow.results$variable.importance)
      rf.shadow.predictions$feature = rownames(rf.shadow.predictions)
      # calculate correlation coefficient (and sign)
      rf.shadow.feature.cor = reshape2::melt(cor(as.matrix(shadow.data), method="spearman")) %>%
        subset(Var1 == "value") %>%
        subset(Var2 != "value") %>%
        dplyr::select(Var2, value)
      colnames(rf.shadow.feature.cor) = c("feature", "cor")
      rf.shadow.predictions = merge(rf.shadow.predictions, rf.shadow.feature.cor, type="feature")
      
      # merge
      predictions = rbind(rf.predictions %>% mutate(type = "real"),
                          rf.shadow.predictions%>% mutate(type = "shadow")) %>% as.data.frame()
      
      ## collate data
      print(paste(x, " ", y, " ", z))
      predictions = predictions %>%
        mutate(RS_Name = x,
               iter = z,
               target = y)
      predictions
    }))}))
t2 <- Sys.time()
t2 - t1 # ~1 mins; 15x faster with ranger


# calculate mean per RS (real and shadow)
ml.rf.regression.oob.oars.ph.pbs.fs.importances.mean = ml.rf.regression.oob.oars.ph.pbs.fs.importances %>%
  group_by(type, feature, RS_Name, iter) %>%
  group_by(type, RS_Name, feature) %>%
  mutate(mean.cor = mean(cor)) %>%
  mutate(mean.imp = mean(imp)) %>%
  mutate(imp.low = mean(imp) - (sd(imp)/sqrt(n()) * 1.96),
         imp.high = mean(imp) + (sd(imp)/sqrt(n()) * 1.96))%>%
  dplyr::select(type, RS_Name, feature, iter, mean.cor, mean.imp, imp, imp.low, imp.high) %>% distinct() %>% data.frame()

# run stats
ml.rf.regression.oob.oars.ph.pbs.fs.importances.pval = 
  do.call(rbind, lapply(rs.names.pbs, function(rs){
    do.call(rbind, lapply(unique(subset(ml.rf.regression.oob.oars.ph.pbs.fs.importances.mean, RS_Name == rs)$feature), function(z){
      data.subset = subset(ml.rf.regression.oob.oars.ph.pbs.fs.importances.mean, feature == z & RS_Name == rs)
      ttest.results = wilcox.test(subset(data.subset, type=="real")$imp,
                                  subset(data.subset, type=="shadow")$imp)$p.value
      ttest.results = data.frame(pval = ttest.results,
                                 RS_Name = rs,
                                 feature = z)
      ttest.results
    }))}))

ml.rf.regression.oob.oars.ph.pbs.fs.importances.pval
ml.rf.regression.oob.oars.ph.pbs.fs.importances.pval$padj = p.adjust(ml.rf.regression.oob.oars.ph.pbs.fs.importances.pval$pval, method="BH")
ml.rf.regression.oob.oars.ph.pbs.fs.importances.pval$sig = ifelse(ml.rf.regression.oob.oars.ph.pbs.fs.importances.pval$padj < 0.05, "*", "")

# merge
ml.rf.regression.oob.oars.ph.pbs.fs.importances.mean = merge(ml.rf.regression.oob.oars.ph.pbs.fs.importances.mean,
                                                          ml.rf.regression.oob.oars.ph.pbs.fs.importances.pval, by=c("RS_Name", "feature"))
# delete if mean real < mean null
ml.rf.regression.oob.oars.ph.pbs.fs.importances.mean = ml.rf.regression.oob.oars.ph.pbs.fs.importances.mean %>%
  group_by(RS_Name, feature) %>%
  mutate(sig = ifelse(mean.imp > mean.imp[type == "shadow"], sig, ""))

# subset to sig 
ml.rf.regression.oob.oars.ph.pbs.fs.importances.sig = subset(ml.rf.regression.oob.oars.ph.pbs.fs.importances,
                                                          feature %in% subset(ml.rf.regression.oob.oars.ph.pbs.fs.importances.mean, sig == "*" & mean.imp > 0.005)$feature)
ml.rf.regression.oob.oars.ph.pbs.fs.importances.mean.sig = subset(ml.rf.regression.oob.oars.ph.pbs.fs.importances.mean,
                                                               feature %in% subset(ml.rf.regression.oob.oars.ph.pbs.fs.importances.mean, sig == "*" & mean.imp > 0.005)$feature) %>% dplyr::select(-iter, -imp) %>% distinct()

unique(ml.rf.regression.oob.oars.ph.pbs.fs.importances.sig$feature) # 11

ml.rf.regression.oob.oars.ph.pbs.fs.importances.mean.sig$cor.sign = sign(ml.rf.regression.oob.oars.ph.pbs.fs.importances.mean.sig$mean.cor)

# adjust importance for sign
ml.rf.regression.oob.oars.ph.pbs.fs.importances.mean.sig$adj.imp = ml.rf.regression.oob.oars.ph.pbs.fs.importances.mean.sig$mean.imp * ml.rf.regression.oob.oars.ph.pbs.fs.importances.mean.sig$cor.sign*100

# :: :: PLOT: feature importance ------------------------------------------

reshape2::acast(subset(ml.rf.regression.oob.oars.ph.pbs.fs.importances.mean.sig, type == "real"),
                RS_Name ~ feature, value.var="adj.imp") %>%
  t() %>%
  pheatmap::pheatmap(color = colorRampPalette(c("red", "white", "blue"))(100),
                     annotation_col = data.frame(RS = rs.names.pbs)%>% `rownames<-`(rs.names.pbs),
                     annotation_colors = list(RS = c(labelcolors.rs, "PBS" = "black")),
                     annotation_legend = F,
                     # display_numbers = ml.rapidaim.data.ph.lfc.sig.mat.stars[,rs.names],
                     fontsize_number = 10,number_color="white",
                     breaks=c(seq(min(na.omit(ml.rf.regression.oob.oars.ph.pbs.fs.importances.mean.sig$adj.imp)), 0, length.out=ceiling(100/2) + 1), 
                              seq(max(na.omit(ml.rf.regression.oob.oars.ph.pbs.fs.importances.mean.sig$adj.imp))/100, max(na.omit(ml.rf.regression.oob.oars.ph.pbs.fs.importances.mean.sig$adj.imp)), length.out=floor(100/2))),
                     border_color = "white") 
# value = importance (% decrease in accuracy, e.g. 1 = 1%) * correlation coefficient sign

# Blue = negatively correlated with fermentation (higher baseline = higher pH)
# Red = negatively correlated with fermentation (higher baseline = lower pH)



# 4. LOOCV ----------------------------------------------------------------

#setwd("./ml_git_data/ml_loocv_data/")
ml_fs_array_files <- list.files()
ml_fs_array_files <- ml_fs_array_files[grepl("regression_results", ml_fs_array_files)]
ml_fs_array_files <- ml_fs_array_files[grepl("2025_11_07", ml_fs_array_files)]
ml_fs_array_files <- ml_fs_array_files[grepl("_fs_", ml_fs_array_files)]

# remove leading text
ml_fs_algorithms <- gsub("2025_11_07_ml_loocv_regression_results_fs_", "", gsub(".rds", "", ml_fs_array_files))
ml_fs_algorithms <- unique(data.frame(ml_fs_algorithms))
ml_fs_algorithms.df <- tidyr::separate(ml_fs_algorithms, col=ml_fs_algorithms, into=c("algo", "iter"))
table(ml_fs_algorithms.df$algo)
ml_fs_algorithms <- unique(ml_fs_algorithms.df$algo)
# list files (takes a while)
ml.fs.loocv.regression.results <- do.call(rbind, lapply(ml_fs_array_files, readRDS))

subset(ml.fs.loocv.regression.results, algo == "nnet" & iter == 15)$HM %>% table()

# now predict change in pH using predictions
ml.fs.loocv.regression.results.ph.pbs = ml.fs.loocv.regression.results %>%
  #subset(algo != "nnet") %>%
  group_by(HM, iter, data.type, algo) %>%
  mutate(true.delta = true.response - true.response[RS_Name == "PBS"],
         pred.delta = prediction - prediction[RS_Name == "PBS"]) %>%
  subset(RS_Name != "PBS")

ggplot(ml.fs.loocv.regression.results.ph.pbs %>%
         subset(data.type == "real")%>%
         mutate(RS_Name = factor(RS_Name, levels=rs.names)),
       aes(x=true.delta, y=pred.delta))+
  geom_point(aes(fill=RS_Name), shape=21)+
  scale_fill_manual(values=labelcolors.rs)+
  geom_smooth(method="lm", color="black")+
  ggpubr::stat_cor(method="spearman")+
  theme_classic()+theme(legend.position="none")+
  facet_grid(algo ~ RS_Name, scales="free")
  

# :: Evaluate AUC -------------------------------------------------------------

# calculate AUC
ml.fs.loocv.regression.results.ph.pbs.auc = ml.fs.loocv.regression.results.ph.pbs %>%
  subset(RS_Name != "PBS") %>%
  #subset(RS_Name %in% c("LetsDoOrganic", "Versafibe1490")) %>%
  mutate(true.response = ifelse(true.delta < -1.27, "strong", "weak")) %>%
  group_by(RS_Name, iter, algo, data.type) %>%
  mutate(auc = pROC::auc(as.factor(true.response), pred.delta, 
                         direction = ">", levels=c("weak", "strong")) %>% mean()) %>%
  dplyr::select(RS_Name, iter, auc) %>% distinct() %>%
  mutate(mean.auc = (auc)) %>%
  dplyr::select(RS_Name, mean.auc) %>% distinct()  %>%
  mutate(RS_Name = factor(RS_Name, levels=rs.names))

ml.fs.loocv.regression.results.ph.pbs.auc.p = 
  do.call(rbind, lapply(unique(ml.fs.loocv.regression.results.ph.pbs.auc$algo), function(algorithm){
  do.call(rbind, lapply(rs.names, function(rs){
  data.subset = subset(ml.fs.loocv.regression.results.ph.pbs.auc, RS_Name == rs & algo == algorithm)
  data.frame(RS_Name = rs,
             auc.pval = wilcox.test(subset(data.subset, data.type == "real")$mean.auc,
                                    subset(data.subset, data.type == "shadow")$mean.auc)$p.value,
             algo = algorithm)
}))})) %>%
  mutate(auc.padj = p.adjust(auc.pval, method = "bonferroni")) %>%
  mutate(auc.sig = ifelse(auc.padj < 0.05, "*", "")) %>%
  mutate(RS_Name = factor(RS_Name, levels=rs.names))

ml.fs.loocv.regression.results.ph.pbs.auc.mean = ml.fs.loocv.regression.results.ph.pbs.auc %>%
  mutate(RS_Name = factor(RS_Name, levels=rs.names))%>%
  group_by(RS_Name, algo, data.type) %>%
  mutate(mean.auc = mean(mean.auc)) %>%
  dplyr::select(RS_Name, mean.auc) %>% distinct() 

ml.fs.loocv.regression.results.ph.pbs.auc.mean = merge(ml.fs.loocv.regression.results.ph.pbs.auc.mean,
                                                       ml.fs.loocv.regression.results.ph.pbs.auc.p, by=c("RS_Name", "algo")) %>%
  mutate(RS_Name = factor(RS_Name, levels=rs.names))

# clean algo names
algo.names = data.frame(algo = c("rpart", "knn", "glmnet","pls","nnet", "svmLinear", "svmRadial", "ranger", "gbm", "xgbTree", "stack"),
                        clean = c("Decision tree", "k-Nearest neighbours", "Elastic net","Partial least squares", "Neural network", "SVM (Linear)", "SVM (Radial)", "Random forest", "Gradient boosting", "XGBoost", "Stacked ensemble"))
algo.names$clean = factor(algo.names$clean, levels=rev(algo.names$clean))

# bold the best
algo.best = ml.fs.loocv.regression.results.ph.pbs.auc.mean %>%
  subset(data.type == "real") %>%
  group_by(RS_Name) %>%
  mutate(best = ifelse(mean.auc == max(mean.auc), "bold", "plain")) %>%
  dplyr::select(RS_Name, algo, best)


# :: :: PLOT: LOOCV AUC ---------------------------------------------------

ml.fs.loocv.regression.results.ph.pbs.auc.mean.plot = ggplot(ml.fs.loocv.regression.results.ph.pbs.auc.mean %>%
         merge(algo.names, by="algo") %>%
         merge(algo.best, by=c("RS_Name", "algo")) %>%
         mutate(clean = factor(clean, levels=rev(algo.names$clean)))%>%
         subset(data.type == "real") %>%
         mutate(norm.auc = (mean.auc - min(mean.auc)) / (max(mean.auc) - min(mean.auc))),
       aes(x=RS_Name, y=clean))+
  geom_tile(aes(fill=(norm.auc)), color="white")+
  geom_text(aes(label = paste(round(mean.auc, digits=2), auc.sig, sep=""),
                fontface=best, color = ifelse(mean.auc < 0.25 | mean.auc > 0.6, 
                                              "black", "black")),
             size=3)+
  scale_y_discrete(position="left")+
  scale_fill_gradientn(colors=c("blue", "white", "red"))+
  scale_color_manual(values=c("black", "white"))+
  theme_classic(base_size = 16)+theme(legend.position="none",
                        #axis.text.x = element_text(angle=45, hjust=1),
                        axis.text.x = element_blank(),
                        strip.text = element_text(size=16))+
  facet_wrap(~"LOOCV AUC with Feature Selection")+
  labs(x=NULL, y=NULL)
ml.fs.loocv.regression.results.ph.pbs.auc.mean.plot


# :: Final predictions ----------------------------------------------------

# take optimal models
top.models = ml.fs.loocv.regression.results.ph.pbs.auc.mean %>%
  subset(data.type == "real") %>%
  group_by(RS_Name) %>%
  slice_max(mean.auc)

# isolate median prediction
top.model.preds = do.call(rbind, lapply(rs.names, function(rs){
  top.model.rs = subset(top.models, RS_Name == rs)
  data.subset = subset(ml.fs.loocv.regression.results.ph.pbs, 
                       RS_Name == rs &
                         algo == top.model.rs$algo &
                         data.type == "real")
  # calculate median prediction
  data.subset = data.subset %>%
    group_by(HM) %>%
    mutate(med.pred = median(pred.delta)) %>%
    dplyr::select(HM, RS_Name, algo, med.pred, true.delta) %>% distinct()
}))


# :: :: PLOT: Predictions -------------------------------------------------


# plot delta-pH 
top.model.preds.plot = ggplot(top.model.preds %>%
                                #group_by(HM) %>%
                                #slice_min(med.pred)%>%
                                mutate(RS_Name = factor(RS_Name, levels=rs.names))%>%
         mutate(target = "Predicted Δ pH"),
       aes(x=(RS_Name), y=(med.pred), fill=(RS_Name)))+
  ggbeeswarm::geom_beeswarm(shape=21, cex=0.75,
                            aes(fill=RS_Name), size=2)+
  geom_boxplot(notch = TRUE, outlier.shape=NA, aes(fill=RS_Name), alpha=0.3)+
  scale_fill_manual(values=labelcolors$cols[c(1:9)])+
  scale_color_manual(values=labelcolors$cols[c(1:9)])+
  scale_alpha_manual(values=c(0.3,1))+
  scale_y_continuous(limits=c(-2.9, 0.6))+
  geom_hline(yintercept = -1.27, linetype=2, linewidth = 0.2)+
  #geom_text(data=ml.rapidaim.ph.stats, aes(x=RS_Name, 
  #                                         label=ifelse(sig == "***", "*", "")), y=Inf, vjust=1.4, size=6)+
  theme_classic(base_size = 16)+theme(legend.position="none",
                        axis.text.x=element_text(angle=45,vjust=1, hjust=1),
                        legend.title = element_text(hjust=0.5),
                        strip.text = element_text(size=16),
                        strip.background = element_rect(
                          color="black"))+
  facet_wrap(~"Predicted Δ pH")+
  labs(x=NULL, y="Predicted Δ pH versus PBS")
top.model.preds.plot

# plot optimal RS (lowest pH)

# use raw real data
top.model.preds.optimal = top.model.preds %>%
  group_by(RS_Name) %>%
  mutate(scaled.pred = scale(med.pred)) %>%
  group_by(HM) %>%
  slice_min(med.pred) %>%
  mutate(is.strong = ifelse(true.delta < -1.27, "strong", "not_strong"),
         pred.strong = ifelse(med.pred < -1.27, "strong", "not_strong"))
top.model.preds.optimal$RS_Name %>% table() %>% as.data.frame() %>% mutate(perc = Freq / (nrow(subset(top.model.preds.optimal, is.strong == "strong"))))
top.model.preds.optimal[,c("is.strong", "pred.strong")] %>% table() %>%
  as.data.frame() %>%
  mutate(correct = sum(Freq[is.strong == pred.strong]),
         incorrect = sum(Freq[is.strong != pred.strong])) %>%
  dplyr::select(correct, incorrect) %>% distinct() %>%
  mutate(acc = correct / 117)
# ~72.6% accuracy

top.model.preds.cor.plot = ggplot(top.model.preds %>%
         group_by(HM) %>%
         slice_min(med.pred)%>%
         mutate(RS_Name = factor(RS_Name, levels=rs.names))%>%
         mutate(target = "Predicted Δ pH"),
       aes(x=(true.delta), y=(med.pred), fill=(RS_Name)))+
  geom_vline(xintercept = -1.27, linetype=2, linewidth = 0.2)+
  geom_hline(yintercept = -1.27, linetype=2, linewidth = 0.2)+
  geom_point(shape=21, aes(fill=RS_Name), size=2)+
  scale_fill_manual(values=labelcolors.rs)+
  geom_smooth(method="lm", color="black", aes(group=1))+
  ggpubr::stat_cor(method="spearman", aes(group=1))+
  scale_y_continuous(limits=c(-2.2, -0.9))+
  theme_classic(base_size = 16)+theme(legend.position="none",
                        #axis.text.x=element_text(angle=45,vjust=1, hjust=1),
                        legend.title = element_text(hjust=0.5),
                        strip.text = element_text(size=16),
                        strip.background = element_rect(
                          color="black"))+
  facet_wrap(~"Minimum Δ pH RS")+
  labs(x="Measured Δ pH versus PBS", y="Predicted Δ pH versus PBS")
top.model.preds.cor.plot

top.model.preds.plot+
  patchwork::free(top.model.preds.cor.plot, type="label")

# Matthew Correlation Coefficient (pearson correlation of diagonal values in a contingency table)
top.model.preds.mcc = data.frame(top.model.preds %>%
  group_by(HM) %>%
  slice_min(med.pred)%>%
  mutate(pred.response = ifelse(med.pred < -1.27, 1, 0),
         true.response = ifelse(true.delta < -1.27, 1, 0)))[,c("pred.response", "true.response")] 

caret::confusionMatrix(as.factor(top.model.preds.mcc$pred.response),
                       as.factor(top.model.preds.mcc$true.response)) 
# MCC
cor.test((top.model.preds.mcc$pred.response), 
         (top.model.preds.mcc$true.response))
# 0.397, p = 9.343e-06

(72+13) / 117

# 6. SAMPLE SIZE EXTRAPOLATION -----------------------------------------

sample.curve.size = seq(from=20, to=115, by=5)

t1 <- Sys.time()
ml.rf.regression.oob.oars.ph.pbs.fs.sample.sizes <- 
  do.call(rbind, lapply(sample.curve.size, function(sample.size){
    print(sample.size)
  do.call(rbind, parallel::mclapply(rs.names.pbs, function(x){
    # y = unique(ml.target.data$variable)[1]
    # x = rs.names[6]
    y = "pH"
    ml.target.data.subset = subset(ph.data, RS_Name == x)
    rownames(ml.target.data.subset) = ml.target.data.subset$HM
    # do not scale delta, we need absolute measurement
    ml.target.data.subset$value = (ml.target.data.subset$med.ph)
    # repeat 15 iterations
    do.call(rbind, lapply(1:100, function(z){
      # z = 1
      
      if(x != "PBS"){
        ml.input.data = ml.input.data[,colnames(ml.input.data) %in% lasso_coefs$feature]
      }
      
      full.data = merge(data.frame(ml.input.data), ml.target.data.subset[,c("HM", "value")], by="row.names")
      rownames(full.data) <- full.data$HM
      full.data$Row.names <- NULL
      full.data$HM <- NULL
      full.data
      
      # subset samples
      set.seed(z)
      samples.to.keep = sample(1:nrow(full.data), size=sample.size)
      full.data = full.data[sort(samples.to.keep),]
      ml.target.data.subset = subset(ml.target.data.subset, HM %in% rownames(full.data))
      
      
      ## Main RF
      set.seed(z)
      rf.results = ranger::ranger(value ~., full.data, importance="permutation")
      rf.predictions = data.frame(pred = as.vector(rf.results$predictions)) %>%
        mutate(HM = rownames(full.data))
      
      ## collate data
      print(paste(x, " ", y, " ", z, " ", sample.size))
      rf.predictions = rf.predictions %>%
        mutate(RS_Name = x,
               iter = z,
               size = sample.size)
      rf.predictions
    }))}))}))
t2 <- Sys.time()
t2 - t1 # 15 min


# now predict change in pH using predictions
ml.sample.size.curve = ml.rf.regression.oob.oars.ph.pbs.fs.sample.sizes %>%
  merge(ph.data, by=c("HM", "RS_Name")) %>%
  mutate(true = med.ph) %>%
  group_by(HM, iter, size) %>%
  mutate(true.delta = true - true[RS_Name == "PBS"],
         pred.delta = pred - pred[RS_Name == "PBS"]) %>%
  subset(RS_Name != "PBS")

subset(ml.sample.size.curve, iter == 1 & RS_Name == "Novelose330" & HM == "HM0239.00")


# :: :: Extrapolate AUC ---------------------------------------------------


# calculate AUCs
ml.sample.size.curve.auc = 
  do.call(rbind, lapply(rs.names, function(rs){
  do.call(rbind, lapply(sample.curve.size, function(sample.size){
    print(paste0(rs, " ", sample.size))
    do.call(rbind, lapply(1:100, function(seed){
      # subset
      data.subset = ml.sample.size.curve %>%
        subset(RS_Name == rs & size == sample.size & iter == seed) %>%
        mutate(response = ifelse(true.delta < -1.27, "strong", "weak")) %>%
        mutate(response = factor(response, levels=c("strong", "weak")))
      
      # calculate MSE (for Back-calculate section)
      mse.rs =  mean((data.subset$true.delta - data.subset$pred.delta)^2)
      
      # skip 1 class situations
      if(length(table(as.character(data.subset$response)))!=2){
        auc = NA
        cor = NA
      }else{
      # calculate AUC
      auc = pROC::auc(data.subset$response, 
                      data.subset$pred.delta, 
                      direction = "<") %>% suppressMessages()
      cor = cor.test(data.subset$pred.delta,
                     data.subset$true.delta)$estimate
      }
      # export
       data.frame(RS_Name = rs,
                  size = sample.size,
                  iter = seed,
                  auc = auc,
                  cor = cor,
                  mse = mse.rs) 
    }))}))}))
        
ml.sample.size.curve.lm = 
  do.call(rbind, lapply(rs.names, function(rs){
  #do.call(rbind, lapply(1:100, function(seed){
    print(paste0(rs))
    data.subset = subset(ml.sample.size.curve.auc, RS_Name == rs)
    
    
    do.call(rbind, lapply(1:15, function(seed){
    # sample 15 data points per sample size
    set.seed(seed)
    data.subset = data.subset %>%
      group_by(RS_Name, size) %>%
      slice_sample(n = 15)
    
    
    if(sum(!is.na(data.subset$auc))>=3){
    # build LINEAR model
    linear.model = lm(auc ~ size, na.omit(data.subset))
    # build EXPONENTIAL model
    nonlinear.model = lm(auc ~ log10(size),na.omit(data.subset))

    # apply models
    sample.sizes = seq(from = 20, to = 1000, by=10)
    do.call(rbind, lapply(sample.sizes, function(x){
      lm.pred = predict(linear.model, data.frame(size = x))
      nlm.pred = (predict(nonlinear.model, data.frame(size = x)))
      data.frame(RS_Name = rs,
                 iter = seed,
                 lm.prediction = lm.pred,
                 nlm.prediction = nlm.pred,
                 sample.size = x)
    }))
    }else{
    data.frame(RS_Name = rs,
               iter = seed,
               lm.prediction = NA,
               nlm.prediction = NA,
               sample.size = NA)
    }      
    }))
    }))


ml.sample.size.curve.lm$RS_Name = factor(ml.sample.size.curve.lm$RS_Name, levels = rs.names)
ml.sample.size.curve.auc$RS_Name = factor(ml.sample.size.curve.auc$RS_Name, levels = rs.names)

ml.sample.size.curve.linear <- ggplot()+
  # linear points
  geom_point(data=ml.sample.size.curve.lm,
             aes(x=(sample.size), y=lm.prediction, color=RS_Name), alpha=0.1)+
  # original points
  geom_point(data=ml.sample.size.curve.auc,
             aes(x=(size), y=auc, color=RS_Name), alpha=0.1)+
  # lines
  geom_smooth(data=ml.sample.size.curve.lm%>% subset(RS_Name != "Versafibe1490"), 
              aes(x=(sample.size), y=lm.prediction), color="black", linetype = 1, method="loess", alpha=1, linewidth=0.5)+
  #geom_smooth(data=ml.sample.size.curve.data.lm, 
  #            aes(x=(size), y=auc), color="black", linetype = 1, method="lm", alpha=1, linewidth=0.5)+
  
  # other
  geom_hline(yintercept=0.90, linewidth=0.2, linetype=2)+
  scale_color_manual(values=scales::muted(labelcolors$cols[c(1:9)], l=70))+
  xlim(1,(1000))+
  ylim(0,1)+
  scale_x_log10()+
  facet_wrap(~RS_Name, nrow=1)+
  geom_text(data = subset(ml.sample.size.curve.lm, lm.prediction < 0.90) %>% 
              subset(RS_Name != "Versafibe1490") %>%
              group_by(RS_Name, iter) %>% 
              slice_max(lm.prediction) %>%
              group_by(RS_Name) %>%
              mutate(mean.size = mean(sample.size)) %>% dplyr::select(RS_Name, mean.size) %>% distinct() ,
            aes(x=75, y=0.25, label=round(mean.size, digits=0)), size=4)+
  geom_segment(data = subset(ml.sample.size.curve.lm, lm.prediction < 0.90) %>% 
                 subset(RS_Name != "Versafibe1490") %>%
                 group_by(RS_Name, iter) %>% 
                 slice_max(lm.prediction) %>%
                 group_by(RS_Name) %>%
                 mutate(mean.size = mean(sample.size)) %>% dplyr::select(RS_Name, mean.size) %>% distinct() ,
               aes(x=mean.size, xend=mean.size, y=0, yend=1),
               linewidth=0.2, linetype=2)+
  geom_segment(data = subset(ml.sample.size.curve.lm, lm.prediction < 0.90) %>% 
                 subset(RS_Name != "Versafibe1490") %>%
                 group_by(RS_Name, iter) %>% 
                 slice_max(lm.prediction) %>%
                 group_by(RS_Name) %>%
                 mutate(mean.size = mean(sample.size)) %>% dplyr::select(RS_Name, mean.size) %>% distinct() ,
               aes(x=100, xend=mean.size*0.85, y=0.15, yend=0),
               arrow = arrow(length=unit(0.15,"cm"), ends="last", type = "closed"))+
  labs(x="Sample Size Assumption (Linear Model)", y="Extrapolated AUC")+
  theme_classic()+theme(legend.position="none",
                        panel.grid.minor = element_blank(),
                        legend.title = element_text(hjust=0.5),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
ml.sample.size.curve.linear


ml.sample.size.curve.nonlinear <- ggplot()+
  # linear points
  geom_point(data=ml.sample.size.curve.lm,
             aes(x=(sample.size), y=nlm.prediction, color=RS_Name), alpha=0.1)+
  # original points
  geom_point(data=ml.sample.size.curve.auc,
             aes(x=(size), y=auc, color=RS_Name), alpha=0.1)+
  # lines
  geom_smooth(data=ml.sample.size.curve.lm %>% subset(RS_Name != "Versafibe1490"), 
              aes(x=(sample.size), y=nlm.prediction), color="black", linetype = 1, method="loess", alpha=1, linewidth=0.5)+
  #geom_smooth(data=ml.sample.size.curve.data.lm, 
  #            aes(x=(size), y=auc), color="black", linetype = 1, method="lm", alpha=1, linewidth=0.5)+
  
  # other
  geom_hline(yintercept=0.90, linewidth=0.2, linetype=2)+
  scale_color_manual(values=scales::muted(labelcolors$cols[c(1:9)], l=70))+
  xlim(1,(1000))+
  ylim(0,1)+
  scale_x_log10()+
  facet_wrap(~RS_Name, nrow=1)+
  geom_text(data = subset(ml.sample.size.curve.lm, nlm.prediction < 0.90) %>% 
              subset(RS_Name != "Versafibe1490") %>%
              group_by(RS_Name, iter) %>% 
              slice_max(nlm.prediction) %>%
              group_by(RS_Name) %>%
              mutate(mean.size = mean(sample.size)) %>% dplyr::select(RS_Name, mean.size) %>% distinct() ,
            aes(x=75, y=0.25, label=round(mean.size, digits=0)), size=4)+
  geom_segment(data = subset(ml.sample.size.curve.lm, nlm.prediction < 0.90) %>% 
                 subset(RS_Name != "Versafibe1490") %>%
                 group_by(RS_Name, iter) %>% 
                 slice_max(nlm.prediction) %>%
                 group_by(RS_Name) %>%
                 mutate(mean.size = mean(sample.size)) %>% dplyr::select(RS_Name, mean.size) %>% distinct() ,
               aes(x=mean.size, xend=mean.size, y=0, yend=1),
               linewidth=0.2, linetype=2)+
  geom_segment(data = subset(ml.sample.size.curve.lm, nlm.prediction < 0.90) %>% 
                 subset(RS_Name != "Versafibe1490") %>%
                 group_by(RS_Name, iter) %>% 
                 slice_max(nlm.prediction) %>%
                 group_by(RS_Name) %>%
                 mutate(mean.size = mean(sample.size)) %>% dplyr::select(RS_Name, mean.size) %>% distinct() ,
               aes(x=100, xend=mean.size*0.85, y=0.15, yend=0),
               arrow = arrow(length=unit(0.15,"cm"), ends="last", type = "closed"))+
  labs(x="Sample Size Assumption (Non-linear Model)", y="Extrapolated AUC")+
  theme_classic()+theme(legend.position="none",
                        panel.grid.minor = element_blank(),
                        legend.title = element_text(hjust=0.5),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
ml.sample.size.curve.nonlinear

ml.sample.size.curve.linear/
  ml.sample.size.curve.nonlinear


# when do RS hit AUC = 0.90?
ml.sample.size.curve.lm.values = subset(ml.sample.size.curve.lm, lm.prediction < 0.90) %>% 
  group_by(RS_Name, iter) %>% 
  slice_max(lm.prediction) %>%
  group_by(RS_Name) %>%
  mutate(mean.size = mean(sample.size),
         imp.low = mean(sample.size) - (sd(sample.size)/sqrt(n()) * 1.96),
         imp.high = mean(sample.size) + (sd(sample.size)/sqrt(n()) * 1.96))%>%
  dplyr::select(RS_Name, mean.size, imp.low,imp.high) %>% distinct() 

subset(ml.sample.size.curve.lm, (nlm.prediction) < 0.90) %>% 
  group_by(RS_Name, iter) %>% 
  slice_max((nlm.prediction)) %>%
  group_by(RS_Name) %>%
  mutate(mean.size = mean(sample.size),
         imp.low = mean(sample.size) - (sd(sample.size)/sqrt(n()) * 1.96),
         imp.high = mean(sample.size) + (sd(sample.size)/sqrt(n()) * 1.96))%>%
  dplyr::select(RS_Name, mean.size, imp.low,imp.high) %>% distinct() 

# subsampled 15 AUC values per sample.size (20 - 115)
# build a LM with those n=15 per size, predicted (120 - 1000); 15 times
# selected mean of 15 predictions that reach AUC 0.90


# :: :: PLOT: Extrapolations ----------------------------------------------

# Just show LDO's curve, then show others as heatmap or something

ml.sample.size.curve.linear.ldo <- ggplot()+
  # linear points
  geom_point(data=ml.sample.size.curve.lm %>% subset(RS_Name == "LetsDoOrganic"),
             aes(x=(sample.size), y=lm.prediction, color=RS_Name), alpha=0.1)+
  # original points
  geom_point(data=ml.sample.size.curve.auc%>% subset(RS_Name == "LetsDoOrganic"),
             aes(x=(size), y=auc, color=RS_Name), alpha=0.1)+
  # lines
  geom_smooth(data=ml.sample.size.curve.lm%>% 
                subset(RS_Name == "LetsDoOrganic"), 
              aes(x=(sample.size), y=lm.prediction), color="black", linetype = 1, method="loess", alpha=1, linewidth=0.5)+
  #geom_smooth(data=ml.sample.size.curve.data.lm, 
  #            aes(x=(size), y=auc), color="black", linetype = 1, method="lm", alpha=1, linewidth=0.5)+
  
  # other
  geom_hline(yintercept=0.90, linewidth=0.2, linetype=2)+
  scale_color_manual(values=labelcolors.rs)+
  xlim(1,(1000))+
  ylim(0,1)+
  scale_x_log10()+
  facet_wrap(~RS_Name, nrow=1)+
  geom_text(data = subset(ml.sample.size.curve.lm, lm.prediction < 0.90) %>% 
              subset(RS_Name == "LetsDoOrganic") %>%
              group_by(RS_Name, iter) %>% 
              slice_max(lm.prediction) %>%
              group_by(RS_Name) %>%
              mutate(mean.size = mean(sample.size)) %>% dplyr::select(RS_Name, mean.size) %>% distinct() ,
            aes(x=75, y=0.25, label=round(mean.size, digits=0)), size=4)+
  geom_segment(data = subset(ml.sample.size.curve.lm, lm.prediction < 0.90) %>% 
                 subset(RS_Name == "LetsDoOrganic") %>%
                 group_by(RS_Name, iter) %>% 
                 slice_max(lm.prediction) %>%
                 group_by(RS_Name) %>%
                 mutate(mean.size = mean(sample.size)) %>% dplyr::select(RS_Name, mean.size) %>% distinct() ,
               aes(x=mean.size, xend=mean.size, y=0, yend=1),
               linewidth=0.2, linetype=2)+
  geom_segment(data = subset(ml.sample.size.curve.lm, lm.prediction < 0.90) %>% 
                 subset(RS_Name == "LetsDoOrganic") %>%
                 group_by(RS_Name, iter) %>% 
                 slice_max(lm.prediction) %>%
                 group_by(RS_Name) %>%
                 mutate(mean.size = mean(sample.size)) %>% dplyr::select(RS_Name, mean.size) %>% distinct() ,
               aes(x=100, xend=mean.size*0.85, y=0.15, yend=0),
               arrow = arrow(length=unit(0.15,"cm"), ends="last", type = "closed"))+
  labs(x="Sample Size Assumption", y="Extrapolated AUC")+
  theme_classic(base_size = 16)+theme(legend.position="none",
                        panel.grid.minor = element_blank(),
                        legend.title = element_text(hjust=0.5),
                        strip.text = element_text(size=16),
                        strip.background = element_rect(
                          color="black"))
ml.sample.size.curve.linear.ldo

ml.sample.size.curve.lm.values.plot = ggplot(ml.sample.size.curve.lm.values %>%
         mutate(plot.version = ifelse(RS_Name !="Versafibe1490", "B", "A")),
       aes(x=RS_Name, y=mean.size))+
  geom_segment(aes(y=imp.low, yend=imp.high, color=RS_Name), 
               linewidth=1, alpha=0.75)+
  geom_point(aes(fill=RS_Name), shape=21, size=2.5)+
  scale_color_manual(values=labelcolors.rs)+
  scale_fill_manual(values=labelcolors.rs)+
  facet_grid(plot.version~., scales="free")+
  theme_classic(base_size=16)+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        legend.position="none",
        strip.background = element_blank(),
        strip.text = element_blank())+
  labs(x=NULL, y="Extrapolated Sample Size")

ml.sample.size.curve.linear.ldo+
ml.sample.size.curve.lm.values.plot

# :: Back-calculate --------------------------------------------------------

# 1. calculate MSE of predictions for subsamples [done with AUC]
# 2. extrapolate MSE of larger sample sizes
# 3. calibrate MSE penalization using true values
# 4. calculate sample size needed to outperform LDO default (10%, 20%, 30% higher)

# :: :: Extrapolate MSE ---------------------------------------------------

# extrapolate MSE of larger sample sizes

ml.sample.size.curve.mse.lm = 
  do.call(rbind, lapply(rs.names, function(rs){
    #do.call(rbind, lapply(1:100, function(seed){
    print(paste0(rs))
    data.subset = subset(ml.sample.size.curve.auc, RS_Name == rs)
    
    
    do.call(rbind, lapply(1:15, function(seed){
      # sample 15 data points per sample size
      set.seed(seed)
      data.subset = data.subset %>%
        group_by(RS_Name, size) %>%
        slice_sample(n = 15)
      
      
      if(sum(!is.na(data.subset$auc))>=3){
        # build LINEAR model
        linear.model = lm(mse ~ size, na.omit(data.subset))
        # build LOG-LOG model
        nonlinear.model = lm(log10(mse) ~ log10(size), na.omit(data.subset))
        
        # apply models
        sample.sizes = seq(from = 20, to = 1000, by=10)
        do.call(rbind, lapply(sample.sizes, function(x){
          lm.pred = predict(linear.model, data.frame(size = x))
          nlm.pred = predict(nonlinear.model, data.frame(size = x))
          data.frame(RS_Name = rs,
                     iter = seed,
                     lm.prediction = lm.pred,
                     nlm.prediction = nlm.pred,
                     lm.learn.rate = coef(linear.model)[2],
                     nlm.learn.rate = coef(nonlinear.model)[2],
                     sample.size = x)
        }))
      }else{
        data.frame(RS_Name = rs,
                   iter = seed,
                   lm.prediction = NA,
                   nlm.prediction = NA,
                   sample.size = NA)
      }      
    }))
  }))


ml.sample.size.curve.linear.mse <- ggplot()+
  # linear points
  geom_point(data=ml.sample.size.curve.mse.lm %>% subset(lm.prediction > 0) %>%
               mutate(RS_Name = factor(RS_Name, levels=rs.names)),
             aes(x=(sample.size), y=lm.prediction, color=RS_Name), alpha=0.1)+
  # original points
  geom_point(data=ml.sample.size.curve.auc %>%
               mutate(RS_Name = factor(RS_Name, levels=rs.names)),
             aes(x=(size), y=mse, color=RS_Name), alpha=0.1)+
  # lines
  geom_smooth(data=ml.sample.size.curve.mse.lm  %>% subset(lm.prediction > 0) %>%
                mutate(RS_Name = factor(RS_Name, levels=rs.names)),
              method="lm",
              aes(x=(sample.size), y=lm.prediction), color="black", linetype = 1,alpha=1, linewidth=0.5)+
  #geom_smooth(data=ml.sample.size.curve.data.lm, 
  #            aes(x=(size), y=auc), color="black", linetype = 1, method="lm", alpha=1, linewidth=0.5)+
  
  # other
 # geom_hline(yintercept=0.90, linewidth=0.2, linetype=2)+
  scale_color_manual(values=scales::muted(labelcolors$cols[c(1:9)], l=70))+
  #xlim(1,(1000))+
  ylim(c(0,0.76))+
  #scale_x_log10()+
  facet_wrap(~RS_Name, nrow=1)+
  labs(x="Sample Size Assumption (Linear Model)", y="Extrapolated MSE")+
  theme_classic()+theme(legend.position="none",
                        panel.grid.minor = element_blank(),
                        legend.title = element_text(hjust=0.5),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
ml.sample.size.curve.linear.mse


# :: :: Calibrate parameters -------------------------------------------------

# calibrate MSE penalization using true values
ml.sample.size.curve.mse.lm

# first, generate raw values that can be transformed outside
ml.sample.size.curve.mse.penalty = do.call(rbind, lapply(1:15, function(seed){
  do.call(rbind, lapply(unique(subset(ml.sample.size.curve.mse.lm, sample.size < 120)$sample.size), function(s){
    print(paste0(seed, " ", s))
    data.subset.rs = do.call(rbind, parallel::mclapply(rs.names, function(rs){
      # subset
      data.subset.rs = subset(ml.sample.size.curve.mse.lm, RS_Name == rs & iter == seed & sample.size == s)
      # 
      # shrink noise i.e. difference(true - pred)
      # based on MSE coefficient
      data.subset.true = subset(ml.sample.size.curve, RS_Name == rs & iter == seed & size == s)
      cor.test(data.subset.true$true.delta, data.subset.true$pred.delta)
      set.seed(seed)
      # noise
      mse.val = data.subset.rs$lm.prediction
      # if mse is negative, transform to shrink it
      if(sign(mse.val) == -1){
        mse.val = 0
      }
      
      disrupt.coef = sample(rnorm(n = nrow(data.subset.true), 
                                  mean = 0, 
                                  sd = sqrt(mse.val)  / (1 + exp(-sqrt((mse.val))))))
      # record noise factor 
      data.subset.true$disrupt.true = disrupt.coef

      #data.subset.true$true.noise = data.subset.true$true.delta - data.subset.true$disrupt.true

      #cor.test(data.subset.true$true.noise,
      #         data.subset.true$true.delta)
      
      data.subset.true %>%
        mutate(RS_Name = rs,
               iter = seed,
               sample.size=s,
               lm.learn.rate = data.subset.rs$lm.learn.rate)
      
    }))

  }))
}))

penalties = seq(from=0.5,
                to=10,
                by=0.5)

ml.sample.size.curve.mse.penalty[,c("sample.size", "disrupt.true", "RS_Name")] %>% 
  group_by(sample.size, RS_Name) %>% 
  mutate(mean.d = mean(disrupt.true)) %>% 
  dplyr::select(sample.size, mean.d) %>% distinct() %>%
  ggplot(aes(x=sample.size, y=mean.d))+
  geom_point(shape=21)+
  geom_smooth(method="lm", color="black")+
  theme_classic()+
  labs(x="Sample Size", y="Mean Noise")
# noise gets smaller as N increases (good)

# now we want to calibrate the noise so that
# "synthesized predictions" (i.e. true(noise) )
# correlate with the truth to a similar degree as real predictions

penalties = seq(from=0, to=10, by=0.01)

#penalties = 5
ml.sample.size.curve.mse.penalty.plot = do.call(rbind, lapply(penalties, function(p){
ml.sample.size.curve.mse.penalty %>%
  group_by(iter, sample.size) %>%
    mutate(adj.pred = (true.delta - (true.delta * disrupt.true * p/(sample.size^(1/p*0.56))))) %>%
    #  mutate(adj.pred = (true.delta + disrupt.true * (p / log10(sample.size*p))) %>%
    mutate(cor.real = cor.test(true.delta, pred.delta)$estimate,
         cor.pred = cor.test(true.delta, adj.pred)$estimate,
         cor.delta = cor.real - cor.pred,
         penalty = p)
  })) %>% #as.data.frame()
  dplyr::select(sample.size, cor.delta, cor.real, cor.pred, penalty) %>% distinct() %>%
  ungroup() %>%
  dplyr::select(-iter) %>%
  group_by(penalty)%>%
  # calculate y-intercept and record slope
  mutate(yint = coef(lm(cor.delta ~ sample.size))[1],
         slope = cor.test(cor.delta, sample.size)$estimate) %>%
  dplyr::select(yint, slope, penalty) %>% distinct() %>%
  mutate(score = yint-slope)%>%
  ggplot(aes(x=yint, y=slope))+
  #geom_line(linetype=2, linewidth=0.2)+
  geom_point(shape=21,  aes(color=penalty))+
  scale_color_gradient2(low="blue", mid="white", high="red")+
  #ggpubr::stat_cor(method="pearson")+
  geom_hline(yintercept=0, color="red",linetype=2)+
  geom_vline(xintercept=0, color="red",linetype=2)+
  theme_classic()+theme(legend.position="inside",
                        legend.position.inside=c(0.8,0.8))+
  labs(x="Y-Intercept", y="Slope (Δ R)", color="Penalty")
ml.sample.size.curve.mse.penalty.plot

ml.sample.size.curve.mse.penalty.plot.2 = do.call(rbind, lapply(c(2,3,4,4.07, 5,6,7,8,9), function(p){
  ml.sample.size.curve.mse.penalty %>%
    group_by(iter, sample.size) %>%
    mutate(adj.pred = (true.delta - (true.delta * disrupt.true * p/(sample.size^(1/p*0.56))))) %>%
    mutate(cor.real = cor.test(true.delta, pred.delta)$estimate,
           cor.pred = cor.test(true.delta, adj.pred)$estimate,
           cor.delta = cor.real - cor.pred,
           penalty = p)
})) %>%
  dplyr::select(RS_Name, sample.size, cor.delta, cor.real, cor.pred, penalty) %>% distinct() %>%
  group_by(penalty) %>%
  mutate(cor.lambda = cor.test(sample.size, cor.delta)$estimate,
         p.lambda = cor.test(sample.size, cor.delta)$p.value) %>%
  dplyr::select(cor.lambda, p.lambda, penalty, sample.size, cor.delta) %>% distinct() %>%
  ggplot(aes(x=sample.size, y=cor.delta))+
  geom_point(shape=21)+
  geom_hline(yintercept=0, linetype=2, color="red")+
  geom_smooth(method="lm", color="orange")+
  ggpubr::stat_cor()+
  #geom_vline(xintercept=16.75, linetype=2, color="red")+
  theme_classic()+
  facet_wrap(~penalty)+
  labs(x="Sample Size", y="Δ R")
# Parameters identified!
ml.sample.size.curve.mse.penalty.plot.2

# now apply calibration and predict

t1 = Sys.time()
ml.sample.size.curve.mse.scores = do.call(rbind, lapply(1:15, function(seed){
  do.call(rbind, lapply(unique(ml.sample.size.curve.mse.lm$sample.size), function(s){
    print(paste0(seed, " ", s))
    data.subset.rs = do.call(rbind, parallel::mclapply(rs.names, function(rs){
      # subset
      data.subset.rs = subset(ml.sample.size.curve.mse.lm, RS_Name == rs & iter == seed & sample.size == s)
      # 
      # shrink noise i.e. difference(true - pred)
      # based on MSE coefficient
      data.subset.true = subset(ml.sample.size.curve, RS_Name == rs & iter == seed & size == 115)

      set.seed(seed)
      # noise
      mse.val = data.subset.rs$lm.prediction
      # if mse is negative, shrink to 0
      if(sign(mse.val) == -1){
        mse.val = 0
      }
      
      
      disrupt.coef = sample(rnorm(n = nrow(data.subset.true), 
                                  mean = 0, 
                                  sd = sqrt(abs(mse.val)) / (1 + exp(-sqrt(abs(mse.val))))))
      
      p = 4.07 # empirically calibrated
      
      data.subset.true$disrupt.coef = disrupt.coef
      
      data.subset.true = data.subset.true %>%
          mutate(disrupt.true = (true.delta - (true.delta * disrupt.coef * p/(s^(1/p*0.56)))))

      data.subset.true
    }))
    # select lowest prediction
    rs.score = data.subset.rs %>% 
      group_by(HM) %>%
      slice_min(disrupt.true) %>%
      mutate(true.response = ifelse(true.delta < -1.27, "strong", "weak"))
    # LDO default
    ldo.score = data.subset.rs %>%
      subset(RS_Name == "LetsDoOrganic") %>%
      mutate(true.response = ifelse(true.delta < -1.27, "strong", "weak"))
    # compare accuracy rates
    data.frame(pred.rate = nrow(rbind(subset(rs.score, true.response == "strong" & disrupt.true < -1.27),
                                      subset(rs.score, true.response == "weak" & disrupt.true > -1.27)))/
                 nrow(rs.score), # calculate total correct, not just TP (i.e. TN, too)
               # because LDO will only be TP + FP (no TN)
               ldo.rate = nrow(subset(ldo.score, true.response == "strong" & true.delta < -1.27))/nrow(rs.score),
               sample.size = s,
               iter = seed)
  }))
}))
Sys.time() - t1 # ~20 min

# ldo only
ldo.default = ml.target.data %>%
  subset(variable == "delta_pH")%>%
  group_by(HM) %>%
  subset(RS_Name == "LetsDoOrganic") %>%
  mutate(response = ifelse(value < -1.27, "strong", "weak"))
ldo.default = table(ldo.default$response)[1] / 117

ml.sample.size.curve.mse.scores.plot = ggplot(ml.sample.size.curve.mse.scores %>%
         mutate(ldo.rate = ldo.default),
       aes(x = sample.size))+
  geom_line(aes(y=pred.rate*100, color="Prediction", group=iter), alpha=0.2)+
  #geom_line(aes(y=ldo.rate, color="LetsDoOrganic",group=iter))+
  geom_smooth(aes(y=ldo.rate*100, color="LetsDoOrganic"), method="lm")+
  geom_smooth(aes(y=pred.rate*100, color="Prediction"))+
  scale_color_manual(values=c(labelcolors.rs, "Prediction" = "black"))+
  scale_y_continuous(limits=c(65, 100), breaks=seq(from=65, to=100, by=5))+
  scale_x_continuous(breaks=seq(from=0, to=1000, by=100))+
  #scale_x_log10(breaks=c(seq(from=100, to=700, by=100)))+
  theme_classic(base_size = 16)+theme(#axis.text.x=element_text(angle=45, hjust=1),
    
                        legend.position = c(0.8,0.4),
                        panel.grid.major = element_line(linewidth=0.2, color="grey", linetype=2))+
  labs(x="Sample Size Assumption", 
       y="Correct Prediction Rate",
       color = NULL)
ml.sample.size.curve.mse.scores.plot

# decelerates after ~n=400


# :: Tabulate performance -------------------------------------------------

# for results section, what are Accuracy TP, TP, Sens, Spec, etc

# LDO Default "one-size-fits-some"
ldo.default = ml.target.data %>%
  subset(variable == "delta_pH")%>%
  group_by(HM) %>%
  subset(RS_Name == "LetsDoOrganic") %>%
  mutate(response = ifelse(value < -1.27, "strong", "weak"))
ldo.default = table(ldo.default$response)[1] / 117

# Sens (TP / TP+FN) = 100 / 100

# Spec (TN / TN+FP) = 0 / 100


# RF OOB model performances
rf.model.preds.optimal = ml.rf.regression.oob.oars.ph.pbs.fs %>% 
  group_by(RS_Name, HM) %>%
  mutate(med.pred = mean(pred.delta)) %>%
  dplyr::select(HM, true.delta, med.pred, RS_Name) %>% distinct() %>%
  group_by(HM) %>%
  slice_min(med.pred) %>%
  mutate(is.strong = ifelse(true.delta < -1.27, "strong", "not_strong"),
         pred.strong = ifelse(med.pred < -1.27, "strong", "not_strong"))
rf.model.preds.optimal$RS_Name %>% table() %>% as.data.frame() %>% mutate(perc = Freq / (nrow((rf.model.preds.optimal))))
rf.model.preds.optimal[,c("is.strong", "pred.strong")] %>% table() %>%
  as.data.frame() %>%
  mutate(correct = sum(Freq[is.strong == pred.strong]),
         incorrect = sum(Freq[is.strong != pred.strong])) %>%
  dplyr::select(correct, incorrect) %>% distinct() %>%
  mutate(acc = correct / 117)
# ~70.1% accuracy

caret::confusionMatrix(table(rf.model.preds.optimal[,c("is.strong", "pred.strong")]),
                       positive = "strong")
# Sensitivity : 0.6903         
# Specificity : 1.0000

# optimal model performances
top.model.preds.optimal = top.model.preds %>%
  group_by(RS_Name) %>%
  mutate(scaled.pred = scale(med.pred)) %>%
  group_by(HM) %>%
  slice_min(med.pred) %>%
  mutate(is.strong = ifelse(true.delta < -1.27, "strong", "not_strong"),
         pred.strong = ifelse(med.pred < -1.27, "strong", "not_strong"))
top.model.preds.optimal$RS_Name %>% table() %>% as.data.frame() %>% mutate(perc = Freq / (nrow(top.model.preds.optimal)))
top.model.preds.optimal[,c("is.strong", "pred.strong")] %>% table() %>%
  as.data.frame() %>%
  mutate(correct = sum(Freq[is.strong == pred.strong]),
         incorrect = sum(Freq[is.strong != pred.strong])) %>%
  dplyr::select(correct, incorrect) %>% distinct() %>%
  mutate(acc = correct / 117)
# ~72.6% accuracy

caret::confusionMatrix(table(top.model.preds.optimal[,c("is.strong", "pred.strong")]),
                       positive = "strong")
# Sensitivity : 0.7059          
# Specificity : 0.8667

# visualize conf matrix
table(top.model.preds.optimal[,c("is.strong", "pred.strong")]) %>%
  as.data.frame() %>%
  mutate(perc = round(Freq / sum(Freq), digits=3)*100)%>%
  mutate(is.strong = factor(is.strong, levels=c("strong", "not_strong")),
         pred.strong = factor(pred.strong, levels=c("strong", "not_strong"))) %>%
  ggplot(aes(x=is.strong, y=pred.strong))+
  geom_tile(aes(fill=perc), color="white")+
  geom_text(aes(label=perc), color="white")+
  scale_fill_gradient2(low="blue", mid="white", midpoint = 35, high="red")+
  theme_classic()+theme(legend.position="none")+
  labs(x="True", y="Predicted")

(61.5 + 11.1) / 100

# FN rate = 2.7%
(1.7) / (1.7 + 61.5)
# FP rate = 29%
25.6 / (61.5 + 25.6)
# TP rate = 70.6%
61.5 / (61.5 + 25.6)
# TN rate = 30.2%
11.1 / (11.1 + 25.6)

# Sens (TP / TP+FN) = 0.706
61.5 / (61.5+25.6)
# Spec (TN / TN+FP) = 0.867
11.1 / (11.1+1.7)


# visualize perturbation value

perturbation.function = function(x){
  if(x < 0){
    x = 0
  }
  sqrt(x) / (1 + exp(-sqrt(x)))
}

data.frame(x = seq(-1, 1, by=0.01)) %>%
  group_by(x)%>%
  mutate(y = perturbation.function(x)) %>%
  ggplot(aes(x=x, y=y))+
  geom_point()+
  theme_classic()

# >> FIGURES --------------------------------------------------------------


## TOP
(ml.rapidaim.ph.plot+
   patchwork::free(dbrda_model_ph_plot, type="label")+
   patchwork::free(ml.rf.regression.oob.oars.ph.pbs.plot.neat, type="label")) %>%
  ggsave(
    filename="../../ml_plots/oars_ml_plot_top.pdf",
    width=14,
    height=5,
    units="in",
    device=cairo_pdf)

reshape2::acast(subset(ml.rf.regression.oob.oars.ph.pbs.importances.mean.sig, type == "real"),
                RS_Name ~ feature, value.var="adj.imp") %>%
  t() %>%
  pheatmap::pheatmap(color = colorRampPalette(c("red", "white", "blue"))(100),
                     annotation_col = data.frame(RS = rs.names.pbs)%>% `rownames<-`(rs.names.pbs),
                     annotation_colors = list(RS = c(labelcolors.rs, "PBS" = "black")),
                     annotation_legend = F,
                     # display_numbers = ml.rapidaim.data.ph.lfc.sig.mat.stars[,rs.names],
                     fontsize_number = 10,number_color="white",
                     breaks=c(seq(min(na.omit(ml.rf.regression.oob.oars.ph.pbs.importances.mean.sig$adj.imp)), 0, length.out=ceiling(100/2) + 1), 
                              seq(max(na.omit(ml.rf.regression.oob.oars.ph.pbs.importances.mean.sig$adj.imp))/100, max(na.omit(ml.rf.regression.oob.oars.ph.pbs.importances.mean.sig$adj.imp)), length.out=floor(100/2))),
                     border_color = "white") %>%
  ggsave(
    filename="../../ml_plots/oars_ml_plot_heatmap.pdf",
    width=6,
    height=5,
    units="in",
    device=cairo_pdf)




bottom.1 = lasso_coefs_plot+
  lasso_scores_plot + patchwork::inset_element(lasso_scores_wilcox_plot, 0.65, 0.05, 1, 0.35)+
  patchwork::plot_layout(widths=c(1,1.5))

bottom.3 = 
  patchwork::free(ml.rf.regression.oob.oars.ph.pbs.fs.cor.plot, type="label")+
  patchwork::free(ml.fs.loocv.regression.results.ph.pbs.auc.mean.plot, type="label")+
  patchwork::free(ml.rf.regression.oob.oars.ph.pbs.fs.auc.plot, type="label")+
  patchwork::free(top.model.preds.plot, type="label")+
  patchwork::free(ml.rf.oob.auc.ph.pbs.fs.roc.plot, type="label")+
  patchwork::free(top.model.preds.cor.plot, type="label")+
  patchwork::plot_layout(nrow=3)
bottom.3



bottom.4 = ml.sample.size.curve.linear.ldo+
  ml.sample.size.curve.lm.values.plot+
  patchwork::free(ml.sample.size.curve.mse.scores.plot, type="label")+
  patchwork::plot_layout(heights=c(1,1,1))



cowplot::plot_grid(bottom.1, 
                   bottom.3, 
                   bottom.4,
                   nrow=1, rel_widths=c(3,2,1)) %>%
  ggsave(
    filename="./ml_plots/oars_ml_plot_bottom.pdf",
    width=33,
    height=10,
    units="in",
    device=cairo_pdf)

## 

r.callidus.plot %>%
  ggsave(
    filename="../../ml_plots/oars_ml_plot_taxa.pdf",
    width=18,
    height=3,
    units="in",
    device=cairo_pdf)


(ml.sample.size.curve.linear/
    ml.sample.size.curve.nonlinear)%>%
  ggsave(
    filename="../../ml_plots/oars_ml_plot_sample_size.pdf",
    width=18,
    height=5,
    units="in",
    device=cairo_pdf)



# >> SUPPLEMENTS ----------------------------------------------------------

(ml.sample.size.curve.linear.mse/
  (ml.sample.size.curve.mse.penalty.plot+
  ml.sample.size.curve.mse.penalty.plot.2+ patchwork::plot_layout(widths=c(1,1)))+
  patchwork::plot_layout(nrow=2, heights=c(1,3))) %>%
  ggsave(
    filename="../../ml_plots/oars_ml_supp_calibration.pdf",
    width=14,
    height=10,
    units="in",
    device=cairo_pdf)


# >> Graphical Abstract ---------------------------------------------------

# use ROC curves

ml.rf.oob.auc.ph.pbs.fs.roc.ga = ggplot(ml.rf.oob.auc.ph.pbs.fs.roc %>%
         subset(RS_Name != "Versafibe1490")%>%
         mutate(RS_Name = factor(RS_Name, levels=rs.names)),
       aes(x=FPR, y=TPR))+
  geom_line(aes(color=RS_Name, linetype=RS_Name), linewidth = 0.75)+
  scale_color_manual(values=labelcolors.rs)+
  scale_linetype_manual(values=c(4,8,7,6,5,1,3,2,9))+
  geom_abline(slope=1, linetype=2, linewidth=0.2)+
  theme_classic(base_size = 10)+
  theme(strip.text = element_text(size=12),
        #legend.position = c(0.85, 0.35),
        legend.position = "none",
        legend.text = element_text(size=7),
        legend.title=element_blank(),
        legend.key.spacing.y = unit(-0.2, "cm"))+
  facet_wrap(~"Random Forest")+
  labs(x="False Positive Rate", y="True Positive Rate") 

ml.rf.oob.auc.ph.pbs.fs.roc.ga%>%
  ggsave(
    filename="./ml_plots/oars_ml_graphical_abstract_roc.pdf",
    width=2.5,
    height=2.3,
    units="in",
    device=cairo_pdf)

# . -----------------------------------------------------------------------


# // defunct --------------------------------------------------------------


# :: Re-Scale Preds ----------------------------------------------------------

# What if we rescale the predictions to the min/max of true?

ml.rf.regression.oob.oars.ph.pbs.scaled = 
  do.call(rbind, lapply(1:15, function(seed){
    do.call(rbind, lapply(rs.names, function(rs){
      data.subset = ml.rf.regression.oob.oars.ph.pbs %>%
        subset(RS_Name == rs & iter == seed)
      
      # outlier predictions skew results; take percentiles
      pred_low  <- quantile(data.subset$pred.delta, probs = 0.01,  na.rm = TRUE)
      pred_high  <- quantile(data.subset$pred.delta, probs = 0.99,  na.rm = TRUE)
      meas_low  <- min(data.subset$true.delta)
      meas_high  <- max((data.subset$true.delta))
      
      data.subset$scaled.pred = scales::rescale(
        data.subset$pred.delta,
        from = c(pred_low, pred_high),     # target range
        to = c(meas_low, meas_high))        # source range
      data.subset
    }))}))
# use SCALED feature selected
ml.rf.regression.oob.oars.ph.pbs.scaled.optimal = ml.rf.regression.oob.oars.ph.pbs.scaled %>%
  group_by(HM, iter) %>%
  slice_min(scaled.pred) %>%
  mutate(is.strong = ifelse(true.delta < -1.27, "strong", "not_strong"))
ml.rf.regression.oob.oars.ph.pbs.scaled.optimal$RS_Name %>% table() %>% as.data.frame() %>% mutate(perc = Freq / (117*15))
ml.rf.regression.oob.oars.ph.pbs.scaled.optimal$is.strong %>% table()%>% as.data.frame() %>% mutate(perc = Freq / (117*15))
# 57%


# Overall correlation: 0.57
cor.test(ml.rf.regression.oob.oars.ph.pbs.scaled$scaled.pred,
         ml.rf.regression.oob.oars.ph.pbs.scaled$true.delta, method="spearman")
ggplot(ml.rf.regression.oob.oars.ph.pbs.scaled %>%
         mutate(RS_Name = factor(RS_Name, levels=rs.names)),
       aes(x=true.delta, y=scaled.pred))+
  geom_point(shape=21, aes(fill=RS_Name))+
  geom_smooth(method="lm", color="black")+
  geom_hline(yintercept=-1.27, linetype=2)+
  geom_vline(xintercept=-1.27, linetype=2)+
  ggpubr::stat_cor(method="spearman")+
  scale_fill_manual(values=labelcolors$cols)+
  theme_classic()+
  labs(title = "FS + PBS",
       x="Measured Δ pH",
       y="Predicted Δ pH")

# RS-specific correlations
do.call(rbind, lapply(rs.names, function(x){data.frame(RS_Name = x, cor = cor.test(subset(ml.rf.regression.oob.oars.ph.pbs.scaled, RS_Name == x)$scaled.pred,
                                                                                   subset(ml.rf.regression.oob.oars.ph.pbs.scaled, RS_Name == x)$true.delta, method="spearman")$estimate)}))

ggplot(ml.rf.regression.oob.oars.ph.pbs.scaled %>%
         subset(RS_Name != "PBS") %>%
         mutate(RS_Name = factor(RS_Name, levels=rs.names)),
       aes(x=true.delta, y=scaled.pred))+
  geom_point(shape=21, aes(fill=RS_Name))+
  geom_smooth(method="lm", color="black")+
  geom_hline(yintercept=-1.27, linetype=2)+
  geom_vline(xintercept=-1.27, linetype=2)+
  ggpubr::stat_cor(method="spearman")+
  scale_fill_manual(values=labelcolors$cols)+
  theme_classic()+
  facet_wrap(~RS_Name)+
  labs(x="Measured Δ pH",
       y="Predicted Δ pH")



# X. PRACTICAL UTILITY -------------------------------------------------------------


# :: Evaluate Predictors ------------------------------------------------------

# Select an optimal RS (based on pH) and evaluate rate of Strong response

# use raw real data
ml.rf.regression.oob.oars.ph.optimal = ml.rf.regression.oob.oars.ph %>%
  group_by(HM, iter) %>%
  slice_min(pred) %>%
  mutate(is.strong = ifelse(true < -1.27, "strong", "not_strong"))
ml.rf.regression.oob.oars.ph.optimal$RS_Name %>% table() %>% as.data.frame() %>% mutate(perc = Freq / (117*15))
ml.rf.regression.oob.oars.ph.optimal$is.strong %>% table()%>% as.data.frame() %>% mutate(perc = Freq / (117*15))
# 68.3%

# use pbs real data
ml.rf.regression.oob.oars.ph.pbs.optimal = ml.rf.regression.oob.oars.ph.pbs %>%
  group_by(HM, iter) %>%
  slice_min(pred) %>%
  mutate(is.strong = ifelse(true.delta < -1.27, "strong", "not_strong"))
ml.rf.regression.oob.oars.ph.pbs.optimal$RS_Name %>% table() %>% as.data.frame() %>% mutate(perc = Freq / (117*15))
ml.rf.regression.oob.oars.ph.pbs.optimal$is.strong %>% table()%>% as.data.frame() %>% mutate(perc = Freq / (117*15))
# 68.3%

# use feature selected data
ml.rf.regression.oob.oars.ph.fs.optimal = ml.rf.regression.oob.oars.ph.fs %>%
  group_by(HM, iter) %>%
  slice_min(pred) %>%
  mutate(is.strong = ifelse(true < -1.27, "strong", "not_strong"))
ml.rf.regression.oob.oars.ph.fs.optimal$RS_Name %>% table() %>% as.data.frame() %>% mutate(perc = Freq / (117*15))
ml.rf.regression.oob.oars.ph.fs.optimal$is.strong %>% table()%>% as.data.frame() %>% mutate(perc = Freq / (117*15))
# 67.5%

# use SCALED feature selected
ml.rf.regression.oob.oars.ph.pbs.scaled.optimal = ml.rf.regression.oob.oars.ph.pbs.scaled %>%
  group_by(HM, iter) %>%
  slice_min(scaled.pred) %>%
  mutate(is.strong = ifelse(true.delta < -1.27, "strong", "not_strong"))
ml.rf.regression.oob.oars.ph.pbs.scaled.optimal$RS_Name %>% table() %>% as.data.frame() %>% mutate(perc = Freq / (117*15))
ml.rf.regression.oob.oars.ph.pbs.scaled.optimal$is.strong %>% table()%>% as.data.frame() %>% mutate(perc = Freq / (117*15))
# 59%

# use shadow data
ml.rf.regression.oob.oars.ph.fs.optimal.null = ml.rf.regression.oob.oars.ph.fs %>%
  group_by(HM, iter) %>%
  slice_min(shadow) %>%
  mutate(is.strong = ifelse(true < -1.27, "strong", "not_strong"))
ml.rf.regression.oob.oars.ph.fs.optimal.null$RS_Name %>% table() %>% as.data.frame() %>% mutate(perc = Freq / (117*15))
ml.rf.regression.oob.oars.ph.fs.optimal.null$is.strong %>% table()%>% as.data.frame() %>% mutate(perc = Freq / (117*15))
# 67%

# compare to LDO by default
ml.rf.regression.oob.oars.ph.fs.optimal.ldo = ml.rf.regression.oob.oars.ph.fs %>%
  group_by(HM, iter) %>%
  subset(RS_Name == "LetsDoOrganic") %>%
  mutate(is.strong = ifelse(true < -1.27, "strong", "not_strong"))
ml.rf.regression.oob.oars.ph.fs.optimal.ldo$RS_Name %>% table() %>% as.data.frame() %>% mutate(perc = Freq / (117*15))
ml.rf.regression.oob.oars.ph.fs.optimal.ldo$is.strong %>% table()%>% as.data.frame() %>% mutate(perc = Freq / (117*15))
# 68.4%

# LDO = 68.4% by default
# NULL = 67% using nonsense microbiome data
# MLFS = 67.5% using real data
# ML = 68%

# are these significant (repeat calculations for 15 iterations)

ml.rf.regression.selection.accuracy = do.call(rbind, lapply(1:15, function(seed){
  do.call(rbind, lapply(c("LetsDoOrganic", "Null","Null (Feature Select)", "ML", "ML (Feature Select)"), function(type){
    if(type == "LetsDoOrganic"){
      # compare to LDO by default
      data.to.evaluate = ml.rf.regression.oob.oars.ph %>%
        subset(iter == seed) %>%
        group_by(HM) %>%
        subset(RS_Name == "LetsDoOrganic") %>%
        mutate(is.strong = ifelse(true < -1.27, "strong", "not_strong"))
    }
    if(type == "Null"){
      data.to.evaluate = ml.rf.regression.oob.oars.ph %>%
        subset(iter == seed) %>%
        group_by(HM) %>%
        slice_min(shadow) %>%
        mutate(is.strong = ifelse(true < -1.27, "strong", "not_strong"))
    }
    if(type == "Null (Feature Select)"){
      data.to.evaluate = ml.rf.regression.oob.oars.ph.fs %>%
        subset(iter == seed) %>%
        group_by(HM) %>%
        slice_min(shadow) %>%
        mutate(is.strong = ifelse(true < -1.27, "strong", "not_strong"))
    }
    if(type == "ML"){
      data.to.evaluate = ml.rf.regression.oob.oars.ph %>%
        subset(iter == seed) %>%
        group_by(HM) %>%
        slice_min(pred) %>%
        mutate(is.strong = ifelse(true < -1.27, "strong", "not_strong"))
    }
    if(type == "ML (Feature Select)"){
      data.to.evaluate = ml.rf.regression.oob.oars.ph.fs %>%
        subset(iter == seed) %>%
        group_by(HM) %>%
        slice_min(pred) %>%
        mutate(is.strong = ifelse(true < -1.27, "strong", "not_strong"))
    }
    
    data.output = data.to.evaluate$is.strong %>% table()%>% as.data.frame() %>% mutate(perc = Freq / (117))
    colnames(data.output)[1] = c("strong")
    data.frame(accuracy = subset(data.output, strong == "strong")$perc,
               type = type,
               iter = seed)
  }))})) %>%
  mutate(type = factor(type, levels=c("LetsDoOrganic",
                                      "Null", "ML",
                                      "Null (Feature Select)","ML (Feature Select)")))

# stats?
kruskal.test(accuracy ~ type, ml.rf.regression.selection.accuracy)
dunn.test::dunn.test(ml.rf.regression.selection.accuracy$accuracy, ml.rf.regression.selection.accuracy$type)

ml.rf.regression.selection.accuracy %>%
  group_by(type) %>%
  mutate(mean.acc = mean(accuracy)) %>%
  dplyr::select(type, mean.acc) %>% distinct()


ggplot(ml.rf.regression.selection.accuracy,
       aes(x=type, y=accuracy))+
  geom_violin()+
  ggbeeswarm::geom_beeswarm()+
  theme_classic()+
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  labs(x=NULL, y="Prediction Accuracy")

## evaluate RS selections

ml.rf.regression.selections = do.call(rbind, lapply(1:15, function(seed){
  do.call(rbind, lapply(c("LetsDoOrganic", "Null","Null (Feature Select)", "ML", "ML (Feature Select)"), function(type){
    if(type == "LetsDoOrganic"){
      # compare to LDO by default
      data.to.evaluate = ml.rf.regression.oob.oars.ph %>%
        subset(iter == seed) %>%
        group_by(HM) %>%
        subset(RS_Name == "LetsDoOrganic") %>%
        mutate(is.strong = ifelse(true < -1.27, "strong", "not_strong"))
    }
    if(type == "Null"){
      data.to.evaluate = ml.rf.regression.oob.oars.ph %>%
        subset(iter == seed) %>%
        group_by(HM) %>%
        slice_min(shadow) %>%
        mutate(is.strong = ifelse(true < -1.27, "strong", "not_strong"))
    }
    if(type == "Null (Feature Select)"){
      data.to.evaluate = ml.rf.regression.oob.oars.ph.fs %>%
        subset(iter == seed) %>%
        group_by(HM) %>%
        slice_min(shadow) %>%
        mutate(is.strong = ifelse(true < -1.27, "strong", "not_strong"))
    }
    if(type == "ML"){
      data.to.evaluate = ml.rf.regression.oob.oars.ph %>%
        subset(iter == seed) %>%
        group_by(HM) %>%
        slice_min(pred) %>%
        mutate(is.strong = ifelse(true < -1.27, "strong", "not_strong"))
    }
    if(type == "ML (Feature Select)"){
      data.to.evaluate = ml.rf.regression.oob.oars.ph.fs %>%
        subset(iter == seed) %>%
        group_by(HM) %>%
        slice_min(pred) %>%
        mutate(is.strong = ifelse(true < -1.27, "strong", "not_strong"))
    }
    
    data.output = data.to.evaluate$RS_Name %>% table()%>% as.data.frame() %>% mutate(perc = Freq / (117)) %>% as.data.frame()
    colnames(data.output)[1] = c("RS_Name")
    data.output = data.output %>%
      mutate(type = type,
             iter = seed)
  }))}))%>%
  mutate(type = factor(type, levels=c("LetsDoOrganic",
                                      "Null", "ML",
                                      "Null (Feature Select)","ML (Feature Select)")))

ml.rf.regression.selections = ml.rf.regression.selections %>%
  mutate(RS_Name = factor(RS_Name, levels=rs.names))

ggplot(ml.rf.regression.selections,
       aes(x=type, y=Freq/15/117 * 100))+
  geom_bar(stat="identity", position="stack", aes(fill=RS_Name))+
  scale_fill_manual(values=labelcolors.rs)+
  theme_classic()+
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  labs(x=NULL, y="Frequency Selected")


# :: BEST RAW pH ---------------------------------------------------------------


ml.rf.regression.raw.selection.accuracy = 
  do.call(rbind, lapply(c(-1.27), function(thresh){
    do.call(rbind, lapply(1:15, function(seed){
      do.call(rbind, lapply(c("LetsDoOrganic", "Null", "ML"), function(type){
        if(type == "LetsDoOrganic"){
          # compare to LDO by default
          data.to.evaluate = ml.rf.regression.oob.oars.raw.ph %>%
            subset(iter == seed) %>%
            group_by(HM) %>%
            subset(RS_Name == "LetsDoOrganic") %>%
            mutate(is.strong = ifelse(true.delta < thresh, "strong", "not_strong"))
        }
        if(type == "Null"){
          data.to.evaluate = ml.rf.regression.oob.oars.raw.ph %>%
            subset(iter == seed) %>%
            group_by(HM) %>%
            slice_min(shadow.delta) %>%
            mutate(is.strong = ifelse(true.delta < thresh, "strong", "not_strong"))
        }
        
        if(type == "ML"){
          data.to.evaluate = ml.rf.regression.oob.oars.raw.ph %>%
            subset(iter == seed) %>%
            group_by(HM) %>%
            slice_min(pred.delta) %>%
            mutate(is.strong = ifelse(true.delta < thresh, "strong", "not_strong"))
        }
        
        data.output = data.to.evaluate$is.strong %>% table()%>% as.data.frame() %>% mutate(perc = Freq / (117))
        colnames(data.output)[1] = c("strong")
        data.frame(accuracy = subset(data.output, strong == "strong")$perc,
                   type = type,
                   iter = seed,
                   thresh = thresh)
      }))}))})) %>%
  mutate(type = factor(type, levels=c("LetsDoOrganic",
                                      "Null", "ML")))

# stats?
kruskal.test(accuracy ~ type, ml.rf.regression.raw.selection.accuracy)
dunn.test::dunn.test(ml.rf.regression.raw.selection.accuracy$accuracy, 
                     ml.rf.regression.raw.selection.accuracy$type)

ml.rf.regression.raw.selection.accuracy %>%
  group_by(type, thresh) %>%
  mutate(mean.acc = mean(accuracy)) %>%
  dplyr::select(type, mean.acc, thresh) %>% distinct()


ggplot(ml.rf.regression.raw.selection.accuracy,
       aes(x=type, y=accuracy))+
  geom_boxplot()+
  ggbeeswarm::geom_beeswarm()+
  theme_classic()+
  #scale_y_continuous(limits=c(0.6,0.75))+
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  facet_wrap(~thresh, scales="free")+
  labs(x=NULL, y="Prediction Accuracy")



# :: Supp. Multivariate ---------------------------------------------------------------

# Results are comparable to univariate, but generally inferior

t1 = Sys.time()
ml.multivariate.raw.rf.results <- 
  do.call(rbind, lapply(c("delta_pH"), function(y){
    # y = "delta_pH"
    # x = rs.names[6]
    
    ml.target.data.subset = ph.data
    
    # apply feature selection
    ml.input.data = ml.input.data[,colnames(ml.input.data) %in% rev(ml_rs_features)[c(1:n.features)]]
    
    # DONT scale responses
    ml.target.data.subset = ml.target.data.subset %>% group_by(RS_Name) %>% 
      mutate(value = (med.ph))
    # cast responses
    ml.target.data.subset = reshape2::acast(ml.target.data.subset, HM ~ RS_Name, value.var = "value")
    # repeat 15 iterations
    do.call(rbind, lapply(1:15, function(z){
      print(z)
      # z = 1
      full.data = merge(data.frame(ml.input.data), ml.target.data.subset, by="row.names")
      rownames(full.data) <- full.data$HM
      full.data$Row.names <- NULL
      full.data$HM <- NULL
      full.data
      
      ## Main RF
      set.seed(z)
      ml.mtl.subset = randomForestSRC::rfsrc(get.mv.formula(rs.names.pbs),
                                             data.frame(ml.target.data.subset, ml.input.data),
                                             importance=TRUE, nsplit = 10, splitrule = "mahalanobis")
      # need to loop through to extract
      mtl.output = do.call(rbind, lapply(1:10, function (x){
        mtl.output.rs = data.frame(RS_Name = (rs.names.pbs)[x],
                                   ml.mtl.subset$regrOutput[[x]][2]) # 2 = out of bag
      }))
      mtl.output$HM = rep(unique(ml.target.data$HM), times=10)
      mtl.output = reshape2::acast(mtl.output, HM ~ RS_Name, value.var = "predicted.oob" )
      # melt again, now that HMs are added
      mtl.output = reshape2::melt(mtl.output)
      colnames(mtl.output) <- c("HM", "RS_Name", "predicted")
      mtl.output = merge(mtl.output, ph.data, by=c("HM", "RS_Name"))
      
      ## Shadow RF
      set.seed(z)
      # shuffle values
      ml.target.data.subset.shuffled = do.call(cbind, lapply(rs.names.pbs, function(x){
        value.vector = ml.target.data.subset[,colnames(ml.target.data.subset) == x]
        value.vector = sample(value.vector, size=length(value.vector))
        value.df = data.frame(value.vector)
        value.df
      }))
      colnames(ml.target.data.subset.shuffled) = rs.names.pbs
      ml.target.data.subset.shuffled = ml.target.data.subset.shuffled[,rs.names.pbs]
      # run model
      ml.mtl.subset.shadow = randomForestSRC::rfsrc(get.mv.formula(rs.names.pbs),
                                                    data.frame(ml.target.data.subset.shuffled, ml.input.data),
                                                    importance=TRUE, nsplit = 10, splitrule = "mahalanobis")
      # need to loop through to extract
      mtl.output.shadow = do.call(rbind, lapply(1:10, function (x){
        mtl.output.rs = data.frame(RS_Name = (rs.names.pbs)[x],
                                   ml.mtl.subset.shadow$regrOutput[[x]][2]) # 2 = out of bag
      }))
      mtl.output.shadow$HM = rep(unique(ml.target.data$HM), times=10)
      mtl.output.shadow = reshape2::acast(mtl.output.shadow, HM ~ RS_Name, value.var = "predicted.oob" )
      # melt again, now that HMs are added
      mtl.output.shadow = reshape2::melt(mtl.output.shadow)
      colnames(mtl.output.shadow) <- c("HM", "RS_Name", "shadow")
      
      # merge
      mtl.output = merge(mtl.output, mtl.output.shadow, by=c("HM", "RS_Name"))
      
      mtl.output$iter =  z
      return(mtl.output)
    }))
  }))

t2 <- Sys.time()
t2 - t1 # ~ 5 min


# now predict change in pH using predictions
ml.multivariate.raw.rf.results = ml.multivariate.raw.rf.results %>%
  group_by(HM, iter) %>%
  mutate(true.delta = med.ph - med.ph[RS_Name == "PBS"],
         pred.delta = predicted - predicted[RS_Name == "PBS"],
         shadow.delta = shadow - shadow[RS_Name == "PBS"])


# Overall correlation
ggplot(ml.multivariate.raw.rf.results %>%
         subset(RS_Name != "PBS")%>%
         mutate(RS_Name = factor(RS_Name, levels=rs.names)),
       aes(x=true.delta, y=pred.delta))+
  geom_point(shape=21, aes(fill=RS_Name))+
  geom_smooth(method="lm", color="black")+
  geom_hline(yintercept=-1.27, linetype=2)+
  geom_vline(xintercept=-1.27, linetype=2)+
  ggpubr::stat_cor(method="spearman")+
  scale_fill_manual(values=labelcolors$cols)+
  theme_classic()+
  labs(title = "Overall",
       x="Measured Δ pH",
       y="Predicted Δ pH")

# RS-specific correlations
ggplot(ml.multivariate.raw.rf.results %>%
         subset(RS_Name != "PBS") %>%
         mutate(RS_Name = factor(RS_Name, levels=rs.names)),
       aes(x=true.delta, y=pred.delta))+
  geom_point(shape=21, aes(fill=RS_Name))+
  geom_smooth(method="lm", color="black")+
  geom_hline(yintercept=-1.27, linetype=2)+
  geom_vline(xintercept=-1.27, linetype=2)+
  ggpubr::stat_cor(method="spearman")+
  scale_fill_manual(values=labelcolors$cols)+
  theme_classic()+
  facet_wrap(~RS_Name)+
  labs(title = "RS-Specific",
       x="Measured Δ pH",
       y="Predicted Δ pH")
# not good!

# calculate AUC
ml.multivariate.raw.rf.auc = ml.multivariate.raw.rf.results %>%
  subset(RS_Name != "PBS") %>%
  #subset(RS_Name %in% c("LetsDoOrganic", "Versafibe1490")) %>%
  mutate(true.response = ifelse(true.delta < -1.27, "strong", "weak")) %>%
  group_by(RS_Name, iter) %>%
  mutate(cor = cor.test(true.delta, pred.delta)$estimate,
         cor.shadow = cor.test(true.delta, shadow.delta)$estimate) %>%
  mutate(auc = pROC::auc(as.factor(true.response), pred.delta, 
                         direction = ">", levels=c("weak", "strong")) %>% mean()) %>%
  mutate(auc.shadow = pROC::auc(as.factor(true.response), shadow.delta, 
                                direction = ">", levels=c("weak", "strong")) %>% mean()) %>%
  dplyr::select(RS_Name, iter, auc, auc.shadow, cor, cor.shadow) %>% distinct() %>%
  mutate(mean.auc = mean(auc),
         mean.auc.shadow = mean(auc.shadow),
         mean.cor = mean(cor),
         mean.cor.shadow = mean(cor.shadow)) %>%
  dplyr::select(RS_Name, mean.auc, mean.auc.shadow, mean.cor, mean.cor.shadow) %>% distinct()  %>%
  mutate(RS_Name = factor(RS_Name, levels=rs.names))

ml.multivariate.raw.rf.auc.p = do.call(rbind, lapply(rs.names, function(rs){
  data.subset = subset(ml.multivariate.raw.rf.auc, RS_Name == rs)
  data.frame(RS_Name = rs,
             pval = wilcox.test(data.subset$mean.auc,
                                data.subset$mean.auc.shadow)$p.value)
})) %>%
  mutate(padj = p.adjust(pval, method = "bonferroni")) %>%
  mutate(sig = ifelse(padj < 0.05, "*", "")) %>%
  mutate(RS_Name = factor(RS_Name, levels=rs.names))

ml.multivariate.raw.rf.auc.mean = ml.multivariate.raw.rf.auc %>%
  mutate(RS_Name = factor(RS_Name, levels=rs.names))%>%
  group_by(RS_Name) %>%
  mutate(mean.auc = mean(mean.auc),
         mean.shadow = mean(mean.auc.shadow)) %>%
  dplyr::select(RS_Name, mean.auc, mean.shadow) %>% distinct() 

ml.multivariate.raw.rf.auc.mean = merge(ml.multivariate.raw.rf.auc.mean,
                                        ml.multivariate.raw.rf.auc.p, by="RS_Name") %>%
  mutate(RS_Name = factor(RS_Name, levels=rs.names))


ml.multivariate.raw.rf.auc.plot <- ggplot()+
  # plot shadow
  geom_violin(data=ml.multivariate.raw.rf.auc, aes(x=RS_Name, y=mean.auc.shadow, fill=RS_Name),fill="grey", alpha=0.4, color=NA)+
  geom_point(data=ml.multivariate.raw.rf.auc.mean, aes(x=RS_Name, y=mean.shadow), color="grey", size=3, alpha=0.8)+
  # plot real
  geom_violin(data=ml.multivariate.raw.rf.auc, aes(x=RS_Name, y=mean.auc, fill=RS_Name), color="black")+
  geom_point(data=ml.multivariate.raw.rf.auc.mean, aes(x=RS_Name, y=mean.auc), shape=21, fill="white", size=2.5)+
  # plot median value
  geom_text(data=ml.multivariate.raw.rf.auc.mean, aes(x=RS_Name, y=(mean.auc+0.05), label=round(mean.auc, digits=2)), nudge_y=0.04, size=3)+
  # add sig vs shadow
  geom_text(data=ml.multivariate.raw.rf.auc.mean, aes(x=RS_Name, y=Inf, label=sig), vjust=1.5, size=6)+
  geom_hline(yintercept=0.5, linetype=2, alpha=0.5)+
  scale_y_continuous(breaks=seq(from=0.1, to=0.9, by=0.1), limits=c(0.1,0.9))+
  scale_fill_manual(values=labelcolors$cols[c(1:9)])+
  scale_color_manual(values=labelcolors$cols[c(1:9)])+
  theme_classic()+theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1),
                        legend.position="none",
                        panel.grid.major.y = element_line(color="grey", linewidth=0.2),
                        legend.title = element_text(hjust=0.5),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  facet_wrap(~"Random Forest Out-of-bag AUC")+
  labs(x=NULL, y="AUC", title = "Default")
ml.multivariate.raw.rf.auc.plot


ml.multivariate.raw.selection.accuracy = 
  do.call(rbind, lapply(c(-1.27), function(thresh){
    do.call(rbind, lapply(1:15, function(seed){
      do.call(rbind, lapply(c("LetsDoOrganic", "Null", "ML"), function(type){
        if(type == "LetsDoOrganic"){
          # compare to LDO by default
          data.to.evaluate = ml.multivariate.raw.rf.results %>%
            subset(iter == seed) %>%
            group_by(HM) %>%
            subset(RS_Name == "LetsDoOrganic") %>%
            mutate(is.strong = ifelse(true.delta < thresh, "strong", "not_strong"))
        }
        if(type == "Null"){
          data.to.evaluate = ml.multivariate.raw.rf.results %>%
            subset(iter == seed) %>%
            group_by(HM) %>%
            slice_min(shadow.delta) %>%
            mutate(is.strong = ifelse(true.delta < thresh, "strong", "not_strong"))
        }
        
        if(type == "ML"){
          data.to.evaluate = ml.multivariate.raw.rf.results %>%
            subset(iter == seed) %>%
            group_by(HM) %>%
            slice_min(pred.delta) %>%
            mutate(is.strong = ifelse(true.delta < thresh, "strong", "not_strong"))
        }
        
        data.output = data.to.evaluate$is.strong %>% table()%>% as.data.frame() %>% mutate(perc = Freq / (117))
        colnames(data.output)[1] = c("strong")
        data.frame(accuracy = subset(data.output, strong == "strong")$perc,
                   type = type,
                   iter = seed,
                   thresh = thresh)
      }))}))})) %>%
  mutate(type = factor(type, levels=c("LetsDoOrganic",
                                      "Null", "ML")))

# stats?
kruskal.test(accuracy ~ type, ml.multivariate.raw.selection.accuracy)
dunn.test::dunn.test(ml.multivariate.raw.selection.accuracy$accuracy, 
                     ml.multivariate.raw.selection.accuracy$type)

ml.multivariate.raw.selection.accuracy %>%
  group_by(type, thresh) %>%
  mutate(mean.acc = mean(accuracy)) %>%
  dplyr::select(type, mean.acc, thresh) %>% distinct()


ggplot(ml.multivariate.raw.selection.accuracy,
       aes(x=type, y=accuracy))+
  geom_boxplot()+
  ggbeeswarm::geom_beeswarm()+
  theme_classic()+
  #scale_y_continuous(limits=c(0.6,0.75))+
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  facet_wrap(~thresh, scales="free")+
  labs(x=NULL, y="Prediction Accuracy")





# >> LOCO -----------------------------------------------------------------


# :: Maaslin2 -------------------------------------------------------------

# Apply multiple linear regression with random effects to identify bacteria significantly associated with RS intake.


# run maaslin2
rs_maaslin_output = Maaslin2::Maaslin2(
  input_data = unified.df.filt.badj$feature_abd_adj,
  input_metadata = unified.meta,
  normalization = "NONE",
  transform = "LOG",
  fixed_effects = c("timepoint"),
  random_effects = c("group","subject"),
  output = "maaslin_output", 
  plot_heatmap = F,
  plot_scatter = F)
rs_maaslin_df = data.frame(rs_maaslin_output$results)
rs_maaslin_df$study = "all studies"


# and loop through each study
rs_maaslin_output_loop = 
  do.call(rbind, lapply(rs.groups, function(x){
    print(x)
    data.subset = subset(unified.meta, group == x)
    rs_maaslin_output = Maaslin2::Maaslin2(
      input_data = unified.df.filt,
      input_metadata = data.subset,
      normalization = "NONE",
      transform = "LOG",
      fixed_effects = c("timepoint"), # so adjusted LogFC are compatible with group-level
      random_effects = c("subject"),
      output = "maaslin_output", 
      plot_heatmap = F,
      plot_scatter = F)
    rs_maaslin_df = data.frame(rs_maaslin_output$results)
    rs_maaslin_df$study = x
    rs_maaslin_df
  }))



# :: evaluate -------------------------------------------------------------


# Now we can plot heatmaps like Chen in Fig 4 and trend plots like in Fig 2 a.


rs_maaslin_df = rs_maaslin_df %>%
  mutate(feature = gsub("g__CAG.", "g__CAG-", feature))
rs_maaslin_output_loop = rs_maaslin_output_loop %>%
  mutate(feature = gsub("g__CAG.", "g__CAG-", feature))

# combine maaslin results
rs_maaslin_all = rbind(rs_maaslin_df,
                       rs_maaslin_output_loop) %>% data.frame()

# select sig taxa (at all study level)
rs_maaslin_all = subset(rs_maaslin_all, feature %in%
                          subset(rs_maaslin_df, qval < 0.05)$feature)
unique(rs_maaslin_all$feature) # 31 sig

unique(rs_maaslin_df$feature) %>% length() # 270 eligible
unique(subset(rs_maaslin_all, qval < 0.05)$feature) %>% length() # 31 sig

loco.names = data.frame(study = unique(rs_maaslin_all$study),
                        clean = factor(c("All studies", "Deehan (Cross-linked corn)", 
                                         "Deehan (Cross-linked potato)", "Deehan (Cross-linked tapioca)",
                                         "DeMartino (Potato)","Flowers (Potato)", 
                                         "Hanes (Potato-Banana-ApplePectin_A)","Hanes (Potato-Banana-ApplePectin_B)",
                                         "Hanes (Potato-Banana-ApplePectin_C)", "Hanes (Potato)",
                                         "Hughes (Wheat)", "Maier (Corn_A)", "Maier (Corn_B)",
                                         "Upadhyaya (Cross-linked wheat)", "Venkataraman (Potato)"),
                                       levels=c("All studies", "Deehan (Cross-linked corn)", 
                                                "Deehan (Cross-linked potato)", "Deehan (Cross-linked tapioca)",
                                                "DeMartino (Potato)","Flowers (Potato)", 
                                                "Hanes (Potato-Banana-ApplePectin_A)","Hanes (Potato-Banana-ApplePectin_B)",
                                                "Hanes (Potato-Banana-ApplePectin_C)", "Hanes (Potato)",
                                                "Hughes (Wheat)", "Maier (Corn_A)", "Maier (Corn_B)",
                                                "Upadhyaya (Cross-linked wheat)", "Venkataraman (Potato)")))
# add names
rs_maaslin_all = merge(rs_maaslin_all,
                       loco.names, by="study")

# order taxa
rs_maaslin_all$feature = factor(rs_maaslin_all$feature, arrange(subset(rs_maaslin_df, qval < 0.05), coef)$feature)


# :: :: PLOT: Maaslin2 ----------------------------------------------------


rs_maaslin_all_plot = ggplot(rs_maaslin_all,
                             aes(x=clean, y=feature))+
  geom_tile(aes(fill=coef), color="white")+
  geom_text(aes(label=ifelse(qval < 0.05, "*", NA)),
            vjust=0.75, size=5, color="black")+
  theme_classic()+
  scale_fill_gradient2(low = ("blue"),
                       high = ("red"))+
  geom_vline(xintercept=1.5,color="black")+
  theme(#axis.text.x=element_text(angle=45, hjust=1, vjust=1),
    axis.text.x = element_blank(),
    strip.text = element_text(size=10))+
  facet_wrap(~"Random Effects Regression")+
  labs(x=NULL, y=NULL, fill=expression(Log[2]*FC))
rs_maaslin_all_plot


# :: :: Response ----------------------------------------------------------

# Response = Log ratio of RS-Enriched / RS-Depleted

all.tax.table = readRDS("../../chen_validation/processed_files/all.tax.table.gg2.Rds")

delta.rs.response = function(otu_table = otu_table){
  # create list of taxa that include butyrogens
  up.taxa = subset(rs_maaslin_all, clean == "All studies" & coef > 0)$feature
  down.taxa =  subset(rs_maaslin_all, clean == "All studies" & coef < 0)$feature
  
  # add pseudo
  pseudo = min(otu_table[otu_table!=0])/2
  
  # subset otu_table to up taxa and down taxa
  up_table = otu_table[,colnames(otu_table) %in% up.taxa]+pseudo
  down_table = otu_table[,colnames(otu_table) %in% down.taxa]+pseudo
  
  # calculate sums of all butyrogens per sample
  up_sums = rowSums(up_table)
  down_sums = rowSums(down_table)
  
  # log ratio
  ratio = up_sums / down_sums
  
  # save as dataframe
  ratio_df = data.frame(sample = rownames(otu_table),
                        rs_ratio = ratio)
  # calculate change
  ratio_df = merge(ratio_df,
                   unified.meta, by="sample")
  
  ratio_df = ratio_df %>%
    group_by(subject) %>%
    mutate(delta = log(rs_ratio / rs_ratio[timepoint == "baseline"], base=2)) %>%
    subset(timepoint == "RS") %>%
    dplyr::select(subject, timepoint, delta) %>% distinct() %>% data.frame()
  
  return(ratio_df)
}


# apply butyrogen function once
rs_targets = delta.rs.response(otu_table = unified.df.filt)

# add study ID
rs_targets = merge(rs_targets,
                   unified.meta[,c("subject","RS_Name", "study", "group")] %>% distinct(), by="subject")

# Now we can plot these values

rs_target_rs_plot = ggplot(lasso_ratio,
                           aes(x=group, y=ratio))+
  geom_violin(aes(fill=study), alpha=0.3,
              draw_quantiles = 0.5)+
  geom_point(shape=21, aes(fill=study))+
  theme_classic()+theme(legend.position="none")+
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1),
        panel.grid.major.y = element_line(color="grey", linewidth=0.2))+
  labs(x=NULL, y=expression(Log[2]*FC~RS~Response~Ratio))
rs_target_rs_plot

rs_target_rs_notscaled = lasso_ratio %>%
  group_by(study, RS_Name) %>%
  mutate(delta = (ratio)) %>% data.frame() %>%
  # add "responder" or not, based on median
  mutate(delta.resp = ifelse(delta > median(delta), "responder", "nonresponder"))

# color based on RS type
rs_target_rs_notscaled = merge(rs_target_rs_notscaled,
                               subset(rs.types, group != "all studies")[,c("group", "cols", "RS_type")], by="group") 


rs.types.list = rs.types$cols
names(rs.types.list) = rs.types$group

loco.names.simple = data.frame(group = unique(rs_target_rs_notscaled$group),
                               clean = c("Deehan (Cross-linked corn)", "Deehan (Cross-linked potato)", "Deehan (Cross-linked tapioca)",
                                         "DeMartino (Potato)", "Flowers (Potato)", 
                                         "Hanes (Potato-Banana-ApplePectin A)", "Hanes (Potato-Banana-ApplePectin B)", "Hanes (Potato-Banana-ApplePectin C)",
                                         "Hanes (Potato)", "Hughes (Wheat)", "Maier (Corn A)", "Maier (Corn B)", "Upadhyaya (Cross-linked wheat)",
                                         "Venkataraman (Potato)"))

rs_target_rs_notscaled = merge(rs_target_rs_notscaled, loco.names.simple, by="group")
rs_target_rs_notscaled$clean


# ADD ALL STUDIES

# position beneath Maaslin2 heatmap
rs_target_rs_notscaled = rbind(rs_target_rs_notscaled,
                               rs_target_rs_notscaled %>% mutate(clean = "All studies", study = "All studies")) %>%
  as.data.frame() %>%
  mutate(clean = factor(clean, levels=c("All studies", loco.names.simple$clean)))


# :: :: PLOT: Response ----------------------------------------------------

rs_target_rs_notscaled_plot = ggplot(rs_target_rs_notscaled,
                                     aes(x=clean, y=delta))+
  geom_violin(aes(fill=study), alpha=0.3)+
  ggbeeswarm::geom_beeswarm(shape=21, aes(fill=study), cex=0.5)+
  scale_fill_manual(values=c("grey", RColorBrewer::brewer.pal(name="Set3", n=length(unique(rs_target_rs_notscaled$study)))))+
  theme_classic()+theme(legend.position="none")+
  facet_wrap(~"RS-Response Signature")+
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1),
        panel.grid.major.y = element_line(color="grey", linewidth=0.2),
        strip.text = element_text(size=10))+
  labs(x=NULL, y=expression(Log[2]*FC~to~Baseline))
rs_target_rs_notscaled_plot


# :: :: LOCO function --------------------------------------------------------

# use this as target/metadata
rs_target_rs_notscaled
# use this as ASV data
rs_ml_data = unified.df.filt[subset(unified.meta, timepoint=="baseline")$sample,]
# replace rowname with subject
rownames(rs_ml_data) = unified.meta[match(rownames(rs_ml_data), unified.meta$sample),]$subject

# log-transform predictors + pseudocount
ml_pseudo = min(rs_ml_data[rs_ml_data!=0])/2
rs_ml_data = log(rs_ml_data+ml_pseudo, base=2)

# build function
ml_loco = function(otu_table = rs_ml_data,
                   metadata = rs_target_rs_notscaled,
                   heldout = study,
                   target = "delta",
                   algo = "oob", # or oob
                   output = "auc", # AUC or importance
                   auc = TRUE, # only calculate AUC, otherwise return predictions
                   scale = T,
                   iter = iter){
  # set holdout meta data
  training_metadata = subset(metadata, group != heldout)
  heldout_metadata = subset(metadata, group == heldout)
  # set holdout ASV data
  training_data = otu_table[training_metadata$subject,] %>% as.data.frame()
  heldout_data = otu_table[heldout_metadata$subject,] %>% as.data.frame()
  # append target
  
  if(scale == T){
    # scale training data WITHOUT test data
    training_metadata$delta = scale(training_metadata$delta)
    # scale test data WITH training data
    heldout_metadata = metadata %>%
      mutate(delta = scale(delta)) %>%
      subset(group == heldout)
  }
  
  if(target %in% c("butyrogen", "but", "butyrogens", "delta")){
    training_data$subject = rownames(training_data)
    training_data = merge(training_data, 
                          training_metadata[,c("subject", "delta")], by="subject")
    training_data$subject = NULL
    # rename target column (we'll keep track later)
    colnames(training_data)[colnames(training_data)=="delta"] <- "target"
  }
  
  # build ranger
  if(algo %in% c("oob", "ranger_oob", "OOB")){
    
    set.seed(iter)
    ml_model = ranger::ranger(
      (target) ~ .,
      data.frame(training_data),
      importance = "permutation") %>% suppressWarnings() %>% suppressMessages()
    
  } 
  
  if(output %in% c("auc", "AUC")){
    # apply predictions
    ml_predictions = predict(ml_model, data.frame(heldout_data)) # %>% scale()
    if(algo %in% c("oob", "ranger_oob", "OOB")){
      ml_predictions = ml_predictions$predictions
    }
    # collate data
    
    if(target %in% c("butyrogen", "but", "butyrogens", "delta")){
      ml_output = data.frame(pred = ml_predictions,
                             true = heldout_metadata$delta,
                             class = heldout_metadata$delta.resp,
                             target = "delta")
    }
    # add other vars
    ml_output$group = heldout
    ml_output$model = algo
    ml_output$iter = iter
    # make binary pred (via sigmoid transformation)
    ml_output$pred.class = ifelse(sigmoid::sigmoid(ml_output$pred)>0.5, "responder", "nonresponder")
    
    if(auc == TRUE){
      ml_auc = pROC::auc(ml_output$class, ml_output$pred, 
                         # must specify these:
                         levels=c("responder", "nonresponder"),
                         direction = ">")
      ml_output = data.frame(target = target,
                             model = algo,
                             iter = iter,
                             group = heldout,
                             cor = cor.test(ml_output$pred, ml_output$true, method="spearman")$estimate,
                             auc = ml_auc)
      return(ml_output)
    }
  }
  
  # if not auc, then importance:
  
  if(output %in% c("importance", "importances")){
    
    ml_importances = data.frame(importance = ml_model$variable.importance) %>%
      mutate(feature = rownames(.)) %>%
      arrange(-importance) %>%
      mutate(target = target,
             model = algo,
             iter = iter,
             group = heldout)
    
    return(ml_importances)
    
  }
}

# :: :: LOCO -----------------------------------------------------------------


# check number of participants (n=208)
rs_target_rs_notscaled$subject %>% unique() %>% length()

ml.models = "oob"

# data
# rs_targets_scaled
# rs_ml_data

iters = c(1:15)
data.types = c("real", "null")
studies = rs.groups
algos = ml.models

t0 <- Sys.time()
ml_loop_output = 
  do.call(rbind, lapply(iters, function(iter){
    do.call(rbind, lapply(data.types, function(data.type){
      
      # conditionally scramble target values, per study
      if(data.type == "null"){
        set.seed(iter)
        rs_ml_input = rs_target_rs_notscaled %>% 
          group_by(study) %>%
          mutate(delta = sample(delta)) %>%
          # redefine response
          mutate(delta.resp = ifelse(delta > median(delta), "responder", "nonresponder"))
      } else {
        rs_ml_input = rs_target_rs_notscaled
      }
      
      # otherwise, proceed
      do.call(rbind, lapply(studies, function(study){
        # print status outside of mclapply applied on models
        t1 <- Sys.time()
        print(paste0(iter, " ", data.type, " ", study, " ", round(t1-t0, digits=1), sep=""))
        do.call(rbind, parallel::mclapply(algos, function(algo){
          # apply function
          ml_loco_output = ml_loco(otu_table = rs_ml_data,
                                   metadata = rs_ml_input,
                                   heldout = study,
                                   target="delta",
                                   algo = algo,
                                   auc = TRUE, # only calculate AUC, otherwise return predictions
                                   scale = T,
                                   iter = iter)
          ml_loco_output$data.type = data.type
          ml_loco_output
        }))}))}))}))
t2 = Sys.time() # 5.5 hours; 2 min for ranger oob

t2 - t0

# calculate median
ml_loop_output_mean = ml_loop_output %>% as.data.frame() %>%
  dplyr::group_by(data.type, group) %>%
  mutate(mean.auc = mean(auc)) %>%
  dplyr::select(data.type, model, target, group, mean.auc) %>% distinct() %>% data.frame()

# run stats
ml_loop_output_pvalue = 
  do.call(rbind, lapply(rs.groups, function(x){
    data.subset = subset(ml_loop_output, group == x)
    ttest.results = wilcox.test(subset(data.subset, data.type=="real")$auc,
                                subset(data.subset, data.type=="null")$auc)$p.value
    ttest.results = data.frame(group = x,
                               pval = ttest.results)
    ttest.results
  }))
ml_loop_output_pvalue
ml_loop_output_pvalue$padj = p.adjust(ml_loop_output_pvalue$pval, method="BH")
ml_loop_output_pvalue$sig = ifelse(ml_loop_output_pvalue$padj < 0.05, "*", "")

# merge
ml_loop_output_mean = merge(ml_loop_output_mean,
                            ml_loop_output_pvalue, by=c("group"))
# delete if mean real < mean null
ml_loop_output_mean = ml_loop_output_mean %>%
  group_by(group) %>%
  mutate(sig = ifelse(mean.auc > mean.auc[data.type == "null"], sig, ""))


# forest plot
ml_loop_output_forest <- ml_loop_output %>%
  group_by(data.type, model) %>%
  summarise(
    mean_auc = median(auc),
    lower_ci = quantile(auc, 0.025),
    upper_ci = quantile(auc, 0.975),
    .groups = "drop"
  ) %>% data.frame() %>% mutate(group = "all studies")
ml_loop_output_forest

# add dummies
ml_loop_output = rbind(ml_loop_output,
                       data.frame(target = "delta", model = "oob", iter = NA, group = "all studies",cor = -1,  auc = -1, data.type = "real"),
                       data.frame(target = "delta", model = "oob", iter = NA, group = "all studies",cor = -1,  auc = -1, data.type = "null")) %>% data.frame()

ml_loop_output_mean = rbind(ml_loop_output_mean,
                            data.frame(group = "all studies",  target = "delta", data.type = "real", model = "oob", mean.auc = 0, pval = NA, padj = NA, sig = NA), 
                            data.frame(group = "all studies",  target = "delta", data.type = "null", model = "oob", mean.auc = 0, pval = NA, padj = NA, sig = NA)) %>% data.frame()

# overall lm
ml_loop_output_mean_sig = data.frame(lmerTest::lmer(auc ~ data.type + (1|group), subset(ml_loop_output, target == "delta")) %>% summary() %>% coef())[2,]
ml_loop_output_mean_sig$target = c("delta")
ml_loop_output_mean_sig$group = "all studies"
ml_loop_output_mean_sig$sig = ifelse(ml_loop_output_mean_sig$Estimate > 0 & ml_loop_output_mean_sig$`Pr...t..` < 0.05, "*", "")

# plot
ml_loop_output_rs = subset(ml_loop_output, target == "delta")
ml_loop_output_mean_rs = subset(ml_loop_output_mean, target == "delta")
# reorder studies based on mean AUC
ml_loop_output_rs = ml_loop_output_rs %>% 
  mutate(group = factor(group, levels = c("all studies", arrange(subset(ml_loop_output_mean, target == "delta" & data.type == "real" & group != "all studies"), mean.auc)$group)))
ml_loop_output_mean_rs = ml_loop_output_mean_rs %>% 
  mutate(group = factor(group, levels = c("all studies", arrange(subset(ml_loop_output_mean, target == "delta" & data.type == "real" & group != "all studies"), mean.auc)$group)))

ml_loop_output_rs.custom = ml_loop_output_rs%>%
  merge(rbind(loco.names.simple, data.frame(group = "all studies", clean = "All studies"))) %>%
  mutate(study = gsub(" .*", "", clean))
ml_loop_output_mean_rs.custom = ml_loop_output_mean_rs %>%
  merge(rbind(loco.names.simple, data.frame(group = "all studies", clean = "All studies"))) %>%
  mutate(study = gsub(" .*", "", clean))
ml_loop_output_rs.custom.forest = ml_loop_output_forest %>%
  merge(rbind(loco.names.simple, data.frame(group = "all studies", clean = "All studies")))%>%
  mutate(study = gsub(" .*", "", clean))
ml_loop_output_mean_sig.custom = ml_loop_output_mean_sig %>% 
  mutate(study = group) %>%
  merge(rbind(loco.names.simple, data.frame(group = "all studies", clean = "All studies")))%>%
  mutate(study = gsub(" .*", "", clean))


# :: :: PLOT: LOCO AUC ----------------------------------------------------

ml_loop_output_plot_rs =  ggplot(ml_loop_output_rs.custom)+
  #geom_hline(yintercept=0.5, linetype=2, color="grey")+
  geom_violin(data=subset(ml_loop_output_rs.custom, data.type == "null"), 
              aes(x=clean, y=auc), fill="grey", color=NA, alpha=0.4)+
  geom_point(data=subset(ml_loop_output_mean_rs.custom, data.type == "null"& group != "all studies"),
             aes(x=clean, y=mean.auc), shape=21, color="grey", fill="grey", size=2.5)+
  geom_violin(data=subset(ml_loop_output_rs.custom, data.type == "real"), 
              aes(x=clean, y=auc, fill=study), width=0.4, alpha=0.4)+
  geom_point(data=subset(ml_loop_output_mean_rs.custom, data.type == "real" & group != "all studies"),
             aes(x=clean, y=mean.auc, fill=study), shape=21, size=2.5)+
  geom_text(data=subset(ml_loop_output_mean_rs.custom, data.type == "real"),
            aes(x=clean, y=1.04, label=sig), vjust=0.75, size=6)+
  # overall mean
  geom_point(data=subset(ml_loop_output_rs.custom.forest, model == "oob" & data.type == "null"),
             aes(x=clean, y=mean_auc), shape=23, size=3, fill="grey")+
  geom_point(data=subset(ml_loop_output_rs.custom.forest,model == "oob" &  data.type == "real"),
             aes(x=clean, y=mean_auc), fill="black", shape=23, size=3)+
  geom_text(data=subset(ml_loop_output_mean_sig.custom, target == "butyrogens"),
            aes(x=clean, y=1.04, label=sig), vjust=0.75, size=6)+
  scale_fill_manual(values=c(RColorBrewer::brewer.pal(name="Set3",
                                                      n=length(unique(rs_target_rs_notscaled$study)))))+
  #facet_grid(group~target, scales="free_y")+
  scale_y_continuous(
    limits=c(0,1.05),
    breaks=seq(0, 1, by=0.1))+
  scale_x_discrete(limits = c("All studies", arrange(subset(ml_loop_output_mean_rs.custom, data.type == "real" & clean != "All studies"), mean.auc)$clean))+
  facet_wrap(~"Leave-one-cohort-out CV")+
  theme_classic()+theme(#axis.text.x=element_text(angle=45, hjust=1, vjust=1),
    #strip.background = element_blank(),
    strip.text = element_text(size=10),
    #panel.grid.major.x = element_line(color="lightgrey", linewidth=0.2),
    strip.text.y = element_blank(),
    panel.grid.major.x = element_line(linewidth=0.2, linetype=2, color="grey"),
    legend.position="none")+
  labs(x=NULL, y="AUC")+
  coord_flip()
ml_loop_output_plot_rs


# :: :: LOCO feature importance ----------------------------------------------

t0 <- Sys.time()
ml_loop_output_importances_rs = 
  do.call(rbind, lapply(iters, function(iter){
    do.call(rbind, lapply(data.types, function(data.type){
      
      # conditionally scramble target values, per study
      if(data.type == "null"){
        set.seed(iter)
        rs_ml_input = rs_target_rs_notscaled %>% 
          group_by(study) %>%
          mutate(delta = sample(delta)) %>%
          # redefine response
          mutate(delta.resp = ifelse(delta > median(delta), "responder", "nonresponder"))
      } else {
        rs_ml_input = rs_target_rs_notscaled
      }
      
      
      # otherwise, proceed
      do.call(rbind, lapply("delta", function(target){
        do.call(rbind, lapply(studies, function(study){
          # print status outside of mclapply applied on models
          t1 <- Sys.time()
          print(paste0(iter, " ", data.type, " ", target, " ", study, " ", round(t1-t0, digits=1), sep=""))
          do.call(rbind, parallel::mclapply(algos, function(algo){
            # apply function
            ml_loco_output = ml_loco(otu_table = rs_ml_data,
                                     metadata = rs_ml_input,
                                     heldout = study,
                                     target = target,
                                     output = "importance",
                                     algo = algo,
                                     auc = TRUE, # only calculate AUC, otherwise return predictions
                                     scale = T,
                                     iter = iter)
            ml_loco_output$data.type = data.type
            ml_loco_output
          }))}))}))}))}))
t2 = Sys.time() # 5.5 hours; 2 min for ranger oob
ml_loop_output_importances_rs

t2 - t0

# fix names
ml_loop_output_importances_rs = ml_loop_output_importances_rs %>%
  mutate(feature = gsub("g__CAG.", "g__CAG-",
                        gsub("NA.", "NA", feature)))

# calculate mean
ml_loop_output_importances_rs_mean = ml_loop_output_importances_rs %>%
  # subset to sig studies
  subset(group %in% c(subset(ml_loop_output_mean_rs.custom, sig == "*")$group))%>%
  group_by(data.type, model, target, feature, iter) %>%
  mutate(mean.imp.sig = mean(importance)) %>%
  group_by(data.type, model, target, feature) %>%
  mutate(mean.imp = mean(importance)) %>%
  mutate(imp.low = mean(importance) - (sd(importance)/sqrt(n()) * 1.96),
         imp.high = mean(importance) + (sd(importance)/sqrt(n()) * 1.96))%>%
  dplyr::select(data.type, model, feature, target, iter, mean.imp, mean.imp.sig, imp.low, imp.high) %>% distinct() %>% data.frame()

# run stats
ml_loop_output_importances_rs_pval = 
  do.call(rbind, lapply(unique(ml_loop_output_importances_rs_mean$feature), function(z){
    do.call(rbind, lapply("delta", function(y){
      data.subset = subset(ml_loop_output_importances_rs_mean, target == y & feature == z)
      ttest.results = wilcox.test(subset(data.subset, data.type=="real")$mean.imp.sig,
                                  subset(data.subset, data.type=="null")$mean.imp.sig)$p.value
      ttest.results = data.frame(target = y,
                                 pval = ttest.results,
                                 feature = z)
      ttest.results
    }))}))

ml_loop_output_importances_rs_pval
ml_loop_output_importances_rs_pval$padj = p.adjust(ml_loop_output_importances_rs_pval$pval, method="BH")
ml_loop_output_importances_rs_pval$sig = ifelse(ml_loop_output_importances_rs_pval$padj < 0.05, "*", "")

# merge
ml_loop_output_importances_rs_mean = merge(ml_loop_output_importances_rs_mean,
                                           ml_loop_output_importances_rs_pval, by=c("target", "feature"))
# delete if mean real < mean null
ml_loop_output_importances_rs_mean = ml_loop_output_importances_rs_mean %>%
  group_by(target, feature) %>%
  mutate(sig = ifelse(mean.imp > mean.imp[data.type == "null"], sig, ""))

# subset to sig 
ml_loop_output_importances_rs_sig = subset(ml_loop_output_importances_rs, feature %in% subset(ml_loop_output_importances_rs_mean, sig == "*")$feature)
ml_loop_output_importances_rs_mean_sig = subset(ml_loop_output_importances_rs_mean, feature %in% subset(ml_loop_output_importances_rs_mean, sig == "*")$feature)
# subset to overall mean
ml_loop_output_importances_mean_rs_sig = ml_loop_output_importances_rs_mean_sig[,c("target", "feature", "data.type", "model", "mean.imp", "sig", "imp.low", "imp.high")] %>% distinct()

# save for ex vivo
ml_rs_features = subset(ml_loop_output_importances_rs_mean_sig, target == "delta" & data.type == "real" & sig == "*")$feature %>% unique()
length(ml_rs_features)

saveRDS(ml_rs_features, "../2025_11_07_meta_rs_features.Rds")

# order features
ml_rs_features = factor(ml_rs_features,
                        levels = arrange(subset(ml_loop_output_importances_mean_rs_sig, data.type == "real"), -mean.imp)$feature)
ml_rs_features = as.character(levels(ml_rs_features))
ml_loop_output_importances_rs_sig$feature = (factor(ml_loop_output_importances_rs_sig$feature, levels=ml_rs_features))
ml_loop_output_importances_rs_mean_sig$feature = (factor(ml_loop_output_importances_rs_mean_sig$feature, levels=ml_rs_features))

ml_loop_output_importances_rs_sig$feature.order = as.numeric(factor(ml_loop_output_importances_rs_sig$feature, levels=ml_rs_features))
ml_loop_output_importances_rs_mean_sig$feature.order = as.numeric(factor(ml_loop_output_importances_rs_mean_sig$feature, levels=ml_rs_features))

n.features = 30

ml_rs_features_top = ml_rs_features[1:n.features]
lasso_coefs$feature
# indicate (in bold) whether taxa are detected in our data (and are in top 15)
ml_loop_output_importances_mean_rs_sig$detected = ifelse(ml_loop_output_importances_mean_rs_sig$feature %in%
                                                           colnames(ml.input.data[,colnames(ml.input.data) %in% ml_rs_features_top]), "bold", "plain")
ml_loop_output_importances_rs_mean_sig$top15 = ifelse(ml_loop_output_importances_rs_mean_sig$feature %in% colnames(ml.input.data[,colnames(ml.input.data) %in% ml_rs_features_top]), "Yes", "No")
ml_loop_output_importances_mean_rs_sig$top15 = ifelse(ml_loop_output_importances_mean_rs_sig$feature %in% colnames(ml.input.data[,colnames(ml.input.data) %in% ml_rs_features_top]), "Yes", "No")

# :: :: PLOT: feature importance ------------------------------------------


ml_loop_output_rs_importance_plot = ggplot(subset(ml_loop_output_importances_rs_sig, target=="delta"))+
  geom_hline(yintercept=0, linetype=2, linewidth=0.2, color="black")+
  #geom_violin(data=subset(subset(ml_loop_output_importances_sig, target=="butyrogens" & feature %in% ml_loop_butyrogens), data.type == "null"), 
  #            aes(x=feature, y=importance), fill="grey", color=NA, alpha=0.4)+
  #geom_violin(data=subset(subset(ml_loop_output_importances_sig, target=="butyrogens" & feature %in% ml_loop_butyrogens), data.type == "real"), 
  #            aes(x=feature, y=importance), width=0.4, fill="red", alpha=0.4)+
  geom_segment(data=subset(subset(ml_loop_output_importances_rs_mean_sig, target=="delta"), data.type == "null"),
               aes(x=reorder(feature,feature.order), xend=reorder(feature,feature.order), y=imp.low, yend=imp.high), color="grey")+
  geom_segment(data=subset(subset(ml_loop_output_importances_rs_mean_sig, target=="delta"), data.type == "real"),
               aes(x=reorder(feature,feature.order), xend=reorder(feature,feature.order), y=imp.low, yend=imp.high,
                   alpha=top15), color="red")+
  geom_point(data=subset(subset(ml_loop_output_importances_mean_rs_sig, target=="delta"), data.type == "null"),
             aes(x=reorder(feature,mean.imp), y=mean.imp), color="grey")+
  geom_point(data=subset(subset(ml_loop_output_importances_mean_rs_sig, target=="delta"), data.type == "real"),
             aes(x=reorder(feature,mean.imp), y=mean.imp), shape=21, fill="white")+
  geom_point(data=subset(subset(ml_loop_output_importances_mean_rs_sig, target=="delta"), data.type == "real"),
             aes(x=reorder(feature,mean.imp), y=mean.imp, alpha=top15), shape=21, fill="red")+
  # indicate whether detected in our data
  #geom_text(data = subset(ml_loop_output_importances_mean_rs_sig, data.type == "real" & mean.imp > 0.04),
  #          aes(y=mean.imp-0.005, x=reorder(feature,mean.imp), label=feature, fontface=detected), hjust=1, size=3)+
  #geom_text(data = subset(ml_loop_output_importances_mean_rs_sig, data.type == "real" & mean.imp < 0.04),
  #          aes(y=mean.imp+0.005, x=reorder(feature,mean.imp), label=feature, fontface=detected), hjust=0, size=3)+
  theme_classic()+theme(#axis.text.x=element_text(angle=45, hjust=1, vjust=1),
    strip.text.y = element_blank(),
    #axis.text.y = element_blank(),
    #axis.ticks.y = element_blank(),
    strip.text=element_text(size=10),
    legend.position = c(0.85, 0.2))+
  scale_x_discrete(limits=rev(unique(arrange(subset(ml_loop_output_importances_rs_sig, target=="delta" & feature.order <= n.features),feature.order)$feature)))+
  guides(alpha = guide_legend(reverse = TRUE))+
  labs(x="", y="Mean Decrease in Accuracy",
       alpha="Detected")+
  facet_wrap(~"Feature Importances")+
  coord_flip()
ml_loop_output_rs_importance_plot


# 5. Primary Degraders ----------------------------------------------------

# visualize whether Primary degraders are differentially abundant
# in fermenters vs non-fermenters

# loop through and perform wilcox on strong vs weak

ml.input.data.for.wilcox = ml.input.data %>% as.data.frame() %>%
  mutate(HM = rownames(.)) %>%
  merge(subset(ml.target.data, variable=="delta_pH"), by="HM") %>%
  mutate(response = ifelse(value < -1.27, "strong", "weak")) 

ml.input.data.for.wilcox.output = do.call(rbind, lapply(colnames(ml.input.data), function(taxa){
  do.call(rbind, lapply(rs.names, function(rs){
    data.subset = ml.input.data.for.wilcox[,colnames(ml.input.data.for.wilcox) %in%
                                             c(taxa, "RS_Name", "value", "response")]
    colnames(data.subset)[1] = "taxa"
    data.subset = subset(data.subset, RS_Name == rs)
    test.output = wilcox.test(subset(data.subset, response == "strong")$taxa,
                              subset(data.subset, response == "weak")$taxa)
    data.frame(RS_Name = rs,
               taxa = taxa,
               pval = test.output$p.value)
  }))}))
ml.input.data.for.wilcox.output$padj = p.adjust(ml.input.data.for.wilcox.output$pval, method="BH")
ml.input.data.for.wilcox.output$sig = ifelse(ml.input.data.for.wilcox.output$padj < 0.25, "*", "")

subset(ml.input.data.for.wilcox.output, sig == "*")


# :: :: PLOT: Primary Degraders -------------------------------------------

subset(ml.input.data.for.wilcox.output, sig == "*")
r.callidus.p = subset(ml.input.data.for.wilcox.output, taxa == "g__Ruminococcus_C_58660_s__callidus") %>%
  mutate(RS_Name = factor(RS_Name, levels=rs.names)) 


r.callidus.plot = ml.input.data %>% as.data.frame()%>%
  dplyr::select(g__Bifidobacterium_388775_s__faecale, 
                g__Ruminococcus_C_58660_s__callidus,
                g__Ruminococcus_E_s__bromii_B,
                g__Ruminococcus_E) %>%
  mutate(HM = rownames(.)) %>%
  merge(subset(ml.target.data, variable=="delta_pH"), by="HM") %>%
  mutate(response = ifelse(value < -1.27, "strong", "weak")) %>%
  mutate(RS_Name = factor(RS_Name, levels=rs.names)) %>%
  mutate(response = ifelse(response=="strong", "Strong", "Weak")) %>%
  ggplot(aes(x=response, y=(g__Ruminococcus_E_s__bromii_B+0.0005)))+
  geom_text(data=r.callidus.p,
            aes(label=paste0("FDR = ", round(padj, digits=3)),
                fontface = ifelse(padj < 0.20, "bold", "plain")), x=1.5, y=log10(6.8), size=3,
            color="black")+
  scale_y_log10(limits =c(0.0005, 7))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(shape=21, aes(fill=response), width=0.2)+
  theme_classic()+theme(legend.position="none")+
  labs(x=NULL, y="R. callidus %")+
  facet_wrap(~RS_Name, nrow=1)
r.callidus.plot

r.bromii.p = subset(ml.input.data.for.wilcox.output, taxa == "g__Ruminococcus_E_s__bromii_B") %>%
  mutate(RS_Name = factor(RS_Name, levels=rs.names)) 

r.bromii.plot = ml.input.data %>% as.data.frame()%>%
  dplyr::select(g__Ruminococcus_E_s__bromii_B) %>%
  mutate(HM = rownames(.)) %>%
  merge(subset(ml.target.data, variable=="delta_pH"), by="HM") %>%
  mutate(response = ifelse(value < -1.27, "strong", "weak")) %>%
  mutate(RS_Name = factor(RS_Name, levels=rs.names)) %>%
  mutate(response = ifelse(response=="strong", "Strong", "Weak")) %>%
  ggplot(aes(x=response, y=(g__Ruminococcus_E_s__bromii_B+0.0005)))+
  geom_text(data=r.bromii.p,
            aes(label=paste0("FDR = ", round(padj, digits=3)),
                fontface = ifelse(padj < 0.20, "bold", "plain")), x=1.5, y=log10(20), size=3,
            color="black")+
  scale_y_log10(limits =c(0.0005, 20))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(shape=21, aes(fill=response), width=0.2)+
  theme_classic()+theme(legend.position="none")+
  labs(x=NULL, y="R. bromii %")+
  facet_wrap(~RS_Name, nrow=1)
r.bromii.plot

b.adol.p = subset(ml.input.data.for.wilcox.output, taxa == "g__Bifidobacterium_388775_s__faecale") %>%
  mutate(RS_Name = factor(RS_Name, levels=rs.names)) 
b.adol.plot = ml.input.data %>% as.data.frame()%>%
  dplyr::select(g__Bifidobacterium_388775_s__faecale) %>%
  mutate(HM = rownames(.)) %>%
  merge(subset(ml.target.data, variable=="delta_pH"), by="HM") %>%
  mutate(response = ifelse(value < -1.27, "strong", "weak")) %>%
  mutate(RS_Name = factor(RS_Name, levels=rs.names)) %>%
  mutate(response = ifelse(response=="strong", "Strong", "Weak")) %>%
  ggplot(aes(x=response, y=(g__Bifidobacterium_388775_s__faecale+0.0005)))+
  geom_text(data=b.adol.p,
            aes(label=paste0("FDR = ", round(padj, digits=3)),
                fontface = ifelse(padj < 0.20, "bold", "plain")), x=1.5, y=log10(1), size=3,
            color="black")+
  scale_y_log10(limits =c(0.0005, 1))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(shape=21, aes(fill=response), width=0.2)+
  theme_classic()+theme(legend.position="none")+
  labs(x=NULL, y="B. faecale %")+
  facet_wrap(~RS_Name, nrow=1)
b.adol.plot

r.callidus.plot/
  r.bromii.plot/
  b.adol.plot


# // old FIGURES --------------------------------------------------------------


## TOP
(ml.rapidaim.ph.plot+
   patchwork::free(dbrda_model_ph_plot, type="label")+
   patchwork::free(ml.rf.regression.oob.oars.ph.pbs.plot.neat, type="label")) %>%
  ggsave(
    #filename="../../ml_plots/oars_ml_plot_top.pdf",
    width=14,
    height=5,
    units="in",
    device=cairo_pdf)

reshape2::acast(subset(ml.rf.regression.oob.oars.ph.pbs.importances.mean.sig, type == "real"),
                RS_Name ~ feature, value.var="adj.imp") %>%
  t() %>%
  pheatmap::pheatmap(color = colorRampPalette(c("red", "white", "blue"))(100),
                     annotation_col = data.frame(RS = rs.names.pbs)%>% `rownames<-`(rs.names.pbs),
                     annotation_colors = list(RS = c(labelcolors.rs, "PBS" = "black")),
                     annotation_legend = F,
                     # display_numbers = ml.rapidaim.data.ph.lfc.sig.mat.stars[,rs.names],
                     fontsize_number = 10,number_color="white",
                     breaks=c(seq(min(na.omit(ml.rf.regression.oob.oars.ph.pbs.importances.mean.sig$adj.imp)), 0, length.out=ceiling(100/2) + 1), 
                              seq(max(na.omit(ml.rf.regression.oob.oars.ph.pbs.importances.mean.sig$adj.imp))/100, max(na.omit(ml.rf.regression.oob.oars.ph.pbs.importances.mean.sig$adj.imp)), length.out=floor(100/2))),
                     border_color = "white") %>%
  ggsave(
    #filename="../../ml_plots/oars_ml_plot_heatmap.pdf",
    width=5,
    height=5,
    units="in",
    device=cairo_pdf)

bottom.1 = (rs_maaslin_all_plot/patchwork::free(rs_target_rs_notscaled_plot, type="label"))+
  patchwork::plot_layout(heights=c(1.5,1))

## BOTTOM rs_target_rs_notscaled_plot
bottom.2a = (patchwork::free(ml_loop_output_plot_rs, type="label"))
bottom.2b = (ml_loop_output_rs_importance_plot)

bottom.3 = 
  patchwork::free(ml.rf.regression.oob.oars.ph.pbs.fs.cor.plot, type="label")+
  patchwork::free(ml.fs.loocv.regression.results.ph.pbs.auc.mean.plot, type="label")+
  patchwork::free(ml.rf.regression.oob.oars.ph.pbs.fs.auc.plot, type="label")+
  patchwork::free(top.model.preds.plot, type="label")+
  patchwork::free(ml.rf.oob.auc.ph.pbs.fs.roc.plot, type="label")+
  patchwork::free(top.model.preds.cor.plot, type="label")+
  patchwork::plot_layout(nrow=3)
bottom.3



bottom.4 = patchwork::free(ml.sample.size.curve.linear+facet_wrap(~RS_Name, nrow=3), type="label")


cowplot::plot_grid(bottom.1, 
                   cowplot::plot_grid(bottom.2a, bottom.2b, nrow=2), 
                   bottom.3, 
                   bottom.4,
                   nrow=1, rel_widths=c(1.5,1,1.5,1.5)) %>%
  ggsave(
    #filename="../../ml_plots/oars_ml_plot_bottom.pdf",
    width=30,
    height=9,
    units="in",
    device=cairo_pdf)

r.callidus.plot %>%
  ggsave(
    #filename="../../ml_plots/oars_ml_plot_taxa.pdf",
    width=18,
    height=3,
    units="in",
    device=cairo_pdf)


(ml.sample.size.curve.linear/
    ml.sample.size.curve.nonlinear)%>%
  ggsave(
    #filename="../../ml_plots/oars_ml_plot_sample_size.pdf",
    width=18,
    height=5,
    units="in",
    device=cairo_pdf)

