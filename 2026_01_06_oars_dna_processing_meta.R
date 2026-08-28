### 2025_05_24  OARS DNA Processing (including metadata)

# :: load packages --------------------------------------------------------
library("ggplot2"); library("dplyr"); library("tidyverse"); library("patchwork")

rs.names <- c("Authentic", "BobsRedMill", "MSPrebiotic", "LetsDoOrganic", "HiMaize260", "Novelose330", "ActistarRT", "FibersymRW", "Versafibe1490")

rs.names.pbs <- c("PBS", "Authentic", "BobsRedMill", "MSPrebiotic", "LetsDoOrganic", "HiMaize260", "Novelose330", "ActistarRT", "FibersymRW", "Versafibe1490")



# :: necessary files ------------------------------------------------------

# OARS patient list
"2023_12_04_oars_patient_list.csv"

# Master mapping file (for MBX)
"2024_11_08_metabolomics_samples_prioritized_PD.csv"

# James' ASV table
"2025_09_08_full_mapping_rapidaim.jb_gg2species_physeq_pooled_wTree_250909.Rds"

# James' MGX manifest
"humann3_outputs_summary_250502.csv"

# James' MGX files (in several directories)
"metaphlan_bugs_list"

# RapidAIM data
"2025_06_09_oars_scores.Rds"




# >> 16S ------------------------------------------------------------------


# :: process ALL 16S data --------------------------------------------------------
process.16s.data = FALSE
if(process.16s.data == TRUE){
  # load updated data
  load("./all_trials_16S_min50k_mappingFile_250530_gg2species_allData_250910.Rdata")
  amplicon.data.gg <- physeq.pooled.wTree
  # 
  amplicon.data.gg.tax = phyloseq::tax_table(amplicon.data.gg)
  # sum(row.names(data.frame(phyloseq::refseq(amplicon.data.gg))) == colnames(data.frame(phyloseq::otu_table(amplicon.data.gg))))
  # add LCA
  amplicon.data.gg.tax <- data.frame(amplicon.data.gg.tax) %>%
    mutate(LCA = 
             ifelse(!is.na(Species)&Species!="s__", paste(as.character(Genus), as.character(Species), sep="_"), 
                    ifelse(!is.na(Genus)&Genus!="g__", paste(Genus), 
                           ifelse(!is.na(Family)&Family!="f__", paste(Family),
                                  ifelse(!is.na(Order)&Order!="o__", paste(Order),
                                         ifelse(!is.na(Class)&Class!="c__", paste(Class),
                                                ifelse(!is.na(Phylum)&Phylum!="p__", paste(Phylum),
                                                       ifelse(is.na(Phylum), "undefined", paste(Phylum)))))))))
  # call this "Taxa", so that LCA remained the "glom'ed" name
  amplicon.data.gg.tax$Taxa <- make.unique(as.character(amplicon.data.gg.tax$LCA), sep = "_")
  # replace DNA sequence with ASV# to enable merging with phyloseq object
  rownames(amplicon.data.gg.tax) = row.names(data.frame(phyloseq::refseq(amplicon.data.gg)))
  # merge updated tax table
  amplicon.data.gg = phyloseq::merge_phyloseq(otu_table = phyloseq::otu_table(amplicon.data.gg),
                                              sample_data = phyloseq::sample_data(amplicon.data.gg),
                                              tax_table = phyloseq::tax_table(as.matrix(amplicon.data.gg.tax)),
                                              refseq = phyloseq::refseq(amplicon.data.gg))
  # identify ASVs annotated as mitochondria /chloroplast
  mitochondria.indices = grep("itochondria", data.frame(phyloseq::tax_table(amplicon.data.gg))$Family) # 12 ASVs = Mitochondria
  chloroplast.indices = grep("hloroplast", data.frame(phyloseq::tax_table(amplicon.data.gg))$Class) # 3 ASVs = Chloroplast
  # remove these ASVs
  asvs.to.keep = data.frame(phyloseq::tax_table(amplicon.data.gg)[-c(mitochondria.indices,chloroplast.indices)]) %>% row.names()
  amplicon.data.gg = phyloseq::prune_taxa(asvs.to.keep, amplicon.data.gg)
  # rarefy
  set.seed(25)
  seq.depth = 50000 # note: range is 47,931 to >4 million (but need to check which samples, may not be OARS)
  amplicon.data.gg.rare <- phyloseq::rarefy_even_depth(amplicon.data.gg, 
                                                       sample.size = seq.depth,
                                                       rngseed = 25, replace = FALSE)
  
  # redo tax assignment using gg13.8 (for Butyrogen calculation) (used for LSARP, too)
  amplicon.data.gg138.tax = dada2::assignTaxonomy((phyloseq::refseq(amplicon.data.gg.rare)), "~/Documents/PhD/16S_databases/gg_13_8_train_set_97.fa.gz")
  amplicon.data.gg138.tax.df = data.frame(amplicon.data.gg138.tax)
  amplicon.data.gg138.tax.df$ASV = names(phyloseq::refseq(amplicon.data.gg.rare)) 
  # and redo tax assignment using RDP (for 16S MLP) (used for LSARP, too)
  amplicon.data.rdp.tax = dada2::assignTaxonomy((phyloseq::refseq(amplicon.data.gg.rare)), "~/Documents/PhD/16S_databases/rdp_train_set_16.fa.gz")
  amplicon.data.rdp.tax.df = data.frame(amplicon.data.rdp.tax)
  amplicon.data.rdp.tax.df$ASV = names(phyloseq::refseq(amplicon.data.gg.rare)) 
  
  # fix standard.name column
  #phyloseq::sample_data(amplicon.data.gg.rare)$standard.name = phyloseq::sample_data(amplicon.data.gg.rare)$Standard_name
  #phyloseq::sample_data(amplicon.data.gg.rare)$standard.name = phyloseq::sample_data(amplicon.data.gg.rare)$Standard_name
  
  # save objects
  saveRDS(amplicon.data.gg.rare, "./2025_09_10_rs_trial_16s_data_rarefied_gg2.Rds")
  saveRDS(amplicon.data.gg, "./2025_09_10_rs_trial_16s_data_unrarefied_gg2.Rds") # note: this doesn't get used anymore
  saveRDS(amplicon.data.gg138.tax.df, "./2025_09_10_rs_trial_16s_data_tax_table_gg2.Rds") # note: misnamed "gg2", should be "gg138"
  saveRDS(amplicon.data.rdp.tax.df, "./2025_09_10_rs_trial_16s_data_tax_table_rdp.Rds")
  # amplicon.data.gg
}
# load

# :: prepare OARS mapping ----------------------------------------------------

# obtain list of OARS patients for analysis
oars.patient.list <- read.csv("./oars_git_mapping/2023_12_04_oars_patient_list.csv")
# create variable to link to REDCap-derived data (study_id)
oars.patient.list$study_id = paste(oars.patient.list$HM, ".00", sep="")
# check that stools are OARS and not other study (this information is in MBX mapping file)
metadata.oars.stool <- read.csv("./2024_11_08_metabolomics_samples_prioritized_PD.csv")
metadata.oars.stool = subset(metadata.oars.stool, grepl("OARS", stl_study))
metadata.oars.stool[,c("standard.name")] %>% unique() %>% length() # 71
# subset to those intended for analysis (i.e. not all OARS participants are to be analyzed)
metadata.oars.stool = subset(metadata.oars.stool, HM %in% oars.patient.list$HM)
metadata.oars.stool$study_id = paste(metadata.oars.stool$HM, ".00", sep="")
nrow(metadata.oars.stool) # 160 samples to analyze
# keep clean
metadata.oars.stool = metadata.oars.stool[,c("HM", "standard.name",
                                             "stool_date_rec_v2",
                                             "oars_start_date", "oars_product_end_date",
                                             "oars.timing","oars.days","oars.off",
                                             "oars.rs.1", "oars.rs.2")]

# add "phase"; phase 1 = baseline sample + on RS; phase 2 = off RS (washout)
# note: for statistics, baseline sample is added to phase 2 as comparison group
metadata.oars.stool$phase = ifelse(grepl(paste(c("baseline", "3", "6"), collapse="|"),metadata.oars.stool$oars.timing), "treatment", "washout")
# identify baseline sample
metadata.oars.stool$baseline  =ifelse(grepl("baseline",metadata.oars.stool$oars.timing), "baseline", "not_baseline")

# good
metadata.oars.stool = metadata.oars.stool %>% distinct() %>% data.frame()
dim(metadata.oars.stool)
# 68 unique samples

## add other notes
metadata.oars.stool = merge(metadata.oars.stool,
                            oars.patient.list[,c("HM", "diagnosis", "sex",
                                                 "dose_altered","doses_patient",
                                                 "missed_patient", "rs_notes",
                                                 "rs_1_compliance", "rs_2_compliance",
                                                 "clinic_notes","medication_notes")], by="HM")
# clean timing

# add simpler timing variable
metadata.oars.stool$timing = 
  ifelse(grepl("baseline", metadata.oars.stool$oars.timing), "0M",
         ifelse(grepl("3", metadata.oars.stool$oars.timing), "3M",
                ifelse(grepl("6", metadata.oars.stool$oars.timing), "6M",
                       ifelse(grepl("9", metadata.oars.stool$oars.timing), "9M",
                              ifelse(grepl("12", metadata.oars.stool$oars.timing), "12M", NA)))))

# add RS_Name for RS at that time point (on RS; for coloring plots)
metadata.oars.stool$RS_Name = ifelse(metadata.oars.stool$timing == "3M", metadata.oars.stool$oars.rs.1,
                                     ifelse(metadata.oars.stool$timing == "6M", metadata.oars.stool$oars.rs.2, NA))
# add "oars.on.rs"
metadata.oars.stool$oars.on.rs = ifelse(metadata.oars.stool$timing == "0M", "preRS",
                                        ifelse(metadata.oars.stool$timing %in% c("3M", "6M"), "onRS", "postRS"))
# fix outlier
metadata.oars.stool$oars.on.rs = ifelse(metadata.oars.stool$standard.name == "HM0902-STL-06", "onRS",
                                       metadata.oars.stool$oars.on.rs)

metadata.oars.stool$oars.on.rs = factor(metadata.oars.stool$oars.on.rs, levels=c("preRS", "onRS", "postRS"))


# >> Load REDCAP ------------------------------------------------------

# 2025_07_08  Goal: add supplemented RS to RS intake

# load REDCap for RS dose
redcap.data.loaded <- read.csv("~/Downloads/redCap_sampleData_uniques_240724.tsv", sep="\t")

# subset to OARS samples
redcap.data.loaded = subset(redcap.data.loaded, study_id %in% paste(metadata.oars.stool$HM, ".00", sep=""))
# generate clinical table


# 1. take diagnosis values
subset(redcap.data.loaded, redcap.event.name.simple == "diagnosis")

redcap.data.loaded[,c("study_id")]
redcap.data.loaded$age_calculated # age at diagnosis
redcap.data.loaded$paris.current.cd.calc # location & behaviour

redcap.data.loaded$wpcdai_category # severity
redcap.data.loaded$pucai_sev # severity

# time from diagnosis to study start
redcap.data.loaded$oars_start_date

# need to do current drug manually
run.this = F
if(run.this == T){
metadata.oars.stool.clinical.data = do.call(rbind, lapply(unique(paste(metadata.oars.stool$HM, ".00", sep="")), function(x){
  print(x)
  # first, subset to HM
  redcap.data.hm = subset(redcap.data.loaded, study_id == x)
  # second, separate "diagnosis" from "enrollment in oars"
  redcap.data.hm.diag = subset(redcap.data.hm, redcap.event.name.simple == "diagnosis")
  redcap.data.hm.oars = subset(redcap.data.hm, stl_oars_time %in% c("OARS baseline") | mend_stl_time == "12 +/- 2 months end of trial")
  # > sex
  hm.sex = unique(redcap.data.hm$gender)
  # > age at entry
  hm.age.at.entry = unique(redcap.data.hm.oars$current.age)
  # > age at diagnosis
  hm.age.at.diagnosis = unique(redcap.data.hm.diag$current.age)
  # > clinical severity at entry
  hm.clinical.severity.at.entry = unique(na.omit(redcap.data.hm.oars$wpcdai_category,
                                                 redcap.data.hm.oars$pucai_sev))
  # > clinical severity at diagnosis
  hm.clinical.severity.at.diagnosis = unique(na.omit(redcap.data.hm.diag$wpcdai_category,
                                                     redcap.data.hm.diag$pucai_sev))
  # > site of disease at diagnosis
  hm.clinical.site.at.diagnosis = unique(redcap.data.hm.diag$paris.current.cd.calc)
  # > time from diagnosis to study start
  hm.time.to.start =  as.numeric(as.Date(min(na.omit(redcap.data.hm.oars$sample.collection.date)), "%Y-%m-%d") - as.Date(min(na.omit(redcap.data.hm.diag$sample.collection.date)), "%Y-%m-%d")) / 30.5
  
  # save as dataframe
  data.frame(HM = x,
             sex = hm.sex,
             age_at_entry = hm.age.at.entry,
             age_at_diagnosis = hm.age.at.diagnosis,
             clinical_severity_at_entry = hm.clinical.severity.at.entry,
             clinical_severity_at_diagnosis = hm.clinical.severity.at.diagnosis,
             clinical_site_at_diagnosis = hm.clinical.site.at.diagnosis,
             months_to_start = hm.time.to.start)
}))
  
}
# >> FFQ ------------------------------------------------------------------


# subset to OARS samples
redcap.data.loaded = subset(redcap.data.loaded, standard.name %in% metadata.oars.stool$standard.name)
nrow(redcap.data.loaded) # 68 redcap FFQ entries

# note: we can use FFQ as predictors for ML models,
# so let's prepare that data here:
redcap.data.ffq = redcap.data.loaded %>% select(apple:other_take_away_meals)
# change to ordinal (note: I end up subtracting one later for filtering purposes)
redcap.data.ffq.map = data.frame(entry = unique(unlist(redcap.data.ffq)),
           value = rank(c(2/7, 5/7, 1/7, 1/1, 1/30, 2/30, 0, 2/1, 4/1, 6/1))-1) %>%
  arrange(value)
# apply replacement
redcap.data.ffq = do.call(cbind, lapply(1:ncol(redcap.data.ffq), function(x){
  data.subset = data.frame(feature = redcap.data.ffq[,x])
  # apply replacement
  data.subset$feature = redcap.data.ffq.map$value[match(data.subset$feature,redcap.data.ffq.map$entry)]
  colnames(data.subset) = colnames(redcap.data.ffq)[x]
  data.subset
}))
# add standard.name
redcap.data.ffq$standard.name = redcap.data.loaded$standard.name
redcap.data.ffq
# save
rownames(redcap.data.ffq) = redcap.data.ffq$standard.name
redcap.data.ffq$standard.name = NULL
saveRDS(redcap.data.ffq, "./2025_08_16_oars_ffq_table.Rds")


# add nutrient data (hypothesis: Mg helps RS fermentation)
# so let's prepare that data here:
redcap.data.nutrients = redcap.data.loaded %>% select(energy:b_carotene_equivalents)
# change to ordinal (note: I end up subtracting one later for filtering purposes)

# add standard.name
redcap.data.nutrients$standard.name = redcap.data.loaded$standard.name
redcap.data.nutrients
# save
rownames(redcap.data.nutrients) = redcap.data.nutrients$standard.name
redcap.data.nutrients$standard.name = NULL
saveRDS(redcap.data.nutrients, "./2025_08_16_oars_nutrient_table.Rds")

# >> Calculate Fibre + RS ------------------------------------------------------

# clean to important columns
redcap.data.loaded = subset(redcap.data.loaded, !is.na(oars_dose))[,c("standard.name", "fecalcal_res",
                                                                      "energy", "resistant_starch","dietary_fibre", "oars_dose", "oars_dose_2", "oars_dose_3",
                                                                      "oars_dose_ro", "oars_dose_ro_2", "oars_dose_ro_3",
                                                                      "oars_percent_starch", "oars_percent_starch_2", "oars_percent_starch_3",
                                                                      "oars_percent_starch_ro", "oars_percent_starch_ro_2", "oars_percent_starch_ro_3",
                                                                      "oars_lot", "oars_lot_2", "oars_lot_3")]

## calculate adj.fibre and adj.rs
# adjust g for kcal energy intake
redcap.data.loaded$fibre_adj = redcap.data.loaded$dietary_fibre / redcap.data.loaded$energy * 4.184 * 1000
redcap.data.loaded$rs_adj = redcap.data.loaded$resistant_starch / redcap.data.loaded$energy * 4.184 * 1000
# merge with metadata
redcap.data.loaded = merge(redcap.data.loaded,
                           metadata.oars.stool[,c("timing", "standard.name")], by="standard.name")

# replace Inf with NA
redcap.data.loaded$oars_dose_ro_2[is.infinite(redcap.data.loaded$oars_dose_ro_2)] = NA
redcap.data.loaded$oars_dose_ro_3[is.infinite(redcap.data.loaded$oars_dose_ro_3)] = NA
redcap.data.loaded$oars_percent_starch_ro_2[(redcap.data.loaded$oars_percent_starch_ro_2)==0] = NA
redcap.data.loaded$oars_percent_starch_ro_3[(redcap.data.loaded$oars_percent_starch_ro_3)==0] = NA

# take average RS% of lots per RS period
redcap.data.loaded = redcap.data.loaded %>%
  group_by(standard.name) %>%
  mutate(rs_dose = ifelse(timing == "3M", mean(na.omit(c(oars_dose,oars_dose_2,oars_dose_3))),
                          ifelse(timing == "6M", mean(na.omit(c(oars_dose_ro,oars_dose_ro_2,oars_dose_ro_3))),
                                 0))) %>%
  mutate(rs_perc = ifelse(timing == "3M", mean(na.omit(c(oars_percent_starch,oars_percent_starch_2,oars_percent_starch_3))),
                          ifelse(timing == "6M", mean(na.omit(c(oars_percent_starch_ro,oars_percent_starch_ro_2,oars_percent_starch_ro_3))),
                                 0))) %>%
  mutate(rs_intake = rs_dose * rs_perc/100) %>% data.frame()

# add calories from starch (note: RS calories != starch calories)
redcap.data.loaded$kcal = redcap.data.loaded$energy / 4.184
redcap.data.loaded$rs.kcal = 
  # add calories from RS fraction
  (redcap.data.loaded$rs_dose * redcap.data.loaded$rs_perc/100 * 2.7) + # cite: Miketinas, 2020 https://pubmed.ncbi.nlm.nih.gov/32840627/
  # plus calories from non-RS fraction
  (redcap.data.loaded$rs_dose * (1-redcap.data.loaded$rs_perc/100) * 4) # 4 kcal/g carbohydrates

# add energy from diet + RS supplement
redcap.data.loaded$energy = redcap.data.loaded$kcal + redcap.data.loaded$rs.kcal
# add RS from diet + RS supplement
redcap.data.loaded$total_rs = redcap.data.loaded$resistant_starch + redcap.data.loaded$rs_intake
# re-normalize RS intake to energy
redcap.data.loaded$total_rs_adj = redcap.data.loaded$total_rs / redcap.data.loaded$energy * 1000

# clean and merge with other data
oars.nutrient.data = merge(redcap.data.loaded[,c("standard.name", "energy","dietary_fibre", "resistant_starch", "total_rs","fibre_adj","rs_adj", "total_rs_adj")],
                           metadata.oars.stool, by="standard.name", all.y=T) # keep all entries

# :: ``` longitudinal fiber plot -----------------------------------------------------------

metadata.oars.fiber.plot = ggplot(oars.nutrient.data,
                                  aes(x=oars.days, 
                                      y=fibre_adj))+
  annotate("rect", xmin=0, xmax=max(subset(oars.nutrient.data, oars.on.rs=="onRS")$oars.days),
           ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(color="black",  se=T, linewidth=0.5)+
  geom_point(aes(fill=oars.on.rs, shape=ifelse(oars.on.rs=="preRS", "1", "2")), size=2.5)+
  scale_fill_manual(values=c("grey",2,"grey"))+
  scale_alpha_manual(values=c(1,0.2))+
  scale_shape_manual(values=c(23,21,21))+
  theme_classic()+theme(legend.position="none")+
  labs(x="Days since starting RS", y="Fiber Intake (g per day / 1000 kcal)")
metadata.oars.fiber.plot

# :: ``` average fiber plot -----------------------------------------------------------

oars.nutrient.data.average.plot = ggplot()+
  geom_segment(data=oars.nutrient.data %>%
                 group_by(HM) %>%
                 mutate(min.fiber = min(na.omit(fibre_adj)),
                        max.fiber = max(na.omit(fibre_adj)),
                        mean.fiber = mean(na.omit(fibre_adj))), 
               aes(x=min.fiber, y=reorder(HM, mean.fiber),
                   xend=max.fiber, yend=reorder(HM, mean.fiber)),
               linewidth=0.3)+
  geom_point(data=oars.nutrient.data %>% group_by(HM) %>% 
               mutate(mean.fib.adj = mean(na.omit(fibre_adj))) %>%
               dplyr::select(HM, mean.fib.adj) %>% distinct(),
             aes(x=mean.fib.adj,
                 y=reorder(HM, mean.fib.adj), 
                 fill=scale(mean.fib.adj)), shape=21, size=3)+
  coord_flip()+
  scale_fill_gradient2(low="blue", high="red")+
  theme_classic()+theme(legend.position="none",
                        axis.text.x=element_blank())+
  labs(x="Average Fiber Intake (g per day / 1000 kcal)", y="Participant", fill = "")
oars.nutrient.data.average.plot


# >> Fiber vs RS ----------------------------------------------------------

oars.nutrient.data.fiber.rs.plot = ggplot(oars.nutrient.data %>%
                                            subset(!is.na(fibre_adj)),
                                          aes(x=fibre_adj, rs_adj))+
  geom_path(aes(group=HM), linetype=2, alpha=0.5)+
  geom_point(shape=21, color="white", fill="black")+
  geom_smooth(method="lm", color="black", fill="black", se=T)+
  ggpubr::stat_cor(method="spearman")+
  theme_classic()+
  labs(x="Fiber Intake g / 1000 kcal per day",
       y="RS Intake g / 1000 kcal per day")
oars.nutrient.data.fiber.rs.plot



# :: ``` longitudinal rs plot ------------------------------------------------------------

metadata.oars.rs.plot = ggplot(oars.nutrient.data,
                                  aes(x=oars.days, 
                                      y=resistant_starch))+
  annotate("rect", xmin=0, xmax=max(subset(oars.nutrient.data, oars.on.rs=="onRS")$oars.days),
           ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(color="white",  se=T)+
  geom_point(aes(fill=oars.on.rs, 
                 shape=ifelse(oars.on.rs=="preRS", "1", "2")), 
             color="white", size=3)+
  scale_fill_manual(values=c("grey",2,"grey"))+
  scale_alpha_manual(values=c(1,0.2))+
  scale_shape_manual(values=c(23,21,21))+
  theme_classic()+theme(legend.position="none")+
  labs(x="Days since starting RS", y="Dietary Resistant Starch Intake (g per day)")
metadata.oars.rs.plot


metadata.oars.rs.supp.plot = ggplot(oars.nutrient.data,
                                    aes(x=oars.days, 
                                        y=total_rs))+
  annotate("rect", xmin=0, xmax=max(subset(oars.nutrient.data, oars.on.rs=="onRS")$oars.days),
           ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(color="white",  se=T)+
  geom_point(aes(fill=oars.on.rs, shape=ifelse(oars.on.rs=="preRS", "1", "2")), 
             color="white", size=3)+
  scale_fill_manual(values=c("grey",2,"grey"))+
  scale_alpha_manual(values=c(1,0.2))+
  scale_shape_manual(values=c(23,21,21))+
  theme_classic()+theme(legend.position="none")+
  labs(x="Days since starting RS", y="Total Resistant Starch Intake (g per day)")
metadata.oars.rs.supp.plot

metadata.oars.rs.plot+metadata.oars.rs.supp.plot


# :: mean intake ----------------------------------------------------------

mean(na.omit(oars.nutrient.data$fibre_adj))
# 14.02605

sd(na.omit(oars.nutrient.data$fibre_adj))
# 3.090316

# :: impute average fiber -------------------------------------------------

# Since dietary fiber intake is so variable (within- and between- individuals)
# we'll include as a covariate in certain analyses
# However, we need to impute missing values (mostly the latter values)
# So we'll calculate the average for each patient and use that

# use average fiber intake to impute NA values
oars.nutrient.data = oars.nutrient.data %>%
  group_by(HM) %>%
  mutate(adj.fiber = ifelse(is.na(fibre_adj), mean(na.omit(fibre_adj)), fibre_adj))
# double check

ggplot(oars.nutrient.data,
       aes(x=oars.days, 
           y=adj.fiber))+
  annotate("rect", xmin=0, xmax=max(subset(oars.nutrient.data, oars.on.rs=="onRS")$oars.days),
           ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_smooth(color="white",  se=T)+
  geom_point(aes(fill=oars.on.rs, shape=ifelse(oars.on.rs=="preRS", "1", "2")), 
             color="white", size=3)+
  scale_fill_manual(values=c("grey",2,"grey"))+
  scale_alpha_manual(values=c(1,0.2))+
  scale_shape_manual(values=c(23,21,21))+
  theme_classic()+theme(legend.position="none")+
  labs(x="Days since starting RS", 
       y="Fiber Intake (g per day)")


# add fiber intake to main mapping
metadata.oars.stool = merge(metadata.oars.stool, 
                            oars.nutrient.data[,c("standard.name", "adj.fiber")], 
                            by="standard.name", all.x = T)
# Good; done.

# :: save -----------------------------------------------------------------

write.csv(oars.nutrient.data, "./2025_07_24_oars_fiber.csv")

# And plots:

(metadata.oars.rs.plot+metadata.oars.rs.supp.plot+
    metadata.oars.fiber.plot+oars.nutrient.data.average.plot
) %>%
  ggsave(filename="./oars_plots/2026_01_13_oars_1_ffq_fiber_rs.pdf",
         width=10, height=6.5)

oars.nutrient.data.fiber.rs.plot%>%
  ggsave(filename="./oars_plots/2026_01_13_oars_1b_fiber_vs_rs.pdf",
         width=4, height=4)

# Brief methods write up:
# Fiber intakes were calculated from Monash FFQ and adjusted for calculated energy intake. 
# Values were calculated per stool. Missing values were imputed using the average of 2+ FFQs.

# Resistant starch doses & percentages were averaged per RS phase.
# Total energy was adjusted for the RS (2.7 kcal/g) and digestible starch (4 kcal/g)


# :: link to STL ----------------------------------------------------------

# check plot
ggplot()+
  # FFQ = circle
  geom_point(data = subset(oars.nutrient.data,!is.na(energy)), aes(x=oars.days, y=HM), shape=15, size=3, alpha=0.8, color="blue")+
  # STL = triangle
  geom_point(data = oars.nutrient.data, aes(x=oars.days, y=HM), shape=17, size=3, alpha=0.8, color="brown")+
  geom_vline(xintercept = 0)+
  theme_bw()+
  labs(x="Days since RS start", 
       y="")

# Goal: load 16S data and replace metadata
amplicon.data.gg.rare <- readRDS("./2025_09_10_rs_trial_16s_data_rarefied_gg2.Rds")
# subset to samples in OARS metadata
amplicon.data.gg.rare.oars = phyloseq::subset_samples(amplicon.data.gg.rare, 
                                                      standard.name %in% metadata.oars.stool$standard.name)
amplicon.data.gg.rare.oars.new.meta = merge(data.frame(phyloseq::sample_data(amplicon.data.gg.rare.oars)),
                                            metadata.oars.stool, by="standard.name")
amplicon.data.gg.rare.oars.new.meta$SampleID = gsub("\\.", "_", gsub("gg13v5.", "", amplicon.data.gg.rare.oars.new.meta$SampleID))
rownames(amplicon.data.gg.rare.oars.new.meta) = amplicon.data.gg.rare.oars.new.meta$SampleID
# merge in new metadata
amplicon.data.gg.rare.oars = phyloseq::merge_phyloseq(phyloseq::otu_table(amplicon.data.gg.rare.oars),
                                                      phyloseq::tax_table(amplicon.data.gg.rare.oars),
                                                      phyloseq::sample_data(amplicon.data.gg.rare.oars.new.meta))
# Good; 138 samples
oars.asv.meta = phyloseq::sample_data(amplicon.data.gg.rare.oars) %>% data.frame()
nrow(oars.asv.meta)

# collapse replicates to median (at ASV-level, for Alpha- and Beta-diversity)
oars.asv.data = speedyseq::psmelt(amplicon.data.gg.rare.oars)
oars.asv.data.median = reshape2::acast(oars.asv.data,
                                       standard.name ~ Taxa, value.var="Abundance", fun.aggregate=median)

# collapse replicates (and glom) to median (at LCA-level, for fold-change)
oars.asv.data.glom = oars.asv.data %>%
  group_by(Sample, LCA) %>%
  mutate(Abundance = sum(Abundance)) %>%
  dplyr::select(Sample, standard.name, LCA, Abundance) %>% distinct()
# then cast while taking median per standard.name
oars.asv.data.glom = reshape2::acast(oars.asv.data.glom,
                                     standard.name ~ LCA, value.var="Abundance", fun.aggregate=median)
# good

# examine missing stools
metadata.oars.stool$missing.16s = ifelse(!metadata.oars.stool$standard.name %in% data.frame(phyloseq::sample_data(amplicon.data.gg.rare.oars))$standard.name, "missing", "")
sum(metadata.oars.stool$missing.16s == "missing")
subset(metadata.oars.stool, missing.16s == "missing")
# HM0932-STL-08 is missing
# HM0924 provided 2 stools at 3M, the first of which is sequenced


# :: add stool water ----------------------------------------------------------

water.oars.stool <- read.csv("./2024_11_08_metabolomics_samples_prioritized_PD.csv")[,c("standard.name", "HM", "oars.days","stool_date_rec_v2", "stool_water_perc", "oars.timing")] %>% distinct()
# take unique
water.oars.stool = water.oars.stool[,c("standard.name", "stool_water_perc")] %>% distinct()

# append to original mapping, to ensure correct phase label
metadata.oars.stool = merge(metadata.oars.stool,
                            water.oars.stool, by="standard.name", all.x = T)


# :: add fecalcal ---------------------------------------------------------

# take from redcap (to ensure completeness)
redcap.data.loaded

oars.fcal = redcap.data.loaded[,c("standard.name", "fecalcal_res")]
colnames(oars.fcal)[2] = "fcal"

metadata.oars.stool = merge(metadata.oars.stool,
                            oars.fcal, by="standard.name", all.x=T) %>% distinct()


# >> Compliance -----------------------------------------------------------

metadata.oars.stool[,c("rs_notes","HM")] %>% distinct()

# these samples are "low-compliant" or were collected AFTER low-compliance
# (including their washout would not make sense)
non.compliant = c(
  # All of HM0639 (non-compliant overall)
  "HM0639-STL-05", # 0M; # recount indicates 50 / 138 missed doses; ~65%
  "HM0639-STL-06", # 3M; # poor compliance
  "HM0639-STL-07", # 6M; # poor compliance
  "HM0639-STL-08", # 9M; washout
  "HM0639-STL-09", # 12M; washout
  
  # After 3M for HM0924 (half-dose for re-opt)
  "HM0924-STL-13", # 6M; half doses for re-opt
  "HM0924-STL-14", # 9M; washout
  "HM0924-STL-15", # 12M; washout
  
  # After 3M for HM0932 (low compliance for last month)
  "HM0932-STL-07", # 6M; missed 50-75% of last month
  "HM0932-STL-08", # 9M; washout
  "HM0932-STL-09", # 12M; washout
  
  # All of HM0940 (low compliance)
  "HM0940-STL-12", # 0M; # recount indicates 60 / 120 missed doses; ~66%
  "HM0940-STL-13", # 3M; # poor compliance
  "HM0940-STL-14", # 6M; # poor compliance
  "HM0940-STL-15", # 9M; # washout
  "HM0940-STL-16", # 12M; # washout
  
  # After 3M for HM0759 (antibiotics right before 6M; then excluded)
  "HM0759-STL-07" # 6M; compliant but started antibiotics shortly before
)

metadata.oars.stool$compliant = ifelse(metadata.oars.stool$standard.name %in% non.compliant, FALSE, TRUE)

# add missing values, provided by Nathan
#HM0899-STL-03       Fecal Calprotectin result (ug/g) = 37.468
#HM0899-STL-05        Fecal Calprotectin result (ug/g) =   497
#HM0903-STL-10       Fecal Calprotectin result (ug/g) = 70
metadata.oars.stool$fcal = ifelse(metadata.oars.stool$standard.name == "HM0899-STL-03", "37",
                                      ifelse(metadata.oars.stool$standard.name == "HM0899-STL-05", "497",
                                             ifelse(metadata.oars.stool$standard.name == "HM0903-STL-10", "70",metadata.oars.stool$fcal)))
metadata.oars.stool$fcal = as.numeric(metadata.oars.stool$fcal)


# save this
saveRDS(metadata.oars.stool, "~/Documents/PhD/git_oars_archfolder/2025_07_22_oars_mapping_2.Rds")


# visualize stools + timings
ggplot(metadata.oars.stool,
       aes(x=oars.days, y=reorder(HM, compliant)))+
  annotate("rect", xmin=0, xmax=max(subset(metadata.oars.stool, oars.on.rs=="onRS")$oars.days),
           ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_point(aes(fill=oars.on.rs, 
                 alpha = ifelse(compliant == T, "1", "0")),
             shape=21, size=3, color="white")+
  scale_fill_manual(values=c("black", 2, "grey"))+
  scale_alpha_manual(values=c(0.25,1))+
  guides(alpha= "none")+
  #geom_text(aes(label=missing.asv), nudge_y=-0.15, size=6)+
  #geom_text(aes(label=substr(standard.name, nchar(standard.name)-2, nchar(standard.name))), size=2)+
  theme_classic()+theme(panel.grid.major.y=element_line(color="grey", linewidth=0.2))+
  labs(x="Day since starting RS", y="", fill="Phase", alpha="")
# Second HM0924 gets handled later

# :: 16S Alpha ------------------------------------------------------------

# calculate alpha measures
oars.asv.data.median.richness = vegan::specnumber(oars.asv.data.median)
oars.asv.data.median.shannon = vegan::diversity(oars.asv.data.median, index="shannon")

oars.asv.data.median.alpha = data.frame(
  standard.name = names(oars.asv.data.median.richness),
  richness = oars.asv.data.median.richness,
  shannon = oars.asv.data.median.shannon
)
# create new mapping file for Omics data
metadata.oars.stool.asv = merge(metadata.oars.stool,
                                oars.asv.data.median.alpha, by="standard.name")
# make rownames standard.name on both metadata files
rownames(metadata.oars.stool.asv) = metadata.oars.stool.asv$standard.name
rownames(metadata.oars.stool) = metadata.oars.stool$standard.name

nrow(metadata.oars.stool)
# 68 total samples
nrow(metadata.oars.stool.asv)
# 66 samples with ASV (2 samples missing 16S)

table(metadata.oars.stool.asv$HM)

# :: 16S Butyrogens I -----------------------------------------------------

# Use classic definition of butyrogens (Mottawea et al)
# Identify ASV's annotated as these butyrogens, and sum them in gg2 data
# reclassify using rarefied data and gg13.8

amplicon.data.gg138.tax.df = readRDS("./2025_09_10_rs_trial_16s_data_tax_table_gg2.Rds") # note: gg2 is a mis-type; it's actually 138 (see line 95)
amplicon.data.gg.rare <- readRDS("./2025_09_10_rs_trial_16s_data_rarefied_gg2.Rds") %>%
  phyloseq::subset_samples(standard.name %in% metadata.oars.stool$standard.name)
# subset to butyrogens
amplicon.data.gg138.tax.df = subset(data.frame(amplicon.data.gg138.tax.df), Family=="f__Lachnospiraceae" | Genus=="g__Blautia" | Genus=="g__Roseburia" | Genus=="g__Eubacterium" | Genus=="g__Ruminococcus" | Genus=="g__Clostridium" | Genus=="g__Faecalibacterium") # Note, Lachnospiraceae already includes several genera listed; listed again for clarity
# apply this filter to phyloseq data
oars.phyloseq.butyrogens.i = speedyseq::psmelt(amplicon.data.gg.rare)
oars.phyloseq.butyrogens.i = oars.phyloseq.butyrogens.i %>%
  # subset OTUs to butyrogens
  subset(OTU %in% amplicon.data.gg138.tax.df$ASV) %>%
  # for each sample, sum up all butyrogen reads (and divide by rarefaction depth for %)
  group_by(Sample) %>%
  mutate(but = sum(Abundance/50000)) %>%
  # then take median of replicates
  group_by(standard.name) %>%
  mutate(but = median(but)) %>%
  dplyr::select(study_id, standard.name, trial.stool.timing, fecalcal_res, but) %>% distinct() %>% data.frame()
oars.phyloseq.butyrogens.i %>% arrange(standard.name)

# add to mapping later
oars.phyloseq.butyrogens.i = oars.phyloseq.butyrogens.i[,c("standard.name", "but")]
colnames(oars.phyloseq.butyrogens.i) = c("standard.name", "but.i")


# :: 16S PICRUSt2 -------------------------------------------------------

# next, perform PICRUSt2 for Functional Redundancy and Kircher butyrogens

run.picrust = FALSE
if(run.picrust == T) {
# note: make new picrust2 instance
# use this data:
amplicon.data.gg.rare <- readRDS("./2025_09_10_rs_trial_16s_data_rarefied_gg2.Rds") %>%
  phyloseq::subset_samples(standard.name %in% metadata.oars.stool$standard.name)
# prepare seq.table for picrust2
fasta_seqs <- Biostrings::DNAStringSet(data.frame(phyloseq::refseq(amplicon.data.gg.rare))[,1])
names(fasta_seqs) <- rownames(data.frame(phyloseq::refseq(amplicon.data.gg.rare)))  # assign ASV IDs as sequence names
# prepare asv.table for picrust2
oars.picrust2.abuntable = data.frame(phyloseq::otu_table(amplicon.data.gg.rare))
# export
Biostrings::writeXStringSet(fasta_seqs, filepath = "oars_picrust2/oars.picrust2.seqtable.fasta")
Biostrings::writeXStringSet(fasta_seqs, filepath = "oars_picrust2/predict_SCFA_producers/oars.picrust2.seqtable.fasta")
write.table(t(oars.picrust2.abuntable), "oars_picrust2/oars.picrust2.abuntable.tsv", 
            sep="\t", quote = F, col.names = NA)

# STEP 1: run through default picrust2 (to ensure picrust2 works)
cd ~/Documents/PhD/git_oars_archfolder/oars_picrust2
conda activate oars_picrust2
  # In R, may need to install Rcpp, jsonlite, lattice, Matrix, RSpectra, castor
picrust2_pipeline.py \
-s oars.picrust2.seqtable.fasta \
-i oars.picrust2.abuntable.tsv \
-o oars_picrust2_out \
-p 1

# optionally perform stratification to source functions to taxa
# use EC's (not preferable; 32 million rows)
metagenome_pipeline.py -i oars.picrust2.abuntable.tsv -m oars_picrust2_out/marker_predicted_and_nsti.tsv.gz -f oars_picrust2_out/EC_predicted.tsv.gz \
-o oars_picrust2_out/EC_metagenome_out --strat_out

# :: 16S Butyrogens II (Vital) -------------------------------------------------------

# from: 

# STEP 2: run through Vitals' predict_SCFA_producers to identify butyrogens
# source: https://github.com/ag-vital/predict_SCFA_producers/tree/master
# first: git clone https://github.com/ag-vital/predict_SCFA_producers.git
# then, rename "picrust" folder as "SCFA"
cd ~/Documents/PhD/git_oars_archfolder/oars_picrust2/predict_SCFA_producers
conda activate oars_picrust2
# when rerunning, delete placement_working first
rm -r placement_working
place_seqs.py -s oars.picrust2.seqtable.fasta -o placed_seqs.tre -p 1 --intermediate placement_working --ref_dir SCFA
# 97 sequences failed to align
hsp.py -t placed_seqs.tre --observed_trait_table SCFA/SCFA_pathwaydata.txt -o SCFA_predicted.tsv -p 1 -m emp_prob -n

# save files to new folder
cp ./SCFA_predicted.tsv ../picrust2_saved/oars_SCFA_predicted.tsv
cp ../oars_picrust2_out/EC_metagenome_out/pred_metagenome_contrib.tsv.gz ../picrust2_saved/oars_pred_metagenome_contrib.tsv

}
# import predicted SCFA data
oars.vital.butyrogens = read.csv("./oars_picrust2/picrust2_saved/oars_SCFA_predicted.tsv", sep="\t")
# Butyrogens = "Possesses a) AcetylCoA pathway, and at least one of: b) but or c) buk"
oars.vital.butyrogens.acetylcoa = subset(oars.vital.butyrogens, acetylcoa == 1)
oars.vital.butyrogens.acetylcoa = subset(oars.vital.butyrogens.acetylcoa, but == 1 | buk == 1)
# these will be "butyrogens"
oars.vital.butyrogens.acetylcoa

# Subset phyloseq object to only include these ASVs
oars.phyloseq.butyrogens <- phyloseq::otu_table(amplicon.data.gg.rare)[,oars.vital.butyrogens.acetylcoa$sequence]
oars.phyloseq.butyrogens <- phyloseq::merge_phyloseq(oars.phyloseq.butyrogens, 
                                                     phyloseq::tax_table(amplicon.data.gg.rare), 
                                                     phyloseq::sample_data(amplicon.data.gg.rare))
# format to dataframe
oars.phyloseq.butyrogens = speedyseq::psmelt(oars.phyloseq.butyrogens)
# visualize to ensure these make sense
kircher.butyrogens.plot = oars.phyloseq.butyrogens %>%
  group_by(LCA) %>%
  mutate(sum.abun = sum(Abundance/50000)) %>% dplyr::select(LCA, sum.abun) %>% distinct() %>% data.frame() %>%
  slice_max(n=20, sum.abun) %>%
  ggplot(aes(x=reorder(LCA, sum.abun), y=sum.abun)) +
  coord_flip()+
  scale_y_log10()+
  geom_point(shape=21, aes(fill=scale(sum.abun)), size=2.5)+
  scale_fill_gradient2(low="blue", high="red")+
  facet_wrap(~"Kircher Butyrogens")+
  theme_classic()+theme(legend.position="none",
                        strip.text=element_text(size=10),
                        panel.grid.major.y = element_line(color="grey", linewidth=0.2),
                        strip.background = element_rect(color="black"))+
  labs(x="", y="Σ Relative Abundance")
kircher.butyrogens.plot
# top are Faecalibacterium, Gemmiger, Agathobacter, Anaerostipes, Anaerobutyricum, Roseburia
# save this plot in full data

# repeat But.i calculation (using new Butyrogens definition)
oars.phyloseq.butyrogens.df = oars.phyloseq.butyrogens %>%
  group_by(Sample) %>%
  mutate(but = sum(Abundance/50000)) %>%
  group_by(standard.name) %>%
  mutate(but.ii = median(but)) %>%
  dplyr::select(study_id, standard.name, trial.stool.timing, fecalcal_res, but.ii) %>% distinct() %>% data.frame()
# take median
oars.phyloseq.butyrogens.df %>% arrange(standard.name)

oars.phyloseq.butyrogens.ii = oars.phyloseq.butyrogens.df[,c("standard.name", "but.ii")]
colnames(oars.phyloseq.butyrogens.ii) = c("standard.name", "but.ii")


# :: 16S Functional Diversity ---------------------------------------------

run.fd = F
if(run.fd == T){
# using PICRUSt2
amplicon.data.gg.rare <- readRDS("./2025_09_10_rs_trial_16s_data_rarefied_gg2.Rds") # it shouldn't matter which version is used (gg2 or gg138), since we're taking ASV sequences

picrust_data <- read.csv("./oars_picrust2/picrust2_saved/oars_pred_metagenome_contrib.tsv", sep="\t")

dim(picrust_data)

# calculate redundancy per function
picrust_data_redun = picrust_data %>%
  group_by(sample, `function.`) %>%
  # ntaxa = how many ASVs possess this function
  mutate(ntaxa = length(unique(taxon))) %>%
  dplyr::select(sample, `function.`, ntaxa) %>% distinct()
  
# average functional redundancy (across functions) per sample
picrust_data_redun_mean <- picrust_data_redun %>%
  group_by(sample) %>%
  # nfunctions = count of unique functions per sample
  mutate(nfunctions = length(unique(`function.`))) %>%
  # meantaxa = mean number of taxa possessing this function
  mutate(meantaxa = mean(ntaxa)) %>%
  dplyr::select(sample, meantaxa, nfunctions) %>% distinct() %>%
  # fd = functional redundancy = ratio of taxa(function)
  # fr = functional richness
  mutate(fd = meantaxa) %>%
  mutate(fr = nfunctions) %>%
  dplyr::select(sample, fd, fr)

colnames(picrust_data_redun_mean)[1] = "dada2.sampleNames"

# link to STL and take median of replicates
oars.asv.functionalredundancy = merge(picrust_data_redun_mean,
                                data.frame(phyloseq::sample_data(amplicon.data.gg.rare))[,c("standard.name", "dada2.sampleNames")],
                                by="dada2.sampleNames")
oars.asv.functionalredundancy = oars.asv.functionalredundancy %>%
  group_by(standard.name) %>%
  mutate(fd = median(fd)) %>%
  mutate(fr = median(fr)) %>%
  dplyr::select(standard.name, fd, fr) %>% distinct() %>% data.frame()

saveRDS(oars.asv.functionalredundancy, "./2026_01_13_oars_functional_redundancy.Rds")

}

oars.asv.functionalredundancy = readRDS("./2026_01_13_oars_functional_redundancy.Rds")

# quickly check correlation
ggplot(oars.asv.functionalredundancy,
       aes(x=fd, y=fr))+
  geom_point(shape=21, color="white", fill="black", size=2.5)+
  geom_smooth(method="lm", color="white", fill="black")+
  ggpubr::stat_cor(method="spearman")+
  theme_classic()+
  labs(x="Functional Redundancy",
       y="Functional Richness")

# fd = functional redundancy
# fr = functional richness

# :: 16S MLP --------------------------------------------------------------

# Apply (but also, validate) Bork's Microbial Load Predictor

# load RDP-annotated ASV data
amplicon.data.rdp.tax.df = readRDS("./2025_09_10_rs_trial_16s_data_tax_table_rdp.Rds")
amplicon.data.rdp.tax.df$OTU = amplicon.data.rdp.tax.df$ASV

# collapse to median
oars.asv.data.for.mlp = speedyseq::psmelt(phyloseq::subset_samples(amplicon.data.gg.rare, standard.name %in%metadata.oars.stool$standard.name))
oars.asv.data.for.mlp = oars.asv.data.for.mlp[,c("OTU", "Abundance", "standard.name")]
oars.asv.data.for.mlp$Abundance = oars.asv.data.for.mlp$Abundance / 50000
# take median % of ASVs
oars.asv.data.for.mlp = oars.asv.data.for.mlp %>%
  group_by(standard.name, OTU) %>%
  mutate(Abundance = median(Abundance)) %>% distinct()
oars.asv.data.for.mlp = merge(oars.asv.data.for.mlp[,c("OTU", "standard.name", "Abundance")],
                              amplicon.data.rdp.tax.df, by="OTU") #

# format OARS taxa (e.g. all genera, or uc_o/c/f) to match MLP
oars.asv.data.for.mlp = oars.asv.data.for.mlp %>%
  mutate(RDP_taxa = ifelse(!is.na(Genus), paste("g_", Genus,sep=""),
                           ifelse(!is.na(Family), paste("uc_f_", Family, sep=""),
                                  ifelse(!is.na(Order), paste("uc_o_", Order, sep=""),
                                         ifelse(!is.na(Class), paste("uc_c_", Class, sep=""),
                                                ifelse(!is.na(Phylum), paste("uc_p_", Phylum, sep=""),
                                                       "Bacteria"))))))

# make matrix; sum up glom'ed taxa using RDP annotation (range = 0-73.6%)
oars.asv.data.for.mlp = reshape2::acast(oars.asv.data.for.mlp,
                                        standard.name ~ RDP_taxa, value.var="Abundance",
                                        fun.aggregate=sum) %>% as.data.frame()
# Note: I'm doing this in a different order than other 16S processing
# e.g. here: median --> sum (if median is 0, lose that taxa)
# elsewhere: sum --> median (if median is 0, sum may rescue detection)

# add Shannon
oars.asv.data.for.mlp.shannon = vegan::diversity(oars.asv.data.for.mlp)
# pseudocount
asv.mlp.min = min(oars.asv.data.for.mlp[oars.asv.data.for.mlp!=0])/2

# log transform with pseudocount later* (because we need to use the smallest pseudo)

# re-train MLP model
mlp.rdp = readRDS("./oars_mlp/model.16S_rRNA.rds")
# extract training data
mlp.rdp.original.data = mlp.rdp$trainingData
# and, we'll use the source data to reprocess
mlp.rdp.retrain.data = read.csv("./oars_mlp/2025_03_05_vandeputte_otu_unrarefied.csv")
mlp.rdp.retrain.data.load = read.csv("./oars_mlp/2025_03_05_vandeputte_load.csv")
# 1. match samples (add rowname, so we can match later, after rarefying)
rownames(mlp.rdp.retrain.data) = paste("sample", mlp.rdp.retrain.data$X)
rownames(mlp.rdp.retrain.data.load) = paste("sample", mlp.rdp.retrain.data.load$X)
# 2. rarefy data
mlp.rdp.retrain.data = mlp.rdp.retrain.data[!is.na(rowSums(mlp.rdp.retrain.data)),]
set.seed(25)
mlp.rdp.retrain.data = phyloseq::rarefy_even_depth(phyloseq::otu_table(mlp.rdp.retrain.data,taxa_are_rows=F), 
                                                   sample.size=10000, replace=F)
mlp.rdp.retrain.data = mlp.rdp.retrain.data / 10000
mlp.rdp.retrain.data = as.data.frame(mlp.rdp.retrain.data)
# add Shannon
mlp.rdp.retrain.data.shannon = vegan::diversity(mlp.rdp.retrain.data)

# note: my data has 212 taxa; theirs is 194; comparable

# 3. log2 transform + pseudo
mlp.rdp.retrain.data = as.data.frame(mlp.rdp.retrain.data)
mlp.rdp.pseudo = min(mlp.rdp.retrain.data[mlp.rdp.retrain.data!=0])/2
if(asv.mlp.min < mlp.rdp.pseudo){
  asv.mlp.min.to.use = asv.mlp.min
  print(paste("Use OARS pseudocount:", asv.mlp.min.to.use))
}else{
  asv.mlp.min.to.use = mlp.rdp.pseudo
  print(paste("Use MLP pseudocount:", asv.mlp.min.to.use))
}
# since my pseudo is smaller, we'll use it instead
# --
# briefly return to OARS data
oars.asv.data.for.mlp = log2(oars.asv.data.for.mlp+asv.mlp.min.to.use)
oars.asv.data.for.mlp$Shannon = oars.asv.data.for.mlp.shannon
# --
# now back to MLP data
mlp.rdp.retrain.data = log2(mlp.rdp.retrain.data+asv.mlp.min.to.use)
mlp.rdp.retrain.data$Shannon = mlp.rdp.retrain.data.shannon
# 4. add load data
mlp.rdp.retrain.data = merge(mlp.rdp.retrain.data,
                             mlp.rdp.retrain.data.load[,c("Cell_count_per_gram", "X")], by="row.names")
mlp.rdp.retrain.data$X.x = NULL
mlp.rdp.retrain.data$X.y = NULL
rownames(mlp.rdp.retrain.data) = mlp.rdp.retrain.data$Row.names
mlp.rdp.retrain.data$Row.names = NULL
mlp.rdp.retrain.data.load = log10(mlp.rdp.retrain.data$Cell_count_per_gram) # Note: Log10
mlp.rdp.retrain.data$Cell_count_per_gram = NULL
# 5. take overlapping features
mlp.rdp.retrain.data = mlp.rdp.retrain.data[,colnames(mlp.rdp.retrain.data) %in% colnames(oars.asv.data.for.mlp)]
mlp.rdp.retrain.data$load = mlp.rdp.retrain.data.load
dim(mlp.rdp.retrain.data) # 133 taxa
mlp.rdp.retrain.data = subset(mlp.rdp.retrain.data, !is.na(load))

# 6. build/validate model

# I will re-validate their results using their validation data.
# Perhaps try different models (glmnet seems to work better, for instance)

t1 = Sys.time()
set.seed(25)
mlp.rdp.model.xgb = caret::train(load ~.,
                             mlp.rdp.retrain.data,
                             method = "xgbTree",
                             trControl = caret::trainControl(method = "cv", number = 10,
                                                             savePredictions=TRUE),
                             verbose = FALSE
)
t2 = Sys.time()
t2-t1 # 2 min

# compare correlation with original
mlp.rdp.model.xgb.preds = subset(mlp.rdp.model.xgb$pred, 
                             eta == mlp.rdp.model.xgb$finalModel$tuneValue$eta &
                               gamma  ==  mlp.rdp.model.xgb$finalModel$tuneValue$gamma &
                               nrounds  ==  mlp.rdp.model.xgb$finalModel$tuneValue$nrounds &
                               max_depth  ==  mlp.rdp.model.xgb$finalModel$tuneValue$max_depth &
                               colsample_bytree ==  mlp.rdp.model.xgb$finalModel$tuneValue$colsample_bytree &
                               min_child_weight == mlp.rdp.model.xgb$finalModel$tuneValue$min_child_weight &
                               subsample == mlp.rdp.model.xgb$finalModel$tuneValue$subsample)[,c("pred", "obs")]
cor.test(mlp.rdp.model.xgb.preds$pred, mlp.rdp.model.xgb.preds$obs, method="pearson")
mlp.rdp.retrain.xgb.plot = ggplot(mlp.rdp.model.xgb.preds,
       aes(x=pred, y=obs))+
  geom_point(shape=21)+geom_smooth(method="lm")+
  ggpubr::stat_cor(method="pearson", size=3)+theme_classic()+
  facet_wrap(~"Retrained")+labs(x="Prediction", y="Observed")
mlp.rdp.retrain.xgb.plot
# Pearson Cor = 0.79

mlp.rdp.original.preds = subset(mlp.rdp$pred, 
                                eta == mlp.rdp$finalModel$tuneValue$eta &
                                  gamma  ==  mlp.rdp$finalModel$tuneValue$gamma &
                                  nrounds  ==  mlp.rdp$finalModel$tuneValue$nrounds &
                                  max_depth  ==  mlp.rdp$finalModel$tuneValue$max_depth &
                                  colsample_bytree ==  mlp.rdp$finalModel$tuneValue$colsample_bytree &
                                  min_child_weight == mlp.rdp$finalModel$tuneValue$min_child_weight &
                                  subsample == mlp.rdp$finalModel$tuneValue$subsample)[,c("pred", "obs", "rowIndex")] %>%
  # need to take average of preds across 10x
  group_by(rowIndex) %>%
  mutate(pred = mean(pred)) %>%
  distinct() %>% data.frame()
cor.test(mlp.rdp.original.preds$pred, mlp.rdp.original.preds$obs, method="pearson")
mlp.rdp.original.plot = ggplot(mlp.rdp.original.preds,
       aes(x=pred, y=obs))+
  geom_point(shape=21, color="white", fill="black")+geom_smooth(method="lm", color="black")+
  ggpubr::stat_cor(method="pearson", size=3)+theme_classic()+
  facet_wrap(~"Original")+labs(x="Prediction", y="Observed")
mlp.rdp.original.plot
# Cor = 0.79
# Retrained is the same! Let's use it


# apply to OARS
oars.asv.mlp.preds.xgb = predict(mlp.rdp.model.xgb, oars.asv.data.for.mlp)
oars.asv.mlp.preds.xgb = data.frame(standard.name = gsub("\\.", "-", rownames(oars.asv.data.for.mlp)),
                                load = oars.asv.mlp.preds.xgb)

oars.asv.mlp.preds.xgb.plot = ggplot(oars.asv.mlp.preds.xgb,
       aes(x=load, y=1))+
  ggridges::geom_density_ridges2()+
  theme_classic()+
  facet_wrap(~"Retrained XGBoost")+
  labs(x="Predicted Load", y="Density")
oars.asv.mlp.preds.xgb.plot
# 
# let's assess other models 

# :: 16S MLP Validation ---------------------------------------------------

# select models
models = c("glmnet", "svmRadial", "xgbTree", "ranger", "pls")

mlp.rdp.retrain.models = do.call(rbind, lapply(models, function(model){
  print(model)
  t1 = Sys.time()
  set.seed(25)
  mlp.rdp.model = caret::train(load ~.,
                               mlp.rdp.retrain.data,
                               method = model,
                               trControl = caret::trainControl(method = "cv", number = 10,
                                                               savePredictions=TRUE),
                               verbose = FALSE,
                               preProcess = c("center", "scale")
  )
  t2 = Sys.time()
  t2-t1 # 2 min
  
  # extract preds per model
  if(model == "glmnet"){
    mlp.rdp.model.preds = subset(mlp.rdp.model$pred, alpha == mlp.rdp.model$bestTune$alpha &
                                   lambda == mlp.rdp.model$bestTune$lambda)[,c("rowIndex", "pred", "obs")] %>% data.frame()
  }
  if(model == "svmRadial"){
    mlp.rdp.model.preds = subset(mlp.rdp.model$pred, sigma == mlp.rdp.model$bestTune$sigma &
                                   C == mlp.rdp.model$bestTune$C)[,c("rowIndex", "pred", "obs")] %>% data.frame()
  }
  if(model == "xgbTree"){
    mlp.rdp.model.preds = subset(mlp.rdp.model$pred, eta == mlp.rdp.model$bestTune$eta &
                                   gamma  ==  mlp.rdp.model$bestTune$gamma &
                                   nrounds  ==  mlp.rdp.model$bestTune$nrounds &
                                   max_depth  ==  mlp.rdp.model$bestTune$max_depth &
                                   colsample_bytree ==  mlp.rdp.model$bestTune$colsample_bytree &
                                   min_child_weight == mlp.rdp.model$bestTune$min_child_weight &
                                   subsample == mlp.rdp.model$bestTune$subsample)[,c("rowIndex", "pred", "obs")] %>% data.frame()
  }
  if(model == "ranger"){
    mlp.rdp.model.preds = subset(mlp.rdp.model$pred, mtry == mlp.rdp.model$bestTune$mtry &
                                   splitrule == mlp.rdp.model$bestTune$splitrule &
                                   min.node.size == mlp.rdp.model$bestTune$min.node.size)[,c("rowIndex", "pred", "obs")] %>% data.frame()
  }
  if(model == "pls"){
    mlp.rdp.model.preds = subset(mlp.rdp.model$pred, ncomp == mlp.rdp.model$bestTune$ncomp)[,c("rowIndex", "pred", "obs")] %>% data.frame()
  }
  # need to take average of preds across 10x
  
  mlp.rdp.model.preds = mlp.rdp.model.preds %>%
    group_by(rowIndex) %>%
    mutate(pred = mean(pred)) %>%
    distinct() %>% data.frame()
  mlp.rdp.model.preds$model = model
  mlp.rdp.model.preds$time = t2-t1
  mlp.rdp.model.preds
}))

mlp.rdp.retrain.models.plots = ggplot(mlp.rdp.retrain.models %>%
                                        rbind(data.frame(mlp.rdp.original.preds[,c("rowIndex", "pred", "obs")]) %>%
                                              mutate(model = "Original", time = NA)) %>% data.frame() %>%
                                        mutate(model = factor(model, levels=c("Original", "glmnet","pls","svmRadial","ranger","xgbTree"))),
       aes(x=pred, y=obs))+
  geom_point(shape=21, color="white", fill="black")+geom_smooth(method="lm", color="black")+
  ggpubr::stat_cor(method="pearson", size=3)+theme_classic()+
  facet_wrap(~model)+labs(x="Prediction", y="Observed")
mlp.rdp.retrain.models.plots
# RandomForest is comparable

# Correlate predictions
Hmisc::rcorr(reshape2::acast(mlp.rdp.retrain.models,
                               model ~ rowIndex, value.var="pred")%>%t() ,
               type="pearson")$r %>% pheatmap::pheatmap(color=colorRampPalette(c("blue","white", "red"))(100))
# indeed, the most similar to xgBoost predictions

# build final model (RandomForest) and apply
set.seed(25)
mlp.rdp.model.final = caret::train(load ~.,
                             mlp.rdp.retrain.data,
                             method = "ranger",
                             trControl = caret::trainControl(method = "cv", number = 10,
                                                             savePredictions=TRUE),
                             verbose = FALSE,
                             preProcess = c("center", "scale")
)

# apply to OARS
oars.asv.mlp.preds = predict(mlp.rdp.model.final, oars.asv.data.for.mlp)
oars.asv.mlp.preds = data.frame(standard.name = gsub("\\.", "-", rownames(oars.asv.data.for.mlp)),
                                load.asv = oars.asv.mlp.preds)

oars.asv.mlp.preds.rf.plot = ggplot(oars.asv.mlp.preds,
       aes(x=load.asv, y=1))+
  ggridges::geom_density_ridges2()+
  theme_classic()+
  facet_wrap(~"Retrained Ranger")+
  labs(x="Predicted Load", y="Density")
oars.asv.mlp.preds.rf.plot
# 



# >> MGX ------------------------------------------------------------------

# make new meta
metadata.oars.stool.mgx = metadata.oars.stool # 68 total samples

# :: load MGX ----------------------------------------------------------

# load manifest (list of files and their MGX completeness)
mgx.oars.manifest = read.csv("~/Downloads/humann3_outputs_summary_250502.csv")

# subset manifest to OARS data samples
mgx.oars.manifest = subset(mgx.oars.manifest, standard.name %in% metadata.oars.stool.mgx$standard.name)

# create variable based on MGX completion
mgx.oars.manifest$mgx.tax.complete = ifelse(mgx.oars.manifest$standard.name %in% subset(mgx.oars.manifest, metaphlan.bugs == "completed")$standard.name, "yes", "")
mgx.oars.manifest$mgx.path.complete = ifelse(mgx.oars.manifest$standard.name %in% subset(mgx.oars.manifest, pathabundance == "completed")$standard.name, "yes", "")
mgx.oars.manifest$mgx.sequenced = ifelse(mgx.oars.manifest$standard.name %in% subset(mgx.oars.manifest, !grepl("failed", sample.id))$standard.name, "yes", "")

library("tidyverse")
# list all files
all_files <- list.files(path = "~/Downloads/humann3_main_outputs/metaphlan_bugs_list", full.names = TRUE)
# reduce to desired samples
desired_samples = subset(mgx.oars.manifest, mgx.tax.complete == "yes")$standard.name
file_list <- all_files[stringr::str_detect(basename(all_files), 
                                           paste(desired_samples, collapse = "|"))]
# write function to read and process MGX files
read_metaphlan_file <- function(file_path) {
  # read the file
  data <- read.csv(file_path, sep="\t", header=F) # 
  # remove first few rows (blank + titles); and select columns (NCBI tax ID + "additional species")
  data = data[-c(1:5),c(1,3)]
  # add colnames
  colnames(data) = c("taxa", "abundance")
  data$sample = paste("HM", gsub(".*HM", "", gsub("\\_merged.*", "", file_path)), sep="")
  # optionally filter out to just species level
  data = subset(data, grepl("s__", taxa))
  data = subset(data, !grepl("t__", taxa))
  data$taxa = gsub(".*s__", "", data$taxa) # clean taxa name to only list species
  # convert abundance to numeric
  data$abundance = as.numeric(data$abundance)
  
  return(data)
}
# apply function to data
oars.mgx.taxa <- do.call(rbind, lapply(file_list, read_metaphlan_file))
# make matrix (and take average of repeat sequenced samples)
oars.mgx.taxa = reshape2::acast(oars.mgx.taxa, sample ~ taxa, 
                                value.var="abundance", fun.aggregate=mean)
# replace NA with 0
oars.mgx.taxa[is.na(oars.mgx.taxa)] = 0

dim(oars.mgx.taxa) # ~1000 species; 66 samples

# :: MGX MLP --------------------------------------------------------------

# now apply MLP to unfiltered MGX data
oars.mgx.data.for.mlp = as.data.frame(oars.mgx.taxa)
# rescale to %
oars.mgx.data.for.mlp <- sweep(oars.mgx.data.for.mlp, 1, rowSums(oars.mgx.data.for.mlp), FUN = "/") * 1
range(oars.mgx.data.for.mlp)

# calculate Shannon diversity
oars.mgx.data.for.mlp.shannon = vegan::diversity((oars.mgx.data.for.mlp), index="shannon")

#library("MLP")
# extract MGX MLP data
mlp.mp3 = readRDS("./oars_mlp/model.metaphlan3.rds")
mlp.mp3.data = mlp.mp3$trainingData
# undo log transform (they say log10 in their methods)
mlp.mp3.data = 10^mlp.mp3.data
mlp.mp3.data.load = log10(mlp.mp3.data$.outcome) # undo the log on load column
mlp.mp3.data.shannon = log10(mlp.mp3.data$`Shannon diversity`) # undo the log on shannon column and save
mlp.mp3.data$`.outcome` = NULL
mlp.mp3.data$`Shannon diversity` = NULL

# replace min value with 0 (per taxa); but ignore Shannon and Load columns
mlp.mp3.data <- do.call(cbind, lapply(1:ncol(mlp.mp3.data), function(x) {
  col = data.frame(mlp.mp3.data[,x])
  col[col == min(col, na.rm = TRUE)] <- 0
  # and convert to %
  col = col/100
  col.df = data.frame(feature = col)
  colnames(col.df) = colnames(mlp.mp3.data)[x]
  col.df
})) %>% data.frame()
rowSums(mlp.mp3.data)
# these values are way off of what'd be expected; unlikely to want to use MGX

# rescale to %
mlp.mp3.data <- sweep(mlp.mp3.data, 1, rowSums(mlp.mp3.data), FUN = "/") * 1

# clean up taxa names
colnames(mlp.mp3.data) = gsub(".*s__", "", colnames(mlp.mp3.data))

# take only overlapping taxa
dim(mlp.mp3.data) # 156 taxa
mlp.mp3.data = mlp.mp3.data[,colnames(mlp.mp3.data) %in% c(colnames(oars.mgx.data.for.mlp), ".outcome")] 
dim(mlp.mp3.data) # 90 taxa (92 - 2)

# log2 transform + pseudocount
mgx.mlp.min = min(mlp.mp3.data[mlp.mp3.data!=0])/2
mgx.oars.mlp.min = min(oars.mgx.data.for.mlp[oars.mgx.data.for.mlp!=0])/2
if(mgx.oars.mlp.min < mgx.mlp.min){
  mgx.mlp.min.to.use = mgx.oars.mlp.min
  print(paste("Use OARS pseudocount:", mgx.mlp.min.to.use))
}else{
  mgx.mlp.min.to.use = mgx.mlp.min
  print(paste("Use MLP pseudocount:", mgx.mlp.min.to.use))
}
mlp.mp3.data = log2(mlp.mp3.data+mgx.mlp.min.to.use)
# add Shannon and Load cols back
mlp.mp3.data$load = mlp.mp3.data.load
mlp.mp3.data$Shannon = mlp.mp3.data.shannon
# and finish off OARS data
oars.mgx.data.for.mlp = log2(oars.mgx.data.for.mlp+mgx.mlp.min.to.use)
oars.mgx.data.for.mlp$Shannon = oars.mgx.data.for.mlp.shannon

# train model
t1 = Sys.time()
mlp.mp3.model.xgb = caret::train(load ~.,
             mlp.mp3.data,
             method = "xgbTree",
             trControl = caret::trainControl(method = "cv", number = 10, # note: they perform 5x 10x
                                             savePredictions=TRUE),
             verbose = FALSE
             )
t2 = Sys.time()
t2-t1 # 5 min

# compare correlation with original
mlp.mp3.model.xgb.preds = subset(mlp.mp3.model.xgb$pred, 
                             eta == mlp.mp3.model.xgb$finalModel$tuneValue$eta &
                               gamma  ==  mlp.mp3.model.xgb$finalModel$tuneValue$gamma &
                               nrounds  ==  mlp.mp3.model.xgb$finalModel$tuneValue$nrounds &
                               max_depth  ==  mlp.mp3.model.xgb$finalModel$tuneValue$max_depth &
                               colsample_bytree ==  mlp.mp3.model.xgb$finalModel$tuneValue$colsample_bytree &
                               min_child_weight == mlp.mp3.model.xgb$finalModel$tuneValue$min_child_weight &
                               subsample == mlp.mp3.model.xgb$finalModel$tuneValue$subsample)[,c("pred", "obs")]
cor.test(mlp.mp3.model.xgb.preds$pred, mlp.mp3.model.xgb.preds$obs, method="pearson")

mlp.mgx.retrain.xgb.plot = ggplot(mlp.mp3.model.xgb.preds,
                               aes(x=pred, y=obs))+
  geom_point(shape=21, color="white", fill="black")+geom_smooth(method="lm", color="black")+
  ggpubr::stat_cor(method="pearson", size=3)+theme_classic()+
  facet_wrap(~"Retrained")+labs(x="Prediction", y="Observed")
mlp.mgx.retrain.xgb.plot# from 0.61 to 0.54 after re-processing (including reducing features)

mlp.mp3.original.preds = subset(mlp.mp3$pred, 
                                eta == mlp.mp3$finalModel$tuneValue$eta &
                                  gamma  ==  mlp.mp3$finalModel$tuneValue$gamma &
                                  nrounds  ==  mlp.mp3$finalModel$tuneValue$nrounds &
                                  max_depth  ==  mlp.mp3$finalModel$tuneValue$max_depth &
                                  colsample_bytree ==  mlp.mp3$finalModel$tuneValue$colsample_bytree &
                                  min_child_weight == mlp.mp3$finalModel$tuneValue$min_child_weight &
                                  subsample == mlp.mp3$finalModel$tuneValue$subsample)[,c("pred", "obs", "rowIndex")] %>%
  # need to take average of preds across 10x
  group_by(rowIndex) %>%
  mutate(pred = mean(pred)) %>%
  distinct() %>% data.frame()
cor.test(mlp.mp3.original.preds$pred, mlp.mp3.original.preds$obs, method="pearson")
mlp.mgx.original.plot = ggplot(mlp.mp3.original.preds,
                               aes(x=pred, y=obs))+
  geom_point(shape=21, color="white", fill="black")+geom_smooth(method="lm", color="black")+
  ggpubr::stat_cor(method="pearson", size=3)+theme_classic()+
  facet_wrap(~"Original")+labs(x="Prediction", y="Observed")
mlp.mgx.original.plot

# apply to OARS

oars.mgx.mlp.preds = predict(mlp.mp3.model.xgb, oars.mgx.data.for.mlp)
oars.mgx.mlp.preds = data.frame(standard.name = gsub("\\.", "-", rownames(oars.mgx.data.for.mlp)),
           load.mgx = oars.mgx.mlp.preds)

hist(oars.mgx.mlp.preds$load.mgx)
# these look more continuous than the ASV predictions; perhaps xgBoost is fine for MGX

# compare to 16S MLP results


# :: MGX MLP Validation ---------------------------------------------------

# Like with ASVs, validate other models for MGX MLP

# select models
models = c("glmnet", "svmRadial", "xgbTree", "ranger", "pls")

mlp.mgx.retrain.models = do.call(rbind, lapply(models, function(model){
  print(model)
  t1 = Sys.time()
  set.seed(25)
  mlp.mgx.model = caret::train(load ~.,
                               mlp.mp3.data %>% as.data.frame(),
                               method = model,
                               trControl = caret::trainControl(method = "cv", number = 10,
                                                               savePredictions=TRUE),
                               verbose = FALSE,
                               preProcess = c("center", "scale")
  )
  t2 = Sys.time()
  t2-t1 # 2 min
  
  # extract preds per model
  if(model == "glmnet"){
    mlp.mgx.model.preds = subset(mlp.mgx.model$pred, alpha == mlp.mgx.model$bestTune$alpha &
                                   lambda == mlp.mgx.model$bestTune$lambda)[,c("rowIndex", "pred", "obs")] %>% data.frame()
  }
  if(model == "svmRadial"){
    mlp.mgx.model.preds = subset(mlp.mgx.model$pred, sigma == mlp.mgx.model$bestTune$sigma &
                                   C == mlp.mgx.model$bestTune$C)[,c("rowIndex", "pred", "obs")] %>% data.frame()
  }
  if(model == "xgbTree"){
    mlp.mgx.model.preds = subset(mlp.mgx.model$pred, eta == mlp.mgx.model$finalModel$tuneValue$eta &
                                   gamma  ==  mlp.mgx.model$finalModel$tuneValue$gamma &
                                   nrounds  ==  mlp.mgx.model$finalModel$tuneValue$nrounds &
                                   max_depth  ==  mlp.mgx.model$finalModel$tuneValue$max_depth &
                                   colsample_bytree ==  mlp.mgx.model$finalModel$tuneValue$colsample_bytree &
                                   min_child_weight == mlp.mgx.model$finalModel$tuneValue$min_child_weight &
                                   subsample == mlp.mgx.model$finalModel$tuneValue$subsample)[,c("rowIndex", "pred", "obs")] %>% data.frame()
  }
  if(model == "ranger"){
    mlp.mgx.model.preds = subset(mlp.mgx.model$pred, mtry == mlp.mgx.model$bestTune$mtry &
                                   splitrule == mlp.mgx.model$bestTune$splitrule &
                                   min.node.size == mlp.mgx.model$bestTune$min.node.size)[,c("rowIndex", "pred", "obs")] %>% data.frame()
  }
  if(model == "pls"){
    mlp.mgx.model.preds = subset(mlp.mgx.model$pred, ncomp == mlp.mgx.model$bestTune$ncomp)[,c("rowIndex", "pred", "obs")] %>% data.frame()
  }
  # need to take average of preds across 10x
  
  mlp.mgx.model.preds = mlp.mgx.model.preds %>%
    group_by(rowIndex) %>%
    mutate(pred = mean(pred)) %>%
    distinct() %>% data.frame()
  mlp.mgx.model.preds$model = model
  mlp.mgx.model.preds$time = t2-t1
  mlp.mgx.model.preds
}))

mlp.mgx.retrain.models.plots = ggplot(mlp.mgx.retrain.models %>%
                                        rbind(data.frame(mlp.mp3.original.preds[,c("rowIndex", "pred", "obs")]) %>%
                                                mutate(model = "Original", time = NA)) %>% data.frame() %>%
                                        mutate(model = factor(model, levels=c("Original", "glmnet","pls","svmRadial","ranger","xgbTree"))),
                                      aes(x=pred, y=obs))+
  geom_point(shape=21, color="white", fill="black")+geom_smooth(method="lm", color="black")+
  ggpubr::stat_cor(method="pearson", size=3)+theme_classic()+
  facet_wrap(~model)+labs(x="Prediction", y="Observed")
mlp.mgx.retrain.models.plots

# Correlate predictions
Hmisc::rcorr(reshape2::acast(mlp.mgx.retrain.models,
                             model ~ rowIndex, value.var="pred")%>%t() ,
             type="pearson")$r %>% pheatmap::pheatmap(color=colorRampPalette(c("blue","white", "red"))(100))
# ranger again is most comparable


# :: ```plot: MLP comparison -------------------------------------------------

# comparison of predictions
mlp.rdp.mgx.comparison.plot = ggplot(merge(oars.asv.mlp.preds,
             oars.mgx.mlp.preds, by="standard.name"),
       aes(x=load.asv, y=load.mgx))+
  geom_point(shape=21)+geom_smooth(method="lm")+
  ggpubr::stat_cor(method="pearson", size=3)+theme_classic()+
  facet_wrap(~"ASV vs MGX")+labs(x="ASV Load", y="MGX Load")

# ASV performs better, but they do not correlate
library("patchwork")
(mlp.rdp.retrain.models.plots / 
                      mlp.mgx.retrain.models.plots) %>%
  ggsave(filename="./oars_plots/oars_1b_mlp_asv_mgx.pdf",
         width=8, height=8,device = cairo_pdf)

# plus
(oars.asv.mlp.preds.xgb.plot/
  oars.asv.mlp.preds.rf.plot/
  mlp.rdp.mgx.comparison.plot)%>%
  ggsave(filename="./oars_plots/oars_1b_mlp_retrain.pdf",
         width=4, height=8,device = cairo_pdf)
# plus
Hmisc::rcorr(reshape2::acast(mlp.mgx.retrain.models,
                             model ~ rowIndex, value.var="pred")%>%t() ,
             type="pearson")$r %>% pheatmap::pheatmap(color=colorRampPalette(c("blue","white", "red"))(100))




oars.asv.mlp.preds
oars.mgx.mlp.preds

# >> MPX ----------------------------------------------------------------

oars.mpx = read.csv("./metaproteomics/MetaLab_iterative/final_proteins.tsv", sep="\t")
dim(oars.mpx)
# note: these include LSARP samples

# add plate number to samples in pro data
# 1), remove Intensity
colnames(oars.mpx) = gsub("Intensity\\.", "", colnames(oars.mpx))
# 2) remove Sara_Stintzi
colnames(oars.mpx) = gsub("Sara_Stintzi", "", colnames(oars.mpx))
colnames(oars.mpx) = gsub("Lab", "", colnames(oars.mpx))
# 3) subset to important samples
oars.mpx = oars.mpx[,grepl(paste(c("Protein.IDs", "HM", "QC"), collapse="|"),colnames(oars.mpx))]


# save mapping file (with plate locations), but needs cleaning
oars.mpx.map = read.csv("./metaproteomics/metaproteomics_mapping.csv")
# do the same for mapping
oars.mpx.map = oars.mpx.map %>%
  mutate(clean = gsub("Intensity\\.", "",
                      gsub("Sara_Stintzi", "", gsub("Lab", "", sample)))) 
oars.mpx.map = oars.mpx.map %>%
  tidyr::separate(clean, sep="__", into=c("date", "standard.name"))
# replace _ with -
oars.mpx.map$standard.name = gsub("_", "-", oars.mpx.map$standard.name)
# use this map to test in analysis
# and add extra consideration for seq replicates
oars.mpx.map.plates = oars.mpx.map %>%
  group_by(standard.name) %>%
  #subset(standard.name == "HM0878-STL-11")%>%
  mutate(plate = stringr::str_c(plate, collapse = "_")) %>%
  dplyr::select(plate, standard.name) %>% distinct()

# >> MPX QC ---------------------------------------------------------------

# check QC samples
oars.mpx.qc = oars.mpx[,!grepl("HM", colnames(oars.mpx))]
colnames(oars.mpx.qc)
rownames(oars.mpx.qc) = oars.mpx.qc[,1]
oars.mpx.qc[,1] = NULL

# (log2 + pseudo)
oars.mpx.qc.log = t(log2(oars.mpx.qc+(min(oars.mpx.qc[oars.mpx.qc!=0])/2)))
# filter out 0 variance
oars.mpx.qc.log = oars.mpx.qc.log[,apply(oars.mpx.qc.log, 2, sd)>0]
# pca
oars.mpx.qc.pca = prcomp(oars.mpx.qc.log,
                          scale=T, center=T)
oars.mpx.qc.pca.df = data.frame(oars.mpx.qc.pca$x[,c(1:2)])
oars.mpx.qc.pca.var = (oars.mpx.qc.pca$sdev)^2 / sum(oars.mpx.qc.pca$sdev^2) * 100

# plot
ggplot(oars.mpx.qc.pca.df%>%
         mutate(sample = rownames(.)),
       aes(x=PC1, y=PC2))+
  geom_point(shape=21, fill="black", color="white", size=3)+
  ggnetwork::geom_nodetext_repel(aes(label=sample))+
  theme_classic()+
  theme(strip.text = element_text(size=10))+
  labs(x=paste("PC1 (", round(oars.mpx.qc.pca.var[1], digits=1), "%)", sep=""),
       y=paste("PC2 (", round(oars.mpx.qc.pca.var[2], digits=1), "%)", sep=""))+
  facet_wrap(~"PCA of QC samples (Protein)")
# appears like #.2 has concerningly higher variance

# >> MPX samples ----------------------------------------------------------

# clean colnames (to be consistent with one another)
colnames(oars.mpx) = gsub(".*__", "", colnames(oars.mpx))

# filter to oars samples
oars.mpx.oars.only = oars.mpx[,colnames(oars.mpx) %in% c("Protein.IDs", "QC",
                                                             gsub("-", "_", metadata.oars.stool$standard.name))]
colnames(oars.mpx.oars.only) = gsub("\\.1", "",colnames(oars.mpx.oars.only))

# load functional annotations
mpx.functions = read.csv("./metaproteomics/MetaLab_iterative/functional_annotation/functions.tsv", sep="\t", header=F)
colnames(mpx.functions) = mpx.functions[1,] # make first rows into colnames
mpx.functions = mpx.functions[-1,] # and then delete them (there's a more efficient way of doing this)
mpx.functions = mpx.functions[,colnames(mpx.functions) %in% c("Group_ID", "Name", "Protein name",
                                              "Description", "Taxonomy Id", "Taxonomy name",
                                              "Preferred name", "Gene_Ontology_id", "Gene_Ontology_name",
                                              "Gene_Ontology_namespace", "EC_id", "EC_de", "EC_an",
                                              "EC_ca", "KEGG_ko", "KEGG_Pathway_Entry", "KEGG_Pathway_Name",
                                              "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC",
                                              "CAZy", "BiGG_Reaction", "PFAMs", "COG accession", "COG category",
                                              "COG name", "NOG accession", "NOG category", "NOG name")] %>% distinct()


colnames(mpx.functions) = gsub(" ", "_", colnames(mpx.functions))

# create function to collapse proteins by annotation
protein.collapser = function(data,
                             protein.information,
                             annotation,
                             add_tax = FALSE,
                             delimiter = ";"){
  # annotation = bquote(annotation) # dplyr cannot handle quotes
  data.melt <- reshape2::melt(data)
  # remove rows with 0 abundance
  data.melt <- subset(data.melt, value != 0)
  # replace sample names
  colnames(data.melt)[2] <- "code"
  
  # select annotation of choice, and lengthen annotation file, such that X;Y;Z becomes X, Y, Z on separate rows
  # optionally append Taxonomy_name
  if(add_tax == T & annotation == "CAZy"){
    # paste Taxonomy_name to each annotation
    protein.information <- protein.information %>%
      dplyr::select(Name, CAZy, Taxonomy_name)%>%
      rowwise() %>%
      mutate(
        CAZy = paste0(str_split(CAZy, ",")[[1]], "_", Taxonomy_name) %>%
          paste(collapse = ",")
      ) %>%
      ungroup()
  }
  if(add_tax == T & annotation == "Preferred_name"){
    # paste Taxonomy_name to each annotation
    protein.information <- protein.information %>%
      dplyr::select(Name, Preferred_name, Taxonomy_name)%>%
      rowwise() %>%
      mutate(
        Preferred_name = paste0(str_split(Preferred_name, ",")[[1]], "_", Taxonomy_name) %>%
          paste(collapse = ",")
      ) %>%
      ungroup()
  }
  protein.information.long = protein.information[,colnames(protein.information) %in% c("Name", annotation)]
  
  # identify the maximum number of splits (optional, for explicit control)
  max_splits <- paste("split", 1:(max(str_count(protein.information.long[,2], delimiter)+1)))
  
  protein.information.long = tidyr::separate(protein.information.long, col={{annotation}}, into=max_splits, sep=delimiter, remove=T)
  
  protein.information.long = reshape2::melt(protein.information.long, id="Name")
  protein.information.long$value = ifelse(protein.information.long$value == "", NA, protein.information.long$value)
  protein.information.long$variable = NULL
  colnames(protein.information.long)[2] <- "annotation"
  protein.information.long = subset(protein.information.long, !is.na(annotation))
  # append to data
  colnames(data.melt)[1] <- "Name"
  data.melt <- merge(data.melt, 
                     protein.information.long, by="Name")
  # EDGECASE: fix Purine"-"/" "nucleoside (sum them together)
  data.melt$annotation <- ifelse(data.melt$annotation == "Purine-nucleoside phosphorylase",
                                             "Purine nucleoside phosphorylase", data.melt$annotation)
  
  colnames(data.melt)[2] = "code"
  # sum proteins by annotation per sample (code)
  data.melt.meta.summed <- data.melt %>%
    group_by(code, annotation) %>%
    mutate(sum.intensity = sum(value)) %>%
    dplyr::select(code, annotation, sum.intensity) %>% distinct() %>% data.frame()
  
  # done; ready for downstream analyses
  data.melt.meta.summed
}



# :: MPX KEGG PATHWAY --------------------------------------------------------------

oars.mpx.kegg.mat = protein.collapser(data=oars.mpx.oars.only, 
                                               protein.information = mpx.functions,
                                               annotation="KEGG_Pathway_Name",
                                               delimiter = ";")
oars.mpx.kegg.mat$code = gsub("\\.1", "", oars.mpx.kegg.mat$code)

# matrix
oars.mpx.kegg.mat = reshape2::acast(oars.mpx.kegg.mat %>% distinct(),
                                                 code ~ annotation, value.var="sum.intensity",
                                                 fun.aggregate = mean)
rownames(oars.mpx.kegg.mat) = gsub("_", "-", rownames(oars.mpx.kegg.mat))
# save
saveRDS(oars.mpx.kegg.mat, "./metaproteomics/2026_01_15_oars_mpx_kegg.Rds")

dim(oars.mpx.kegg.mat) # 67 samples x 180 KEGG Pathways

oars.mpx.kegg.mat = readRDS("./metaproteomics/2026_01_15_oars_mpx_kegg.Rds")


# :: MPX COG --------------------------------------------------------------

oars.mpx.cog.mat = protein.collapser(data=oars.mpx.oars.only, 
                                               protein.information = mpx.functions,
                                               annotation="COG_name")
oars.mpx.cog.mat$code = gsub("\\.1", "", oars.mpx.cog.mat$code)

# matrix
oars.mpx.cog.mat = reshape2::acast(oars.mpx.cog.mat %>% distinct(),
                                                 code ~ annotation, value.var="sum.intensity",
                                                 fun.aggregate = mean)
rownames(oars.mpx.cog.mat) = gsub("_", "-", rownames(oars.mpx.cog.mat))
# save
saveRDS(oars.mpx.cog.mat, "./metaproteomics/2026_01_15_oars_mpx_cog.Rds")

dim(oars.mpx.cog.mat) # 67 samples x 2653 COGs

oars.mpx.cog.mat = readRDS("./metaproteomics/2026_01_15_oars_mpx_cog.Rds")

# :: MPX CAZy -------------------------------------------------------------

oars.mpx.cazy.mat = protein.collapser(data = oars.mpx.oars.only, 
                                  protein.information = mpx.functions,
                                  annotation = "CAZy",
                                  delimiter = ",")
oars.mpx.cazy.mat$code = gsub("\\.1", "", oars.mpx.cazy.mat$code)
# remove CAZy titled "-"
oars.mpx.cazy.mat = subset(oars.mpx.cazy.mat, annotation != "-")

# matrix
oars.mpx.cazy.mat = reshape2::acast(oars.mpx.cazy.mat %>% distinct() %>%
                                      mutate(code = as.character(code)),
                                                 code ~ annotation, value.var="sum.intensity",
                                                 fun.aggregate = mean)
rownames(oars.mpx.cazy.mat) = gsub("_", "-", rownames(oars.mpx.cazy.mat))


saveRDS(oars.mpx.cazy.mat, "./metaproteomics/2026_01_15_oars_mpx_cazy.Rds")

dim(oars.mpx.cazy.mat) # 67 samples x 66 CAZy

oars.mpx.cazy.mat = readRDS("./metaproteomics/2026_01_15_oars_mpx_cazy.Rds")




# :: MPX CAZy Starch:Mucin ----------------------------------------------------

oars.mpx.cazy.mat = readRDS("./metaproteomics/2026_01_15_oars_mpx_cazy.Rds")

# https://pmc.ncbi.nlm.nih.gov/articles/PMC9120202/ for GH's involved in mucin degradation

mucin.cazy = c("GH101", "GH20", "GH29", "GH33", "GH84", "GH95")
starch.cazy = c("GH13", "CBM48", "CBM20", "GH77")

# calculate Starch:Mucin ratio per sample
oars.mpx.starch.mucin = oars.mpx.cazy.mat %>%
  reshape2::melt() %>%
  mutate(value = ifelse(is.na(value), 0, value)) %>%
  mutate(value = ifelse(value == 0, min(value[value!=0])/2, value)) %>%
  mutate(type = ifelse(Var2 %in% starch.cazy, "starch",
                       ifelse(Var2 %in% mucin.cazy, "mucin", "other"))) %>%
  subset(type != "other") %>%
  group_by(Var1, type) %>%
  mutate(sum.cazy = sum(value)) %>%
  dplyr::select(Var1, sum.cazy, type) %>% distinct() %>%
  reshape2::acast(Var1 ~ type, value.var="sum.cazy") %>%
  data.frame() %>%
  mutate(starch.mucin = log2(starch / mucin)) %>% 
  mutate(standard.name = rownames(.)) %>%
  # keep starch, mucin, and starch:mucin ratio
  mutate(starch = (starch)) %>%
  mutate(mucin = (mucin)) %>%
  dplyr::select(standard.name, starch.mucin, starch, mucin)

ggplot(oars.mpx.starch.mucin,
       aes(x=mucin, y=starch))+
  geom_point(shape=21)+
  geom_smooth(method="lm", color="black")+
  ggpubr::stat_cor(method="spearman")+
  scale_y_log10()+
  scale_x_log10()+
  theme_classic()+
  labs(x=c("Σ Mucin CAZy Intensities"),
       y=c("Σ Starch CAZy Intensities"))


# >> MBX ------------------------------------------------------------------


# :: load data ------------------------------------------------------------

# feature table
oars.mbx.data = read.csv("./metabolomics/oars_mbx_data_raw.csv")

# feature map
oars.mbx.feature.data = read.csv("./metabolomics/oars_mbx_feature_data.csv")

# mbx map
oars.mbx.metadata <- read.csv("./2024_11_08_metabolomics_samples_prioritized_PD.csv")

# mbx sample map
oars.mbx.samples = read.csv("~/Documents/PhD/git_oars_archfolder/metabolomics/Naming-STIN-HM.csv")

# main map
metadata.oars.stool = readRDS("~/Documents/PhD/git_oars_archfolder/2025_07_22_oars_mapping_2.Rds")



# :: process --------------------------------------------------------------

# how many samples
ncol(oars.mbx.data) - 6
# n = 70

# add rownames (bucket = features)
rownames(oars.mbx.data) = oars.mbx.data$Bucket
# remove feature ID data
oars.mbx.data[,c(1:6)]
oars.mbx.data = oars.mbx.data[,-c(1:6)]

# fix sample names
colnames(oars.mbx.data) = oars.mbx.samples$standard.name[match(colnames(oars.mbx.data), make.names(oars.mbx.samples[,1]))]

# subset to OARS
oars.mbx.metadata = subset(oars.mbx.metadata, standard.name %in% metadata.oars.stool$standard.name)
length(unique(oars.mbx.metadata$standard.name))
# n=68 samples

# annotations
oars.mbx.feature.data %>% nrow()
subset(oars.mbx.feature.data, Annotation.choisie != "")$Annotation.choisie %>% unique() %>% length()
# 3668 total features
subset(oars.mbx.feature.data, Annotation.choisie != "")$Annotation.choisie %>% unique() %>% sort()
# 1360 unique annotated features


# Save 2 versions:
# 1: annotated (DO NOT sum up features annotated the same)
# 2: raw (annotate if possible, otherwise use spectra)

oars.mbx.raw = oars.mbx.data %>% as.matrix() %>% reshape2::melt()
colnames(oars.mbx.raw) = c("Bucket", "standard.name", "value")
oars.mbx.raw = merge(oars.mbx.raw,
                     oars.mbx.feature.data[,c("Bucket", "Annotation.choisie")], by="Bucket")

oars.mbx.raw$Annotation = ifelse(as.character(oars.mbx.raw$Annotation.choisie) != "", 
                                 paste(as.character(oars.mbx.raw$Annotation.choisie), as.character(oars.mbx.raw$Bucket), sep=" | "),
                                 as.character(oars.mbx.raw$Bucket))
# good

unique(oars.mbx.raw$Annotation) %>% length()
# 3668 unique 

# now, slice off unannotated
oars.mbx.annotated = subset(oars.mbx.raw, grepl(" \\| ", Annotation))
unique(oars.mbx.annotated$Annotation) %>% length()
# 2012 features

# make matrix
oars.mbx.raw.mat = reshape2::acast(oars.mbx.raw,
                                   standard.name ~ Annotation, value="value") %>% as.data.frame()
dim(oars.mbx.raw.mat) # 70 x 3668

oars.mbx.annotated.mat = reshape2::acast(oars.mbx.annotated,
                                         standard.name ~ Annotation, value="value") %>% as.data.frame()
dim(oars.mbx.annotated.mat) # 70 x 2012

saveRDS(oars.mbx.raw.mat, "./metabolomics/2026_01_15_oars_mbx_raw.Rds")
saveRDS(oars.mbx.annotated.mat, "./metabolomics/2026_01_15_oars_mbx_annotated.Rds")

oars.mbx.raw.mat = readRDS("./metabolomics/2026_01_15_oars_mbx_raw.Rds")
oars.mbx.annotated.mat = readRDS("./metabolomics/2026_01_15_oars_mbx_annotated.Rds")


# :: filter ---------------------------------------------------------------

# note: filtration is applied in analysis script
# the reason for this is to add flexibility for removing non-compliant samples first (or not)


# >> RapidAIM pH ----------------------------------------------------------

# :: load pH data ---------------------------------------------------------

# select OARS microbiomes
oars.rapidaim.list <- read.csv("./oars_performance_characteristics_data/2024_11_22_oars_rapidaim_list.csv")

# process
oars.rapidaim.scores <- readRDS("./2025_06_09_oars_scores.Rds") %>%
  # calculate LFC to PBS
  group_by(HM, timing) %>%
  mutate(but.lfc = log(med.but / med.but[RS_Name == "PBS"], base=2)) %>%
  # indicate if + of -
  mutate(but.dir = ifelse(but.lfc > 0, "response", "nonresponse")) %>%
  # also calculate delta.pH
  mutate(delta.ph = med.ph - med.ph[RS_Name == "PBS"]) %>%
  # subset to just the selected RS
  subset(RS_Name == ifelse(HM == "HM0618" & timing %in% c("0M", "3M"), "ActistarRT", selected))
# good


# >> RapidAIM 16S ---------------------------------------------------------

# :: Process 16S data -----------------------------------------------------

# goal: process data to enable in vitro vs in vivo comparisons (i.e. same taxonomy needed)

# adjust mapping file
oars.mapping.rapidaim$barcode = gsub("\\.", "_", gsub("\\.gg13v5", "", oars.mapping.rapidaim$SampleID))
# I'm not sure where this came from!

# load and merge rds objects
oars.rapidaim.ps = list.files(path = "./oars_performance_characteristics_data", pattern = "pooled", full.names = TRUE)
oars.rapidaim.ps = lapply(oars.rapidaim.ps, readRDS)
oars.rapidaim.ps = lapply(oars.rapidaim.ps, function(x){phyloseq::otu_table(x, taxa_are_rows=F)})

# merge all phyloseqs into 1
oars.rapidaim.ps = do.call(phyloseq::merge_phyloseq, oars.rapidaim.ps)
ncol(oars.rapidaim.ps) # 71665 seq

# remove chimeras
t1 <- Sys.time()
oars.rapidaim.ps.nochimera = dada2::removeBimeraDenovo(oars.rapidaim.ps, method="pooled")
t2 <- Sys.time()
t2-t1 # 1.5 h elapsed for chimera removal
ncol(oars.rapidaim.ps.nochimera) # 8602 ASVs

# visualize seq sums
rownames(oars.rapidaim.ps.nochimera) = gsub("_filtered.fastq.gz", "", rownames(oars.rapidaim.ps.nochimera))
oars.rapidaim.ps.nochimera.df = data.frame(barcode = rownames(oars.rapidaim.ps.nochimera),
                                           sums = rowSums(oars.rapidaim.ps.nochimera))
oars.rapidaim.ps.nochimera.df = merge(oars.rapidaim.ps.nochimera.df,
                                      oars.mapping.rapidaim, by="barcode") %>%
  mutate(sample = paste(HM, RS_Name, Replicate, sep="_"))
# plot
ggplot(oars.rapidaim.ps.nochimera.df,
       aes(x=sums, y=reorder(HM, sums)))+
  geom_point(shape=21, aes(fill=RS_Name))+
  geom_vline(xintercept=50000)+
  scale_x_log10()+
  theme_minimal()

# rarefy
t1 = Sys.time() 
set.seed(25)
oars.rapidaim.ps.nochimera.rare = phyloseq::rarefy_even_depth(oars.rapidaim.ps.nochimera, sample.size=50000, replace=F)
t2 = Sys.time()
t2 - t1 # 5 min elapsed for rarefaction
# lost ~ half the samples

# assign taxonomy
t1 = Sys.time() # 9:45pm
oars.rapidaim.taxa = dada2::assignTaxonomy(phyloseq::taxa_names(oars.rapidaim.ps.nochimera.rare), "~/Documents/PhD/16s_databases/gg2_2024_09_toSpecies_trainset.fa.gz", multithread=TRUE)
t2 = Sys.time()
t2 - t1 # 20 min elapsed for taxonomy assignment

# add taxa and lca groups
oars.rapidaim.taxa <- data.frame(oars.rapidaim.taxa) %>%
  mutate(LCA = 
           ifelse(!is.na(Species)&Species!="s__", paste(as.character(Genus), as.character(Species), sep="_"), 
                  ifelse(!is.na(Genus)&Genus!="g__", paste(Genus), 
                         ifelse(!is.na(Family)&Family!="f__", paste(Family),
                                ifelse(!is.na(Order)&Order!="o__", paste(Order),
                                       ifelse(!is.na(Class)&Class!="c__", paste(Class),
                                              ifelse(!is.na(Phylum)&Phylum!="p__", paste(Phylum),
                                                     ifelse(is.na(Phylum), "undefined", paste(Phylum)))))))))
# append unique number to duplicated ASVs
# Call this "Taxa"
oars.rapidaim.taxa$Taxa <- make.unique(as.character(oars.rapidaim.taxa$LCA), sep = "_")

# make phyloseq
rownames(oars.mapping.rapidaim) = oars.mapping.rapidaim$barcode
rownames(oars.rapidaim.ps.nochimera.rare) = gsub("_filtered.fastq.gz", "", rownames(oars.rapidaim.ps.nochimera.rare))

oars.rapidaim.ps.final = phyloseq::merge_phyloseq(phyloseq::otu_table(oars.rapidaim.ps.nochimera.rare, taxa_are_rows=F),
                                                  phyloseq::tax_table(as.matrix(oars.rapidaim.taxa)),
                                                  phyloseq::sample_data(oars.mapping.rapidaim))

saveRDS(oars.rapidaim.ps.final, "./2025_07_24_oars.rapidaim.ps.final.Rds")
oars.rapidaim.ps.final = readRDS("./2025_07_24_oars.rapidaim.ps.final.Rds")




# >> final data process -----------------------------------------------------------

metadata.oars.stool %>% nrow()
metadata.oars.stool.asv %>% nrow()

# add HM0932-STL-08 to .asv file (missing ASV, but has MGX, clinical, FFQ, etc)
metadata.oars.stool.asv = rbind(subset(metadata.oars.stool.asv, standard.name != "HM0932-STL-08"),
                                data.frame(metadata.oars.stool[rownames(metadata.oars.stool)=="HM0932-STL-08",] %>%
                                mutate(richness = NA, shannon = NA)))%>%data.frame() %>%
  group_by(HM) %>% arrange(standard.name)
metadata.oars.stool.asv %>% nrow()

# add DNA variables
metadata.oars.stool.asv = merge(metadata.oars.stool.asv,
                                oars.phyloseq.butyrogens.i, by="standard.name", all=T)
metadata.oars.stool.asv = merge(metadata.oars.stool.asv,
                                oars.phyloseq.butyrogens.ii, by="standard.name", all=T)
metadata.oars.stool.asv = merge(metadata.oars.stool.asv,
                                oars.asv.functionalredundancy, by="standard.name", all=T)
metadata.oars.stool.asv = merge(metadata.oars.stool.asv,
                                oars.asv.mlp.preds, by="standard.name", all=T)
metadata.oars.stool.asv = merge(metadata.oars.stool.asv,
                                oars.mgx.mlp.preds, by="standard.name", all=T)
metadata.oars.stool.asv = merge(metadata.oars.stool.asv,
                                oars.mpx.starch.mucin, by="standard.name", all=T)

metadata.oars.stool.asv %>% nrow()

# check dataset completion
metadata.oars.stool.asv$standard.name %>% unique() %>% length()
# HM0924-STL-12 and HM0932-STL-08  are missing values
# fcal
subset(metadata.oars.stool.asv, is.na(fcal)) # HM0924-STL-12 missing fcal
# water
subset(metadata.oars.stool.asv, is.na(stool_water_perc)) # HM0924-STL-12  missing stool water
# load
subset(metadata.oars.stool.asv, is.na(load.asv)) # HM0924-STL-12 missing stool water
# Note: HM0924-STL-12 can be removed from many analyses because it was a repeat collection at timepoint 3M
# so, let's remove it from ASV data
metadata.oars.stool = subset(metadata.oars.stool, standard.name != "HM0924-STL-12")
metadata.oars.stool.asv = subset(metadata.oars.stool.asv, standard.name != "HM0924-STL-12")
# and remove from mgx
oars.mgx.taxa = oars.mgx.taxa[rownames(oars.mgx.taxa)!= "HM0924-STL-12",]
# down to 67 samples, tops

metadata.oars.stool$timing = factor(metadata.oars.stool$timing, levels=c("0M", "3M", "6M", "9M", "12M"))
metadata.oars.stool.asv$timing = factor(metadata.oars.stool.asv$timing, levels=c("0M", "3M", "6M", "9M", "12M"))

metadata.oars.stool$rs.col = ifelse(!is.na(metadata.oars.stool$RS_Name), metadata.oars.stool$RS_Name, "grey")
metadata.oars.stool$rs.col = factor(metadata.oars.stool$rs.col, levels=c("grey", rs.names[rs.names %in% unique(metadata.oars.stool$RS_Name)]))
metadata.oars.stool.asv$rs.col = ifelse(!is.na(metadata.oars.stool.asv$RS_Name), metadata.oars.stool.asv$RS_Name, "grey")
metadata.oars.stool.asv$rs.col = factor(metadata.oars.stool.asv$rs.col, levels=c("grey", rs.names[rs.names %in% unique(metadata.oars.stool.asv$RS_Name)]))

# fix up scores names
oars.rapidaim.scores = merge(oars.rapidaim.scores,
                             metadata.oars.stool[,c("HM", "timing", "standard.name")], by=c("HM", "timing"))

# add MPX plates
metadata.oars.stool = merge(metadata.oars.stool,
                                oars.mpx.map.plates[,c("standard.name", "plate")],
                                by="standard.name", all.x=T, all.y=F) %>%
  mutate(plate = as.factor(plate))

metadata.oars.stool.asv = merge(metadata.oars.stool.asv,
                                oars.mpx.map.plates[,c("standard.name", "plate")],
                                by="standard.name", all.x=T, all.y=F) %>%
  mutate(plate = as.factor(plate))

table(metadata.oars.stool.asv$standard.name) %>% range()

# >> save -----------------------------------------------------------

# save independently
saveRDS(metadata.oars.stool, "~/Documents/PhD/git_oars_archfolder/2026_01_15_oars_mapping_2.Rds")
saveRDS(metadata.oars.stool.asv, "~/Documents/PhD/git_oars_archfolder/2026_01_15_oars_mapping_asv_2.Rds")


# so, we have:
# general metadata
metadata.oars.stool
# 16S metadata, with alpha diversity
metadata.oars.stool.asv
# ASV table (median)
oars.asv.data.median
# Glom table (median)
oars.asv.data.glom

# save these, load into analysis script
save(
  # non-compliant
  non.compliant,
  # 16S data
  metadata.oars.stool, 
  metadata.oars.stool.asv, 
  oars.asv.data.median,
  oars.asv.data.glom,
  oars.asv.functionalredundancy,
  oars.phyloseq.butyrogens.i,
  oars.phyloseq.butyrogens.ii,
  oars.asv.mlp.preds,
  # MGX data
  oars.mgx.taxa,
  oars.mgx.mlp.preds,
  # MPX data
  oars.mpx.kegg.mat,
  oars.mpx.cog.mat,
  oars.mpx.cazy.mat,
  # MBX data
  oars.mbx.raw.mat,
  oars.mbx.annotated.mat,
  # RapidAIM pH
  oars.rapidaim.scores,
  # RapidAIM 16S
  #oars.rapidaim.ps.final,
  # Supplemental data
  kircher.butyrogens.plot, # butyrogen contributors
  # destination
  file = "./2026_01_13_oars_16s_data_meta.Renv")


# :: ----------------------------------------------------------------------


# :: ----------------------------------------------------------------------


# >> old ------------------------------------------------------------------


# :: MPX CAZy w Taxa -------------------------------------------------------------

oars.mpx.cazy.tax = protein.collapser(data = oars.mpx.oars.only, 
                                      protein.information = mpx.functions,
                                      add_tax = T,
                                      annotation = "CAZy",
                                      delimiter = ",")
# collapse replicates
oars.mpx.cazy.tax$code = gsub("\\.1", "", oars.mpx.cazy.tax$code)


# matrix
oars.mpx.cazy.tax.mat = reshape2::acast(oars.mpx.cazy.tax %>% data.frame() %>% distinct(),
                                        code ~ annotation, value.var="sum.intensity",
                                        fun.aggregate = mean)
# delete unannotated
oars.mpx.cazy.tax.mat = oars.mpx.cazy.tax.mat[,!grepl("-", colnames(oars.mpx.cazy.tax.mat))]
# replace NA with 0
oars.mpx.cazy.tax.mat[is.na(oars.mpx.cazy.tax.mat)] = 0

saveRDS(oars.mpx.cazy.tax.mat, "./metaproteomics/2025_06_28_oars_mpx_cazy_tax.Rds")

# :: MPX Butyrate -------------------------------------------------------------

# use Preferred_name
# select only the following: thl, bhbd, cro, bcd, but, buk

oars.mpx.but.mat = protein.collapser(data = oars.mpx.oars.only, 
                                     protein.information = mpx.functions,
                                     add_tax = F,
                                     annotation = "Preferred_name",
                                     delimiter = ",")
# collapse replicates
oars.mpx.but.mat$code = gsub("\\.1", "", oars.mpx.but.mat$code)

# matrix
oars.mpx.but.mat = reshape2::acast(oars.mpx.but %>% data.frame() %>% distinct(),
                                   code ~ annotation, value.var="sum.intensity",
                                   fun.aggregate = mean)
# delete unannotated
oars.mpx.but.mat = oars.mpx.but.mat[,!grepl("-", colnames(oars.mpx.but.mat))]
# replace NA with 0
oars.mpx.but.mat[is.na(oars.mpx.but.mat)] = 0

saveRDS(oars.mpx.but.mat, "./metaproteomics/2025_06_28_oars_mpx_but.Rds")

# :: MPX Butyrate w Taxa-------------------------------------------------------------

# use Preferred_name
# select only the following: thl, bhbd, cro, bcd, but, buk

oars.mpx.but.mat = protein.collapser(data = oars.mpx.oars.only, 
                                     protein.information = mpx.functions,
                                     add_tax = T,
                                     annotation = "Preferred_name",
                                     delimiter = ",")
# collapse replicates
oars.mpx.but.mat$code = gsub("\\.1", "", oars.mpx.but.mat$code)


# matrix
oars.mpx.but.tax.mat = reshape2::acast(oars.mpx.but.tax %>% data.frame() %>% distinct(),
                                       code ~ annotation, value.var="sum.intensity",
                                       fun.aggregate = mean)
# delete unannotated
oars.mpx.but.tax.mat = oars.mpx.but.tax.mat[,!grepl("-", colnames(oars.mpx.but.tax.mat))]
# replace NA with 0
oars.mpx.but.tax.mat[is.na(oars.mpx.but.tax.mat)] = 0

saveRDS(oars.mpx.but.tax.mat, "./metaproteomics/2025_06_28_oars_mpx_but_tax.Rds")

mpx.functions = NULL
oars.mpx = NULL




# :: 16S Phylotree --------------------------------------------------------


# build a phylogenetic tree using ASV alignments

amplicon.data.gg.rare <- readRDS("./2025_09_10_rs_trial_16s_data_rarefied_gg2.Rds")

# need to select representative ASV for each taxa

# original phyloseq object
tree_tax_sums = data.frame(sum_abun = phyloseq::taxa_sums(amplicon.data.gg.rare))
tree_tax_sums$OTU = rownames(tree_tax_sums)
colnames(tree_tax_sums)
# tax_table
asv.oars.rare.tax.df = phyloseq::tax_table(amplicon.data.gg.rare) %>% data.frame()
asv.oars.rare.tax.df$OTU = rownames(asv.oars.rare.tax.df)

# merge
tree_tax_sums = merge(tree_tax_sums,
                      asv.oars.rare.tax.df[,c("OTU","Phylum","Class","Order", "Family","Genus","Species", "Taxa", "LCA")], by="OTU") %>% distinct()
# add ASV sequence
tree_tax_sums$ASV = data.frame(phyloseq::refseq(amplicon.data.gg.rare))[,1]

# subset to sig taxa
# tree_tax_sums = subset(tree_tax_sums, LCA %in% subset(mend.asv.glom.lm.results.all, sig == "*")$feature)

# take largest sum of each LCA
tree_tax_sums = tree_tax_sums %>%
  group_by(LCA) %>%
  slice_max(sum_abun, n = 1, with_ties=FALSE)

hist(tree_tax_sums$sum_abun)
tree_tax_sums %>% arrange(-sum_abun)

## build tree
library("DECIPHER"); packageVersion("DECIPHER")
library("phangorn"); packageVersion("phangorn")

amplicon.data.gg.rare.seqs = Biostrings::DNAStringSet(tree_tax_sums$ASV)
names(amplicon.data.gg.rare.seqs) = tree_tax_sums$LCA

alignment = AlignSeqs(amplicon.data.gg.rare.seqs, anchor=NA, processors=4)

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm)

fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2) # faster
t1 <- Sys.time()
#fitGTR <- optim.pml(fit, model = "GTR", optInv = TRUE, optGamma = TRUE, k = 8)
t2 <- Sys.time() # TOO SLOW and no progress bar 

tree = phyloseq::phy_tree(fitGTR$tree)

# save
saveRDS(tree, "./2025_06_29_oars_tree.Rds")
saveRDS(tree_tax_sums, "./2025_06_29_oars_tree_meta.Rds")


# :: MPX Reaction Network -------------------------------------------------

# Goal: establish a full reaction network (link all possible reactants and products inferrable from protein annotations 'EC_an')

mpx.functions[,13] %>% head(n=30)
colnames(mpx.functions)[13]

mpx.reaction.network = mpx.functions %>% dplyr::select(EC_an) %>% distinct() # %>%unique() %>% length()
# 1799 unique entries

# split multiple (e.g. EC_1;EC_2 into EC_1 on one row, and EC_2 on another)
max_splits_reactions <- paste("split", 1:(max(str_count(mpx.reaction.network[,1], ";")+1)))
max_splits_reactions = tidyr::separate(mpx.reaction.network, col=EC_an, into=max_splits_reactions, sep=";", remove=T)
# make long
max_splits_reactions = reshape2::melt(as.matrix(max_splits_reactions))
max_splits_reactions[,1] = NULL
max_splits_reactions = distinct(max_splits_reactions)
nrow(max_splits_reactions)
# 1933 unique reactions
colnames(max_splits_reactions)

# make dataframe
max_splits_reactions = data.frame(max_splits_reactions) %>%
  tidyr::separate(col=value, into=c("reactants", "products"), sep= " = ", remove=T)
max_splits_reactions
max_splits_reactions = subset(max_splits_reactions, reactants != "")
max_splits_reactions = na.omit(max_splits_reactions)

max_splits_reactions %>% head(n = 20)




