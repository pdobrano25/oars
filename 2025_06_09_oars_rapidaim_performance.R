### 2025_06_09  OARS RS Selections

# Goal: Re-perform RS selections using OTUs for OARS patients; plus add 6M, 12M where available
# Assess stability compared to permuted selections
# Assess outcome of pH-based selection


library("ggplot2"); library("dplyr"); library("tidyverse")
setwd("~/Documents/PhD/git_oars_archfolder")

# :: load -----------------------------------------------------------------

rs.names <- c("Authentic", "BobsRedMill", "MSPrebiotic", "LetsDoOrganic", "HiMaize260", "Novelose330", "ActistarRT", "FibersymRW", "Versafibe1490")

rs.names.pbs <- c("PBS", "Authentic", "BobsRedMill", "MSPrebiotic", "LetsDoOrganic", "HiMaize260", "Novelose330", "ActistarRT", "FibersymRW", "Versafibe1490")

# give 9 RS consistent colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
labelcolors <- data.frame(cols = c(gg_color_hue(9), "#000000"),
                          label = c(seq(1:9), "P"))
cols = c(gg_color_hue(5), "#000000")
labelcolors = data.frame(cols = c(
  cols[1], # Authentic
  cols[1], # BobsRedMill
  cols[1], # MSPrebiotic
  cols[2], # LetsDoOrganic
  cols[3], # HiMaize260
  cols[3], # Novelose330
  cols[4], # ActistarRT
  cols[4], # FibersymRW
  cols[5], # Versafibe1490
  "#000000"# PBS
))

# :: load data --------------------------------------------------------------

# load OARS patient list
oars.patient.list <- read.csv("./oars_git_mapping/2023_12_04_oars_patient_list.csv")

# load OARS RapidAIM (by Peter) list
oars.rapidaim.list <- read.csv("./oars_performance_characteristics_data/2024_11_22_oars_rapidaim_list.csv")

# subset to samples that have been sequenced
oars.rapidaim.list = subset(oars.rapidaim.list, grepl("AS16S", seq_run))

## Process Peter's mapping
# list mapping files
oars_files_to_read <- list.files(
  path = "./oars_performance_characteristics_data", # directories to search within (but ignore MBX folder, which has pseudoduplicates of n=15 Mapping files )
  pattern = "*HM.*.csv", # regex pattern, some explanation below
  recursive = TRUE,          # search subdirectories
  full.names = TRUE          # return the full path
)
# remove first HM0924.03
oars_files_to_read = oars_files_to_read[oars_files_to_read!="./oars_performance_characteristics_data/2022_11_29_HM0924.03_Map.csv"]
# and from oars.rapidaim.list
oars.rapidaim.list = subset(oars.rapidaim.list, name != "HM0924.03")

# load files
oars.mapping.files.rs.list <- lapply(oars_files_to_read, read_csv)  # read all the matching files
length(oars.mapping.files.rs.list) # check number of files processed: 34
# rename first column to "barcode" and create new column for RapidAIM run
oars.mapping.files.rs.list <- lapply(1:length(oars.mapping.files.rs.list), function(x){
  data.subset = oars.mapping.files.rs.list[[x]]
  colnames(data.subset)[1] = "Barcode"
  data.subset$source = oars_files_to_read[x]
  data.subset
})
# in order to merge into a single dataframe, must take intersect of colnames
mincol <- Reduce(intersect, lapply(oars.mapping.files.rs.list, colnames))
oars.mapping.files.rs.list <- lapply(oars.mapping.files.rs.list, function(x) x[mincol])
# now rbind
oars.mapping.df.rs = do.call(rbind, oars.mapping.files.rs.list) %>% distinct()
# clean up dataframe for downstream processing
colnames(oars.mapping.df.rs) = make.names(colnames(oars.mapping.df.rs))
# fix Post_PCR column
oars.mapping.df.rs$Post.PCR.DNA.Yield..ng.uL. <- as.numeric(oars.mapping.df.rs$Post.PCR.DNA.Yield..ng.uL.)
oars.mapping.df.rs <- oars.mapping.df.rs[,c("Barcode", "Treatment","Replicate", "HM", "RS_Name",
                                            "RS_Type", "RS_Source", "RS_Polymorph", "pH", "Qubit",
                                            "Post.PCR.DNA.Yield..ng.uL.", "source")]
oars.mapping.df.rs <- oars.mapping.df.rs %>% data.frame()
# good

## Process Mais' mapping
# load Mais' mapping files
oars.mais.mapping <- read.csv("./oars_git_mapping/2025_01_09_mais_oars_mapping_250107.csv")
# clean Post_PCR name
oars.mais.mapping$Post.PCR.DNA.Yield..ng.uL. <- as.numeric(oars.mais.mapping$amplicon.conc)
# make Barcode from SampleID
oars.mais.mapping$Barcode = oars.mais.mapping$SampleID
# remove HM from HM
oars.mais.mapping$HM <- gsub("HM", "", oars.mais.mapping$sample)
# add pseudo-source
oars.mais.mapping$source = paste("Mais", oars.mais.mapping$sample, sep="_")

## Unify both mapping files
# merge with larger df
all_columns <- union(names(oars.mapping.df.rs), names(oars.mais.mapping)) # Get all unique column names
# Function to add missing columns (from ChatGPT)
add_missing_columns <- function(df, all_cols) {
  missing_cols <- setdiff(all_cols, names(df))
  df[missing_cols] <- NA
  df <- df[, all_cols] # Reorder columns
  return(df)
}
# Standardize both data frames
oars.mapping.df.rs <- add_missing_columns(oars.mapping.df.rs, all_columns)
oars.mais.mapping <- add_missing_columns(oars.mais.mapping, all_columns)
# rbind the standardized data frames
oars.mapping.df.rs.mais <- rbind(oars.mapping.df.rs, oars.mais.mapping)
# remove NA rows
oars.mapping.df.rs.mais = subset(oars.mapping.df.rs.mais, !is.na(Barcode))
tail(oars.mapping.df.rs.mais)
# 2291 unique barcodes

# which RapidAIMs were performed by Mais?
oars.mais.mapping[,"sample"] %>% unique()
mais.samples = c("HM0844.05", "HM0924.04", "HM0940.04", "HM0639.05","HM0902.03", "HM0903.06","HM0883.04","HM0874.04","HM0899.04","HM0906.04","HM0932.04","HM0844.06")

# subset(oars.rapidaim.scores.pc, timing == "6M" & sample %in% mais.samples)[,c("sample")] %>%unique() %>% sort()
# "HM0639.05" "HM0844.05" "HM0902.03" "HM0903.06" "HM0924.04" "HM0940.04"


# >> run through RS Selection Algorithm -----------------------------------

# Note: copied from original algorithm implementation

controls <- readRDS("./oars_performance_characteristics_data/controls_nonibd-dc-cloud_n83.rds")
tree <- phyloseq::read_tree( "~/Documents/PhD/Resistant_Starch/97_otus_gg13v5.tree")

rerun=F
if(rerun==T){
t0 <- Sys.time()

oars.rapidaim.data <- lapply(1:nrow(oars.rapidaim.list), function(run) {
  
  # Troubleshoot 
  # run = 7
  
  # select microbiome
  sample = oars.rapidaim.list[run,]$name
  HM = oars.rapidaim.list[run,]$HM
  timing = oars.rapidaim.list[run,]$timing
  seq_run = oars.rapidaim.list[run,]$seq_run
  
  print(paste(sample, " ", run))
  

  # :: load data -------------------------------------------------------------
  
  # mapping file; subset to microbiome of interest
  oars.mapping.df.rs.hm <- oars.mapping.df.rs.mais[grepl(sample, oars.mapping.df.rs.mais$source),]
  # add timing
  oars.mapping.df.rs.hm$timing = timing
  # remove non-existant rows
  oars.mapping.df.rs.hm = oars.mapping.df.rs.hm[!is.na(oars.mapping.df.rs.hm$Barcode),]
  # replace Post.PCR value == 0 with a half/min
  oars.mapping.df.rs.hm$Post.PCR.DNA.Yield..ng.uL. <- as.numeric(lapply(oars.mapping.df.rs.hm$Post.PCR.DNA.Yield..ng.uL., function(x) replace(x, x == 0, min(oars.mapping.df.rs.hm$Post.PCR.DNA.Yield..ng.uL.[oars.mapping.df.rs.hm$Post.PCR.DNA.Yield..ng.uL.>0], na.rm = TRUE)/2))) # replace 0's with half the smallest value
  
  
  # load 16S data; pre-processed QIIME 97% OTU + GreenGenes13.5
  oars.rapidaim.biom <- phyloseq::import_biom(paste("./oars_performance_characteristics_data/", seq_run, "_gg13v5_no_doubletons.biom", sep=""), parseFunction = phyloseq::parse_taxonomy_greengenes) #import biom file
  row.names(oars.mapping.df.rs.hm) <- oars.mapping.df.rs.hm$Barcode
  oars.mapping.df.rs.hm <- phyloseq::sample_data(oars.mapping.df.rs.hm, errorIfNULL=T) #create sample data phyloseq object from metadata file
  oars.mapping.df.rs.hm$Sample <- rownames(oars.mapping.df.rs.hm) #add a new column with sample names (for downstream analyses)
  RapidAIM.ps <- phyloseq::merge_phyloseq(oars.rapidaim.biom, oars.mapping.df.rs.hm, tree) #Create phyloseq object (with tree)
  RapidAIM.ps #check phyloseq object  
  
  # HM0819.03 = 5908 taxa
  
  # :: decontam --------------------------------------------------------------
  
  # remove taxa that aren't actually present
  RapidAIM.ps <- phyloseq::filter_taxa(RapidAIM.ps, function(x) mean(x) != 0, TRUE)
  
  # HM0819.03 = 4596 taxa
  
  
  # decontam using frequency method (note: elsewhere, I tend to use concentration method)
  contamdf.freq <- decontam::isContaminant(RapidAIM.ps, method="frequency", conc="Post.PCR.DNA.Yield..ng.uL.")
  RapidAIM.ps.noncontam <- phyloseq::prune_taxa(!contamdf.freq$contaminant, RapidAIM.ps) # 162
  # RapidAIM.ps.noncontam <- RapidAIM.ps # this would override the decontam!
  RapidAIM.counts <- data.frame(Depth=phyloseq::sample_sums(RapidAIM.ps.noncontam),
                                Sample=phyloseq::sample_data(RapidAIM.ps.noncontam)$Sample,
                                RS_Name=phyloseq::sample_data(RapidAIM.ps.noncontam)$RS_Name,
                                Replicate=phyloseq::sample_data(RapidAIM.ps.noncontam)$Replicate,
                                stringsAsFactors = FALSE)
  RapidAIM.seq.depth.plot <- ggplot(RapidAIM.counts, aes(x=reorder(RS_Name, Depth), Depth)) + 
    geom_boxplot() +
    geom_point(aes(color=Replicate))+
    geom_hline(yintercept=120000, linetype="solid", color="red",size=0.5)+
    coord_flip()+
    theme_bw()+
    labs(x="", y="Depth (Reads)", title="Sequencing Depth")
  RapidAIM.seq.depth.plot
  
  # HM0819.03 = 4434 taxa
  
  # :: read depth filter -----------------------------------------------------
  
  if(sample == "HM0618.01"){
    depth = 75000
  }
  if(sample == "HM0940.02"){
    depth = 100000
  }
  if(!sample %in% c("HM0618.01", "HM0940.02")){
    depth = 120000
  }
  RapidAIM.filtered <- phyloseq::prune_samples(phyloseq::sample_sums(RapidAIM.ps.noncontam)>=depth, RapidAIM.ps.noncontam) # Include only samples with "sufficient" depth (120k reads)
  
  # :: prevalence filter -----------------------------------------------------
  
  # Execute filter: Keep OTUs with a minimum read >=2 in at least 4% of the samples (ie 2, to keep the Zymo controls intact)
  RapidAIM.filtered.2 <- phyloseq::filter_taxa(RapidAIM.filtered,function(x) sum(x>=2)>=length(x)*.04,prune = TRUE)
  RapidAIM.filtered.3 <- phyloseq::subset_taxa(RapidAIM.filtered.2, Family!="mitochondria") # remove reads classified as mitochondria
  RapidAIM.filtered.4 <- phyloseq::subset_taxa(RapidAIM.filtered.3, Phylum!="Chloroplast") # remove reads classified as mitochondria

  # EDGE CASES:
  if(sample == "HM0874.01"){
  RapidAIM.filtered.4 = phyloseq::subset_samples(RapidAIM.filtered.4,
                                                 Barcode != "AS16S.101.gg13v5.BC015")
  }
  if(sample == "HM0883.01"){
  outlier <- phyloseq::sample_data(RapidAIM.filtered.4)[which(sample_data(RapidAIM.filtered.4)$RS_Name=="MED" & sample_data(RapidAIM.filtered.4)$Replicate=="C"),]$Sample
  RapidAIM.filtered.4 <- subset_samples(RapidAIM.filtered.4, Sample!=outlier)
  }
  if(sample == "HM0906.02"){
  RapidAIM.filtered.4 <- subset_samples(RapidAIM.filtered.4, Barcode!="AS16S.170.gg13v5.BC001")
  }

  # :: export phyloseq ------------------------------------------------------
  
  RapidAIM.filtered.4 # for exporting
  
  # HM0819.03 = 2403 taxa
  
  
  # :: merge with controls ------------------------------
  controls <- readRDS("./oars_performance_characteristics_data/controls_nonibd-dc-cloud_n83.rds")

  RapidAIM.controls.ps <- phyloseq::merge_phyloseq(RapidAIM.filtered.4, controls)  # Combine  RapidAIM data with reference cohort
  # NOTE: includes only taxa found in RapidAIM.filtered.4, not controls (probably overlap)
  RapidAIM.controls.ps <- phyloseq::transform_sample_counts(RapidAIM.controls.ps, function(x) 100 * x/sum(x)) # TSS scale library sizes to 100%
  
  # :: distance --------------------------------------------------------------
  
  dist_methods <- unlist(phyloseq::distanceMethodList)
  RapidAIM.controls.distancematrix <- phyloseq::distance(RapidAIM.controls.ps, method="bray", type="samples") # Calculate Bray-Curtis dissimilarity matrix
  id.control <- as.character(phyloseq::sample_data(phyloseq::subset_samples(RapidAIM.controls.ps, Location=="Distal Colon" & Diagnosis=="Control"))$Sample) # Vector of ONLY reference samples (comparisons will be made to this group)
  id.all <- as.character(phyloseq::sample_data(phyloseq::subset_samples(RapidAIM.controls.ps))$Sample) # Vector of ALL samples
  RapidAIM.controls.distancematrix.matrix <- as.matrix(RapidAIM.controls.distancematrix) # Matrix needed for median-distance code, but not distance-to-centroid
  
  # Calculate Median Distance to Reference Cohort for all samples
  median.distance.df <- do.call(rbind.data.frame,
                                lapply(id.all, function (x) median(RapidAIM.controls.distancematrix.matrix[x, id.control])))
  names(median.distance.df)[1] <- "Distance" 
  
  RapidAIM.distance.df <- cbind(phyloseq::sample_data(RapidAIM.controls.ps), distance=median.distance.df)
  RapidAIM.distance.df <- RapidAIM.distance.df[grepl(HM, RapidAIM.distance.df$source),]
  RapidAIM.distance.df = RapidAIM.distance.df[, colnames(RapidAIM.distance.df) %in% c("RS_Name", "Replicate", "pH", "Qubit", "source", "timing", "Distance")]
  # Remove scores that will not be included in RS selection analysis
  RapidAIM.distance.df <- subset(RapidAIM.distance.df, !RS_Name %in%
                                   c("Stool", "WaterB", "RawStool", "MED", "ZymoC", "ZymoD","Amioca",
                                     "FOS", "ExtractC", "Spike-in", "FOS-100", "PBS-100", "CLS"))
  RapidAIM.distance.df$Sample = rownames(RapidAIM.distance.df)
  
  ggplot(RapidAIM.distance.df %>% mutate(RS_Name = factor(RS_Name, levels=rs.names.pbs)),
         aes(x=RS_Name, y=Distance))+
    geom_boxplot()+
    geom_point()+
    theme_minimal()

  # :: butyrogens ------------------------------------------------------------

  
  # Convert reads to relative abundances
  RapidAIM.tss  <- phyloseq::transform_sample_counts(RapidAIM.filtered.4, function(x) x / sum(x) )
  # Remove unimportant samples
  RapidAIM.tss <- phyloseq::subset_samples(RapidAIM.tss, !RS_Name %in%
                                             c("ZymoC", "ZymoD", "NA", "ExtractC", "MED", "Amioca",
                                               "Spike-in", "WaterB", "FOS-100", "PBS-100", "Stool", "RawStool"))
  
  # Create dataframes for butyrogens derived from Mottawea et al., 2016)
  RapidAIM.butyrogens <- phyloseq::subset_taxa(RapidAIM.tss, Family=="Lachnospiraceae" | Genus=="Blautia" | Genus=="Roseburia" | Genus=="Eubacterium" | Genus=="Ruminococcus" | Genus=="Clostridium" | Genus=="Faecalibacterium") # Note, Lachnospiraceae already includes several genera listed; listed again for clarity
  
  RapidAIM.butyrogens.ps.placeholder <- microbiome::aggregate_taxa(RapidAIM.butyrogens, 'Kingdom', verbose=FALSE) #collapse all butyrogens to one group (kingdom)
  RapidAIM.butyrogens.counts <- as.data.frame(colSums((phyloseq::otu_table(RapidAIM.butyrogens.ps.placeholder))))
  RapidAIM.butyrogens.counts <- RapidAIM.butyrogens.counts %>% tibble::rownames_to_column(var="Sample")
  names(RapidAIM.butyrogens.counts)[2] <- "But_Abundance"

  # :: finalize --------------------------------------------------------------
  
  RapidAIM.final = merge(RapidAIM.distance.df, RapidAIM.butyrogens.counts, by = "Sample") %>% data.frame()
  RapidAIM.final$HM = HM
  RapidAIM.final$sample = sample
  
  
  
  ggplot(RapidAIM.final %>% mutate(RS_Name = factor(RS_Name, levels=rs.names.pbs)),
         aes(x=RS_Name, y=But_Abundance))+
    geom_boxplot()+
    geom_point()+
    theme_minimal()
  
  # calculate median per replicate
  
  RapidAIM.final.median = RapidAIM.final %>%
    group_by(RS_Name) %>%
    mutate(med.dis = median(Distance)) %>%
    mutate(med.but = median(But_Abundance)) %>% 
    mutate(med.ph = median(pH)) %>% # note: median pH may be missing a replicate if barcode did not sequence
    dplyr::select(RS_Name, med.ph, med.dis, med.but) %>% distinct() %>% data.frame()
  
  RapidAIM.final.median$HM = HM
  RapidAIM.final.median$timing = timing
  RapidAIM.final.median$sample = sample
  
  # select RS
  # first, remove PBS and re-add
  RapidAIM.final.median.pbs = subset(RapidAIM.final.median, RS_Name == "PBS")
  RapidAIM.final.median.pbs$Z_score = NA
  
  RapidAIM.final.median = subset(RapidAIM.final.median, RS_Name != "PBS")
  
  # Convert  values to Z-scores to identify top RS
  RapidAIM.final.median$Z_score <- scale(RapidAIM.final.median$med.but, center = TRUE, scale = TRUE)
  candidate.RS <- subset(RapidAIM.final.median, Z_score > 1)
  if(nrow(candidate.RS) >0 ){ # select lowest distance among top butyrogen RS
    top.RS <- subset(candidate.RS, med.dis == min(candidate.RS$med.dis)) # Determine optimal RS by which one has the lowest Distance value
  } else { # or else select top butyrogen RS
    top.RS = subset(RapidAIM.final.median, med.but == max(med.but))
  }
  # re-add PBS
  RapidAIM.final.median = rbind(RapidAIM.final.median, RapidAIM.final.median.pbs) %>% data.frame()
  RapidAIM.final.median$selected = top.RS$RS_Name
  
  
  # :: export all ----------------------------------------------------------------
  
  list(
    # phyloseq object
    RapidAIM.filtered.4,
    # median scores
    RapidAIM.final.median,
    # scores
    RapidAIM.final
  )
  
}) # END LOOP
t1 <- Sys.time()

t1 - t0 # 20 min

saveRDS(oars.rapidaim.data, "2025_08_08_oars.rapidaim.data.Rds")
}

oars.rapidaim.data = readRDS("2025_08_08_oars.rapidaim.data.Rds")

oars.rapidaim.data

# :: analysis -------------------------------------------------------------


## 1. merge phyloseq object
oars.rapidaim.phyloseq = do.call(phyloseq::merge_phyloseq, lapply(1:nrow(oars.rapidaim.list), function(x){
  #x = 1
  oars.rapidaim.data[[x]][[1]]@phy_tree <- NULL
  # need tree for identifying relationship of important OTUs
  oars.rapidaim.data[[x]][[1]]
}))
oars.rapidaim.phyloseq # 1891 samples

# add tree
oars.rapidaim.phyloseq = phyloseq::merge_phyloseq(
  oars.rapidaim.phyloseq,
  phyloseq::read_tree("../16s_databases/97_otus_gg13v5.tree"))


saveRDS(oars.rapidaim.phyloseq, "./2025_06_09_oars_phyloseq.Rds")



## 2. median scores
oars.rapidaim.scores.pc = do.call(rbind, lapply(1:nrow(oars.rapidaim.list), function(x){
  oars.rapidaim.data[[x]][[2]]
}))
saveRDS(oars.rapidaim.scores.pc, "./2025_06_09_oars_scores.Rds")
# note: some median pH scores are based on conditions with missing replicates (because barcodes did not sequence)

# Replace pH values with full data
oars.rapidaim.ph = oars.mapping.df.rs.mais %>%
  subset(RS_Name %in% rs.names.pbs) %>%
  group_by(HM) %>%
  # add replicate if missing
  mutate(Replicate = ifelse(is.na(Replicate), c("A", "B", "C"), Replicate))%>%
  group_by(HM, RS_Name)%>%
  mutate(med.ph = median(pH)) %>%
  # calculate delta
  group_by(HM) %>%
  mutate(delta.ph = med.ph - med.ph[RS_Name == "PBS"]) %>%
  dplyr::select(HM, RS_Name, Replicate, pH, med.ph, delta.ph)
colnames(oars.rapidaim.ph)[1] = "sample"
oars.rapidaim.ph$sample = paste("HM", oars.rapidaim.ph$sample, sep="")
# merge with patient data
oars.rapidaim.ph = merge(oars.rapidaim.ph,
                         oars.rapidaim.scores.pc[,c("HM","RS_Name", "sample", "timing", "selected", "Z_score", "med.but", "med.dis")] %>% distinct(), by=c("sample", "RS_Name"))
colnames(oars.rapidaim.ph)[colnames(oars.rapidaim.ph)=="selected"] = "rs.selected"
# fix edge cases
oars.rapidaim.ph = oars.rapidaim.ph %>%
  mutate(rs.selected = ifelse(sample %in% c("HM0618.00", "HM0618.01"), "ActistarRT", rs.selected))

saveRDS(oars.rapidaim.ph, "./2025_06_09_oars_scores_ph.Rds")

write.csv(oars.rapidaim.ph, "./2025_06_09_oars_scores_ph.csv")



# ::RESET -----------------------------------------------------------------


###
oars.rapidaim.phyloseq = readRDS("./2025_06_09_oars_phyloseq.Rds")

oars.rapidaim.scores.pc = readRDS("./2025_06_09_oars_scores.Rds")

# add HMno
oars.rapidaim.scores.pc$HMno = substr(oars.rapidaim.scores.pc$HM, 1, 6)
# factor RS
oars.rapidaim.scores.pc$RS_Name <- factor(oars.rapidaim.scores.pc$RS_Name, levels=rs.names)
# factor timing (not needed if only 0-6 are plotted)
#oars.rapidaim.scores.pc$timing <- factor(oars.rapidaim.scores.pc$timing, levels=c("0M", "3M", "6M", "12M"))

# :: rs selections --------------------------------------------------

oars.deid.list <- setNames(oars.deid$code, oars.deid$HM)

oars.rapidaim.scores.pc.selection.plot <- 
  subset(oars.rapidaim.scores.pc, RS_Name != "PBS" & timing %in% c("0M", "3M", "6M"))%>%
  #mutate(Timing = as.character(timing)) %>%
  ggplot(aes(x=timing, y=Z_score))+
  # add vertical bars
  geom_vline(xintercept=c(1:3), color="black")+
  # geom_vline(xintercept=4, color="black")+
  # add points
  geom_point(aes(fill=RS_Name), shape=21, size=2, alpha=1)+
  # overlay selected RS (for clarity)
  geom_point(data = subset(oars.rapidaim.scores.pc, RS_Name == selected & timing %in% c("0M", "3M", "6M")), 
             aes(fill=RS_Name), size=2.5, shape=21, alpha=1)+
  # add lines
  geom_path(aes(group = RS_Name, color=RS_Name), alpha=0.6, linewidth=0.3, linetype=1)+
  # label selected RS
  geom_label(data = subset(oars.rapidaim.scores.pc, RS_Name == selected & timing %in% c("0M", "3M", "6M"))%>%
               mutate(Timing = timing), 
             aes(label=RS_Name, color=RS_Name, x=timing, y=3, vjust=0.55),
             size=2)+
  scale_y_continuous(limits=c(-2.2,3.2))+
  geom_hline(yintercept=1, linetype=2, color="red", alpha=0.5)+
  scale_color_manual(values = c(labelcolors$cols[c(1:9)]))+
  scale_fill_manual(values = c(labelcolors$cols[c(1:9)]))+
  theme_minimal()+theme(legend.position="none")+
  labs(x="", y="Butyrogen Z-Score")+
  facet_wrap(~HMno, nrow=3,
             # deid
             labeller = labeller(HMno = oars.deid.list))+
  theme_minimal()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))
oars.rapidaim.scores.pc.selection.plot

# break down frequency of RS's being selected, per timepoint (Stacked barplot)
oars.rapidaim.scores.pc.frequencies = oars.rapidaim.scores.pc %>%
  subset(timing %in% c("0M", "3M", "6M")) %>%
  subset(RS_Name != "PBS") %>%
  dplyr::select(HM, timing, selected) %>% distinct() %>% data.frame() %>%
  dplyr::select(timing, selected) %>%
  table() %>% data.frame() %>%
  group_by(timing) %>%
  mutate(perc = Freq / sum(Freq)) %>% data.frame() 

oars.rapidaim.scores.pc.frequencies$selected = factor(oars.rapidaim.scores.pc.frequencies$selected, levels=rs.names)

oars.rapidaim.scores.pc.frequencies.plot = ggplot(subset(oars.rapidaim.scores.pc.frequencies, timing != "12M") %>%
                                                    mutate(Timing = timing),
       aes(x=selected, y=Freq))+
  geom_bar(stat="identity", fill="white", alpha=1)+
  geom_bar(stat="identity", position="stack", 
           aes(fill=selected), alpha=1,
           color="black")+
  scale_fill_manual(values= labelcolors$cols[c(1:9)])+
  scale_y_continuous(breaks=seq(from=0, to=13, by=1))+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=45, hjust=1),
        plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.x=element_blank(),

        strip.text = element_text(size=10),
        strip.background = element_rect(
          color="black"),
        legend.position="none")+
  labs(x="", y="Frequency Selected", fill="")+
  facet_wrap(~Timing, labeller = label_both, nrow=3)
oars.rapidaim.scores.pc.frequencies.plot

oars.rapidaim.scores.pc.selection.plot+ 
  oars.rapidaim.scores.pc.frequencies.plot+
  patchwork::plot_layout(widths=c(4,1))

# :: ph selections --------------------------------------------------

oars.rapidaim.scores.pc = oars.rapidaim.scores.pc %>%
  group_by(sample) %>%
  mutate(ph.selected = RS_Name[med.ph == min(med.ph)])

oars.rapidaim.scores.pc.ph.selection.plot <- ggplot(subset(oars.rapidaim.scores.pc, RS_Name != "PBS" & timing != "12M"),
                                              #     aes(x=timing, y=log(med.but,base=2)))+
                                              aes(x=timing, y=med.ph))+
  
  # add vertical bars
  geom_vline(xintercept=1, color="black")+
  geom_vline(xintercept=2, color="black")+
  geom_vline(xintercept=3, color="black")+
  #geom_vline(xintercept=4, color="black")+

  # add lines
  geom_path(aes(group = RS_Name, color=RS_Name), alpha=0.6)+
  # label selected RS
  geom_label(data = subset(oars.rapidaim.scores.pc, RS_Name == ph.selected & timing != "12M"), 
             aes(label=RS_Name, color=RS_Name, x=timing, y=8.5, vjust=0.95),
             size=2.5)+
  # overlay selected RS (for clarity)
  geom_point(aes(fill=RS_Name), shape=21, size=1.5)+
  # geom_hline(yintercept=1, color="red", alpha=0.5)+
  scale_color_manual(values = c(labelcolors$cols[c(1:9)]))+
  scale_fill_manual(values = c(labelcolors$cols[c(1:9)]))+
  theme_minimal()+theme(legend.position="none",
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="", y="Median pH")+
  facet_wrap(~HM, nrow=4)

oars.rapidaim.scores.pc.ph.selection.plot

# break down frequency of RS's being selected, per timepoint (Stacked barplot)
oars.rapidaim.scores.pc.ph.frequencies = oars.rapidaim.scores.pc %>%
  subset(RS_Name != "PBS") %>%
  dplyr::select(HM, timing, ph.selected) %>% distinct() %>% data.frame() %>%
  dplyr::select(timing, ph.selected) %>%
  table() %>% data.frame() %>%
  group_by(timing) %>%
  mutate(perc = Freq / sum(Freq)) %>% data.frame() 

oars.rapidaim.scores.pc.ph.frequencies$ph.selected = factor(oars.rapidaim.scores.pc.ph.frequencies$ph.selected, levels=rs.names)

oars.rapidaim.scores.pc.ph.frequencies.plot = ggplot(
  subset(oars.rapidaim.scores.pc.ph.frequencies, timing != "12M") %>%
    mutate(Timing = timing),
       aes(x=ph.selected, y=Freq))+
  geom_bar(stat="identity", fill="white", alpha=1)+
  geom_bar(stat="identity", position="stack", 
           aes(fill=ph.selected), alpha=1,
           color="black")+
  scale_fill_manual(values= labelcolors$cols[c(1:9)])+
  scale_y_continuous(breaks=seq(from=0, to=13, by=1))+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=45, hjust=1),
        plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y=element_blank(),
        strip.text = element_text(size=12),
        strip.background = element_rect(
          color="black"),
        legend.position="none")+
  labs(x="", y="Frequency Selected", fill="")+
  facet_wrap(~Timing, labeller = label_both, nrow=3)
oars.rapidaim.scores.pc.ph.frequencies.plot

oars.rapidaim.scores.pc.ph.selection.plot+ 
  oars.rapidaim.scores.pc.ph.frequencies.plot+
  patchwork::plot_layout(widths=c(4,1))


# >> vignette (UPDATED HM's) -------------------------------------------------------------

# select 2 microbiomes to illustrate selection algorithm
# HM0899 0M (Scenario A) and HM0639 3M (Scenario B) look like solid examples

# :: Z-scores -------------------------------------------------------------

oars.rapidaim.scores.pc.899 = readRDS("./2025_06_09_oars_scores.Rds")
oars.rapidaim.scores.pc.899 = subset(oars.rapidaim.scores.pc.899, HM == "HM0899")
oars.rapidaim.scores.pc.899 = subset(oars.rapidaim.scores.pc.899, timing %in% c("0M"))
oars.rapidaim.scores.pc.899 = subset(oars.rapidaim.scores.pc.899, RS_Name != "PBS")
oars.rapidaim.scores.pc.899$RS_Name = factor(oars.rapidaim.scores.pc.899$RS_Name, levels=rs.names)

# create dataframe for labels
oars.rapidaim.scores.pc.899.1 = subset(oars.rapidaim.scores.pc.899, RS_Name != "PBS" & timing == "0M")
oars.rapidaim.scores.pc.899.1 = oars.rapidaim.scores.pc.899.1 %>% arrange(-med.but)
oars.rapidaim.scores.pc.899.1$y.pos = seq(from = 0.07, 
                                       to = min(oars.rapidaim.scores.pc.899.1$med.but),length=9)
oars.rapidaim.scores.pc.899.1$y.pos = ifelse(oars.rapidaim.scores.pc.899.1$RS_Name == "LetsDoOrganic", 
                                             subset(oars.rapidaim.scores.pc.899.1, RS_Name == "LetsDoOrganic")$med.but,
                                             oars.rapidaim.scores.pc.899.1$y.pos)
# shift sec_axis calculated values
scale_factor.899 = max(scale(range(oars.rapidaim.scores.pc.899.1$med.but)))

oars.rapidaim.scores.pc.899.1.selection <- ggplot(oars.rapidaim.scores.pc.899.1,
                                              #     aes(x=timing, y=log(med.but,base=2)))+
                                              aes(x=1, y=med.but))+
  # add vertical bars
  geom_vline(xintercept=1, color="black")+
 
  # add label
  geom_text(data=oars.rapidaim.scores.pc.899.1,
            x=1.01, aes(y=y.pos, label=RS_Name, color=RS_Name), hjust=0)+
  geom_segment(data=oars.rapidaim.scores.pc.899.1,
               x=1.0025, xend=1, aes(y=y.pos, yend=med.but, 
                            color=RS_Name))+
  geom_segment(data=oars.rapidaim.scores.pc.899.1,
               x=1.0025, xend=1.0095, aes(y=y.pos, yend=y.pos,
                                    color=RS_Name))+
  geom_label(data=subset(oars.rapidaim.scores.pc.899.1, RS_Name == selected),
            x=1.01, aes(y=y.pos, label=RS_Name, color=RS_Name), hjust=0)+
  # add points
  geom_point(aes(fill=RS_Name), size=3, shape=21)+
  # add others (including Z-score = 1)
  geom_hline(yintercept=(mean(oars.rapidaim.scores.pc.899.1$med.but) + 1 * sd(oars.rapidaim.scores.pc.899.1$med.but)), 
             color="red", linetype=2, alpha=0.5)+
  scale_fill_manual(values = c(labelcolors$cols[c(1:9)]))+
  scale_color_manual(values = c(labelcolors$cols[c(1:9)]))+
  scale_x_continuous(limits=c(0.995, 1.03))+
  scale_y_continuous(
    # Features of the first axis
    name = " ",
    labels = scales::label_percent(),
    # Add a second axis and specify its features
    sec.axis = sec_axis(transform=~scale(.)+scale_factor.899, 
                        name=" "))+
  labs(x="")+
  facet_wrap(~"Butyrogens % (Z-score)")+
  theme_minimal()+theme(legend.position="none",
                        panel.grid.minor.x  = element_blank(),
                        panel.grid.major.x  = element_blank(),
                        axis.text.x=element_blank(),
                        axis.ticks.x=element_blank(),
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="white"))
oars.rapidaim.scores.pc.899.1.selection


# create dataframe for labels
oars.rapidaim.scores.pc.883 = readRDS("./2025_06_09_oars_scores.Rds")
oars.rapidaim.scores.pc.883 = subset(oars.rapidaim.scores.pc.883, HM == "HM0883")
oars.rapidaim.scores.pc.883.2 = subset(oars.rapidaim.scores.pc.883, RS_Name != "PBS" & timing == "3M")
oars.rapidaim.scores.pc.883.2 = oars.rapidaim.scores.pc.883.2 %>% arrange(-med.but)
oars.rapidaim.scores.pc.883.2$RS_Name = factor(oars.rapidaim.scores.pc.883.2$RS_Name, levels=rs.names)
# make non-top labels below Z-score 1

oars.rapidaim.scores.pc.883.2$y.pos =c(subset(oars.rapidaim.scores.pc.883.2, RS_Name %in% c("BobsRedMill", "MSPrebiotic"))$med.but,
                                       seq(from = max(subset(oars.rapidaim.scores.pc.883.2, !RS_Name %in% c("BobsRedMill", "MSPrebiotic"))$med.but), 
                                           to = min(oars.rapidaim.scores.pc.883.2$med.but), length=7)) 
                                             

# shift sec_axis calculated values
scale_factor.883 = ((mean(oars.rapidaim.scores.pc.883.2$med.but)))


oars.rapidaim.scores.pc.883.2.selection <- ggplot(oars.rapidaim.scores.pc.883.2,
                                                  #     aes(x=timing, y=log(med.but,base=2)))+
                                                  aes(x=1, y=med.but))+
  # add vertical bars
  geom_vline(xintercept=1, color="black")+
  # add others (including Z-score = 1)
  geom_hline(yintercept=(mean(oars.rapidaim.scores.pc.883.2$med.but) + 1 * sd(oars.rapidaim.scores.pc.883.2$med.but)), 
             color="red", linetype=2, alpha=0.5)+
  # add label
  geom_text(data=oars.rapidaim.scores.pc.883.2,
            x=1.01, aes(y=y.pos, label=RS_Name, color=RS_Name), hjust=0)+
  geom_segment(data=oars.rapidaim.scores.pc.883.2,
               x=1.0025, xend=1, aes(y=y.pos, yend=med.but, 
                                     color=RS_Name))+
  geom_segment(data=oars.rapidaim.scores.pc.883.2,
               x=1.0025, xend=1.0095, aes(y=y.pos, yend=y.pos,
                                          color=RS_Name))+
  geom_label(data=subset(oars.rapidaim.scores.pc.883.2, RS_Name == selected),
             x=1.01, aes(y=y.pos, label=RS_Name, color=RS_Name), hjust=0)+
  # add points
  geom_point(aes(fill=RS_Name), size=3, shape=21)+

  scale_fill_manual(values = c(labelcolors$cols[c(1:9)]))+
  scale_color_manual(values = c(labelcolors$cols[c(1:9)]))+
  scale_x_continuous(limits=c(0.995, 1.03))+
  scale_y_continuous(
    # Features of the first axis
    name = " ",
    labels = scales::label_percent(),
    # Add a second axis and specify its features
    sec.axis = sec_axis(transform= ~scale(., center=T)*0.77,
                        name=" "))+
  labs(x="")+
  facet_wrap(~"Butyrogens % (Z-score)")+
  theme_minimal()+theme(legend.position="none",
                        panel.grid.minor.x  = element_blank(),
                        panel.grid.major.x  = element_blank(),
                        axis.text.x=element_blank(),
                        axis.ticks.x=element_blank(),
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="white"))
oars.rapidaim.scores.pc.883.2.selection

# and without selecting

oars.rapidaim.scores.pc.883.2.tie <- ggplot(oars.rapidaim.scores.pc.883.2,
                                                  #     aes(x=timing, y=log(med.but,base=2)))+
                                                  aes(x=1, y=med.but))+
  # add vertical bars
  geom_vline(xintercept=1, color="black")+
  # add others (including Z-score = 1)
  geom_hline(yintercept=(mean(oars.rapidaim.scores.pc.883.2$med.but) + 1 * sd(oars.rapidaim.scores.pc.883.2$med.but)), 
             color="red", linetype=2, alpha=0.5)+
  # add label
  geom_text(data=oars.rapidaim.scores.pc.883.2,
            x=1.01, aes(y=y.pos, label=RS_Name, color=RS_Name), hjust=0)+
  geom_segment(data=oars.rapidaim.scores.pc.883.2,
               x=1.0025, xend=1, aes(y=y.pos, yend=med.but, 
                                     color=RS_Name))+
  geom_segment(data=oars.rapidaim.scores.pc.883.2,
               x=1.0025, xend=1.0095, aes(y=y.pos, yend=y.pos,
                                          color=RS_Name))+
  #geom_label(data=subset(oars.rapidaim.scores.pc.883.2, RS_Name == selected),
  #           x=1.01, aes(y=y.pos*100, label=RS_Name, color=RS_Name), hjust=0)+
  # add points
  geom_point(aes(fill=RS_Name), size=3, shape=21)+
  scale_fill_manual(values = c(labelcolors$cols[c(1:9)]))+
  scale_color_manual(values = c(labelcolors$cols[c(1:9)]))+
  scale_x_continuous(limits=c(0.995, 1.03))+
  scale_y_continuous(
    # Features of the first axis
    name = " ",
    labels = scales::label_percent(),
    # Add a second axis and specify its features
    sec.axis = sec_axis(transform= ~scale(., center=T)*0.77,
                        name=" "))+
  labs(x="")+
  facet_wrap(~"Butyrogens % (Z-score)")+
  theme_minimal()+theme(legend.position="none",
                        panel.grid.minor.x  = element_blank(),
                        panel.grid.major.x  = element_blank(),
                        axis.text.x=element_blank(),
                        axis.ticks.x=element_blank(),
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="white"))
oars.rapidaim.scores.pc.883.2.tie

# :: PCoA -------------------------------------------------------------

# and a PCoA for HM0883 3M
controls <- readRDS("./oars_performance_characteristics_data/controls_nonibd-dc-cloud_n83.rds")

oars.rapidaim.scores.pc.883 = readRDS("./2025_06_09_oars_phyloseq.Rds")
oars.rapidaim.scores.pc.883 = phyloseq::subset_samples(oars.rapidaim.scores.pc.883, HM == "0883.02")
oars.rapidaim.scores.pc.883 = phyloseq::merge_phyloseq(oars.rapidaim.scores.pc.883,
                                                    controls)
oars.rapidaim.scores.pc.883 <- phyloseq::transform_sample_counts(oars.rapidaim.scores.pc.883, function(x) 100 * x/sum(x)) # TSS scale library sizes to 100%
oars.rapidaim.scores.pc.883.bray = phyloseq::ordinate(oars.rapidaim.scores.pc.883, "PCoA", "bray")
oars.rapidaim.scores.pc.883.bray = oars.rapidaim.scores.pc.883.bray$vectors[,c(1,2)] %>% data.frame() %>%
  mutate(Sample = rownames( oars.rapidaim.scores.pc.883.bray$vectors))
oars.rapidaim.scores.pc.883.bray = merge(oars.rapidaim.scores.pc.883.bray,
                                      phyloseq::sample_data(oars.rapidaim.scores.pc.883) %>% data.frame(),
                                      by="Sample")
oars.rapidaim.scores.pc.883.bray$RS_Name = factor(oars.rapidaim.scores.pc.883.bray$RS_Name, levels=rs.names.pbs)
oars.rapidaim.scores.pc.883.bray$type = ifelse(!is.na(oars.rapidaim.scores.pc.883.bray$source), "RapidAIM", "Controls")
# take median
oars.rapidaim.scores.pc.883.bray = oars.rapidaim.scores.pc.883.bray %>%
  group_by(type, RS_Name) %>%
  mutate(med.1 = ifelse(type == "RapidAIM", median(Axis.1), Axis.1)) %>%
  mutate(med.2 = ifelse(type == "RapidAIM", median(Axis.2), Axis.2)) %>%
  dplyr::select(type, RS_Name, med.1, med.2) %>% distinct()
# add medoid of controls
oars.rapidaim.scores.pc.883.bray.medoid = oars.rapidaim.scores.pc.883.bray %>%
  subset(type == "Controls") %>%
  mutate(medoid.1 = median(med.1)) %>%
  mutate(medoid.2 = median(med.2)) %>% dplyr::select(medoid.1, medoid.2) %>% distinct()

oars.rapidaim.scores.pc.883.1.pcoa = ggplot()+
  geom_point(data = subset(oars.rapidaim.scores.pc.883.bray, type == "RapidAIM" & RS_Name %in% rs.names),
             aes(x=med.1, y=med.2, fill = RS_Name), shape=21, size=2.5)+
    # draw distance
  geom_segment(aes(
    x=subset(oars.rapidaim.scores.pc.883.bray, RS_Name == "BobsRedMill")$med.1,
    xend=oars.rapidaim.scores.pc.883.bray.medoid$medoid.1,
    y=subset(oars.rapidaim.scores.pc.883.bray, RS_Name == "BobsRedMill")$med.2,
    yend=oars.rapidaim.scores.pc.883.bray.medoid$medoid.2),
    color = labelcolors$cols[1],
    linetype=2, alpha=0.9)+
  geom_segment(aes(
    x=subset(oars.rapidaim.scores.pc.883.bray, RS_Name == "MSPrebiotic")$med.1,
    xend=oars.rapidaim.scores.pc.883.bray.medoid$medoid.1,
    y=subset(oars.rapidaim.scores.pc.883.bray, RS_Name == "MSPrebiotic")$med.2,
    yend=oars.rapidaim.scores.pc.883.bray.medoid$medoid.2),
    color = labelcolors$cols[1],
    linetype=2, alpha=0.9)+
  # draw controls
  stat_density_2d(data = subset(oars.rapidaim.scores.pc.883.bray, type == "Controls"),
                  geom = "polygon", aes(x=med.1, y=med.2,alpha = ..level..), fill = "grey",
                  adjust=1.1)+
  geom_point(data = subset(oars.rapidaim.scores.pc.883.bray, type == "Controls"),
             aes(x=med.1, y=med.2), shape=21, fill = "black", color="white", size=2)+
  scale_fill_manual(values=labelcolors$cols[c(1:9)])+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=12),
                        strip.background = element_rect(color="white"))+
  facet_wrap(~"PCoA + non-IBD microbiomes (Bray-Curtis)")+
  labs(x="Axis 1", y="\nAxis 2")
oars.rapidaim.scores.pc.883.1.pcoa
# remember: it's MEDIAN DISTANCE to controls; not distance to medoid
oars.rapidaim.scores.pc.883.2.selection
oars.rapidaim.scores.pc.899.1.selection

# plot bar/stacked plot for distances
oars.rapidaim.scores.pc.883.2 = oars.rapidaim.scores.pc.883.2 %>% arrange(med.dis)
oars.rapidaim.scores.pc.883.2$y.pos.dis = seq(from = min(oars.rapidaim.scores.pc.883.2$med.dis), 
                                       to = max(oars.rapidaim.scores.pc.883.2$med.dis),length=9)

oars.rapidaim.scores.pc.883.2.distance.plot <- ggplot(subset(oars.rapidaim.scores.pc.883.2, RS_Name != "PBS"),
                                               #     aes(x=timing, y=log(med.but,base=2)))+
                                               aes(x=1, y=med.dis))+
  # add vertical bars
  geom_vline(xintercept=1, color="black")+
  
  # add label
  geom_text(data=oars.rapidaim.scores.pc.883.2,
            x=1.01, aes(y=y.pos.dis, label=RS_Name, color=RS_Name), hjust=0,
            size=3)+
  geom_segment(data=oars.rapidaim.scores.pc.883.2,
               x=1.0025, xend=1, aes(y=y.pos.dis, yend=med.dis, 
                                     color=RS_Name))+
  geom_segment(data=oars.rapidaim.scores.pc.883.2,
               x=1.0025, xend=1.0095, aes(y=y.pos.dis, yend=y.pos.dis,
                                          color=RS_Name))+
  geom_label(data=subset(oars.rapidaim.scores.pc.883.2, RS_Name == selected),
             x=1.01, aes(y=y.pos.dis, label=RS_Name, color=RS_Name), hjust=0,
             size=3)+
  # add points
  geom_point(aes(fill=RS_Name), size=2, shape=21)+
  # add others
  scale_fill_manual(values = c(labelcolors$cols[c(1:9)]))+
  scale_color_manual(values = c(labelcolors$cols[c(1:9)]))+
  scale_x_continuous(limits=c(0.995, 1.03))+
  labs(x="", y=" ")+
  facet_wrap(~"Dissimilarity to non-IBD")+
  theme_minimal()+theme(legend.position="none",
                        panel.grid.minor.x  = element_blank(),
                        panel.grid.major.x  = element_blank(),
                        axis.text.x=element_blank(),
                        axis.text.y=element_text(size=8),
                        axis.ticks.x=element_blank(),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(color="white"))
oars.rapidaim.scores.pc.883.2.distance.plot


(oars.rapidaim.scores.pc.899.1.selection+
  (oars.rapidaim.scores.pc.883.1.pcoa + 
  inset_element(oars.rapidaim.scores.pc.883.2.distance.plot,
                left = 0.003, bottom = 0.005, 
                right = 0.52, top = 0.55))+
  oars.rapidaim.scores.pc.883.2.selection+
    patchwork::plot_layout(nrow=1, widths=c(1,1.6,1))) %>%
  ggsave(filename="./oars_plots/oars_rs_selection.pdf",
         width=12, height=4.5, device = cairo_pdf)

# :: pH -------------------------------------------------------------------

oars.rapidaim.scores.pc.819.1 = oars.rapidaim.scores.pc.819.1 %>% arrange(med.ph)
oars.rapidaim.scores.pc.819.1$y.ph = seq(from = min(oars.rapidaim.scores.pc.819.1$med.ph), 
                                          to = max(oars.rapidaim.scores.pc.819.1$med.ph),length=9) 

oars.rapidaim.scores.pc.819.1.selection.ph <- ggplot(subset(oars.rapidaim.scores.pc.819.1, RS_Name != "PBS") %>%
                                                       mutate(ph.select = ifelse(med.ph == min(med.ph), RS_Name, NA)),
                                                  #     aes(x=timing, y=log(med.but,base=2)))+
                                                  aes(x=1, y=med.ph))+
  # add vertical bars
  geom_vline(xintercept=1, color="black")+
  # add label
  geom_text(data=oars.rapidaim.scores.pc.819.1,
            x=1.01, aes(y=y.ph, label=RS_Name, color=RS_Name), hjust=0)+
  geom_segment(data=oars.rapidaim.scores.pc.819.1,
               x=1.0025, xend=1, aes(y=y.ph, yend=med.ph, 
                                     color=RS_Name))+
  geom_segment(data=oars.rapidaim.scores.pc.819.1,
               x=1.0025, xend=1.0095, aes(y=y.ph, yend=y.ph,
                                          color=RS_Name))+
  # add points
  geom_point(aes(fill=RS_Name), size=2.5, shape=21)+
  # add others
  scale_fill_manual(values = c(labelcolors$cols[c(1:9)]))+
  scale_color_manual(values = c(labelcolors$cols[c(1:9)]))+
  scale_x_continuous(limits=c(0.995, 1.03))+
  scale_y_continuous(
    # Features of the first axis
    name = "pH"
  )+
  labs(x="")+
  theme_minimal()+theme(legend.position="none",
                        panel.grid.minor.x  = element_blank(),
                        panel.grid.major.x  = element_blank(),
                        axis.text.x=element_blank(),
                        axis.ticks.x=element_blank(),
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))
oars.rapidaim.scores.pc.819.1.selection.ph

oars.rapidaim.scores.pc.819.2.selection+
oars.rapidaim.scores.pc.819.2.selection.ph

# save
oars.rapidaim.scores.pc.819.1.selection.ph %>%
  ggsave(filename="./oars_plots/oars_supp_ph_selections_819.pdf",
         width=3, height=4, device = cairo_pdf)

# >> stability ------------------------------------------------------------

# does instability of RS rankings correlate with instability of microbiome composition?

# first, calculate Bray-Curtis of the microbiomes
# focus on 0M and 3M
load(file = "./2025_05_24_oars_16s_data_meta.Renv")
oars.stool.stability.data = oars.asv.data.median
# subset to 0M and 3M stools
oars.stool.stability.data = oars.stool.stability.data[subset(metadata.oars.stool.asv, timing %in% c("0M", "3M", "6M"))$standard.name,]
# calculate Bray-Curtis and PCoA
oars.stool.stability.bray = vegan::vegdist(oars.stool.stability.data, method="bray") 
# perform PCoA
oars.stool.stability.pcoa = ape::pcoa(oars.stool.stability.bray)
# extract data from pcoa
oars.stool.stability.pcoa.df = data.frame(oars.stool.stability.pcoa$vectors[,c(1:2)])
oars.stool.stability.pcoa.df$standard.name = rownames(oars.stool.stability.pcoa.df)
# add metadata
oars.stool.stability.pcoa.df = merge(oars.stool.stability.pcoa.df,
                                     metadata.oars.stool.asv, by="standard.name")
# PCoA is now plottable
oars.stool.stability.pcoa.df
oars.stool.stability.pcoa.plot = ggplot(oars.stool.stability.pcoa.df,
       aes(x=Axis.1, y=Axis.2))+
  geom_line(aes(group=HM), linetype=2, alpha=0.5)+
  geom_point(shape=21, size=2.5, fill="white")+
  facet_wrap(~"Bray-Curtis PCoA")+
  theme_classic()+theme(strip.text=element_text(size=10))+
  labs(x=paste("Axis.1 (", round(oars.stool.stability.pcoa$values$Relative_eig[1],digits=3)*100, "%)", sep=""),
       y=paste("Axis.2 (", round(oars.stool.stability.pcoa$values$Relative_eig[2],digits=3)*100, "%)", sep=""))
oars.stool.stability.pcoa.plot

# now calculate distances
# Subset to samples and comparisons of interest
oars.stool.stability.bray.df = reshape2::melt(as.matrix(oars.stool.stability.bray))
# configure as a "mapping file" so we can use other functions
oars.stool.stability.bray.df$standard.name = oars.stool.stability.bray.df$Var1
# subset to self-comparisons
oars.stool.stability.bray.df = subset(oars.stool.stability.bray.df, substr(Var1, 1, 6) == substr(Var2, 1, 6))
# delete same-sample comparison
oars.stool.stability.bray.df = subset(oars.stool.stability.bray.df, Var1 != Var2)
# add timepoints
oars.stool.stability.bray.df$time1 = metadata.oars.stool.asv[match(oars.stool.stability.bray.df$Var1, metadata.oars.stool.asv$standard.name),]$timing
oars.stool.stability.bray.df$time2 = metadata.oars.stool.asv[match(oars.stool.stability.bray.df$Var2, metadata.oars.stool.asv$standard.name),]$timing
# subset to 0M_3M and 3M_6M
oars.stool.stability.bray.df$time = paste(oars.stool.stability.bray.df$time1, 
                                          oars.stool.stability.bray.df$time2, sep="_")
oars.stool.stability.bray.df = subset(oars.stool.stability.bray.df, time %in% c("0M_3M", "3M_6M"))
oars.stool.stability.bray.df$HM = substr(oars.stool.stability.bray.df$Var1, 1, 6)
colnames(oars.stool.stability.bray.df)[colnames(oars.stool.stability.bray.df)=="value"] = "dissimilarity"
# good

# now ICC of RS rankings
oars.rapidaim.stability = readRDS("./2025_06_09_oars_scores.Rds")
oars.rapidaim.stability = subset(oars.rapidaim.stability, timing %in% c("0M", "3M", "6M"))
oars.rapidaim.stability$microbiome = paste(oars.rapidaim.stability$HM,
                                           oars.rapidaim.stability$timing, sep="_")
oars.rapidaim.stability = subset(oars.rapidaim.stability, RS_Name != "PBS")
# Remove HM0759 because it is missing 6M sample
oars.rapidaim.stability = subset(oars.rapidaim.stability, HM != "HM0759")

oars.rapidaim.stability.icc = do.call(rbind, lapply(unique(oars.rapidaim.stability$HM), function(x){
  data.subset = subset(oars.rapidaim.stability, HM == x)
  # cast
  do.call(rbind, lapply(c("0M_3M", "3M_6M"), function(y){
  if(y == "0M_3M"){ 
  data.subset = merge(subset(data.subset, timing == "0M")[,c("RS_Name", "Z_score")],
                      subset(data.subset, timing == "3M")[,c("RS_Name", "Z_score")], by="RS_Name")
  }
  if(y == "3M_6M"){
    data.subset = merge(subset(data.subset, timing == "3M")[,c("RS_Name", "Z_score")],
                        subset(data.subset, timing == "6M")[,c("RS_Name", "Z_score")], by="RS_Name")
  }
  icc_result = irr::icc(data.subset[,-1], model = "oneway", type = "agreement", unit = "single")
  kendall_result = cor.test(data.subset$Z_score.x,
                            data.subset$Z_score.y, method="kendall")
  # save
  data.frame(HM = x,
             icc = icc_result$value,
             tau = kendall_result$estimate,
             time = y)
  }))
}))
plot(oars.rapidaim.stability.icc[,c(2:3)])

# now correlate Microbiome Stability with RS Ranking Stability
oars.rapidaim.stability = merge(oars.rapidaim.stability.icc,
                                oars.stool.stability.bray.df, by=c("HM", "time"))
# add metadata
oars.rapidaim.stability = merge(oars.rapidaim.stability,
                                metadata.oars.stool.asv[,c("standard.name","diagnosis", "compliant", "richness","shannon","fcal","adj.fiber")], by="standard.name")
# lm
oars.rapidaim.stability.stats = lmerTest::lmer(scale(icc) ~ scale(dissimilarity) + diagnosis + adj.fiber + (1|HM), subset(oars.rapidaim.stability, compliant==T)) %>%
  summary() %>% coef() %>% data.frame()
# p = 0.045, significant after adjusting for diagnosis and fiber.intake

# plot
oars.rapidaim.stability.plot = ggplot(subset(oars.rapidaim.stability,compliant==T),
       aes(x=dissimilarity, y=icc))+
  geom_line(aes(group=HM), alpha=0.5, linetype=2)+
  geom_point(shape=21, size=2.5, fill="white")+
  geom_smooth(method="lm", se=F, color="black")+
  #ggpubr::stat_cor(method="pearson")+
  geom_polygon(data=data.frame(x=c(0.35, 1, 1),
                               y=c(-0.5, -0.5, 0.45)), aes(x=x, y=y),fill="salmon", alpha=0.3)+
  annotate(geom="text", label = paste("β =", round(oars.rapidaim.stability.stats[2,1], digits=2),
                          ", p =", round(oars.rapidaim.stability.stats[2,5], digits=3)),
            x=0.65, y=0.85)+
  facet_wrap(~"Instability Correlation")+
  theme_classic()+theme(legend.position="none",
                        strip.text=element_text(size=10))+
  labs(x="Bray-Curtis Dissimilarity\n(Time 1 vs Time 2)",
       y="ICC RS Z-Scores\n(Time 1 vs Time 2)")
oars.rapidaim.stability.plot

# plot distribution of ICCs
oars.rapidaim.stability.ridges = ggplot(subset(oars.rapidaim.stability,compliant==T),
                                      aes(x=icc, y=1))+
  ggridges::geom_density_ridges(fill="white")+
  coord_flip()+
  theme_minimal()+theme(legend.position="none")
oars.rapidaim.stability.ridges
# plot distribution of Distances
oars.rapidaim.distances.ridges = ggplot(subset(oars.rapidaim.stability,compliant==T),
                                        aes(x=dissimilarity, y=1))+
  ggridges::geom_density_ridges(fill="white")+
  theme_void()+theme(legend.position="none")
oars.rapidaim.distances.ridges

(oars.stool.stability.pcoa.plot+
     oars.rapidaim.stability.plot) %>%
  ggsave(filename="./oars_plots/oars_1b_stability_zone.pdf",
         width=8, height=3.5, device = cairo_pdf)

# The more different 2 microbiomes are (within a patient), the more similar the RS butyrogen Z-score rankings are
# Surprising!



# >> sensitivity ----------------------------------------------------------

# What I'm doing:
# Perturbing (+/-) each taxa % by some % (i.e. change is proportional to the abundance of the taxa)
# X + X * %

upset.values = c(0, 0.001, 0.0025, 0.005, 0.0075, 0.01, 0.025, 0.05, 0.075, 0.1)
iter = 15

skip =T
t0 <- Sys.time()
if(skip == F){
oars.sensitivity.results = 
    do.call(rbind, lapply(upset.values, function(y1){
      do.call(rbind, lapply(1:nrow(oars.rapidaim.list), function(x1){
        print(paste0("disruption = ", y1, ". sample = ", x1))
        
        do.call(rbind, parallel::mclapply(1:iter, function(z1){
          
        # Troubleshoot HM0844.03, HM0874.01, HM0874.04 (row 12 and 41)
        # x1 = 12
        # set variables
        sample = oars.rapidaim.list[x1,]$name
        HM = oars.rapidaim.list[x1,]$HM
        timing = oars.rapidaim.list[x1,]$timing
        seq_run = oars.rapidaim.list[x1,]$seq_run
        
        # > Extract data
        data.subset = oars.rapidaim.data[[x1]][[1]]
        
        # > Upset values by y1 % (each cell independently)
        set.seed(z1)
        
        data.subset.df = phyloseq::otu_table(data.subset)
        
        # disrupt function (multiplies the value by 1 * X %)
        disrupt.fn = function(x) if (x != 0) abs(x + x*as.numeric(sample(c(-1,1)*y1, size=1)) ) else (x)
        # apply function (rescaling is performed later)
        data.subset.otu <- apply(data.subset.df, 1:2, disrupt.fn)
        
        
        # re-make phyloseq object
        data.subset <- phyloseq::merge_phyloseq(phyloseq::otu_table(data.subset.otu, taxa_are_rows = T),
                                                phyloseq::sample_data(data.subset),
                                                phyloseq::tax_table(data.subset),
                                                tree) 
        
        RapidAIM.controls.ps <- phyloseq::merge_phyloseq(data.subset, controls)  # Combine  RapidAIM data with reference cohort
        # includes ALL taxa by default; issue is Tree got dropped when dissecting otu table
        RapidAIM.controls.ps <- phyloseq::transform_sample_counts(RapidAIM.controls.ps, function(x) 100 * x/sum(x)) # TSS scale library sizes to 100%
        
        RapidAIM.controls.ps
        
        # > DISTANCE
        
        RapidAIM.controls.distancematrix <- phyloseq::distance(RapidAIM.controls.ps, method="bray", type="samples") # Calculate Bray-Curtis dissimilarity matrix
        id.control <- as.character(phyloseq::sample_data(phyloseq::subset_samples(RapidAIM.controls.ps, Location=="Distal Colon" & Diagnosis=="Control"))$Sample) # Vector of ONLY reference samples (comparisons will be made to this group)
        id.all <- as.character(phyloseq::sample_data(phyloseq::subset_samples(RapidAIM.controls.ps))$Sample) # Vector of ALL samples
        RapidAIM.controls.distancematrix.matrix <- as.matrix(RapidAIM.controls.distancematrix) # Matrix needed for median-distance code, but not distance-to-centroid
        
        # Calculate Median Distance to Reference Cohort for all samples
        median.distance.df <- do.call(rbind.data.frame,
                                      lapply(id.all, function (x) median(RapidAIM.controls.distancematrix.matrix[x, id.control])))
        names(median.distance.df)[1] <- "Distance" 
        
        RapidAIM.distance.df <- cbind(phyloseq::sample_data(RapidAIM.controls.ps), distance=median.distance.df)
        RapidAIM.distance.df <- RapidAIM.distance.df[grepl(HM, RapidAIM.distance.df$source),]
        RapidAIM.distance.df = RapidAIM.distance.df[, colnames(RapidAIM.distance.df) %in% c("RS_Name", "Replicate", "pH", "Qubit", "source", "timing", "Distance")]
        # Remove scores that will not be included in RS selection analysis
        RapidAIM.distance.df <- subset(RapidAIM.distance.df, !RS_Name %in%
                                         c("Stool", "WaterB", "RawStool", "MED", "ZymoC", "ZymoD","Amioca",
                                           "FOS", "ExtractC", "Spike-in", "FOS-100", "PBS-100", "CLS"))
        RapidAIM.distance.df$Sample = rownames(RapidAIM.distance.df)
        
        # > BUTYROGENS
        
        # EDGE CASE: remove outlier for HM0874.01
        RapidAIM.filtered.4 = phyloseq::subset_samples(data.subset,
                                                       Barcode != "AS16S.101.gg13v5.BC015")
        
        # Convert reads to relative abundances
        RapidAIM.tss  <- phyloseq::transform_sample_counts(RapidAIM.filtered.4, function(x) x / sum(x) )
        # Remove unimportant samples
        RapidAIM.tss <- phyloseq::subset_samples(RapidAIM.tss, !RS_Name %in%
                                                   c("ZymoC", "ZymoD", "NA", "ExtractC", "MED", "Amioca",
                                                     "Spike-in", "WaterB", "FOS-100", "PBS-100", "Stool", "RawStool"))
        
        # Create dataframes for butyrogens derived from Mottawea et al., 2016)
        RapidAIM.butyrogens <- phyloseq::subset_taxa(RapidAIM.tss, Family=="Lachnospiraceae" | Genus=="Blautia" | Genus=="Roseburia" | Genus=="Eubacterium" | Genus=="Ruminococcus" | Genus=="Clostridium" | Genus=="Faecalibacterium") # Note, Lachnospiraceae already includes several genera listed; listed again for clarity
        
        RapidAIM.butyrogens.ps.placeholder <- microbiome::aggregate_taxa(RapidAIM.butyrogens, 'Kingdom', verbose=FALSE) #collapse all butyrogens to one group (kingdom)
        RapidAIM.butyrogens.counts <- as.data.frame(colSums((phyloseq::otu_table(RapidAIM.butyrogens.ps.placeholder))))
        RapidAIM.butyrogens.counts <- RapidAIM.butyrogens.counts %>% rownames_to_column(var="Sample")
        names(RapidAIM.butyrogens.counts)[2] <- "But_Abundance"
        
        # > FINALISE
        
        RapidAIM.final = merge(RapidAIM.distance.df, RapidAIM.butyrogens.counts, by = "Sample") %>% data.frame()
        
        # calculate median per replicate
        
        RapidAIM.final.median = RapidAIM.final %>%
          group_by(RS_Name) %>%
          mutate(med.dis = median(Distance)) %>%
          mutate(med.but = median(But_Abundance)) %>% 
          mutate(med.ph = median(pH)) %>%
          dplyr::select(RS_Name, med.ph, med.dis, med.but) %>% distinct() %>% data.frame()
        
        RapidAIM.final.median$HM = HM
        RapidAIM.final.median$timing = timing
        RapidAIM.final.median$sample = sample
        
        # select RS
        # first, remove PBS and re-add
        RapidAIM.final.median.pbs = subset(RapidAIM.final.median, RS_Name == "PBS")
        RapidAIM.final.median.pbs$Z_score = NA
        
        RapidAIM.final.median = subset(RapidAIM.final.median, RS_Name != "PBS")
        
        # Convert  values to Z-scores to identify top RS
        RapidAIM.final.median$Z_score <- scale(RapidAIM.final.median$med.but, center = TRUE, scale = TRUE)
        candidate.RS <- subset(RapidAIM.final.median, Z_score > 1)
        if(nrow(candidate.RS) >0 ){ # select lowest distance among top butyrogen RS
          top.RS <- subset(candidate.RS, med.dis == min(candidate.RS$med.dis)) # Determine optimal RS by which one has the lowest Distance value
        } else { # or else select top butyrogen RS
          top.RS = subset(RapidAIM.final.median, med.but == max(med.but))
        }
        # re-add PBS
        RapidAIM.final.median = rbind(RapidAIM.final.median, RapidAIM.final.median.pbs) %>% data.frame()
        RapidAIM.final.median$selected = top.RS$RS_Name
        
        
        # > CALCULATE RANKINGS
        RapidAIM.final.median$disrupt = y1
        RapidAIM.final.median$seed = z1
        
        # OUTPUT
        RapidAIM.final.median
        
      }))}))})) # END LOOP
t1 <- Sys.time() # 2 h with parallel
t1 - t0
# Save
saveRDS(oars.sensitivity.results, "./2025_06_23_oars_rapidaim_perturbation_results.Rds")
}
oars.sensitivity.results = readRDS("./2025_06_23_oars_rapidaim_perturbation_results.Rds")

oars.sensitivity.results.df = oars.sensitivity.results

# First, identify the rate by which the RS selected remains the same
empirical.results = oars.rapidaim.scores.pc %>%
  subset(RS_Name == selected) %>%
  dplyr::select(selected, sample, HM, timing) %>% distinct() %>% data.frame()
colnames(empirical.results)[1] = "empirical"

oars.sensitivity.results.df = merge(oars.sensitivity.results.df, 
                                    empirical.results,
                                    by=c("HM", "timing"))

oars.sensitivity.results.rates = oars.sensitivity.results.df %>%
  subset(RS_Name == selected) %>%
  group_by(disrupt, seed, HM) %>%
  mutate(same = ifelse(selected == empirical, TRUE, FALSE)) %>%
  group_by(disrupt, seed) %>%
  summarize(
    same.rate = sum(same == TRUE)) %>%
  mutate(same.rate.perc = same.rate / nrow(empirical.results))


# stats
oars.sens.rates.main <- subset(oars.sensitivity.results.rates, disrupt <= 0.1) # remove values greater than 10% disruption
oars.sens.rates.main.lm <- do.call(rbind, lapply(unique(oars.sens.rates.main$disrupt)[-1], function(x) {
  # x = unique(oars.sens.rates.main$disrupt)[2]
  # lm
  lm.results = glm(same.rate.perc ~ disrupt,
                   #binomial(link = "logit"),
                   data=subset(oars.sens.rates.main, disrupt %in% c(min(oars.sens.rates.main$disrupt), x))) %>% summary() %>% coef()
  data.frame(disrupt = x,
             pval = lm.results[2,4],
             coef = lm.results[2,1])
}))
oars.sens.rates.main.lm$padj <- p.adjust(oars.sens.rates.main.lm$pval, method="bonferroni")
oars.sens.rates.main.lm$sig = ifelse(oars.sens.rates.main.lm$padj < 0.05, "*", "")
oars.sens.rates.main.lm

oars.sensitivity.results.rates.plot = 
  ggplot(subset(oars.sensitivity.results.rates, disrupt <= 0.1),
         aes(x=as.factor(disrupt*100), y=same.rate.perc*100))+
  geom_point(position=position_jitter(width=0.1, height = 0.0, seed=1), 
             color="black",
             aes(fill=same.rate.perc),
             size=2, shape=21)+
  #geom_violin(draw_quantiles=c(0.05, 0.50, 0.95),alpha=0.5)+
  geom_boxplot(alpha=0.5, outlier.shape=NA, width=0.5)+
  geom_hline(yintercept=95, color="red", linetype=2, alpha=0.5)+
  scale_y_continuous(breaks=seq(from=70, to=100, by=5))+
  #geom_text(data=oars.sens.rates.main.lm,
  #          aes(x=as.factor(disrupt*100), y=101, label=sig),
  #          size=7)+
  scale_fill_gradient2(low = "red", high="blue", midpoint = 0.9)+
  facet_wrap(~"RS Selection Sensitivity")+
  theme_classic()+theme(legend.position="none",
                        strip.text=element_text(size=10))+
  labs(x="OTU % Perturbation", y="Unchanged RS Selection %")
oars.sensitivity.results.rates.plot



# >> butyrogens vs pH -----------------------------------------------------

# note: much of this is hypothesis generating and negative results

# correlate change in pH with butyrogens (only)
oars.rapidaim.scores.pc.full = readRDS("./2025_06_09_oars_scores.Rds")
# calculate delta pH
oars.rapidaim.scores.pc.full = oars.rapidaim.scores.pc.full %>%
  group_by(sample) %>%
  mutate(delta.ph = med.ph - med.ph[RS_Name == "PBS"])

# identify microbiomes with + or - correlation between delta pH and butyrogens
oars.rapidaim.ph.but.cor = do.call(rbind, lapply(unique(oars.rapidaim.scores.pc.full$sample), function(hm){
  data.subset = subset(oars.rapidaim.scores.pc.full, sample == hm)
  data.frame(HM = unique(data.subset$HM),
             timing = unique(data.subset$timing),
             sample = hm,
             cor = cor.test(data.subset$med.but, data.subset$med.ph, method="spearman")$estimate,
             cor.sign = sign(cor.test(data.subset$med.but, data.subset$med.ph, method="spearman")$estimate),
             pval = cor.test(data.subset$med.but, data.subset$med.ph, method="spearman")$p.value)
}))

oars.rapidaim.ph.but.cor %>% arrange(pval)

oars.rapidaim.scores.pc.full$RS_Name = factor(oars.rapidaim.scores.pc.full$RS_Name, levels=rs.names.pbs)
oars.rapidaim.ph.butyrogens.plot = ggplot(subset(oars.rapidaim.scores.pc.full, timing %in% c("0M", "3M", "6M")),
       aes(x=med.but*100, y=med.ph))+
  # make transparent if below R value
  geom_point(shape=21, aes(fill=RS_Name))+
  geom_smooth(method="lm", color="black", se=F)+
  ggpubr::stat_cor(aes(label = ..r.label..), 
                   method="spearman", 
                   label.y=8.3, vjust=0.75,
                   label.x.npc=0.5, hjust=0.5)+
  scale_x_log10()+
  scale_fill_manual(values=labelcolors$cols[c(10,1:9)])+
  facet_grid(timing~HM,scales="free")+
  theme_minimal()+theme(legend.position="none",
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="Median Butyrogens %",
       y="Median pH")
oars.rapidaim.ph.butyrogens.plot

# Compress this into a heatmap;
# select 2 microbiomes to show examples

oars.rapidaim.ph.but.cor.data = subset(oars.rapidaim.scores.pc.full, timing %in% c("0M", "3M", "6M"))

oars.rapidaim.ph.but.cor.data = do.call(rbind, lapply(unique(oars.rapidaim.ph.but.cor.data$sample), function(hm){
  data.subset = subset(oars.rapidaim.ph.but.cor.data, sample == hm) 
  # run cor
  data.frame(cor = cor.test(data.subset$med.but, data.subset$med.ph, method="spearman")$estimate,
             sample = hm)
}))
oars.rapidaim.ph.but.cor.data = merge(oars.rapidaim.ph.but.cor.data,
                                      subset(oars.rapidaim.scores.pc.full, timing %in% c("0M", "3M", "6M"))[,c("sample", "timing", "HM")]%>%distinct,
                                      by="sample")
oars.rapidaim.ph.but.cor.plot = ggplot(oars.rapidaim.ph.but.cor.data %>% mutate(HM = factor(HM, levels=rev(sort(unique(HM))))),
       aes(x=timing, y=HM))+
  geom_tile(aes(fill=cor), color="white")+
  theme_minimal()+
  geom_text(aes(label=round(cor, digits=2), color=ifelse(abs(cor) > 0.5, "black","white")), size=3)+
  scale_fill_gradient2(low="blue", high="red")+
  scale_color_manual(values=c("white", "black"))+
  guides(color="none")+
  theme(panel.grid=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        legend.position="right",
        strip.text = element_text(size=10),
        strip.background = element_rect(
          color="black"))+
  facet_wrap(~timing, scales="free_x")+
  labs(fill = "Spearman ρ", color=NA)
oars.rapidaim.ph.but.cor.plot

oars.rapidaim.ph.butyrogens.plot.819 = ggplot(subset(oars.rapidaim.scores.pc.full, sample %in% c("HM0819.03", "HM0819.04")),
                                          aes(x=med.but*100, y=med.ph))+
  # make transparent if below R value
  geom_point(shape=21, aes(fill=RS_Name), size=3)+
  geom_smooth(method="lm", color="black", se=F)+
  ggpubr::stat_cor(aes(label = ..r.label..), 
                   method="spearman", 
                   label.y=7.6, vjust=0.75,
                   label.x.npc=0.5, hjust=0.5, size=3)+
  scale_x_log10()+
  scale_fill_manual(values=labelcolors$cols[c(10,1:9)])+
  facet_wrap(HM~timing,scales="free")+
  theme_minimal()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="Median Butyrogens %",
       y="Median pH")
oars.rapidaim.ph.butyrogens.plot.819


oars.rapidaim.ph.butyrogens.plot.819+oars.rapidaim.ph.but.cor.plot


table(subset(oars.rapidaim.ph.but.cor, timing %in% c("0M", "3M", "6M"))$cor.sign)

# rarefy
redo = F
if(redo == T){
set.seed(25)
oars.rapidaim.phyloseq.rare = phyloseq::rarefy_even_depth(oars.rapidaim.phyloseq, sample.size=50000, replace=F)
# save
saveRDS(oars.rapidaim.phyloseq.rare, "./oars.rapidaim.phyloseq.rare.Rds")
}
oars.rapidaim.phyloseq.rare = readRDS("./oars.rapidaim.phyloseq.rare.Rds")


# butyrogens
oars.rapidaim.phyloseq.butyrogens = speedyseq::psmelt(oars.rapidaim.phyloseq.rare)
oars.rapidaim.phyloseq.butyrogens = subset(oars.rapidaim.phyloseq.butyrogens, Family=="Lachnospiraceae" | Genus=="Blautia" | Genus=="Roseburia" | Genus=="Eubacterium" | Genus=="Ruminococcus" | Genus=="Clostridium" | Genus=="Faecalibacterium")
# make taxa names
oars.rapidaim.phyloseq.butyrogens$taxa = gsub("_NA", "", ifelse(!is.na(oars.rapidaim.phyloseq.butyrogens$Genus), 
                                                paste(oars.rapidaim.phyloseq.butyrogens$Genus, oars.rapidaim.phyloseq.butyrogens$Species, sep="_"),
                                                oars.rapidaim.phyloseq.butyrogens$Family))
# make matrix
oars.rapidaim.phyloseq.butyrogens$Abundance = oars.rapidaim.phyloseq.butyrogens$Abundance / 50000
length(unique(oars.rapidaim.phyloseq.butyrogens$taxa))
# make the dataframe a bit smaller
oars.rapidaim.phyloseq.butyrogens = subset(oars.rapidaim.phyloseq.butyrogens, RS_Name %in% rs.names.pbs)
oars.rapidaim.phyloseq.butyrogens.mat = reshape2::acast(oars.rapidaim.phyloseq.butyrogens,
                                                        Sample ~ OTU, value.var="Abundance")

# create tax_table with OTU and taxa name
# first, print out full GG name
oars.rapidaim.phyloseq.butyrogens = oars.rapidaim.phyloseq.butyrogens %>%
  mutate(full_name = 
  paste(paste("k__", Kingdom, sep=""), paste("p__", Phylum, sep=""), paste("c__", Class, sep=""), 
        paste("o__", Order, sep=""), paste("f__", Family, sep=""), paste("g__", Genus, sep=""), paste("s__", Species, sep=""), sep=";"))
# replace OTU with taxa name
oars.rapidaim.otu.names = data.frame(
  full_name = oars.rapidaim.phyloseq.butyrogens[match(colnames(oars.rapidaim.phyloseq.butyrogens.mat), oars.rapidaim.phyloseq.butyrogens$OTU),]$full_name,
  taxa = oars.rapidaim.phyloseq.butyrogens[match(colnames(oars.rapidaim.phyloseq.butyrogens.mat), oars.rapidaim.phyloseq.butyrogens$OTU),]$taxa) %>%
  mutate(OTU = colnames(oars.rapidaim.phyloseq.butyrogens.mat)) %>% distinct() %>%
  mutate(name = make.unique(taxa))
# make names unique

# replace OTU with taxa name
colnames(oars.rapidaim.phyloseq.butyrogens.mat) = oars.rapidaim.otu.names[match(colnames(oars.rapidaim.phyloseq.butyrogens.mat), oars.rapidaim.otu.names$OTU),]$name

# make mapping file
oars.rapidaim.phyloseq.butyrogens.map = oars.rapidaim.phyloseq.butyrogens[,c("Barcode", "RS_Name", "HM", "pH")] %>% distinct() %>% data.frame()
oars.rapidaim.phyloseq.butyrogens.map$sample = paste("HM", oars.rapidaim.phyloseq.butyrogens.map$HM, sep="")
oars.rapidaim.phyloseq.butyrogens.map = merge(oars.rapidaim.phyloseq.butyrogens.map,
                                              oars.rapidaim.scores.pc.full[,c("sample", "RS_Name","timing", "delta.ph")], by=c("sample", "RS_Name"))
rownames(oars.rapidaim.phyloseq.butyrogens.map) = oars.rapidaim.phyloseq.butyrogens.map$Barcode
# run maaslin2

oars.rapidaim.phyloseq.butyrogens.maaslin = Maaslin2::Maaslin2(input_data = oars.rapidaim.phyloseq.butyrogens.mat,
                                          input_metadata = subset(oars.rapidaim.phyloseq.butyrogens.map, timing %in% c("0M", "3M", "6M")),
                                          output = "~/Downloads",
                                          fixed_effects = c("delta.ph"),  # Example fixed effects
                                          random_effects = c("HM", "RS_Name"),       # Example random effects
                                          normalization = "TSS",                       # Total Sum Scaling normalization
                                          transform = "LOG",                           # Log transformation
                                          analysis_method = "LM",                      # Linear model
                                          plot_scatter = FALSE,                        # Disable scatterplot generation
                                          plot_heatmap = FALSE,                         # Keep heatmap generation (optional)
                                          max_significance = 0.05,                     # Significance threshold for q-values
                                          standardize = TRUE                           # Disable standardization (optional)
)
oars.rapidaim.phyloseq.butyrogens.maaslin.df = oars.rapidaim.phyloseq.butyrogens.maaslin$results
oars.rapidaim.phyloseq.butyrogens.maaslin.df$feature = ifelse(grepl("X.", oars.rapidaim.phyloseq.butyrogens.maaslin.df$feature),
                                                              gsub("X.Ruminococcus.", "Ruminococcus", oars.rapidaim.phyloseq.butyrogens.maaslin.df$feature),
                                                              oars.rapidaim.phyloseq.butyrogens.maaslin.df$feature)


oars.rapidaim.phyloseq.butyrogens.maaslin.volcano = ggplot(oars.rapidaim.phyloseq.butyrogens.maaslin.df,
       aes(x=coef, y=-log10(qval)))+
  geom_point(shape=21, aes(fill=coef))+
  scale_fill_gradient2(low="blue", high="red")+
  ggrepel::geom_text_repel(aes(label = ifelse(qval < 0.05, gsub("\\..*", "", feature), NA)), size=3)+
  theme_minimal()+theme(legend.position="none")+
  labs(x="Adjusted Coefficient",
       y="-log10(padj)")
oars.rapidaim.phyloseq.butyrogens.maaslin.volcano

# loop through butyrogens and see ratios of up vs down with delta ph
oars.rapidaim.phyloseq.butyrogens.maaslin.sign = do.call(rbind, lapply(unique(gsub("\\..*", "", oars.rapidaim.phyloseq.butyrogens.maaslin.df$feature)), function(taxa){
  print(taxa)
  # taxa = unique(gsub("\\..*", "", oars.rapidaim.phyloseq.butyrogens.maaslin.df$feature))[1]
  # taxa = "Oribacterium"
  data.subset = subset(oars.rapidaim.phyloseq.butyrogens.maaslin.df, gsub("\\..*", "", feature) == taxa ) %>% 
    subset(qval < 0.05) %>%
    mutate(coef.sign = sign(coef)) %>% dplyr::select(coef.sign) %>% table() %>% data.frame()
  if(nrow(data.subset) > 0){
  data.subset$taxa = taxa
  data.subset
  }else{
    data.frame(coef.sign = NA,
               Freq = NA,
               taxa = taxa)
  }
}))
oars.rapidaim.phyloseq.butyrogens.maaslin.sign$sign = ifelse(oars.rapidaim.phyloseq.butyrogens.maaslin.sign$coef.sign == 1, 1, -1)

oars.rapidaim.phyloseq.butyrogens.maaslin.sign.plot = ggplot(na.omit(oars.rapidaim.phyloseq.butyrogens.maaslin.sign),
       aes(x=Freq*sign, y=reorder(taxa, Freq)))+
  geom_col(aes(group=taxa), fill="white", alpha=1, linetype=2)+
  geom_col(aes(group=taxa, fill=sign), alpha=0.5, linetype=2)+
  geom_text(aes(label = Freq, hjust=ifelse(sign > 0, -0.5, 1.5)), size=2.5)+
  geom_vline(xintercept=0, linetype=1, alpha=1)+
  scale_fill_gradient(low=("blue"), high="red")+
  scale_x_continuous(limits=c(min(na.omit(oars.rapidaim.phyloseq.butyrogens.maaslin.sign)$Freq*-1)-10,
                              max(na.omit(oars.rapidaim.phyloseq.butyrogens.maaslin.sign)$Freq)+10))+
  theme_minimal()+theme(legend.position="none")+
  labs(x="Number of Significant OTU-pH Associations",
       y="",
       )
oars.rapidaim.phyloseq.butyrogens.maaslin.sign.plot

oars.rapidaim.ph.butyrogens.plot
oars.rapidaim.phyloseq.butyrogens.maaslin.volcano+oars.rapidaim.phyloseq.butyrogens.maaslin.sign.plot


# :: per correlation group: negative ------------------------------------------------

# subset microbiomes to sig (pval < 0.05) + and - correlations
# repeat maaslin per group

oars.rapidaim.ph.but.cor


oars.rapidaim.ph.negative.maaslin = Maaslin2::Maaslin2(input_data = oars.rapidaim.phyloseq.butyrogens.mat,
                                                               input_metadata = subset(oars.rapidaim.phyloseq.butyrogens.map, 
                                                                                       timing %in% c("0M", "3M", "6M") & 
                                                                                         sample %in% subset(oars.rapidaim.ph.but.cor, cor > 0 & pval < 0.05)$sample),
                                                               output = "~/Downloads",
                                                               fixed_effects = c("delta.ph"),  # Example fixed effects
                                                               random_effects = c("HM", "RS_Name"),       # Example random effects
                                                               normalization = "TSS",                       # Total Sum Scaling normalization
                                                               transform = "LOG",                           # Log transformation
                                                               analysis_method = "LM",                      # Linear model
                                                               plot_scatter = FALSE,                        # Disable scatterplot generation
                                                               plot_heatmap = FALSE,                         # Keep heatmap generation (optional)
                                                               max_significance = 0.05,                     # Significance threshold for q-values
                                                               standardize = TRUE                           # Disable standardization (optional)
)
oars.rapidaim.ph.negative.maaslin.df = oars.rapidaim.ph.negative.maaslin$results
oars.rapidaim.ph.negative.maaslin.df$feature = ifelse(grepl("X.", oars.rapidaim.ph.negative.maaslin.df$feature),
                                                              gsub("X.Ruminococcus.", "Ruminococcus", oars.rapidaim.ph.negative.maaslin.df$feature),
                                                              oars.rapidaim.ph.negative.maaslin.df$feature)


oars.rapidaim.ph.negative.maaslin.volcano = ggplot(oars.rapidaim.ph.negative.maaslin.df,
                                                           aes(x=coef, y=-log10(qval)))+
  geom_point(shape=21, aes(fill=coef))+
  scale_fill_gradient2(low="blue", high="red")+
  ggrepel::geom_text_repel(aes(label = ifelse(qval < 0.05, gsub("\\..*", "", feature), NA)))+
  theme_minimal()+theme(legend.position="none")+
  labs(x="Adjusted Coefficient",
       y="-log10(padj)")
oars.rapidaim.ph.negative.maaslin.volcano

# loop through butyrogens and see ratios of up vs down with delta ph
oars.rapidaim.ph.negative.maaslin.sign = do.call(rbind, lapply(unique(gsub("\\..*", "", oars.rapidaim.ph.negative.maaslin.df$feature)), function(taxa){
  print(taxa)
  # taxa = unique(gsub("\\..*", "", oars.rapidaim.ph.negative.maaslin.df$feature))[1]
  # taxa = "Oribacterium"
  data.subset = subset(oars.rapidaim.ph.negative.maaslin.df, gsub("\\..*", "", feature) == taxa ) %>% 
    subset(qval < 0.05) %>%
    mutate(coef.sign = sign(coef)) %>% dplyr::select(coef.sign) %>% table() %>% data.frame()
  if(nrow(data.subset) > 0){
    data.subset$taxa = taxa
    data.subset
  }else{
    data.frame(coef.sign = NA,
               Freq = NA,
               taxa = taxa)
  }
}))
oars.rapidaim.ph.negative.maaslin.sign$sign = ifelse(oars.rapidaim.ph.negative.maaslin.sign$coef.sign == 1, 1, -1)

oars.rapidaim.ph.negative.maaslin.sign.plot = ggplot(na.omit(oars.rapidaim.ph.negative.maaslin.sign),
                                                             aes(x=Freq*sign, y=reorder(taxa, Freq)))+
  geom_col(aes(group=taxa), fill="white", alpha=1, linetype=2)+
  geom_col(aes(group=taxa, fill=sign), alpha=0.5, linetype=2)+
  geom_text(aes(label = Freq, hjust=ifelse(sign > 0, -0.5, 1.5)), size=2.5)+
  geom_vline(xintercept=0, linetype=1, alpha=1)+
  scale_fill_gradient(low=("blue"), high="red")+
  scale_x_continuous(limits=c(min(na.omit(oars.rapidaim.ph.negative.maaslin.sign)$Freq*-1)-10,
                              max(na.omit(oars.rapidaim.ph.negative.maaslin.sign)$Freq)+10))+
  theme_minimal()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5))+
  labs(x="Number of Significant OTU-pH Associations (per direction)",
       y="",
       title="Microbiomes with\nNegative pH-Butyrogen Correlation")
oars.rapidaim.ph.negative.maaslin.sign.plot

# :: per correlation group: positive ------------------------------------------------

# subset microbiomes to sig (pval < 0.05) + and - correlations
# repeat maaslin per group

oars.rapidaim.ph.but.cor

oars.rapidaim.ph.positive.maaslin = Maaslin2::Maaslin2(input_data = oars.rapidaim.phyloseq.butyrogens.mat,
                                                       input_metadata = subset(oars.rapidaim.phyloseq.butyrogens.map, 
                                                                               timing %in% c("0M", "3M", "6M") & 
                                                                                 sample %in% subset(oars.rapidaim.ph.but.cor, cor < 0 & pval < 0.05)$sample),
                                                       output = "~/Downloads",
                                                       fixed_effects = c("delta.ph"),  # Example fixed effects
                                                       random_effects = c("HM", "RS_Name"),       # Example random effects
                                                       normalization = "TSS",                       # Total Sum Scaling normalization
                                                       transform = "LOG",                           # Log transformation
                                                       analysis_method = "LM",                      # Linear model
                                                       plot_scatter = FALSE,                        # Disable scatterplot generation
                                                       plot_heatmap = FALSE,                         # Keep heatmap generation (optional)
                                                       max_significance = 0.05,                     # Significance threshold for q-values
                                                       standardize = TRUE                           # Disable standardization (optional)
)
oars.rapidaim.ph.positive.maaslin.df = oars.rapidaim.ph.positive.maaslin$results
oars.rapidaim.ph.positive.maaslin.df$feature = ifelse(grepl("X.", oars.rapidaim.ph.positive.maaslin.df$feature),
                                                      gsub("X.Ruminococcus.", "Ruminococcus", oars.rapidaim.ph.positive.maaslin.df$feature),
                                                      oars.rapidaim.ph.positive.maaslin.df$feature)


oars.rapidaim.ph.positive.maaslin.volcano = ggplot(oars.rapidaim.ph.positive.maaslin.df,
                                                   aes(x=coef, y=-log10(qval)))+
  geom_point(shape=21, aes(fill=coef))+
  scale_fill_gradient2(low="blue", high="red")+
  ggrepel::geom_text_repel(aes(label = ifelse(qval < 0.05, gsub("\\..*", "", feature), NA)))+
  theme_minimal()+theme(legend.position="none")+
  labs(x="Adjusted Coefficient",
       y="-log10(padj)")
oars.rapidaim.ph.positive.maaslin.volcano

# loop through butyrogens and see ratios of up vs down with delta ph
oars.rapidaim.ph.positive.maaslin.sign = do.call(rbind, lapply(unique(gsub("\\..*", "", oars.rapidaim.ph.positive.maaslin.df$feature)), function(taxa){
  print(taxa)
  # taxa = unique(gsub("\\..*", "", oars.rapidaim.ph.positive.maaslin.df$feature))[1]
  # taxa = "Oribacterium"
  data.subset = subset(oars.rapidaim.ph.positive.maaslin.df, gsub("\\..*", "", feature) == taxa ) %>% 
    subset(qval < 0.05) %>%
    mutate(coef.sign = sign(coef)) %>% dplyr::select(coef.sign) %>% table() %>% data.frame()
  if(nrow(data.subset) > 0){
    data.subset$taxa = taxa
    data.subset
  }else{
    data.frame(coef.sign = NA,
               Freq = NA,
               taxa = taxa)
  }
}))
oars.rapidaim.ph.positive.maaslin.sign$sign = ifelse(oars.rapidaim.ph.positive.maaslin.sign$coef.sign == 1, 1, -1)


# :: per correlation group: both  -----------------------------------------

# plot both as facet_Wrap
# combine datasets
oars.rapidaim.ph.both.maaslin.sign = rbind(oars.rapidaim.ph.positive.maaslin.sign %>% mutate(Correlation = "Positive"),
                                           oars.rapidaim.ph.negative.maaslin.sign %>% mutate(Correlation = "Negative")) %>% data.frame()


oars.rapidaim.ph.both.maaslin.sign.plot = ggplot(na.omit(oars.rapidaim.ph.both.maaslin.sign)%>%mutate(Correlation = factor(Correlation, levels=c("Positive", "Negative"))),
                                                     aes(x=Freq*sign, y=reorder(taxa, Freq)))+
  geom_col(aes(group=taxa), fill="white", alpha=1, linetype=2)+
  geom_col(aes(group=taxa, fill=sign), alpha=0.5, linetype=2)+
  geom_text(aes(label = Freq, hjust=ifelse(sign > 0, -0.5, 1.5)), size=2.5)+
  geom_vline(xintercept=0, linetype=1, alpha=1)+
  scale_fill_gradient(low=("blue"), high="red")+
  scale_x_continuous(limits=c(min(na.omit(oars.rapidaim.ph.positive.maaslin.sign)$Freq*-1)-10,
                              max(na.omit(oars.rapidaim.ph.positive.maaslin.sign)$Freq)+10))+
  labs(x="Number of Significant OTU-pH Associations (per direction)",
       y="")+
  facet_wrap(~Correlation,labeller = label_both)+
  theme_minimal()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5),
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))
oars.rapidaim.ph.both.maaslin.sign.plot


# :: Fisher test ----------------------------------------------------------

# are they balanced?
subset(oars.rapidaim.ph.but.cor, pval < 0.05 & timing %in% c("0M", "3M", "6M"))$cor.sign %>% table()
# 8 vs 4

# loop through fisher test
oars.rapidaim.ph.both.maaslin.sign.stats = do.call(rbind, lapply(unique(oars.rapidaim.ph.both.maaslin.sign$taxa), function(x){
  print(x)
  # x = "Faecalibacterium_prausnitzii"
  data.subset = subset(oars.rapidaim.ph.both.maaslin.sign, taxa == x)
  data.subset = na.omit(data.subset)
  # manually reconstruct contingency table, since many values might be missing
  if(sum(data.subset$coef.sign == 1 & data.subset$Correlation == "Positive")==0){
    condition1 = data.frame(coef.sign = 1, Correlation = "Positive", Freq = 0)
  } else {
    condition1 = data.frame(coef.sign = 1, Correlation = "Positive", Freq = subset(data.subset, coef.sign == 1 & Correlation == "Positive")$Freq)
  }
  if(sum(data.subset$coef.sign == -1 & data.subset$Correlation == "Positive")==0){
    condition2 = data.frame(coef.sign = -1, Correlation = "Positive", Freq = 0)
  } else {
    condition2 = data.frame(coef.sign = -1, Correlation = "Positive", Freq = subset(data.subset, coef.sign == -1 & Correlation == "Positive")$Freq)
  }
  if(sum(data.subset$coef.sign == 1 & data.subset$Correlation == "Negative")==0){
    condition3 = data.frame(coef.sign = 1, Correlation = "Negative", Freq = 0)
  } else {
    condition3 = data.frame(coef.sign = 1, Correlation = "Negative", Freq = subset(data.subset, coef.sign == 1 & Correlation == "Negative")$Freq)
  }
  if(sum(data.subset$coef.sign == -1 & data.subset$Correlation == "Negative")==0){
    condition4 = data.frame(coef.sign = -1, Correlation = "Negative", Freq = 0)
  } else {
    condition4 = data.frame(coef.sign = -1, Correlation = "Negative", Freq = subset(data.subset, coef.sign == -1 & Correlation == "Negative")$Freq)
  }
  # merge
  data.subset = rbind(condition1, condition2, condition3, condition4) %>% data.frame() %>% 
    reshape2::acast(coef.sign ~ Correlation, value.var="Freq")
  # run test
  test.result = fisher.test(data.subset)
  data.frame(taxa = x,
             pval =test.result$p.value)
})) %>% mutate(padj = p.adjust(pval, method = "BH")) %>% arrange(padj)

oars.rapidaim.ph.both.maaslin.sign.stats

# F. prausnitzii is significantly negatively associated with pH in negative correlators; 
# positively correlated with ph in positive correlators

# What is special about F. prausnitzii? Absolute pH differences?

# Identify specific ASVs sig in Positive and Negative groups


# :: F. prausnitzii -------------------------------------------------------

oars.fprau.positive = subset(oars.rapidaim.ph.positive.maaslin.df, qval < 0.05 & grepl("Faecalibacterium_prausnitzii", feature)) %>%
  dplyr::select(feature, coef, qval) %>% distinct() %>% mutate(Correlation = "Positive")
oars.fprau.negative = subset(oars.rapidaim.ph.negative.maaslin.df, qval < 0.05 & grepl("Faecalibacterium_prausnitzii", feature)) %>%
  dplyr::select(feature, coef, qval) %>% distinct() %>% mutate(Correlation = "Negative")

oars.fprau.important = rbind(oars.fprau.positive,
                             oars.fprau.negative) %>%
  data.frame() %>%
  group_by(feature) %>% mutate(row.order = prod(abs(coef)*100))

oars.fprau.important.coef.plot = ggplot(oars.fprau.important,
       aes(x=coef, y=reorder(gsub("Faecalibacterium.prausnitzii", "F.p", feature), coef)))+
  geom_vline(xintercept=0, linetype=2, alpha=0.5)+
  geom_path(aes(group=feature), linetype=1, color="black")+
  geom_point(aes(fill=coef, shape=Correlation), size=3)+
  scale_shape_manual(values=c(25,24))+
  scale_fill_gradient2(low="blue", high="red")+
  theme_minimal()+
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5, size=12))+
  labs(x="Adjusted Coefficient",y="")
oars.fprau.important.coef.plot

# phylogenetic tree
# USING OTUs

# isolate F. prausintzii + OTU IDs
oars.rapidaim.otu.names.fprau = subset(oars.rapidaim.otu.names,
                                       name %in% oars.fprau.important$feature)
# prune full gg tree
oars.rapidaim.otu.tree.fprau = phyloseq::prune_taxa(oars.rapidaim.otu.names.fprau$OTU,
                     oars.rapidaim.phyloseq)

tree <- phyloseq::phy_tree(oars.rapidaim.otu.tree.fprau)  # ps is your phyloseq object
tree$node.label = oars.rapidaim.otu.names.fprau$name[match(tree$tip.label, oars.rapidaim.otu.names.fprau$OTU)]
tree$tip.label = oars.rapidaim.otu.names.fprau$name[match(tree$tip.label, oars.rapidaim.otu.names.fprau$OTU)]

# build tree meta
tree_meta = merge(oars.rapidaim.otu.names.fprau[,c("name", "OTU")],
                  oars.fprau.important%>%mutate(name = feature), by="name")

# shorten longest tip
# Identify tip edges
tip_edges <- which(tree$edge[,2] <= ape::Ntip(tree))

# Cap long tip lengths
tree$edge.length[tip_edges] <- pmin(tree$edge.length[tip_edges], 0.02)

library("ggtree");library("ggtreeExtra")


fp.p0 <- ggtree(tree, layout = "fan", open.angle = 45, size=0.5)+
  geom_tiplab(
    aes(label = paste(" ", gsub("Faecalibacterium_prausnitzii", "FP", label), " ")),  # Use feature as label
    offset = 0.02,         # Adjust to place labels beyond bars (0.05 + 0.15 + extra)
    size = 3,              # Match size of previous attempts
    color = "black",       # Consistent color
    align = TRUE,          # Align labels at equal radius
    linetype = 3,          # No connecting lines
    linewidth=5
  )
fp.p0

fp.p0 <- rotate_tree(fp.p0, angle=80)

# Heatmap tiles
fp.p1 <- fp.p0 + geom_fruit(
  data = tree_meta,
  geom = geom_tile,
  mapping = aes(y = feature, x = Correlation, fill = (coef)),
  color = "white",
  offset = -0.05,
  pwidth = 0.25,
  axis.params = list(axis = "none")
) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  #scale_color_manual(values=c("blue", "red", "white"))+
  labs(fill = "Adjusted\nCoefficient")
fp.p1

margins = c(0,0,0,0)

# save all
fp.p2 = fp.p1 + 
  theme(plot.margin = margin(margins, "cm"))+ 
  expand_limits(x = c(0, 0.11))
fp.p2

ggsave(fp.p2,
       filename = "./oars_plots/oars_1_otu_fprau.pdf",
       width=4, height=4, device = cairo_pdf)

# no longer sig (after FDR)


# >> Reproducibility ------------------------------------------------------

# note: moved code from oars_validation_analysis


# :: load original data ---------------------------------------------------

oars.rapidaim.scores.val = readRDS("./2025_06_09_oars_scores.Rds")
# make column names amenable
colnames(oars.rapidaim.scores.val) = c("RS_Name", "med.ph", "med.dis", "med.but", "HMcode", "timing", "HM", "zscore", "z.selected")
# calculate pH selection
oars.rapidaim.scores.val = oars.rapidaim.scores.val %>%
  group_by(HM) %>%
  mutate(ph.selected = RS_Name[med.ph == min(med.ph)]) %>%
  mutate(Attempt = 0)
# clean HM
oars.rapidaim.scores.val$HM = gsub("HM", "", oars.rapidaim.scores.val$HM)
# subset to 819 and 932
oars.rapidaim.scores.val = subset(oars.rapidaim.scores.val,
                                  HM %in% c("0819.03", "0819.04", "0932.01", "0932.02"))


# :: Repro - Minimal Time -------------------------------------------------

# here are the first replicates
plate1_scfa_rep1 = read.csv("~/Documents/PhD/For others/oars_validation/2025_03_09_oars_validation_plate_map_819.03_819.04.csv")
plate2_scfa_rep1 = read.csv("~/Documents/PhD/For others/oars_validation/2025_03_09_oars_validation_plate_map_932.01_932.02.csv")
colnames(plate1_scfa_rep1)[1] <- "Barcode"
colnames(plate2_scfa_rep1)[1] <- "Barcode"

# here are the back-to-back replicates
plate1_scfa_rep4 = read.csv("~/Documents/PhD/For others/oars_validation/2025_06_03_oars_validation_plate_map_819.03_819.04_4_b2b1.csv")
plate1_scfa_rep5 = read.csv("~/Documents/PhD/For others/oars_validation/2025_06_03_oars_validation_plate_map_819.03_819.04_5_b2b2.csv")
colnames(plate1_scfa_rep4)[1] <- "Barcode"
colnames(plate1_scfa_rep5)[1] <- "Barcode"


plates_scfa = rbind(plate1_scfa_rep1,
                    plate2_scfa_rep1,
                    plate1_scfa_rep4, 
                    plate1_scfa_rep5) %>% data.frame()

plates_scfa$RS_Name = factor(plates_scfa$RS_Name, levels=c("PBS", "Authentic", "BobsRedMill","MSPrebiotic", "LetsDoOrganic","HiMaize260","Novelose330","ActistarRT","FibersymRW","Versafibe1490","RS"))

# fix HM names (to be compatible with OTU data)
plates_scfa$HM = ifelse(substr(plates_scfa$HM, 1,1) == "0", 
                        plates_scfa$HM,
                        paste("0", plates_scfa$HM, sep=""))

# :: Load OTUs ------------------------------------------------------------

# bioms
plate1_otu_rep1 <- phyloseq::import_biom("~/Documents/PhD/For others/oars_validation/AS16S_370_gg13v5_no_doubletons.biom", parseFunction = phyloseq::parse_taxonomy_greengenes) #import biom file
plate2_otu_rep1 <- phyloseq::import_biom("~/Documents/PhD/For others/oars_validation/AS16S_371_gg13v5_no_doubletons.biom", parseFunction = phyloseq::parse_taxonomy_greengenes) #import biom file

plate1_otu_rep4 <- phyloseq::import_biom("~/Documents/PhD/For others/oars_validation/AS16S_380_gg13v5_no_doubletons.biom", parseFunction = phyloseq::parse_taxonomy_greengenes) #import biom file
plate1_otu_rep5 <- phyloseq::import_biom("~/Documents/PhD/For others/oars_validation/AS16S_381_gg13v5_no_doubletons.biom", parseFunction = phyloseq::parse_taxonomy_greengenes) #import biom file


# merge
plate_otu = phyloseq::merge_phyloseq(plate1_otu_rep1,
                                     plate2_otu_rep1)

plate_otu = phyloseq::merge_phyloseq(plate_otu,
                                     plate1_otu_rep4,
                                     plate1_otu_rep5)

rownames(plates_scfa) <- make.unique(plates_scfa$Barcode)

plates_scfa$microbiome = paste(plates_scfa$HM, plates_scfa$Attempt, sep="_")
unique(plates_scfa$microbiome)

plate_otu_ps = phyloseq::merge_phyloseq(plate_otu,
                                        phyloseq::sample_data(plates_scfa))
# plot reads
ggplot(merge(data.frame(depth = colSums(phyloseq::otu_table(plate_otu_ps)),
                        Barcode = rownames(data.frame(colSums(phyloseq::otu_table(plate_otu_ps))))),
             plates_scfa, by="Barcode"),
       aes(x=Treatment, y=depth))+
  geom_point(shape=21, aes(fill=Treatment))+
  coord_flip()+
  facet_grid(HM~Attempt)

# convert to %
plate_otu_ps <- phyloseq::transform_sample_counts(plate_otu_ps, function(x) 100 * x/sum(x)) # TSS scale library sizes to 100%

# :: butyrogens -----------------------------------------------------------

# calculate butyrogens

plates_otu_ps_butyrogen = speedyseq::psmelt(plate_otu_ps) %>%
  subset(grepl("Lachnospiraceae", Family) | 
           grepl("Blautia", Genus) | 
           grepl("Roseburia", Genus) | 
           grepl("Eubacterium", Genus) | 
           grepl("Ruminococcus", Genus) | 
           grepl("Clostridium", Genus) | 
           grepl("Faecalibacterium", Genus)) %>%
  group_by(Sample) %>%
  mutate(But = sum(Abundance)) %>%
  dplyr::select(Plate, Attempt, Treatment, Replicate, HM, RS_Name, pH, Acetate, Propionate, Butyrate, But) %>% distinct() %>% data.frame() %>%
  arrange(HM)
# plate_otu_ps = NULL

# merge
plates_scfa_but = merge(plates_scfa, plates_otu_ps_butyrogen[,c("HM", "RS_Name", "Replicate", "Attempt", "But")], by=c("HM", "RS_Name", "Replicate", "Attempt"))
#plates_scfa = NULL
#plates_scfa_but = plates_scfa

unique(plates_scfa_but$Attempt)


# :: distance -------------------------------------------------------------

# load control data
control_data <- readRDS("~/Documents/PhD/Resistant_Starch/controls_nonibd-dc-cloud_n83.rds")

microbiomes_to_keep = unique(plates_scfa$microbiome)[!grepl(paste(c("PR2"), collapse="|"), unique(plates_scfa$microbiome))]

# loop through Microbiomes and calculate distance and score
plates_scfa_dis = do.call(rbind, lapply(microbiomes_to_keep, function(x){
  print(x)
  # subset to microbiome (Plate_Attempt)
  samples_to_keep = data.frame(phyloseq::sample_data(plate_otu_ps)) %>%
    subset(microbiome == x)
  # inelegent solution to subsetting samples:
  data.subset = phyloseq::phyloseq(phyloseq::otu_table(plate_otu_ps), phyloseq::tax_table(plate_otu_ps),
                                   phyloseq::sample_data(plate_otu_ps)[phyloseq::sample_data(plate_otu_ps)$Barcode %in% samples_to_keep$Barcode])
  phyloseq::sample_data(data.subset)$Sample = phyloseq::sample_data(data.subset)$Barcode
  # remove controls
  data.subset = phyloseq::subset_samples(data.subset, RS_Name %in% rs.names.pbs & !grepl("PR2", HM))
  # merge with controls
  microbiome_controls_ps <- phyloseq::merge_phyloseq(data.subset, control_data)
  # ensure 100%
  microbiome_controls_ps <- phyloseq::transform_sample_counts(microbiome_controls_ps, function(n) 100 * n/sum(n)) # TSS scale library sizes to 100%
  # Calculate Median Distance to Reference Cohort for all samples
  microbiome_controls_bray <- phyloseq::distance(microbiome_controls_ps, method="bray", type="samples") # Calculate Bray-Curtis dissimilarity matrix
  id.control <- as.character((phyloseq::sample_data(phyloseq::subset_samples(microbiome_controls_ps, Location=="Distal Colon" & Diagnosis=="Control"))$Sample)) # Vector of ONLY reference samples (comparisons will be made to this group)
  id.all <- as.character((phyloseq::sample_data(phyloseq::subset_samples(microbiome_controls_ps))$Sample)) # Vector of ALL samples
  microbiome_controls_bray_mat <- as.matrix(microbiome_controls_bray) # Matrix needed for median-distance code, but not distance-to-centroid
  microbiome_controls_distances <- do.call(rbind.data.frame,
                                           lapply(id.all, function (i) median(microbiome_controls_bray_mat[i, id.control])))
  names(microbiome_controls_distances)[1] <- "Distance" 
  microbiome_controls_distances$Sample = id.all
  
  # merge with sample data
  microbiome_controls_distances = merge(microbiome_controls_distances, 
                                        data.frame(phyloseq::sample_data(data.subset)), by="Sample")
  # calculate median
  microbiome_controls_distances = microbiome_controls_distances %>%
    group_by(HM, RS_Name) %>%
    mutate(med.dis = median(Distance)) %>%
    dplyr::select(HM, RS_Name, Distance, Replicate, med.dis, Plate, Attempt, microbiome) %>% distinct() %>% data.frame()
}))
plates_scfa_dis

# visualize distances
plates_scfa_dis$RS_Name = factor(plates_scfa_dis$RS_Name, levels=rs.names.pbs)
ggplot(plates_scfa_dis,
       aes(x=RS_Name, y=Distance))+
  geom_boxplot(width=0.2)+
  geom_point(shape=21, aes(fill=RS_Name), size=2.5)+
  scale_fill_manual(values=labelcolors$cols[c(10,1:9)])+
  theme_minimal()+theme(axis.text.x = element_text(angle=45, hjust=1))+
  facet_grid(Attempt~HM)+
  labs(x="", y="Median Distance to Controls")

# :: z-score --------------------------------------------------------------

# obtain median butyrogen, convert to Z-score, per HM and Attempt
plates_otu_ps_butyrogen_median = plates_otu_ps_butyrogen %>%
  subset(RS_Name %in% rs.names) %>%
  group_by(HM, Attempt, RS_Name) %>%
  mutate(med.but = median(But)) %>%
  dplyr::select(HM, RS_Name, Attempt, med.but) %>% distinct() %>%
  group_by(HM, Attempt) %>%
  mutate(zscore = scale(med.but)) %>% data.frame()

plates_otu_ps_butyrogen_median


# :: RS selection ---------------------------------------------------------

# merge in distance data
plates_scfa_but = merge(plates_scfa_but,
                        plates_scfa_dis[,c("HM", "RS_Name", "Replicate", "Attempt", "Distance", "med.dis")], by=c("HM", "RS_Name", "Replicate", "Attempt"))

plates_scfa_but_hm = subset(plates_scfa_but, HM != "0PR2") %>%
  subset(RS_Name %in% rs.names.pbs) %>%
  mutate(RS_Name = factor(RS_Name, levels=rs.names.pbs))%>%
  mutate(Attempt = as.factor(Attempt))


# First, using Hybrid Butyrogen+Distance
plates_scfa_but_hm_selections = do.call(rbind, lapply(microbiomes_to_keep, function(x){
  # subset to microbiome
  data.subset = subset(plates_scfa_but_hm, microbiome == x) %>%
    # calculate median
    subset(RS_Name != "PBS") %>% 
    group_by(RS_Name) %>%
    # calculate median
    mutate(med.but = median(But),
           med.dis = median(Distance),
           med.ph = median(pH)) %>%
    dplyr::select(HM, RS_Name, Attempt, med.but, med.dis, med.ph) %>% distinct() %>% 
    # calculate z-score
    group_by(HM, Attempt) %>%
    mutate(zscore = scale(med.but)) %>%
    # select Z > 1 or highest
    subset(zscore > 1 | zscore == max(zscore)) %>%
    subset(med.dis == max(med.dis))
  # 
  data.subset = data.subset[,c("HM", "RS_Name", "Attempt")] 
  colnames(data.subset)[2] = "z.selected"
  data.subset
}))%>% arrange(HM)
plates_scfa_but_hm_selections

# Second, using pH
plates_scfa_but_hm_selections_ph = do.call(rbind, lapply(microbiomes_to_keep, function(x){
  # subset to microbiome
  data.subset = subset(plates_scfa_but_hm, microbiome == x) %>%
    # calculate median
    subset(RS_Name != "PBS") %>% 
    group_by(RS_Name) %>%
    # calculate median
    mutate(med.but = median(But),
           med.dis = median(Distance),
           med.ph = median(pH)) %>%
    dplyr::select(HM, RS_Name, Attempt, med.but, med.dis, med.ph) %>% distinct() %>% 
    # calculate z-score
    group_by(HM, Attempt) %>%
    mutate(zscore = scale(med.ph)) %>%
    # select Z > 1 or highest
    subset(zscore == min(zscore)) 
  # 
  data.subset = data.subset[,c("HM", "RS_Name", "Attempt")] 
  colnames(data.subset)[2] = "ph.selected"
  data.subset
}))%>% arrange(HM)
# nearly all LDO, but much greater consistency
plates_scfa_but_hm_selections_ph

# :: finalize data frame --------------------------------------------------------------

# add RS selected
plates_otu_ps_butyrogen_median = merge(plates_otu_ps_butyrogen_median,
                                       plates_scfa_but_hm_selections, by=c("HM", "Attempt"))
# add pH selected
plates_scfa_but_ph_scores = merge(plates_scfa_but_hm,
                                  plates_scfa_but_hm_selections_ph, by=c("HM", "Attempt"))


# :: Add Sara-Krystal Replicate -------------------------------------------

# HM0899 was cultured twice, by Sara and Krystal (on the same day)
hm0899.tech1 = readRDS("~/Documents/PhD/For others/oars_validation/oars_validation_experiment_work/HM0899.00_scores_krystal.rds")
hm0899.tech2 = readRDS("~/Documents/PhD/For others/oars_validation/oars_validation_experiment_work/HM0899.00_scores_sara.rds")

hm0899.tech1$Attempt = "A"
hm0899.tech2$Attempt = "B"

hm0899.tech12 = rbind(hm0899.tech1,
                      hm0899.tech2) %>% data.frame()
# select RS
hm0899.tech12 = hm0899.tech12 %>%
  group_by(Attempt) %>%
  mutate(z.selected = ifelse(Z_score > 1, RS, NA))

hm0899.tech12$HM = "HM0899.00"
hm0899.tech12$RS_Name = hm0899.tech12$RS
hm0899.tech12$zscore = hm0899.tech12$Z_score

# merge with HM0819 data (unify replicate names)
similar.replicates.rapidaim = rbind(
  subset(plates_otu_ps_butyrogen_median, Attempt %in% c("4", "5"))[,c("HM", "Attempt", "RS_Name", "zscore", "z.selected")] %>%
    mutate(Attempt = ifelse(Attempt == "4", "A", "B")),
  hm0899.tech12[,c("HM", "Attempt", "RS_Name", "zscore", "z.selected")]
)

similar.replicates.rapidaim$Microbiome = paste("Microbiome", as.numeric(as.factor(similar.replicates.rapidaim$HM)))

similar.replicates.rapidaim.plot = ggplot(similar.replicates.rapidaim,
                                          aes(x=as.factor(Attempt), y=zscore))+
  geom_line(aes(group=paste(HM, RS_Name), color=RS_Name), alpha=0.5)+
  geom_point(shape=21, aes(fill=RS_Name), size=3)+
  geom_hline(yintercept=1, color="red", alpha=0.6)+
  scale_fill_manual(values=labelcolors$cols[c(1:9)])+
  scale_color_manual(values=labelcolors$cols[c(1:9)])+
  ggrepel::geom_text_repel(aes(label = ifelse(z.selected == as.character(RS_Name), as.character(RS_Name), NA)), size=2.5)+
  facet_wrap(~Microbiome,nrow=1)+
  theme_minimal()+theme(legend.position="none",
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="Culture Replicate", y="Butyrogen Z-Score")
similar.replicates.rapidaim.plot
# Microbiomes 1 and 2 were cultured 24 h apart
# Microbiome 3 was cultured at the same time by 2 different technicians

# :: Add Time Difference -------------------------------------------

# Compare Original culture to years later culture
hm0819.03.culture = readRDS("~/Documents/PhD/For others/oars_validation/oars_validation_experiment_work/HM0819.03_scores.rds")
hm0819.04.culture = readRDS("~/Documents/PhD/For others/oars_validation/oars_validation_experiment_work/HM0819.04_scores.rds")
hm0932.01.culture = readRDS("~/Documents/PhD/For others/oars_validation/oars_validation_experiment_work/HM0932.01_scores.rds")
hm0932.02.culture = readRDS("~/Documents/PhD/For others/oars_validation/oars_validation_experiment_work/HM0932.02_scores.rds")

hm0819.03.culture$Attempt = "0"
hm0819.03.culture$HM = "0819.03"

hm0819.04.culture$Attempt = "0"
hm0819.04.culture$HM = "0819.04"

hm0932.01.culture$Attempt = "0"
hm0932.01.culture$HM = "0932.01"

hm0932.02.culture$Attempt = "0"
hm0932.02.culture$HM = "0932.02"


hm.original.cultures = rbind(hm0819.03.culture,
                             hm0819.04.culture,
                             hm0932.01.culture,
                             hm0932.02.culture) %>% data.frame()
# select RS
hm.original.cultures = hm.original.cultures %>%
  group_by(Attempt, HM) %>%
  mutate(z.selected = ifelse(Z_score > 1, RS, NA)) %>% data.frame()

hm.original.cultures$RS_Name = hm.original.cultures$RS
hm.original.cultures$zscore = hm.original.cultures$Z_score

# merge with original data: plates_otu_ps_butyrogen_median
hm.original.cultures.rapidaim = rbind(
  subset(plates_otu_ps_butyrogen_median, Attempt %in% c("1", "2", "3"))[,c("HM", "Attempt", "RS_Name", "zscore", "z.selected")],
  hm.original.cultures[,c("HM", "Attempt", "RS_Name", "zscore", "z.selected")]
)

hm.original.cultures.rapidaim$Microbiome = paste("Microbiome", as.numeric(as.factor(hm.original.cultures.rapidaim$HM)))

# Add date of culture
hm.original.cultures.rapidaim.dates = read.csv("~/Documents/PhD/For others/oars_validation/oars_validation_experiment_work/oars_validation_rapidaim_dates.csv")
hm.original.cultures.rapidaim.dates$HM = paste("0", hm.original.cultures.rapidaim.dates$HM, sep="")
hm.original.cultures.rapidaim = merge(hm.original.cultures.rapidaim,
                                      hm.original.cultures.rapidaim.dates, by=c("HM", "Attempt"))

hm.original.cultures.rapidaim.plot = ggplot(hm.original.cultures.rapidaim %>%
                                              mutate(Microbiome = gsub("Microbiome 3", "Microbiome 4", gsub("Microbiome 4", "Microbiome 5", Microbiome))) %>%
                                              subset(Attempt %in% c(0,1)), 
                                            aes(x=as.character(as.Date(Date, "%d-%b-%y")), y=zscore))+
  geom_line(aes(group=paste(HM, RS_Name), color=RS_Name), alpha=0.5)+
  geom_point(shape=21, aes(fill=RS_Name), size=3)+
  geom_hline(yintercept=1, color="red", alpha=0.6)+
  scale_fill_manual(values=labelcolors$cols[c(1:9)])+
  scale_color_manual(values=labelcolors$cols[c(1:9)])+
  ggrepel::geom_text_repel(aes(label = ifelse(z.selected == as.character(RS_Name), as.character(RS_Name), NA)), size=2.5)+
  facet_wrap(~Microbiome,nrow=1, scales="free_x")+
  theme_minimal()+theme(legend.position="none",
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="Culture Date", y="Butyrogen Z-Score")
hm.original.cultures.rapidaim.plot


# :: Correlate pH~SCFA --------------------------------------------------------
plates_scfa_but_hm$SCFA = plates_scfa_but_hm$Acetate +plates_scfa_but_hm$Propionate+plates_scfa_but_hm$Butyrate

oars.val.cor.plot1 = ggplot(subset(plates_scfa_but_hm, Attempt == 1) %>%
                              mutate(Microbiome = paste("Microbiome", as.numeric(as.factor(HM)), sep=" ")),
                            aes(x=pH, y=SCFA/1000))+
  geom_point(aes(fill=RS_Name), shape=21)+
  geom_smooth(method="lm", color="black")+
  ggpubr::stat_cor(method="spearman", size=3)+
  scale_fill_manual(values=labelcolors$cols[c(10,1:9)])+
  facet_wrap(~Microbiome, scales="free", nrow=1)+
  theme_classic()+theme(legend.position="none")+
  labs(x="pH", y="Total SCFA (mM)")
oars.val.cor.plot1

# :: Correlate pH~Butyrate --------------------------------------------------------

oars.val.cor.plot2 = ggplot(subset(plates_scfa_but_hm, Attempt == 1) %>%
                              mutate(Microbiome = paste("Microbiome", as.numeric(as.factor(HM)), sep=" ")),
                            aes(x=pH, y=Butyrate/1000))+
  geom_point(aes(fill=RS_Name), shape=21)+
  geom_smooth(method="lm", color="black")+
  ggpubr::stat_cor(method="spearman", size=3)+
  scale_fill_manual(values=labelcolors$cols[c(10,1:9)])+
  facet_wrap(~Microbiome, scales="free", nrow=1)+
  theme_classic()+theme(legend.position="none")+
  labs(x="pH", y="Butyrate (mM)")
oars.val.cor.plot2

# :: Correlate Butyrogens~SCFA --------------------------------------------------------

oars.val.cor.plot3 = ggplot(subset(plates_scfa_but_hm, Attempt == 1) %>%
                              mutate(Microbiome = paste("Microbiome", as.numeric(as.factor(HM)), sep=" ")),
                            aes(x=But, y=SCFA/1000))+
  geom_point(aes(fill=RS_Name), shape=21)+
  geom_smooth(method="lm", color="black")+
  ggpubr::stat_cor(method="spearman", size=3)+
  scale_fill_manual(values=labelcolors$cols[c(10,1:9)])+
  facet_wrap(~Microbiome, scales="free", nrow=1)+
  theme_classic()+theme(legend.position="none")+
  labs(x="Butyrogens %", y="Total SCFA (mM)")
oars.val.cor.plot3

# :: Correlate Butyrogens~Butyrate --------------------------------------------------------

oars.val.cor.plot4 = ggplot(subset(plates_scfa_but_hm, Attempt == 1) %>%
                              mutate(Microbiome = paste("Microbiome", as.numeric(as.factor(HM)), sep=" ")),
                            aes(x=But, y=Butyrate/1000))+
  geom_point(aes(fill=RS_Name), shape=21)+
  geom_smooth(method="lm", color="black")+
  ggpubr::stat_cor(method="spearman", size=3)+
  scale_fill_manual(values=labelcolors$cols[c(10,1:9)])+
  facet_wrap(~Microbiome, scales="free", nrow=1)+
  theme_classic()+theme(legend.position="none")+
  labs(x="Butyrogens %", y="Butyrate (mM)")
oars.val.cor.plot4


# :: pH-dependence on SCFA correlation ------------------------------------

# 2025_07_29  Correlations seem to depend on pH; what if we cut off pH values below 5.5
plates_scfa_but_hm_threshold = plates_scfa_but_hm %>%
  subset(HM %in% c("0819.03", "0819.04", "0932.01", "0932.02")) %>%
  subset(RS_Name %in% rs.names.pbs) %>%
  subset(Attempt == 1)

plates_scfa_but_hm_threshold$Microbiome = paste("Microbiome", as.numeric(as.factor(plates_scfa_but_hm_threshold$HM)), sep=" ")


# see how correlations change with threshold
oars_scfa_ph_threshold = do.call(rbind, lapply(seq(from=4.85, to=7, by=0.05), function(ph){
  do.call(rbind, lapply(unique(plates_scfa_but_hm_threshold$Microbiome), function(hm){
    print(paste(hm, ph))
    data.subset.1 = subset(plates_scfa_but_hm_threshold, pH > ph & Microbiome == hm)
    if(nrow(data.subset.1)<3){
      data.right = data.frame(Microbiome = hm,
                              threshold = ph,
                              cor = NA,
                              pval = NA,
                              side = "right")
    }else{
      cor.results.1 = cor.test(data.subset.1$pH,
                               data.subset.1$Butyrate, method="spearman")
      data.right = data.frame(Microbiome = hm,
                              threshold = ph,
                              cor = cor.results.1$estimate,
                              pval = cor.results.1$p.value,
                              side = "right")
    }
    # now do left side
    data.subset.2 = subset(plates_scfa_but_hm_threshold, pH <= ph & Microbiome == hm)
    if(nrow(data.subset.2)<3){
      data.left = data.frame(Microbiome = hm,
                             threshold = ph,
                             cor = NA,
                             pval = NA,
                             side = "left")
    }else{
      cor.results.2 = cor.test(data.subset.2$pH,
                               data.subset.2$Butyrate, method="spearman")
      data.left = data.frame(Microbiome = hm,
                             threshold = ph,
                             cor = cor.results.2$estimate,
                             pval = cor.results.2$p.value,
                             side = "left")
    }
    # merge
    rbind(data.right, data.left) %>% data.frame()
  }))}))

# add minima
plates_scfa_but_hm_threshold = do.call(rbind, lapply(unique(oars_scfa_ph_threshold$Microbiome), function(hm){
  print(hm)
  data.subset = subset(oars_scfa_ph_threshold, Microbiome == hm & !is.na(cor))
  other.data = subset(plates_scfa_but_hm_threshold, Microbiome == hm)
  # build gam
  gam.model.left = mgcv::gam(cor ~ s(threshold), data=subset(data.subset, side == "left"))
  gam.model.right = mgcv::gam(cor ~ s(threshold), data=subset(data.subset, side == "right"))
  
  # narrow to 5 to 5.6
  gam.fitted.left = data.frame(fitted = gam.model.left$fitted.values,
                               threshold = subset(data.subset, side == "left")$threshold) %>%
    subset(threshold > 5 & threshold < 6) %>% mutate(side = "left")
  gam.fitted.right = data.frame(fitted = gam.model.right$fitted.values,
                                threshold = subset(data.subset, side == "right")$threshold) %>%
    subset(threshold > 5 & threshold < 6) %>% mutate(side = "right")
  # max difference
  gam.differences = merge(gam.fitted.left,gam.fitted.right, by="threshold") %>%
    mutate(delta = fitted.x - fitted.y)
  cgpt = subset(gam.differences, abs(delta) == max(abs(delta)))$threshold %>% unique()
  other.data$selected = subset(data.subset, threshold == cgpt)$threshold %>% unique()
  other.data
}))

oars_scfa_ph_threshold = merge(oars_scfa_ph_threshold,
                               plates_scfa_but_hm_threshold[,c("Microbiome", "selected")]%>% distinct(),
                               by="Microbiome")

oars_scfa_ph_threshold.plot = ggplot(oars_scfa_ph_threshold,
                                     aes(x=threshold, y=cor))+
  geom_smooth(aes(group = side, color=side))+
  geom_vline(aes(xintercept = selected), color="red", linetype=2)+
  scale_fill_gradient2(low="blue", high="red")+
  facet_wrap(~Microbiome, nrow=1)+
  xlim(4.85, 7.2)+
  theme_classic()+theme(legend.position="none",
                        strip.text=element_text(size=10))+
  labs(x="pH Threshold", y="Correlation")
oars_scfa_ph_threshold.plot

plates_scfa_but_hm_threshold$threshold = ifelse(plates_scfa_but_hm_threshold$pH < plates_scfa_but_hm_threshold$selected, "left", "right")
# see correlation with 5.5 marker
plates_scfa_but_hm_5.5.plot = ggplot(plates_scfa_but_hm_threshold,
                                     aes(x=pH, y=Butyrate/1000))+
  geom_point(aes(fill=RS_Name), shape=21)+
  geom_vline(aes(xintercept = selected), color="red", linetype=2)+
  geom_smooth(aes(group=threshold, color=threshold), method="lm")+
  ggpubr::stat_cor(aes(group = threshold, color=threshold),
                   method="spearman", size=3,
                   label.x.npc = 0.5,
                   label.y.npc = 0.95)+
  xlim(4.85, 7.2)+
  scale_fill_manual(values=labelcolors$cols[c(10,1:9)])+
  facet_wrap(~Microbiome,  nrow=1)+
  theme_classic()+theme(legend.position="none",
                        strip.text=element_text(size=10))+
  labs(x="pH", y="Butyrate (mM)")
plates_scfa_but_hm_5.5.plot


# see how correlations change 
plates_scfa_but_hm_5.5_butyrogen.plot = ggplot(plates_scfa_but_hm_threshold,
                                               aes(x=But, y=Butyrate/1000))+
  geom_point(aes(fill=RS_Name), shape=21)+
  geom_smooth(aes(group=threshold, color=threshold), method="lm")+
  ggpubr::stat_cor(aes(group = threshold, color=threshold),
                   method="spearman", size=3,
                   label.x.npc = 0.1,
                   label.y.npc = 0.95)+
  scale_fill_manual(values=labelcolors$cols[c(10,1:9)])+
  facet_wrap(~Microbiome, nrow=1)+
  theme_classic()+theme(legend.position="none",
                        strip.text=element_text(size=10))+
  labs(x="Median Butyrogens %", y="Butyrate (mM)")
plates_scfa_but_hm_5.5_butyrogen.plot


# see how correlations change 
plates_scfa_but_hm_5.5_acebut.plot = ggplot(plates_scfa_but_hm_threshold,
                                            aes(x=pH, y=log2(Acetate/Butyrate)))+
  geom_point(aes(fill=RS_Name), shape=21)+
  geom_vline(aes(xintercept = selected), color="red", linetype=2)+
  geom_smooth(aes(group=threshold, color=threshold), method="lm")+
  ggpubr::stat_cor(aes(group = threshold, color=threshold),
                   method="spearman", size=3,
                   label.x.npc = 0.5,
                   label.y.npc = 0.15)+
  xlim(4.85, 7.2)+
  scale_fill_manual(values=labelcolors$cols[c(10,1:9)])+
  facet_wrap(~Microbiome,  nrow=1)+
  theme_classic()+theme(legend.position="none",
                        strip.text=element_text(size=10))+
  labs(x="pH", y="Log2(Acetate/Butyrate)")
plates_scfa_but_hm_5.5_acebut.plot


oars_scfa_ph_threshold.plot/
  plates_scfa_but_hm_5.5.plot/
  plates_scfa_but_hm_5.5_butyrogen.plot/
  plates_scfa_but_hm_5.5_acebut.plot

# > plots -----------------------------------------------------------------

(oars.rapidaim.scores.pc.selection.plot+ 
  oars.rapidaim.scores.pc.frequencies.plot+
  patchwork::plot_layout(widths=c(4,1))) %>%
  ggsave(filename="./oars_plots_deidentified/oars_1_rs_selections_DI.pdf",
         width=15, height=7, device = cairo_pdf)

(oars.rapidaim.scores.pc.ph.selection.plot+ 
    oars.rapidaim.scores.pc.ph.frequencies.plot+
    patchwork::plot_layout(widths=c(4,1))) %>%
  ggsave(filename="./oars_plots/oars_rs_ph_selections.pdf",
         width=15, height=7, device = cairo_pdf)

(oars.stool.stability.pcoa.plot+
    oars.rapidaim.stability.plot) %>%
  ggsave(filename="./oars_plots/oars_1b_rs_selection_stability.pdf",
         width=7, height=3, device = cairo_pdf)

oars.sensitivity.results.rates.plot%>%
  ggsave(filename="./oars_plots/oars_1b_rs_selection_sensitivity.pdf",
         width=4, height=3, device = cairo_pdf)

# pH-Butyrogen correlations
((oars.rapidaim.ph.butyrogens.plot.819+oars.rapidaim.ph.but.cor.plot)+
  patchwork::plot_layout(nrow=1, widths=c(2,1))) %>%
  ggsave(filename="./oars_plots/oars_1b_butyrogen_ph_cor.pdf",
         width=10, height=4, device = cairo_pdf)


(similar.replicates.rapidaim.plot/hm.original.cultures.rapidaim.plot) %>%
  ggsave(filename="./oars_plots/oars_1b_rapidaim_validation_both.pdf",
         width=7, height=7, device = cairo_pdf)

(oars.val.cor.plot1+
    oars.val.cor.plot2+
    oars.val.cor.plot3+
    oars.val.cor.plot4+
    patchwork::plot_layout(nrow=4))  %>%
  ggsave(filename="./oars_plots/oars_1b_rapidaim_validation_correlations.pdf",
         width=10, height=8, device = cairo_pdf)

(oars_scfa_ph_threshold.plot/
    plates_scfa_but_hm_5.5.plot/
    plates_scfa_but_hm_5.5_butyrogen.plot/
    plates_scfa_but_hm_5.5_acebut.plot) %>%
  ggsave(filename="./oars_plots/oars_1b_rapidaim_ph_thresholds.pdf",
         width=12, height=10, device = cairo_pdf)

(plates_scfa_but_hm_5.5.plot) %>%
  ggsave(filename="../For others/2025_08_08_gbm_figure_1_PD.pdf",
         width=12, height=3, device = cairo_pdf)


# >>> Phasco/Dialister -------------------------------------------------------------------

# Goal: Based on LSARP-CD, we need to ask

# For the RS that was selected, does Phasco/Dialister ratio (change to PBS)
# predict clinical response? (in this case, fecal calprotectin)

# Answer: No

# load OTU table
oars.rapidaim.phyloseq = readRDS("./2025_06_09_oars_phyloseq.Rds")

# transform to 100%
oars.rapidaim.phyloseq <- phyloseq::transform_sample_counts(oars.rapidaim.phyloseq, function(x) 100 * x/sum(x)) # TSS scale library sizes to 100%

# convert to data frame
oars.rapidaim.phyloseq.df = speedyseq::psmelt(oars.rapidaim.phyloseq)
# LCA
oars.rapidaim.phyloseq.df = data.frame(oars.rapidaim.phyloseq.df) %>%
  mutate(LCA = 
           ifelse(!is.na(Species)&Species!="", paste(as.character(Genus), as.character(Species), sep="_"), 
                  ifelse(!is.na(Genus)&Genus!="", paste(Genus), 
                         ifelse(!is.na(Family)&Family!="", paste(Family),
                                ifelse(!is.na(Order)&Order!="", paste(Order),
                                       ifelse(!is.na(Class)&Class!="", paste(Class),
                                              ifelse(!is.na(Phylum)&Phylum!="", paste(Phylum),
                                                     ifelse(is.na(Phylum), "undefined", paste(Phylum)))))))))
# collapse to median
oars.rapidaim.phyloseq.df = oars.rapidaim.phyloseq.df %>% dplyr::select(Sample, RS_Name, HM, Replicate, pH, timing, LCA, Abundance)

oars.rapidaim.phyloseq.med = oars.rapidaim.phyloseq.df %>%
  group_by(HM, RS_Name, LCA, Sample) %>%
  # sum up LCA
  mutate(sum.abun = sum(Abundance)) %>%
  # median
  group_by(HM, RS_Name, LCA) %>%
  mutate(med.abun = median(sum.abun)) %>%
  mutate(med.ph = median(pH))%>%
  dplyr::select(HM, RS_Name, med.abun, LCA, med.ph, timing) %>% distinct()
range(oars.rapidaim.phyloseq.med$med.abun)

# logfc to PBS
pseudo = min(oars.rapidaim.phyloseq.med[oars.rapidaim.phyloseq.med$med.abun!=0,]$med.abun)/2

oars.rapidaim.phyloseq.lfc = oars.rapidaim.phyloseq.med %>%
  subset(RS_Name %in% rs.names.pbs) %>%
  group_by(HM, LCA) %>%
  mutate(lfc = log2(med.abun+pseudo) - log2(med.abun[RS_Name == "PBS"]+pseudo)) %>%
  mutate(delta_ph = med.ph - med.ph[RS_Name == "PBS"]) %>%
  mutate(HM = paste("HM", substr(HM, 1, 4), sep=""))%>%
  mutate(sample = paste(HM, "_", timing, sep=""))
 
# add clinical data
oars.rapidaim.responses = subset(metadata.oars.stool.double, reltiming == "pre")[,c("HM", "timing", "response", "rs.selected", "lfc.fcal")] %>% distinct()
oars.rapidaim.responses$fcal.response = ifelse(oars.rapidaim.responses$lfc.fcal > 0, "Increase", "Decrease")

oars.rapidaim.phyloseq.lfc = merge(oars.rapidaim.phyloseq.lfc,
                                   oars.rapidaim.responses, by=c("HM", "timing"))

# subset to selected RS
oars.rapidaim.phyloseq.lfc.rs = subset(oars.rapidaim.phyloseq.lfc, RS_Name == rs.selected)

# subset to Phasco/Dialister
oars.rapidaim.phyloseq.lfc.rs = subset(oars.rapidaim.phyloseq.lfc.rs, LCA %in% c("Phascolarctobacterium", "Dialister")) %>% distinct()
oars.rapidaim.phyloseq.lfc.rs = oars.rapidaim.phyloseq.lfc.rs %>%
  group_by(HM, timing) %>%
  mutate(pd = lfc[LCA == "Phascolarctobacterium"] - lfc[LCA == "Dialister"]) %>%
  dplyr::select(-LCA, -lfc, -med.abun) %>% distinct()

ggplot(oars.rapidaim.phyloseq.lfc.rs %>%
         mutate(response = ifelse(response == "high", "Strong", "Weak")),
       aes(x=response, y=pd))+
  geom_boxplot()+
  geom_point(aes(fill=response), shape=21, size=3)+
  theme_classic()+theme(legend.position="none")+
  labs(x="Response", y="Phascolarctobacterium / Dialister")
# no association with fermentation response

ggplot(oars.rapidaim.phyloseq.lfc.rs,
       aes(x=fcal.response, y=pd))+
  geom_boxplot()+
  geom_point(aes(fill=fcal.response), shape=21, size=3)+
  theme_classic()+theme(legend.position="none")+
  labs(x="Fecal Calprotectin", y="Phascolarctobacterium / Dialister")
# no association with fecal calprotectin response

ggplot(oars.rapidaim.phyloseq.lfc.rs,
       aes(x=lfc.fcal, y=pd))+
  geom_point(aes(fill=fcal.response), shape=21, size=3)+
  geom_smooth(method="lm", color="black")+
  theme_classic()+theme(legend.position="none")+
  labs(x="Log2FC Fecal Calprotectin", y="Phascolarctobacterium / Dialister")
# no association

# Does it need to be fermented
ggplot(oars.rapidaim.phyloseq.lfc.rs,
       aes(x=delta_ph, y=pd))+
  geom_point(aes(fill=delta_ph), shape=21, size=3)+
  geom_smooth(method="lm", color="black")+
  ggpubr::stat_cor(method="spearman")+
  theme_classic()+theme(legend.position="none")+
  labs(x="Δ pH", y="Phascolarctobacterium / Dialister")
# no association

lmerTest::lmer(pd ~ fcal.response + delta_ph + rs.selected + (1|HM),
               oars.rapidaim.phyloseq.lfc.rs) %>%
  summary()

# no associations

# Different patient populations
# LSARP-CD used clinical response, not fecal calprotectin

# :: ----------------------------------------------------------------------


# \\ defunct --------------------------------------------------------------


# \\ vital butyrogens --------------------------------------------------

# goal: assess RS selections when we use Vital's butyrogen algorithm

# here is the ASV phyloseq object
oars.rapidaim.ps.final

# now process through PICRUSt2-Vital
# export
fasta_seqs <- Biostrings::DNAStringSet(rownames(phyloseq::tax_table(oars.rapidaim.ps.final)))
names(fasta_seqs) = fasta_seqs
Biostrings::writeXStringSet(fasta_seqs, filepath = "oars_picrust2/predict_SCFA_producers/oars.rapidaim.picrust2.seqtable.fasta")
write.table(t(phyloseq::otu_table(oars.rapidaim.ps.final)), "oars_picrust2/oars.rapidaim.picrust2.abuntable.tsv", 
            sep="\t", quote = F, col.names = NA)
# in terminal
cd ~/Documents/PhD/git_oars_archfolder/oars_picrust2/predict_SCFA_producers
conda activate oars_picrust2
place_seqs.py -s oars.rapidaim.picrust2.seqtable.fasta -o placed_seqs.tre -p 1 --intermediate placement_working --ref_dir SCFA
hsp.py -t placed_seqs.tre --observed_trait_table SCFA/SCFA_pathwaydata.txt -o SCFA_rapidaim_predicted.tsv -p 1 -m emp_prob -n

oars.rapidaim.vital.butyrogens = read.csv("./oars_picrust2/predict_SCFA_producers/SCFA_rapidaim_predicted.tsv", sep="\t")
# subset to acetylcoa+but or buk alone
oars.rapidaim.vital.butyrogens.acetylcoa = subset(oars.rapidaim.vital.butyrogens, acetylcoa == 1)
oars.rapidaim.vital.butyrogens.acetylcoa = subset(oars.rapidaim.vital.butyrogens.acetylcoa, but == 1 | buk == 1)
# these will be "butyrogens"
oars.rapidaim.vital.butyrogens.acetylcoa

# return to phyloseq object
oars.rapidaim.phyloseq.butyrogens <- phyloseq::otu_table(oars.rapidaim.ps.final)[,oars.rapidaim.vital.butyrogens.acetylcoa$sequence]
oars.rapidaim.phyloseq.butyrogens <- phyloseq::merge_phyloseq(oars.rapidaim.phyloseq.butyrogens, 
                                                              phyloseq::tax_table(oars.rapidaim.ps.final), phyloseq::sample_data(oars.rapidaim.ps.final))
# format to dataframe
oars.rapidaim.phyloseq.butyrogens = speedyseq::psmelt(oars.rapidaim.phyloseq.butyrogens)
oars.rapidaim.phyloseq.butyrogens$Abundance = oars.rapidaim.phyloseq.butyrogens$Abundance / 50000
# visualize
oars.rapidaim.phyloseq.butyrogens %>%
  group_by(LCA) %>%
  mutate(sum.abun = sum(Abundance)) %>% dplyr::select(LCA, sum.abun) %>% distinct() %>% data.frame() %>%
  slice_max(n=20, sum.abun) %>%
  ggplot(aes(x=reorder(LCA, sum.abun), y=sum.abun)) +
  coord_flip()+scale_y_log10()+
  geom_point(shape=21, aes(fill=sum.abun), size=2.5)+
  scale_fill_gradient2(low="blue", high="red")+
  theme_minimal()+theme(legend.position="none")+labs(x="", y="Total Abundance (log10)")
# top are Agathobacter, Faecalibacterium, Gemmiger, Anaerostipes, Anaerobutyricum, Roseburia

# calculate median, and Z-score
oars.rapidaim.phyloseq.butyrogens = oars.rapidaim.phyloseq.butyrogens[,c("Sample", "Abundance", "Replicate", "HM", "RS_Name", "round", "pH")] %>%
  group_by(Sample) %>%
  mutate(but = sum(Abundance)) %>%
  dplyr::select(HM, RS_Name, but, round, pH) %>% distinct() %>%
  # subset to RS
  subset(RS_Name %in% rs.names.pbs) %>%
  # calculate median
  group_by(HM, RS_Name) %>%
  mutate(med.but = median(but)) %>%
  mutate(med.ph = median(pH)) %>%
  dplyr::select(HM, RS_Name, round, med.but, med.ph) %>% distinct()

oars.rapidaim.phyloseq.butyrogens$sample = oars.rapidaim.phyloseq.butyrogens$HM

# correlate with original butyrogens
oars.rapidaim.scores.pc.full = readRDS("./2025_06_09_oars_scores.Rds")

oars.rapidaim.phyloseq.butyrogens = merge(oars.rapidaim.phyloseq.butyrogens,
                                          oars.rapidaim.scores.pc.full[c("sample", "RS_Name", "med.but", "Z_score",  "selected", "med.ph")], by=c("sample", "RS_Name"))
# rename
colnames(oars.rapidaim.phyloseq.butyrogens) = c("sample", "RS_Name", "HM", "round", "but.vital", "ph.x", "but.classic", "z.y", "selected.y", "ph.y")

# ensure pH is the same
cor.test(oars.rapidaim.phyloseq.butyrogens$ph.x,
         oars.rapidaim.phyloseq.butyrogens$ph.y) # good
# correlate but
cor.test(oars.rapidaim.phyloseq.butyrogens$but.classic,
         oars.rapidaim.phyloseq.butyrogens$but.vital, method="spearman") # r = 0.79, sig
# plot
oars.rapidaim.phyloseq.butyrogens$RS_Name = factor(oars.rapidaim.phyloseq.butyrogens$RS_Name, levels=rs.names.pbs)


# :: correlate with classic butyrogens ------------------------------------

oars.rapidaim.phyloseq.butyrogens.plot = ggplot(oars.rapidaim.phyloseq.butyrogens,
                                                aes(x=but.classic, y=but.vital))+
  scale_x_log10()+
  scale_y_log10()+
  geom_point(shape=21, aes(fill=RS_Name), size=2.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_fill_manual(values=labelcolors$cols[c(10,1:9)])+
  ggpubr::stat_cor(method="spearman")+
  theme_minimal()+
  labs(x="Classic Butyrogens %",
       y="Vital Butyrogens %")
oars.rapidaim.phyloseq.butyrogens.plot

# :: compare Z-scores ------------------------------------

oars.rapidaim.phyloseq.butyrogens.z = oars.rapidaim.phyloseq.butyrogens %>%
  subset(RS_Name != "PBS") %>%
  group_by(sample) %>%
  mutate(Z_score.vital = scale(but.vital)) %>%
  mutate(Z_score.classic = scale(but.classic))
# correlate this!
cor.test(oars.rapidaim.phyloseq.butyrogens.z$Z_score.classic,
         oars.rapidaim.phyloseq.butyrogens.z$Z_score.vital, method="spearman") # r = 0.63, sig

oars.rapidaim.phyloseq.butyrogens.z.plot = ggplot(oars.rapidaim.phyloseq.butyrogens.z,
                                                  aes(x=Z_score.classic, y=Z_score.vital))+
  geom_point(shape=21, aes(fill=RS_Name), size=2.5)+
  geom_smooth(method="lm", se=F, color="black")+
  scale_fill_manual(values=labelcolors$cols[c(1:9)])+
  ggpubr::stat_cor(method="spearman")+
  theme_minimal()+
  labs(x="Classic Butyrogens Z-score",
       y="Vital Butyrogens Z-score")
oars.rapidaim.phyloseq.butyrogens.z.plot

# stats
oars.rapidaim.phyloseq.butyrogens.selection = oars.rapidaim.phyloseq.butyrogens.z %>%
  group_by(sample) %>%
  mutate(rs.vital = ifelse(Z_score.classic == max(Z_score.classic), as.character(RS_Name), NA)) %>%
  mutate(rs.classic = ifelse(Z_score.vital == max(Z_score.vital), as.character(RS_Name), NA)) %>%
  # calculate rate that same selection occured
  group_by(sample) %>%
  mutate(same = ifelse(!is.na(rs.vital) == !is.na(rs.classic), "same","different")) %>%
  mutate(same = ifelse(is.na(rs.vital) & is.na(rs.classic), NA, same)) %>%
  subset(!is.na(same)) %>%
  group_by(sample) %>%
  dplyr::select(sample, same) %>% distinct()

oars.rapidaim.phyloseq.butyrogens.selection[,c("same")] %>% table() 
# 22 remain the same; 24 are different



# :: F. prausnitzii in MLI ------------------------------------------------

# need to load Jenn's old UC data
mli_otu_data <- phyloseq::import_biom("./mli_otu_data/otu_otutable_200224.biom", parseFunction = phyloseq::parse_taxonomy_greengenes)
mli_otu_map <- read.csv("./mli_otu_data/glom_map.txt", sep="\t", row.names=1) 
mli_otu_map = phyloseq::sample_data(mli_otu_map, errorIfNULL=T)
mli_otu_map$Sample <- rownames(mli_otu_map)
mli_otu_data <- phyloseq::merge_phyloseq(mli_otu_data, mli_otu_map) #Create phyloseq object (without tree)
# add tree
mli_otu_data = phyloseq::merge_phyloseq(mli_otu_data, phyloseq::read_tree("../16s_databases/97_otus_gg13v5.tree")) 
# good

# among F. prausnitzii, see which are associated with disease
mli_otu_data_fp = speedyseq::psmelt(mli_otu_data)
mli_otu_data_fp = subset(mli_otu_data_fp, Genus == "Faecalibacterium")
# subset to specific OTUs?
#mli_otu_data_fp_otu = subset(mli_otu_data_fp, OTU %in% tree_meta$OTU)
mli_otu_data_fp_otu_mat = reshape2::acast(mli_otu_data_fp,
                                          OTU ~ Sample, value.var="Abundance")
mli_otu_data_fp_otu_mat = mli_otu_data_fp_otu_mat / 150000
# maaslin2

mli_otu_data_fp_otu_maaslin = Maaslin2::Maaslin2(input_data = mli_otu_data_fp_otu_mat,
                                                 input_metadata = subset(mli_otu_map,
                                                                         DiagnosisInflammation %in% c("Control Not Inflamed", "UC Inflamed") & 
                                                                           Location %in% c("Distal Colon")) %>% data.frame(),
                                                 output = "~/Downloads",
                                                 fixed_effects = c("DiagnosisInflammation"),  # Example fixed effects
                                                 #random_effects = c("study_id"),       # Example random effects
                                                 normalization = "TSS",                       # Total Sum Scaling normalization
                                                 transform = "LOG",                           # Log transformation
                                                 analysis_method = "LM",                      # Linear model
                                                 plot_scatter = FALSE,                        # Disable scatterplot generation
                                                 plot_heatmap = FALSE,                         # Keep heatmap generation (optional)
                                                 max_significance = 0.05,                     # Significance threshold for q-values
                                                 standardize = TRUE                           # Disable standardization (optional)
)
mli_otu_data_fp_otu_maaslin = mli_otu_data_fp_otu_maaslin$results %>% data.frame() %>%
  subset(metadata != "Location") %>%
  mutate(padj = p.adjust(pval, method="BH"))
mli_otu_data_fp_otu_maaslin = merge(mli_otu_data_fp_otu_maaslin %>% mutate(OTU = gsub("X", "", feature)),
                                    oars.rapidaim.otu.names.fprau[,c("OTU", "name")], by="OTU")
mli_otu_data_fp_otu_maaslin %>% arrange(padj)
# none are associated with inflammation in UC MLI



# :: F. prausnitzii ConCor ------------------------------------------------

# Use lmtree to find conditions that split F.p. correlation w pH

oars.fp.concor = oars.rapidaim.phyloseq.butyrogens.mat %>% data.frame()
oars.fp.concor = oars.fp.concor[,colnames(oars.fp.concor) %in%
                                  c(subset(oars.rapidaim.ph.positive.maaslin.df, qval < 0.05)$feature,
                                    subset(oars.rapidaim.ph.negative.maaslin.df, qval < 0.05)$feature)]
oars.fp.concor$Barcode = rownames(oars.fp.concor)
oars.fp.concor = merge(oars.fp.concor,
                       subset(oars.rapidaim.phyloseq.butyrogens.map, 
                              timing %in% c("0M", "3M", "6M"))[,c("Barcode", "HM", "delta.ph")],
                       by="Barcode")
rownames(oars.fp.concor) = oars.fp.concor$Barcode
oars.fp.concor$Barcode = NULL
oars.fp.concor$HM = as.factor(oars.fp.concor$HM)

# 177 (no change in sign)
# 173 (changes in sign)
# 154 (changes in sign)
# 162 (no change in sign)


# build lmertree
taxa = "Faecalibacterium_prausnitzii.154"
partition = colnames(oars.fp.concor)[!colnames(oars.fp.concor) %in% c("HM", "delta.ph",taxa)]
formula = as.formula(paste("delta.ph ~ x_norm | HM |", paste(c(partition), collapse="+")))
oars.fp.concor.subset = oars.fp.concor

# transform: target to log; predictors to presence/absence
oars.fp.concor.subset[colnames(oars.fp.concor.subset) == taxa] = log10(oars.fp.concor.subset[colnames(oars.fp.concor.subset) == taxa] + 0.000005)
colnames(oars.fp.concor.subset)[colnames(oars.fp.concor.subset) == taxa] = "taxa"
# scale
oars.fp.concor.subset = oars.fp.concor.subset %>%
  group_by(HM) %>%
  mutate(x_norm = taxa - mean(taxa)) %>%
  ungroup() %>%
  dplyr::select(-taxa)

oars.fp.concor.subset[colnames(oars.fp.concor.subset) %in% partition] = ifelse(oars.fp.concor.subset[colnames(oars.fp.concor.subset) %in% partition]>0, 1, 0)

set.seed(25)
oars.fp.concor.tree = glmertree::lmertree(formula,
                                          oars.fp.concor.subset, minsize = 200)

plot(oars.fp.concor.tree)
# extract nodes
oars.fp.concor.subset$node <- as.factor(predict(oars.fp.concor.tree, type = "node"))

# Plot with ggplot2 + geom_smooth per node
ggplot(oars.fp.concor.subset, 
       aes(x = x_norm, y = delta.ph)) +
  geom_point(shape=21, aes(fill=node)) +
  geom_smooth(method = "lm", se = FALSE, color="black") +
  ggpubr::stat_cor(method="pearson")+
  facet_wrap(~ node, scale="free", nrow=1) +
  theme_minimal() +theme(legend.position = "none")+
  labs(x = "Scaled Taxa %",
       y = "Δ pH", color = "")


# :: All F. prausnitzii  ------------------------------------------------------

# We know that certain butyrogens have sig differential associations with pH
# within the 2 groups
# How about individualized analyses?

# reduce to prevalent F. prausnitzii

oars.rapidaim.phyloseq.butyrogens.mat.pa = oars.rapidaim.phyloseq.butyrogens.mat
oars.rapidaim.phyloseq.butyrogens.mat.pa[oars.rapidaim.phyloseq.butyrogens.mat.pa!=0] = 1
oars.rapidaim.phyloseq.butyrogens.mat.filtered = oars.rapidaim.phyloseq.butyrogens.mat[,colSums(oars.rapidaim.phyloseq.butyrogens.mat.pa)>nrow(oars.rapidaim.phyloseq.butyrogens.mat.pa)*0.2]

oars.rapidaim.all.butyrogens.associations = do.call(rbind, lapply(unique(oars.rapidaim.phyloseq.butyrogens.map$sample), function(hm){
  # subset to f.prau and samples
  samples = subset(oars.rapidaim.phyloseq.butyrogens.map, sample == hm)
  fpau.data = oars.rapidaim.phyloseq.butyrogens.mat.filtered[samples$Barcode,
                                                             grepl("prausnitzii", colnames(oars.rapidaim.phyloseq.butyrogens.mat.filtered))] %>% data.frame()
  # correlation
  fpau.data$pH = samples$pH
  Hmisc::rcorr(as.matrix(fpau.data), type="spearman")$r %>% 
    reshape2::melt() %>% data.frame() %>%
    subset(Var1 == "pH" & !is.na(value)) %>%
    subset(Var1 != Var2) %>%
    mutate(sample = hm)
}))
oars.rapidaim.all.butyrogens.associations = merge(oars.rapidaim.all.butyrogens.associations,
                                                  oars.rapidaim.phyloseq.butyrogens.map[,c("sample", "timing")] %>% distinct(), by="sample")
oars.rapidaim.all.butyrogens.associations = tidyr::separate(oars.rapidaim.all.butyrogens.associations, col=sample, into=c("HM", "other"), sep="\\.", remove=F)
oars.rapidaim.all.butyrogens.associations$hm_time = paste(oars.rapidaim.all.butyrogens.associations$HM,
                                                          oars.rapidaim.all.butyrogens.associations$timing, sep="_")
oars.rapidaim.all.butyrogens.associations = subset(oars.rapidaim.all.butyrogens.associations, !grepl("12M", oars.rapidaim.all.butyrogens.associations$timing))

oars.rapidaim.all.butyrogens.associations = reshape2::acast(oars.rapidaim.all.butyrogens.associations,
                                                            hm_time ~ Var2, value.var="value")
oars.rapidaim.all.butyrogens.associations[is.na(oars.rapidaim.all.butyrogens.associations)] = 0

# make annotation for F. prausnitzii group-level association
oars.fprau.important.mat = reshape2::acast(oars.fprau.important,
                                           feature ~ Correlation, value.var="coef") %>% data.frame()
oars.fprau.important.mat$association = ifelse(!is.na(oars.fprau.important.mat$Negative) & !is.na(oars.fprau.important.mat$Positive), "Variable",
                                              ifelse(!is.na(oars.fprau.important.mat$Negative), "Negative", "Positive"))
oars.fprau.important.mat = data.frame(Association = oars.fprau.important.mat$association)
rownames(oars.fprau.important.mat) = rownames(reshape2::acast(oars.fprau.important,
                                                              feature ~ Correlation, value.var="coef") %>% data.frame())
# make annotation for HM responses
oars.fprau.hm.mat = oars.rapidaim.ph.but.cor %>%
  subset(timing != "12M") %>%
  mutate(hm_timing = paste(HM, timing, sep="_")) %>%
  mutate(Group = ifelse(cor.sign == 1, "Positive", "Negative"))
rownames(oars.fprau.hm.mat) = oars.fprau.hm.mat$hm_timing
oars.fprau.hm.mat = oars.fprau.hm.mat %>% dplyr::select(Group)

pheatmap::pheatmap(t(oars.rapidaim.all.butyrogens.associations),
                   color=colorRampPalette(c("blue","white", "red"))(100),
                   annotation_row=oars.fprau.important.mat,
                   annotation_colors = list(Association = c(Positive = "red", Negative = "blue", Variable="purple"),
                                            Group = c(Positive = "red", Negative = "blue")),
                   annotation_col=oars.fprau.hm.mat)


# :: All Taxa/Butyrogens PCoA -----------------------------------------------------------------

# all taxa
oars.rapidaim.phyloseq.rare.all = speedyseq::psmelt(oars.rapidaim.phyloseq.rare)
# make taxa names
oars.rapidaim.phyloseq.rare.all$taxa = gsub("_NA", "", ifelse(!is.na(oars.rapidaim.phyloseq.rare.all$Genus),  paste(oars.rapidaim.phyloseq.rare.all$Genus, oars.rapidaim.phyloseq.rare.all$Species, sep="_"),
                                                              ifelse(!is.na(oars.rapidaim.phyloseq.rare.all$Family), oars.rapidaim.phyloseq.rare.all$Family,
                                                                     ifelse(!is.na(oars.rapidaim.phyloseq.rare.all$Order), oars.rapidaim.phyloseq.rare.all$Order,
                                                                            ifelse(!is.na(oars.rapidaim.phyloseq.rare.all$Class), oars.rapidaim.phyloseq.rare.all$Class,
                                                                                   ifelse(!is.na(oars.rapidaim.phyloseq.rare.all$Phylum), oars.rapidaim.phyloseq.rare.all$Phylum,
                                                                                          oars.rapidaim.phyloseq.rare.all$Kingdom))))))
# make matrix
oars.rapidaim.phyloseq.rare.all$Abundance = oars.rapidaim.phyloseq.rare.all$Abundance / 50000
length(unique(oars.rapidaim.phyloseq.rare.all$taxa))
# make the dataframe a bit smaller
oars.rapidaim.phyloseq.rare.all = subset(oars.rapidaim.phyloseq.rare.all, RS_Name %in% rs.names.pbs)
oars.rapidaim.phyloseq.rare.all.mat = reshape2::acast(oars.rapidaim.phyloseq.rare.all,
                                                      Sample ~ OTU, value.var="Abundance")


dim(oars.rapidaim.phyloseq.rare.all.mat) # 7164 detected OTUs

# filter (present in 10% of samples)
oars.rapidaim.phyloseq.rare.all.mat.pa = oars.rapidaim.phyloseq.rare.all.mat
oars.rapidaim.phyloseq.rare.all.mat.pa[oars.rapidaim.phyloseq.rare.all.mat.pa !=0] = 1
oars.rapidaim.phyloseq.rare.all.mat.filtered = oars.rapidaim.phyloseq.rare.all.mat[,colSums(oars.rapidaim.phyloseq.rare.all.mat.pa) > nrow(oars.rapidaim.phyloseq.rare.all.mat)*0.1]
dim(butyrogen.data.pcoa) # 1757 OTUs remaining

butyrogen.data.pcoa = oars.rapidaim.phyloseq.rare.all.mat.filtered[rownames(oars.rapidaim.phyloseq.rare.all.mat.filtered) %in% rownames(subset(oars.rapidaim.phyloseq.butyrogens.map, 
                                                                                                                                               timing %in% c("0M", "3M", "6M"))),]
# subset to F.prau
#butyrogen.data.pcoa = butyrogen.data.pcoa[,grepl("prausn", colnames(butyrogen.data.pcoa))]

butyrogen.bray = vegan::vegdist(butyrogen.data.pcoa, method="bray") 
# perform PCoA
butyrogen.pcoa = ape::pcoa(butyrogen.bray)
# extract data from pcoa
butyrogen.pcoa.df = data.frame(butyrogen.pcoa$vectors[,c(1:2)])
# add metadata
butyrogen.pcoa.df = merge(butyrogen.pcoa.df,
                          oars.rapidaim.phyloseq.butyrogens.map, by="row.names")
# extract variance explained
butyrogen.pcoa.var = butyrogen.pcoa$values[c(1:2),2]
butyrogen.pcoa.df$var1 = round(butyrogen.pcoa.var[1]*100, digits=2)
butyrogen.pcoa.df$var2 = round(butyrogen.pcoa.var[2]*100, digits=2)
rownames(butyrogen.pcoa.df) = butyrogen.pcoa.df$Row.names
# clean up "oars.on.rs" variable

# add correlation of HM
butyrogen.pcoa.df = merge(butyrogen.pcoa.df,
                          oars.rapidaim.ph.but.cor.data, by="sample")
# make sign
butyrogen.pcoa.df = butyrogen.pcoa.df %>%
  mutate(cor.sign = as.factor(sign(cor)))

set.seed(25)
t1 = Sys.time()
butyrogen.pcoa.permanova = vegan::adonis2(butyrogen.bray ~ cor.sign,
                                          butyrogen.pcoa.df,
                                          strata = butyrogen.pcoa.df$sample,
                                          by="margin")
t2 = Sys.time()
t2 - t1
butyrogen.pcoa.permanova # not sig

butyrogen.pcoa.plot <- ggplot(
  data=butyrogen.pcoa.df %>% group_by(sample),
  aes(x=Axis.1, y=Axis.2))+
  ggpubr::stat_chull(geom="polygon", aes(group=sample), fill="black", alpha=0.2)+
  #geom_path(aes(group=sample), color="black", linetype=2, alpha=0.5, linewidth=0.3) + 
  geom_point(aes(fill=cor, shape=cor.sign), size=2)+
  #stat_ellipse(aes(group=response, color=response), alpha=0.5)+
  scale_fill_gradient2(low="blue", high="red")+
  scale_shape_manual(values=c(25,24))+
  annotate(geom="text", x=-Inf, y=Inf, hjust=-0.1, vjust=1.2,
           label=paste(paste("R²: ", round(data.frame(butyrogen.pcoa.permanova)[1,3], 3)*100, "%",
                             "  p: ", round(data.frame(butyrogen.pcoa.permanova)[1,5], 3), sep=""), sep=""))+
  facet_wrap(~"Bray-Curtis PCoA")+
  theme_minimal()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(color="black"))+
  labs(x=paste("Axis 1: ", round(unique(butyrogen.pcoa.df$var1), digits=2), "%", sep=""), 
       y=paste("Axis 2: ", round(unique(butyrogen.pcoa.df$var2), digits=2), "%", sep=""))
butyrogen.pcoa.plot
