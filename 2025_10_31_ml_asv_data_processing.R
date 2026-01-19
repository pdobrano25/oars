### 2025_10_31  Process ASV files for ML on ComputeCanada

## from ml_git_archfolder,
#  sbatch ml_git_main/asv_prep

## Goal: prepare asv tables for ML analysis
## Notes: 
# - decontam using amplicon concentrations
# - taxonomy using gg2; no longer gg_13_8
# - rarefy to 50,000 reads / sample
# - save easy-to-use files
# - calculate target variables

steps = c(4:7)

seed = 25

## Note:
# this code will generate larger objects if all HMs available are included; 
# depends on objects in /dada2_pooledSequences/ folder


# upload data via: scp -r Documents/PhD/git_ml_archfolder/ml_git_data/ pdobrano@cedar.computecanada.ca:scratch/ml_git_archfolder
# upload code via: scp -r Documents/PhD/git_ml_archfolder/ml_git_main/ pdobrano@cedar.computecanada.ca:scratch/ml_git_archfolder

# go to: cd scratch/ml_git_archfolder/ml_git_main
#.   setwd("~/Documents/PhD/git_ml_archfolder/ml_git_data/")
# local: setwd("~/Documents/PhD/git_ml_archfolder/ml_git_data/")

date = "2025_11_07" # this checks for mapping.file version + saves file accordingly
print(date)
# steps = c(1, # prepare
#           2, # remove chimeras
#           3, # decontam
#           4, # assign tax
#           5, # rarefy
#           6) # export
#           7) # variables

# packages ----------------------------------------------------------------
library("phyloseq");
print("phyloseq loaded")
library("reshape2"); 
print("reshape2 loaded")
library("dplyr"); 
print("dplyr loaded")
library("tidyr");
print("tidyr loaded")
library("decontam");
print("decontam loaded")
library("parallel");
print("parallel loaded")
library("speedyseq")
print("speedyseq loaded")

print("all packages loaded")

# directory (local)
# setwd("~/Documents/PhD/git_ml_archfolder/")

# directory (remote)
#setwd("~/scratch/ml_git_archfolder/")

# rs names
rs.names <- c("Authentic", "BobsRedMill", "MSPrebiotic", "LetsDoOrganic", "HiMaize260", "Novelose330", "ActistarRT", "FibersymRW", "Versafibe1490")
rs.names.pbs <- c("PBS", "Authentic", "BobsRedMill", "MSPrebiotic", "LetsDoOrganic", "HiMaize260", "Novelose330", "ActistarRT", "FibersymRW", "Versafibe1490")

print("navigated to directory")

t0 <- Sys.time()

# load mapping file -------------------------------------------------------
# directory
final.mapping.file <- read.csv("../ml_git_data/2024_04_15_full_mapping_rapidaim_covariates.csv")
head(final.mapping.file)
# clean mapping file
all.mapping.df <- subset(final.mapping.file, Seq_run != "AS16S..gg") # delete runs not yet completed
all.mapping.df <- all.mapping.df %>% distinct() # remove potential duplicates
row.names(all.mapping.df) <- all.mapping.df$Barcode # add rownames
# Any repeated barcodes?
repeats <- nrow(all.mapping.df[all.mapping.df$Barcode %in% data.frame(table(all.mapping.df$Barcode))$Var1[data.frame(table(all.mapping.df$Barcode))$Freq > 1],])
print(paste0(repeats, " repeated barcodes."))
# clean up barcode name to match Biom file
row.names(all.mapping.df) <- gsub("\\.gg13v5", "", rownames(all.mapping.df))

# subset to HMs for ML analysis (first non-repeat RapidAIM culture)
ml.hms = read.csv("../ml_git_data/2024_03_15_hm_list.csv")
all.mapping.df <- subset(all.mapping.df, HM %in% ml.hms$HM)

print("loaded mapping file")

# extract AS16S runs (to filter next step)
# ALL runs
# asv.runs <- paste("AS16S_", gsub("AS16S.", "", unique(final.mapping.file$Seq_run)), sep="")
# SPECIFIC runs
# find ~/projects/def-astintzi/dada2/dada2_pooledSequences/ |grep -f ./lib_names.txt
# xargs -a <(find ~/projects/def-astintzi/dada2/dada2_pooledSequences/ |grep -f ~/scratch/ml_git_archfolder/ml_git_data/lib_names.txt) cp -t ~/scratch/ml_git_archfolder/ml_git_data/ml_dada2/
asv.runs <- read.csv("../ml_git_data/2023_03_13_ml_list_of_sequence_runs.csv")
# convert to vector
asv.runs <- paste("AS16S_", gsub("AS16S.", "", unique(asv.runs$seq)), sep="")


if (sum(grepl(1, steps)>0)) {
  # load asv files ----------------------------------------------------------
  # directory
  all.asv.list <- list.files("../ml_git_data/ml_dada2/", recursive=TRUE)
  all.asv.list <- all.asv.list[grep(paste(asv.runs, collapse="|"), all.asv.list)]
  all.asv.list <- all.asv.list[!grepl("compressed", all.asv.list)]
  
  print("loaded asv files")
  
  # 1. prepare asv files ------------------------------------------------------------
  
  # load data as list of bioms
  all.asv.ps <- lapply(1:length(all.asv.list), function(x) {
    all.asv.ps.list <- readRDS(file.path("../ml_git_data/ml_dada2/", all.asv.list[x]))
  })
  print("listed asv files")
  
  # convert bioms to phyloseq otu_tables
  all.asv.ps.clean <- lapply(all.asv.ps, function(x){
    convert.to.otutable =  phyloseq::otu_table(x, taxa_are_rows=FALSE)
    return(convert.to.otutable)
  })
  
  # ensure otu_table sample IDs match mapping.file Barcodes
  all.asv.ps.clean.again <- lapply(all.asv.ps.clean, function(x){
    rownames(x) =       gsub("_filtered\\.fastq\\.gz", "", rownames(x))
    rownames(x) =          gsub("_", "\\.", rownames(x))
    return(x)
  })

  t1 <- Sys.time()
  
  print(paste0("1. prepared asv files. ", t1-t0, " elapsed. ", t1-t0, " total."))
  
  # merge into a single otu_table; note: this crashes R locally
  all.asv.ps.merge <- do.call(merge_phyloseq, all.asv.ps.clean.again)
  
  print(paste0("1. and merged. ", t1-t0, " elapsed. ", t1-t0, " total."))
  
  saveRDS(all.asv.ps.merge, "../ml_git_data/asv_prep_step_1.rds")
}
# Works

if (sum(grepl(2, steps)>0)) {
  
  # 2. remove chimeras ------------------------------------------------------
  all.asv.ps.merge <- readRDS("../ml_git_data/asv_prep_step_1.rds")
  
  print("2. begun chimera identification")
  
  # remove chimeras on merged ASV Table
  
  t1 <- Sys.time()
  
  set.seed(seed)
  nonchimeras.to.be.kept = dada2::removeBimeraDenovo(
    dada2::getUniques(all.asv.ps.merge), 
    method="consensus", 
    multithread=F,  # note: must be F, otherwise cores are likely to fail and results don't reproduce!
    verbose=FALSE)
  
  # remove these bimeras
  all.asv.ps.nobim <- all.asv.ps.merge[,colnames(all.asv.ps.merge) %in% names(nonchimeras.to.be.kept)]
  
  # how many sequences + reads kept?
  asvs.kept.number <- length(nonchimeras.to.be.kept) # 
  asvs.kept.percentage <- length(nonchimeras.to.be.kept)  / 
    ncol(phyloseq::otu_table(all.asv.ps.merge)) #
  reads.kept.percentage <- sum(colSums(phyloseq::otu_table(all.asv.ps.nobim)))/
    sum(colSums(phyloseq::otu_table(all.asv.ps.merge))) # 
  
  print(paste("asvs kept (number): ", asvs.kept.number, sep=""))
  print(paste("asvs kept (perc): ", asvs.kept.percentage, sep=""))
  print(paste("reads kept (perc): ", reads.kept.percentage, sep=""))
  
  t2 <- Sys.time() # ~ 2 hours for n=117
  print(paste0("2. removed chimeras. ", t2-t1, " elapsed. ", t2-t0, " total."))
  saveRDS(all.asv.ps.nobim, "../ml_git_data/asv_prep_step_2.rds")
  
}

if (sum(grepl(3, steps)>0)) {
  # 3. decontam  -------------------------------------------------------------
  t2 <- Sys.time()
  
  all.asv.ps.nobim <- readRDS("../ml_git_data/asv_prep_step_2.rds")
  
  print("prepare for decontam")
  
  # prepare phyloseq object with metadata
  all.asv.ps.nobim <- phyloseq::phyloseq(otu_table = otu_table(all.asv.ps.nobim),
                                         sample_data = sample_data(all.mapping.df))
  
  decontam.batches <- as.factor(phyloseq::sample_data(all.asv.ps.nobim)$HM)
  
  # run decontam in manual batches (do not use `batch = ` option)
  contaminants.df = do.call(rbind, lapply(unique(decontam.batches), function(x){
    # x=unique(decontam.batches)[2]
    #print(paste(x))
    # subset to each HM (needs workaround because of phyloseq bug)
    samples.subset = data.frame(phyloseq::sample_data(all.asv.ps.nobim)) %>%
      subset(HM == x)
    ps.subset = phyloseq::prune_samples(rownames(samples.subset), all.asv.ps.nobim)
    # delete taxa not there
    ps.subset = phyloseq::prune_taxa(phyloseq::taxa_sums(ps.subset)>0, ps.subset)
    # run decontam
    contaminants = decontam::isContaminant(ps.subset, 
                                           conc = "Post_PCR", 
                                           detailed = TRUE)
    # record in dataframe
    contaminants.df = data.frame(feature = rownames(contaminants),
                                 contaminant = contaminants$contaminant,
                                 HM = x)
    # return df
    contaminants.df
  }))
  
  contaminants = unique(subset(contaminants.df, contaminant == TRUE)$feature)
  
  print(paste("contaminants = ", length(contaminants), sep="")) # 1955
  
  asvs.before.decontam <- ntaxa(all.asv.ps.nobim)
  
  # apply decontam
  all.asv.ps.decontam = phyloseq::prune_taxa(subset(contaminants.df, contaminant == FALSE)$feature, all.asv.ps.nobim)
  
  t3 <- Sys.time()
  print(paste0("3. list asv tables. ", t3-t2, " elapsed. ", t3-t0, " total."))
  
  saveRDS(all.asv.ps.decontam, "../ml_git_data/asv_prep_step_3.rds")
}
# Works

if (sum(grepl(4, steps)>0)) {
  
  # 4. assign tax --------------------------------------------------------------
  t3 <- Sys.time()
  
  all.asv.ps.decontam = readRDS("../ml_git_data/asv_prep_step_3.rds")
  
  # remove control samples
  all.asv.ps.decontam = phyloseq::subset_samples(all.asv.ps.decontam,
                                                 !RS_Name %in% c("MED", "WaterB", "ZymoC", "ZymoD", "Spike-in","ExtractC","Spike-In"))
  # remove taxa not present
  all.asv.ps.decontam = phyloseq::prune_taxa(names(phyloseq::taxa_sums(all.asv.ps.decontam)[phyloseq::taxa_sums(all.asv.ps.decontam) > 0]), all.asv.ps.decontam)
  # removed ~ 3000 taxa
  
  # using greengenes 13.8
  set.seed(seed)
  asv.taxa <- dada2::assignTaxonomy(colnames(otu_table(all.asv.ps.decontam)), 
                                   #"../ml_git_data/gg_13_8_train_set_97.fa.gz", 
                                   "../ml_git_data/gg2_2024_09_toSpecies_trainset.fa.gz", 
                                    multithread=TRUE,
                                    minBoot = 50) # >50% of bootstrap results support call
  # note: assignTax is not 100% deterministic; 
  # ~1% of ASVs will vary 50% of the time (even with set.seed)
  
  physeq.all.taxa.table <- tax_table(asv.taxa)
  taxa_names(physeq.all.taxa.table) <- rownames(physeq.all.taxa.table)
  colnames(physeq.all.taxa.table) <- c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species")
  all.asv.ps.taxa <- merge_phyloseq(all.asv.ps.decontam, 
                                    tax_table=(physeq.all.taxa.table))
  asvs.after.decontam <- ntaxa(all.asv.ps.taxa)
  all.asv.ps.taxa <- phyloseq::subset_taxa(all.asv.ps.taxa, !Genus %in% c("g__Mitochondria")) # uppercase for gg2, lowercase for gg138
  all.asv.ps.taxa <- phyloseq::subset_taxa(all.asv.ps.taxa, !Family %in% c("f__Mitochondria")) 
  asvs.as.mitochondria <- ntaxa(all.asv.ps.taxa) # tabulate mitochondria
  all.asv.ps.taxa <- phyloseq::subset_taxa(all.asv.ps.taxa, !Class  %in% c("c__Chloroplast")) #
  asvs.as.chloroplast <- ntaxa(all.asv.ps.taxa) # tabulate chloroplast
  all.asv.ps.taxa
  t4 <- Sys.time()
  
  print(paste0("4. tax assignment complete. ", t4-t3, " elapsed. ", t4-t0, " total."))
  
  saveRDS(all.asv.ps.taxa, "../ml_git_data/asv_prep_step_4.rds")
  
  # check sample sums
  check.depths <- data.frame(depth = phyloseq::sample_sums(all.asv.ps.taxa),
                             HM = phyloseq::sample_data(all.asv.ps.taxa)$HM,
                             RS_Name = phyloseq::sample_data(all.asv.ps.taxa)$RS_Name) %>%
    subset(RS_Name %in% c(rs.names, "PBS", "Stool"))
  # good
  
}

if (sum(grepl(5, steps)>0)) {
  # 5. rarefy ------------------------------------------------------------------
  all.asv.ps.taxa <- readRDS("../ml_git_data/asv_prep_step_4.rds")
  t4 <- Sys.time()
  
  # Rarefy to 50,000 reads / sample
  set.seed(seed)
  samples.before.rare <- nsamples(all.asv.ps.taxa)
  physeq.all.taxa.rare <- rarefy_even_depth(all.asv.ps.taxa, replace=FALSE, sample.size=50000)
  samples.after.rare <- nsamples(physeq.all.taxa.rare)
  t5 <- Sys.time()
  print(paste0("5. rarefaction complete. ", t5-t4, " elapsed. ", t5-t0, " total."))
  physeq.all.taxa.tss <- transform_sample_counts(physeq.all.taxa.rare, function(x) 100 * x/sum(x)) 
  t6 <- Sys.time()
  print(paste0("5. TSS complete. ", t6-t5, " elapsed. ", t6-t0, " total."))
  print(physeq.all.taxa.tss)
  
  sample.size = length(unique(sample_data(physeq.all.taxa.tss)$HM))
  paste0("sample size = ", sample.size)
  
  # directory
  saveRDS(physeq.all.taxa.tss, paste0("../ml_git_data/", date, "_n", sample.size,"_s", seed, "_asv_final_50k.rds"))
}

if (sum(grepl(6, steps)>0)) {
  
  # 6. export all -------------------------------------------------------------
  print("Read in file")
  phyloseq.file.name <- list.files(path = "../ml_git_data/", pattern = paste0(date, "_n", sample.size,"_s", seed, "_asv_final_50k.rds", sep=""))
  print(paste0(phyloseq.file.name))
  physeq.all.taxa.tss <- readRDS(paste0("../ml_git_data/", phyloseq.file.name, sep=""))
  sample.size = length(unique(sample_data(physeq.all.taxa.tss)$HM))
  
  t7 <- Sys.time()
  
  # Reduce to rs.names.pbs + stool + FOS
  physeq.all.taxa.tss.rs <- subset_samples(physeq.all.taxa.tss, RS_Name %in% c("Authentic", "BobsRedMill","MSPrebiotic", "LetsDoOrganic","HiMaize260", "Novelose330", "ActistarRT", "FibersymRW", "Versafibe1490", "Stool", "RawStool", "PBS", "FOS"))
  asvs.total = ntaxa(physeq.all.taxa.tss.rs)
  
  print("Add Glom + LCA")
  # Add LCA
  asv.rare.tax.df <- data.frame(tax_table(physeq.all.taxa.tss))
  # remove empty names
  asv.rare.tax.df = asv.rare.tax.df %>%
    mutate(Species = ifelse(asv.rare.tax.df$Species == "s__", NA, asv.rare.tax.df$Species)) %>%
    mutate(Genus = ifelse(asv.rare.tax.df$Genus == "g__", NA, asv.rare.tax.df$Genus)) %>%
    mutate(Family = ifelse(asv.rare.tax.df$Family == "f__", NA, asv.rare.tax.df$Family)) %>%
    mutate(Order = ifelse(asv.rare.tax.df$Order == "o__", NA, asv.rare.tax.df$Order)) %>%
    mutate(Class = ifelse(asv.rare.tax.df$Class == "c__", NA, asv.rare.tax.df$Class)) %>%
    mutate(Phylum = ifelse(asv.rare.tax.df$Phylum == "p__", NA, asv.rare.tax.df$Phylum)) %>% data.frame()
    # assign LCA
  asv.rare.tax.df <- asv.rare.tax.df %>%
    mutate(Taxa = 
             ifelse(!is.na(Species)&Species!="s__", paste(as.character(Genus), as.character(Species), sep="_"), 
                    ifelse(!is.na(Genus)&Genus!="g__", paste(Genus), 
                           ifelse(!is.na(Family)&Family!="f__", paste(Family),
                                  ifelse(!is.na(Order)&Order!="o__", paste(Order),
                                         ifelse(!is.na(Class)&Class!="c__", paste(Class),
                                                paste(Phylum)))))))
  # append unique number to duplicated ASVs
  asv.rare.tax.df$Glom <- asv.rare.tax.df$Taxa
  glom.total = length(unique(asv.rare.tax.df$Glom))
  
  asv.rare.tax.df$Taxa <- make.unique(as.character(asv.rare.tax.df$Taxa), sep = "_")
  asv.rare.tax.df$OTU <- rownames(asv.rare.tax.df) # call them OTUs, not ASVs, so they can merge later
  
  # save this tax file
  saveRDS(asv.rare.tax.df, paste0("../ml_git_data/", date, "_", "n", sample.size,"_s", seed, "_asv_final_tax.rds"))
  
  print("Save Glom and LCA")
  physeq.all.taxa.tss.rs.glom <- physeq.all.taxa.tss.rs
  tax_table(physeq.all.taxa.tss.rs.glom) <- as.matrix(asv.rare.tax.df[,c("Kingdom", "Glom")])
  physeq.all.taxa.tss.rs.glom <- speedyseq::tax_glom(physeq.all.taxa.tss.rs.glom, taxrank = "Glom")
  # Glommed LCA
  physeq.all.taxa.tss.rs.lca <- physeq.all.taxa.tss.rs
  tax_table(physeq.all.taxa.tss.rs.lca) <- as.matrix(asv.rare.tax.df[,c("Kingdom", "Taxa")])
  physeq.all.taxa.tss.rs.lca <- speedyseq::tax_glom(physeq.all.taxa.tss.rs.lca, taxrank = "Taxa")
  # LCA ASV
  saveRDS(physeq.all.taxa.tss.rs.glom, paste0("../ml_git_data/", date, "_", "n",sample.size,"_s", seed, "_asv_final_glom_50k.rds"))
  saveRDS(physeq.all.taxa.tss.rs.lca, paste0("../ml_git_data/", date, "_", "n",sample.size,"_s", seed, "_asv_final_lca_50k.rds"))
  
  print("Save Median Glom and LCA")
  # Median Glom
  physeq.all.taxa.tss.rs.glom.df <- speedyseq::psmelt(physeq.all.taxa.tss.rs.glom)
  physeq.all.taxa.tss.rs.glom.df$Code <- paste(physeq.all.taxa.tss.rs.glom.df$HM, physeq.all.taxa.tss.rs.glom.df$RS_Name, sep="_")
  physeq.all.taxa.tss.rs.glom.median.mat <- reshape2::acast(physeq.all.taxa.tss.rs.glom.df, Code ~ Glom, value.var="Abundance", fun.aggregate=median)
  # Median glommed LCA
  
  # Median LCA
  physeq.all.taxa.tss.rs.lca.df <- speedyseq::psmelt(physeq.all.taxa.tss.rs.lca)
  physeq.all.taxa.tss.rs.lca.df$Code <- paste(physeq.all.taxa.tss.rs.lca.df$HM, physeq.all.taxa.tss.rs.lca.df$RS_Name, sep="_")
  physeq.all.taxa.tss.rs.lca.median.mat <- reshape2::acast(physeq.all.taxa.tss.rs.lca.df, Code ~ Taxa, value.var="Abundance", fun.aggregate=median)
  
  # Save
  saveRDS(physeq.all.taxa.tss.rs.glom.median.mat, paste0("../ml_git_data/", date, "_", "n",sample.size,"_s", seed, "_asv_final_median_glom_50k.rds"))
  saveRDS(physeq.all.taxa.tss.rs.lca.median.mat, paste0("../ml_git_data/", date, "_", "n",sample.size,"_s", seed, "_asv_final_median_lca_50k.rds"))
  
  print("Save Mapping")
  # Lastly, save final mapping.files
  physeq.all.taxa.tss.rs.mapping <- physeq.all.taxa.tss.rs.glom.df[,c("HM","HMcode", "RS_Name", "Replicate", "pH", "Qubit", "Seq_run", "Post_PCR", "TIME","TIME.2", "STUDY","STUDY.2", "REDCap.Stool.ID", "STL_ID", "diagnosis", "gender", "age")] %>% distinct()
  saveRDS(physeq.all.taxa.tss.rs.mapping, paste0("../ml_git_data/", date, "_", "n",sample.size,"_s", seed, "_asv_final_mapping.rds"))
  
  t8 <- Sys.time()
  
  print(paste0("6. Export complete. ", t8-t7, " elapsed. ", t8-t0, " total."))
  
  
}

if (sum(grepl(7, steps)>0)) {
  
  # 7. variables ------------------------------------------------------------
  # load data (remotely)
  sample.size = 117
  t8 <- Sys.time()
  
  # load data (remotely)
  ph.data <- readRDS(paste0("../ml_git_data/", date, "_", "n", sample.size,"_s", seed, "_asv_final_mapping.rds"))
  
  physeq.all.taxa.tss.rs.glom <- readRDS(paste0("../ml_git_data/", date, "_", "n", sample.size,"_s", seed, "_asv_final_glom_50k.rds"))
  phyloseq.file.name <- list.files(path = "../ml_git_data/", pattern = paste0(date, "_", "n", sample.size,"_s", seed, "_asv_final_50k.rds", sep=""))
  physeq.all.taxa.tss <- readRDS(paste0("../ml_git_data/", phyloseq.file.name, sep=""))
  physeq.all.taxa.tss.rs.glom.median.mat <- readRDS(paste0("../ml_git_data/", date, "_", "n", sample.size,"_s", seed, "_asv_final_median_glom_50k.rds"))
  
  sample.size = length(unique(phyloseq::sample_data(physeq.all.taxa.tss)$HM))
  
  # :: pH ----------------------------------------------------------------------
  print("Calculate pH data")
  
  # load data (locally)
  # ph.data <- readRDS("../ml_git_data/2024_03_13_n117_asv_final_mapping.rds")
  
  # subset to important RS
  ph.data <- subset(ph.data, RS_Name %in% c(rs.names, "PBS"))
  
  # calculate median value
  ph.data <- ph.data %>%
    group_by(HM, RS_Name) %>%
    mutate(med.ph = median(na.omit(pH))) %>% 
    dplyr::select(HM, RS_Name, med.ph) %>% distinct() %>% data.frame()
  
  # calculate change from PBS
  ph.data.pbs <- subset(ph.data, RS_Name == "PBS")[,c("HM", "med.ph")] %>% distinct()
  colnames(ph.data.pbs) <- c("HM", "PBS")
  ph.data <- merge(ph.data, ph.data.pbs, by=c("HM"))
  ph.data$delta_ph = ph.data$med.ph - ph.data$PBS
  
  # finalize data
  ph.data <- subset(ph.data, RS_Name %in% c(rs.names, "PBS"))[,c("HM", "RS_Name", "delta_ph")]
  ph.data
  
  
  # > processing report ----------------------------------------------------
  
  processing.report   = list("asvs.kept.number" = print(asvs.kept.number),
                                "asvs.kept.percentage" = print(asvs.kept.percentage),
                                "reads.kept.percentage" = print(reads.kept.percentage),
                                "asvs.before.decontam" = print(asvs.before.decontam),
                                "asvs.after.decontam" = print(asvs.after.decontam),
                                "asvs.as.mitochondria" = print(asvs.as.mitochondria),
                                "asvs.as.chloroplast" = print(asvs.as.chloroplast),
                                "samples.before.rare"= print(samples.before.rare),
                                "samples.after.rare" = print(samples.after.rare),
                                "asvs.total" = print(asvs.total),
                                "glom.total" = print(glom.total))
  t9 <- Sys.time()
  
  print(t9-t8)
  
  saveRDS(processing.report, paste0("../ml_git_data/", date, "_", "n", sample.size,"_s", seed, "_processing_report.rds"))
  
}

# Works!

# Download via
# scp pdobrano@nibi.alliancecan.ca:'~/scratch/rs_ml/ml_git_data/*.rds' Documents/PhD/git_ml_archfolder/ml_git_data/

