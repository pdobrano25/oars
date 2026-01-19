### 2025_05_24  OARS DNA Analysis

# Goal: 
# 1. Standard: Perform standard approach for: water, fcal, 16S, and MGX (alpha, beta, maaslin) (0-6m, 6m-12m)
# 2. PostHoc Response: Perform posthoc analysis on RapidAIM pH responses


# :: save -----------------------------------------------------------------

# Update notes: Split into 2 papers (Standard + Multomic)
# Update 2 notes: Coalesce into dissertation chapter (they can figure it out later)
#  
v.date = "2026_01_19" # this object gets saved to indicate latest version
# save.image("./2026_01_19_oars_analysis.Renv")

# load("./2026_01_19_oars_analysis.Renv")
# note: 2016_01_19 is a bare analysis; only this script and no defunct parts

# :: load packages --------------------------------------------------------
library("ggplot2"); library("dplyr"); library("tidyverse"); library("patchwork")

# load processed R objects
load(file = "./2026_01_13_oars_16s_data_meta.Renv")


# :: color vectors ---------------------------------------------------------------

rs.names.pbs = c("PBS", "Authentic", "BobsRedMill", "MSPrebiotic", "LetsDoOrganic", "HiMaize260", "Novelose330", "ActistarRT","FibersymRW", "Versafibe1490")
rs.names = c("Authentic", "BobsRedMill", "MSPrebiotic", "LetsDoOrganic", "HiMaize260", "Novelose330", "ActistarRT","FibersymRW", "Versafibe1490")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
labelcolors = data.frame(cols = c(gg_color_hue(5)[c(1,1,1,2,3,3,4,4,5)], "#000000"),
                         RS_Name = rs.names.pbs[c(2:10,1)])
# used for nearly all figures
rs.colors = c("Authentic" = "#F8766D", "BobsRedMill" = "#F8766D", "MSPrebiotic" = "#F8766D",
              "LetsDoOrganic" = "#A3A500", "HiMaize260" = "#00BF7D", "Novelose330" = "#00BF7D",
              "ActistarRT" = "#00B0F6", "FibersymRW" = "#00B0F6", "Versafibe1490" = "#E76BF3")
# used for pheatmap, network, and ML
omics.colors = c("ASV" = "#8DD3C7",
                 "Species" = "#FFFFB3",
                 "Pathway" =  "#BEBADA",
                 "COG" = "#FB8072",
                 "CAZy" = "#80B1D3",
                 "Metabolite" = "#FDB462")
# used for heatmap and network
omics.shapes = c("ASV" = 21,
                 "Species" = 22,
                 "Pathway" =  23,
                 "COG" = 23,
                 "CAZy" = 23,
                 "Metabolite" = 24)

# add rownames to data
rownames(metadata.oars.stool) = metadata.oars.stool$standard.name
rownames(metadata.oars.stool.asv) = metadata.oars.stool.asv$standard.name


metadata.oars.stool$plate

# >> Functions ------------------------------------------------------------

# this function clones longitudinal data such that
# a 3M timepoint becomes a baseline timepoint to the 6M
# so we can have up to 2x "before" and "after" comparisons per patient
doubleup = function(data = beta.trajectory.data,
                    delta = FALSE,
                    compliant = TRUE){

  # subset to 0, 3, 6 timepoints
  # separate phases
  data.03 = subset(data, timing %in% c("0M", "3M"))
  data.36 = subset(data, timing %in% c("3M", "6M"))
  # add dummy vars for missing stool (i.e. HM0759, HM0924)
  data.missing.6 = table(data[,c("HM", "timing")]) %>% data.frame() %>% subset(timing == "6M" & Freq == 0)
  data.36 = rbind(data.36,
               data.frame(subset(data.03, HM %in% data.missing.6$HM & timing == "3M") %>% 
                            mutate(timing = "6M") %>%
                            mutate(phase = "rs2"))) %>% data.frame()
  # delta's will be 0, will be removed after
  
  # make new timing name
  data.03$reltiming = ifelse(data.03$timing == "0M", "pre", "post")
  data.03$phase = "rs1"
  data.36$reltiming = ifelse(data.36$timing == "3M", "pre", "post")
  data.36$phase = "rs2"
  # fix target variable
  data.03 = data.03 %>%
    group_by(HM) %>%
    mutate(delta.ph = ifelse(timing == "3M", delta.ph[timing == "0M"],
                             ifelse(timing == "6M", delta.ph[timing == "3M"], NA))) %>%
    mutate(RS_Name = ifelse(timing == "3M", rs.selected[timing == "0M"],
                             ifelse(timing == "6M", rs.selected[timing == "3M"], NA))) 
  data.36 = data.36 %>%
    group_by(HM) %>%
    mutate(delta.ph = ifelse(timing == "3M", delta.ph[timing == "0M"],
                             ifelse(timing == "6M", delta.ph[timing == "3M"], NA))) %>%
    mutate(RS_Name = ifelse(timing == "3M", rs.selected[timing == "0M"],
                            ifelse(timing == "6M", rs.selected[timing == "3M"], NA))) 
    
  # merge
  data.036 = rbind(data.03, data.36) %>% data.frame()
  # clean # data.036 = data.036[,c("HM", "timing", "delta.ph", "rs.selected", "phase", "reltiming", "fcal")]
  data.036$phase = factor(data.036$phase, levels=c("rs1", "rs2"))
  data.036$reltiming = factor(data.036$reltiming, levels=c("pre", "post"))
  # if delta, output delta
  if(delta == TRUE){
    if(sum(colnames(data.036) == "fcal")==1){
      data.036 = data.036 %>%
        group_by(HM, phase) %>%
        mutate(lfc.fcal = log2(fcal[reltiming == "post"]) - (log2(fcal[reltiming == "pre"]))) %>%
        # add baseline value
        mutate(baseline.fcal = ((fcal[reltiming == "pre"])))
    }
    if(sum(colnames(data.036) == "stool_water_perc")==1){
      data.036 = data.036 %>%
        group_by(HM, phase) %>%
        mutate(delta.water = (stool_water_perc)[reltiming == "post"] - ((stool_water_perc))[reltiming == "pre"])%>%
        # add baseline value
        mutate(baseline.water = ((stool_water_perc[reltiming == "pre"])))
    }
    if(sum(colnames(data.036) == "richness")==1){
      data.036 = data.036 %>%
        group_by(HM, phase) %>%
        mutate(delta.richness = log2(richness[reltiming == "post"]) - log2(richness[reltiming=="pre"])) %>%
        # add baseline value
        mutate(baseline.richness = ((richness[reltiming == "pre"])))
    }
    if(sum(colnames(data.036) == "shannon")==1){
      data.036 = data.036 %>%
        group_by(HM, phase) %>%
        mutate(delta.shannon = shannon[reltiming == "post"] - shannon[reltiming == "pre"]) %>%
        # add baseline value
        mutate(baseline.shannon = ((shannon[reltiming == "pre"])))
    }
    if(sum(colnames(data.036) == "fd")==1){
      data.036 = data.036 %>%
        group_by(HM, phase) %>%
        mutate(delta.fd = log2(fd[reltiming == "post"] / fd[reltiming == "pre"])) %>%
        # add baseline value
        mutate(baseline.fd = ((fd[reltiming == "pre"])))
    }
    if(sum(colnames(data.036) == "load.asv")==1){
      data.036 = data.036 %>%
        group_by(HM, phase) %>%
        mutate(delta.load.asv = log2(load.asv[reltiming == "post"] / load.asv[reltiming == "pre"])) %>%
        # add baseline value
        mutate(baseline.load.asv = ((load.asv[reltiming == "pre"])))
    }
    if(sum(colnames(data.036) == "load.mgx")==1){
      data.036 = data.036 %>%
        group_by(HM, phase) %>%
        mutate(delta.load.mgx = log2(load.mgx[reltiming == "post"] / load.mgx[reltiming == "pre"])) %>%
        # add baseline value
        mutate(baseline.load.mgx = ((load.mgx[reltiming == "pre"])))
    }
    if(sum(colnames(data.036) == "med.ph")==1){
      data.036 = data.036 %>%
        group_by(HM, phase) %>%
        mutate(delta.ph = med.ph[reltiming == "post"] - med.ph[reltiming == "pre"]) %>%
        # add baseline value
        mutate(baseline.ph = ((med.ph[reltiming == "pre"])))
    }
    if(sum(colnames(data.036) == "between.beta")==1){
      data.036 = data.036 %>%
        group_by(HM, phase) %>%
        mutate(delta.between = between.beta[reltiming == "post"] - between.beta[reltiming == "pre"]) %>%
        # add baseline value
        mutate(baseline.between = ((between.beta[reltiming == "pre"])))
    }
    if(sum(colnames(data.036) == "within.beta")==1){
      data.036 = data.036 %>%
        dplyr::group_by(HM, phase) %>%
        dplyr::mutate(delta.within = within.beta[reltiming == "post"] - within.beta[reltiming == "pre"]) %>%
        # add baseline value
        dplyr::mutate(baseline.within = ((within.beta[reltiming == "pre"])))
    }
   # data.036 = data.036[,!colnames(data.036) %in% c("standard.name", "stool_date_rec_v2","oars_start_date","oars_product_end_date","oars.timing", "oars.days","oars.off","baseline","stool_water_perc","fcal","timing","oars.on.rs")] %>% distinct()
  }
  # double up delta.ph
  data.036 = data.036 %>%
    group_by(HM, phase) %>%
    mutate(delta.ph = ifelse(is.na(delta.ph), delta.ph[reltiming == "post"], delta.ph)) %>% data.frame()
  # remove non-compliant
  if(compliant == T){
  data.036 = subset(data.036, compliant == TRUE)
  }
  return(data.036)
  
}

# this function applies a prevalence filter to omics data (using only samples of interest)
# applies a log(+pseudo) transform
# and optionally scales to 100%

delta.omic.prepare = function(data = oars.mbx.annotated.mat,
                              prev = 0.20,
                              min.abun = 1,
                              pseudo = T,
                              normalize = F,
                              for_ml = F,
                              double=T){
  # force data.frame
  data = as.data.frame(data)
  # keep samples of interest
  data = data[rownames(data) %in% subset(metadata.oars.stool, timing %in% c("0M", "3M", "6M"))$standard.name,]
  
  # replace NA with 0
  data[is.na(data)] = 0
  
  # scale to 100%
  if(normalize == T){
    data = as.data.frame(((data) / rowSums(data)) * 100)
  }
  if(pseudo==T){
  # add pseudo
  pseudo = min(data[data!=0])/2
  data = data+pseudo
  }
  data$standard.name = rownames(data)
  # merge with meta to add merging variables
  if(double == T){
    data = merge(data,
                        metadata.oars.stool.double[,c("standard.name", "HM", "timing","phase","reltiming", "response", "diagnosis", "adj.fiber")], by=c("standard.name"))
  } else {
    data = merge(data,
                          metadata.oars.stool.double[,c("standard.name", "HM","oars.days", "timing","phase","reltiming", "response", "diagnosis", "adj.fiber")], by=c("standard.name"))
  }  
  # fix colname order
  data.standard.names = data$standard.name
  data$standard.name = NULL
  data$standard.name = data.standard.names
  # apply filter: for ml, do not filter per group; and subset to just baseline samples
  if(for_ml==T){
    data = subset(data, reltiming == "pre")
    data.pa = data[,!colnames(data) %in% colnames(metadata.oars.stool.double)]
    data.pa[data.pa!=pseudo] = 1
    data.pa[data.pa==pseudo] = 0
    data.pa$response = data$response
    data.pa = reshape2::melt(data.pa) %>%
      # DO NOT filter per group
      group_by(variable) %>%
      mutate(prevalence = mean(value)) %>%
      subset(prevalence >= prev) %>%
      mutate(variable = as.character(variable))
    data.filtered = data[,colnames(data) %in% c(unique(data.pa$variable),colnames(metadata.oars.stool.double))]
    dim(data.filtered)
    
    return(data.filtered)
    
  }else{
  data.pa = data[,!colnames(data) %in% colnames(metadata.oars.stool.double)]
  data.pa[data.pa!=pseudo] = 1
  data.pa[data.pa==pseudo] = 0
  data.pa$response = data$response
  data.pa = reshape2::melt(data.pa) %>%
    # FILTER per group
    group_by(response, variable) %>%
    mutate(prevalence = mean(value)) %>%
    subset(prevalence >= prev) %>%
    mutate(variable = as.character(variable))
  data.filtered = data[,colnames(data) %in% c(unique(data.pa$variable),colnames(metadata.oars.stool.double))]
  dim(data.filtered)
  
  return(data.filtered)
  }
}

# this function converts a longitudinal omic table into
# a table of "delta" values; log-fold changes relative to the preceding value

delta.omic = function(data = oars.asv.data.glom.prep) # omics merged with meta
  { 
  
  omic.data = do.call(cbind, lapply(colnames(data)[colnames(data)!="delta.ph"], function(x){
    #x = colnames(data)[colnames(data)!="delta.ph"][1]
    
    data.subset = data.frame(data[,x])
    colnames(data.subset) = "feature"
    data.subset$standard.name = rownames(data)
    data.subset = merge(data.subset, metadata.oars.stool.double[,c("delta.ph","RS_Name", "phase", "HM", "timing","reltiming","diagnosis", "standard.name")],
                        by="standard.name", all.x=T)  
    # remove samples with only 1 timepoint per phase (i.e. non-compliant)
    data.subset = data.subset %>% 
      subset(!is.na(phase)) %>%
      group_by(HM, phase) %>%
      filter(n() > 1) %>% ungroup()
    
   # rownames(data.subset) = data.subset$standard.name
  #  data.subset$standard.name = NULL
    # lm prep
    data.subset$feature = log(data.subset$feature, base=2)
    # calculate delta
    data.subset = data.subset %>%
      group_by(HM, phase) %>%
      arrange(phase) %>%
      mutate(delta.feature = feature[reltiming == "post"] - feature[reltiming == "pre"]) %>% data.frame()
    data.subset = data.subset[,c("delta.feature", "delta.ph","RS_Name", "phase", "HM", "timing","diagnosis", "standard.name", "adj.fiber")] %>% na.omit() %>% distinct()
    colnames(data.subset)[1] = x
    rownames(data.subset) = data.subset$standard.name
    data.subset$standard.name = NULL
    data.subset = data.subset[,c(x,"delta.ph")]
    data.subset
  }))
  omic.data.delta.ph = omic.data$delta.ph
  omic.data[,colnames(omic.data) == "delta.ph"] <- NULL
  omic.data$delta.ph = omic.data.delta.ph
  return(omic.data)
}


# Note: this function is semi-defunct
# this function takes the "delta" data
# and builds an lmer for each feature
# adjusted for HM as random effect (repeat measures)
# as well as the timing (3M or 6M), to adjust for RS sequence

delta.omic.lmer = function(data = oars.asv.data.glom.prep,
                           split = TRUE,
                           mpx = F) {
  
  lm.output = do.call(rbind, lapply((colnames(data)[!colnames(data) %in% c("HM", "timing", "phase","reltiming", "response","diagnosis", "standard.name", "adj.fiber")]), function(x){
    print(x)
    # x = colnames(data)[-length(colnames(data))][3]
    # isolate per taxa and compare to lfc or response
    data.subset = data[,colnames(data) %in% c(x, "HM", "timing", "phase","reltiming", "response","diagnosis", "standard.name", "adj.fiber")] %>% data.frame()
    colnames(data.subset)[1] = "feature"
    data.subset$standard.name = data$standard.name
    data.subset$reltiming = factor(data.subset$reltiming, levels=c("pre", "post"))
    data.subset$response = factor(data.subset$response, levels=c("low", "high"))
    
    # log scale
    data.subset$feature = log2(data.subset$feature)
    # subset based on response level
    if(split == TRUE){
    lmer.output = do.call(rbind, lapply(c("high", "low"), function(res){
       data.response = subset(data.subset, response == res)
      if(sum(data.response[,1]==min(data.response[,1]))<=6){
        #
        data.response$feature = (data.response$feature)
        
        if(mpx == T){
        lmer.output = lmerTest::lmer((feature) ~ reltiming + phase + (1|HM)+(1|plate), data.response) %>% summary() %>% coef()
        lmer.output = data.frame(lmer.output)[2,c(1,5)]
        lmer.output$taxa = x
        colnames(lmer.output)[1] = "estimate"
        colnames(lmer.output)[2] = "pval"
        lmer.output$response = res
        lmer.output
        }else{
          lmer.output = lmerTest::lmer((feature) ~ reltiming + phase + (1|HM), data.response) %>% summary() %>% coef()
          
          lmer.output = data.frame(lmer.output)[2,c(1,5)]
          lmer.output$taxa = x
          colnames(lmer.output)[1] = "estimate"
          colnames(lmer.output)[2] = "pval"
          lmer.output$response = res
          lmer.output
        }
      } else {lmer.output = data.frame(estimate = NA, pval = NA, taxa = x, response = res)}
    }))
    }
    if(split == FALSE){
        if(sum(data.subset[,1]!=0)>6){
          #
          data.subset$feature = (data.subset$feature)
          # conditionally control for baseline abundance (removed: diagnosis+adj.fiber)
          lmer.output = lmerTest::lmer((feature) ~ reltiming*response + phase + (1|HM), data.subset) %>% summary() %>% coef()
          #
          lmer.output = data.frame(lmer.output)[5,c(1,5)]
          lmer.output$taxa = x
          colnames(lmer.output)[1] = "estimate"
          colnames(lmer.output)[2] = "pval"
          lmer.output
        } else {lmer.output = data.frame(estimate = NA, pval = NA, taxa = x, response = res)}
    lmer.output
    }
    return(lmer.output)
  }))
  lm.output$padj = p.adjust(lm.output$pval, method="BH")
  lm.output = lm.output %>% arrange(-pval)
  return(lm.output)
}

omic.gam = function(data = oars.mgx.taxa, # oars.mpx.cog.mat; oars.mgx.taxa; oars.asv.data.glom
                    meta = metadata.oars.stool.asv,
                    transform = "log2",
                    type = "individualized", # | group
                    run.gam = T,
                    prev = 0.10){
  # replace NA with 0
  data[is.na(data)] = 0
  
  # filter to 2+ observations per participant
  # 10% prevalence (or 80% for MBX)
  data.pa = data# *10000 
  data.pa[data.pa != 0] = 1
  data.filt = data[,colSums(data.pa) >= nrow(data)*0.10]
  print(paste(ncol(data.filt), " features remain (out of ", ncol(data), ")", sep =""))
  
  # add pseudo
  data.filt.pseudo = min(data.filt[data.filt!=0])/2
  # log transform

  if(run.gam != TRUE){
    data.filt = log2(data.filt+data.filt.pseudo)
    return(data.filt)
    break
    }else{
  gam.output = do.call(rbind, lapply(colnames(data.filt), function(x){
    print(x)
   # x = colnames(data.filt)[1]
    # select feature and add metadata
    data.subset = data.filt %>% as.data.frame()
    data.subset$standard.name = rownames(data.subset)
    data.subset = data.subset[,colnames(data.subset) %in% c(x, "standard.name")]
    data.subset = merge(data.subset,
                        meta[,c("standard.name", "oars.days", "HM", "compliant")])
    # remove non-compliant samples
    data.subset = subset(data.subset, compliant == T)
    # select individuals with 2+ timepoints
    data.subset = subset(data.subset, HM %in% subset(data.frame(table(data.subset$HM)), Freq >1)$Var1)
    colnames(data.subset)[2] = "feature"
    # fix data types
    data.subset$oars.days = as.numeric(data.subset$oars.days)
    data.subset$HM = as.factor(data.subset$HM)
    
    # transform accordingly
    # ASV, MGX, MPX = log2
    if(transform == "log2"){
      data.subset$feature = log2( data.subset$feature + 1)
    }
    # MBX = sqrt
    if(transform == "sqrt"){
      data.subset$feature = sqrt(data.subset$feature)
    }
    
    # For group-level analyses
    if(type == "group"){
    if(length(unique(subset(data.subset, feature != min(data.subset$feature))$HM))>3){
      
    # perform GAM
      gam.output = mgcv::gam(((feature)) ~ s((oars.days), k=10) + s(HM, bs="re"),
                data = data.subset,
                family = mgcv::Tweedie(p=1.01),
                method = "REML")# %>% summary() 
      # calculate derivative per phase
      smooth_data <- gratia::derivatives(gam.output, term = "s(oars.days)", order=1, type = "central", n = 100)      
      smooth_data$phase = ifelse(smooth_data$oars.days < max(subset(meta, oars.on.rs %in% c("preRS", "onRS"))$oars.days), "treatment", "washout")
      # calculate mean of first derivative per phase (positive = /; negative = \)
      treatment.curvature = mean(subset(smooth_data, phase == "treatment")$`.derivative`)
      washout.curvature = mean(subset(smooth_data, phase == "washout")$`.derivative`)
      
     # extract output
      gam.output.df = summary(gam.output)$s.table[,4] %>% as.data.frame()
      gam.output.df$variable = c("days", "HM")
      colnames(gam.output.df)[1] = "pval"
      gam.output.df$feature = x
      # add curvature
      gam.output.df$treatment =  treatment.curvature
      gam.output.df$washout =  washout.curvature
      # print
      print(gam.output.df)
    }else{
      data.frame(pval = c(NA, NA),
                 variable = c("days", "HM"),
                 feature = x,
                 treatment = NA,
                 washout = NA)
    }
    }
    # for individual level analysis
    if(type == "individualized"){
      if(length(unique(subset(data.subset, feature != min(data.subset$feature))$HM))>3){
        
        # perform GAM
        gam.output = mgcv::gam(((feature)) ~ s((oars.days), k=3) + s(oars.days, by=HM, bs="fs",m=1),
                               data = data.subset,
                               #family = mgcv::Tweedie(p=1.01),
                               method = "REML")# %>% summary() 
        #summary(gam.output) %>% as.data.frame()
        
        # extract
        gam.output.df = data.frame(edf = summary(gam.output)$edf,
                   HM = c("global", unique(as.character(data.subset$HM))))
        gam.output.df = gam.output.df %>%
          mutate(coef = edf - edf[HM == "global"]) %>%
          subset(HM != "global") %>%
          mutate(feature = x)
        
       # print
        print(gam.output.df)
      }else{
        data.frame(edf = NA, HM = NA, coef = NA, feature = x)
      }
    }
  }))
    }
}
    
omic.gam.plot = function(f = feature,
                         yadjust=yadjust){
  
  # identify which data set the feature is in
  data.subset = subset(oars.omics.gam, feature == f)
  omic.type = unique(data.subset$data.type)
  
  # select dataset
  if(omic.type == "ASV"){
    feature.data = oars.asv.data.glom %>% as.data.frame()
  }
  if(omic.type == "Species"){
    feature.data = oars.mgx.taxa %>% as.data.frame() 
  }
  if(omic.type == "Pathway"){
    feature.data = oars.mpx.kegg.mat %>% as.data.frame() 
  }
  if(omic.type == "COG"){
    feature.data = oars.mpx.cog.mat %>% as.data.frame()
  }
  if(omic.type == "CAZy"){
    feature.data = oars.mpx.cazy.mat %>% as.data.frame()
  }
  if(omic.type == "Metabolite"){
    feature.data = oars.mbx.raw.mat.filt 
  }
  # establish pseudo
  feature.data[is.na(feature.data)] = 0
  pseudo = min(feature.data[feature.data!=0])/2
  
  # add standard.name
  feature.data = feature.data %>%
    mutate(standard.name = rownames(.))
  
  # select feature
  feature.data = feature.data[,colnames(feature.data) %in% c("standard.name", f)]
  colnames(feature.data)[1] = "feature"
  # add metadata
  feature.data = merge(feature.data,
                       metadata.oars.stool.asv, by="standard.name")
  
  # plot
  feature.data.plot = ggplot(subset(feature.data, compliant == TRUE),
                                    aes(x=oars.days, 
                                        y=feature+pseudo))+
    annotate("rect", xmin=0, xmax=max(subset(metadata.oars.stool, oars.on.rs=="onRS")$oars.days),
             ymin=0, ymax=Inf, alpha=0.2, fill="salmon") +
    scale_y_log10()+
    geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
    geom_point(aes(fill=rs.col, shape=baseline), size=2)+
    scale_shape_manual(values=c(23,21))+
    scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
    geom_smooth(color="black")+
    theme_classic()+theme(legend.position="none",
                          plot.title = element_text(hjust = 0.5, size=12))+
    # add GAM pvalue
    geom_text(data=pval.positioner("feature",type="omics",
                                   feature.data,
                                   data.subset$padj,
                                   yadjust = 1.01), aes(x=xval, y=yval, label=label), vjust=1,size=4)+
    labs(x="Days since starting RS", y=f)
  feature.data.plot
  
}

    

# this function positions pvalues in a stable location for certain plots
pval.positioner = function(
    variable = "load.mgx",
    data = metadata.oars.stool.asv,
    data1 = stats.fcal.036,
    data2 = stats.fcal.6912,
    type = "regular", # or "omics
    yadjust = 1.1){
  if(is.null(variable)){
    print("Enter variable")
    break
  }
  # 1. calculate x ranges for each phase:
  rs.range = c(0, max(subset(data, compliant == TRUE & oars.on.rs == "onRS")$oars.days))
  wo.range = c(max(subset(data, compliant == TRUE & oars.on.rs == "onRS")$oars.days), max(subset(data, compliant == TRUE & oars.on.rs == "postRS")$oars.days))
  # 2. calculate center value (mean)
  rs.center = mean(rs.range)
  wo.center = mean(wo.range)
  # 3. calculate y value, add 5%
  y.val = max(na.omit(data[,variable]))*yadjust
  if(type == "regular"){
  # 4. save coords as data.frame
  data.frame(phase = c("RS", "Washout"),
             xval = c(rs.center, wo.center),
             yval = y.val,
             label = c(paste("p:", round(data.frame(data1)[2,5], digits=3)),
                       paste("p:", round(data.frame(data2)[2,5], digits=3))))
  }else{
    data.frame(phase = c("RS", "Washout"),
               xval = c(rs.center, wo.center),
               yval = y.val,
               label = c(paste("Time\np:", round(data.frame(data1)[1,], digits=3)),
                         paste("Patient\np:", round(data.frame(data1)[2,], digits=3))))
  }
}

# this function calculates within/between sample beta-diversity
# into a format compatible with standard variable analyses
beta.trajectory = function(data=oars.asv.bray){
  # melt
  data.subset = reshape2::melt(as.matrix(data))
  # note, I'm keeping the duplicate values (eg. X-Y and Y-X)
  # because attempting to remove them could introduce new errors (hard to pin down)
  # and the average will be the same if all values are duplicated
  
  # but we will remove self comparisons
  data.subset = subset(data.subset, Var1 != Var2)
  
  # now calculate between sample beta
  
  # unpack sample names
  data.subset = tidyr::separate(data.subset, col=Var1, into=c("HM.A", "STL.A", "Number.A"), sep="-", remove=F)
  data.subset = tidyr::separate(data.subset, col=Var2, into=c("HM.B", "STL.B", "Number.B"), sep="-", remove=F)
  # calculate between sample diversity
  between.beta = data.subset %>% 
    group_by(Var1) %>%
    mutate(between.beta = mean(value)) %>%
    dplyr::select(Var1, between.beta) %>% distinct()
  # calculate within patient diversity
  within.beta = data.subset %>% 
    subset(HM.A == HM.B) %>%
    group_by(Var1) %>%
    mutate(within.beta = mean(value)) %>%
    dplyr::select(Var1, within.beta) %>% distinct()
  # merge
  beta.data = merge(between.beta,
                    within.beta, by="Var1")
  colnames(beta.data)[1] = "standard.name"
  return(beta.data)
}

# function to make volcano plots have an intelligible y axis
neg_log10_trans <- scales::new_transform(
  name = "neglog10",
  transform = function(x) -log10(x),
  inverse = function(x) 10^(-x),
  format = function(x) format(x, scientific = FALSE)
)


# function to perform RF and associated analyses
rf.function = function(data.types = "All",
                       iters = 15,
                       target = "fermentation", # fermentation | fcal
                       output = "AUC"){ # | ROC | importances | interactions

  # start timer
  t1 = Sys.time()
  
  rf.function.output = 
    
    # select data.type
    do.call(rbind, lapply(data.types, function(data.type){
     
       # select iteration (seed)
      do.call(rbind, lapply(1:iters, function(iter){
        
        # track progress
        print(paste(data.type, iter))
        
        # select data
        if(data.type == "ASV"){
          data.all = oars.asv.baseline
        }
        if(data.type == "Species"){
          data.all = oars.mgx.baseline
        }
        if(data.type == "Pathway"){
          data.all = oars.mpx.kegg.baseline
        }
        if(data.type == "COG"){
          data.all = oars.mpx.cog.baseline
        }
        if(data.type == "CAZy"){
          data.all = oars.mpx.cazy.baseline
        }
        if(data.type == "Metabolite"){
          data.all = oars.mbx.baseline
        }
        if(data.type == "FFQ"){
          data.all = oars.ffq.baseline
        }
        if(data.type == "All"){
          data.all = oars.all.omics.baseline
        }
        
        # subset to samples with matching data across omics (for apples-to-apples comparison across omics)
        rownames(data.all) = paste(data.all$HM, 
                                   data.all$phase, sep="_")
        data.all = data.all[rownames(data.all) %in% rownames(oars.all.omics.baseline),]
        
        # optionally target fcal response instead of fermentation response
        # note, predictions are lower than AUC 0.65, so consider this defunct
        if(target == "fcal"){
          data.all = merge(data.all,
                           oars.fcal.response, by="row.names") %>% as.data.frame()
          rownames(data.all) = data.all$Row.names
          data.all$Row.names = NULL
          data.all$response = data.all$fcal.response
          data.all$fcal.response = NULL
        }
        
        # keep track of HMs and indices (for LOOCV and class balancing)
        hm.list = data.frame(HM = data.all$HM,
                             index = seq(from=1,to=length(data.all$HM)))
        
        # remove HM column
        data.all = data.all[,colnames(data.all) != "HM"]
        
        # factorize features for rfsrc
        data.all = data.all %>%
          mutate(response = as.factor(response),
                 diagnosis = as.factor(diagnosis),
                 phase = as.factor(phase))
          
          
        # for "output == interactions", subset features to important
        if(output == "interactions"){
          data.all = data.all[,colnames(data.all) %in% c("response", oars.rf.loocv.imp.df$feature)]
        }

        # loop through microbiomes (LOOCV); run in parallel
        oars.loocv = do.call(rbind, parallel::mclapply(1:nrow(data.all), function(hm){
          
          print(paste(iter, hm))
          
          # remove potential sample pseudoreplicates
          # (avoids potential overfitting, but cuts n by 1)
          data.all$HM = hm.list$HM
          
          # reduce to train/test
          data.train = data.all[-hm, ]
          # remove other samples from the same HM
          data.train = data.train %>%
            subset(HM != subset(hm.list, index == hm)$HM)
          
          data.train$HM = NULL
          
          # balance data (downsample majority class)
          set.seed(iter)
          data.train = data.train %>%
            group_by(response) %>%
            sample_n(min(table(data.train$response)),replace=F)
          
          # now that we've removed samples, we should remove features with no variance
          no_var_cols = apply(data.train[,!colnames(data.train) %in% c("phase", "response","diagnosis","adj.fiber")], 2, sd) %>% data.frame() %>% subset(. == 0)
          data.train = data.train[,!colnames(data.train) %in% rownames(no_var_cols)]
          
          # additionally, we can identify features with low (non-0) variance
          low_var_cols <- caret::nearZeroVar(data.train[,!colnames(data.train) %in% c("phase", "response","diagnosis","adj.fiber")], 
                                             saveMetrics = TRUE)
          
          # select features with >= 20% unique values to keep consistent with other filtering
          data.train = data.train[,!colnames(data.train) %in% 
                                    rownames(subset(low_var_cols, percentUnique <20))]

          # build RF model
          # use rfsrc so we can extract interaction scores later
          set.seed(iter)
          model.results = randomForestSRC::rfsrc(response~ ., data.train %>%
                                                   as.data.frame(), 
                                                 # importance="permute",
                                                 probability = TRUE)
          
          ## OUTPUT 1 = AUC | ROC (performances of model(s)); continues outside of loop
          if(output %in% c("AUC", "ROC")){
            pred = predict(model.results, 
                           data.all[hm, ])
            
            output = data.frame(pred = pred$predicted[1], # extract prob of "high"
                                high = factor(data.all[hm, ]$response, levels=c("high", "low")), # if "high", true
                                index = hm,
                                iter = iter)
            return(output)
          }
          
          ## OUTPUT 2 = Importances
          if(output == "importances"){
            # calculate importances
            set.seed(iter)
            ints = randomForestSRC::vimp(model.results, method = "permute", joint=F) 
            ints = data.frame(imp = ints$importance[,1]) %>%
              arrange(imp) %>%
              mutate(iter = iter,
                     index = hm) %>%
              mutate(feature = rownames(.))
            rownames(ints) = NULL
            return(ints)
          }
          
          ## OUTPUT 3 = Interactions
          if(output == "interactions"){
            # calculate interactions
            set.seed(iter)
            ints = randomForestSRC::find.interaction(model.results,
                                                     method = "vimp", verbose=F)
            ints.df = ints %>% 
              as.data.frame() %>%
              mutate(features = rownames(.))%>%
              tidyr::separate(features, into=c("var1", "var2"), sep=":", remove=F)%>%
              # save df identifies
              mutate(
                iter = iter,
                hm = hm)
            return(ints.df)
          }
          
        })) # end loop for a single sample
        # for "output==AUC", calculate AUC of that iteration
        if(output == "AUC"){
        # calculate LOOCV AUC
        oars.loocv.auc = oars.loocv %>%
          mutate(auc = pROC::auc(high, pred, 
                                 levels=c("high", "low"),  # define case = "high" delta pH
                                 direction="<")[1], # if probability > 0.5 (e.g.), then it is a case
                 iter = iter,
                 data.type = data.type) %>% suppressWarnings() %>% suppressMessages()
        return(oars.loocv.auc)
        }else{
          return(oars.loocv)
        }
      }))
    }))
  
  print(Sys.time() - t1)
  
  return(rf.function.output)
}
# stop timer outside


# ::``` Stool Plot -----------------------------------------------------------

# n = 50 compliant samples
subset(metadata.oars.stool, compliant == TRUE) %>% nrow()

# visualize stools + timings
metadata.oars.stool.plot = ggplot(metadata.oars.stool %>%
                                    mutate(missing.mpx = ifelse(standard.name %in% rownames(oars.mpx.cog.mat), "","*")) %>%
                                    dplyr::select(HM, oars.days, compliant, oars.on.rs) %>% distinct(),
                                  aes(x=oars.days, y=reorder(HM, compliant)))+
  annotate("rect", xmin=0, xmax=max(subset(metadata.oars.stool, oars.on.rs=="onRS")$oars.days),
           ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_point(aes(shape=oars.on.rs), fill="white",  size=3)+
  geom_point(aes(shape=oars.on.rs, fill=oars.on.rs,
                 alpha = ifelse(compliant == T, "1", "0")),
             size=3)+
  scale_fill_manual(values=c("grey", 2, "grey"))+
  scale_shape_manual(values=c(23, 21,21))+
  scale_alpha_manual(values=c(0.25,1))+
  guides(alpha= "none")+
  #geom_text(aes(label=missing.asv), nudge_y=-0.15, size=6)+
  #geom_text(aes(label=substr(standard.name, nchar(standard.name)-2, nchar(standard.name))), size=2)+
  theme_classic()+theme(panel.grid.major.y=element_line(color="grey", linewidth=0.2),
                        legend.position="left")+
  labs(x="Day since starting RS", y="", shape="Phase", fill="Phase", alpha="")
metadata.oars.stool.plot

# >> 1. STANDARD ANALYSES ----------------------------------------------------

# e.g. fecal calprotectin, stool water, richness, etc
# plus simple ASV analyses (e.g. alpha, beta, functional redundancy, )

# remove non-compliant stools

# fit lmer for each phase, using 0M for treatment baseline, and 6M as washout baseline

# "do features change over treatment, and do they change again over washout"

# :: Fecal Calprotectin ----------------------------------------------------------

# does fecal calprotectin increase or decrease over RS treatment

# stats
stats.fcal.036 = lmerTest::lmer(scale(log10(fcal)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
               subset(metadata.oars.stool, timing %in% c("0M", "3M", "6M") & compliant == TRUE)) %>%
  summary() %>% coef()

stats.fcal.6912 = lmerTest::lmer(scale(log10(fcal)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
               subset(metadata.oars.stool, timing %in% c("6M", "9M", "12M") & compliant == TRUE)) %>%
  summary()  %>% coef()

# GAM
stats.fcal.gam = mgcv::gam(scale(log10(fcal)) ~ s(scale(oars.days), k=10) + diagnosis + adj.fiber + s(HM, bs="re"),
                                data = subset(metadata.oars.stool, compliant == TRUE & !is.na(fcal))%>%
                             mutate(HM = as.factor(HM)),
                                family = gaussian(),
                                method = "REML") %>% summary() 

# plot
metadata.oars.fcal.plot = ggplot(subset(metadata.oars.stool, compliant == TRUE),
       aes(x=oars.days, 
           y=fcal))+
  scale_y_log10()+ # need to put before, so max=Inf works
  annotate("rect", xmin=0, xmax=max(subset(metadata.oars.stool, oars.on.rs=="onRS")$oars.days),
           ymin=0, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(aes(fill=rs.col, shape=baseline), size=3)+
  scale_shape_manual(values=c(23,21))+
  geom_smooth(color="black")+
  geom_hline(yintercept=250, linetype=1, alpha=1, color="red")+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12))+
  geom_text(data=pval.positioner("fcal",
                                        subset(metadata.oars.stool, compliant == TRUE),
                                        stats.fcal.036,
                                        stats.fcal.6912,
                                 yadjust=2), aes(x=xval, y=yval, label=label), vjust=1,size=4)+
  labs(x="Days since starting RS", y="Fecal Calprotectin (μg/g)")

metadata.oars.fcal.plot


# :: Stool Water ----------------------------------------------------------

# does stool water increase or decrease over RS treatment

# stats
stats.water.036 = lmerTest::lmer(scale(log10(stool_water_perc)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                                subset(metadata.oars.stool, timing %in% c("0M", "3M", "6M") & compliant == TRUE)) %>%
  summary() %>% coef()

stats.water.6912 = lmerTest::lmer(scale(log10(stool_water_perc)) ~ scale(oars.days) + diagnosis +  adj.fiber +(1|HM),
                                 subset(metadata.oars.stool, timing %in% c("6M", "9M", "12M") & compliant == TRUE)) %>%
  summary() %>% coef()

# GAM
stats.water.gam = mgcv::gam(scale(log10(stool_water_perc)) ~ s(scale(oars.days), k=10) + diagnosis + adj.fiber + s(HM, bs="re"),
                           data = subset(metadata.oars.stool, compliant == TRUE & !is.na(stool_water_perc))%>%
                             mutate(HM = as.factor(HM)),
                           family = gaussian(),
                           method = "REML") %>% summary() 

# plot
metadata.oars.water.plot = ggplot(subset(metadata.oars.stool, compliant == TRUE),
       aes(x=oars.days, 
           y=stool_water_perc*100))+
  annotate("rect", xmin=0, xmax=max(subset(metadata.oars.stool, oars.on.rs=="onRS")$oars.days),
           ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(aes(fill=rs.col, shape=baseline), size=2)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  geom_smooth(color="black")+
 theme_classic()+theme(legend.position="none",
                       plot.title = element_text(hjust = 0.5, size=12))+
  geom_text(data=pval.positioner("stool_water_perc",
                                 subset(metadata.oars.stool, compliant == TRUE),
                                 stats.water.036,
                                 stats.water.6912,
                                 yadjust = 1.03), aes(x=xval, y=yval*100, label=label), vjust=1,size=4)+
  labs(x="Days since starting RS", y="Stool Moisture (%)")
metadata.oars.water.plot


# :: Microbial Load (ASV) -------------------------------------------------------
# Estimated; predictor had R = 0.79

# stats
stats.load.asv.036 = lmerTest::lmer(scale(load.asv) ~ scale(oars.days) + diagnosis +  adj.fiber +(1|HM),
                                    subset(metadata.oars.stool.asv, timing %in% c("0M", "3M", "6M") & compliant == TRUE)) %>%
  summary() %>% coef()

stats.load.asv.6912 = lmerTest::lmer(scale(load.asv) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                                     subset(metadata.oars.stool.asv, timing %in% c("6M", "9M", "12M") & compliant == TRUE)) %>%
  summary() %>% coef()

# GAM
stats.load.asv.gam = mgcv::gam(scale((load.asv)) ~ s(scale(oars.days), k=10) + diagnosis + adj.fiber + s(HM, bs="re"),
                           data = subset(metadata.oars.stool.asv, compliant == TRUE & !is.na(load.asv))%>%
                             mutate(HM = as.factor(HM)),
                           family = gaussian(),
                           method = "REML") %>% summary() 

# plot
metadata.oars.load.asv.plot = ggplot(subset(metadata.oars.stool.asv, compliant == TRUE),
                                  aes(x=oars.days,
                                      y=load.asv))+
  scale_y_log10()+
  annotate("rect", xmin=0, xmax=max(subset(metadata.oars.stool, oars.on.rs=="onRS")$oars.days),
           ymin=0, ymax=Inf,
           alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(aes(fill=rs.col, shape=baseline), size=2)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  geom_smooth(color="black")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12))+
  geom_text(data=pval.positioner("load.asv",
                                 subset(metadata.oars.stool.asv, compliant == TRUE),
                                 stats.load.asv.036,
                                 stats.load.asv.6912,
                                 yadjust = 1.01), aes(x=xval, y=yval, label=label), vjust=1,size=4)+
  labs(x="Days since starting RS", y="Microbial Load (Log10 Predicted)")
metadata.oars.load.asv.plot


## Do Water and Load correlate? (note: no need to remove non-compliant)
oars.stool.water.load.correlation.plot = ggplot(metadata.oars.stool.asv,
                                                aes(x=stool_water_perc*100, y=load.asv))+
  geom_point(shape=21)+
  geom_smooth(method="lm", color="black")+
  ggpubr::stat_cor(method="spearman", size=4,
                   label.x = 70,
                   label.y=11)+
  theme_classic()+
  theme(strip.text=element_text(size=10))+
  facet_wrap(~"Predicted Bacterial Load vs Stool Moisture")+
  labs(x="Stool Moisture (%)",
       y="Microbial Load (Log10 Predicted)")
oars.stool.water.load.correlation.plot
# Bacterial load decreases as stool water increases
# Note: using ASVs

# :: Microbial Load (MGX) -------------------------------------------------------
# Estimated; predictor had R = 0.54

# stats
stats.load.mgx.036 = lmerTest::lmer(scale(load.mgx) ~ scale(oars.days) + diagnosis +  adj.fiber +(1|HM),
                                    subset(metadata.oars.stool.asv, timing %in% c("0M", "3M", "6M")& compliant == TRUE)) %>%
  summary() %>% coef()

stats.load.mgx.6912 = lmerTest::lmer(scale(load.mgx) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                                     subset(metadata.oars.stool.asv, timing %in% c("6M", "9M", "12M")& compliant == TRUE)) %>%
  summary() %>% coef()

# GAM
stats.load.mgx.gam = mgcv::gam(scale(load.mgx) ~ s(scale(oars.days), k=10) + diagnosis + adj.fiber + s(HM, bs="re"),
                           data = subset(metadata.oars.stool.asv, compliant == TRUE & !is.na(load.mgx))%>%
                             mutate(HM = as.factor(HM)),
                           family = gaussian(),
                           method = "REML") %>% summary() 

# plot
metadata.oars.load.mgx.plot = ggplot(subset(metadata.oars.stool.asv, compliant == TRUE),
                                     aes(x=oars.days,
                                         y=load.mgx))+
  scale_y_log10()+
  annotate("rect", xmin=0, xmax=max(subset(metadata.oars.stool, oars.on.rs=="onRS")$oars.days),
           ymin=0, ymax=Inf,
           alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(aes(fill=rs.col, shape=baseline), size=2)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  geom_smooth(color="black")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12))+
  geom_text(data=pval.positioner("load.mgx",
                                 subset(metadata.oars.stool.asv, compliant == TRUE),
                                 stats.load.mgx.036,
                                 stats.load.mgx.6912,
                                 yadjust = 1.01), aes(x=xval, y=yval, label=label), vjust=1,size=4)+
  labs(x="Days since starting RS", y="Microbial Load (Log10 Predicted)")
metadata.oars.load.mgx.plot



# :: ASV Richness ------------------------------------------------------

# stats
stats.richness.036 = lmerTest::lmer(scale(log10(richness)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                                 subset(metadata.oars.stool.asv, timing %in% c("0M", "3M", "6M") & compliant == TRUE)) %>%
  summary() %>% coef()

stats.richness.6912 = lmerTest::lmer(scale(log10(richness)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                                  subset(metadata.oars.stool.asv, timing %in% c("6M", "9M", "12M") & compliant == TRUE)) %>%
  summary() %>% coef()

# GAM
stats.richness.gam = mgcv::gam(scale(log10(richness)) ~ s(scale(oars.days), k=10) + diagnosis + adj.fiber + s(HM, bs="re"),
                           data = subset(metadata.oars.stool.asv, compliant == TRUE & !is.na(richness))%>%
                             mutate(HM = as.factor(HM)),
                           family = gaussian(),
                           method = "REML") %>% summary() 


# plot
metadata.oars.stool.asv.richness.plot = ggplot(subset(metadata.oars.stool.asv,compliant == TRUE),
       aes(x=oars.days, 
           y=richness))+
  #scale_y_log10()+
  annotate("rect", xmin=0, xmax=max(subset(metadata.oars.stool, oars.on.rs=="onRS")$oars.days),
           ymin=-Inf, ymax=Inf,
           alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(aes(fill=rs.col, shape=baseline), size=2)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  geom_smooth(color="black")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=10))+
  geom_text(data=pval.positioner("richness",
                                 subset(metadata.oars.stool.asv, compliant == TRUE),
                                 stats.richness.036,
                                 stats.richness.6912,
                                 yadjust = 1), aes(x=xval, y=yval, label=label), vjust=1,size=4)+
  labs(x="Days since starting RS", y="ASV Richness")
metadata.oars.stool.asv.richness.plot


# :: ASV Shannon ------------------------------------------------------

# stats
stats.shannon.036 = lmerTest::lmer(scale((shannon)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                                    subset(metadata.oars.stool.asv, timing %in% c("0M", "3M", "6M") & 
                                             compliant == TRUE)) %>%
  summary() %>% coef()

stats.shannon.6912 = lmerTest::lmer(scale((shannon)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                                     subset(metadata.oars.stool.asv, timing %in% c("6M", "9M", "12M") & 
                                              compliant == TRUE)) %>%
  summary() %>% coef()

# GAM
stats.shannon.gam = mgcv::gam(scale(shannon) ~ s(scale(oars.days), k=10) + diagnosis + adj.fiber + s(HM, bs="re"),
                           data = subset(metadata.oars.stool.asv, compliant == TRUE & !is.na(shannon))%>%
                             mutate(HM = as.factor(HM)),
                           family = gaussian(),
                           method = "REML") %>% summary() 

# plot
metadata.oars.stool.asv.shannon.plot = ggplot(subset(metadata.oars.stool.asv,compliant==TRUE),
       aes(x=oars.days, 
           y=shannon))+
  annotate("rect", xmin=0, xmax=max(subset(metadata.oars.stool, oars.on.rs=="onRS")$oars.days),
           ymin=-Inf, ymax=Inf,
           alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(aes(fill=rs.col, shape=baseline), size=2)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  geom_smooth(color="black")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12))+
  geom_text(data=pval.positioner("shannon",
                                 subset(metadata.oars.stool.asv, compliant == TRUE),
                                 stats.shannon.036,
                                 stats.shannon.6912,
                                 yadjust = 1.05), aes(x=xval, y=yval, label=label), vjust=1,size=4)+
  labs(x="Days since starting RS", y="Shannon Diversity")
metadata.oars.stool.asv.shannon.plot


# :: ASV Functional Redundancy --------------------------------------------

# no need to log-transform

# stats
stats.fd.036 = lmerTest::lmer(scale((fd)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                                    subset(metadata.oars.stool.asv, timing %in% c("0M", "3M", "6M")& compliant == TRUE)) %>%
  summary() %>% coef()

stats.fd.6912 = lmerTest::lmer(scale((fd)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                                     subset(metadata.oars.stool.asv, timing %in% c("6M", "9M", "12M")& compliant == TRUE)) %>%
  summary() %>% coef()

# GAM
stats.fd.gam = mgcv::gam(scale((fd)) ~ s(scale(oars.days), k=10) + diagnosis + adj.fiber + s(HM, bs="re"),
                           data = subset(metadata.oars.stool.asv, compliant == TRUE & !is.na(fd))%>%
                             mutate(HM = as.factor(HM)),
                           family = gaussian(),
                           method = "REML") %>% summary() 

# plot
metadata.oars.stool.asv.fd.plot = ggplot(subset(metadata.oars.stool.asv,compliant == TRUE),
                                               aes(x=oars.days, 
                                                   y=fd))+
  annotate("rect", xmin=0, xmax=max(subset(metadata.oars.stool.asv, oars.on.rs=="onRS")$oars.days),
           ymin=-Inf, ymax=Inf,
           alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(aes(fill=rs.col, shape=baseline), size=2)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  geom_smooth(color="black")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=10))+
  geom_text(data=pval.positioner("fd",
                                 subset(metadata.oars.stool.asv, compliant == TRUE),
                                 stats.fd.036,
                                 stats.fd.6912,
                                 yadjust = 1.05), aes(x=xval, y=yval, label=label), vjust=1,size=4)+
  labs(x="Days since starting RS", y="Functional Redundancy")
metadata.oars.stool.asv.fd.plot


# :: ASV Butyrogens I -------------------------------------------------------

# stats
stats.but.i.036 = lmerTest::lmer(scale(log10(but.i)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                                 subset(metadata.oars.stool.asv, timing %in% c("0M", "3M", "6M") & compliant == TRUE)) %>%
  summary() %>% coef()

stats.but.i.6912 = lmerTest::lmer(scale(log10(but.i)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                                  subset(metadata.oars.stool.asv, timing %in% c("6M", "9M", "12M") & compliant == TRUE)) %>%
  summary() %>% coef()

# GAM
stats.but.i.gam = mgcv::gam(scale(log10(but.i)) ~ s(scale(oars.days), k=10) + diagnosis + adj.fiber + s(HM, bs="re"),
                           data = subset(metadata.oars.stool.asv, compliant == TRUE & !is.na(but.i))%>%
                             mutate(HM = as.factor(HM)),
                           family = gaussian(),
                           method = "REML") %>% summary() 


# plot
metadata.oars.but.i.plot = ggplot(subset(metadata.oars.stool.asv, compliant == TRUE),
                                  aes(x=oars.days, 
                                      y=but.i*100))+
  annotate("rect", xmin=0, xmax=max(subset(metadata.oars.stool, oars.on.rs=="onRS")$oars.days),
           ymin=-Inf, ymax=Inf,
           alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(aes(fill=rs.col, shape=baseline), size=3)+
  geom_smooth(color="black")+
  # scale_y_log10()+
  scale_shape_manual(values=c(23,21))+
  #geom_text(aes(label=standard.name), size=3)+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12))+
  geom_text(data=pval.positioner("but.i",
                                 subset(metadata.oars.stool.asv, compliant == TRUE),
                                 stats.but.i.036,
                                 stats.but.i.6912,
                                 yadjust = 1.05), aes(x=xval, y=yval*100, label=label), vjust=1,size=4)+
  labs(x="Days since starting RS", y="Butyrogens (%)")
metadata.oars.but.i.plot


# :: ASV Butyrogens II (Vital) -------------------------------------------------------

# stats
stats.but.ii.036 = lmerTest::lmer(scale(log10(but.ii)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                                  subset(metadata.oars.stool.asv, timing %in% c("0M", "3M", "6M") & compliant == TRUE)) %>%
  summary() %>% coef()

stats.but.ii.6912 = lmerTest::lmer(scale(log10(but.ii)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                                   subset(metadata.oars.stool.asv, timing %in% c("6M", "9M", "12M") & compliant == TRUE)) %>%
  summary() %>% coef()

# GAM
stats.but.ii.gam = mgcv::gam(scale(log10(but.ii)) ~ s(scale(oars.days), k=10) + diagnosis + adj.fiber + s(HM, bs="re"),
                           data = subset(metadata.oars.stool.asv, compliant == TRUE & !is.na(but.ii))%>%
                             mutate(HM = as.factor(HM)),
                           family = gaussian(),
                           method = "REML") %>% summary() 


# plot
metadata.oars.but.ii.plot = ggplot(subset(metadata.oars.stool.asv, compliant == TRUE),
                                   aes(x=oars.days, 
                                       y=but.ii))+
 # scale_y_log10()+
  annotate("rect", xmin=0, xmax=max(subset(metadata.oars.stool, oars.on.rs=="onRS")$oars.days),
           ymin=-Inf, ymax=Inf,
           alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(aes(fill=rs.col, shape=baseline), size=2)+
  scale_shape_manual(values=c(23,21))+
  geom_smooth(color="black")+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12))+
  geom_text(data=pval.positioner("but.ii",
                                 subset(metadata.oars.stool.asv, compliant == TRUE),
                                 stats.but.ii.036,
                                 stats.but.ii.6912,
                                 yadjust = 1.01), aes(x=xval, y=yval, label=label), vjust=1,size=4)+
  labs(x="Days since starting RS", y="Butyrogens (%)")
metadata.oars.but.ii.plot

## correlate
lmerTest::lmer(but.ii ~ but.i + (1|HM), metadata.oars.stool.asv) %>% summary()

# no need to remove non-compliant
metadata.oars.stool.asv.plot = ggplot(metadata.oars.stool.asv,
                                      aes(x=but.i*100, y=but.ii*100))+
  geom_point(shape=21)+
 # scale_x_log10()+
  #scale_y_log10()+
  geom_smooth(method="lm", color="black")+
  ggpubr::stat_cor(method="spearman")+
  theme_classic()+theme(strip.text=element_text(size=10))+
  facet_wrap(~"Butyrogen Comparison")+
  labs(x="Original Butyrogens (%)", y="Kircher Butyrogens (%)")
metadata.oars.stool.asv.plot


# :: ASV Bray-Curtis ------------------------------------------------------

# using: rarefied ASV (ASV-level annotations, annotated, with "make.unique" to give # to shared annotation)
# median value of technical replicates
# permanova sets HM to strata (to control for pseudoreplicates)
# and control for adj.fiber and diagnosis

# remove non-compliant (call this ".1" to indicate `compliant`)
dim(oars.asv.data.median)
oars.asv.data.median.1 = oars.asv.data.median[!rownames(oars.asv.data.median) %in% non.compliant,]
dim(oars.asv.data.median.1)

# filter to 10% prevalence
dim(oars.asv.data.median.1) # 10843 taxa
oars.asv.data.median.1.filt = oars.asv.data.median.1
oars.asv.data.median.1.filt[oars.asv.data.median.1.filt!=0] = 1
oars.asv.data.median.1.filt = oars.asv.data.median.1[,colSums(oars.asv.data.median.1.filt) >= nrow(oars.asv.data.median.1)*0.1]
dim(oars.asv.data.median.1.filt) # 1165


# calculate Bray-Curtis dissimilarities
oars.asv.bray = vegan::vegdist(oars.asv.data.median.1.filt, method="bray") 
# perform PCoA
oars.asv.pcoa = ape::pcoa(oars.asv.bray)
# extract data from pcoa
oars.asv.pcoa.df = data.frame(oars.asv.pcoa$vectors[,c(1:2)])
oars.asv.pcoa.df$standard.name = rownames(oars.asv.pcoa.df)
# add metadata
oars.asv.pcoa.df = merge(oars.asv.pcoa.df,
                         metadata.oars.stool.asv, by="standard.name")
# extract variance explained
oars.asv.pcoa.var_exp = oars.asv.pcoa$values[c(1:2),2]
oars.asv.pcoa.df$var1 = round(oars.asv.pcoa.var_exp[1]*100, digits=2)
oars.asv.pcoa.df$var2 = round(oars.asv.pcoa.var_exp[2]*100, digits=2)

# clean up "oars.on.rs" variable
set.seed(25)
t1 = Sys.time()
oars.asv.permanova = vegan::adonis2(oars.asv.bray ~ on.rs + adj.fiber + HM,
                                                oars.asv.pcoa.df %>% mutate(on.rs = as.factor(ifelse(oars.on.rs=="onRS", "onRS", "offRS"))),
                                                 strata = oars.asv.pcoa.df$HM,
                                    # or, instead of strata:
                                    # with(oars.mgx.pcoa.df, permute::how(nperm = 999, blocks = HM)), 
                                    # from Gavin Simpson: https://github.com/vegandevs/vegan/discussions/600
                                    # note: these give the same results
                                                 by="margin")
t2 = Sys.time()
t2 - t1
oars.asv.permanova # not sig

oars.asv.pcoa.plot <- ggplot(
  data=oars.asv.pcoa.df %>% group_by(HM) %>% arrange(stool_date_rec_v2), 
  aes(x=Axis.1, y=Axis.2))+
  geom_path(aes(group=HM), color="black", linetype=2, alpha=0.5, linewidth=0.3) + 
  geom_point(aes(fill=rs.col, shape=baseline), size=2)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  #geom_text(aes(label=RS_Name))+
  theme_classic()+theme(legend.position="none",
                       plot.title = element_text(hjust = 0.5, size=12),
                       strip.text=element_text(size=10))+
  #annotate(geom="text",
  #         x=-0.27,
  #         y=-0.25,
  #         label = paste(paste("Treatment\nR²: ", round(data.frame(oars.asv.permanova)[1,3], 3)*100, "%",
  #                             "\n p: ", round(data.frame(oars.asv.permanova)[1,5], 3), sep="")),
  #         size=3.3)+
  facet_wrap(~"ASV PCoA")+
  labs(x=paste("Axis 1: ", round(unique(oars.asv.pcoa.df$var1), digits=2), "%", sep=""), 
       y=paste("Axis 2: ", round(unique(oars.asv.pcoa.df$var2), digits=2), "%", sep=""))
oars.asv.pcoa.plot



# :: ASV Beta Diversity ---------------------------------------------------

beta.trajectory.data = beta.trajectory(oars.asv.bray)
  
beta.trajectory.data =   merge(metadata.oars.stool.asv, 
                               beta.trajectory.data, by="standard.name")


# stats
stats.beta.between.036 = lmerTest::lmer(scale((between.beta)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                                   subset(beta.trajectory.data, timing %in% c("0M", "3M", "6M")& 
                                            compliant == TRUE)) %>%
  summary() %>% coef()

stats.beta.between.6912 = lmerTest::lmer(scale((between.beta)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                                    subset(beta.trajectory.data, timing %in% c("6M", "9M", "12M")& 
                                             compliant == TRUE)) %>%
  summary() %>% coef()

# GAM
stats.beta.between.gam = mgcv::gam(scale((between.beta)) ~ s(scale(oars.days), k=10) + diagnosis + adj.fiber + s(HM, bs="re"),
                               data = subset(beta.trajectory.data, compliant == TRUE & !is.na(between.beta))%>%
                                 mutate(HM = as.factor(HM)),
                               family = gaussian(),
                               method = "REML") %>% summary() 

# plot
metadata.oars.stool.asv.beta.between.plot = ggplot(subset(beta.trajectory.data,compliant==TRUE),
                                              aes(x=oars.days, 
                                                  y=(between.beta)))+
  annotate("rect", xmin=0, xmax=max(subset(metadata.oars.stool, oars.on.rs=="onRS")$oars.days),
           ymin=-Inf, ymax=Inf,
           alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(aes(fill=rs.col, shape=baseline), size=2)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  geom_smooth(color="black")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12))+
  geom_text(data=pval.positioner("between.beta",
                                 subset(beta.trajectory.data, compliant == TRUE),
                                 stats.beta.between.036,
                                 stats.beta.between.6912,
                                 yadjust = 1.05), aes(x=xval, y=yval, label=label), vjust=1,size=4)+
  labs(x="Days since starting RS", y="Mean Bray-Curtis Dissimilarity")
metadata.oars.stool.asv.beta.between.plot



# stats
stats.beta.within.036 = lmerTest::lmer(scale((within.beta)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                                        subset(beta.trajectory.data, timing %in% c("0M", "3M", "6M")& 
                                                 compliant == TRUE)) %>%
  summary() %>% coef()

stats.beta.within.6912 = lmerTest::lmer(scale((within.beta)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                                         subset(beta.trajectory.data, timing %in% c("6M", "9M", "12M")& 
                                                  compliant == TRUE)) %>%
  summary() %>% coef()

# GAM
stats.beta.within.gam = mgcv::gam(scale((within.beta)) ~ s(scale(oars.days), k=10) + diagnosis + adj.fiber + s(HM, bs="re"),
                                   data = subset(beta.trajectory.data, compliant == TRUE & !is.na(within.beta))%>%
                                     mutate(HM = as.factor(HM)),
                                   family = gaussian(),
                                   method = "REML") %>% summary() 
# plot
metadata.oars.stool.asv.beta.within.beta.plot = ggplot(subset(beta.trajectory.data,compliant==TRUE),
                                                   aes(x=oars.days, 
                                                       y=within.beta))+
  annotate("rect", xmin=0, xmax=max(subset(metadata.oars.stool, oars.on.rs=="onRS")$oars.days),
           ymin=-Inf, ymax=Inf,
           alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(aes(fill=rs.col, shape=baseline), size=2)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  geom_smooth(color="black")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12))+
  geom_text(data=pval.positioner("within.beta",
                                 subset(beta.trajectory.data, compliant == TRUE),
                                 stats.beta.within.036,
                                 stats.beta.within.6912,
                                 yadjust = 1.05), aes(x=xval, y=yval, label=label), vjust=1,size=4)+
  labs(x="Days since starting RS", y="Intra-individual Dissimilarity")
metadata.oars.stool.asv.beta.within.beta.plot



# :: ASV Maaslin2 ---------------------------------------------------------

# note: remove covariate adjustments (adj.fiber + diagnosis) because we are probably underpowered

# asv or glom: oars.asv.data.median | oars.asv.data.glom

oars.asv.data.glom[rownames(oars.asv.data.glom) %in% subset(metadata.oars.stool.asv, timing %in% c("0M", "3M", "6M")& compliant==TRUE)$standard.name,] %>% nrow()
# n = 33 samples for treatment phase, compliant

oars.asv.maaslin.036 = Maaslin2::Maaslin2(input_data = oars.asv.data.glom,
                                 input_metadata = subset(metadata.oars.stool.asv, timing %in% c("0M", "3M", "6M") & compliant==TRUE),
                                 output = "~/Downloads",         # Discard to downloads
                                 fixed_effects = c("oars.days"), # Note: do not adjust for anything (e.g. diagnosis, adj.fiber)
                                 random_effects = c("HM"),       # Adjust for repeat measures
                                 normalization = "TSS",          # Apply % normalization
                                 transform = "LOG",              # Apply log transformation
                                 analysis_method = "LM",         # Use linear model
                                 plot_scatter = FALSE,           # Disable scatterplot generation
                                 plot_heatmap = FALSE,           # Disable heatmap generation
                                 max_significance = 0.05,        # Ignore; we will re-calculate later
                                 standardize = TRUE              # Standardize coefficients
)

oars.asv.maaslin.036 = oars.asv.maaslin.036$results %>% data.frame() %>% arrange(pval)

# recalculate padj minus diagnosis fixed effect (if they were used)
oars.asv.maaslin.036 = subset(oars.asv.maaslin.036, metadata == "oars.days") %>%
  mutate(padj = p.adjust(pval, method="BH"))

subset(oars.asv.maaslin.036, grepl("bromii", feature))
subset(oars.asv.maaslin.036, grepl("adolescentis", feature))

oars.asv.maaslin.036.volcano = ggplot(subset(oars.asv.maaslin.036, value == "oars.days"),
                                      aes(x=coef, y=padj))+
  geom_hline(yintercept=(0.2), linetype=2, alpha=0.5)+
  geom_vline(xintercept=c(-1,1), linetype=2, alpha=0.5)+
  scale_fill_gradient2(low="blue", high="red")+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.05,0.1, 0.20, 0.5, 1))+
  #ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.2 & abs(coef)>1, gsub("\\..*", "", feature), NA)),size=2.5)+
  geom_point(shape=21, aes(fill=coef))+
  theme_classic()+theme(legend.position="none",
                        strip.text=element_text(size=10))+
  facet_wrap(~"Treatment: ASV")+
  labs(x="Adjusted Coefficient",
       y="FDR")
oars.asv.maaslin.036.volcano
# Christensenellales decreased and Brotophodocola increased over RS

oars.asv.data.glom[rownames(oars.asv.data.glom) %in% subset(metadata.oars.stool.asv, timing %in% c("9M", "12M", "6M")& compliant==TRUE)$standard.name,] %>% nrow()
# n = 26 samples for washout, compliant

oars.asv.maaslin.6912 = Maaslin2::Maaslin2(input_data = oars.asv.data.glom,
                                          input_metadata = subset(metadata.oars.stool.asv, timing %in% c("6M", "9M", "12M")& compliant==TRUE),
                                          output = "~/Downloads",
                                          fixed_effects = c("oars.days"),  # Example fixed effects
                                          random_effects = c("HM"),       # Example random effects
                                          normalization = "TSS",                       # Total Sum Scaling normalization
                                          transform = "LOG",                           # Log transformation
                                          analysis_method = "LM",                      # Linear model
                                          plot_scatter = FALSE,                        # Disable scatterplot generation
                                          plot_heatmap = FALSE,                         # Keep heatmap generation (optional)
                                          max_significance = 0.05,                     # Significance threshold for q-values
                                          standardize = TRUE                           # Disable standardization (optional)
)

oars.asv.maaslin.6912 = oars.asv.maaslin.6912$results %>% data.frame() %>% arrange(pval)

# recalculate padj minus diagnosis fixed effect
oars.asv.maaslin.6912 = subset(oars.asv.maaslin.6912, metadata == "oars.days") %>%
  mutate(padj = p.adjust(pval, method="BH"))

subset(oars.asv.maaslin.6912, grepl("bromii", feature))
subset(oars.asv.maaslin.6912, grepl("adolescentis", feature))

oars.asv.maaslin.6912.volcano = ggplot(oars.asv.maaslin.6912,
                                      aes(x=coef, y=(padj)))+
  geom_hline(yintercept=(0.2), linetype=2, alpha=0.5)+
  geom_vline(xintercept=c(-1,1), linetype=2, alpha=0.5)+
  scale_fill_gradient2(low="blue", high="red")+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.05,0.1, 0.20, 0.5, 1))+
  #ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.2 & abs(coef)>1, gsub("\\..*", "", feature), NA)),size=2.5)+
  geom_point(shape=21, aes(fill=coef))+
  theme_classic()+theme(legend.position="none",
                        strip.text=element_text(size=10))+
  facet_wrap(~"Washout: ASV")+
  labs(x="Adjusted Coefficient",
       y="FDR")
oars.asv.maaslin.6912.volcano
# Negativibacillus, Anaerofustis, Hungatella, and Bacillus increased over washout

oars.asv.maaslin.036.volcano+oars.asv.maaslin.6912.volcano


# :: ASV LFC Heatmaps ---------------------------------------------------------

# Note: July 31, 2025: Only look at butyrogens

# split 0-6 and 6-12 M heatmaps
# use LFC

# remove non-compliant
oars.asv.data.glom.1 = oars.asv.data.glom[!rownames(oars.asv.data.glom) %in% non.compliant,]

# identify most prevalent taxa
oars.asv.data.glom.pa = oars.asv.data.glom.1
oars.asv.data.glom.pa[oars.asv.data.glom.pa!=0] <- 1
# by prevalence (>=80% prev)
oars.asv.selected.taxa = data.frame(prev = colSums(oars.asv.data.glom.pa)) %>% 
  mutate(taxa = rownames(.))%>%
  subset(grepl(paste(c("Bifidobacterium", "Ruminococcus","charta", "Anaerobutyricum", "Anaerostipes","Clostridium","Blautia", "Faecalibacterium", "Eubact", "Roseburia", "Lachnospira"), collapse="|"), taxa)) %>%
  subset(!grepl(paste(c("gnavus", "torques", "innocuum"), collapse="|"), taxa)) %>%
  subset(prev >= nrow(oars.asv.data.glom.pa)*.8 | grepl(paste(c("bromii", "adolescent"), collapse="|"), taxa)) %>% row.names() %>% sort()

# now process for LFC
oars.asv.data.glom.lfc = oars.asv.data.glom.1 
# normalize + pseudo + filter
oars.asv.data.glom.lfc = (oars.asv.data.glom.lfc+min(oars.asv.data.glom.lfc[oars.asv.data.glom.lfc!=0])/2) / 50000
# melt
oars.asv.data.glom.lfc = oars.asv.data.glom.lfc %>% reshape2::melt()
colnames(oars.asv.data.glom.lfc) = c("standard.name", "taxa", "value")
# merge
oars.asv.data.glom.lfc = merge(oars.asv.data.glom.lfc,
                               metadata.oars.stool, by="standard.name")
# lfc from baseline
oars.asv.data.glom.lfc = oars.asv.data.glom.lfc %>%
  group_by(HM, taxa) %>%
  mutate(lfc.rs = log2(value / value[oars.on.rs == "preRS"])) %>%
  # note: made // defunct, because it's more straightforward to plot a single plot
  # i.e., do washout samples cluster with the same patient again
  #mutate(lfc.washout = log2(value / value[timing == "6M"])) %>%  # needs to be 6M because it's 0, 3, 6, 9, 12!
  subset(timing != "0M") %>% data.frame()

# create new name as "RS 1" or "WO 1", etc
oars.asv.data.glom.lfc = oars.asv.data.glom.lfc %>%
  group_by(HM) %>%
  mutate(hm_timing = paste(HM, timing, sep=" ")) 
  

# cast
oars.asv.lfc.treatment = reshape2::acast(oars.asv.data.glom.lfc %>% subset(timing %in% c("3M", "6M", "9M", "12M")),
                                         hm_timing ~ taxa, value.var="lfc.rs")

# or, subset to most prevalent taxa
oars.asv.lfc.treatment.prev = oars.asv.lfc.treatment[,colnames(oars.asv.lfc.treatment) %in% 
                                                       oars.asv.selected.taxa]
# prepare mapping for HM, Timing, and RS_Name
oars.asv.lfc.treatment.prev.mapping = data.frame(
  hm_timing = rownames(oars.asv.lfc.treatment.prev)
) %>% merge(oars.asv.data.glom.lfc[,c("hm_timing", "HM", "timing", "oars.rs.1", "oars.rs.2")] %>% distinct(), by="hm_timing") %>%
  mutate(RS_Name = ifelse(timing == "3M", oars.rs.1, oars.rs.2))
rownames(oars.asv.lfc.treatment.prev.mapping) = oars.asv.lfc.treatment.prev.mapping$hm_timing
oars.asv.lfc.treatment.prev.mapping = oars.asv.lfc.treatment.prev.mapping[,c("HM", "timing", "RS_Name")]
oars.asv.lfc.treatment.prev.mapping$RS_Name = factor(oars.asv.lfc.treatment.prev.mapping$RS_Name, levels=rs.names)
colnames(oars.asv.lfc.treatment.prev.mapping) = c("Patient", "Timing", "RS")
# add colors
oars.asv.lfc.treatment.prev.mapping.colors = list(
  Patient = c(HM0618 = RColorBrewer::brewer.pal(n = 12, name = "Set3")[1],
              HM0819 = RColorBrewer::brewer.pal(n = 12, name = "Set3")[2],
              HM0844 = RColorBrewer::brewer.pal(n = 12, name = "Set3")[3],
              HM0874 = RColorBrewer::brewer.pal(n = 12, name = "Set3")[4],
              HM0883 = RColorBrewer::brewer.pal(n = 12, name = "Set3")[5],
              HM0899 = RColorBrewer::brewer.pal(n = 12, name = "Set3")[6],
              HM0902 = RColorBrewer::brewer.pal(n = 12, name = "Set3")[7],
              HM0903 = RColorBrewer::brewer.pal(n = 12, name = "Set3")[8],
              HM0906 = RColorBrewer::brewer.pal(n = 12, name = "Set3")[9],
              HM0924 = RColorBrewer::brewer.pal(n = 12, name = "Set3")[10],
              HM0932 = RColorBrewer::brewer.pal(n = 12, name = "Set3")[11],
              HM0759 = RColorBrewer::brewer.pal(n = 12, name = "Set3")[12]),
  RS = c(Authentic = labelcolors$cols[c(1)],
         MSPrebiotic = labelcolors$cols[c(1)],
         LetsDoOrganic = labelcolors$cols[c(4)],
         Novelose330 = labelcolors$cols[c(5)],
         ActistarRT = labelcolors$cols[c(7)],
         FibersymRW = labelcolors$cols[c(7)],
         Versafibe1490 = labelcolors$cols[c(9)]),
  Timing = c(`3M` = "black", `6M` = "black", `9M` = "grey", `12M` = "white")
)

# heatmaps
pheatmap::pheatmap(t(oars.asv.lfc.treatment.prev),
                   color=colorRampPalette(c("blue","white", "red"))(100),
                   # note: use correlation when it makes sense
                   # here, the goal is to accentuate taxa that correlate (mostly, but not always, intra-genera correlations)
                   # and samples that correlate (mostly intra-individual correlations)
                   clustering_distance_rows = "correlation", 
                   clustering_distance_cols = "correlation",
                   breaks=c(seq(min(oars.asv.lfc.treatment.prev), 0, length.out=ceiling(100/2) + 1), 
                            seq(max(oars.asv.lfc.treatment.prev)/100, max(oars.asv.lfc.treatment.prev), length.out=floor(100/2))),
                   annotation_col=oars.asv.lfc.treatment.prev.mapping,
                   annotation_colors = oars.asv.lfc.treatment.prev.mapping.colors)


# :: But/Fcal by Time  -----------------------------------------------------

# QQ plots show log10-scaled is more linear
qqnorm(lmerTest::lmer(scale((but.i)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                      subset(metadata.oars.stool.asv, timing %in% c("0M", "3M") & compliant == TRUE))%>%resid())
qqnorm(lmerTest::lmer(scale(log10(but.i)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                      subset(metadata.oars.stool.asv, timing %in% c("0M", "3M") & compliant == TRUE))%>%resid())

# Butyrogens
stats.but.i.12 = lmerTest::lmer(scale(log10(but.i)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                                subset(metadata.oars.stool.asv, timing %in% c("0M", "3M") & compliant == TRUE)) %>%
  summary() %>% coef() %>% data.frame()
stats.but.i.23 = lmerTest::lmer(scale(log10(but.i)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                                subset(metadata.oars.stool.asv, timing %in% c("3M", "6M") & compliant == TRUE)) %>%
  summary() %>% coef()%>% data.frame()

# plot
metadata.oars.but.i.12.plot = ggplot(subset(metadata.oars.stool.asv, compliant == TRUE)%>%
                                       subset(timing %in% c("0M", "3M")),
                                     aes(x=oars.days, 
                                         y=but.i*100))+
  #geom_boxplot(width=0.5)+
  geom_smooth(method="lm", color="black")+
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(aes(fill=rs.col, shape=baseline), size=2.5)+
  # scale_y_log10()+
  scale_shape_manual(values=c(23,21))+
  #geom_text(aes(label=rs.col), size=3)+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  facet_wrap(~"RS 1")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text=element_text(size=10))+
  annotate(geom="text", x=1.5, y=80, 
           label=paste("p:", round(stats.but.i.12[2,5], digits=2), sep=" "))+
  labs(x="Days", y="Butyrogens (%)")
metadata.oars.but.i.12.plot
metadata.oars.but.i.23.plot = ggplot(subset(metadata.oars.stool.asv, compliant == TRUE)%>%
                                       subset(timing %in% c("3M", "6M")),
                                     aes(x=oars.days, 
                                         y=but.i*100))+
  #geom_boxplot(width=0.5)+
  geom_smooth(method="lm", color="black")+
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(aes(fill=rs.col, shape=baseline), size=2.5)+
  # scale_y_log10()+
  scale_shape_manual(values=c(21))+
  # geom_text(aes(label=rs.col), size=3)+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  facet_wrap(~"RS 2")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text=element_text(size=10))+
  annotate(geom="text", x=140, y=80, 
           label=paste("p:", round(stats.but.i.23[2,5], digits=2), sep=" "))+
  labs(x="Days", y="Butyrogens (%)")
metadata.oars.but.i.23.plot

# fecal cal
# QQ plots show log-scaled is more linear
qqnorm(lmerTest::lmer(scale((fcal)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                      subset(metadata.oars.stool, timing %in% c("0M", "3M") & compliant == TRUE))%>%resid())
qqnorm(lmerTest::lmer(scale(log10(fcal)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                      subset(metadata.oars.stool, timing %in% c("0M", "3M") & compliant == TRUE))%>%resid())

stats.fcal.12 = lmerTest::lmer(scale(log10(fcal)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                               subset(metadata.oars.stool, timing %in% c("0M", "3M") & compliant == TRUE)) %>%
  summary() %>% coef() %>% data.frame()
stats.fcal.23 = lmerTest::lmer(scale(log10(fcal)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                               subset(metadata.oars.stool, timing %in% c("3M", "6M") & compliant == TRUE)) %>%
  summary() %>% coef()%>% data.frame()

# plot
metadata.oars.fcal.12.plot = ggplot(subset(metadata.oars.stool.asv, compliant == TRUE)%>%
                                      subset(timing %in% c("0M", "3M")),
                                    aes(x=oars.days, 
                                        y=fcal))+
  #geom_boxplot(width=0.5)+
  geom_smooth(method="lm", color="black")+
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(aes(fill=rs.col, shape=baseline), size=2.5)+
  scale_y_log10()+
  scale_shape_manual(values=c(23,21))+
  geom_hline(yintercept=250, linetype=1, alpha=1, color="red")+
  #geom_text(aes(label=rs.col), size=3)+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  facet_wrap(~"RS 1")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text=element_text(size=10))+
  annotate(geom="text", x=1.5, y=3000, 
           label=paste("p:", round(stats.fcal.12[2,5], digits=2), sep=" "))+
  labs(x="Days", y="Fecal calprotectin (μg/g)")
metadata.oars.fcal.12.plot

metadata.oars.fcal.23.plot = ggplot(subset(metadata.oars.stool, compliant == TRUE)%>%
                                      subset(timing %in% c("3M", "6M")),
                                    aes(x=oars.days, 
                                        y=fcal))+
  #geom_boxplot(width=0.5)+
  geom_smooth(method="lm", color="black")+
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(aes(fill=rs.col, shape=baseline), size=2.5)+
  geom_hline(yintercept=250, linetype=1, alpha=1, color="red")+
  scale_y_log10()+
  scale_shape_manual(values=c(21))+
  #geom_text(aes(label=rs.col), size=3)+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  facet_wrap(~"RS 2")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text=element_text(size=10))+
  annotate(geom="text", x=140, y=3000, 
           label=paste("p:", round(stats.fcal.23[2,5], digits=3), sep=" "))+
  labs(x="Days", y="Fecal calprotectin (μg/g)")
metadata.oars.fcal.23.plot


# Starch:Mucin
stats.starch.mucin.12 = lmerTest::lmer(scale((starch.mucin)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM) + (1|plate),
                                subset(metadata.oars.stool.asv, timing %in% c("0M", "3M") & compliant == TRUE)) %>%
  summary() %>% coef() %>% data.frame()
stats.starch.mucin.23 = lmerTest::lmer(scale((starch.mucin)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM)+ (1|plate),
                                subset(metadata.oars.stool.asv, timing %in% c("3M", "6M") & compliant == TRUE)) %>%
  summary() %>% coef()%>% data.frame()

# plot
metadata.oars.starch.mucin.12.plot = ggplot(subset(metadata.oars.stool.asv, compliant == TRUE)%>%
                                       subset(timing %in% c("0M", "3M")),
                                     aes(x=oars.days, 
                                         y=starch.mucin))+
  #geom_boxplot(width=0.5)+
  geom_smooth(method="lm", color="black")+
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(aes(fill=rs.col, shape=baseline), size=2.5)+
  # scale_y_log10()+
  scale_shape_manual(values=c(23,21))+
  #geom_text(aes(label=rs.col), size=3)+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  facet_wrap(~"RS 1")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text=element_text(size=10))+
  annotate(geom="text", x=1.5, y=3, 
           label=paste("p:", round(stats.starch.mucin.12[2,5], digits=2), sep=" "))+
  labs(x="Days", y="Starch:Mucin CAZy Log2 Ratio")
metadata.oars.starch.mucin.12.plot
metadata.oars.starch.mucin.23.plot = ggplot(subset(metadata.oars.stool.asv, compliant == TRUE)%>%
                                       subset(timing %in% c("3M", "6M")),
                                     aes(x=oars.days, 
                                         y=starch.mucin))+
  #geom_boxplot(width=0.5)+
  geom_smooth(method="lm", color="black")+
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(aes(fill=rs.col, shape=baseline), size=2.5)+
  # scale_y_log10()+
  scale_shape_manual(values=c(21))+
  # geom_text(aes(label=rs.col), size=3)+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  facet_wrap(~"RS 2")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text=element_text(size=10))+
  annotate(geom="text", x=140, y=3, 
           label=paste("p:", round(stats.starch.mucin.23[2,5], digits=2), sep=" "))+
  labs(x="Days", y="Starch:Mucin CAZy Log2 Ratio")
metadata.oars.starch.mucin.23.plot


# plot
metadata.oars.but.i.12.plot+metadata.oars.but.i.23.plot+
  metadata.oars.fcal.12.plot+metadata.oars.fcal.23.plot+
  metadata.oars.starch.mucin.12.plot+metadata.oars.starch.mucin.23.plot+
  patchwork::plot_layout(nrow=1)


# >>> 2. MULTIOMIC ANALYSES -----------------------------------------------------------


# :: MGX Maaslin2 ---------------------------------------------------------

# 
oars.mgx.taxa[rownames(oars.mgx.taxa) %in% subset(metadata.oars.stool.asv, timing %in% c("0M", "3M", "6M")& compliant==TRUE)$standard.name,] %>% nrow()
# n = 32 samples for treatment, compliant

oars.mgx.maaslin.036 = Maaslin2::Maaslin2(input_data = oars.mgx.taxa,
                                          input_metadata = subset(metadata.oars.stool.asv, timing %in% c("0M", "3M", "6M")& compliant==TRUE),
                                          output = "~/Downloads",
                                          fixed_effects = c("oars.days"),  # Example fixed effects
                                          random_effects = c("HM"),       # Example random effects
                                          normalization = "TSS",                       # Total Sum Scaling normalization
                                          transform = "LOG",                           # Log transformation
                                          analysis_method = "LM",                      # Linear model
                                          plot_scatter = FALSE,                        # Disable scatterplot generation
                                          plot_heatmap = FALSE,                         # Keep heatmap generation (optional)
                                          max_significance = 0.05,                     # Significance threshold for q-values
                                          standardize = TRUE                           # Disable standardization (optional)
)

oars.mgx.maaslin.036 = oars.mgx.maaslin.036$results %>% data.frame()

# recalculate padj minus diagnosis fixed effect
oars.mgx.maaslin.036 = subset(oars.mgx.maaslin.036, metadata == "oars.days") %>%
  mutate(padj = p.adjust(pval, method="BH"))

oars.mgx.maaslin.036.volcano = ggplot(oars.mgx.maaslin.036,
                                      aes(x=coef, y=(padj)))+
  geom_point(shape=21, aes(fill=coef))+
  geom_hline(yintercept=(0.2), linetype=2, alpha=0.5)+
  geom_vline(xintercept=c(-1,1), linetype=2, alpha=0.5)+
  scale_fill_gradient2(low="blue", high="red")+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.05,0.1, 0.20, 0.5, 1))+
  # ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20 & abs(coef) > 1, gsub("\\..*", "", feature), NA)), size=2.5)+  
  facet_wrap(~"Treatment: Species")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text=element_text(size=10))+
  labs(x="Adjusted Coefficient",
       y="FDR")
oars.mgx.maaslin.036.volcano
oars.mgx.maaslin.036 %>% head()
# several depleted

oars.mgx.taxa[rownames(oars.mgx.taxa) %in% subset(metadata.oars.stool.asv, timing %in% c("9M", "12M", "6M")& compliant==TRUE)$standard.name,] %>% nrow()
# n = 24 for washout, compliant

oars.mgx.maaslin.6912 = Maaslin2::Maaslin2(input_data = oars.mgx.taxa,
                                           input_metadata = subset(metadata.oars.stool.asv, timing %in% c("6M", "9M", "12M")& compliant==TRUE),
                                           output = "~/Downloads",
                                           fixed_effects = c("oars.days"),  # Example fixed effects
                                           random_effects = c("HM"),       # Example random effects
                                           normalization = "TSS",                       # Total Sum Scaling normalization
                                           transform = "LOG",                           # Log transformation
                                           analysis_method = "LM",                      # Linear model
                                           plot_scatter = FALSE,                        # Disable scatterplot generation
                                           plot_heatmap = FALSE,                         # Keep heatmap generation (optional)
                                           max_significance = 0.05,                     # Significance threshold for q-values
                                           standardize = TRUE                           # Disable standardization (optional)
)

oars.mgx.maaslin.6912 = oars.mgx.maaslin.6912$results %>% data.frame()

# recalculate padj minus diagnosis fixed effect
oars.mgx.maaslin.6912 = subset(oars.mgx.maaslin.6912, metadata == "oars.days") %>%
  mutate(padj = p.adjust(pval, method="BH"))

oars.mgx.maaslin.6912.volcano = ggplot(oars.mgx.maaslin.6912,
                                       aes(x=coef, y=(padj)))+
  geom_point(shape=21, aes(fill=coef))+
  geom_hline(yintercept=(0.2), linetype=2, alpha=0.5)+
  geom_vline(xintercept=c(-1,1), linetype=2, alpha=0.5)+
  scale_fill_gradient2(low="blue", high="red")+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.05,0.1, 0.20, 0.5, 1))+
  #ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20 & abs(coef) > 1, gsub("\\..*", "", feature), NA)), size=2.5)+  
  facet_wrap(~"Washout: Species")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text=element_text(size=10))+
  labs(x="Adjusted Coefficient",
       y="FDR")
oars.mgx.maaslin.6912.volcano
# nothing sig over washout
oars.mgx.maaslin.6912 %>% head()


oars.mgx.maaslin.036.volcano+oars.mgx.maaslin.6912.volcano

# :: MGX Bray-Curtis ------------------------------------------------------

# remove non-compliant before prevalence filter
oars.mgx.taxa.1 = oars.mgx.taxa[!rownames(oars.mgx.taxa) %in% non.compliant,]

# filter to 10% prevalence
dim(oars.mgx.taxa) # 958 taxa
oars.mgx.taxa.filt.1 = oars.mgx.taxa.1
oars.mgx.taxa.filt.1[oars.mgx.taxa.filt.1!=0] = 1
oars.mgx.taxa.filt.1 = oars.mgx.taxa.1[,colSums(oars.mgx.taxa.1) >= nrow(oars.mgx.taxa.filt.1)*0.1]
dim(oars.mgx.taxa.filt.1) # 120


# calculate Bray-Curtis dissimilarities
oars.mgx.bray = vegan::vegdist(oars.mgx.taxa.filt.1, method="bray") 
# perform PCoA
oars.mgx.pcoa = ape::pcoa(oars.mgx.bray)
# extract data from pcoa
oars.mgx.pcoa.df = data.frame(oars.mgx.pcoa$vectors[,c(1:2)])
oars.mgx.pcoa.df$standard.name = rownames(oars.mgx.pcoa.df)
# add metadata
oars.mgx.pcoa.df = merge(oars.mgx.pcoa.df,
                         metadata.oars.stool.asv, by="standard.name")
# extract variance explained
oars.mgx.pcoa.var_exp = oars.mgx.pcoa$values[c(1:2),2]
oars.mgx.pcoa.df$var1 = round(oars.mgx.pcoa.var_exp[1]*100, digits=2)
oars.mgx.pcoa.df$var2 = round(oars.mgx.pcoa.var_exp[2]*100, digits=2)

# clean up "oars.on.rs" variable
rownames(oars.mgx.pcoa.df) = oars.mgx.pcoa.df$standard.name

set.seed(25)
t1 = Sys.time()
oars.mgx.permanova = vegan::adonis2(oars.mgx.bray ~ on.rs + adj.fiber + HM,
                                    oars.mgx.pcoa.df %>% 
                                      mutate(on.rs = as.factor(ifelse(oars.on.rs=="onRS", "onRS", "offRS"))),
                                    strata = oars.mgx.pcoa.df$HM,
                                    by="margin")
t2 = Sys.time()
t2 - t1
oars.mgx.permanova # sig

oars.mgx.pcoa.plot <- ggplot(
  data=oars.mgx.pcoa.df %>% group_by(HM) %>% arrange(stool_date_rec_v2), 
  aes(x=Axis.1, y=Axis.2))+
  geom_path(aes(group=HM), color="black", linetype=2, alpha=0.5, linewidth=0.3) + 
  geom_point(aes(fill=rs.col, shape=baseline), size=2)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  #geom_text(aes(label=RS_Name))+
  #annotate(geom="text",
  #         x=0.25,
  #         y=-0.20,
  #         label = paste(paste("Treatment\nR²: ", round(data.frame(oars.mgx.permanova)[1,3], 3)*100, "%",
  #                             "\n p: ", round(data.frame(oars.mgx.permanova)[1,5], 3), sep="")),
  #         size=3.3)+
  facet_wrap(~"Species PCoA")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text=element_text(size=10))+
  labs(x=paste("Axis 1: ", round(unique(oars.mgx.pcoa.df$var1), digits=2), "%", sep=""), 
       y=paste("Axis 2: ", round(unique(oars.mgx.pcoa.df$var2), digits=2), "%", sep=""))
oars.mgx.pcoa.plot



# :: MPX KEGG -------------------------------------------------------------

# load KEGG Pathway data
oars.mpx.kegg.mat = readRDS("./metaproteomics/2026_01_15_oars_mpx_kegg.Rds")

# create new metadata containing matching samples
metadata.oars.stool.mpx = subset(metadata.oars.stool,
                                 standard.name %in% rownames(oars.mpx.kegg.mat) & !standard.name %in% "HM0924-STL-12")

# create map for protein names
oars.mpx.keggprotein.names = data.frame(good = colnames(oars.mpx.kegg.mat),
                                       feature = make.names(colnames(oars.mpx.kegg.mat)))
length(unique(oars.mpx.keggprotein.names$good))
# 180 KEGG pathways

# removen non-compliant
oars.mpx.kegg.mat.1 = oars.mpx.kegg.mat[!rownames(oars.mpx.kegg.mat) %in% c(non.compliant, "HM0924-STL-12"),] %>% as.data.frame()

# filter to 10% prevalence
dim(oars.mpx.kegg.mat.1) # 180 KEGG
oars.mpx.kegg.mat.filt.1 = oars.mpx.kegg.mat.1
oars.mpx.kegg.mat.filt.1[is.na(oars.mpx.kegg.mat.filt.1)] = 0
oars.mpx.kegg.mat.filt.1[oars.mpx.kegg.mat.filt.1!=0] = 1
oars.mpx.kegg.mat.filt.1 = oars.mpx.kegg.mat.1[,colSums(oars.mpx.kegg.mat.filt.1) >= nrow(oars.mpx.kegg.mat.filt.1)*0.1]
dim(oars.mpx.kegg.mat.filt.1) # 180

# :: MPX KEGG PCA -------------------------------------------------------------

# replace NA with 0
oars.mpx.kegg.mat.filt.1[is.na(oars.mpx.kegg.mat.filt.1)] = 0
# log transform (with pseudo)
oars.mpx.kegg.oars.pca = log2(oars.mpx.kegg.mat.filt.1+(min(oars.mpx.kegg.mat.filt.1[oars.mpx.kegg.mat.filt.1!=0])/2))
# PCA
oars.mpx.kegg.oars.pca = prcomp((oars.mpx.kegg.oars.pca), scale=T)
oars.mpx.kegg.oars.pca.df = oars.mpx.kegg.oars.pca$x[,c(1,2)] %>% data.frame() %>%
  rownames_to_column("standard.name")
oars.mpx.kegg.oars.pca.df = merge(oars.mpx.kegg.oars.pca.df,
                                 metadata.oars.stool.mpx, by="standard.name")
oars.mpx.kegg.oars.pca.var <- (oars.mpx.kegg.oars.pca$sdev^2 / sum(oars.mpx.kegg.oars.pca$sdev^2)) * 100

rownames(oars.mpx.kegg.oars.pca.df) = oars.mpx.kegg.oars.pca.df$standard.name
nrow(oars.mpx.kegg.oars.pca.df)

# permanova
set.seed(25)
t1 = Sys.time()
oars.mpx.kegg.oars.pca.permanova = vegan::adonis2(dist(oars.mpx.kegg.oars.pca$x) ~ on.rs + adj.fiber + HM + plate,
                                                 oars.mpx.kegg.oars.pca.df %>% mutate(on.rs = as.factor(ifelse(oars.on.rs=="onRS", "onRS", "offRS"))),
                                                 strata = oars.mpx.kegg.oars.pca.df$HM,
                                                 by="margin")
t2 = Sys.time()
t2 - t1

oars.mpx.kegg.pca.plot <- ggplot(
  data=oars.mpx.kegg.oars.pca.df %>% group_by(HM) %>% arrange(stool_date_rec_v2), 
  aes(x=PC1, y=PC2))+
  geom_path(aes(group=HM), color="black", linetype=2, alpha=0.5, linewidth=0.3) + 
  geom_point(aes(fill=rs.col, shape=baseline), size=2)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=c("grey", labelcolors$cols[c(1,1,4,5,5,8,8,9)]))+
  #annotate(geom="text",
  #         x=20,
  #         y=9,
  #         label = paste(paste("Treatment\nR²: ", round(data.frame(oars.mpx.kegg.oars.pca.permanova)[1,3], 3)*100, "%",
  #                             "\n p: ", round(data.frame(oars.mpx.kegg.oars.pca.permanova)[1,5], 3), sep="")),
  #         size=3.3)+
  facet_wrap(~"Pathway PCA")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text=element_text(size=10))+
  labs(x=paste("PC1: ", round((oars.mpx.kegg.oars.pca.var)[1], digits=2), "%", sep=""), 
       y=paste("PC2: ", round((oars.mpx.kegg.oars.pca.var)[2], digits=2), "%", sep=""))
oars.mpx.kegg.pca.plot


# :: MPX KEGG Maaslin2 ----------------------------------------------------

oars.mpx.kegg.mat[rownames(oars.mpx.kegg.mat) %in% subset(metadata.oars.stool, timing %in% c("0M", "3M", "6M")& compliant==TRUE)$standard.name,] %>% nrow()
# n = 32 treatment, compliant

oars.mpx.kegg.maaslin.036 = Maaslin2::Maaslin2(input_data = (oars.mpx.kegg.mat),
                                              input_metadata = subset(metadata.oars.stool, timing %in% c("0M", "3M", "6M") & compliant==TRUE),
                                              output = "~/Downloads",
                                              fixed_effects = c("oars.days"),  # Example fixed effects
                                              random_effects = c("HM", "plate"),       # Example random effects
                                              normalization = "NONE",                       # Total Sum Scaling normalization
                                              transform = "LOG",                           # Log transformation
                                              analysis_method = "LM",                      # Linear model
                                              plot_scatter = FALSE,                        # Disable scatterplot generation
                                              plot_heatmap = FALSE,                         # Keep heatmap generation (optional)
                                              max_significance = 0.05,                     # Significance threshold for q-values
                                              standardize = TRUE                           # Disable standardization (optional)
)

oars.mpx.kegg.maaslin.036 = oars.mpx.kegg.maaslin.036$results %>% data.frame() %>% arrange(pval)

# recalculate padj minus diagnosis fixed effect (if used)
oars.mpx.kegg.maaslin.036 = subset(oars.mpx.kegg.maaslin.036, metadata == "oars.days") %>%
  mutate(padj = p.adjust(pval, method="BH"))

# Enriched: Nitrogen reg protein, Isocitrate-isopropylmalate dehydrogenase, tRNA methylthiotransferase, Histidinol dehydrogenase, Pentose-5-phosphate-3-epimerase

# fix feature names
oars.mpx.kegg.maaslin.036$feature = oars.mpx.keggprotein.names$good[match(oars.mpx.kegg.maaslin.036$feature, oars.mpx.keggprotein.names$feature)]

oars.mpx.kegg.maaslin.036.volcano = ggplot(oars.mpx.kegg.maaslin.036,
                                           aes(x=coef, y=(padj)))+
  geom_point(shape=21, aes(fill=coef))+
  geom_hline(yintercept=(0.2), linetype=2, alpha=0.5)+
  geom_vline(xintercept=c(-1,1), linetype=2, alpha=0.5)+
  scale_fill_gradient2(low="blue", high="red")+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.05,0.1, 0.20, 0.5, 1))+
 #  ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, gsub("\\..*", "", feature), NA)), size=2.5)+  
  facet_wrap(~"Treatment: Pathway")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text=element_text(size=10))+
  labs(x="Adjusted Coefficient",
       y="FDR")
oars.mpx.kegg.maaslin.036.volcano
oars.mpx.kegg.maaslin.036 %>% head()
# C5-branched dibasic acid metabolism, and beta-Lactam resistance
subset(oars.mpx.kegg.maaslin.036, padj < 0.20) %>% nrow()


oars.mpx.kegg.mat[rownames(oars.mpx.kegg.mat) %in% subset(metadata.oars.stool, timing %in% c("6M", "9M", "12M")& compliant==TRUE)$standard.name,] %>% nrow()
# n = 26 washout, compliant

oars.mpx.kegg.maaslin.6912 = Maaslin2::Maaslin2(input_data = oars.mpx.kegg.mat,
                                               input_metadata = subset(metadata.oars.stool, timing %in% c("6M", "9M", "12M")& compliant==TRUE),
                                               output = "~/Downloads",
                                               fixed_effects = c("oars.days"),  # Example fixed effects
                                               random_effects = c("HM", "plate"),       # Example random effects
                                               normalization = "NONE",                       # Total Sum Scaling normalization
                                               transform = "LOG",                           # Log transformation
                                               analysis_method = "LM",                      # Linear model
                                               plot_scatter = FALSE,                        # Disable scatterplot generation
                                               plot_heatmap = FALSE,                         # Keep heatmap generation (optional)
                                               max_significance = 0.05,                     # Significance threshold for q-values
                                               standardize = TRUE                           # Disable standardization (optional)
)

oars.mpx.kegg.maaslin.6912 = oars.mpx.kegg.maaslin.6912$results %>% data.frame() %>% arrange(pval)

# recalculate padj minus diagnosis fixed effect (if used)
oars.mpx.kegg.maaslin.6912 = subset(oars.mpx.kegg.maaslin.6912, metadata == "oars.days") %>%
  mutate(padj = p.adjust(pval, method="BH"))

# fix feature names
oars.mpx.kegg.maaslin.6912$feature = oars.mpx.keggprotein.names$good[match(oars.mpx.kegg.maaslin.6912$feature, oars.mpx.keggprotein.names$feature)]

oars.mpx.kegg.maaslin.6912.volcano = ggplot(oars.mpx.kegg.maaslin.6912,
                                            aes(x=coef, y=(padj)))+
  geom_point(shape=21, aes(fill=coef))+
  geom_hline(yintercept=(0.2), linetype=2, alpha=0.5)+
  geom_vline(xintercept=c(-1,1), linetype=2, alpha=0.5)+
  scale_fill_gradient2(low="blue", high="red")+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.05,0.1, 0.20, 0.5, 1))+
 # ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, gsub("\\..*", "", feature), NA)), size=2.5)+  
  facet_wrap(~"Washout: Pathway")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text=element_text(size=10))+
  labs(x="Adjusted Coefficient",
       y="FDR")
oars.mpx.kegg.maaslin.6912.volcano
# nothing

oars.mpx.kegg.maaslin.036.volcano+
  oars.mpx.kegg.maaslin.6912.volcano

# :: MPX COG -------------------------------------------------------------

# load COG data
# oars.mpx.cog.mat = readRDS("./metaproteomics/2026_01_15_oars_mpx_cog.Rds")

# use KEGG metadata
metadata.oars.stool.mpx

# create map for protein names
oars.mpx.cogprotein.names = data.frame(good = colnames(oars.mpx.cog.mat),
                                    feature = make.names(colnames(oars.mpx.cog.mat)))
length(unique(oars.mpx.cogprotein.names$good))
# 2653 unique COGs

# subset to compliant
oars.mpx.cog.mat.1 = oars.mpx.cog.mat[!rownames(oars.mpx.cog.mat) %in% c(non.compliant, "HM0924-STL-12"),]
# n = 48

# filter to 10% prevalence
dim(oars.mpx.cog.mat) # 2653 KEGG
oars.mpx.cog.mat.filt.1 = oars.mpx.cog.mat.1
oars.mpx.cog.mat.filt.1[is.na(oars.mpx.cog.mat.filt.1)] = 0
oars.mpx.cog.mat.filt.1[oars.mpx.cog.mat.filt.1!=0] = 1
oars.mpx.cog.mat.filt.1 = oars.mpx.cog.mat.1[,colSums(oars.mpx.cog.mat.filt.1) >= nrow(oars.mpx.cog.mat.filt.1)*0.1]
dim(oars.mpx.cog.mat.filt.1) # 2481

# :: MPX COG PCA -------------------------------------------------------------

# replace NA with 0
oars.mpx.cog.mat.filt.1[is.na(oars.mpx.cog.mat.filt.1)] = 0
# log transform with pseudo
oars.mpx.cog.oars.pca = log2(oars.mpx.cog.mat.filt.1+(min(oars.mpx.cog.mat.filt.1[oars.mpx.cog.mat.filt.1!=0])/2))
# PCA
oars.mpx.cog.oars.pca = prcomp((oars.mpx.cog.oars.pca), scale=T)
oars.mpx.cog.oars.pca.df = oars.mpx.cog.oars.pca$x[,c(1,2)] %>% data.frame() %>%
  rownames_to_column("standard.name")
oars.mpx.cog.oars.pca.df = merge(oars.mpx.cog.oars.pca.df,
                             metadata.oars.stool.mpx, by="standard.name")
oars.mpx.cog.oars.pca.var <- (oars.mpx.cog.oars.pca$sdev^2 / sum(oars.mpx.cog.oars.pca$sdev^2)) * 100

rownames(oars.mpx.cog.oars.pca.df) = oars.mpx.cog.oars.pca.df$standard.name
nrow(oars.mpx.cog.oars.pca.df)

# permanova
set.seed(25)
t1 = Sys.time()
oars.mpx.cog.oars.pca.permanova = vegan::adonis2(dist(oars.mpx.cog.oars.pca$x) ~ on.rs + adj.fiber + HM + plate,
                                             oars.mpx.cog.oars.pca.df %>% mutate(on.rs = as.factor(ifelse(oars.on.rs=="onRS", "onRS", "offRS"))),
                                             strata = oars.mpx.cog.oars.pca.df$HM,
                                             by="margin")
t2 = Sys.time()
t2 - t1
oars.mpx.cog.oars.pca.permanova

oars.mpx.cog.pca.plot <- ggplot(
  data=oars.mpx.cog.oars.pca.df %>% group_by(HM) %>% arrange(stool_date_rec_v2), 
  aes(x=PC1, y=PC2))+
  geom_path(aes(group=HM), color="black", linetype=2, alpha=0.5, linewidth=0.3) + 
  geom_point(aes(fill=rs.col, shape=baseline), size=2)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  # annotate(geom="text",
  #          x=-40,
  #          y=-30,
  #          label = paste(paste("Treatment\nR²: ", round(data.frame(oars.mpx.cog.oars.pca.permanova)[1,3], 3)*100, "%",
  #                              "\n p: ", round(data.frame(oars.mpx.cog.oars.pca.permanova)[1,5], 3), sep="")),
  #          size=3.3)+
  facet_wrap(~"COG PCA")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text=element_text(size=10))+
  labs(x=paste("PC1: ", round((oars.mpx.cog.oars.pca.var)[1], digits=2), "%", sep=""), 
       y=paste("PC2: ", round((oars.mpx.cog.oars.pca.var)[2], digits=2), "%", sep=""))
oars.mpx.cog.pca.plot



# :: MPX COG Maaslin2 ----------------------------------------------------

oars.mpx.cog.mat[rownames(oars.mpx.cog.mat) %in% subset(metadata.oars.stool, timing %in% c("0M", "3M", "6M")& compliant==TRUE)$standard.name,] %>% nrow()
# n = 32 treatment, compliant

oars.mpx.cog.maaslin.036 = Maaslin2::Maaslin2(input_data = (oars.mpx.cog.mat),
                                          input_metadata = subset(metadata.oars.stool, timing %in% c("0M", "3M", "6M") & compliant==TRUE),
                                          output = "~/Downloads",
                                          fixed_effects = c("oars.days"),  # Example fixed effects
                                          random_effects = c("HM", "plate"),       # Example random effects
                                          normalization = "NONE",                       # Total Sum Scaling normalization
                                          transform = "LOG",                           # Log transformation
                                          analysis_method = "LM",                      # Linear model
                                          plot_scatter = FALSE,                        # Disable scatterplot generation
                                          plot_heatmap = FALSE,                         # Keep heatmap generation (optional)
                                          max_significance = 0.05,                     # Significance threshold for q-values
                                          standardize = TRUE                           # Disable standardization (optional)
)

oars.mpx.cog.maaslin.036 = oars.mpx.cog.maaslin.036$results %>% data.frame() %>% arrange(pval)

# recalculate padj minus diagnosis fixed effect (if used)
oars.mpx.cog.maaslin.036 = subset(oars.mpx.cog.maaslin.036, metadata == "oars.days") %>%
  mutate(padj = p.adjust(pval, method="BH"))

# fix feature names
oars.mpx.cog.maaslin.036$feature = oars.mpx.cogprotein.names$good[match(oars.mpx.cog.maaslin.036$feature, oars.mpx.cogprotein.names$feature)]

oars.mpx.cog.maaslin.036.volcano = ggplot(oars.mpx.cog.maaslin.036,
                                          aes(x=coef, y=(padj)))+
  geom_point(shape=21, aes(fill=coef))+
  geom_hline(yintercept=(0.2), linetype=2, alpha=0.5)+
  geom_vline(xintercept=c(-1,1), linetype=2, alpha=0.5)+
  scale_fill_gradient2(low="blue", high="red")+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.05,0.1, 0.20, 0.5, 1))+
  # ggnetwork::geom_nodetext_repel(aes(label=ifelse(padj < 0.20, feature, NA)), size=2.5)+  
  facet_wrap(~"Treatment: COG")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text=element_text(size=10))+
  labs(x="Adjusted Coefficient",
       y="FDR")
oars.mpx.cog.maaslin.036.volcano
oars.mpx.cog.maaslin.036 %>% head()

oars.mpx.cog.mat[rownames(oars.mpx.cog.mat) %in% subset(metadata.oars.stool, timing %in% c("6M", "9M", "12M")& compliant==TRUE)$standard.name,] %>% nrow()
# n = 26 washout, compliant

oars.mpx.cog.maaslin.6912 = Maaslin2::Maaslin2(input_data = oars.mpx.cog.mat,
                                           input_metadata = subset(metadata.oars.stool, timing %in% c("6M", "9M", "12M")& compliant==TRUE),
                                           output = "~/Downloads",
                                           fixed_effects = c("oars.days"),  # Example fixed effects
                                           random_effects = c("HM", "plate"),       # Example random effects
                                           normalization = "NONE",                       # Total Sum Scaling normalization
                                           transform = "LOG",                           # Log transformation
                                           analysis_method = "LM",                      # Linear model
                                           plot_scatter = FALSE,                        # Disable scatterplot generation
                                           plot_heatmap = FALSE,                         # Keep heatmap generation (optional)
                                           max_significance = 0.05,                     # Significance threshold for q-values
                                           standardize = TRUE                           # Disable standardization (optional)
)

oars.mpx.cog.maaslin.6912 = oars.mpx.cog.maaslin.6912$results %>% data.frame() %>% arrange(pval)

# recalculate padj minus diagnosis fixed effect (if used)
oars.mpx.cog.maaslin.6912 = subset(oars.mpx.cog.maaslin.6912, metadata == "oars.days") %>%
  mutate(padj = p.adjust(pval, method="BH"))

# fix feature names
oars.mpx.cog.maaslin.6912$feature = oars.mpx.cogprotein.names$good[match(oars.mpx.cog.maaslin.6912$feature, oars.mpx.cogprotein.names$feature)]

oars.mpx.cog.maaslin.6912.volcano = ggplot(oars.mpx.cog.maaslin.6912,
                                           aes(x=coef, y=(padj)))+
  geom_point(shape=21, aes(fill=coef))+
  geom_hline(yintercept=(0.2), linetype=2, alpha=0.5)+
  geom_vline(xintercept=c(-1,1), linetype=2, alpha=0.5)+
  scale_fill_gradient2(low="blue", high="red")+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.05,0.1, 0.20, 0.5, 1))+
  # ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, gsub("\\..*", "", feature), NA)), size=2.5)+  
  facet_wrap(~"Washout: COG")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text=element_text(size=10))+
  labs(x="Adjusted Coefficient",
       y="FDR")
oars.mpx.cog.maaslin.6912.volcano
# few
oars.mpx.cog.maaslin.6912 %>% subset(padj < 0.20) %>% head()

oars.mpx.cog.maaslin.036.volcano+
  oars.mpx.cog.maaslin.6912.volcano

# :: MPX CAZy -------------------------------------------------------------

# load CAZy data
# oars.mpx.cazy.mat = readRDS("./metaproteomics/2025_06_28_oars_mpx_cazy.Rds")

# delete "-" (unannotated CAZy)
oars.mpx.cazy.mat = oars.mpx.cazy.mat[,colnames(oars.mpx.cazy.mat) != "-"]

# use KEGG metadata
metadata.oars.stool.mpx

# create map for protein names
oars.mpx.cazy.names = data.frame(good = colnames(oars.mpx.cazy.mat),
                                    feature = make.names(colnames(oars.mpx.cazy.mat)))

# filter to compliant
oars.mpx.cazy.mat.1 = oars.mpx.cazy.mat[!rownames(oars.mpx.cazy.mat) %in% c(non.compliant, "HM0924-STL-12"),]

# filter to 10% prevalence
dim(oars.mpx.cazy.mat) # 66 CAZy
oars.mpx.cazy.mat.filt.1 = oars.mpx.cazy.mat.1
oars.mpx.cazy.mat.filt.1[is.na(oars.mpx.cazy.mat.filt.1)] = 0
oars.mpx.cazy.mat.filt.1[oars.mpx.cazy.mat.filt.1!=0] = 1
oars.mpx.cazy.mat.filt.1 = oars.mpx.cazy.mat.1[,colSums(oars.mpx.cazy.mat.filt.1) >= nrow(oars.mpx.cazy.mat.filt.1)*0.1]
dim(oars.mpx.cazy.mat.filt.1) # 66 

# :: MPX CAZy PCA -------------------------------------------------------------

# replace NA with 0
oars.mpx.cazy.mat.filt.1[is.na(oars.mpx.cazy.mat.filt.1)] = 0
# log transform with pseudo
oars.mpx.cazy.pca = log2(oars.mpx.cazy.mat.filt.1+(min(oars.mpx.cazy.mat.filt.1[oars.mpx.cazy.mat.filt.1!=0])/2))
# PCA
oars.mpx.cazy.pca = prcomp((oars.mpx.cazy.pca), scale=T)
oars.mpx.cazy.pca.df = oars.mpx.cazy.pca$x[,c(1,2)] %>% data.frame() %>%
  rownames_to_column("standard.name")
oars.mpx.cazy.pca.df = merge(oars.mpx.cazy.pca.df,
                             metadata.oars.stool.mpx, by="standard.name")
oars.mpx.cazy.pca.var <- (oars.mpx.cazy.pca$sdev)^2 / sum(oars.mpx.cazy.pca$sdev^2) * 100

rownames(oars.mpx.cazy.pca.df) = oars.mpx.cazy.pca.df$standard.name
nrow(oars.mpx.cazy.pca.df)

# permanova
set.seed(25)
t1 = Sys.time()
oars.mpx.cazy.pca.permanova = vegan::adonis2(dist(oars.mpx.cazy.pca$x) ~ on.rs + adj.fiber + HM + plate,
                                             oars.mpx.cazy.pca.df %>% mutate(on.rs = as.factor(ifelse(oars.on.rs=="onRS", "onRS", "offRS"))),
                                             strata = oars.mpx.cazy.pca.df$HM,
                                             by="margin")
t2 = Sys.time()
t2 - t1
oars.mpx.cazy.pca.permanova

oars.mpx.cazy.pca.plot <- ggplot(
  data=oars.mpx.cazy.pca.df %>% group_by(HM) %>% arrange(stool_date_rec_v2), 
  aes(x=PC1, y=PC2))+
  geom_path(aes(group=HM), color="black", linetype=2, alpha=0.5, linewidth=0.3) + 
  geom_point(aes(fill=rs.col, shape=baseline), size=2)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  # annotate(geom="text",
  #          x=-7,
  #          y=-3,
  #          label = paste(paste("Treatment\nR²: ", round(data.frame(oars.mpx.cazy.pca.permanova)[1,3], 3)*100, "%",
  #                              "\n p: ", round(data.frame(oars.mpx.cazy.pca.permanova)[1,5], 3), sep="")),
  #          size=3.3)+
  facet_wrap(~"CAZy PCA")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text=element_text(size=10))+
  labs(x=paste("PC1: ", round((oars.mpx.cazy.pca.var)[1], digits=2), "%", sep=""), 
       y=paste("PC2: ", round((oars.mpx.cazy.pca.var)[2], digits=2), "%", sep=""))
oars.mpx.cazy.pca.plot

# very sig

# :: MPX CAZy Maaslin2 ----------------------------------------------------

oars.mpx.cazy.mat[rownames(oars.mpx.cazy.mat) %in% subset(metadata.oars.stool, timing %in% c("0M", "3M", "6M")& compliant==TRUE)$standard.name,] %>% nrow()
# n = 32 treatment, compliant

oars.mpx.cazy.maaslin.036 = Maaslin2::Maaslin2(input_data = (oars.mpx.cazy.mat),
                                          input_metadata = subset(metadata.oars.stool, timing %in% c("0M", "3M", "6M") & compliant==TRUE),
                                          output = "~/Downloads",
                                          fixed_effects = c("oars.days"),  # Example fixed effects
                                          random_effects = c("HM", "plate"),       # Example random effects
                                          normalization = "NONE",                       # Total Sum Scaling normalization
                                          transform = "LOG",                           # Log transformation
                                          analysis_method = "LM",                      # Linear model
                                          plot_scatter = FALSE,                        # Disable scatterplot generation
                                          plot_heatmap = FALSE,                         # Keep heatmap generation (optional)
                                          max_significance = 0.05,                     # Significance threshold for q-values
                                          standardize = TRUE                           # Disable standardization (optional)
)

oars.mpx.cazy.maaslin.036 = oars.mpx.cazy.maaslin.036$results %>% data.frame() %>% arrange(pval)


# recalculate padj minus diagnosis fixed effect
oars.mpx.cazy.maaslin.036 = subset(oars.mpx.cazy.maaslin.036, metadata == "oars.days") %>%
  mutate(padj = p.adjust(pval, method="BH"))

# fix feature names
oars.mpx.cazy.maaslin.036$feature = oars.mpx.cazy.names$good[match(oars.mpx.cazy.maaslin.036$feature, oars.mpx.cazy.names$feature)]

oars.mpx.cazy.maaslin.036.volcano = ggplot(oars.mpx.cazy.maaslin.036,
                                           aes(x=coef, y=(padj)))+
  geom_point(shape=21, aes(fill=coef))+
  geom_hline(yintercept=(0.2), linetype=2, alpha=0.5)+
  geom_vline(xintercept=c(-1,1), linetype=2, alpha=0.5)+
  scale_fill_gradient2(low="blue", high="red")+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.05,0.1, 0.20, 0.5, 1))+
  # ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, feature, NA)), size=2.5)+  
  facet_wrap(~"Treatment: CAZy")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text=element_text(size=10))+
  labs(x="Adjusted Coefficient",
       y="FDR")
oars.mpx.cazy.maaslin.036.volcano
subset(oars.mpx.cazy.maaslin.036, padj < 0.20) %>% nrow()



oars.mpx.cazy.mat[rownames(oars.mpx.cazy.mat) %in% subset(metadata.oars.stool, timing %in% c("6M", "9M", "12M")& compliant==TRUE)$standard.name,] %>% nrow()
# n = 26 treatment, compliant

oars.mpx.cazy.maaslin.6912 = Maaslin2::Maaslin2(input_data = oars.mpx.cazy.mat,
                                           input_metadata = subset(metadata.oars.stool, timing %in% c("6M", "9M", "12M")& compliant==TRUE),
                                           output = "~/Downloads",
                                           fixed_effects = c("oars.days"),  # Example fixed effects
                                           random_effects = c("HM", "plate"),       # Example random effects
                                           normalization = "NONE",                       # Total Sum Scaling normalization
                                           transform = "LOG",                           # Log transformation
                                           analysis_method = "LM",                      # Linear model
                                           plot_scatter = FALSE,                        # Disable scatterplot generation
                                           plot_heatmap = FALSE,                         # Keep heatmap generation (optional)
                                           max_significance = 0.05,                     # Significance threshold for q-values
                                           standardize = TRUE                           # Disable standardization (optional)
)

oars.mpx.cazy.maaslin.6912 = oars.mpx.cazy.maaslin.6912$results %>% data.frame() %>% arrange(pval)

# recalculate padj minus diagnosis fixed effect
oars.mpx.cazy.maaslin.6912 = subset(oars.mpx.cazy.maaslin.6912, metadata == "oars.days") %>%
  mutate(padj = p.adjust(pval, method="BH"))

# fix feature names
oars.mpx.cazy.maaslin.6912$feature = oars.mpx.cazy.names$good[match(oars.mpx.cazy.maaslin.6912$feature, oars.mpx.cazy.names$feature)]

oars.mpx.cazy.maaslin.6912.volcano = ggplot(oars.mpx.cazy.maaslin.6912,
                                            aes(x=coef, y=(padj)))+
  geom_point(shape=21, aes(fill=coef))+
  geom_hline(yintercept=(0.2), linetype=2, alpha=0.5)+
  geom_vline(xintercept=c(-1,1), linetype=2, alpha=0.5)+
  scale_fill_gradient2(low="blue", high="red")+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.05,0.1, 0.20, 0.5, 1))+
  # ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, feature, NA)), size=2.5)+  
  facet_wrap(~"Washout: CAZy")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text=element_text(size=10))+
  labs(x="Adjusted Coefficient",
       y="FDR")
oars.mpx.cazy.maaslin.6912.volcano
# 

oars.mpx.cazy.maaslin.036.volcano+
  oars.mpx.cazy.maaslin.6912.volcano

# :: MPX Starch -----------------------------------------------------

# stats
stats.starch.036 = lmerTest::lmer(scale(log2(starch)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM) + (1|plate),
                                        subset(metadata.oars.stool.asv, timing %in% c("0M", "3M", "6M")& 
                                                 compliant == TRUE)) %>%
  summary() %>% coef()

stats.starch.6912 = lmerTest::lmer(scale(log2(starch)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM) + (1|plate),
                                         subset(metadata.oars.stool.asv, timing %in% c("6M", "9M", "12M")& 
                                                  compliant == TRUE)) %>%
  summary() %>% coef()

# GAM
stats.starch.gam = mgcv::gam(scale(log2(starch)) ~ s(scale(oars.days), k=10) + diagnosis + adj.fiber + s(HM, bs="re") + s(plate, bs="re"),
                            data = subset(metadata.oars.stool.asv, compliant == TRUE & !is.na(starch))%>%
                              mutate(HM = as.factor(HM)),
                            family = gaussian(),
                            method = "REML") %>% summary() 
# very sig

# plot
metadata.oars.stool.starch.plot = ggplot(subset(metadata.oars.stool.asv,compliant==TRUE),
                                               aes(x=oars.days, 
                                                   y=starch))+
  annotate("rect", xmin=0, xmax=max(subset(metadata.oars.stool, oars.on.rs=="onRS")$oars.days),
           ymin=0, ymax=Inf,
           alpha=0.2, fill="salmon") +
  scale_y_log10()+
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(aes(fill=rs.col, shape=baseline), size=2)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  geom_smooth(color="black")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12))+
  geom_text(data=pval.positioner("starch",
                                 subset(metadata.oars.stool.asv, compliant == TRUE),
                                 stats.starch.036,
                                 stats.starch.6912,
                                 yadjust = 1.2), aes(x=xval, y=yval, label=label), vjust=1,size=4)+
  labs(x="Days since starting RS", y="Starch CAZy Intensity")
metadata.oars.stool.starch.plot

# :: MPX Mucin -----------------------------------------------------

# stats
stats.mucin.036 = lmerTest::lmer(scale(log2(mucin)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM) + (1|plate),
                                        subset(metadata.oars.stool.asv, timing %in% c("0M", "3M", "6M")& 
                                                 compliant == TRUE)) %>%
  summary() %>% coef()

stats.mucin.6912 = lmerTest::lmer(scale(log2(mucin)) ~ scale(oars.days) + diagnosis + adj.fiber+ (1|HM) + (1|plate),
                                         subset(metadata.oars.stool.asv, timing %in% c("6M", "9M", "12M")& 
                                                  compliant == TRUE)) %>%
  summary() %>% coef()

# GAM
stats.mucin.gam = mgcv::gam(scale(log10(mucin)) ~ s(scale(oars.days), k=10) + diagnosis + adj.fiber + s(HM, bs="re") + s(plate, bs="re"),
                            data = subset(metadata.oars.stool.asv, compliant == TRUE & !is.na(mucin))%>%
                              mutate(HM = as.factor(HM)),
                            family = gaussian(),
                            method = "REML") %>% summary() 
# highly individual specific

# plot
metadata.oars.stool.mucin.plot = ggplot(subset(metadata.oars.stool.asv,compliant==TRUE),
                                               aes(x=oars.days, 
                                                   y=mucin))+
  annotate("rect", xmin=0, xmax=max(subset(metadata.oars.stool, oars.on.rs=="onRS")$oars.days),
           ymin=-Inf, ymax=Inf,
           alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(aes(fill=rs.col, shape=baseline), size=2)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  geom_smooth(color="black")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12))+
  geom_text(data=pval.positioner("mucin",
                                 subset(metadata.oars.stool.asv, compliant == TRUE),
                                 stats.mucin.036,
                                 stats.mucin.6912,
                                 yadjust = 1.15), aes(x=xval, y=yval, label=label), vjust=1,size=4)+
  labs(x="Days since starting RS", y="Mucin CAZy Intensity")
metadata.oars.stool.mucin.plot

# :: MPX Starch:Mucin -----------------------------------------------------

# stats
stats.starch.mucin.036 = lmerTest::lmer(scale((starch.mucin)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM) + (1|plate),
                                   subset(metadata.oars.stool.asv, timing %in% c("0M", "3M", "6M")& 
                                            compliant == TRUE)) %>%
  summary() %>% coef()

stats.starch.mucin.6912 = lmerTest::lmer(scale((starch.mucin)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM) + (1|plate),
                                    subset(metadata.oars.stool.asv, timing %in% c("6M", "9M", "12M")& 
                                             compliant == TRUE)) %>%
  summary() %>% coef()

# GAM
stats.starch.mucin.gam = mgcv::gam(scale((starch.mucin)) ~ s(scale(oars.days), k=10) + diagnosis + adj.fiber + s(HM, bs="re") + s(plate, bs="re"),
                            data = subset(metadata.oars.stool.asv, compliant == TRUE & !is.na(starch.mucin))%>%
                              mutate(HM = as.factor(HM)),
                            family = gaussian(),
                            method = "REML") %>% summary() 
# sig

# plot
metadata.oars.stool.starch.mucin.plot = ggplot(subset(metadata.oars.stool.asv,compliant==TRUE),
                                              aes(x=oars.days, 
                                                  y=starch.mucin))+
  annotate("rect", xmin=0, xmax=max(subset(metadata.oars.stool, oars.on.rs=="onRS")$oars.days),
           ymin=-Inf, ymax=Inf,
           alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(aes(fill=rs.col, shape=baseline), size=2)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  geom_smooth(color="black")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12))+
  geom_text(data=pval.positioner("starch.mucin",
                                 subset(metadata.oars.stool.asv, compliant == TRUE),
                                 stats.starch.mucin.036,
                                 stats.starch.mucin.6912,
                                 yadjust = 1.15), aes(x=xval, y=yval, label=label), vjust=1,size=4)+
  labs(x="Days since starting RS", y="Starch:Mucin CAZy Log2 Ratio")
metadata.oars.stool.starch.mucin.plot


# :: MBX Filter --------------------------------------------------------------

# for the sake of interpretability, only examine annotated features here

# load annotated
oars.mbx.raw.mat = readRDS("./metabolomics/2026_01_15_oars_mbx_annotated.Rds")

# filter to compliant
oars.mbx.raw.mat.1 = oars.mbx.raw.mat[rownames(oars.mbx.raw.mat) %in% subset(metadata.oars.stool, compliant==T)$standard.name,]

# 80% prevalence; PER TIME PHASE (preRS, onRS, postRS)
dim(oars.mbx.raw.mat.1)
oars.mbx.raw.mat.filt.1 = oars.mbx.raw.mat.1 %>% as.matrix() %>%
  reshape2::melt()
colnames(oars.mbx.raw.mat.filt.1)[1] = "standard.name"
oars.mbx.raw.mat.filt.1 = merge(oars.mbx.raw.mat.filt.1,
                              metadata.oars.stool[,c("standard.name", "oars.on.rs")], by="standard.name")
oars.mbx.raw.mat.filt.1$value = ifelse(oars.mbx.raw.mat.filt.1$value >0, 1, 0)
# calc prevalence per phase (pre, on, post)
oars.mbx.raw.mat.filt.1 = oars.mbx.raw.mat.filt.1 %>%
  group_by(Var2, oars.on.rs) %>%
  summarise(prevalence = mean(value, na.rm = TRUE)) %>% # cool, use mean of presence to calculate prevalence
  ungroup() %>%
  subset(prevalence >= .80)
# keep features with 80% prevalence in at least one of the groups
length(unique(oars.mbx.raw.mat.filt.1$Var2))
# n = 203 features
# apply filter
oars.mbx.raw.mat.filt.1 = oars.mbx.raw.mat.1[,colnames(oars.mbx.raw.mat.1) %in% unique(oars.mbx.raw.mat.filt.1$Var2)]
dim(oars.mbx.raw.mat.filt.1)

# create feature map
# fix feature names
oars.mbx.names = data.frame(feature = make.names(colnames(oars.mbx.raw.mat.filt.1)),
                            good = colnames(oars.mbx.raw.mat.filt.1))

# :: MBX PCA --------------------------------------------------------------

# log transform with pseudo
oars.mbx.pca = log2(oars.mbx.raw.mat.filt.1+(min(oars.mbx.raw.mat.filt.1[oars.mbx.raw.mat.filt.1!=0])/2))
# PCA
oars.mbx.pca = prcomp((oars.mbx.pca), scale=T)
oars.mbx.pca.df = oars.mbx.pca$x[,c(1,2)] %>% data.frame() %>%
  rownames_to_column("standard.name")
oars.mbx.pca.df = merge(oars.mbx.pca.df,
                                  metadata.oars.stool, by="standard.name")
oars.mbx.pca.var <- (oars.mbx.pca$sdev^2 / sum(oars.mbx.pca$sdev^2)) * 100

rownames(oars.mbx.pca.df) = oars.mbx.pca.df$standard.name
nrow(oars.mbx.pca.df)

# permanova
set.seed(25)
t1 = Sys.time()
oars.mbx.pca.permanova = vegan::adonis2(dist(oars.mbx.pca$x) ~ on.rs + adj.fiber + HM,
                                        oars.mbx.pca.df %>% mutate(on.rs = as.factor(ifelse(oars.on.rs=="onRS", "onRS", "offRS"))),
                                                  strata = oars.mbx.pca.df$HM,
                                                  by="margin")
t2 = Sys.time()
t2 - t1

oars.mbx.pca.plot <- ggplot(
  data=oars.mbx.pca.df %>% group_by(HM) %>% arrange(stool_date_rec_v2) %>%
    mutate(rs.col = factor(RS_Name, levels=rs.names)), 
  aes(x=PC1, y=PC2))+
  geom_path(aes(group=HM), color="black", linetype=2, alpha=0.5, linewidth=0.3) + 
  geom_point(aes(fill=rs.col, shape=baseline), size=2)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  # annotate(geom="text",
  #          x=-9,
  #          y=-6,
  #          label = paste(paste("Treatment\nR²: ", round(data.frame(oars.mbx.pca.permanova)[1,3], 3)*100, "%",
  #                              "\n p: ", round(data.frame(oars.mbx.pca.permanova)[1,5], 3), sep="")),
  #          size=3.3)+
  facet_wrap(~"Metabolite PCA")+
  #geom_text(aes(label=rs.col, color=rs.col))+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text=element_text(size=10))+
  labs(x=paste("PC1: ", round((oars.mbx.pca.var)[1], digits=2), "%", sep=""), 
       y=paste("PC2: ", round((oars.mbx.pca.var)[2], digits=2), "%", sep=""))
oars.mbx.pca.plot


# :: MBX Maaslin2 ----------------------------------------------------

oars.mbx.raw.mat.filt.1[rownames(oars.mbx.raw.mat.filt.1) %in% subset(metadata.oars.stool, timing %in% c("0M", "3M", "6M")& compliant==TRUE)$standard.name,] %>% nrow()
# n = 32 treatment, compliant

oars.mbx.maaslin.036 = Maaslin2::Maaslin2(input_data = (oars.mbx.raw.mat.filt.1),
                                               input_metadata = subset(metadata.oars.stool, timing %in% c("0M", "3M", "6M") & compliant==TRUE),
                                               output = "~/Downloads",
                                               fixed_effects = c("oars.days"),  # Example fixed effects
                                               random_effects = c("HM"),       # Example random effects
                                               normalization = "NONE",                       # Total Sum Scaling normalization
                                               transform = "LOG",                           # Log transformation
                                               analysis_method = "LM",                      # Linear model
                                               min_prevalence = 0, # use custom filter
                                               plot_scatter = FALSE,                        # Disable scatterplot generation
                                               plot_heatmap = FALSE,                         # Keep heatmap generation (optional)
                                               max_significance = 0.05,                     # Significance threshold for q-values
                                               standardize = TRUE                           # Disable standardization (optional)
)

oars.mbx.maaslin.036 = oars.mbx.maaslin.036$results %>% data.frame() %>% arrange(pval)

# recalculate padj minus diagnosis fixed effect
oars.mbx.maaslin.036 = subset(oars.mbx.maaslin.036, metadata == "oars.days") %>%
  mutate(padj = p.adjust(pval, method="BH"))

# fix feature names
oars.mbx.maaslin.036$feature = oars.mbx.names$good[match(oars.mbx.maaslin.036$feature, oars.mbx.names$feature)]

oars.mbx.maaslin.036 %>% head(n=10)

oars.mbx.maaslin.036.volcano = ggplot(oars.mbx.maaslin.036,
                                           aes(x=coef, y=(padj)))+
  geom_point(shape=21, aes(fill=coef))+
  geom_hline(yintercept=(0.2), linetype=2, alpha=0.5)+
  geom_vline(xintercept=c(-1,1), linetype=2, alpha=0.5)+
  scale_fill_gradient2(low="blue", high="red")+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.05,0.1, 0.20, 0.5, 1))+
  # ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, feature, NA)), size=2.5)+  
  facet_wrap(~"Treatment: Metabolite")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text=element_text(size=10))+
  labs(x="Adjusted Coefficient",
       y="FDR")
oars.mbx.maaslin.036.volcano
subset(oars.mbx.maaslin.036, padj < 0.20) %>% nrow()


oars.mbx.raw.mat.filt.1[rownames(oars.mbx.raw.mat.filt.1) %in% subset(metadata.oars.stool, timing %in% c("6M", "9M", "12M")& compliant==TRUE)$standard.name,] %>% nrow()
# n = 22 compliant

oars.mbx.maaslin.6912 = Maaslin2::Maaslin2(input_data = oars.mbx.raw.mat.filt.1,
                                                input_metadata = subset(metadata.oars.stool, timing %in% c("6M", "9M", "12M")& compliant==TRUE),
                                                output = "~/Downloads",
                                                fixed_effects = c("oars.days"),  # Example fixed effects
                                                random_effects = c("HM"),       # Example random effects
                                                normalization = "NONE",                       # Total Sum Scaling normalization
                                                transform = "LOG",                           # Log transformation
                                                analysis_method = "LM",                      # Linear model
                                                min_prevalence = 0, # use custom filter
                                                plot_scatter = FALSE,                        # Disable scatterplot generation
                                                plot_heatmap = FALSE,                         # Keep heatmap generation (optional)
                                                max_significance = 0.05,                     # Significance threshold for q-values
                                                standardize = TRUE                           # Disable standardization (optional)
)

oars.mbx.maaslin.6912 = oars.mbx.maaslin.6912$results %>% data.frame() %>% arrange(pval)

# recalculate padj minus diagnosis fixed effect
oars.mbx.maaslin.6912 = subset(oars.mbx.maaslin.6912, metadata == "oars.days") %>%
  mutate(padj = p.adjust(pval, method="BH"))

# fix feature names
oars.mbx.maaslin.6912$feature = oars.mbx.names$good[match(oars.mbx.maaslin.6912$feature, oars.mbx.names$feature)]

oars.mbx.maaslin.6912.volcano = ggplot(oars.mbx.maaslin.6912,
                                            aes(x=coef, y=(padj)))+
  geom_point(shape=21, aes(fill=coef))+
  geom_hline(yintercept=(0.2), linetype=2, alpha=0.5)+
  geom_vline(xintercept=c(-1,1), linetype=2, alpha=0.5)+
  scale_fill_gradient2(low="blue", high="red")+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.05,0.1, 0.20, 0.5, 1))+
 #  ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, feature, NA)), size=2.5)+  
  facet_wrap(~"Washout: Metabolite")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text=element_text(size=10))+
  labs(x="Adjusted Coefficient",
       y="FDR")
oars.mbx.maaslin.6912.volcano
# 

oars.mbx.maaslin.036.volcano+
  oars.mbx.maaslin.6912.volcano


# :: PCA  R2 ------------------------------------------------------------------

# goal: plot Var-Exp and p values for PCA PERMANOVAs
oars.pca.permanova = rbind(
  oars.asv.permanova %>% as.data.frame() %>% mutate(data.type = "ASV") %>% mutate(var = rownames(.)),
  oars.mgx.permanova %>% as.data.frame() %>% mutate(data.type = "Species") %>% mutate(var = rownames(.)),
  oars.mpx.kegg.oars.pca.permanova %>% as.data.frame() %>% mutate(data.type = "Pathway") %>% mutate(var = rownames(.)),
  oars.mpx.cog.oars.pca.permanova%>% as.data.frame() %>% mutate(data.type = "COG")  %>% mutate(var = rownames(.)),
  oars.mpx.cazy.pca.permanova%>% as.data.frame() %>% mutate(data.type = "CAZy") %>% mutate(var = rownames(.)),
  oars.mbx.pca.permanova%>% as.data.frame() %>% mutate(data.type = "Metabolite") %>% mutate(var = rownames(.)))
# remove extra var
oars.pca.permanova = subset(oars.pca.permanova, var %in% c("on.rs", "adj.fiber", "HM"))
oars.pca.permanova$Variable = ifelse(oars.pca.permanova$var == "on.rs", "Treatment",
                                     ifelse(oars.pca.permanova$var == "adj.fiber", "Adjusted Fiber",
                                            "Individual"))
# plot
oars.pca.permanova.plot = ggplot(oars.pca.permanova %>% 
         mutate(data.type = factor(data.type, levels=c("ASV", "Species", "Pathway", "COG", "CAZy", "Metabolite"))) %>%
         mutate(Variable = factor(Variable, levels=rev(c("Individual", "Treatment", "Adjusted Fiber")))) %>%
         mutate(padj = p.adjust(`Pr(>F)`, method="BH")),
       aes(x=R2*100, y=Variable))+
  geom_bar(stat="identity", aes(fill=data.type), color="black")+
  #geom_text(aes(label=ifelse(`Pr(>F)` < 0.05, "*", ""), x=R2*100), size=7, vjust=0.8, nudge_x=2)+
  geom_text(aes(label=ifelse(padj < 0.20, "*", ""), x=R2*100), size=7, alpha=0.5, vjust=0.8, nudge_x=2)+
  geom_text(aes(label=ifelse(padj < 0.05, "*", ""), x=R2*100), size=7, vjust=0.8, nudge_x=2)+
  scale_x_continuous(breaks=c(0, 10, 20, 30, 40, 50, 60))+
  scale_fill_manual(values=omics.colors)+
  theme_classic()+theme(legend.position = "none",
                        strip.text=element_text(size=10),
                        panel.grid.major.x = element_line(color="grey", linewidth=0.2))+
  facet_wrap(~data.type, ncol=1)+
  labs(x="Variance Explained", y="")
oars.pca.permanova.plot
# adjusted p values < 0.20

# does feature space correlate with individualized effect?

oars.pca.permanova.vs.features = subset(oars.pca.permanova, Variable == "Individual") %>% 
  merge(data.frame(data.type = c("ASV", "Species", "Pathway", "COG", "CAZy", "Metabolite"),
           features = c(ncol(oars.asv.data.median.1.filt), # 1165
                        ncol(oars.mgx.taxa.filt.1), # 120
                        ncol(oars.mpx.kegg.mat.filt.1), # 180
                        ncol(oars.mpx.cog.mat.filt.1), #2481
                        ncol(oars.mpx.cazy.mat.filt.1), # 66
                        ncol(oars.mbx.raw.mat.filt.1))) # 203
        , by="data.type") 

cor.test(oars.pca.permanova.vs.features$R2,
         oars.pca.permanova.vs.features$features, method="spearman")
# p = 0.42
ggplot(oars.pca.permanova.vs.features,
       aes(x=features, y=R2))+
  geom_point(aes(shape=data.type, fill=data.type), size=3)+
  geom_smooth(method="lm", color="black")+
  ggpubr::stat_cor(method="spearman")+
  scale_fill_manual(values=omics.colors)+
  scale_shape_manual(values=omics.shapes)+
  theme_classic()+
  labs(x="Feature Space", y="Individual R2", fill="Data type", shape="Data type")
# not at all

# :: Hysteresis? ----------------------------------------------------------

# Make a heatmap splitting main variables into RS1 vs RS2

oars.hysteresis = rbind(
  lmerTest::lmer(scale(log10(but.i)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                 subset(metadata.oars.stool.asv, timing %in% c("0M", "3M") & compliant == TRUE)) %>% summary() %>% coef() %>% as.data.frame() %>% dplyr::select(Estimate, `Pr(>|t|)`) %>% mutate(feature = rownames(.))%>% mutate(time = "RS1", variable = "Butyrogens (%)") %>%  subset(feature == "scale(oars.days)"),
  lmerTest::lmer(scale(log10(but.i)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                 subset(metadata.oars.stool.asv, timing %in% c("3M", "6M") & compliant == TRUE)) %>% summary() %>% coef() %>% as.data.frame() %>% dplyr::select(Estimate, `Pr(>|t|)`) %>% mutate(feature = rownames(.))%>% mutate(time = "RS2", variable = "Butyrogens (%)") %>%  subset(feature == "scale(oars.days)"),
  
  lmerTest::lmer(scale(log10(fcal)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                 subset(metadata.oars.stool, timing %in% c("0M", "3M") & compliant == TRUE)) %>% summary() %>% coef() %>% as.data.frame() %>% dplyr::select(Estimate, `Pr(>|t|)`) %>% mutate(feature = rownames(.))%>% mutate(time = "RS1", variable = "Fecal calprotectin (μg/g)") %>%  subset(feature == "scale(oars.days)"),
  lmerTest::lmer(scale(log10(fcal)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                 subset(metadata.oars.stool, timing %in% c("3M", "6M") & compliant == TRUE)) %>% summary() %>% coef() %>% as.data.frame() %>% dplyr::select(Estimate, `Pr(>|t|)`) %>% mutate(feature = rownames(.))%>% mutate(time = "RS2", variable = "Fecal calprotectin (μg/g)") %>%  subset(feature == "scale(oars.days)"),
  
  lmerTest::lmer(scale(log10(richness)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                 subset(metadata.oars.stool.asv, timing %in% c("0M", "3M") & compliant == TRUE)) %>% summary() %>% coef() %>% as.data.frame() %>% dplyr::select(Estimate, `Pr(>|t|)`) %>% mutate(feature = rownames(.))%>% mutate(time = "RS1", variable = "ASV Richness") %>%  subset(feature == "scale(oars.days)"),
  lmerTest::lmer(scale(log10(richness)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                 subset(metadata.oars.stool.asv, timing %in% c("3M", "6M") & compliant == TRUE)) %>% summary() %>% coef() %>% as.data.frame() %>% dplyr::select(Estimate, `Pr(>|t|)`) %>% mutate(feature = rownames(.))%>% mutate(time = "RS2", variable = "ASV Richness") %>%  subset(feature == "scale(oars.days)"),
  
  lmerTest::lmer(scale((shannon)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                 subset(metadata.oars.stool.asv, timing %in% c("0M", "3M") & compliant == TRUE)) %>% summary() %>% coef() %>% as.data.frame() %>% dplyr::select(Estimate, `Pr(>|t|)`) %>% mutate(feature = rownames(.)) %>% mutate(time = "RS1", variable = "Shannon Diversity") %>%subset(feature == "scale(oars.days)"),
  lmerTest::lmer(scale((shannon)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                 subset(metadata.oars.stool.asv, timing %in% c("3M", "6M") & compliant == TRUE)) %>% summary() %>% coef() %>% as.data.frame() %>% dplyr::select(Estimate, `Pr(>|t|)`) %>% mutate(feature = rownames(.)) %>% mutate(time = "RS2", variable = "Shannon Diversity") %>% subset(feature == "scale(oars.days)"),
  
  lmerTest::lmer(scale((fd)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                 subset(metadata.oars.stool.asv, timing %in% c("0M", "3M") & compliant == TRUE)) %>% summary() %>% coef() %>% as.data.frame() %>% dplyr::select(Estimate, `Pr(>|t|)`) %>% mutate(feature = rownames(.)) %>% mutate(time = "RS1", variable = "Functional Redundancy") %>% subset(feature == "scale(oars.days)"),
  lmerTest::lmer(scale((fd)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                 subset(metadata.oars.stool.asv, timing %in% c("3M", "6M") & compliant == TRUE)) %>% summary() %>% coef() %>% as.data.frame() %>% dplyr::select(Estimate, `Pr(>|t|)`) %>% mutate(feature = rownames(.)) %>% mutate(time = "RS2", variable = "Functional Redundancy") %>% subset(feature == "scale(oars.days)"),
  
  lmerTest::lmer(scale((stool_water_perc)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                 subset(metadata.oars.stool, timing %in% c("0M", "3M") & compliant == TRUE)) %>% summary() %>% coef() %>% as.data.frame() %>% dplyr::select(Estimate, `Pr(>|t|)`) %>% mutate(feature = rownames(.)) %>% mutate(time = "RS1", variable = "Stool Moisture (%)") %>% subset(feature == "scale(oars.days)"),
  lmerTest::lmer(scale((stool_water_perc)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                 subset(metadata.oars.stool, timing %in% c("3M", "6M") & compliant == TRUE)) %>% summary() %>% coef() %>% as.data.frame() %>% dplyr::select(Estimate, `Pr(>|t|)`) %>% mutate(feature = rownames(.)) %>% mutate(time = "RS2", variable = "Stool Moisture (%)") %>% subset(feature == "scale(oars.days)"),
  
  lmerTest::lmer(scale((between.beta)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                 subset(beta.trajectory.data, timing %in% c("0M", "3M") & compliant == TRUE)) %>% summary() %>% coef() %>% as.data.frame() %>% dplyr::select(Estimate, `Pr(>|t|)`) %>% mutate(feature = rownames(.))%>% mutate(time = "RS1", variable = "Mean Bray-Curtis Dissimilarity") %>% subset(feature == "scale(oars.days)"),
  lmerTest::lmer(scale((between.beta)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                 subset(beta.trajectory.data, timing %in% c("3M", "6M") & compliant == TRUE)) %>% summary() %>% coef() %>% as.data.frame() %>% dplyr::select(Estimate, `Pr(>|t|)`) %>% mutate(feature = rownames(.))%>% mutate(time = "RS2", variable = "Mean Bray-Curtis Dissimilarity") %>% subset(feature == "scale(oars.days)"),
  
  lmerTest::lmer(scale((load.asv)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                 subset(metadata.oars.stool.asv, timing %in% c("0M", "3M") & compliant == TRUE)) %>% summary() %>% coef() %>% as.data.frame() %>% dplyr::select(Estimate, `Pr(>|t|)`) %>% mutate(feature = rownames(.))%>% mutate(time = "RS1", variable = "Microbial Load (Log10 Predicted)") %>% subset(feature == "scale(oars.days)"),
  lmerTest::lmer(scale((load.asv)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                 subset(metadata.oars.stool.asv, timing %in% c("3M", "6M") & compliant == TRUE)) %>% summary() %>% coef() %>% as.data.frame() %>% dplyr::select(Estimate, `Pr(>|t|)`) %>% mutate(feature = rownames(.))%>% mutate(time = "RS2", variable = "Microbial Load (Log10 Predicted)") %>% subset(feature == "scale(oars.days)"),
  
  lmerTest::lmer(scale(log10(starch)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                 subset(metadata.oars.stool.asv, timing %in% c("0M", "3M") & compliant == TRUE)) %>% summary() %>% coef() %>% as.data.frame() %>% dplyr::select(Estimate, `Pr(>|t|)`) %>% mutate(feature = rownames(.))%>% mutate(time = "RS1", variable = "Starch CAZy Intensity") %>%  subset(feature == "scale(oars.days)"),
  lmerTest::lmer(scale(log10(starch)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                 subset(metadata.oars.stool.asv, timing %in% c("3M", "6M") & compliant == TRUE)) %>% summary() %>% coef() %>% as.data.frame() %>% dplyr::select(Estimate, `Pr(>|t|)`) %>% mutate(feature = rownames(.)) %>% mutate(time = "RS2", variable = "Starch CAZy Intensity") %>% subset(feature == "scale(oars.days)"),
  
  lmerTest::lmer(scale(log10(mucin)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                 subset(metadata.oars.stool.asv, timing %in% c("0M", "3M") & compliant == TRUE)) %>% summary() %>% coef() %>% as.data.frame() %>% dplyr::select(Estimate, `Pr(>|t|)`) %>% mutate(feature = rownames(.)) %>% mutate(time = "RS1", variable = "Mucin CAZy Intensity") %>% subset(feature == "scale(oars.days)"),
  lmerTest::lmer(scale(log10(mucin)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                 subset(metadata.oars.stool.asv, timing %in% c("3M", "6M") & compliant == TRUE)) %>% summary() %>% coef() %>% as.data.frame() %>% dplyr::select(Estimate, `Pr(>|t|)`) %>% mutate(feature = rownames(.)) %>% mutate(time = "RS2", variable = "Mucin CAZy Intensity") %>% subset(feature == "scale(oars.days)"),
  
  lmerTest::lmer(scale((starch.mucin)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                 subset(metadata.oars.stool.asv, timing %in% c("0M", "3M") & compliant == TRUE)) %>% summary() %>% coef() %>% as.data.frame() %>% dplyr::select(Estimate, `Pr(>|t|)`) %>% mutate(feature = rownames(.)) %>% mutate(time = "RS1", variable = "Starch:Mucin CAZy Log2 Ratio") %>% subset(feature == "scale(oars.days)"),
  lmerTest::lmer(scale((starch.mucin)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                 subset(metadata.oars.stool.asv, timing %in% c("3M", "6M") & compliant == TRUE)) %>% summary() %>% coef() %>% as.data.frame() %>% dplyr::select(Estimate, `Pr(>|t|)`) %>% mutate(feature = rownames(.)) %>% mutate(time = "RS2", variable = "Starch:Mucin CAZy Log2 Ratio") %>% subset(feature == "scale(oars.days)")
)


oars.hysteresis.plot = ggplot(oars.hysteresis %>% 
         mutate(time = factor(time, levels=c("RS2", "RS1"))) %>%
         mutate(variable = factor(variable, levels=c("Butyrogens (%)", "Fecal calprotectin (μg/g)", "ASV Richness", "Shannon Diversity",
                                                     "Functional Redundancy", "Mean Bray-Curtis Dissimilarity", "Stool Moisture (%)",
                                                     "Microbial Load (Log10 Predicted)", "Starch CAZy Intensity", "Mucin CAZy Intensity", "Starch:Mucin CAZy Log2 Ratio"))),
       aes(x=variable, y=time))+
  geom_tile(aes(fill=scale(Estimate)), color="white")+
  geom_text(aes(label=ifelse(`Pr(>|t|)` < 0.05, "*", NA)), size=8, vjust=0.8, color="white")+
  scale_fill_gradient2(low="blue", mid="white", high="red")+
  theme_classic()+theme(#legend.position="none",
                        axis.text.x=element_text(angle=45, hjust=1))+
  labs(x="", y="RS Phase", fill="Coefficient")
oars.hysteresis.plot
  

# :: GAM Variables // defunct ------------------------------------------------------------

# print table of GAM results

oars.gam.results = rbind(
  stats.but.i.gam$s.table[,4],
  stats.fcal.gam$s.table[,4],
  stats.water.gam$s.table[,4],
  stats.load.asv.gam$s.table[,4],
  stats.fd.gam$s.table[,4],
  stats.beta.between.gam$s.table[,4],
  stats.starch.gam$s.table[,4],
  stats.mucin.gam$s.table[,4],
  stats.starch.mucin.gam$s.table[,4]) %>% as.data.frame() %>%
  mutate(feature = factor(c("Butyrogens", "Fecal Calprotectin",
                      "Stool Moisture", "Microbial Load","Functional Redundancy",
                      "Bray-Curtis Dissimilarity", "Starch CAZy", "Mucin CAZy", "Starch:Mucin"), levels=rev(c("Butyrogens", "Fecal Calprotectin",
                                                                                                           "Stool Moisture", "Microbial Load","Functional Redundancy",
                                                                                                           "Bray-Curtis Dissimilarity", "Starch CAZy", "Mucin CAZy", "Starch:Mucin")))) %>%
  rename("Days since starting RS" = "s(scale(oars.days))", 
         "Patient" = "s(HM)") %>%
  reshape2::melt() %>%
  mutate(padj = p.adjust(value, method="BH"))%>%
  ggplot(aes(x=variable, y=feature))+
  geom_tile(aes(fill=scale(padj)), color="white")+
  geom_text(aes(label=round(padj, digits=3), color=ifelse(padj < 0.20, "1", "2")))+
  scale_fill_gradient2(low="blue", high="red")+
  scale_color_manual(values=c("1" = "white", "2" = "black"))+
  theme_classic()+
  guides(fill="none", color="none")+
  facet_wrap(~"GAM FDR")+
  labs(x="", y="", fill="p-value")
oars.gam.results
# fill and color = adjusted p value (blue = sig, padj < 0.20)
# left tiles = time (RS)
# right tiles = individualized effect




# >> 2.5. PostHoc pH Response ----------------------------------------------------------

oars.rapidaim.scores = readRDS("./2025_06_09_oars_scores_ph.Rds")
oars.rapidaim.scores = oars.rapidaim.scores[,!colnames(oars.rapidaim.scores) %in% c("Replicate", "pH")] %>% distinct()

# subset RS to the one that was selected for the microbiome
oars.rapidaim.scores = subset(oars.rapidaim.scores, rs.selected == RS_Name)

# append ph scores to metadata
metadata.oars.stool.ph = merge(metadata.oars.stool, oars.rapidaim.scores[,c("HM","rs.selected", "timing","delta.ph")], by=c("HM", "timing"))
metadata.oars.stool.asv.ph = merge(metadata.oars.stool.asv, oars.rapidaim.scores[,c("HM","rs.selected", "timing","delta.ph")], by=c("HM", "timing"))

# double up variables
metadata.oars.stool.double = doubleup(metadata.oars.stool.ph, delta=T)
metadata.oars.stool.asv.double = doubleup(metadata.oars.stool.asv.ph, delta=T)

# noncompliant are removed within doubleup function

# add signifier for removing NA data
metadata.oars.stool.double$hm_phase = paste(metadata.oars.stool.double$HM, 
                                            metadata.oars.stool.double$phase, sep="_")
metadata.oars.stool.asv.double$hm_phase = paste(metadata.oars.stool.asv.double$HM, 
                                                metadata.oars.stool.asv.double$phase, sep="_")
# remove samples missing paired value
metadata.oars.stool.double = subset(metadata.oars.stool.double, !hm_phase %in% c("HM0924_rs2", "HM0932_rs2", "HM0759_rs2"))
metadata.oars.stool.asv.double = subset(metadata.oars.stool.asv.double, !hm_phase %in% c("HM0924_rs2", "HM0932_rs2", "HM0759_rs2"))

# factor RS
oars.rapidaim.scores$RS_Name = factor(oars.rapidaim.scores$RS_Name, levels=rs.names)

oars.ph.zscore.scatterplot = ggplot(subset(oars.rapidaim.scores, timing %in% c("0M", "3M")),
       aes(x=Z_score, y=delta.ph))+
  geom_path(aes(group=HM), linetype=2, alpha=0.5)+
  geom_smooth(method="lm", color="black")+
  geom_point(shape=21, 
             aes(fill=RS_Name), size=2.5)+
  ggpubr::stat_cor(method="spearman", label.x.npc="center")+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12))+
  labs(y="Δ pH of Selected RS",
       x="Butyrogen Z-Score")
oars.ph.zscore.scatterplot

oars.ph.rs.boxplot = ggplot(subset(oars.rapidaim.scores, timing %in% c("0M", "3M")) %>%
                              # remove unused HM0759
                              subset(!sample %in% c("HM0759.04", "HM0932.02", "HM0924.03R")) %>%
                              mutate(HM_time = paste(HM, timing, sep="_")) %>%
                              subset(HM_time %in% paste(subset(metadata.oars.stool, compliant == T)$HM, 
                                                        subset(metadata.oars.stool, compliant == T)$timing, sep="_"))%>%
                              mutate(rs.selected = factor(rs.selected, levels=rs.names))%>%
                              dplyr::select(HM_time, rs.selected, delta.ph)%>%
                              # add blanks
                              rbind(data.frame(HM_time = NA, rs.selected = "BobsRedMill", delta.ph=NA),
                                    data.frame(HM_time = NA, rs.selected = "HiMaize260", delta.ph=NA),
                                    data.frame(HM_time = NA, rs.selected = "Versafibe1490", delta.ph=NA)),
       aes(x=rs.selected, y=delta.ph))+
  geom_boxplot(width=0.5, alpha=0.5, color="black")+
  geom_point(shape=21, 
             aes(fill=rs.selected), size=3)+
  geom_hline(yintercept=mean(metadata.oars.stool.double$delta.ph), 
             linetype=2, color="red")+
  #ggrepel::geom_text_repel(aes(label=HM_time),size=3)+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  ylim(c(min(metadata.oars.stool.double$delta.ph)-0.25,
         max(metadata.oars.stool.double$delta.ph)+0.25))+
  theme_classic()+theme(legend.position="none",
                        panel.grid.minor=element_blank(),
                        panel.grid.major.x=element_blank(),
                        axis.text.x=element_text(angle=45, hjust=1))+
  labs(y="Δ pH of Selected RS",
       x="")
oars.ph.rs.boxplot


# :: Response Threshold --------------------------------------------------------

# assign response group
response = mean(distinct(metadata.oars.stool.double[,c("HM", "timing", "delta.ph")])$delta.ph)
# these are doubled, but the mean is the same

metadata.oars.stool.double$response = ifelse(metadata.oars.stool.double$delta.ph > mean(metadata.oars.stool.double$delta.ph), "low", "high")
metadata.oars.stool.asv.double$response = ifelse(metadata.oars.stool.asv.double$delta.ph > mean(metadata.oars.stool.double$delta.ph), "low", "high")

table(distinct(metadata.oars.stool.double[,c("HM", "phase","response")])$response)
# 8 and 13

# what is the sensitivity range? (used for LSARP-CD validation)
response.sens.1 = subset(metadata.oars.stool.double, response == "low")$delta.ph %>% min()
response.sens.2 = subset(metadata.oars.stool.double, response == "high")$delta.ph %>% max()
mean(metadata.oars.stool.double$delta.ph)
# -1.18 to -1.52 (mean -1.27)


oars.ph.bidist.plot.2 = ggplot((metadata.oars.stool.double[,c("delta.ph", "response", "HM")] %>% distinct()),
                               aes(x = delta.ph, y=0.1))+
  coord_flip()+
  annotate(geom="rect", xmin=-Inf, xmax=Inf,
           ymin=-Inf, ymax=Inf,
           alpha=1, fill="white") +
  ggridges::geom_density_ridges(fill="white",scale = 0.5, alpha=0)+
  geom_boxplot(data=(metadata.oars.stool.double[,c("delta.ph", "response", "HM")] %>% distinct())%>%
                 subset(response == "high"),aes(x = delta.ph, y=0), width=0.2, alpha=0.5)+
  ggbeeswarm::geom_beeswarm(data=(metadata.oars.stool.double[,c("delta.ph", "response", "HM")] %>% distinct())%>%
               subset(response == "high"), aes(x = delta.ph, y=0, fill=response), shape=21, size=3)+
  geom_boxplot(data=(metadata.oars.stool.double[,c("delta.ph", "response", "HM")] %>% distinct())%>%
               subset(response == "low"),aes(x = delta.ph, y=0), width=0.2, alpha=0.5)+
  ggbeeswarm::geom_beeswarm(data=(metadata.oars.stool.double[,c("delta.ph", "response", "HM")] %>% distinct())%>%
                              subset(response == "low"), aes(x = delta.ph, y=0, fill=response), shape=21, size=3)+
  geom_vline(xintercept=mean(metadata.oars.stool.double$delta.ph), 
             linetype=2, color="red")+
  scale_x_continuous(position = "top",
                     limits = c(min(metadata.oars.stool.double$delta.ph)-0.25,
                                max(metadata.oars.stool.double$delta.ph)+0.25))+
  #ggrepel::geom_text_repel(aes(label=HM), size=2)+
  theme_classic()+theme(legend.position="none",
                        panel.grid.minor=element_blank(),
                        panel.grid.major.x=element_blank(),
                        axis.text.x=element_blank(),
                        axis.ticks.x=element_blank())+
  labs(#x="Δ pH of Selected RS",
    x="",
       y="")
oars.ph.bidist.plot.2

# show groups per HM
oars.ph.score.change.plot = ggplot(metadata.oars.stool.double %>%
         mutate(Phase = ifelse(phase == "rs1", "RS 1", "RS 2")) %>%
         mutate(ordering = ifelse(phase == "rs1" & response == "high", 2,
                                  ifelse(phase == "rs1" & response == "low", 1, 
                                         ifelse(phase == "rs2" & response == "high", 1.5,
                                                ifelse(phase == "rs2" & response == "low", 1, 0))))) %>%
           mutate(Response = ifelse(response == "low", "Weak", "Strong"))%>%
         group_by(HM) %>%
         mutate(ordering = sum(ordering)),
       aes(x=Phase, y=reorder(HM, ordering)))+
  geom_path(aes(group=HM), linetype=2, alpha=0.5)+
  geom_point(shape=21, aes(fill = Response), size=4)+
  theme_minimal()+
  labs(x="", y="", fill="Response")
oars.ph.score.change.plot

# V2
(oars.ph.rs.boxplot+
    oars.ph.bidist.plot.2+oars.ph.score.change.plot)+
  patchwork::plot_layout(widths=c(2,1,1))

((oars.ph.rs.boxplot+
    oars.ph.bidist.plot.2)+
  patchwork::plot_layout(widths=c(2,1))) %>%
  ggsave(filename="../For others/2025_08_08_gbm_figure_2_PD.pdf",
         width=8, height=4, device = cairo_pdf)



# >>> 3. RESPONDERS -----------------------------------------------------------

# :: Fecal Calprotectin ---------------------------------------------------

# note: keep only paired datapoints 
stats.fcal.ph.low = lmerTest::lmer(scale(log10(fcal)) ~ reltiming + phase + diagnosis + adj.fiber + (1|HM), 
                                   subset(metadata.oars.stool.double, response == "low")) %>% summary() %>% coef()
stats.fcal.ph.high = lmerTest::lmer(scale(log10(fcal)) ~ reltiming + phase + diagnosis + adj.fiber + (1|HM), 
                                   subset(metadata.oars.stool.double, response == "high")) %>% summary() %>% coef()
stats.fcal.ph.interact = lmerTest::lmer(scale(log10(fcal)) ~ reltiming*response + phase + diagnosis+adj.fiber + (1|HM),
                                        metadata.oars.stool.double) %>% summary() %>% coef()
# fcal is not significantly reduced in Weak response Pre
lmerTest::lmer(scale(log10(fcal)) ~ response + phase + diagnosis + adj.fiber + (1|HM), 
               subset(metadata.oars.stool.double, reltiming == "pre")) %>% summary() %>% coef()

# fcal remains significant after adjusting for stool_water_perc
lmerTest::lmer(scale(log10(fcal)) ~ reltiming + phase + (stool_water_perc) + diagnosis + adj.fiber + (1|HM), 
               subset(metadata.oars.stool.double, response == "high")) %>% summary() %>% coef()
# p = 0.0450

# and when explicitly normalizing value to dry weight
lmerTest::lmer(scale(log10(fcal / (1-stool_water_perc))) ~ reltiming + phase + diagnosis + adj.fiber + (1|HM), 
               subset(metadata.oars.stool.double, response == "high")) %>% summary() %>% coef()
# p = 0.0150

# fcal ug / g stool
# moisture g / g stool
# goal: fecal cal / g dry stool

# 2025_05_28 Note: I spent an hour reanalyzing this data manually 
# (starting a new spreadsheet, collecting HM, RS, delta.ph + response, and fecal cal)
# and the results are consistent

# make data.frame
stats.fcal.ph = data.frame(Response = c("Strong Response", "Weak Response"),
                           pval = c(stats.fcal.ph.high[2,5], stats.fcal.ph.low[2,5]))

metadata.oars.stool.double$RS_Name = factor(metadata.oars.stool.double$RS_Name, levels=rs.names)

metadata.oars.fcal.double.plot = ggplot(subset(metadata.oars.stool.double, !hm_phase %in% c("HM0924_rs2", "HM0932_rs2", "HM0759_rs2")) %>%
                                          mutate(reltiming = factor(ifelse(reltiming == "pre", "Pre", "Post"), levels=c("Pre", "Post"))) %>%
                                          mutate(Response = ifelse(response == "low", "Weak Response", "Strong Response")),
       aes(x=reltiming, y=fcal))+
  scale_y_log10()+
  geom_boxplot(width=0.5)+
  geom_line(aes(group=hm_phase), linetype=2, alpha=0.5, linewidth=0.5)+
  geom_point(shape=21, aes(fill=RS_Name), size=3)+
  geom_hline(yintercept=250, linetype=1, alpha=1, color="red")+
  #ggrepel::geom_label_repel(aes(label=RS_Name, fill=RS_Name),size=4)+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  geom_text(data = stats.fcal.ph%>%
              mutate(Response = ifelse(Response == "Weak Response", "Weak Response", "Strong Response")), 
            x=1.5, y=log10(max(na.omit(metadata.oars.stool.double$fcal))), 
            aes(label=paste("p:", round(pval, digits=3))), size=4.5)+
  facet_wrap(~Response)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="",
       y="Fecal Calprotectin (μg/g)")
metadata.oars.fcal.double.plot
stats.fcal.ph.interact # sig interaction


metadata.oars.fcal.double.plot%>%
  ggsave(filename="../For others/2025_08_08_gbm_figure_3_PD.pdf",
         width=8, height=4, device = cairo_pdf)

# :: Stool Water ----------------------------------------------------------

# delta calculated with fcal

# note: keep only paired datapoints (e.g. remove HM924-High)
stats.water.ph.low = lmerTest::lmer(scale(log10(stool_water_perc)) ~ reltiming + phase + diagnosis+adj.fiber + (1|HM), 
                                   subset(metadata.oars.stool.double, response == "low")) %>% summary() %>% coef()
stats.water.ph.high = lmerTest::lmer(scale(log10(stool_water_perc)) ~ reltiming + phase + diagnosis+adj.fiber + (1|HM), 
                                    subset(metadata.oars.stool.double, response == "high")) %>% summary() %>% coef()
stats.water.ph.interact = lmerTest::lmer(scale(log10(stool_water_perc)) ~ reltiming*response + phase + diagnosis+adj.fiber + (1|HM),
                                         metadata.oars.stool.double) %>% summary() %>% coef()

# make data.frame
stats.water.ph = data.frame(Response = c("Strong Response", "Weak Response"),
                           pval = c(stats.water.ph.high[2,5], stats.water.ph.low[2,5]))

metadata.oars.stool.double$RS_Name = factor(metadata.oars.stool.double$RS_Name, levels=rs.names)

metadata.oars.water.double.plot = ggplot(subset(metadata.oars.stool.double, !hm_phase %in% c("HM0924_rs2", "HM0932_rs2", "HM0759_rs2")) %>%
                                           mutate(reltiming = factor(ifelse(reltiming == "pre", "Pre", "Post"), levels=c("Pre", "Post"))) %>%
                                           mutate(Response = ifelse(response == "low", "Weak Response", "Strong Response")),
       aes(x=reltiming, y=stool_water_perc*100))+
  geom_boxplot(width=0.5)+
  geom_line(aes(group=hm_phase), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(shape=21, aes(fill=RS_Name), size=2)+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  geom_text(data = stats.water.ph, x=1.5, y=(max(na.omit(metadata.oars.stool.double$stool_water_perc)))*100, 
            aes(label=paste("p:", round(pval, digits=3))), size=3.5)+
  facet_wrap(~Response)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="",
       y="Stool Moisture (%)")
metadata.oars.water.double.plot
stats.water.ph.interact # no interaction


# :: Microbial Load (ASV) -------------------------------------------------------

# run stats (already log transformed)
stats.load.asv.ph.low = lmerTest::lmer(scale(load.asv) ~ reltiming + phase + diagnosis+adj.fiber + (1|HM), 
                                       subset(metadata.oars.stool.asv.double, response == "low")) %>% summary() %>% coef()
stats.load.asv.ph.high = lmerTest::lmer(scale(load.asv) ~ reltiming + phase + diagnosis+ adj.fiber + (1|HM), 
                                        subset(metadata.oars.stool.asv.double, response == "high")) %>% summary() %>% coef()
stats.load.asv.ph.interact = lmerTest::lmer(scale(load.asv) ~ reltiming*response + phase +diagnosis+ adj.fiber + (1|HM), metadata.oars.stool.asv.double) %>% summary() %>% coef()
# make data.frame
stats.load.asv.ph = data.frame(Response = c("Strong Response", "Weak Response"),
                            pval = c(stats.load.asv.ph.high[2,5], stats.load.asv.ph.low[2,5]))
metadata.oars.stool.asv.double$RS_Name = factor(metadata.oars.stool.asv.double$RS_Name, levels=rs.names)

metadata.oars.load.asv.double.plot = ggplot(subset(metadata.oars.stool.asv.double, !hm_phase %in% c("HM0924_rs2", "HM0932_rs2", "HM0759_rs2")) %>%
                                          mutate(reltiming = factor(ifelse(reltiming == "pre", "Pre", "Post"), levels=c("Pre", "Post"))) %>%
                                            mutate(Response = ifelse(response == "low", "Weak Response", "Strong Response")),
                                         aes(x=reltiming, y=load.asv))+
  geom_boxplot(width=0.5)+
  geom_line(aes(group=hm_phase), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(shape=21, aes(fill=RS_Name), size=2)+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  geom_text(data = stats.load.asv.ph, x=1.5, y=(max(na.omit(metadata.oars.stool.asv.double$load.asv))), 
            aes(label=paste("p =", round(pval, digits=3))),size=3.5)+
  facet_wrap(~Response,)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="",
       y="Microbial Load") # shortening title here so it fits in the plot
metadata.oars.load.asv.double.plot
stats.load.asv.ph.interact # sig interaction


# :: Microbial Load (MGX) -------------------------------------------------------

# run stats
stats.load.mgx.ph.low = lmerTest::lmer(scale(load.mgx) ~ reltiming + phase + diagnosis+adj.fiber + (1|HM), 
                                       subset(metadata.oars.stool.asv.double, response == "low")) %>% summary() %>% coef()
stats.load.mgx.ph.high = lmerTest::lmer(scale(load.mgx) ~ reltiming + phase + diagnosis+ adj.fiber + (1|HM), 
                                        subset(metadata.oars.stool.asv.double, response == "high")) %>% summary() %>% coef()
stats.load.mgx.ph.interact = lmerTest::lmer(scale(load.mgx) ~ reltiming*response + phase +diagnosis+ adj.fiber + (1|HM), metadata.oars.stool.asv.double) %>% summary() %>% coef()
# make data.frame
stats.load.mgx.ph = data.frame(Response = c("Strong Response", "Weak Response"),
                               pval = c(stats.load.mgx.ph.high[2,5], stats.load.mgx.ph.low[2,5]))
metadata.oars.stool.asv.double$RS_Name = factor(metadata.oars.stool.asv.double$RS_Name, levels=rs.names)

metadata.oars.load.mgx.double.plot = ggplot(subset(metadata.oars.stool.asv.double, !hm_phase %in% c("HM0924_rs2", "HM0932_rs2", "HM0759_rs2")) %>%
                                              mutate(reltiming = factor(ifelse(reltiming == "pre", "Pre", "Post"), levels=c("Pre", "Post")))%>%
                                              mutate(Response = ifelse(response == "low", "Weak Response", "Strong Response")),
                                            aes(x=reltiming, y=load.mgx))+
  geom_boxplot(width=0.5)+
  geom_line(aes(group=hm_phase), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(shape=21, aes(fill=RS_Name), size=2)+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  geom_text(data = stats.load.mgx.ph, x=1.5, y=(max(na.omit(metadata.oars.stool.asv.double$load.mgx))), 
            aes(label=paste("p:", round(pval, digits=3))), size=3.5)+
  facet_wrap(~Response)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="",
       y="Microbial Load (Log10 Predicted)")
metadata.oars.load.mgx.double.plot
stats.load.mgx.ph.interact # no interaction

# :: ASV Richness ---------------------------------------------------------

# delta calculated with richness

stats.richness.ph.low = lmerTest::lmer(scale(log10(richness)) ~ reltiming + phase + diagnosis+adj.fiber + (1|HM), 
                                       subset(metadata.oars.stool.asv.double, response == "low")) %>% summary() %>% coef()
stats.richness.ph.high = lmerTest::lmer(scale(log10(richness)) ~ reltiming + phase +diagnosis+adj.fiber +  (1|HM), 
                                        subset(metadata.oars.stool.asv.double, response == "high")) %>% summary() %>% coef()
stats.richness.ph.interact = lmerTest::lmer(scale(log10(richness)) ~ reltiming*response + phase + diagnosis+adj.fiber + (1|HM), metadata.oars.stool.asv.double) %>% summary() %>% coef()

# make data.frame
stats.richness.ph = data.frame(Response = c("Strong Response", "Weak Response"),
                               pval = c(stats.richness.ph.high[2,5], stats.richness.ph.low[2,5]))

metadata.oars.stool.asv.double$RS_Name = factor(metadata.oars.stool.asv.double$RS_Name, levels=rs.names)

metadata.oars.stool.asv.double.richness.plot = ggplot(subset(metadata.oars.stool.asv.double, !hm_phase %in% c("HM0924_rs2", "HM0932_rs2", "HM0759_rs2")) %>% 
                                                        mutate(HM_phase = paste(HM, phase, sep="_")) %>%
                                                        mutate(reltiming = factor(ifelse(reltiming == "pre", "Pre", "Post"), levels=c("Pre", "Post"))) %>%
                                                        mutate(Response = ifelse(response == "low", "Weak Response", "Strong Response")),
                                                            aes(x=reltiming, y=richness))+
  geom_boxplot(width=0.5)+
  geom_line(aes(group=HM_phase), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(shape=21, aes(fill=RS_Name), size=2)+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  geom_text(data = stats.richness.ph, x=1.5, y=(max(na.omit(metadata.oars.stool.asv.double$richness))), 
            aes(label=paste("p:", round(pval, digits=3))), size=3.5)+
  facet_wrap(~Response)+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="",
       y="ASV Richness")
metadata.oars.stool.asv.double.richness.plot
# funneling effect!
stats.richness.ph.interact # sig interaction


# :: ASV Shannon ----------------------------------------------------------

# delta calculated with shannon

stats.shannon.ph.low = lmerTest::lmer(scale(shannon) ~ reltiming + phase + diagnosis+adj.fiber + (1|HM), 
                                       subset(metadata.oars.stool.asv.double, response == "low")) %>% summary() %>% coef()
stats.shannon.ph.high = lmerTest::lmer(scale(shannon) ~ reltiming + phase + diagnosis+adj.fiber + (1|HM), 
                                        subset(metadata.oars.stool.asv.double, response == "high")) %>% summary() %>% coef()
stats.shannon.ph.interact = lmerTest::lmer(scale(shannon) ~ reltiming*response + phase + diagnosis+adj.fiber + (1|HM), metadata.oars.stool.asv.double) %>% summary() %>% coef()

# make data.frame
stats.shannon.ph = data.frame(Response = c("Strong Response", "Weak Response"),
                               pval = c(stats.shannon.ph.high[2,5], stats.shannon.ph.low[2,5]))

metadata.oars.stool.asv.double$RS_Name = factor(metadata.oars.stool.asv.double$RS_Name, levels=rs.names)

metadata.oars.stool.asv.double.shannon.plot = ggplot(subset(metadata.oars.stool.asv.double, !hm_phase %in% c("HM0924_rs2", "HM0932_rs2", "HM0759_rs2")) %>%
                                                       mutate(HM_phase = paste(HM, phase, sep="_")) %>%
                                                       mutate(reltiming = factor(ifelse(reltiming == "pre", "Pre", "Post"), levels=c("Pre", "Post"))) %>%
                                                       mutate(Response = ifelse(response == "low", "Weak Response", "Strong Response")),
                                                      aes(x=reltiming, y=shannon))+
  geom_boxplot(width=0.5)+
  geom_line(aes(group=HM_phase), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(shape=21, aes(fill=RS_Name), size=2)+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  geom_text(data = stats.shannon.ph, x=1.5, y=(max(na.omit(metadata.oars.stool.asv.double$shannon))), 
            aes(label=paste("p:", round(pval, digits=3))), size=3.5)+
  facet_wrap(~Response)+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="",
       y="Shannon Diversity")
metadata.oars.stool.asv.double.shannon.plot
stats.shannon.ph.interact # no interaction

# :: ASV Functional Redundancy --------------------------------------------

stats.fd.ph.low = lmerTest::lmer(scale(fd) ~ reltiming + phase + diagnosis+adj.fiber + (1|HM), 
                                      subset(metadata.oars.stool.asv.double, response == "low")) %>% summary() %>% coef()
stats.fd.ph.high = lmerTest::lmer(scale(fd) ~ reltiming + phase + diagnosis+adj.fiber + (1|HM), 
                                       subset(metadata.oars.stool.asv.double, response == "high")) %>% summary() %>% coef()
stats.fd.ph.interact = lmerTest::lmer(scale(fd) ~ reltiming*response + phase + diagnosis+adj.fiber + (1|HM), metadata.oars.stool.asv.double) %>% summary() %>% coef()

# make data.frame
stats.fd.ph = data.frame(Response = c("Strong Response", "Weak Response"),
                              pval = c(stats.fd.ph.high[2,5], stats.fd.ph.low[2,5]))

metadata.oars.stool.asv.double$RS_Name = factor(metadata.oars.stool.asv.double$RS_Name, levels=rs.names)

metadata.oars.stool.asv.double.fd.plot = ggplot(subset(metadata.oars.stool.asv.double, !hm_phase %in% c("HM0924_rs2", "HM0932_rs2", "HM0759_rs2")) %>%
                                                  mutate(HM_phase = paste(HM, phase, sep="_")) %>%
                                                  mutate(reltiming = factor(ifelse(reltiming == "pre", "Pre", "Post"), levels=c("Pre", "Post"))) %>%
                                                  mutate(Response = ifelse(response == "low", "Weak Response", "Strong Response")),
                                                     aes(x=reltiming, y=fd))+
  geom_boxplot(width=0.5)+
  geom_line(aes(group=HM_phase), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(shape=21, aes(fill=RS_Name), size=2)+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  geom_text(data = stats.fd.ph, x=1.5, y=(max(na.omit(metadata.oars.stool.asv.double$fd))), 
            aes(label=paste("p =", round(pval, digits=3))), size=3.5)+
  facet_wrap(~Response)+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="",
       y="Functional Redundancy")
metadata.oars.stool.asv.double.fd.plot
stats.fd.ph.interact # sig interaction


# :: ASV Butyrogens I -------------------------------------------------------

# delta calculated with butyrogens

stats.but.i.ph.low = lmerTest::lmer(scale(log10(but.i)) ~ reltiming + phase + diagnosis+adj.fiber + (1|HM), 
                                  subset(metadata.oars.stool.asv.double, response == "low")) %>% summary() %>% coef()
stats.but.i.ph.high = lmerTest::lmer(scale(log10(but.i)) ~ reltiming + phase +diagnosis+ adj.fiber + (1|HM), 
                                   subset(metadata.oars.stool.asv.double, response == "high")) %>% summary() %>% coef()
stats.but.i.ph.interact = lmerTest::lmer(scale(log10(but.i)) ~ reltiming*response + phase + diagnosis+adj.fiber + (1|HM), metadata.oars.stool.asv.double) %>% summary() %>% coef()

# make data.frame
stats.but.i.ph = data.frame(Response = c("Strong Response", "Weak Response"),
                          pval = c(stats.but.i.ph.high[2,5], stats.but.i.ph.low[2,5]))

metadata.oars.stool.asv.double$RS_Name = factor(metadata.oars.stool.asv.double$RS_Name, levels=rs.names)

metadata.oars.stool.asv.double.but.i.plot = ggplot(subset(metadata.oars.stool.asv.double, !hm_phase %in% c("HM0924_rs2", "HM0932_rs2", "HM0759_rs2")) %>% 
                                                     mutate(HM_phase = paste(HM, phase, sep="_")) %>%
                                                     mutate(reltiming = factor(ifelse(reltiming == "pre", "Pre", "Post"), levels=c("Pre", "Post")))%>%
                                                     mutate(Response = ifelse(response == "low", "Weak Response", "Strong Response")),
                                                 aes(x=reltiming, y=but.i*100))+
  geom_boxplot(width=0.5)+
  geom_line(aes(group=HM_phase), linetype=2, alpha=0.5, linewidth=0.5)+
  geom_point(shape=21, aes(fill=RS_Name), size=3)+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  geom_text(data = stats.but.i.ph, x=1.5, y=(max(na.omit(metadata.oars.stool.asv.double$but.i)*100)), 
            aes(label=paste("p =", round(pval, digits=3))), size=4.5)+
  facet_wrap(~Response)+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="",
       y="Butyrogens (%)")
metadata.oars.stool.asv.double.but.i.plot
# not significant


# :: ASV Butyrogens II (Vital) -------------------------------------------------------


# delta calculated with butyrogens

stats.but.ii.ph.low = lmerTest::lmer(scale(log10(but.ii)) ~ reltiming + phase + diagnosis+adj.fiber + (1|HM), 
                                       subset(metadata.oars.stool.asv.double, response == "low")) %>% summary() %>% coef()
stats.but.ii.ph.high = lmerTest::lmer(scale(log10(but.ii)) ~ reltiming + phase +diagnosis+ adj.fiber + (1|HM), 
                                        subset(metadata.oars.stool.asv.double, response == "high")) %>% summary() %>% coef()
stats.but.ii.ph.interact = lmerTest::lmer(scale(log10(but.ii)) ~ reltiming*response + phase + diagnosis+adj.fiber + (1|HM), metadata.oars.stool.asv.double) %>% summary() %>% coef()

# make data.frame
stats.but.ii.ph = data.frame(Response = c("Strong Response", "Weak Response"),
                               pval = c(stats.but.ii.ph.high[2,5], stats.but.ii.ph.low[2,5]))

metadata.oars.stool.asv.double$RS_Name = factor(metadata.oars.stool.asv.double$RS_Name, levels=rs.names)

metadata.oars.stool.asv.double.but.ii.plot = ggplot(subset(metadata.oars.stool.asv.double, !hm_phase %in% c("HM0924_rs2", "HM0932_rs2", "HM0759_rs2")) %>% 
                                                      mutate(HM_phase = paste(HM, phase, sep="_")) %>%
                                                      mutate(reltiming = factor(ifelse(reltiming == "pre", "Pre", "Post"), levels=c("Pre", "Post")))%>%
                                                      mutate(Response = ifelse(response == "low", "Weak Response", "Strong Response")),
                                                      aes(x=reltiming, y=but.ii*100))+
  geom_boxplot(width=0.5)+
  geom_line(aes(group=HM_phase), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(shape=21, aes(fill=RS_Name), size=2)+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  geom_text(data = stats.but.ii.ph, x=1.5, y=(max(na.omit(metadata.oars.stool.asv.double$but.ii)*100)), 
            aes(label=paste("p =", round(pval, digits=3))), size=3.5)+
  facet_wrap(~Response)+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="",
       y="ASV Butyrogens %")
metadata.oars.stool.asv.double.but.ii.plot
# not significant



# :: ASV Bray-Curtis ------------------------------------------------------

# Calculate distances (use 10% prevalence filtered with non-compliant removed)
dim(oars.asv.data.median.1.filt)
oars.asv.bray = vegan::vegdist(oars.asv.data.median.1.filt, method="bray") 

# Subset to samples and comparisons of interest
oars.asv.bray.df = reshape2::melt(as.matrix(oars.asv.bray))
# configure as a "mapping file" so we can use other functions
oars.asv.bray.df$standard.name = oars.asv.bray.df$Var1
# subset to self-comparisons
oars.asv.bray.df = subset(oars.asv.bray.df, substr(Var1, 1, 6) == substr(Var2, 1, 6))
# delete same-sample comparisons
oars.asv.bray.df = subset(oars.asv.bray.df, Var1 != Var2)
# replace Var with timing
oars.asv.bray.df$Var1 = metadata.oars.stool.asv$timing[match(oars.asv.bray.df$Var1, metadata.oars.stool.asv$standard.name)]
oars.asv.bray.df$Var2 = metadata.oars.stool.asv$timing[match(oars.asv.bray.df$Var2, metadata.oars.stool.asv$standard.name)]
# calculate distances to preceding / baseline
oars.asv.bray.df = oars.asv.bray.df %>%
  mutate(vars = paste(Var1, Var2, sep="_")) %>%
  subset(vars %in% c("0M_3M", "3M_6M", "0M_6M", "6M_9M", "9M_12M", "6M_12M"))
# add other metadata
oars.asv.bray.df = merge(oars.asv.bray.df, metadata.oars.stool.asv, by="standard.name")
lmerTest::lmer(value ~ vars + timing  + diagnosis+adj.fiber + (1|HM), subset(oars.asv.bray.df, vars %in% c("0M_3M", "3M_6M", "6M_9M", "9M_12M"))) %>% summary()
# none are significant at group level
ggplot(subset(oars.asv.bray.df, vars %in% c("0M_3M", "3M_6M", "6M_9M", "9M_12M")),
       aes(x=vars, y=value))+
  geom_boxplot(outlier.shape=NA, width=0.3)+
  geom_path(aes(group=HM), alpha=0.5, linetype=2, linewidth=0.2)+
  geom_point(shape=21, aes(fill=(timing)), size=2)+
  scale_fill_manual(values=c(2,2, "grey", "grey"))+
  theme_classic()+theme(legend.position="none")+
  labs(x="", y="Distance ")
unique(oars.asv.bray.df$HM)

# OPTIONAL: and responses
oars.asv.bray.responses.df = merge(oars.asv.bray.df, subset(metadata.oars.stool.double, reltiming == "pre")[,c("standard.name", "response")]) %>%
  subset(vars %in% c("0M_3M", "3M_6M"))
# stats
lmerTest::lmer(value ~ response + timing  + diagnosis+adj.fiber + (1|HM), oars.asv.bray.responses.df) %>% summary()
ggplot(oars.asv.bray.responses.df,
       aes(x=response, y=value))+
  geom_boxplot(outlier.shape=NA, width=0.3)+
  ggbeeswarm::geom_beeswarm(shape=21, aes(fill=response), size=2)+
  ggnetwork::geom_nodetext_repel(aes(label=paste(HM, vars, sep=" ")),size=3)+
  theme_classic()+theme(legend.position="none")+
  labs(x="", y="Distance to Preceding Stool")
# not significant


# :: ASV Beta Diversty -------------------------------------------------------------
# remake data
beta.trajectory.data = beta.trajectory(oars.asv.bray)

beta.trajectory.data =   merge(metadata.oars.stool.asv, 
                               beta.trajectory.data, by="standard.name")

beta.trajectory.data = merge(beta.trajectory.data, oars.rapidaim.scores[,c("HM", "timing","rs.selected", "delta.ph")], by=c("HM", "timing"))
# remove extra columns
beta.trajectory.data = beta.trajectory.data[,!colnames(beta.trajectory.data) %in% c("fcal", "stool_water_perc", "richness", "shannon", "fd", "but.i", "but.ii", "load.asv", "load.mgx")]

beta.trajectory.data.double = doubleup(beta.trajectory.data, delta=T)
beta.trajectory.data.double = subset(beta.trajectory.data.double, delta.between != 0)
beta.trajectory.data.double$response = ifelse(beta.trajectory.data.double$delta.ph > mean(metadata.oars.stool.double$delta.ph), "low", "high")


# INTER
stats.beta.between.ph.low = lmerTest::lmer(scale((between.beta)) ~ reltiming + phase + diagnosis+adj.fiber + (1|HM), 
                                       subset(beta.trajectory.data.double, response == "low")) %>% summary() %>% coef()
stats.beta.between.ph.high = lmerTest::lmer(scale(between.beta) ~ reltiming + phase +diagnosis+adj.fiber +  (1|HM), 
                                        subset(beta.trajectory.data.double, response == "high")) %>% summary() %>% coef()
stats.beta.between.ph.interact = lmerTest::lmer(scale(between.beta) ~ reltiming*response + phase + diagnosis+adj.fiber + (1|HM), beta.trajectory.data.double) %>% summary() %>% coef()

lmerTest::lmer(scale(between.beta) ~ response + phase + diagnosis+adj.fiber + (1|HM), 
                                                subset(beta.trajectory.data.double, reltiming == "pre")) %>% summary() %>% coef()
# note: Weak response has greater inter-individual variability than Strong

# make data.frame
stats.beta.between.ph = data.frame(Response = c("Strong Response", "Weak Response"),
                               pval = c(stats.beta.between.ph.high[2,5], stats.beta.between.ph.low[2,5]))

beta.trajectory.data.double$RS_Name = factor(beta.trajectory.data.double$RS_Name, levels=rs.names)

metadata.oars.stool.asv.double.beta.between.plot = ggplot(beta.trajectory.data.double %>% 
                                                        mutate(HM_phase = paste(HM, phase, sep="_")) %>%
                                                        subset( !HM_phase %in% c("HM0924_rs2", "HM0932_rs2", "HM0759_rs2")) %>%
                                                        mutate(reltiming = factor(ifelse(reltiming == "pre", "Pre", "Post"), levels=c("Pre", "Post"))) %>%
                                                        mutate(Response = ifelse(response == "low", "Weak Response", "Strong Response")),
                                                      aes(x=reltiming, y=between.beta))+
  geom_boxplot(width=0.5)+
  geom_line(aes(group=HM_phase), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(shape=21, aes(fill=RS_Name), size=2)+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  geom_text(data = stats.beta.between.ph, x=1.5, y=(max(na.omit(beta.trajectory.data.double$between.beta))), 
            aes(label=paste("p:", round(pval, digits=3))), size=3.5)+
  facet_wrap(~Response)+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="",
       y="Mean Bray-Curtis Dissimilarity")
metadata.oars.stool.asv.double.beta.between.plot

stats.beta.between.ph.interact # nearly p < 0.05

# INTRA
stats.beta.within.ph.low = lmerTest::lmer(scale(within.beta) ~ reltiming + phase + diagnosis+adj.fiber + (1|HM), 
                                           subset(beta.trajectory.data.double, response == "low")) %>% summary() %>% coef()
stats.beta.within.ph.high = lmerTest::lmer(scale(within.beta) ~ reltiming + phase +diagnosis+adj.fiber +  (1|HM), 
                                            subset(beta.trajectory.data.double, response == "high")) %>% summary() %>% coef()
stats.beta.within.ph.interact = lmerTest::lmer(scale(within.beta) ~ reltiming*response + phase + diagnosis+adj.fiber + (1|HM), beta.trajectory.data.double) %>% summary() %>% coef()

lmerTest::lmer(scale(within.beta) ~ response + phase + diagnosis+adj.fiber + (1|HM), 
               subset(beta.trajectory.data.double, reltiming == "pre")) %>% summary() %>% coef()
# note: No difference at baseline

# make data.frame
stats.beta.within.ph = data.frame(Response = c("Strong Response", "Weak Response"),
                                   pval = c(stats.beta.within.ph.high[2,5], stats.beta.within.ph.low[2,5]))

beta.trajectory.data.double$RS_Name = factor(beta.trajectory.data.double$RS_Name, levels=rs.names)

metadata.oars.stool.asv.double.beta.within.plot = ggplot(beta.trajectory.data.double %>% 
                                                            mutate(HM_phase = paste(HM, phase, sep="_")) %>%
                                                            subset( !HM_phase %in% c("HM0924_rs2", "HM0932_rs2", "HM0759_rs2")) %>%
                                                            mutate(reltiming = factor(ifelse(reltiming == "pre", "Pre", "Post"), levels=c("Pre", "Post"))) %>%
                                                            mutate(Response = ifelse(response == "low", "Weak Response", "Strong Response")),
                                                          aes(x=reltiming, y=within.beta))+
  geom_boxplot(width=0.5)+
  geom_line(aes(group=HM_phase), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(shape=21, aes(fill=RS_Name), size=2)+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  geom_text(data = stats.beta.within.ph, x=1.5, y=(max(na.omit(beta.trajectory.data.double$within.beta))), 
            aes(label=paste("p:", round(pval, digits=3))), size=3.5)+
  facet_wrap(~Response)+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="",
       y="Inter-individual Dissimilarity")
metadata.oars.stool.asv.double.beta.within.plot

stats.beta.within.ph.interact # no interaction



# [] Prep data for lmer ----------------------------------------------------

# used for ML, too

oars.asv.data.glom.prep = delta.omic.prepare(oars.asv.data.glom,  # rarefied reads (50,000)
                                             normalize=T, # normalize to 100%
                                             min.abun = 1, # 1 read
                                             prev=0.1) # 10% prev per group
dim(oars.asv.data.glom.prep) # 495

oars.mgx.prep = delta.omic.prepare(oars.mgx.taxa, # sum to ~100% data (~ because of median)
                                   normalize=T, # re-normalize to 100%
                                   min.abun = 0.01, # 0.01% min abun for prev
                                   prev = 0.1) # 10% prev per group
dim(oars.mgx.prep) # 627

oars.mpx.kegg.prep = delta.omic.prepare(oars.mpx.kegg.mat, # norm intensity
                                        normalize=F, # don't normalize
                                        min.abun = 1, # detected at all
                                        prev=0.1) # 10% prev per group
dim(oars.mpx.kegg.prep) # 184

oars.mpx.cog.prep = delta.omic.prepare(oars.mpx.cog.mat, 
                                       normalize=F, # don't normalize
                                       min.abun = 1,# detected
                                       prev=0.1)# 10% prev per group
dim(oars.mpx.cog.prep) # 2530

oars.mpx.cazy.prep = delta.omic.prepare(oars.mpx.cazy.mat, 
                                        normalize=F, # don't normalize
                                        min.abun = 1,# detected
                                        prev=0.1)# 10% prev per group
dim(oars.mpx.cazy.prep) # 74

oars.mbx.prep = delta.omic.prepare(oars.mbx.annotated.mat, 
                                        normalize=F, # don't normalize
                                        min.abun = 1, # detected
                                        pseudo=T,
                                        prev=0.8) # 80% prev per group
dim(oars.mbx.prep) # 197



# :: ASV lmer ---------------------------------------------------------

# 2025_05_27 NOTES ON LMER:
# Interaction-based model preserves power
# while no features are sig after FDR,
# those nominally sig seem to agree between 16S and MGX 
# (e.g. F. prausnitzii, Bifidobacterium, Bacteroides, Blautia, Clostridium)


## LMER
oars.asv.data.glom.lmer = delta.omic.lmer(
  split = FALSE,
  data = oars.asv.data.glom.prep,
  mpx = F
  )

oars.asv.data.glom.lmer %>% arrange(pval)

# volcano plot
oars.asv.data.glom.lmer.volcano = ggplot(oars.asv.data.glom.lmer,
                                         aes(x=estimate, y=padj))+
  geom_point(shape=21, aes(fill=estimate), size=2)+
  geom_hline(yintercept=(0.20), linetype=2, alpha=0.5)+
  geom_vline(xintercept=c(-1,1), linetype=2, alpha=0.5)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.05,0.1, 0.20, 0.5, 1))+
  #ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, taxa, NA)), size=3)+
  scale_fill_gradient2(low="blue", high="red")+
  facet_wrap(~"Interaction: ASV")+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10))+
  labs(x="Interaction Coefficient", y="FDR")
oars.asv.data.glom.lmer.volcano
# nothing sig

## Repeat with split // defunct
oars.asv.data.glom.lmer.split = delta.omic.lmer(
  split = TRUE,
  data = oars.asv.data.glom.prep
)

oars.asv.data.glom.lmer.split.volcano = ggplot(oars.asv.data.glom.lmer.split %>%
                                                 mutate(Response = ifelse(response == "high", "Strong", "Weak")),
                                         aes(x=estimate, y=padj))+
  geom_point(shape=21, aes(fill=estimate), size=2)+
  geom_hline(yintercept=0.20, linetype=2, alpha=0.5)+
  geom_vline(xintercept=c(-1,1), linetype=2, alpha=0.5)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.05,0.1, 0.20, 0.5, 1))+
  ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, taxa, NA)), size=3)+
  scale_fill_gradient2(low="blue", high="red")+
  facet_wrap(~Response)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="Adjusted Coefficient", y="FDR", title="ASV")
oars.asv.data.glom.lmer.split.volcano
# nothing sig

# :: MGX lmer -------------------------------------------------------------

## LMER
oars.mgx.lmer = delta.omic.lmer(
  data = oars.mgx.prep,
  split=F,
  mpx=F
)
oars.mgx.lmer %>% arrange(pval)

# volcano plot (FDR)
oars.mgx.lmer.volcano = ggplot(oars.mgx.lmer,
                               aes(x=estimate, y=padj))+
  geom_point(shape=21, aes(fill=estimate), size=2)+
  geom_hline(yintercept=(0.20), linetype=2, alpha=0.5)+
  geom_vline(xintercept=c(-1,1), linetype=2, alpha=0.5)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.05,0.1, 0.20, 0.5, 1))+
  #ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, taxa, NA)), size=3)+
  scale_fill_gradient2(low="blue", high="red")+
  facet_wrap(~"Interaction: Species")+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10))+
  labs(x="Interaction Coefficient", y="FDR")
oars.mgx.lmer.volcano
# nothing sig

# volcano plot (unadjusted)
oars.mgx.lmer.volcano.unadj = ggplot(oars.mgx.lmer,
                               aes(x=estimate, y=pval))+
  geom_point(shape=21, aes(fill=estimate), size=2)+
  geom_hline(yintercept=(0.05), linetype=2, alpha=0.5)+
  geom_vline(xintercept=c(-1,1), linetype=2, alpha=0.5)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05, 1))+
  ggrepel::geom_text_repel(aes(label=ifelse(pval < 0.05, taxa, NA)), size=3)+
  scale_fill_gradient2(low="blue", high="red")+
  facet_wrap(~"Species")+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10))+
  labs(x="Interaction Coefficient", y="p-value")

## Repeat with split // defunct
oars.mgx.lmer.split = delta.omic.lmer(
  split = TRUE,
  data = oars.mgx.prep
)
oars.mgx.lmer.split.volcano = ggplot(oars.mgx.lmer.split %>%
                                     mutate(Response = ifelse(response == "high", "Strong", "Weak")),
                                     aes(x=estimate, y=padj))+
  geom_point(shape=21, aes(fill=estimate), size=2)+
  geom_hline(yintercept=0.20, linetype=2, alpha=0.5)+
  geom_vline(xintercept=c(-1,1), linetype=2, alpha=0.5)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.05,0.1, 0.20, 0.5, 1))+
  ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, taxa, NA)), size=3)+
  scale_fill_gradient2(low="blue", high="red")+
  facet_wrap(~Response)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="Adjusted Coefficient", y="FDR", title="Species")
oars.mgx.lmer.split.volcano
# nothing sig

oars.mgx.lmer.split.volcano.unadj = ggplot(oars.mgx.lmer.split %>%
                                                  mutate(Response = ifelse(response == "high", "Strong", "Weak")),
                                                aes(x=estimate, y=padj))+
  geom_point(shape=21, aes(fill=estimate), size=2)+
  geom_hline(yintercept=0.20, linetype=2, alpha=0.5)+
  geom_vline(xintercept=c(-1,1), linetype=2, alpha=0.5)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05,1))+
  ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20 & estimate>4, taxa, NA)), size=3)+
  scale_fill_gradient2(low="blue", high="red")+
  facet_wrap(~Response)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="Adjusted Coefficient", y="FDR", title="Species")
oars.mgx.lmer.split.volcano.unadj


# histogram of abs(lfc)
lm(abs(estimate) ~ response, oars.mgx.lmer.split) %>% summary()
ggplot(oars.mgx.lmer.split)+
  geom_violin(aes(x=abs(estimate), y=response, fill=response),
              draw_quantiles=0.5)+
  coord_flip()+
  theme_minimal()+theme(legend.position="none")+
  labs(y="", x="|Standardized Coefficient|")
# stronger impacts in strong responders

# :: MPX KEGG lmer -------------------------------------------------------------

dim(oars.mpx.kegg.mat)

# create map for protein names
oars.mpx.kegg.names = data.frame(good = colnames(oars.mpx.kegg.mat),
                                    feature = make.names(colnames(oars.mpx.kegg.mat))) %>%
  subset(feature %in% colnames(oars.mpx.kegg.prep))

## LMER
oars.mpx.kegg.lmer = delta.omic.lmer(
  split = FALSE,
  data = oars.mpx.kegg.prep,
  mpx = T
)
oars.mpx.kegg.lmer %>% arrange(pval)
oars.mpx.kegg.lmer$feature = oars.mpx.kegg.names$good[match(oars.mpx.kegg.lmer$taxa, oars.mpx.kegg.names$feature)]


# volcano plot
oars.mpx.kegg.lmer.volcano = ggplot(oars.mpx.kegg.lmer,
                                    aes(x=estimate, y=padj))+
  geom_point(shape=21, aes(fill=estimate), size=2)+
  geom_hline(yintercept=(0.20), linetype=2, alpha=0.5)+
  geom_vline(xintercept=c(-1,1), linetype=2, alpha=0.5)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.05,0.1, 0.20, 0.5, 1))+
  #ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.05, gsub("\\.", " ", gsub("\\.\\.", ", ", taxa)), NA)), size=3)+
  scale_fill_gradient2(low="blue", high="red")+
  facet_wrap(~"Interaction: Pathway")+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10))+
  labs(x="Interaction Coefficient", y="FDR")
oars.mpx.kegg.lmer.volcano


## Repeat with split // defunct
oars.mpx.kegg.lmer.split = delta.omic.lmer(
  split = TRUE,
  data = oars.mpx.kegg.prep
)
oars.mpx.kegg.lmer.split$feature = oars.mpx.kegg.names$good[match(oars.mpx.kegg.lmer.split$taxa, oars.mpx.kegg.names$feature)]


oars.mpx.kegg.lmer.split.volcano = ggplot(oars.mpx.kegg.lmer.split  %>%
                                            mutate(Response = ifelse(response == "high", "Strong", "Weak")),
                                          aes(x=estimate, y=padj))+
  geom_point(shape=21, aes(fill=estimate), size=2)+
  geom_hline(yintercept=0.20, linetype=2, alpha=0.5)+
  geom_vline(xintercept=c(-1,1), linetype=2, alpha=0.5)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.05,0.1, 0.20, 0.5, 1))+
  ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, gsub("\\.", " ", gsub("\\.\\.", ", ", taxa)), NA)), size=3)+
  scale_fill_gradient2(low="blue", high="red")+
  facet_wrap(~Response)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="Adjusted Coefficient", y="FDR", title="Pathway")
oars.mpx.kegg.lmer.split.volcano
# nothing sig

# :: MPX COG lmer -------------------------------------------------------------

dim(oars.mpx.cog.mat)

# create map for protein names
oars.mpx.cog.protein.names = data.frame(good = colnames(oars.mpx.cog.mat),
                                    feature = make.names(colnames(oars.mpx.cog.mat))) %>%
  subset(feature %in% colnames(oars.mpx.cog.prep))

## LMER
oars.mpx.cog.lmer = delta.omic.lmer(
  split = FALSE,
  data = oars.mpx.cog.prep,
  mpx = T
)
oars.mpx.cog.lmer %>% arrange(padj)
oars.mpx.cog.lmer$feature = oars.mpx.cog.protein.names$good[match(oars.mpx.cog.lmer$taxa, oars.mpx.cog.protein.names$feature)]


# volcano plot
oars.mpx.cog.lmer.volcano = ggplot(oars.mpx.cog.lmer%>%
                                     mutate(taxa = gsub("\\."," ", taxa)),
                                   aes(x=estimate, y=padj))+
  geom_point(shape=21, aes(fill=estimate), size=2)+
  geom_hline(yintercept=(0.20), linetype=2, alpha=0.5)+
  geom_vline(xintercept=c(-1,1), linetype=2, alpha=0.5)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.05,0.1, 0.20, 0.5, 1))+
  #ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.2, taxa, NA)), size=3)+
  scale_fill_gradient2(low="blue", high="red")+
  facet_wrap(~"Interaction: COG")+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10))+
  labs(x="Interaction Coefficient", y="FDR")
oars.mpx.cog.lmer.volcano
# a few sig

## Repeat with split // defunct
oars.mpx.cog.lmer.split = delta.omic.lmer(
  split = TRUE,
  data = oars.mpx.cog.prep
)
oars.mpx.cog.lmer.split$feature = oars.mpx.cog.protein.names$good[match(oars.mpx.cog.lmer.split$taxa, oars.mpx.cog.protein.names$feature)]


oars.mpx.cog.lmer.split.volcano = ggplot(oars.mpx.cog.lmer.split %>%
                                           mutate(taxa = gsub("  ", " ", gsub("\\.", " ", taxa))) %>%
                                           mutate(Response = ifelse(response == "high", "Strong", "Weak")),
                                         aes(x=estimate, y=padj))+
  geom_point(shape=21, aes(fill=estimate), size=2)+
  geom_hline(yintercept=0.20, linetype=2, alpha=0.5)+
  geom_vline(xintercept=c(-1,1), linetype=2, alpha=0.5)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.05,0.1, 0.20, 0.5, 1))+
  ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, taxa, NA)), size=3)+
  scale_fill_gradient2(low="blue", high="red")+
  facet_wrap(~Response)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="Adjusted Coefficient", y="FDR", title="COG")
oars.mpx.cog.lmer.split.volcano
# 1 sig


# :: MPX CAZy lmer -------------------------------------------------------------

dim(oars.mpx.cazy.mat)

# create map for cazy names
oars.mpx.cazy.names = data.frame(good = colnames(oars.mpx.cazy.mat),
                                    feature = make.names(colnames(oars.mpx.cazy.mat))) %>%
  subset(feature %in% colnames(oars.mpx.cazy.prep))

## LMER
oars.mpx.cazy.lmer = delta.omic.lmer(
  split = FALSE,
  data = oars.mpx.cazy.prep,
  mpx=T
)
oars.mpx.cazy.lmer %>% arrange(pval)
oars.mpx.cazy.lmer$feature = oars.mpx.cazy.names$good[match(oars.mpx.cazy.lmer$taxa, oars.mpx.cazy.names$feature)]

# volcano plot
oars.mpx.cazy.lmer.volcano = ggplot(oars.mpx.cazy.lmer,
                                    aes(x=estimate, y=padj))+
  geom_point(shape=21, aes(fill=estimate), size=2)+
  geom_hline(yintercept=(0.20), linetype=2, alpha=0.5)+
  geom_vline(xintercept=c(-1,1), linetype=2, alpha=0.5)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.05,0.1, 0.20, 0.5, 1))+
  #ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.2, taxa, NA)), size=3)+
  scale_fill_gradient2(low="blue", high="red")+
  facet_wrap(~"Interaction: CAZy")+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10))+
  labs(x="Interaction Coefficient", y="FDR")
oars.mpx.cazy.lmer.volcano
# none sig

oars.mpx.cazy.lmer.volcano.unadj = ggplot(oars.mpx.cazy.lmer,
                                     aes(x=estimate, y=pval))+
  geom_point(shape=21, aes(fill=estimate), size=2)+
  geom_hline(yintercept=(0.05), linetype=2, alpha=0.5)+
  geom_vline(xintercept=c(-1,1), linetype=2, alpha=0.5)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05, 1))+
  ggrepel::geom_text_repel(aes(label=ifelse(pval < 0.05, taxa, NA)), size=3)+
  scale_fill_gradient2(low="blue", high="red")+
  facet_wrap(~"CAZymes")+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10))+
  labs(x="Interaction Coefficient", y="p-value")

## Repeat with split
oars.mpx.cazy.lmer.split = delta.omic.lmer(
  split = TRUE,
  data = oars.mpx.cazy.prep
)
oars.mpx.cazy.lmer.split$feature = oars.mpx.cazy.names$good[match(oars.mpx.cazy.lmer.split$taxa, oars.mpx.cazy.names$feature)]

oars.mpx.cazy.lmer.split.volcano = ggplot(oars.mpx.cazy.lmer.split %>%
                                            mutate(Response = ifelse(response == "high", "Strong", "Weak")),
                                          aes(x=estimate, y=padj))+
  geom_point(shape=21, aes(fill=estimate), size=2)+
  geom_hline(yintercept=0.20, linetype=2, alpha=0.5)+
  geom_vline(xintercept=c(-1,1), linetype=2, alpha=0.5)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.05,0.1, 0.20, 0.5, 1))+
  ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, taxa, NA)), size=3)+
  scale_fill_gradient2(low="blue", high="red")+
  facet_wrap(~Response)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="Adjusted Coefficient", y="FDR", title="CAZy")
oars.mpx.cazy.lmer.split.volcano
# 1 sig

oars.mpx.cazy.lmer.split.volcano.unadj = ggplot(oars.mpx.cazy.lmer.split %>%
                                            mutate(Response = ifelse(response == "high", "Intermediate", "Minor")),
                                          aes(x=estimate, y=pval))+
  geom_point(shape=21, aes(fill=estimate), size=2)+
  geom_hline(yintercept=0.05, linetype=2, alpha=0.5)+
  geom_vline(xintercept=c(-1,1), linetype=2, alpha=0.5)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05,1))+
  ggrepel::geom_text_repel(aes(label=ifelse(pval < 0.01 & estimate>0.7, taxa, NA)), size=3)+
  scale_fill_gradient2(low="blue", high="red")+
  facet_wrap(~Response)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="Adjusted Coefficient", y="p-value", title="CAZy")
oars.mpx.cazy.lmer.split.volcano.unadj


# plot all
oars.mpx.kegg.lmer.split.volcano+
  oars.mpx.cog.lmer.split.volcano+
  oars.mpx.cazy.lmer.split.volcano


# :: MPX CAZy Starch ------------------------------------------------


# delta calculated with starch

stats.starch.ph.low = lmerTest::lmer(scale(log10(starch)) ~ reltiming + phase + diagnosis+adj.fiber + (1|HM) + (1|plate), 
                                           subset(metadata.oars.stool.asv.double, response == "low")) %>% summary() %>% coef()
stats.starch.ph.high = lmerTest::lmer(scale(log10(starch)) ~ reltiming + phase + diagnosis+adj.fiber +(1|HM) + (1|plate),
                                            subset(metadata.oars.stool.asv.double, response == "high")) %>% summary() %>% coef()
stats.starch.ph.interact = lmerTest::lmer(scale(log10(starch)) ~ reltiming*response + phase + diagnosis+adj.fiber +(1|HM)+(1|plate), metadata.oars.stool.asv.double) %>% summary() %>% coef()

# make data.frame
stats.starch.ph = data.frame(Response = c("Strong Response", "Weak Response"),
                                   pval = c(stats.starch.ph.high[2,5], stats.starch.ph.low[2,5]))

metadata.oars.stool.asv.double$RS_Name = factor(metadata.oars.stool.asv.double$RS_Name, levels=rs.names)

metadata.oars.stool.asv.double.starch.plot = ggplot(subset(metadata.oars.stool.asv.double, !hm_phase %in% c("HM0924_rs2", "HM0932_rs2", "HM0759_rs2")) %>%
                                                            mutate(HM_phase = paste(HM, phase, sep="_")) %>%
                                                            mutate(reltiming = factor(ifelse(reltiming == "pre", "Pre", "Post"), levels=c("Pre", "Post"))) %>%
                                                            mutate(Response = ifelse(response == "low", "Weak Response", "Strong Response")),
                                                          aes(x=reltiming, y=starch))+
  geom_boxplot(width=0.5)+
  scale_y_log10()+
  geom_line(aes(group=HM_phase), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(shape=21, aes(fill=RS_Name), size=2)+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  geom_text(data = stats.starch.ph, x=1.5, y=log10(max(na.omit(metadata.oars.stool.asv.double$starch))), 
            aes(label=paste("p:", round(pval, digits=3))), size=3.5)+
  facet_wrap(~Response)+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="",
       y="Starch CAZy Intensity")
metadata.oars.stool.asv.double.starch.plot
stats.starch.ph.interact # no interaction


# :: MPX CAZy Mucin ------------------------------------------------

# delta calculated with mucin

stats.mucin.ph.low = lmerTest::lmer(scale(log10(mucin)) ~ reltiming + phase + diagnosis+adj.fiber + (1|HM) + (1|plate), 
                                           subset(metadata.oars.stool.asv.double, response == "low")) %>% summary() %>% coef()
stats.mucin.ph.high = lmerTest::lmer(scale(log10(mucin)) ~ reltiming + phase + diagnosis+adj.fiber  + (1|HM) + (1|plate), 
                                            subset(metadata.oars.stool.asv.double, response == "high")) %>% summary() %>% coef()
stats.mucin.ph.interact = lmerTest::lmer(scale(log10(mucin)) ~ reltiming*response + phase + diagnosis+adj.fiber  + (1|HM) + (1|plate), metadata.oars.stool.asv.double) %>% summary() %>% coef()

# make data.frame
stats.mucin.ph = data.frame(Response = c("Strong Response", "Weak Response"),
                                   pval = c(stats.mucin.ph.high[2,5], stats.mucin.ph.low[2,5]))

metadata.oars.stool.asv.double$RS_Name = factor(metadata.oars.stool.asv.double$RS_Name, levels=rs.names)

metadata.oars.stool.asv.double.mucin.plot = ggplot(subset(metadata.oars.stool.asv.double, !hm_phase %in% c("HM0924_rs2", "HM0932_rs2", "HM0759_rs2")) %>%
                                                            mutate(HM_phase = paste(HM, phase, sep="_")) %>%
                                                            mutate(reltiming = factor(ifelse(reltiming == "pre", "Pre", "Post"), levels=c("Pre", "Post"))) %>%
                                                            mutate(Response = ifelse(response == "low", "Weak Response", "Strong Response")),
                                                          aes(x=reltiming, y=mucin))+
  geom_boxplot(width=0.5)+
  scale_y_log10()+
  geom_line(aes(group=HM_phase), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(shape=21, aes(fill=RS_Name), size=2)+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  geom_text(data = stats.mucin.ph, x=1.5, y=log10(max(na.omit(metadata.oars.stool.asv.double$mucin))), 
            aes(label=paste("p:", round(pval, digits=3))), size=3.5)+
  facet_wrap(~Response)+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="",
       y="Mucin CAZy Intensity")
metadata.oars.stool.asv.double.mucin.plot
stats.mucin.ph.interact # no interaction


# :: MPX CAZy Starch:Mucin ------------------------------------------------


# delta calculated with starch.mucin

stats.starch.mucin.ph.low = lmerTest::lmer(scale(starch.mucin) ~ reltiming + phase + diagnosis+adj.fiber + (1|HM) + (1|plate),  
                                      subset(metadata.oars.stool.asv.double, response == "low")) %>% summary() %>% coef()
stats.starch.mucin.ph.high = lmerTest::lmer(scale(starch.mucin) ~ reltiming + phase + diagnosis+adj.fiber  + (1|HM) + (1|plate), 
                                       subset(metadata.oars.stool.asv.double, response == "high")) %>% summary() %>% coef()
stats.starch.mucin.ph.interact = lmerTest::lmer(scale(starch.mucin) ~ reltiming*response + phase + diagnosis+adj.fiber + (1|HM) + (1|plate), metadata.oars.stool.asv.double) %>% summary() %>% coef()

# make data.frame
stats.starch.mucin.ph = data.frame(Response = c("Strong Response", "Weak Response"),
                              pval = c(stats.starch.mucin.ph.high[2,5], stats.starch.mucin.ph.low[2,5]))

metadata.oars.stool.asv.double$RS_Name = factor(metadata.oars.stool.asv.double$RS_Name, levels=rs.names)

metadata.oars.stool.asv.double.starch.mucin.plot = ggplot(subset(metadata.oars.stool.asv.double, !hm_phase %in% c("HM0924_rs2", "HM0932_rs2", "HM0759_rs2")) %>%
                                                       mutate(HM_phase = paste(HM, phase, sep="_")) %>%
                                                       mutate(reltiming = factor(ifelse(reltiming == "pre", "Pre", "Post"), levels=c("Pre", "Post"))) %>%
                                                       mutate(Response = ifelse(response == "low", "Weak Response", "Strong Response")),
                                                     aes(x=reltiming, y=starch.mucin))+
  geom_boxplot(width=0.5)+
  geom_line(aes(group=HM_phase), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(shape=21, aes(fill=RS_Name), size=2)+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  geom_text(data = stats.starch.mucin.ph, x=1.5, y=(max(na.omit(metadata.oars.stool.asv.double$starch.mucin))), 
            aes(label=paste("p:", round(pval, digits=3))), size=3.5)+
  facet_wrap(~Response)+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="",
       y="Starch:Mucin CAZy Log2 Ratio")
metadata.oars.stool.asv.double.starch.mucin.plot
stats.starch.mucin.ph.interact # no interaction


# :: MBX lmer ------------------------------------------------------------------

dim(oars.mbx.prep) # 80% per group prevalence

## LMER
oars.mbx.lmer = delta.omic.lmer(
  data = oars.mbx.prep,
  split=F,
  mpx=F
)
oars.mbx.lmer %>% arrange(pval) %>% head()

oars.mbx.lmer$feature = oars.mbx.lmer$taxa

# volcano plot (FDR)
oars.mbx.lmer.volcano = ggplot(oars.mbx.lmer,
                               aes(x=estimate, y=padj))+
  geom_point(shape=21, aes(fill=estimate), size=2)+
  geom_hline(yintercept=(0.20), linetype=2, alpha=0.5)+
  geom_vline(xintercept=c(-1,1), linetype=2, alpha=0.5)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.05,0.1, 0.20, 0.5, 1))+
  #ggrepel::geom_text_repel(aes(label=ifelse(pval < 0.05, feature, NA)), size=3)+
  scale_fill_gradient2(low="blue", high="red")+
  facet_wrap(~"Interaction: Metabolite")+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10))+
  labs(x="Interaction Coefficient", y="FDR")
oars.mbx.lmer.volcano
# none

# volcano plot (unadjusted)
oars.mbx.lmer.volcano.unadj = ggplot(oars.mbx.lmer,
                                     aes(x=estimate, y=pval))+
  geom_point(shape=21, aes(fill=estimate), size=2)+
  geom_hline(yintercept=(0.05), linetype=2, alpha=0.5)+
  geom_vline(xintercept=c(-1,1), linetype=2, alpha=0.5)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05, 1))+
  ggrepel::geom_text_repel(aes(label=ifelse(pval < 0.05, feature, NA)), size=3)+
  scale_fill_gradient2(low="blue", high="red")+
  facet_wrap(~"Metabolite")+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10))+
  labs(x="Interaction Coefficient", y="p-value")
oars.mbx.lmer.volcano.unadj

oars.mbx.lmer.heatmap = ggplot(subset(oars.mbx.lmer,pval<0.05),
                               aes(x=1, y=reorder(taxa, estimate)))+
  geom_tile(aes(fill=estimate), color="black")+
  scale_fill_gradient2(low="blue", high="red")+
  theme_minimal()+theme(axis.text.x = element_blank(),
                        plot.title = element_text(hjust = 0.5, size=12))+#theme(legend.position="none")+
  labs(x="", y="", title="Metabolite", fill="adjusted\ninteraction\ncoefficient")
oars.mbx.lmer.heatmap


## Repeat with split
oars.mbx.lmer.split = delta.omic.lmer(
  split = TRUE,
  data = oars.mbx.prep
)

oars.mbx.lmer.split$feature = oars.mbx.lmer.split$taxa

oars.mbx.lmer.split.volcano = ggplot(oars.mbx.lmer.split %>%
                                       mutate(Response = ifelse(response == "high", "Strong", "Weak")),
                                     aes(x=estimate, y=padj))+
  geom_point(shape=21, aes(fill=estimate), size=2)+
  geom_hline(yintercept=0.20, linetype=2, alpha=0.5)+
  geom_vline(xintercept=c(-1,1), linetype=2, alpha=0.5)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.05,0.1, 0.20, 0.5, 1))+
  ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, feature, NA)), size=3)+
  scale_fill_gradient2(low="blue", high="red")+
  facet_wrap(~Response)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="Adjusted Coefficient", y="FDR", title="Metabolite")
oars.mbx.lmer.split.volcano
# none

oars.mbx.lmer.split.volcano.unadj = ggplot(oars.mbx.lmer.split %>%
                                             mutate(Response = ifelse(response == "high", "Strong", "Weak")),
                                           aes(x=estimate, y=pval))+
  geom_point(shape=21, aes(fill=estimate), size=2)+
  geom_hline(yintercept=0.05, linetype=2, alpha=0.5)+
  geom_vline(xintercept=c(-1,1), linetype=2, alpha=0.5)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.001, 0.01, 0.05,1))+
  ggrepel::geom_text_repel(aes(label=ifelse(pval < 0.05, feature, NA)), size=3)+
  scale_fill_gradient2(low="blue", high="red")+
  facet_wrap(~Response)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="Adjusted Coefficient", y="p-value", title="Metabolite")
oars.mbx.lmer.split.volcano.unadj



# >> Heatmaps (Group) ----------------------------------------------------

# first, tabulate # sig features over treatment
subset(oars.asv.maaslin.036, padj < 0.20)
subset(oars.mgx.maaslin.036, padj < 0.20)
subset(oars.mpx.kegg.maaslin.036, padj < 0.20)
subset(oars.mpx.cog.maaslin.036, padj < 0.20)
subset(oars.mpx.cazy.maaslin.036, padj < 0.20)
subset(oars.mbx.maaslin.036, padj < 0.20)

# GOAL: Make Log2FC heatmap of FDR < 0.20 features

# Subset data to FDR < 0.20 features
oars.group.asv.lfc = oars.asv.data.glom.prep %>% reshape2::melt() %>%
  group_by(HM, phase, variable) %>%
  mutate(lfc = log2(value) - lag(log2(value))) %>%
  subset(reltiming == "post") %>%
  mutate(sample = paste(HM, timing, sep="_")) %>%
  subset(variable %in% subset(oars.asv.maaslin.036, padj < 0.20)$feature)

oars.group.mgx.lfc = oars.mgx.prep %>% reshape2::melt() %>%
  group_by(HM, phase, variable) %>%
  mutate(lfc = log2(value) - lag(log2(value))) %>%
  subset(reltiming == "post") %>%
  mutate(sample = paste(HM, timing, sep="_")) %>%
  subset(variable %in% subset(oars.mgx.maaslin.036, padj < 0.20)$feature)

oars.group.mpx.cog.lfc = oars.mpx.cog.prep %>% reshape2::melt() %>%
  group_by(HM, phase, variable) %>%
  mutate(lfc = log2(value) - lag(log2(value))) %>%
  mutate(lfc = ifelse(is.na(lfc), 0, lfc))%>%
  subset(reltiming == "post") %>%
  mutate(sample = paste(HM, timing, sep="_"))  %>%
  subset(variable %in% make.names(subset(oars.mpx.cog.maaslin.036, padj < 0.20)$feature))

oars.group.mpx.kegg.lfc = oars.mpx.kegg.prep %>% reshape2::melt() %>%
  group_by(HM, phase, variable) %>%
  mutate(lfc = log2(value) - lag(log2(value))) %>%
  mutate(lfc = ifelse(is.na(lfc), 0, lfc))%>%
  subset(reltiming == "post") %>%
  mutate(sample = paste(HM, timing, sep="_"))  %>%
  subset(variable %in% subset(oars.mpx.kegg.maaslin.036, padj < 0.20)$feature)

oars.group.mpx.cazy.lfc = oars.mpx.cazy.prep %>% reshape2::melt() %>%
  group_by(HM, phase, variable) %>%
  mutate(lfc = log2(value) - lag(log2(value))) %>%
  mutate(lfc = ifelse(is.na(lfc), 0, lfc))%>%
  subset(reltiming == "post") %>%
  mutate(sample = paste(HM, timing, sep="_"))  %>%
  subset(variable %in% subset(oars.mpx.cazy.maaslin.036, padj < 0.20)$feature)

oars.group.mbx.lfc = oars.mbx.prep %>% reshape2::melt() %>%
  group_by(HM, phase, variable) %>%
  mutate(lfc = log2(value) - lag(log2(value))) %>%
  mutate(lfc = ifelse(is.na(lfc), 0, lfc))%>%
  subset(reltiming == "post") %>%
  mutate(sample = paste(HM, timing, sep="_"))  %>%
  subset(variable %in% subset(oars.mbx.maaslin.036, padj < 0.20)$feature)

# merge (keep matching samples)
# skip ASV because no padj < 0.20

oars.group.mgx.mpx.mbx.lfc = rbind(
  oars.group.asv.lfc,
  oars.group.mgx.lfc,
  oars.group.mpx.kegg.lfc,
  oars.group.mpx.cog.lfc,
  oars.group.mpx.cazy.lfc,
  oars.group.mbx.lfc) %>% as.data.frame()
unique(oars.group.mgx.mpx.mbx.lfc$variable) %>% length()

# make mapping object
oars.group.omics.lfc.map = oars.mbx.prep %>% reshape2::melt() %>%
  mutate(sample = paste(HM, timing, sep="_")) %>%
  dplyr::select(c(sample, HM, timing)) %>% distinct()

rownames(oars.group.omics.lfc.map) = oars.group.omics.lfc.map$sample
oars.group.omics.lfc.map$sample = NULL

# optional: keep samples with overlapping omics features

# scale per variable
oars.group.mgx.mpx.mbx.lfc = oars.group.mgx.mpx.mbx.lfc %>% 
  group_by(variable) %>% 
  mutate(slfc = scale(lfc))
# add feature map
oars.group.mgx.mpx.mbx.lfc.map = data.frame(variable = oars.group.mgx.mpx.mbx.lfc[,c("variable")]) %>% distinct() %>%
  mutate(Datatype = ifelse(variable %in% oars.group.mgx.lfc$variable, "Species",
                           ifelse(variable %in% oars.group.mpx.kegg.lfc$variable, "Pathway", 
                                  ifelse(variable %in% oars.group.mpx.cog.lfc$variable, "COG", 
                                         ifelse(variable %in% oars.group.mpx.cazy.lfc$variable, "CAZy", 
                                                ifelse(variable %in% oars.group.mbx.lfc$variable, "Metabolite", "ASV")))))) %>% data.frame()
rownames(oars.group.mgx.mpx.mbx.lfc.map) = oars.group.mgx.mpx.mbx.lfc.map$variable
colnames(oars.group.mgx.mpx.mbx.lfc.map)[1] = "feature"


# fix names in map and data
rownames(oars.group.mgx.mpx.mbx.lfc.map) = oars.group.mgx.mpx.mbx.lfc.map$feature

oars.group.mgx.mpx.mbx.lfc.map.map = merge(oars.group.mgx.mpx.mbx.lfc.map,
                                       data.frame(rbind(oars.asv.maaslin.036,
                                                        oars.mgx.maaslin.036,
                                                        oars.mpx.kegg.maaslin.036,
                                                        oars.mpx.cog.maaslin.036,
                                                        oars.mpx.cazy.maaslin.036,
                                                        oars.mbx.maaslin.036))[,c("feature", "coef")], by="feature")
rownames(oars.group.mgx.mpx.mbx.lfc.map.map) = oars.group.mgx.mpx.mbx.lfc.map.map$feature
oars.group.mgx.mpx.mbx.lfc.map.map$feature = NULL
colnames(oars.group.mgx.mpx.mbx.lfc.map.map)[2] = "Coefficient"

# clean feature names (if needed)
rownames(oars.group.mgx.mpx.mbx.lfc.map.map) = gsub(" involved.*", "", rownames(oars.group.mgx.mpx.mbx.lfc.map.map))
oars.group.mgx.mpx.mbx.lfc$variable = gsub(" involved.*", "", oars.group.mgx.mpx.mbx.lfc$variable)

# add HM map (RS_Name) (use from ASV)
oars.omics.heatmap.mapping = oars.asv.lfc.treatment.prev.mapping
oars.omics.heatmap.mapping = oars.omics.heatmap.mapping[rownames(oars.omics.heatmap.mapping) %in%
                                                          paste(oars.group.mgx.mpx.mbx.lfc$HM, 
                                                                oars.group.mgx.mpx.mbx.lfc$timing),]
#rownames(oars.omics.heatmap.mapping) = gsub(" ", "_", rownames(oars.omics.heatmap.mapping))

pheatmap::pheatmap(reshape2::acast(oars.group.mgx.mpx.mbx.lfc %>% mutate(sample2 = paste(HM, timing)),
                                   variable ~ sample2, value.var="slfc"),
                   color=colorRampPalette(c("blue","white", "red"))(100),
                   #angle_col = 45,
                   fontsize_row = 9,
                   fontsize_col = 8,
                   annotation_row = oars.group.mgx.mpx.mbx.lfc.map.map,
                   annotation_col = oars.omics.heatmap.mapping[,colnames(oars.omics.heatmap.mapping)!="Patient"],
                   annotation_colors = list(Coefficient = colorRampPalette(c("blue","white", "red"))(100),
                                            Datatype = omics.colors,
                                            Timing = c("3M" = "white", "6M" = "lightgrey"),
                                            RS = rs.colors),
                   breaks=c(seq(min(oars.group.mgx.mpx.mbx.lfc$slfc), 0, length.out=ceiling(100/2) + 1), 
                            seq(max(oars.group.mgx.mpx.mbx.lfc$slfc)/100, max(oars.group.mgx.mpx.mbx.lfc$slfc), length.out=floor(100/2))))

# LFC to previous timepoint


# >> MOFA // defunct -----------------------------------------------------------------

library("MOFA2")

# create list of data with shared samples
oars.asv.data.glom
oars.mgx.taxa.1
oars.mpx.kegg.mat.filt.1
oars.mpx.cog.mat.filt.1
oars.mpx.cazy.mat.filt.1
oars.mbx.raw.mat.filt.1

#
oars.mofa.data = 
  rbind(reshape2::melt(oars.asv.data.glom) %>% 
          mutate(value = log2(value + min(oars.asv.data.median.1[oars.asv.data.median.1!=0])/2)) %>%
          mutate(view = "ASV"),
        reshape2::melt(oars.mgx.taxa.1) %>% 
          mutate(value = log2(value + min(oars.mgx.taxa.1[oars.mgx.taxa.1!=0])/2)) %>%
          mutate(view = "Species"),
        reshape2::melt(as.matrix(oars.mpx.kegg.mat.filt.1)) %>% 
          mutate(value = log2(value + min(oars.mpx.kegg.mat.filt.1[oars.mpx.kegg.mat.filt.1!=0])/2)) %>%
          mutate(view = "Pathway"),
        reshape2::melt(as.matrix(oars.mpx.cog.mat.filt.1)) %>% 
          mutate(value = log2(value + min(oars.mpx.cog.mat.filt.1[oars.mpx.cog.mat.filt.1!=0])/2)) %>%
          mutate(view = "COG"),
        reshape2::melt(as.matrix(oars.mpx.cazy.mat.filt.1)) %>% 
          mutate(value = log2(value + min(oars.mpx.cazy.mat.filt.1[oars.mpx.cazy.mat.filt.1!=0])/2)) %>%
          mutate(view = "CAZy"),
        reshape2::melt(as.matrix((oars.mbx.raw.mat.filt.1))) %>% 
          mutate(Var2 = make.names(Var2)) %>%
          mutate(value = log2(value + min(oars.mbx.raw.mat.filt.1[oars.mbx.raw.mat.filt.1!=0])/2)) %>%
          mutate(view = "Metabolite")) %>%
  as.data.frame()
colnames(oars.mofa.data) = c("sample", "feature", "value", "view")      

oars.mofa.data = oars.mofa.data %>%
  mutate(feature = gsub("Â", "A", feature))

# make sample data
oars.mofa.meta = metadata.oars.stool.asv
oars.mofa.meta$sample = rownames(oars.mofa.meta)

# Follow tutorial: https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/microbiome_vignette.html#create-mofa-object

oars.mofa.object <- MOFA2::create_mofa(oars.mofa.data)

MOFA2::plot_data_overview(oars.mofa.object)

model_opts <- MOFA2::get_default_model_options(oars.mofa.object)
model_opts$num_factors <- 10

oars.mofa.object <- MOFA2::prepare_mofa(oars.mofa.object, model_options = model_opts)

set.seed(25)
oars.mofa.trained <- MOFA2::run_mofa(oars.mofa.object, use_basilisk = F)

MOFA2::plot_variance_explained(oars.mofa.trained, plot_total = T)[[2]]

MOFA2::plot_variance_explained(oars.mofa.trained, max_r2=15)
# 1 and 2 look good

# add metadata
MOFA2::samples_metadata(oars.mofa.trained) <- subset(oars.mofa.meta, sample %in% oars.mofa.trained@samples_metadata$sample)
MOFA2::samples_metadata(oars.mofa.trained)$RS_Name = factor(MOFA2::samples_metadata(oars.mofa.trained)$RS_Name, levels=rs.names)

# screeplot
MOFA2::calculate_variance_explained(oars.mofa.trained) %>% 
  .$r2_per_factor %>% 
  as.data.frame() %>% 
  rownames_to_column("factor") %>% 
  pivot_longer(-factor, names_to = "view", values_to = "R2") %>%
  group_by(factor) %>%
  summarise(total_R2 = sum(R2)) %>%
  mutate(factor = factor(factor, levels=c("Factor1", "Factor2", "Factor3", "Factor4", "Factor5", "Factor6", "Factor7", "Factor8", "Factor9", "Factor10")))%>%
  ggplot(aes(factor, total_R2)) +
  geom_col(fill = "steelblue", color="black") +
  theme_classic() +
  labs(y = "Total variance explained (%)", x = NULL)

# create dataframe for plotting
oars.mofa.trained.df <- MOFA2::get_factors(oars.mofa.trained, factors = "all") %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>% 
  merge(oars.mofa.meta, by="sample")

colnames(oars.mofa.trained.df) = gsub("single_group.Factor", "fac", colnames(oars.mofa.trained.df))

oars.mofa.trained.plot = ggplot(oars.mofa.trained.df%>%
         mutate(RS_Name = factor(RS_Name, levels=rs.names))%>%
           arrange(oars.days),
       aes(x=fac1,
           y=fac2))+
  geom_line(aes(group=HM), linewidth=0.2, linetype=2)+
  geom_point(aes(fill=RS_Name, shape=baseline),
             size=2.5)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values = c(rs.colors), na.value="lightgrey")+
  guides(shape = "none",
         fill = guide_legend(override.aes = list(shape=21)))+
  theme_classic()+
  facet_wrap(~"MOFA")+
  labs(x = "Factor 1", y = "Factor 2")
oars.mofa.trained.plot

# compare euclidean distance from baseline
oars.mofa.trained.df %>%
  group_by(HM) %>%
  mutate()

# permanova
set.seed(25)
t1 = Sys.time()
oars.mofa.permanova = vegan::adonis2(dist(oars.mofa.trained.df[,c("fac1", "fac2")]) ~ on.rs + adj.fiber + HM,
                                        oars.mofa.trained.df %>% mutate(on.rs = as.factor(ifelse(oars.on.rs=="onRS", "onRS", "offRS"))),
                                        strata = oars.mofa.trained.df$HM,
                                        by="margin")
t2 = Sys.time()
t2 - t1

oars.mofa.permanova

# MOFA does the poorest job of partitioning variance to treatment



# >>> 4. Response Tree lmer --------------------------------------------------

# Premise: RS must first be degraded (A) before it can be fermented (B).
# This (A) requires CAZymes.
# So, to control for the fact that some patients might not
# degrade the RS (A), it's unlikely they could ferment it (B).
# So, let's control for these non-digesters and ask:
# "which features correlate with starch CAZy"

# Main method: MLR with random slopes per microbiome (use maaslin/gam data)

# ** use compliant samples **

# :: filter ---------------------------------------------------------------

# Here, we'll apply a 10% prevalence filter (and 80% for MBX)

# ASV
dim(oars.asv.data.glom) # 1554
oars.asv.data.glom.filt.3 = oars.asv.data.glom[rownames(oars.asv.data.glom) %in% subset(metadata.oars.stool, compliant==T)$standard.name,]
oars.asv.data.glom.filt.3.pa = oars.asv.data.glom.filt.3
oars.asv.data.glom.filt.3.pa[oars.asv.data.glom.filt.3.pa!=0] = 1
oars.asv.data.glom.filt.3 = oars.asv.data.glom.filt.3[,colSums(oars.asv.data.glom.filt.3.pa) >= nrow(oars.asv.data.glom.filt.3)*0.1]
dim(oars.asv.data.glom.filt.3) # 443
oars.asv.data.glom.filt.3 = oars.asv.data.glom.filt.3/50000

# MGX
dim(oars.mgx.taxa) # 958
oars.mgx.taxa.filt.3 = oars.mgx.taxa[rownames(oars.mgx.taxa) %in% subset(metadata.oars.stool, compliant==T)$standard.name,]
oars.mgx.taxa.filt.3.pa = oars.mgx.taxa.filt.3
oars.mgx.taxa.filt.3.pa[oars.mgx.taxa.filt.3.pa!=0] = 1
oars.mgx.taxa.filt.3 = oars.mgx.taxa.filt.3[,colSums(oars.mgx.taxa.filt.3.pa) >= nrow(oars.mgx.taxa.filt.3)*0.1]
dim(oars.mgx.taxa.filt.3) # 542

# KEGG
dim(oars.mpx.kegg.mat) # 180
oars.mpx.kegg.filt.3 = oars.mpx.kegg.mat[rownames(oars.mpx.kegg.mat) %in% subset(metadata.oars.stool, compliant==T)$standard.name,]
oars.mpx.kegg.filt.3.pa = oars.mpx.kegg.filt.3
oars.mpx.kegg.filt.3.pa[is.na(oars.mpx.kegg.filt.3.pa)] = 0
oars.mpx.kegg.filt.3.pa[oars.mpx.kegg.filt.3.pa!=0] = 1
oars.mpx.kegg.filt.3=oars.mpx.kegg.filt.3
oars.mpx.kegg.filt.3[is.na(oars.mpx.kegg.filt.3)] = 0
oars.mpx.kegg.filt.3 = oars.mpx.kegg.filt.3[,colSums(oars.mpx.kegg.filt.3.pa) >= nrow(oars.mpx.kegg.filt.3)*0.1]
dim(oars.mpx.kegg.filt.3) # 180

# COG
dim(oars.mpx.cog.mat) # 2653
oars.mpx.cog.filt.3 = oars.mpx.cog.mat[rownames(oars.mpx.cog.mat) %in% subset(metadata.oars.stool, compliant==T)$standard.name,]
oars.mpx.cog.filt.3.pa = oars.mpx.cog.filt.3
oars.mpx.cog.filt.3.pa[is.na(oars.mpx.cog.filt.3.pa)] = 0
oars.mpx.cog.filt.3.pa[oars.mpx.cog.filt.3.pa!=0] = 1
oars.mpx.cog.filt.3=oars.mpx.cog.filt.3
oars.mpx.cog.filt.3[is.na(oars.mpx.cog.filt.3)] = 0
oars.mpx.cog.filt.3 = oars.mpx.cog.filt.3[,colSums(oars.mpx.cog.filt.3.pa) >= nrow(oars.mpx.cog.filt.3)*0.1]
dim(oars.mpx.cog.filt.3) # 2484

# CAZy
dim(oars.mpx.cazy.mat) # 66
oars.mpx.cazy.filt.3 = oars.mpx.cazy.mat[rownames(oars.mpx.cazy.mat) %in% subset(metadata.oars.stool, compliant==T)$standard.name,]
oars.mpx.cazy.filt.3.pa = oars.mpx.cazy.filt.3
oars.mpx.cazy.filt.3.pa[is.na(oars.mpx.cazy.filt.3.pa)] = 0
oars.mpx.cazy.filt.3.pa[oars.mpx.cazy.filt.3.pa!=0] = 1
oars.mpx.cazy.filt.3=oars.mpx.cazy.filt.3
oars.mpx.cazy.filt.3[is.na(oars.mpx.cazy.filt.3)] = 0
oars.mpx.cazy.filt.3 = oars.mpx.cazy.filt.3[,colSums(oars.mpx.cazy.filt.3.pa) >= nrow(oars.mpx.cazy.filt.3)*0.1]
dim(oars.mpx.cazy.filt.3) # 66

# MBX
dim(oars.mbx.raw.mat) # 2012
oars.mbx.filt.3 = oars.mbx.raw.mat[rownames(oars.mbx.raw.mat) %in% subset(metadata.oars.stool, compliant==T)$standard.name,]
oars.mbx.filt.3.pa = oars.mbx.filt.3
oars.mbx.filt.3.pa[oars.mbx.filt.3.pa!=0] = 1
oars.mbx.filt.3 = oars.mbx.filt.3[,colSums(oars.mbx.filt.3.pa) >= nrow(oars.mbx.filt.3)*0.8]
dim(oars.mbx.filt.3) # 153

oars.asv.data.glom.filt.3
oars.mgx.taxa.filt.3
oars.mpx.kegg.filt.3
oars.mpx.cog.filt.3
oars.mpx.cazy.filt.3
oars.mbx.filt.3

# :: process --------------------------------------------------------------

# add min!=0 / 2 pseudocount

oars.double.response.tree = merge(as.data.frame(oars.asv.data.glom.filt.3+min(oars.asv.data.glom.filt.3[oars.asv.data.glom.filt.3!=0])/2) %>%
                                    mutate(standard.name = rownames(.)),
                                  metadata.oars.stool.asv%>%subset(timing %in% c("0M", "3M", "6M")) %>% select(standard.name, HM, oars.days, starch, compliant, fcal, plate),
                                  by="standard.name")
oars.double.response.tree = merge(as.data.frame(oars.mgx.taxa.filt.3+min(oars.mgx.taxa.filt.3[oars.mgx.taxa.filt.3!=0])/2) %>%
                                    mutate(standard.name = rownames(.)),
                                  oars.double.response.tree,
                                  by="standard.name")
oars.double.response.tree = merge(as.data.frame(oars.mpx.kegg.filt.3+min(oars.mpx.kegg.filt.3[oars.mpx.kegg.filt.3!=0])/2) %>%
                                    mutate(standard.name = rownames(.)),
                                  oars.double.response.tree,
                                  by="standard.name")
oars.double.response.tree = merge(as.data.frame(oars.mpx.cog.filt.3+min(oars.mpx.cog.filt.3[oars.mpx.cog.filt.3!=0])/2) %>%
                                    mutate(standard.name = rownames(.)),
                                  oars.double.response.tree,
                                  by="standard.name")
oars.double.response.tree = merge(as.data.frame(oars.mpx.cazy.filt.3+min(oars.mpx.cazy.filt.3[oars.mpx.cazy.filt.3!=0])/2) %>%
                                    mutate(standard.name = rownames(.)),
                                  oars.double.response.tree,
                                  by="standard.name")
oars.double.response.tree = merge(as.data.frame(oars.mbx.filt.3+min(oars.mbx.filt.3[oars.mbx.filt.3!=0])/2) %>%
                                    mutate(standard.name = rownames(.)),
                                  oars.double.response.tree,
                                  by="standard.name")

dim(oars.double.response.tree)
# 3961 features

# :: mlr ------------------------------------------------------------------

oars.double.response.tree.lm = do.call(rbind, lapply(colnames(oars.double.response.tree)[!colnames(oars.double.response.tree)%in%c("standard.name", "oars.days", "HM", "starch", "plate")], function(x){
  print(x)
  data.subset = oars.double.response.tree[colnames(oars.double.response.tree)%in%c("standard.name", "oars.days", "HM", "starch", "plate")]
  data.subset$feature = oars.double.response.tree[,x]
  # process (log2 transform)
  data.subset$starch = log2(data.subset$starch)
  data.subset$feature = log2(data.subset$feature)
  # skip if feature is not unique in < 3 samples
  if(sum(data.subset$feature!=min(data.subset$feature))>=3){
    # build model (predict starch with feature, random SLOPES)
    # add plate random effect for MPX
    if(x %in% c(colnames(oars.mpx.kegg.filt.3), colnames(oars.mpx.cog.filt.3), colnames(oars.mpx.cazy.filt.3))){
    lm.output = lmerTest::lmer(scale(starch) ~ scale(feature) + (feature | HM) + (1|plate), data.subset)
    }else{
      # otherwise no additional covariate
    lm.output = lmerTest::lmer(scale(starch) ~ scale(feature) + (feature | HM), data.subset)
    }
    # extract slope coefficients (which are ~slopes relative to global slope)
    data.frame(slope = lme4::ranef(lm.output)$HM[,2],
               HM = rownames(lme4::ranef(lm.output)$HM)) %>%
      mutate(feature = x) %>%
      mutate(coef = coef(summary(lm.output))[2,1]) %>%
      mutate(pval = coef(summary(lm.output))[2,5])
  } else{
    data.frame(slope = NA,
               HM = NA,
               feature = x,
               coef = NA,
               pval = NA)
  }
}))

oars.double.response.tree.lm = subset(oars.double.response.tree.lm,
                                      !is.na(HM))

# sig models:  
oars.response.features = oars.double.response.tree.lm %>%
  select(feature, pval, coef) %>%
  distinct() %>%
  mutate(padj = p.adjust(pval, method="BH")) %>%
  arrange(padj) %>% subset(padj < 0.05) %>% select(feature) 

# variation

oars.response.features.slopes = oars.double.response.tree.lm %>%
  #select(feature, HM, slope) %>%
  #subset(feature %in% oars.response.features$feature) %>%
  group_by(feature) %>%
  mutate(slope.ratio = abs(table(sign(slope))[1]/table(sign(slope))[2]))%>%
  mutate(slope.sd = sd(slope),
         slope.min = range(slope)[1],
         slope.max = range(slope)[2],
         slope.range = range(slope)[2] - range(slope)[1],
         slope.mean = mean(slope)) %>%
  select(feature, pval, slope.sd,slope.min, slope.max, slope.range, slope.mean, slope.ratio, coef) %>% distinct() %>%
  mutate(padj = p.adjust(pval, method="BH"))%>%
  arrange(pval) 
oars.response.features.slopes %>% head(n=20)

# to plot a response tree you need to:
# 1. show strongest associations with Starch CAZy at root (e.g. GH13)
# 2. show correlated (low pval) features ("conserved co-responses", e.g. Nicotinate)
# 3. these transition to individualized impacts (i.e. variable slopes)
# Careful, because #3 could be noise, or strongly impacted by noise

oars.double.response.tree.lm.features = oars.double.response.tree.lm %>%
  subset(slope != 0)%>%
  group_by(feature) %>%
  mutate(slope.sd = sd(slope))%>%
  arrange(slope.sd) %>%
  dplyr::select(feature, slope.sd, pval) %>% distinct() %>%
  arrange(slope.sd)

# add omic type
oars.double.response.tree.lm = oars.double.response.tree.lm %>%
  mutate(Datatype = ifelse(feature %in% colnames(oars.asv.data.glom.filt.3), "ASV",
                           ifelse(feature %in% colnames(oars.mgx.taxa.filt.3), "Species",
                                  ifelse(feature %in% colnames(oars.mpx.kegg.filt.3), "Pathway",
                                         ifelse(feature %in% colnames(oars.mpx.cog.filt.3), "COG", 
                                                ifelse(feature %in% colnames(oars.mpx.cazy.filt.3), "CAZy", 
                                                       ifelse(feature %in% colnames(oars.mbx.filt.3), "Metabolite", "other")))))))

# first, pick features with most conserved associations with starch
oars.response.tree.conserved = oars.double.response.tree.lm.features %>%
  arrange(pval)%>%
  head(n=20)
# second, pick features with least conserved associations with starch
oars.response.tree.variable = oars.double.response.tree.lm.features %>%
  arrange(slope.sd) %>%
  tail(n=20)

table(oars.double.response.tree.lm$Datatype)


# :: tree plot -----------------------------------------------------------------


# show bottom 20 and top 20
oars.double.response.tree.lm.plot = ggplot(oars.double.response.tree.lm %>%
                                             subset(slope != 0)%>%
                                             group_by(feature) %>%
                                             mutate(slope.sd = sd(slope))%>%
                                             arrange(slope.sd) %>%
                                             subset(feature %in% c(oars.response.tree.conserved$feature,
                                                                   oars.response.tree.variable$feature))%>%
                                             mutate(stage = ifelse(feature %in% oars.response.tree.conserved$feature,
                                                                   "Conserved", "Variable"))%>%
                                             mutate(stage = factor(stage, levels=rev(c("Variable", "Conserved")))) %>%
                                             mutate(Datatype = factor(Datatype, levels=(names(omics.colors))))%>%
                                             #subset(feature %in% oars.response.features$feature)%>%
                                             arrange(-slope.sd),
                                           aes(x=slope, y=reorder(feature, -rank(slope.sd))))+
  geom_line(aes(group=feature))+
  geom_path(aes(group=HM), linetype=2, linewidth=0.2, color="black", alpha=0.5)+
  geom_point(aes(shape=Datatype, fill=Datatype), size=2.5)+
  #geom_point(aes(x=coef), shape=10)+ # plot overall model coef
  scale_fill_manual(values=omics.colors)+
  scale_shape_manual(values=omics.shapes)+
  theme_classic()+theme(legend.position.inside = c(0.9, 0.75),
                        legend.position="inside",
                        axis.title.y=element_blank())+
  facet_wrap(~ stage, nrow=2, scales="free_y")+
  labs(x="Patient Slope")
oars.double.response.tree.lm.plot
# this is it; just need to check for bugs

# filtered features
# built lmer with random slopes (treatment period, starch ~ feature + (feature | HM)
# "conserved" = top 15 lowest pvals (most consistent)
# "variable" = top 15 highest slope variance


# :: omics ----------------------------------------------------------------

# now we show that metabolites have more individualized responses

oars.double.response.tree.lm.sd = oars.double.response.tree.lm %>%
  subset(feature != "fcal") %>%
  subset(slope != 0)%>%
  group_by(feature) %>%
  mutate(feature.sd = sd(slope))%>%
  arrange(feature.sd) %>% 
  group_by(feature) %>%
  dplyr::select(feature, feature.sd, `Datatype`) %>% distinct() %>%
  # select top X based on parameter
  arrange(-feature.sd)

# stats
oars.double.response.tree.lm.sd


oars.double.response.tree.lm.sd.plot = ggplot(oars.double.response.tree.lm.sd %>%
         subset(feature != "fcal") %>%
         mutate(`Datatype` = factor(`Datatype`, levels=names(omics.colors))) %>%
           arrange(feature.sd),
       aes(x=`Datatype`, y=sqrt(feature.sd)))+
  geom_jitter(aes(fill=Datatype, shape=Datatype), color="black", 
              alpha=0.25, width=0.1, size=2)+
  geom_violin(aes(fill=`Datatype`), draw_quantiles = 0.5)+
  scale_y_sqrt(trans="reverse", position="right")+
  scale_x_discrete(position="top")+
  scale_fill_manual(values=omics.colors)+
  scale_color_manual(values=omics.colors)+
  scale_shape_manual(values=omics.shapes)+
  theme_classic()+theme(legend.position="none",
                        axis.text.x = element_text(angle=315, hjust=1),
                        axis.title.x=element_blank())+
  labs(y="√ Patient Slope SD") # √
oars.double.response.tree.lm.sd.plot
# don't need stats

# Let's plot Niconitate and Pantothenate
# and some highly variable ones


stats.tree.features.data = oars.double.response.tree[,c("Nicotinate and nicotinamide metabolism",
                                                        "Pantothenate and CoA biosynthesis",
                                                        "Vancomycin resistance",
                                                        tail(oars.response.tree.variable$feature, 2), "standard.name")]

# :: :: Nicotinamide ------------------------------------------------------

stats.tree.features.kegg.data = merge(oars.mpx.kegg.mat.filt.1[,c("Nicotinate and nicotinamide metabolism",
                                                          "Pantothenate and CoA biosynthesis",
                                                          "Vancomycin resistance",
                                                          "Glycerophospholipid metabolism")] %>% mutate(standard.name = rownames(.)),
                                 metadata.oars.stool.asv, by="standard.name")


# stats
stats.tree.feature.niconitate.036 = lmerTest::lmer(scale(log2(`Nicotinate and nicotinamide metabolism`)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM)+(1|plate),
                                        subset(stats.tree.features.kegg.data, timing %in% c("0M", "3M", "6M")& 
                                                 compliant == TRUE)) %>%
  summary() %>% coef()

stats.tree.feature.niconitate.6912 = lmerTest::lmer(scale(log2(`Nicotinate and nicotinamide metabolism`)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM)+(1|plate),
                                                   subset(stats.tree.features.kegg.data, timing %in% c("6M", "9M", "12M")& 
                                                            compliant == TRUE)) %>%
  summary() %>% coef()

# plot
stats.tree.features.kegg.data.niconitate.plot = ggplot(subset(stats.tree.features.kegg.data,compliant==TRUE),
                                               aes(x=oars.days, 
                                                   y=`Nicotinate and nicotinamide metabolism`))+
  annotate("rect", xmin=0, xmax=max(subset(metadata.oars.stool, oars.on.rs=="onRS")$oars.days),
           ymin=-Inf, ymax=Inf,
           alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(aes(fill=rs.col, shape=baseline), size=2)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  geom_smooth(color="black")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12))+
  geom_text(data=pval.positioner("Nicotinate and nicotinamide metabolism",
                                 subset(stats.tree.features.kegg.data, compliant == TRUE),
                                 stats.tree.feature.niconitate.036,
                                 stats.tree.feature.niconitate.6912,
                                 yadjust = 1.15), aes(x=xval, y=yval, label=label), vjust=1,size=4)+
  labs(x="Days since starting RS", y="Nicotinate and nicotinamide metabolism")
stats.tree.features.kegg.data.niconitate.plot


# :: :: Pantothenate ------------------------------------------------------


# stats
stats.tree.feature.pantothenate.036 = lmerTest::lmer(scale(log2(`Pantothenate and CoA biosynthesis`)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM)+(1|plate),
                                                   subset(stats.tree.features.kegg.data, timing %in% c("0M", "3M", "6M")& 
                                                            compliant == TRUE)) %>%
  summary() %>% coef()

stats.tree.feature.pantothenate.6912 = lmerTest::lmer(scale(log2(`Pantothenate and CoA biosynthesis`)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM)+(1|plate),
                                                    subset(stats.tree.features.kegg.data, timing %in% c("6M", "9M", "12M")& 
                                                             compliant == TRUE)) %>%
  summary() %>% coef()

# plot
stats.tree.features.kegg.data.pantothenate.plot = ggplot(subset(stats.tree.features.kegg.data,compliant==TRUE),
                                                       aes(x=oars.days, 
                                                           y=`Pantothenate and CoA biosynthesis`))+
  annotate("rect", xmin=0, xmax=max(subset(metadata.oars.stool, oars.on.rs=="onRS")$oars.days),
           ymin=-Inf, ymax=Inf,
           alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(aes(fill=rs.col, shape=baseline), size=2)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  geom_smooth(color="black")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12))+
  geom_text(data=pval.positioner("Pantothenate and CoA biosynthesis",
                                 subset(stats.tree.features.kegg.data, compliant == TRUE),
                                 stats.tree.feature.pantothenate.036,
                                 stats.tree.feature.pantothenate.6912,
                                 yadjust = 1.15), aes(x=xval, y=yval, label=label), vjust=1,size=4)+
  labs(x="Days since starting RS", y="Pantothenate and CoA biosynthesis")
stats.tree.features.kegg.data.pantothenate.plot




# :: :: Vancomycin resistance ------------------------------------------------------

# stats
stats.tree.feature.vancomycin.036 = lmerTest::lmer(scale(log2(`Vancomycin resistance`)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM)+(1|plate),
                                                     subset(stats.tree.features.kegg.data, timing %in% c("0M", "3M", "6M")& 
                                                              compliant == TRUE)) %>%
  summary() %>% coef()

stats.tree.feature.vancomycin.6912 = lmerTest::lmer(scale(log2(`Vancomycin resistance`)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM)+(1|plate),
                                                      subset(stats.tree.features.kegg.data, timing %in% c("6M", "9M", "12M")& 
                                                               compliant == TRUE)) %>%
  summary() %>% coef()

# plot
stats.tree.features.kegg.data.vancomycin.plot = ggplot(subset(stats.tree.features.kegg.data,compliant==TRUE),
                                                         aes(x=oars.days, 
                                                             y=`Vancomycin resistance`))+
  annotate("rect", xmin=0, xmax=max(subset(metadata.oars.stool, oars.on.rs=="onRS")$oars.days),
           ymin=-Inf, ymax=Inf,
           alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(aes(fill=rs.col, shape=baseline), size=2)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  geom_smooth(color="black")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12))+
  geom_text(data=pval.positioner("Vancomycin resistance",
                                 subset(stats.tree.features.kegg.data, compliant == TRUE),
                                 stats.tree.feature.vancomycin.036,
                                 stats.tree.feature.vancomycin.6912,
                                 yadjust = 1.15), aes(x=xval, y=yval, label=label), vjust=1,size=4)+
  labs(x="Days since starting RS", y="Vancomycin resistance")
stats.tree.features.kegg.data.vancomycin.plot



# :: :: Hexatheyleneglycol ------------------------------------------------------


stats.tree.features.metabolite.data = merge(oars.mbx.raw.mat.filt.1[,c("Hexaethyleneglycol mono-n-tetradecyl ether | 11.51min : 478.387667m/z",
                                                                       "Octaethyleneglycol monododecyl ether | 11.13min : 538.408992m/z",
                                                                       "13-Docosenamide | 11.69min : 337.334677m/z",
                                                                       "N-Stearoyl Valine | 11.69min : 383.340342m/z")] %>% mutate(standard.name = rownames(.)),
                                      metadata.oars.stool.asv, by="standard.name")

# stats
stats.tree.feature.hexatheyleneglycol.036 = lmerTest::lmer(scale(log2(`Hexaethyleneglycol mono-n-tetradecyl ether | 11.51min : 478.387667m/z`)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                                                           subset(stats.tree.features.metabolite.data, timing %in% c("0M", "3M", "6M")& 
                                                                    compliant == TRUE)) %>%
  summary() %>% coef()

stats.tree.feature.hexatheyleneglycol.6912 = lmerTest::lmer(scale(log2(`Hexaethyleneglycol mono-n-tetradecyl ether | 11.51min : 478.387667m/z`)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                                                            subset(stats.tree.features.metabolite.data, timing %in% c("6M", "9M", "12M")& 
                                                                     compliant == TRUE)) %>%
  summary() %>% coef()

# plot
stats.tree.features.kegg.data.hexatheyleneglycol.plot = ggplot(subset(stats.tree.features.metabolite.data,compliant==TRUE),
                                                               aes(x=oars.days, 
                                                                   y=`Hexaethyleneglycol mono-n-tetradecyl ether | 11.51min : 478.387667m/z`))+
  annotate("rect", xmin=0, xmax=max(subset(metadata.oars.stool, oars.on.rs=="onRS")$oars.days),
           ymin=-Inf, ymax=Inf,
           alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(aes(fill=rs.col, shape=baseline), size=2)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  geom_smooth(color="black")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12))+
  geom_text(data=pval.positioner("Hexaethyleneglycol mono-n-tetradecyl ether | 11.51min : 478.387667m/z",
                                 subset(stats.tree.features.metabolite.data, compliant == TRUE),
                                 stats.tree.feature.hexatheyleneglycol.036,
                                 stats.tree.feature.hexatheyleneglycol.6912,
                                 yadjust = 1.15), aes(x=xval, y=yval, label=label), vjust=1,size=4)+
  labs(x="Days since starting RS", y="Hexaethyleneglycol...")
stats.tree.features.kegg.data.hexatheyleneglycol.plot


# :: :: Stearoyl Valine ------------------------------------------------------


# stats
stats.tree.feature.stearoylvaline.036 = lmerTest::lmer(scale(log2(`N-Stearoyl Valine | 11.69min : 383.340342m/z`)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                                                            subset(stats.tree.features.metabolite.data, timing %in% c("0M", "3M", "6M")& 
                                                                     compliant == TRUE)) %>%
  summary() %>% coef()

stats.tree.feature.stearoylvaline.6912 = lmerTest::lmer(scale(log2(`N-Stearoyl Valine | 11.69min : 383.340342m/z`)) ~ scale(oars.days) + diagnosis + adj.fiber + (1|HM),
                                                             subset(stats.tree.features.metabolite.data, timing %in% c("6M", "9M", "12M")& 
                                                                      compliant == TRUE)) %>%
  summary() %>% coef()

# plot
stats.tree.features.kegg.data.stearoylvaline.plot = ggplot(subset(stats.tree.features.metabolite.data,compliant==TRUE),
                                                                aes(x=oars.days, 
                                                                    y=`N-Stearoyl Valine | 11.69min : 383.340342m/z`))+
  annotate("rect", xmin=0, xmax=max(subset(metadata.oars.stool, oars.on.rs=="onRS")$oars.days),
           ymin=-Inf, ymax=Inf,
           alpha=0.2, fill="salmon") +
  geom_line(aes(group=HM), linetype=2, alpha=0.5, linewidth=0.3)+
  geom_point(aes(fill=rs.col, shape=baseline), size=2)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  geom_smooth(color="black")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12))+
  geom_text(data=pval.positioner("N-Stearoyl Valine | 11.69min : 383.340342m/z",
                                 subset(stats.tree.features.metabolite.data, compliant == TRUE),
                                 stats.tree.feature.stearoylvaline.036,
                                 stats.tree.feature.stearoylvaline.6912,
                                 yadjust = 1.15), aes(x=xval, y=yval, label=label), vjust=1,size=4)+
  labs(x="Days since starting RS", y="N-Stearoyl Valine")
stats.tree.features.kegg.data.stearoylvaline.plot



# plots
stats.tree.features.plots = (stats.tree.features.kegg.data.niconitate.plot+
                               stats.tree.features.kegg.data.pantothenate.plot+
                               stats.tree.features.kegg.data.stearoylvaline.plot+
                               stats.tree.features.kegg.data.vancomycin.plot+
                               patchwork::plot_layout(nrow=1))


# :: plots ----------------------------------------------------------------

stats.tree.features.plots.top = oars.double.response.tree.lm.plot+
  patchwork::plot_spacer()+
  patchwork::free(oars.double.response.tree.lm.sd.plot,type="label")+
  patchwork::plot_layout(widths=c(1,0.3,1))

cowplot::plot_grid(stats.tree.features.plots.top, stats.tree.features.plots,
                   nrow=2, rel_heights=c(3,1))


# >>> 4 Dendrograms // defunct--------------------------------------

# Goal: Identify variables strongly affected per individual
# Method: Apply MLR with random slopes for all variables ~ oars.days
# (do not chase statistical significance with MLR, since underpowered)


# :: filter ---------------------------------------------------------------

# but we'll use ALL samples from compliant patients
# so we need to refilter
# we'll apply a 10% prevalence filter (and 80% for MBX)

phases = c("preRS", "onRS", "postRS")

# ASV
dim(oars.asv.data.glom) # 1554
oars.asv.data.glom.filt.4 = oars.asv.data.glom[rownames(oars.asv.data.glom) %in% subset(metadata.oars.stool, oars.on.rs %in% phases & compliant==T)$standard.name,]
oars.asv.data.glom.filt.4.pa = oars.asv.data.glom.filt.4
oars.asv.data.glom.filt.4.pa[oars.asv.data.glom.filt.4.pa!=0] = 1
oars.asv.data.glom.filt.4 = oars.asv.data.glom.filt.4[,colSums(oars.asv.data.glom.filt.4.pa) >= nrow(oars.asv.data.glom.filt.4)*0.1]
dim(oars.asv.data.glom.filt.4) # 406 x 50 samples
oars.asv.data.glom.filt.4 = oars.asv.data.glom.filt.4/50000

# MGX
dim(oars.mgx.taxa) # 958
oars.mgx.taxa.filt.4 = oars.mgx.taxa[rownames(oars.mgx.taxa) %in% subset(metadata.oars.stool,oars.on.rs %in% phases & compliant==T)$standard.name,]
oars.mgx.taxa.filt.4.pa = oars.mgx.taxa.filt.4
oars.mgx.taxa.filt.4.pa[oars.mgx.taxa.filt.4.pa!=0] = 1
oars.mgx.taxa.filt.4 = oars.mgx.taxa.filt.4[,colSums(oars.mgx.taxa.filt.4.pa) >= nrow(oars.mgx.taxa.filt.4)*0.1]
dim(oars.mgx.taxa.filt.4) # 489 x 48 samples

# KEGG
dim(oars.mpx.kegg.mat) # 180
oars.mpx.kegg.filt.4 = oars.mpx.kegg.mat[rownames(oars.mpx.kegg.mat) %in% subset(metadata.oars.stool, oars.on.rs %in%phases & compliant==T)$standard.name,]
oars.mpx.kegg.filt.4.pa = oars.mpx.kegg.filt.4
oars.mpx.kegg.filt.4.pa[is.na(oars.mpx.kegg.filt.4.pa)] = 0
oars.mpx.kegg.filt.4.pa[oars.mpx.kegg.filt.4.pa!=0] = 1
oars.mpx.kegg.filt.4=oars.mpx.kegg.filt.4
oars.mpx.kegg.filt.4[is.na(oars.mpx.kegg.filt.4)] = 0
oars.mpx.kegg.filt.4 = oars.mpx.kegg.filt.4[,colSums(oars.mpx.kegg.filt.4.pa) >= nrow(oars.mpx.kegg.filt.4)*0.1]
dim(oars.mpx.kegg.filt.4) # 174 x 48 samples

# COG
dim(oars.mpx.cog.mat) # 2653
oars.mpx.cog.filt.4 = oars.mpx.cog.mat[rownames(oars.mpx.cog.mat) %in% subset(metadata.oars.stool, oars.on.rs %in% phases& compliant==T)$standard.name,]
oars.mpx.cog.filt.4.pa = oars.mpx.cog.filt.4
oars.mpx.cog.filt.4.pa[is.na(oars.mpx.cog.filt.4.pa)] = 0
oars.mpx.cog.filt.4.pa[oars.mpx.cog.filt.4.pa!=0] = 1
oars.mpx.cog.filt.4=oars.mpx.cog.filt.4
oars.mpx.cog.filt.4[is.na(oars.mpx.cog.filt.4)] = 0
oars.mpx.cog.filt.4 = oars.mpx.cog.filt.4[,colSums(oars.mpx.cog.filt.4.pa) >= nrow(oars.mpx.cog.filt.4)*0.1]
dim(oars.mpx.cog.filt.4) # 2454 x 48 samples

# CAZy
dim(oars.mpx.cazy.mat) # 66
oars.mpx.cazy.filt.4 = oars.mpx.cazy.mat[rownames(oars.mpx.cazy.mat) %in% subset(metadata.oars.stool, oars.on.rs %in% phases & compliant==T)$standard.name,]
oars.mpx.cazy.filt.4.pa = oars.mpx.cazy.filt.4
oars.mpx.cazy.filt.4.pa[is.na(oars.mpx.cazy.filt.4.pa)] = 0
oars.mpx.cazy.filt.4.pa[oars.mpx.cazy.filt.4.pa!=0] = 1
oars.mpx.cazy.filt.4=oars.mpx.cazy.filt.4
oars.mpx.cazy.filt.4[is.na(oars.mpx.cazy.filt.4)] = 0
oars.mpx.cazy.filt.4 = oars.mpx.cazy.filt.4[,colSums(oars.mpx.cazy.filt.4.pa) >= nrow(oars.mpx.cazy.filt.4)*0.1]
dim(oars.mpx.cazy.filt.4) # 65 x 48 samples

# MBX
dim(oars.mbx.raw.mat) # 2012
oars.mbx.filt.4 = oars.mbx.raw.mat[rownames(oars.mbx.raw.mat) %in% subset(metadata.oars.stool, oars.on.rs %in% phases & compliant==T)$standard.name,]
oars.mbx.filt.4.pa = oars.mbx.filt.4
oars.mbx.filt.4.pa[oars.mbx.filt.4.pa!=0] = 1
oars.mbx.filt.4 = oars.mbx.filt.4[,colSums(oars.mbx.filt.4.pa) >= nrow(oars.mbx.filt.4)*0.8]
dim(oars.mbx.filt.4) # 151 x 46 samples

oars.asv.data.glom.filt.4
oars.mgx.taxa.filt.4
oars.mpx.kegg.filt.4
oars.mpx.cog.filt.4
oars.mpx.cazy.filt.4
oars.mbx.filt.4

# :: process --------------------------------------------------------------


oars.individualized.lm.data = merge(as.data.frame(oars.asv.data.glom.filt.4+min(oars.asv.data.glom.filt.4[oars.asv.data.glom.filt.4!=0])/2) %>%
                                    mutate(standard.name = rownames(.)),
                                  metadata.oars.stool.asv%>%subset(oars.on.rs %in% phases) %>% select(standard.name, HM, oars.days, starch, compliant, fcal, oars.on.rs),
                                  by="standard.name")
oars.individualized.lm.data = merge(as.data.frame(oars.mgx.taxa.filt.4+min(oars.mgx.taxa.filt.4[oars.mgx.taxa.filt.4!=0])/2) %>%
                                    mutate(standard.name = rownames(.)),
                                  oars.individualized.lm.data,
                                  by="standard.name")
oars.individualized.lm.data = merge(as.data.frame(oars.mpx.kegg.filt.4+min(oars.mpx.kegg.filt.4[oars.mpx.kegg.filt.4!=0])/2) %>%
                                    mutate(standard.name = rownames(.)),
                                  oars.individualized.lm.data,
                                  by="standard.name")
oars.individualized.lm.data = merge(as.data.frame(oars.mpx.cog.filt.4+min(oars.mpx.cog.filt.4[oars.mpx.cog.filt.4!=0])/2) %>%
                                    mutate(standard.name = rownames(.)),
                                  oars.individualized.lm.data,
                                  by="standard.name")
oars.individualized.lm.data = merge(as.data.frame(oars.mpx.cazy.filt.4+min(oars.mpx.cazy.filt.4[oars.mpx.cazy.filt.4!=0])/2) %>%
                                    mutate(standard.name = rownames(.)),
                                  oars.individualized.lm.data,
                                  by="standard.name")
oars.individualized.lm.data = merge(as.data.frame(oars.mbx.filt.4+min(oars.mbx.filt.4[oars.mbx.filt.4!=0])/2) %>%
                                    mutate(standard.name = rownames(.)),
                                  oars.individualized.lm.data,
                                  by="standard.name")

dim(oars.individualized.lm.data)
# 3872 features x 42 samples


# :: dendrograms ----------------------------------------------------------

# use spearman distance for dendrograms

oars.spear.dist.asv = 1 - Hmisc::rcorr(t(as.data.frame(oars.asv.data.glom.filt.4)))$r
oars.spear.dist.mgx = 1 - Hmisc::rcorr(t(as.data.frame(oars.mgx.taxa.filt.4)))$r
oars.spear.dist.pathway = 1 - Hmisc::rcorr(t(as.data.frame(oars.mpx.kegg.filt.4)))$r
oars.spear.dist.cog = 1 - Hmisc::rcorr(t(as.data.frame(oars.mpx.cog.filt.4)))$r
oars.spear.dist.cazy = 1 - Hmisc::rcorr(t(as.data.frame(oars.mpx.cazy.filt.4)))$r
oars.spear.dist.mbx = 1 - Hmisc::rcorr(t(as.data.frame(oars.mbx.filt.4)))$r

# make meta
oars.spear.dist.meta = subset(metadata.oars.stool, oars.on.rs %in% phases & compliant == T)
oars.spear.dist.meta = oars.spear.dist.meta[,c("RS_Name", "timing","fcal", "adj.fiber")]
oars.spear.dist.meta$label = rownames(oars.spear.dist.meta)
oars.spear.dist.meta$RS_Name = factor(oars.spear.dist.meta$RS_Name, levels=rs.names)

# prepare data
oars.spear.dist.asv.dendro = ggdendro::dendro_data(as.dendrogram(hclust(dist(oars.spear.dist.asv, method="euclidean"), method = "ward.D2")), type = "rectangle")
oars.spear.dist.asv.dendro$labels = merge(oars.spear.dist.asv.dendro$labels, oars.spear.dist.meta)
# rescale
oars.spear.dist.asv.dendro.scaling.factor = max(oars.spear.dist.asv.dendro$segments$y, oars.spear.dist.asv.dendro$segments$yend, na.rm = TRUE)
oars.spear.dist.asv.dendro$segments$y <- oars.spear.dist.asv.dendro$segments$y / oars.spear.dist.asv.dendro.scaling.factor  # Normalize to [0, 1]
oars.spear.dist.asv.dendro$segments$yend <- oars.spear.dist.asv.dendro$segments$yend /oars.spear.dist.asv.dendro.scaling.factor  # Normalize endpoints
# plot
oars.spear.dist.asv.dendro.plot = ggplot() +
  geom_segment(data = (oars.spear.dist.asv.dendro$segments), aes(x = x, y = y, xend = xend, yend = yend)) +  # Dendrogram segments
  geom_text(data = (oars.spear.dist.asv.dendro$labels), aes(x = x, y = -0.45, label = substr(label, 1,6)),  # Labels with color by slope
            hjust = 1, size = 3) +
  geom_tile(data=(oars.spear.dist.asv.dendro$labels), aes(x = x, y = -0.05, fill = as.factor(RS_Name)),height=0.1,color="white")+
  geom_text(data=(oars.spear.dist.asv.dendro$labels), aes(x = x, y = -0.5, label=timing), hjust=0, size=3)+
  scale_fill_manual(values=rs.colors)+
  coord_flip() +  # Flip for vertical dendrogram
  scale_y_reverse(limits=c(1,-0.5), expand = c(0.2, 0)) +  # Reverse y-axis for traditional dendrogram look
  labs(title = "", x = "", y = "") +
  theme_void()+
  facet_wrap(~"ASV")+
  theme(legend.position="none",
        strip.text=element_text(size=10),
        axis.text.y = element_blank(),  # Remove y-axis labels
        axis.ticks.y = element_blank())
oars.spear.dist.asv.dendro.plot

# prepare data
oars.spear.dist.mgx.dendro = ggdendro::dendro_data(as.dendrogram(hclust(dist(oars.spear.dist.mgx, method="euclidean"), method = "ward.D2")), type = "rectangle")
oars.spear.dist.mgx.dendro$labels = merge(oars.spear.dist.mgx.dendro$labels, oars.spear.dist.meta)
# rescale
oars.spear.dist.mgx.dendro.scaling.factor = max(oars.spear.dist.mgx.dendro$segments$y, oars.spear.dist.mgx.dendro$segments$yend, na.rm = TRUE)
oars.spear.dist.mgx.dendro$segments$y <- oars.spear.dist.mgx.dendro$segments$y / oars.spear.dist.mgx.dendro.scaling.factor  # Normalize to [0, 1]
oars.spear.dist.mgx.dendro$segments$yend <- oars.spear.dist.mgx.dendro$segments$yend /oars.spear.dist.mgx.dendro.scaling.factor  # Normalize endpoints
# plot
oars.spear.dist.mgx.dendro.plot = ggplot() +
  geom_segment(data = (oars.spear.dist.mgx.dendro$segments), aes(x = x, y = y, xend = xend, yend = yend)) +  # Dendrogram segments
  geom_text(data = (oars.spear.dist.mgx.dendro$labels), aes(x = x, y = -0.45, label = substr(label, 1,6)),  # Labels with color by slope
            hjust = 1, size = 3) +
  geom_tile(data=(oars.spear.dist.mgx.dendro$labels), aes(x = x, y = -0.05, fill = as.factor(RS_Name)),height=0.1,color="white")+
  geom_text(data=(oars.spear.dist.mgx.dendro$labels), aes(x = x,  y = -0.5, label=timing), hjust=0, size=3)+
  scale_fill_manual(values=rs.colors)+
  coord_flip() +  # Flip for vertical dendrogram
  scale_y_reverse(limits=c(1,-0.5), expand = c(0.2, 0)) +  # Reverse y-axis for traditional dendrogram look
  labs(title = "", x = "", y = "") +
  theme_void()+
  facet_wrap(~"Species")+
  theme(legend.position="none",
        strip.text=element_text(size=10),
        axis.text.y = element_blank(),  # Remove y-axis labels
        axis.ticks.y = element_blank())
oars.spear.dist.mgx.dendro.plot

# prepare data
oars.spear.dist.pathway.dendro = ggdendro::dendro_data(as.dendrogram(hclust(dist(oars.spear.dist.pathway, method="euclidean"), method = "ward.D2")), type = "rectangle")
oars.spear.dist.pathway.dendro$labels = merge(oars.spear.dist.pathway.dendro$labels, oars.spear.dist.meta)
# rescale
oars.spear.dist.pathway.dendro.scaling.factor = max(oars.spear.dist.pathway.dendro$segments$y, oars.spear.dist.pathway.dendro$segments$yend, na.rm = TRUE)
oars.spear.dist.pathway.dendro$segments$y <- oars.spear.dist.pathway.dendro$segments$y / oars.spear.dist.pathway.dendro.scaling.factor  # Normalize to [0, 1]
oars.spear.dist.pathway.dendro$segments$yend <- oars.spear.dist.pathway.dendro$segments$yend /oars.spear.dist.pathway.dendro.scaling.factor  # Normalize endpoints
# plot
oars.spear.dist.pathway.dendro.plot = ggplot() +
  geom_segment(data = (oars.spear.dist.pathway.dendro$segments), aes(x = x, y = y, xend = xend, yend = yend)) +  # Dendrogram segments
  geom_text(data = (oars.spear.dist.pathway.dendro$labels), aes(x = x, y = -0.45, label = substr(label, 1,6)),  # Labels with color by slope
            hjust = 1, size = 3) +
  geom_tile(data=(oars.spear.dist.pathway.dendro$labels), aes(x = x, y = -0.05, fill = as.factor(RS_Name)),height=0.1,color="white")+
  geom_text(data=(oars.spear.dist.pathway.dendro$labels), aes(x = x,  y = -0.5, label=timing), hjust=0, size=3)+
  scale_fill_manual(values=rs.colors)+
  coord_flip() +  # Flip for vertical dendrogram
  scale_y_reverse(limits=c(1,-0.5), expand = c(0.2, 0)) +  # Reverse y-axis for traditional dendrogram look
  labs(title = "", x = "", y = "") +
  theme_void()+
  facet_wrap(~"Pathway")+
  theme(legend.position="none",
        strip.text=element_text(size=10),
        axis.text.y = element_blank(),  # Remove y-axis labels
        axis.ticks.y = element_blank())
oars.spear.dist.pathway.dendro.plot


# prepare data
oars.spear.dist.cog.dendro = ggdendro::dendro_data(as.dendrogram(hclust(dist(oars.spear.dist.cog, method="euclidean"), method = "ward.D2")), type = "rectangle")
oars.spear.dist.cog.dendro$labels = merge(oars.spear.dist.cog.dendro$labels, oars.spear.dist.meta)
# rescale
oars.spear.dist.cog.dendro.scaling.factor = max(oars.spear.dist.cog.dendro$segments$y, oars.spear.dist.cog.dendro$segments$yend, na.rm = TRUE)
oars.spear.dist.cog.dendro$segments$y <- oars.spear.dist.cog.dendro$segments$y / oars.spear.dist.cog.dendro.scaling.factor  # Normalize to [0, 1]
oars.spear.dist.cog.dendro$segments$yend <- oars.spear.dist.cog.dendro$segments$yend /oars.spear.dist.cog.dendro.scaling.factor  # Normalize endpoints
# plot
oars.spear.dist.cog.dendro.plot = ggplot() +
  geom_segment(data = (oars.spear.dist.cog.dendro$segments), aes(x = x, y = y, xend = xend, yend = yend)) +  # Dendrogram segments
  geom_text(data = (oars.spear.dist.cog.dendro$labels), aes(x = x, y = -0.45, label = substr(label, 1,6)),  # Labels with color by slope
            hjust = 1, size = 3) +
  geom_tile(data=(oars.spear.dist.cog.dendro$labels), aes(x = x, y = -0.05, fill = as.factor(RS_Name)),height=0.1,color="white")+
  geom_text(data=(oars.spear.dist.cog.dendro$labels), aes(x = x,  y = -0.5, label=timing), hjust=0,size=3)+
  scale_fill_manual(values=rs.colors)+
  coord_flip() +  # Flip for vertical dendrogram
  scale_y_reverse(limits=c(1,-0.5), expand = c(0.2, 0)) +  # Reverse y-axis for traditional dendrogram look
  labs(title = "", x = "", y = "") +
  theme_void()+
  facet_wrap(~"COG")+
  theme(legend.position="none",
        strip.text=element_text(size=10),
        axis.text.y = element_blank(),  # Remove y-axis labels
        axis.ticks.y = element_blank())
oars.spear.dist.cog.dendro.plot

# prepare data
oars.spear.dist.cazy.dendro = ggdendro::dendro_data(as.dendrogram(hclust(dist(oars.spear.dist.cazy, method="euclidean"), method = "ward.D2")), type = "rectangle")
oars.spear.dist.cazy.dendro$labels = merge(oars.spear.dist.cazy.dendro$labels, oars.spear.dist.meta)
# rescale
oars.spear.dist.cazy.dendro.scaling.factor = max(oars.spear.dist.cazy.dendro$segments$y, oars.spear.dist.cazy.dendro$segments$yend, na.rm = TRUE)
oars.spear.dist.cazy.dendro$segments$y <- oars.spear.dist.cazy.dendro$segments$y / oars.spear.dist.cazy.dendro.scaling.factor  # Normalize to [0, 1]
oars.spear.dist.cazy.dendro$segments$yend <- oars.spear.dist.cazy.dendro$segments$yend /oars.spear.dist.cazy.dendro.scaling.factor  # Normalize endpoints
# plot
oars.spear.dist.cazy.dendro.plot = ggplot() +
  geom_segment(data = (oars.spear.dist.cazy.dendro$segments), aes(x = x, y = y, xend = xend, yend = yend)) +  # Dendrogram segments
  geom_text(data = (oars.spear.dist.cazy.dendro$labels), aes(x = x, y = -0.45, label = substr(label, 1,6)),  # Labels with color by slope
            hjust = 1, size = 3) +
  geom_tile(data=(oars.spear.dist.cazy.dendro$labels), aes(x = x, y = -0.05, fill = as.factor(RS_Name)),height=0.1,color="white")+
  geom_text(data=(oars.spear.dist.cazy.dendro$labels), aes(x = x, y = -0.5, label=timing), hjust=0, size=3)+
  scale_fill_manual(values=rs.colors)+
  coord_flip() +  # Flip for vertical dendrogram
  scale_y_reverse(limits=c(1,-0.5), expand = c(0.2, 0)) +  # Reverse y-axis for traditional dendrogram look
  labs(title = "", x = "", y = "") +
  theme_void()+
  facet_wrap(~"CAZy")+
  theme(legend.position="none",
        strip.text=element_text(size=10),
        axis.text.y = element_blank(),  # Remove y-axis labels
        axis.ticks.y = element_blank())
oars.spear.dist.cazy.dendro.plot


# prepare data
oars.spear.dist.mbx.dendro = ggdendro::dendro_data(as.dendrogram(hclust(dist(oars.spear.dist.mbx, method="euclidean"), method = "ward.D2")), type = "rectangle")
oars.spear.dist.mbx.dendro$labels = merge(oars.spear.dist.mbx.dendro$labels, oars.spear.dist.meta)
# rescale
oars.spear.dist.mbx.dendro.scaling.factor = max(oars.spear.dist.mbx.dendro$segments$y, oars.spear.dist.mbx.dendro$segments$yend, na.rm = TRUE)
oars.spear.dist.mbx.dendro$segments$y <- oars.spear.dist.mbx.dendro$segments$y / oars.spear.dist.mbx.dendro.scaling.factor  # Normalize to [0, 1]
oars.spear.dist.mbx.dendro$segments$yend <- oars.spear.dist.mbx.dendro$segments$yend /oars.spear.dist.mbx.dendro.scaling.factor  # Normalize endpoints
# plot
oars.spear.dist.mbx.dendro.plot = ggplot() +
  geom_segment(data = (oars.spear.dist.mbx.dendro$segments), aes(x = x, y = y, xend = xend, yend = yend)) +  # Dendrogram segments
  geom_text(data = (oars.spear.dist.mbx.dendro$labels), aes(x = x, y = -0.45, label = substr(label, 1,6)),  # Labels with color by slope
            hjust = 1, size = 3) +
  geom_tile(data=(oars.spear.dist.mbx.dendro$labels), aes(x = x, y = -0.05, fill = as.factor(RS_Name)),height=0.1,color="white")+
  geom_text(data=(oars.spear.dist.mbx.dendro$labels), aes(x = x,  y = -0.5, label=timing), hjust=0,size=3)+
  scale_fill_manual(values=rs.colors)+
  coord_flip() +  # Flip for vertical dendrogram
  scale_y_reverse(limits=c(1,-0.5), expand = c(0.2, 0)) +  # Reverse y-axis for traditional dendrogram look
  labs(title = "", x = "", y = "") +
  theme_void()+
  facet_wrap(~"Metabolite")+
  theme(legend.position="none",
        strip.text=element_text(size=10),
        axis.text.y = element_blank(),  # Remove y-axis labels
        axis.ticks.y = element_blank())
oars.spear.dist.mbx.dendro.plot

oars.spear.dist.asv.dendro.plot|
oars.spear.dist.mgx.dendro.plot|
oars.spear.dist.pathway.dendro.plot|
oars.spear.dist.cog.dendro.plot|
oars.spear.dist.cazy.dendro.plot|
oars.spear.dist.mbx.dendro.plot

# what's the point of this?
# what does it tell us that's different from the PCoA/PCA

# >>> 5. FCAL-OMICS ---------------------------------------------------

# For fecal cal data, we will use non-compliant samples, too
# which means we need to filter data again

# ** use non-compliant, too **

# :: filter ---------------------------------------------------------------

# since we're including non-compliant,
# we need to filter
# so, apply 10% prevalence filter (or 80% for MBX)

# ASV
dim(oars.asv.data.glom) # 1554
oars.asv.data.glom.filt.2.pa = oars.asv.data.glom
oars.asv.data.glom.filt.2.pa[oars.asv.data.glom.filt.2.pa!=0] = 1
oars.asv.data.glom.filt.2 = oars.asv.data.glom[,colSums(oars.asv.data.glom.filt.2.pa) >= nrow(oars.asv.data.glom)*0.1]
dim(oars.asv.data.glom.filt.2) # 340
oars.asv.data.glom.filt.2 = oars.asv.data.glom.filt.2/50000

# MGX
dim(oars.mgx.taxa) # 958
oars.mgx.taxa.filt.2.pa = oars.mgx.taxa
oars.mgx.taxa.filt.2.pa[oars.mgx.taxa.filt.2.pa!=0] = 1
oars.mgx.taxa.filt.2 = oars.mgx.taxa[,colSums(oars.mgx.taxa.filt.2.pa) >= nrow(oars.mgx.taxa)*0.1]
dim(oars.mgx.taxa.filt.2) # 514

# KEGG
dim(oars.mpx.kegg.mat) # 180
oars.mpx.kegg.filt.2.pa = oars.mpx.kegg.mat
oars.mpx.kegg.filt.2.pa[is.na(oars.mpx.kegg.filt.2.pa)] = 0
oars.mpx.kegg.filt.2.pa[oars.mpx.kegg.filt.2.pa!=0] = 1
oars.mpx.kegg.filt.2=oars.mpx.kegg.mat
oars.mpx.kegg.filt.2[is.na(oars.mpx.kegg.filt.2)] = 0
oars.mpx.kegg.filt.2 = oars.mpx.kegg.filt.2[,colSums(oars.mpx.kegg.filt.2.pa) >= nrow(oars.mpx.kegg.mat)*0.1]
dim(oars.mpx.kegg.filt.2) # 178

# COG
dim(oars.mpx.cog.mat) # 2653
oars.mpx.cog.filt.2.pa = oars.mpx.cog.mat
oars.mpx.cog.filt.2.pa[is.na(oars.mpx.cog.filt.2.pa)] = 0
oars.mpx.cog.filt.2.pa[oars.mpx.cog.filt.2.pa!=0] = 1
oars.mpx.cog.filt.2=oars.mpx.cog.mat
oars.mpx.cog.filt.2[is.na(oars.mpx.cog.filt.2)] = 0
oars.mpx.cog.filt.2 = oars.mpx.cog.filt.2[,colSums(oars.mpx.cog.filt.2.pa) >= nrow(oars.mpx.cog.mat)*0.1]
dim(oars.mpx.cog.filt.2) # 2478

# CAZy
dim(oars.mpx.cazy.mat) # 66
oars.mpx.cazy.filt.2.pa = oars.mpx.cazy.mat
oars.mpx.cazy.filt.2.pa[is.na(oars.mpx.cazy.filt.2.pa)] = 0
oars.mpx.cazy.filt.2.pa[oars.mpx.cazy.filt.2.pa!=0] = 1
oars.mpx.cazy.filt.2=oars.mpx.cazy.mat
oars.mpx.cazy.filt.2[is.na(oars.mpx.cazy.filt.2)] = 0
oars.mpx.cazy.filt.2 = oars.mpx.cazy.filt.2[,colSums(oars.mpx.cazy.filt.2.pa) >= nrow(oars.mpx.cazy.mat)*0.1]
dim(oars.mpx.cazy.filt.2) # 65

# MBX
dim(oars.mbx.raw.mat) # 2012
oars.mbx.filt.2.pa = oars.mbx.raw.mat
oars.mbx.filt.2.pa[oars.mbx.filt.2.pa!=0] = 1
oars.mbx.filt.2 = oars.mbx.raw.mat[,colSums(oars.mbx.filt.2.pa) >= nrow(oars.mbx.raw.mat)*0.8]
dim(oars.mbx.filt.2) # 150


# :: process --------------------------------------------------------------


# merge all data sets
oars.fecalcal.omics.data = merge(oars.asv.data.glom.filt.2,
                                 oars.mgx.taxa.filt.2,
                                 by="row.names")
oars.fecalcal.omics.data = merge(oars.fecalcal.omics.data,
                                 oars.mpx.kegg.filt.2%>%as.data.frame()%>%
                                   mutate(`Row.names` = rownames(.)),
                                 by="Row.names")
oars.fecalcal.omics.data = merge(oars.fecalcal.omics.data,
                                 oars.mpx.cog.filt.2%>%as.data.frame()%>%
                                   mutate(`Row.names` = rownames(.)),
                                 by="Row.names")
oars.fecalcal.omics.data = merge(oars.fecalcal.omics.data,
                                 oars.mpx.cazy.filt.2%>%as.data.frame()%>%
                                   mutate(`Row.names` = rownames(.)),
                                 by="Row.names")
oars.fecalcal.omics.data = merge(oars.fecalcal.omics.data,
                                 oars.mbx.filt.2%>%as.data.frame()%>%
                                   mutate(`Row.names` = rownames(.)),
                                 by="Row.names")

oars.fecalcal.omics.data$standard.name = oars.fecalcal.omics.data$`Row.names`
# add fecal cal data
oars.fecalcal.omics.data = merge(oars.fecalcal.omics.data,
                                 metadata.oars.stool[,c("fcal", "standard.name")],
                                 by="standard.name")
rownames(oars.fecalcal.omics.data) = oars.fecalcal.omics.data$Row.names
oars.fecalcal.omics.data$Row.names=NULL
dim(oars.fecalcal.omics.data) 


rownames(oars.fecalcal.omics.data) = oars.fecalcal.omics.data$standard.name
oars.fecalcal.omics.data$standard.name = NULL
# replace NA with 0
oars.fecalcal.omics.data[is.na(oars.fecalcal.omics.data)] = 0

# apply additional filter: remove low prevalent (< 50%)
oars.fecalcal.omics.data.pa = oars.fecalcal.omics.data
oars.fecalcal.omics.data.pa[oars.fecalcal.omics.data.pa!=0] = 1
cols.to.keep = data.frame(keep = colSums(oars.fecalcal.omics.data.pa)>=
                            (nrow(oars.fecalcal.omics.data)*0.50))
oars.fecalcal.omics.data = oars.fecalcal.omics.data[,rownames(subset(cols.to.keep, keep == T))]
dim(oars.fecalcal.omics.data)
# 2712 features
# 59 samples with matching omics data

# convert to lfc (quite slow)
oars.fecalcal.omics.lfc = do.call(cbind, lapply(1:ncol(oars.fecalcal.omics.data), function(col){
  print(col)
  data.subset = data.frame(sample = rownames(oars.fecalcal.omics.data),
                           feature = oars.fecalcal.omics.data[,col])
  data.subset = tidyr::separate(data.subset, col="sample", into=c("HM", "stl", "no"), sep="-", remove=F)
  data.subset = data.subset %>%
    group_by(HM) %>%
    mutate(feature = feature + (min(feature[feature!=0])/2)) %>%
    arrange(no) %>% # ensure stools are ordered by collection time
    mutate(lfc = log2(feature / lag(feature))) # lag takes the preceding value
  new.data = data.subset[,c("sample", "lfc")] %>% as.data.frame()
  colnames(new.data) = c("standard.name", colnames(oars.fecalcal.omics.data)[col])
  rownames(new.data) = new.data$standard.name
  new.data = new.data %>% select(-standard.name)
  return(new.data)
})) %>% as.data.frame()

# save
oars.fecalcal.omics.lfc.saved = oars.fecalcal.omics.lfc
oars.fecalcal.omics.lfc = oars.fecalcal.omics.lfc.saved

# optional: remove non-compliant samples
# note: since we're trying to find associations between
# the microbiome and fecal calprotectin, compliance to RS doesn't matter as much
# oars.fecalcal.omics.lfc = oars.fecalcal.omics.lfc[rownames(oars.fecalcal.omics.lfc) %in% subset(metadata.oars.stool.asv, compliant==T)$standard.name,]

# replace NaN with 0
#oars.fecalcal.omics.lfc[is.na(oars.fecalcal.omics.lfc)] = 0
rownames(oars.fecalcal.omics.lfc) = oars.fecalcal.omics.lfc$standard.name
oars.fecalcal.omics.lfc$standard.name = NULL

# spearman correlation (takes a min)
metadata.oars.stool.omics.cor = Hmisc::rcorr(as.matrix(oars.fecalcal.omics.lfc[,colnames(oars.fecalcal.omics.lfc)!="standard.name"]), type="spearman")

metadata.oars.stool.omics.cor.df = 
  data.frame(reshape2::melt(metadata.oars.stool.omics.cor$r)) %>%
  mutate(pval = reshape2::melt(metadata.oars.stool.omics.cor$P)$value) %>%
  subset(!is.na(pval)) %>%
  arrange(pval)

metadata.oars.stool.omics.cor.bh = subset(metadata.oars.stool.omics.cor.df, Var1 == "fcal") %>%
  mutate(feature = as.character(Var2))%>%
  mutate(padj = p.adjust(pval, method="BH")) %>%
  arrange(padj)

subset(metadata.oars.stool.omics.cor.bh, padj < 0.20)$feature %>% unique() %>% length()

oars.fecalcal.omics.cor.plot = ggplot(metadata.oars.stool.omics.cor.bh,
                                      aes(x=value, y=(padj)))+
  geom_point(shape=21, size=2.5, aes(fill = ifelse(padj < 0.20, "sig", "notsig")))+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, gsub("\\..*", "", gsub("\\(.*", "", gsub("\\/.*", "", gsub("\\,.*", "", feature)))), NA)),
                           size=3)+  
  scale_fill_manual(values=c("black", "white"))+
  geom_hline(yintercept=0.20, linetype=2, alpha=0.5)+
  facet_wrap(~"Fecal Calprotectin Correlation")+
  theme_classic()+theme(legend.position="none",
                        strip.text=element_text(size=10))+
  labs(x="Spearman ρ",
       y="FDR")
oars.fecalcal.omics.cor.plot
# Note: UDP-N-acetyl... is co-annotated as UDP-N-acetylglucosamine-1-phosphate transferase
# suggesting a closer link to Heather's finding (could say, e.g., UDP-N-acetylglucosamine metabolism == IBD)

oars.fecalcal.omics.cor.df = reshape2::melt(oars.fecalcal.omics.lfc[,colnames(oars.fecalcal.omics.lfc)%in%c("fcal", subset(metadata.oars.stool.omics.cor.bh, padj < 0.20)$feature)],
               id.vars="fcal")

oars.fecalcal.omics.cor.df$clean.feature = ifelse(nchar(as.character(oars.fecalcal.omics.cor.df$variable))>30, paste(substr(oars.fecalcal.omics.cor.df$variable, 1, 30), "...", sep=""), as.character(oars.fecalcal.omics.cor.df$variable))

# add omics shape
oars.fecalcal.omics.cor.df = oars.fecalcal.omics.cor.df %>%
  mutate(data.type = ifelse(variable %in% colnames(oars.mgx.taxa.filt.2), "Species",
                            ifelse(variable %in% colnames(oars.mpx.kegg.filt.2), "Pathway", 
                                   ifelse(variable %in% colnames(oars.mpx.cog.filt.2), "COG", 
                                          ifelse(variable %in% colnames(oars.mpx.cazy.filt.2), "CAZy", 
                                                 ifelse(variable %in% colnames(oars.mbx.filt.2), "Metabolite", 
                                                 ifelse(variable %in% colnames(oars.asv.data.glom.filt.2), "ASV", "Fecal calprotectin")))))))
oars.fecalcal.omics.cor.df$omic.type = ifelse(grepl(paste(c("Pathway", "COG", "CAZy"),collapse="|"), oars.fecalcal.omics.cor.df$data.type), "MPX", oars.fecalcal.omics.cor.df$data.type)

oars.fecalcal.omics.cor.df$data.type = factor(oars.fecalcal.omics.cor.df$data.type, levels=c("ASV", "Species", "COG", "Pathway", "CAZy", "Metabolite", "Fecal calprotectin"))
oars.fecalcal.omics.cor.df$omic.type = factor(oars.fecalcal.omics.cor.df$omic.type, levels=c("ASV", "Species", "MPX", "Metabolite", "Fecal calprotectin"))

oars.fecalcal.omics.cor.data = oars.fecalcal.omics.cor.df %>%
  subset(variable %in% slice_min(metadata.oars.stool.omics.cor.bh, order_by=abs(pval), n=8)$feature) %>%
  arrange(order_by=value)  %>%
  subset(!is.na(value)) %>%
  mutate(Var2 = variable)

oars.fecalcal.omics.cor.data.order = merge(oars.fecalcal.omics.cor.data[,c("Var2", "clean.feature")],
                                     metadata.oars.stool.omics.cor.bh[,c("Var2", "value")], by="Var2") %>%
  dplyr::select(clean.feature, value) %>% distinct() %>%
  arrange(value)
oars.fecalcal.omics.cor.data = oars.fecalcal.omics.cor.data %>%
  mutate(clean.feature = factor(clean.feature, levels=oars.fecalcal.omics.cor.data.order$clean.feature))

oars.fecalcal.omics.cor.data = merge(oars.fecalcal.omics.cor.data,
                                     data.frame(omics.colors) %>% mutate(data.type = rownames(.)), by="data.type")

oars.fecalcal.omics.cor.plots = ggplot(oars.fecalcal.omics.cor.data,
                                       aes(x=value,
                                           y=fcal))+
  geom_point(size = 2, aes(shape=data.type, fill=data.type), color="black", alpha=0.9) +
  scale_shape_manual(values=omics.shapes)+
  guides(fill=FALSE, shape=FALSE, size = FALSE)+
  # colors: 
  scale_fill_manual(values=omics.colors)+
  geom_smooth(method="lm", color="black")+
  ggpubr::stat_cor(aes(label = ..r.label..), method="spearman", size=3, 
                   label.y=5.5, #vjust=1.5,
                   label.x.npc=0.5, hjust=0.5)+
  #scale_y_log10()+
  #scale_x_log10()+
  theme_classic()+theme(strip.text=element_text(size=8))+
  facet_wrap(~clean.feature, scales="free", nrow=2)+
  labs(y="Calprotectin Log2FC",
       x="Feature Log2FC")
oars.fecalcal.omics.cor.plots

# draw network of these features // added to Targeted Analysis (Responders)

# 50% prevalence across ALL samples with matching omics (n=58)


# now prepare to make Log2FC heatmap of features most correlated with fecal cal (padj < 0.2) in Resp vs Non-resp

# Subset data to nominally sig taxa (== fcal)
oars.asv.data.glom.filt.2
oars.mgx.taxa.filt.2
oars.mpx.kegg.filt.2
oars.mpx.cog.filt.2
oars.mpx.cazy.filt.2
oars.mbx.filt.2

# features to include
features.to.include = metadata.oars.stool.omics.cor.bh %>%
  subset(padj < 0.20) %>%
  arrange(abs(value)) # do not filter beyond this; filter later
features.to.include = features.to.include$feature

# ASV Taxa
oars.asv.lfc = oars.asv.data.glom.filt.2 %>% reshape2::melt() %>%
  mutate(value = value+min(value[value!=0])/2)%>%
  rename("standard.name" = "Var1", "variable" = "Var2") %>%
  merge(metadata.oars.stool.double[,c("standard.name", "phase", "HM", "reltiming", "timing", "response")])%>%
  group_by(HM, phase, variable) %>%
  mutate(lfc = log2(value) - lag(log2(value))) %>%
  subset(reltiming == "post") %>%
  mutate(sample = paste(HM, timing, sep="_")) %>%
  subset(variable %in% features.to.include)

# MGX Taxa
oars.mgx.lfc = oars.mgx.taxa.filt.2 %>% reshape2::melt() %>%
  mutate(value = value+min(value[value!=0])/2)%>%
  rename("standard.name" = "Var1", "variable" = "Var2") %>%
  merge(metadata.oars.stool.double[,c("standard.name", "phase", "HM", "reltiming", "timing", "response")])%>%
  group_by(HM, phase, variable) %>%
  mutate(lfc = log2(value) - lag(log2(value))) %>%
  subset(reltiming == "post") %>%
  mutate(sample = paste(HM, timing, sep="_")) %>%
  #subset(variable %in% (subset(oars.asv.data.glom.lmer, pval < 0.05) %>%
  subset(variable %in% features.to.include)

# MPX KEGG
oars.mpx.kegg.lfc = oars.mpx.kegg.filt.2 %>% reshape2::melt() %>%
  mutate(value = value+min(value[value!=0])/2)%>%
  rename("standard.name" = "Var1", "variable" = "Var2") %>%
  merge(metadata.oars.stool.double[,c("standard.name", "phase", "HM", "reltiming", "timing", "response")])%>%
  group_by(HM, phase, variable) %>%
  mutate(lfc = log2(value) - lag(log2(value))) %>%
  mutate(lfc = ifelse(is.na(lfc), 0, lfc))%>%
  subset(reltiming == "post") %>%
  mutate(sample = paste(HM, timing, sep="_"))  %>%
  #subset(variable %in% (subset(oars.asv.data.glom.lmer, pval < 0.05) %>%
  subset(variable %in% (features.to.include))

# MPX COG
oars.mpx.cog.lfc = oars.mpx.cog.filt.2 %>% reshape2::melt() %>%
  mutate(value = value+min(value[value!=0])/2)%>%
  rename("standard.name" = "Var1", "variable" = "Var2") %>%
  merge(metadata.oars.stool.double[,c("standard.name", "phase", "HM", "reltiming", "timing", "response")])%>%
  group_by(HM, phase, variable) %>%
  mutate(lfc = log2(value) - lag(log2(value))) %>%
  mutate(lfc = ifelse(is.na(lfc), 0, lfc))%>%
  subset(reltiming == "post") %>%
  mutate(sample = paste(HM, timing, sep="_"))  %>%
  #subset(variable %in% (subset(oars.asv.data.glom.lmer, pval < 0.05) %>%
  subset(variable %in% (features.to.include))

# MPX CAZy
oars.mpx.cazy.lfc = oars.mpx.cazy.filt.2 %>% reshape2::melt() %>%
  mutate(value = value+min(value[value!=0])/2)%>%
  rename("standard.name" = "Var1", "variable" = "Var2") %>%
  merge(metadata.oars.stool.double[,c("standard.name", "phase", "HM", "reltiming", "timing", "response")])%>%
  group_by(HM, phase, variable) %>%
  mutate(lfc = log2(value) - lag(log2(value))) %>%
  mutate(lfc = ifelse(is.na(lfc), 0, lfc))%>%
  subset(reltiming == "post") %>%
  mutate(sample = paste(HM, timing, sep="_"))  %>%
  #subset(variable %in% (subset(oars.asv.data.glom.lmer, pval < 0.05) %>%
  subset(variable %in% (features.to.include))

# MBX
oars.mbx.lfc = oars.mbx.filt.2 %>% as.matrix() %>% reshape2::melt() %>%
  mutate(value = value+min(value[value!=0])/2)%>%
  rename("standard.name" = "Var1", "variable" = "Var2") %>%
  merge(metadata.oars.stool.double[,c("standard.name", "phase", "HM", "reltiming", "timing", "response")])%>%
  group_by(HM, phase, variable) %>%
  mutate(lfc = log2(value) - lag(log2(value))) %>%
  mutate(lfc = ifelse(is.na(lfc), 0, lfc))%>%
  subset(reltiming == "post") %>%
  mutate(sample = paste(HM, timing, sep="_"))  %>%
  subset(variable %in% (features.to.include))

# merge all omics lfc data
oars.asv.mgx.mpx.mbx.lfc = rbind(
  oars.asv.lfc,
  oars.mgx.lfc,
  oars.mpx.kegg.lfc,
  oars.mpx.cog.lfc,
  oars.mpx.cazy.lfc,
  oars.mbx.lfc) %>% data.frame() 

# keep samples with matching omics
oars.overlap.samples = Reduce(intersect, list(oars.asv.lfc$sample,
                                 oars.mgx.lfc$sample,
                                 oars.mpx.kegg.lfc$sample,
                                 oars.mbx.lfc$sample))
oars.asv.mgx.mpx.mbx.lfc = subset(oars.asv.mgx.mpx.mbx.lfc, sample %in% oars.overlap.samples)

# create mapping for patients
oars.omics.lfc.map = oars.mbx.lfc[,c("sample", "response")] %>% data.frame() %>%distinct()
rownames(oars.omics.lfc.map) = oars.omics.lfc.map$sample
oars.omics.lfc.map$sample = NULL

oars.omics.lfc.map = data.frame(sample = oars.overlap.samples,
                                response = oars.omics.lfc.map[rownames(oars.omics.lfc.map) %in% oars.overlap.samples,])
# n = 18

# scale per variable (across samples)
oars.asv.mgx.mpx.mbx.lfc = oars.asv.mgx.mpx.mbx.lfc %>% 
  group_by(variable) %>% 
  subset(!is.na(lfc))%>%
  mutate(slfc = scale(lfc))
# add feature map
oars.asv.mgx.mpx.mbx.lfc.map = data.frame(variable = oars.asv.mgx.mpx.mbx.lfc[,c("variable")]) %>% distinct() %>%
  mutate(datatype = ifelse(variable %in% oars.mgx.lfc$variable, "Species",
                                  ifelse(variable %in% oars.mpx.kegg.lfc$variable, "Pathway", 
                                         ifelse(variable %in% oars.mpx.cog.lfc$variable, "COG", 
                                                ifelse(variable %in% oars.mpx.cazy.lfc$variable, "CAZy",
                                                       ifelse(variable %in% oars.mbx.lfc$variable, "Metabolite",
                                                              "ASV")))))) %>% data.frame()
rownames(oars.asv.mgx.mpx.mbx.lfc.map) = oars.asv.mgx.mpx.mbx.lfc.map$variable
colnames(oars.asv.mgx.mpx.mbx.lfc.map)[1] = "feature"
nrow(oars.asv.mgx.mpx.mbx.lfc.map)

rownames(oars.asv.mgx.mpx.mbx.lfc.map) = oars.asv.mgx.mpx.mbx.lfc.map$feature
oars.asv.mgx.mpx.mbx.lfc.map$feature = NULL
colnames(oars.asv.mgx.mpx.mbx.lfc.map)[1] = "Data type"
dim(oars.asv.mgx.mpx.mbx.lfc.map)

# refresh response names
oars.omics.lfc.map$response = ifelse(oars.omics.lfc.map$response == "low", "Weak Response", "Strong Response")
colnames(oars.omics.lfc.map) = c("sample", "Response")

# refresh variable names
#oars.asv.mgx.mpx.mbx.lfc.feature.map = data.frame(variable = make.names(features.to.include),
#                                              clean.name = features.to.include)
#oars.asv.mgx.mpx.mbx.lfc$feature = oars.asv.mgx.mpx.mbx.lfc.feature.map[match(oars.asv.mgx.mpx.mbx.lfc$variable,
#                                                                        oars.asv.mgx.mpx.mbx.lfc.feature.map$variable),]$clean.name
#rownames(oars.asv.mgx.mpx.mbx.lfc.map) = oars.asv.mgx.mpx.mbx.lfc.feature.map[match(rownames(oars.asv.mgx.mpx.mbx.lfc.map),
#                                                                            oars.asv.mgx.mpx.mbx.lfc.feature.map$variable),]$clean.name
oars.asv.mgx.mpx.mbx.lfc.map$`Fcal Correlation` = metadata.oars.stool.omics.cor.bh[match(rownames(oars.asv.mgx.mpx.mbx.lfc.map),
                                                                              metadata.oars.stool.omics.cor.bh$Var2),]$value
# shorten
rownames(oars.asv.mgx.mpx.mbx.lfc.map) = ifelse(nchar(as.character(rownames(oars.asv.mgx.mpx.mbx.lfc.map)))>40, paste(substr(rownames(oars.asv.mgx.mpx.mbx.lfc.map), 1, 40), "...", sep=""), as.character(rownames(oars.asv.mgx.mpx.mbx.lfc.map)))
oars.asv.mgx.mpx.mbx.lfc$feature = oars.asv.mgx.mpx.mbx.lfc$variable
oars.asv.mgx.mpx.mbx.lfc$feature = ifelse(nchar(as.character(oars.asv.mgx.mpx.mbx.lfc$feature))>40, paste(substr(oars.asv.mgx.mpx.mbx.lfc$feature, 1, 40), "...", sep=""), as.character(oars.asv.mgx.mpx.mbx.lfc$feature))

# clean
oars.asv.mgx.mpx.mbx.lfc.map$feature = NULL


# :: PLS-DA -----------------------------------------------------------------


# Step 1: Prepare data
# Assume lfc_matrix is your data (rows = samples, columns = features)
# Example: lfc_matrix <- metadata.oars.stool.omics.cor[, numeric_cols]
# labels is a factor vector of "strong" vs "weak"
# Replace with your actual data

lfc_matrix <- reshape2::acast(oars.asv.mgx.mpx.mbx.lfc,
                              sample ~ feature, value.var="lfc") %>% as.matrix()

if(all(oars.omics.lfc.map$sample == rownames(lfc_matrix))){
labels <- as.factor(oars.omics.lfc.map$Response)  #
}

# Check for NAs and impute if needed
lfc_matrix[is.na(lfc_matrix)] <- 0  # Simple imputation; consider mixOmics::nipals for better handling

# Step 2: Run PLS-DA
plsda_model <- mixOmics::plsda(X = lfc_matrix, 
                               Y = labels, 
                               ncomp = 10)  # Try 5 components

# Tune number of components using cross-validation
set.seed(25)
perf_plsda <- mixOmics::perf(plsda_model, 
                             #validation = "Mfold", 
                             nrepeat = 10,
                             folds = 5, progressBar = TRUE)

ggplot(perf_plsda$error.rate.class$mahalanobis.dist%>%
         data.frame() %>%
         summarize(total.err = colSums(.)) %>%
         mutate(ncomp = seq(1:10)))+
  aes(x=ncomp, y=total.err*100)+
  scale_x_continuous(breaks=seq(1:10))+
  geom_line(color="red", linewidth=1)+
  theme_classic()+theme(strip.text=element_text(size=10))+
  facet_wrap(~"Method: Mahalanobis Distance")+
  labs(x="Number of Components",
       y="Total Error (%)")

ncomp_opt = 5

# Refit with optimal components
plsda_model <- mixOmics::plsda(X = lfc_matrix, Y = labels, ncomp = ncomp_opt)

# Step 4: Extract discriminative features (VIP scores)
vip_scores <- mixOmics::vip(plsda_model)  # VIP for each component
vip_df <- data.frame(
  Feature = colnames(lfc_matrix),
  VIP = vip_scores[, ncomp_opt]  # Use last component
) %>%
  arrange(desc(VIP))

ggplot(vip_df,
       aes(x=reorder(Feature, VIP), y=VIP))+
  coord_flip()+
  geom_point(shape=21, aes(fill=scale(VIP)))+
  scale_fill_gradient2(low="blue", high="red")+
  theme_classic()+theme(legend.position="none")

# Filter top discriminative features (e.g., VIP > 1)
top_features <- vip_df %>%
  arrange(-VIP) %>%
 # filter(VIP > 1) %>%
  head(10)  # Top 20 features
nrow(top_features)

# Step 5: Visualize
# Score plot (sample separation)
mixOmics::plotIndiv(plsda_model, comp = c(1, 2), group = labels, legend = TRUE,
          title = "PLS-DA Score Plot (Strong vs Weak)", ellipse = TRUE)
# plot manually

rownames(oars.omics.lfc.map) = oars.omics.lfc.map$sample
plsda_model_plot = plsda_model$variates$X[,c(1,2)] %>% data.frame() %>%
  merge(oars.omics.lfc.map, by="row.names") %>%
  ggplot(aes(x=comp1, y=comp2))+
  stat_ellipse(aes(color=Response), linewidth=1.5, alpha=0.5)+
  scale_color_discrete(guide = "none")+
  geom_point(aes(fill=Response), shape=21, size=3)+
  theme_classic()+theme(strip.text=element_text(size=10),
                        legend.position=c(0.8, 0.15),
                        legend.title = element_blank(),
                        legend.background = element_rect(color="black"))+
  facet_wrap(~"PLS-DA")+
  labs(x=paste("Comp 1: ", round(plsda_model$prop_expl_var$X[1]*100, digits=1), "%", sep=""),
       y=paste("Comp 2: ", round(plsda_model$prop_expl_var$X[2]*100, digits=1), "%", sep=""))
plsda_model_plot

# select top 20 features and add VIMP score
top_features # done
nrow(top_features)

# :: Heatmap --------------------------------------------------------------

oars.asv.mgx.mpx.mbx.lfc.map.vimp = merge(oars.asv.mgx.mpx.mbx.lfc.map,
                                      top_features, by="row.names")
rownames(oars.asv.mgx.mpx.mbx.lfc.map.vimp) = oars.asv.mgx.mpx.mbx.lfc.map.vimp$Row.names
oars.asv.mgx.mpx.mbx.lfc.map.vimp$Row.names = NULL
oars.asv.mgx.mpx.mbx.lfc.map.vimp$Feature = NULL
colnames(oars.asv.mgx.mpx.mbx.lfc.map.vimp)[1] = "Data type"

pheatmap::pheatmap(reshape2::acast(subset(oars.asv.mgx.mpx.mbx.lfc, feature %in% rownames(oars.asv.mgx.mpx.mbx.lfc.map.vimp)),
                                   feature ~ sample, value.var="slfc"),
                   color=colorRampPalette(c("blue","white", "red"))(100),
                   #angle_col = 45,
                   clustering_distance_rows = "correlation",
                   clustering_distance_cols = "correlation",
                   fontsize_row = 8,
                   fontsize_col = 8,
                   annotation_col = oars.omics.lfc.map %>% dplyr::select(-sample),
                   annotation_colors = list(Response = c(`Strong Response` = gg_color_hue(2)[1],
                                                         `Weak Response` = gg_color_hue(2)[2]),
                                            #`Coefficient` = colorRampPalette(c("blue", "red"))(2),
                                            #Coefficient = c(`Positive` = "red",
                                            #                `Negative` = "blue"),
                                            `Fcal Correlation` = colorRampPalette(c("blue","white", "red"))(100),
                                            `VIP` = colorRampPalette(c("white", "red"))(100),
                                            `Data type` = c(ASV = "#8DD3C7",
                                                            Species = "#FFFFB3",
                                                            Pathway = "#BEBADA",
                                                            COG = "#FB8072",
                                                            CAZy = "#80B1D3",
                                                            Metabolite = "#FDB462")),
                   annotation_row = oars.asv.mgx.mpx.mbx.lfc.map.vimp,
                   breaks=c(seq(min(na.omit(oars.asv.mgx.mpx.mbx.lfc$slfc)), 0, length.out=ceiling(100/2) + 1), 
                            seq(max(na.omit(oars.asv.mgx.mpx.mbx.lfc$slfc))/100, max(na.omit(oars.asv.mgx.mpx.mbx.lfc$slfc)), length.out=floor(100/2))))
# correlation clustering
# scaled Log2FC

# :: Network --------------------------------------------------------------

# these are the features sig associated with fcal (padj < 0.20)
features.cor.fcal.df = subset(metadata.oars.stool.omics.cor.bh, padj < 0.20)
# features
features.cor.fcal = c(features.cor.fcal.df$Var1, features.cor.fcal.df$Var2) %>% unique()

# all feature cors
all.features.cor.all = metadata.oars.stool.omics.cor.df %>%
  mutate(padj = p.adjust(pval, method="BH")) 
# apply filter for plot
all.features.cor = all.features.cor.all %>%
  subset(padj < 0.05) %>%
  #subset(abs(value) > 0.65 | Var1 == "fcal" | Var2 == "fcal") %>%
  subset(abs(value) > 0.5) %>%
  subset(Var1 %in% features.cor.fcal & Var2 %in% features.cor.fcal)
c(unique(all.features.cor$Var1),unique(all.features.cor$Var2))
# Calculate adjusted P values of features sig associated with fecal calprotectin

# build with UMAP instead (regularized spearman)
set.seed(25)
adj_matrix.umap = #metadata.oars.stool.omics.cor.sig %>% # or regularized with adj_matrix.df
  all.features.cor %>%
  # subset to strong
  #subset(abs(value) > 0.5) %>%
  # cycle through melt/cast to impute NA as 0
  reshape2::acast(Var1 ~ Var2, value.var="value") %>%
  reshape2::melt() %>%
  mutate(value = ifelse(is.na(value), 0, value)) %>% 
  #mutate(value = 1-value) %>% 
  reshape2::acast(Var1 ~ Var2, value.var="value") %>% as.matrix() 
adj_matrix.umap[is.na(adj_matrix.umap)]=0

# UMAP is a bit strong
# adj_matrix.umap = umap::umap(t(adj_matrix.umap))
# adj_matrix.umap.df = adj_matrix.umap$layout %>% as.data.frame()   

# use PCoA
adj_matrix.umap.df = ape::pcoa(as.dist(1-(adj_matrix.umap)))$vectors[,c(1:2)] %>% as.data.frame()
colnames(adj_matrix.umap.df) = c("V1", "V2")

# extract coords
#rownames(adj_matrix.umap.df) = colnames(oars.fecalcal.omics.data[,colnames(oars.fecalcal.omics.data) %in% c("fcal", features.to.include)])
adj_matrix.umap.df$feature = rownames(adj_matrix.umap.df)

# add edge data
adj_matrix.umap.df = do.call(rbind, lapply(1:nrow(adj_matrix.umap.df),function(x){
  data.subset = adj_matrix.umap.df[x,]
  # which features is this feature correlated to?
  cors = subset(all.features.cor, Var1 == data.subset$feature)
  colnames(cors)[1:3] = c("feature", "Var2", "cor")
  data.subset = merge(data.subset, cors, by="feature")
  # add coords
  cors.coords = adj_matrix.umap.df
  colnames(cors.coords) = c("V1.B", "V2.B", "Var2")
  data.subset = merge(data.subset,
                      cors.coords, by="Var2")
  # good
  return(data.subset)
}))


adj_matrix.umap.df$feature = gsub("fcal", "Fecal calprotectin", adj_matrix.umap.df$feature)
adj_matrix.umap.df = adj_matrix.umap.df %>%
  mutate(data.type = ifelse(feature %in% colnames(oars.mgx.taxa.filt.2), "Species",
                            ifelse(feature %in% colnames(oars.mpx.kegg.filt.2), "Pathway", 
                                   ifelse(feature %in% colnames(oars.mpx.cog.filt.2), "COG", 
                                          ifelse(feature %in% colnames(oars.mpx.cazy.filt.2), "CAZy", 
                                                 ifelse(feature %in% colnames(oars.mbx.filt.2), "Metabolite",
                                                        ifelse(feature %in% colnames(oars.asv.data.glom.filt.2), "ASV", "Fecal calprotectin")))))))
adj_matrix.umap.df$omic.type = ifelse(grepl(paste(c("COG", "Pathway", "CAZy"), collapse="|"), adj_matrix.umap.df$data.type), "MPX", adj_matrix.umap.df$data.type)

adj_matrix.umap.df$data.type = factor(adj_matrix.umap.df$data.type, levels=c("Fecal calprotectin", "ASV", "Species", "Pathway", "COG", "CAZy", "Metabolite")) 

adj_matrix.umap.df$short = ifelse(nchar(as.character(adj_matrix.umap.df$feature))>40, paste(substr(adj_matrix.umap.df$feature, 1, 40), "...", sep=""), as.character(adj_matrix.umap.df$feature))

# only label VIP
adj_matrix.umap.df$vip = ifelse(adj_matrix.umap.df$short %in% c(top_features$Feature, "Fecal calprotectin"), adj_matrix.umap.df$short, NA)

metadata.oars.stool.omics.cor.network.pcoa = ggplot(adj_matrix.umap.df, 
                                                    aes(V1, V2))+
  geom_segment(aes(x=V1, xend=V1.B, y=V2, yend=V2.B, color = (cor), 
                   size = abs(cor)), alpha = 0.1) +  # Edge color by correlation, size by |corr|
  geom_point(size = 4, aes(shape=data.type), fill="white", color = "white", alpha=1) +
  geom_point(size = 4, aes(shape=data.type, fill=data.type), color = "black", alpha=0.9) +
  scale_shape_manual(values=c("Fecal calprotectin" = 21, omics.shapes))+
  guides(shape = FALSE, fill=FALSE, size = FALSE)+
  scale_fill_manual(values= c("Fecal calprotectin" = "white", omics.colors))+
  ggnetwork::geom_nodetext_repel(aes(label = short), 
                      size = 3) +  # Node labels
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +  # Color for pos/neg correlations
  scale_size(range = c(0.5, 2)) +  # Edge thickness
  ggnetwork::theme_blank() +  # ggnetwork's blank theme
 # theme(legend.position="none")+
  labs(color = "Correlation",
       fill="Data type")
metadata.oars.stool.omics.cor.network.pcoa
# maybe keep as LARGE supplemental



# :: Response Plot --------------------------------------------------------

# Goal: plot B. adol, F.cal, CAZyme change, per individual, over time

# need: HM, time, response, B. adol, fcal, starch CAZyme


oars.badol.starch.fcal.plot.data = metadata.oars.stool.asv[,c("HM","timing", "phase","compliant","oars.timing",
                              "fcal", "starch", "oars.days", "rs.col")] %>%
  mutate(standard.name = rownames(.)) %>%
  merge(data.frame(standard.name = names(oars.mgx.taxa[,"Bifidobacterium_adolescentis"]),
                   badol = oars.mgx.taxa[,"Bifidobacterium_adolescentis"]),
        by="standard.name") %>%
  subset(!is.na(fcal)) %>% 
  subset(HM != "HM0618") %>%
  group_by(HM) %>%
  mutate(scaled.badol = scale(log10(badol+0.5e-05)),
         scaled.starch = scale(log10(na.omit(starch))),
         scaled.fcal = scale(log10(fcal))) %>%
  merge(metadata.oars.stool.asv.double[,c("HM","timing","reltiming", "response")] %>%
          subset(reltiming == "pre"), by=c("HM", "timing"), all.x = T) %>%
  arrange(HM, oars.days) %>%
  subset(compliant == T) %>%
  # add max days and response
  group_by(HM) %>%
  mutate(max.days.1 = ifelse(!is.na(max(oars.days[timing=="3M"])), max(oars.days[timing=="3M"]), NA),
         response.1 = ifelse(!is.na(response[reltiming == "pre" & timing=="0M"]),response[reltiming == "pre" & timing=="0M"], NA),
         max.days.2 = ifelse(!HM %in% c("HM0759", "HM0924", "HM0932"), max(oars.days[timing=="6M"]), NA),
         response.2 = ifelse(!HM %in% c("HM0759", "HM0924", "HM0932"),response[reltiming == "pre" & timing=="3M"], NA)) %>%
  # and now make NA for repeated
  mutate(response.1 = ifelse(timing == "0M", response.1, NA),
         response.2 = ifelse(timing == "3M", response.2, NA)) %>%
  # reorder HM manually
  mutate(HM = factor(HM, levels=c("HM0883", "HM0906","HM0903", "HM0874", "HM0819", "HM0618",
                                  "HM0932", "HM0902", "HM0899", "HM0844", "HM0924", "HM0759"))) %>%
  as.data.frame() 
  #dplyr::select(HM, timing, response, max.days.1, response.1, max.days.2, response.2)

# plot

oars.badol.starch.fcal.plot = ggplot(data=oars.badol.starch.fcal.plot.data,
       aes(x=oars.days))+
  # add rectangles to indicate phase + response
  #annotate("rect", xmin=0, 
  #         xmax=max(subset(metadata.oars.stool.asv, oars.on.rs == "onRS")$oars.days), 
  #         ymin=-Inf, ymax=Inf, alpha=0.2) +
  geom_rect(data=subset(oars.badol.starch.fcal.plot.data, timing == "0M"), 
            aes(xmin=0, xmax=max.days.1, 
                    fill=(response.1)),na.rm = TRUE,
           ymin=-Inf, ymax=Inf, alpha=0.5, color="white") +
  geom_rect(data=subset(oars.badol.starch.fcal.plot.data, timing == "3M"), 
            aes(xmin=max.days.1, xmax=max.days.2, 
                fill=(response.2)),na.rm = TRUE,
            ymin=-Inf, ymax=Inf, alpha=0.5, color="white") +
  # add lines for variables
  geom_line(aes(y=scaled.badol), linetype=1, linewidth=0.5, color="black")+
  geom_line(aes(y=scaled.starch), linetype=2, linewidth=0.5, color="black")+
  geom_line(aes(y=scaled.fcal), linetype=3, linewidth=0.75, color="red")+
  scale_fill_manual(values=c("salmon", "darkturquoise"), na.value="white")+
  theme_classic()+
  theme(strip.text.y.right = element_text(angle = 0),
        strip.background = element_blank(),
        legend.position="none")+
  labs(x="Days since starting RS",
       y="Scaled Log Abundance")+
  facet_grid(HM~.)
oars.badol.starch.fcal.plot

# cool


# :: B. adol --------------------------------------------------------------

oars.badol.starch.fcal.plot.data

#oars.badol.simple.data = oars.asv.data.glom.prep[,colnames(oars.asv.data.glom.prep) %in% c("g__Bifidobacterium_388775_s__faecale", "HM", "timing", "phase", "reltiming", "response","diagnosis","adj.fiber","standard.name")]

oars.badol.simple.data = oars.mgx.prep[,colnames(oars.mgx.prep) %in% c("Bifidobacterium_adolescentis", "HM", "timing", "phase", "reltiming", "response","diagnosis","adj.fiber","standard.name")]

# note: keep only paired datapoints 
stats.badol.ph.low = lmerTest::lmer(scale(log10(Bifidobacterium_adolescentis+0.000005)) ~ reltiming + phase + diagnosis + adj.fiber + (1|HM), 
                                   subset(oars.badol.simple.data, response =="low")) %>% summary() %>% coef()
stats.badol.ph.high = lmerTest::lmer(scale(log10(Bifidobacterium_adolescentis+0.000005)) ~ reltiming + phase + diagnosis + adj.fiber + (1|HM), 
                                    subset(oars.badol.simple.data, response == "high")) %>% summary() %>% coef()
stats.badol.ph.interact = lmerTest::lmer(scale(log10(Bifidobacterium_adolescentis+0.000005)) ~ reltiming*response + phase + diagnosis+adj.fiber + (1|HM),
                                          oars.badol.simple.data) %>% summary() %>% coef()

# make data.frame
stats.badol.ph = data.frame(Response = c("Strong Response", "Weak Response"),
                           pval = c(stats.badol.ph.high[2,5], stats.badol.ph.low[2,5]))

oars.badol.simple.data$hm_phase = paste(oars.badol.simple.data$HM, oars.badol.simple.data$phase, sep="_")

oars.badol.simple.plot = ggplot(subset(oars.badol.simple.data, !hm_phase %in% c("HM0924_rs2", "HM0932_rs2", "HM0759_rs2")) %>%
                                          mutate(Response = ifelse(response == "low", "Weak Response", "Strong Response")) %>%
                                  merge(metadata.oars.stool.double[,c("hm_phase","reltiming", "RS_Name")],
                                        by=c("hm_phase", "reltiming")) %>% distinct() %>%
                                  mutate(reltiming = factor(ifelse(reltiming == "pre", "Pre", "Post"), levels=c("Pre", "Post"))),
                                        aes(x=reltiming, y=Bifidobacterium_adolescentis))+
  scale_y_log10(limits = c(0.0000005, 100))+
  geom_boxplot(width=0.5)+
  geom_line(aes(group=hm_phase), linetype=2, alpha=0.5, linewidth=0.5)+
  geom_point(shape=21, aes(fill=RS_Name), size=3)+
  #geom_hline(yintercept=250, linetype=1, alpha=1, color="red")+
 # ggrepel::geom_label_repel(aes(label=RS_Name),size=2)+
  scale_fill_manual(values=rs.colors, na.value = "lightgrey")+
  geom_text(data = stats.badol.ph %>%
              mutate(Response = ifelse(Response == "Weak Response", "Weak Response", "Strong Response")), 
            x=1.5, y=Inf, vjust=1.3,
            aes(label=paste("p:", round(pval, digits=3))), size=3)+
  facet_wrap(~Response)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=10),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="",
       y="B. adolescentis (%)")
oars.badol.simple.plot
oars.badol.starch.fcal.plot
# Compelling (but not supported by ASV data)


# >>> 6. MACHINE LEARNING -----------------------------------------------------------

# goal: see if features predict RS fermentation response

# :: process ---------------------------------------------------------

# prepare data; use omic.prep function, but use 20% prevalence filter

# ASV
oars.asv.data.glom.prep.ml = delta.omic.prepare(oars.asv.data.glom,  # rarefied reads (50,000)
                                             normalize=T, # normalize to 100%
                                             min.abun = 1, # 1 read
                                             for_ml = T,
                                             prev=0.1) #
ncol(oars.asv.data.glom.prep.ml)-8 # 406

# MGX
oars.mgx.prep.ml = delta.omic.prepare(oars.mgx.taxa, # sum to ~100% data (~ because of median)
                                   normalize=T, # re-normalize to 100%
                                   min.abun = 0.01, # 0.01% min abun for prev
                                   for_ml=T,
                                   prev = 0.1) #
ncol(oars.mgx.prep.ml)-8 # 480

# KEGG
oars.mpx.kegg.prep.ml = delta.omic.prepare(oars.mpx.kegg.mat, # norm intensity
                                        normalize=F, # don't normalize
                                        min.abun = 1, # detected at all
                                        for_ml=T,
                                        prev=0.1) # 10% prev per group
ncol(oars.mpx.kegg.prep.ml)-8 # 176

# COG
oars.mpx.cog.prep.ml = delta.omic.prepare(oars.mpx.cog.mat, 
                                       normalize=F, # don't normalize
                                       min.abun = 1,# detected
                                       for_ml=T,
                                       prev=0.1)# 10% prev per group
ncol(oars.mpx.cog.prep.ml)-8 # 2483

# CAZy
oars.mpx.cazy.prep.ml = delta.omic.prepare(oars.mpx.cazy.mat, 
                                        normalize=F, # don't normalize
                                        min.abun = 1,# detected
                                        for_ml=T,
                                        prev=0.1)# 10% prev per group
ncol(oars.mpx.cazy.prep.ml)-8 # 66

# MBX
oars.mbx.prep.ml = delta.omic.prepare(oars.mbx.annotated.mat, 
                                   normalize=F, # don't normalize
                                   min.abun = 1, # detected
                                   pseudo=T,
                                   for_ml=T,
                                   prev=0.8) # 80% prev per group
ncol(oars.mbx.prep.ml)-8 # 165

# FFQ
redcap.data.ffq = readRDS("./2025_08_16_oars_ffq_table.Rds")
redcap.data.ffq = redcap.data.ffq - 1
redcap.data.ffq.prep.ml = delta.omic.prepare(data = redcap.data.ffq,
                                             prev = 0.1,
                                             min.abun = 1,
                                             pseudo = F,
                                             for_ml = T,
                                             normalize = F)
ncol(redcap.data.ffq.prep.ml)-8 # 233



# ASV
oars.asv.baseline = oars.asv.data.glom.prep.ml[,!colnames(oars.asv.data.glom.prep.ml) %in% c("timing","reltiming","standard.name")] %>%
  mutate(response = as.factor(response))
# MGX
oars.mgx.baseline = oars.mgx.prep.ml[,!colnames(oars.mgx.prep.ml) %in% c("timing","reltiming","standard.name")] %>%
  mutate(response = as.factor(response))
# MPX (KEGG)
oars.mpx.kegg.baseline = oars.mpx.kegg.prep.ml[,!colnames(oars.mpx.kegg.prep.ml) %in% c("timing","reltiming","standard.name")] %>%
  mutate(response = as.factor(response))
# MPX (COG)
oars.mpx.cog.baseline = oars.mpx.cog.prep.ml[,!colnames(oars.mpx.cog.prep.ml) %in% c("timing","reltiming","standard.name")] %>%
  mutate(response = as.factor(response))
# MPX (CAZy)
oars.mpx.cazy.baseline = oars.mpx.cazy.prep.ml[,!colnames(oars.mpx.cazy.prep.ml) %in% c("timing","reltiming","standard.name")] %>%
  mutate(response = as.factor(response))
# MBX
oars.mbx.baseline = oars.mbx.prep.ml[,!colnames(oars.mbx.prep.ml) %in% c("timing","reltiming","standard.name")] %>%
  mutate(response = as.factor(response))
# FFQ
oars.ffq.baseline = redcap.data.ffq.prep.ml[,!colnames(redcap.data.ffq.prep.ml) %in% c("timing","reltiming","standard.name")] %>%
  mutate(response = as.factor(response))

## merge datasets for RF

# All Omics
# add ASVs (n=21 samples, 12 HM)
# to MGX (n=21 samples, 12 HM)
oars.all.omics.baseline = merge(oars.asv.baseline,
                                oars.mgx.baseline, by=c("HM", "phase", "response","diagnosis","adj.fiber"))
# add COGs (n=19 samples, 11 HM)
oars.all.omics.baseline = merge(oars.all.omics.baseline,
                                oars.mpx.cog.baseline,  by=c("HM", "phase", "response","diagnosis","adj.fiber"))
# add KEGG (n=19 samples, 11 HM)
oars.all.omics.baseline = merge(oars.all.omics.baseline,
                                oars.mpx.kegg.baseline,  by=c("HM", "phase", "response","diagnosis","adj.fiber"))
# add CAZy (n=19 samples, 11 HM)
oars.all.omics.baseline = merge(oars.all.omics.baseline,
                                oars.mpx.cazy.baseline,  by=c("HM", "phase", "response","diagnosis","adj.fiber"))
# add MBX (n=21 samples, 12 HM)
oars.all.omics.baseline = merge(oars.all.omics.baseline,
                                oars.mbx.baseline,  by=c("HM", "phase", "response","diagnosis","adj.fiber"))
# add FFQ (n=21 samples, 12 HM)
oars.all.omics.baseline = merge(oars.all.omics.baseline,
                                oars.ffq.baseline,  by=c("HM", "phase", "response","diagnosis","adj.fiber"))

# add rownames
rownames(oars.all.omics.baseline) = paste(oars.all.omics.baseline$HM, 
                                          oars.all.omics.baseline$phase, sep="_")

dim(oars.all.omics.baseline) # 19 samples and 4014 features

# data is ready for random forest

# Key question: Should we subset the datasets to only include matching samples
# Option 1: For an apples-to-apples comparison (e.g. ranking omics data), yes, we should prune to matching!
# Option 2: To maximize sample size (i.e. generalizability), we should include all samples where possible

# To compare importance of omics data, we will prune (Option 1)

# here is fcal response:
oars.fcal.response = metadata.oars.stool.double[,c("hm_phase", "lfc.fcal")] %>% unique() %>%
  mutate(fcal.response = ifelse(lfc.fcal < 0, "low", "high")) %>% dplyr::select(hm_phase, fcal.response)
rownames(oars.fcal.response) = oars.fcal.response$hm_phase
oars.fcal.response$hm_phase = NULL
oars.fcal.response
# high = "increase in fcal"

# :: ML Omics -------------------------------------------------

# Goal: compare AUC of omics
data.types = c("ASV", "Species", "Pathway", "COG", "CAZy","Metabolite", "FFQ", "All")
iters = c(1:15)

oars.loocv.models.ml.omics = rf.function(data.types = c("ASV", "Species", "Pathway", "COG", "CAZy","Metabolite", "FFQ", "All"),
                                         iters = 15,
                                         target = "fermentation",
                                         output = "AUC")


# summarize
oars.loocv.models.ml.omics.df = oars.loocv.models.ml.omics %>%
  group_by(data.type) %>%
  dplyr::select(data.type, auc) %>% distinct() %>%
  mutate(mean.auc = mean(auc),
         median.auc = median(auc),
         auc.low = mean(auc) - (sd(auc)/sqrt(n()) * 1.96), # 95% CI
         auc.high = mean(auc) + (sd(auc)/sqrt(n()) * 1.96)) %>%
  dplyr::select(mean.auc, median.auc, auc.low, auc.high, data.type) %>% distinct() %>%
  arrange(-mean.auc)
oars.loocv.models.ml.omics.df
# Nice, reproducible

oars.loocv.models.ml.omics.rf.plot = ggplot(oars.loocv.models.ml.omics.df %>%
                                              mutate(data.type = ifelse(data.type == "All", "Multi-Omic", data.type))%>%
         mutate(data.type = factor(data.type, levels= c("FFQ", "ASV", "Species", "COG","Pathway", "CAZy","Metabolite", "Multi-Omic"))),
       aes(y=reorder(data.type, mean.auc), x=mean.auc))+
  geom_bar(stat="identity", width=0.75,
           aes(fill=data.type),color="black")+
  geom_segment(aes(x=auc.low, xend=auc.high))+
  scale_x_continuous(breaks=seq(0,1, by=0.1),
                     limits=c(0,1))+
  geom_text(data=oars.loocv.models.ml.omics.df %>%
              mutate(data.type = ifelse(data.type == "All", "Multi-Omic", data.type))%>%
              mutate(data.type = factor(data.type, levels= c("FFQ", "ASV", "Species", "COG","Pathway", "CAZy","Metabolite", "Multi-Omic"))), 
            aes(x=auc.high, y=data.type, 
                label=paste(" ", round(mean.auc, digits=3))), hjust=0, size=3)+
  scale_fill_manual(values=c(omics.colors, "Multi-Omic" = "black", "FFQ" = "white"))+
  #geom_segment()+
  facet_wrap(~"RandomForest")+
  theme_classic()+theme(legend.position="none",
                        strip.text=element_text(size=10))+
  labs(x="Mean AUC", y="")
oars.loocv.models.ml.omics.rf.plot

# :: RF LOOCV ----------------------------------------------------------------

# Now that we know Multi-Omic has highest AUC, let's plot it's ROC curve

oars.loocv.models.rf = rf.function(data.types = "All",
                                   iters = 15,
                                   output = "ROC")

oars.loocv.models.ml.omics.fcal.df

# calculate values and extract best
oars.loocv.models.df = oars.loocv.models.rf %>%
  group_by(iter) %>%
  summarize(auc = pROC::auc(high, pred, 
                            levels=c("high", "low"),  # define case = "high" delta pH
                            direction="<")[1])%>%  
  mutate(mean.auc = mean(auc),
         median.auc = median(auc),
         auc.low = mean(auc) - (sd(auc)/sqrt(n()) * 1.96), # 95% CI
         auc.high = mean(auc) + (sd(auc)/sqrt(n()) * 1.96)) %>%
  dplyr::select(mean.auc, median.auc, auc.low, auc.high) %>% distinct()

# extract best
oars.loocv.models.best.stats = subset(oars.loocv.models.df, mean.auc == max(oars.loocv.models.df$mean.auc))

# plot
oars.loocv.roc = do.call(rbind, lapply(1:15, function(seed){
  data.subset = subset(oars.loocv.models.rf, iter == seed)
  # calibrate pred; not necessary
  #data.subset$pred = 1 / (1 + exp(-data.subset$pred))
  data.frame(sens = pROC::roc(data.subset$high, data.subset$pred)$sensitivities,
             spec = pROC::roc(data.subset$high, data.subset$pred)$specificities,
             iter = seed)
}))
  

# Step 2: Compute mean and standard error for sensitivities across iterations
roc_summary <- oars.loocv.roc %>%
  group_by(spec) %>%  # Group by specificity (or alternatively by sens)
  summarise(
    mean_sens = mean(sens, na.rm = TRUE),
    se_sens = sd(sens, na.rm = TRUE) / sqrt(n()),  # Standard error
    lower = mean_sens - 1.96 * se_sens,           # 95% CI lower bound
    upper = mean_sens + 1.96 * se_sens            # 95% CI upper bound
  ) %>%
  mutate(fpr = 1 - spec) %>%  # False Positive Rate (1 - specificity)
  filter(!is.na(mean_sens) & !is.na(se_sens))  # Remove any NA values


oars.rf.loocv.roc.plot = ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  # Individual ROC curves for each iteration
  #geom_line(data = oars.rf.loocv.roc %>% arrange(sens),
  #          aes(x=1-spec, y=sens, group=iter))+
  # Error ribbon (95% CI)
  geom_ribbon(data = roc_summary, 
              aes(x = fpr, ymin = lower, ymax = upper), 
              #fill = RColorBrewer::brewer.pal(n=5, "Set3")[1], alpha = 0.2) +
              fill = "black", alpha=0.2)+
  # Mean ROC curve
  geom_path(data = roc_summary, 
            aes(x = fpr, y = mean_sens), 
            #color = RColorBrewer::brewer.pal(n=5, "Set3")[1], size = 1) +
            color="black")+
  # add label
  annotate(geom="text", x=0.75, y=0.25,
           label=paste("15x LOOCV\n",
                       "AUC: ", round(oars.loocv.models.best.stats$mean.auc, digits=2),
                 "\n(", round(oars.loocv.models.best.stats$auc.low, digits=2), 
                 ", ", round(oars.loocv.models.best.stats$auc.high,digits=2), ")", sep=""),
           size=4)+
  theme_classic()+theme(strip.text=element_text(size=10))+
  facet_wrap(~"Multi-Omic RandomForest")+
  xlim(0,1)+
  ylim(0,1)+
  labs(x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)") 
oars.rf.loocv.roc.plot



# :: RF Importances -------------------------------------------------------


oars.loocv.models.importances = rf.function(data.types = "All",
                                            iters = 15,
                                            output = "importances")
# check for HM0924-STL-12

oars.rf.loocv.imp.df = oars.loocv.models.importances %>%
  #subset(imp != 0) %>%
  group_by(feature) %>%
  mutate(n.present = n()) %>%
  mutate(mean.imp = mean(na.omit(imp))) %>%
  mutate(imp.low = mean(na.omit(imp)) - (sd(na.omit(imp))/sqrt(n()) * 1.96)) %>% # 95% CI
  mutate(imp.high = mean(na.omit(imp))+ (sd(na.omit(imp))/sqrt(n()) * 1.96)) %>% # 95% CI
  subset(mean.imp != 0) %>%
  dplyr::select(feature, mean.imp, imp.low, imp.high, n.present) %>% distinct() %>%
  arrange(-mean.imp) %>% data.frame()

# Note: need to be careful selecting which features to highlight,
# as their frequency of being selected impacts results considerably
ggplot(oars.rf.loocv.imp.df,
       aes(x=n.present, y=mean.imp))+
  scale_x_log10()+geom_smooth()+
  geom_point()+theme_classic()


# log2fc
oars.rf.loocv.imp.df = oars.rf.loocv.imp.df %>% 
  subset(!is.na(mean.imp))%>%
  # select the top 10 mean importance
  arrange(mean.imp) %>% slice_max(mean.imp, n=10) 

# calculate wilcox test pvalue of these features between response groups
oars.rf.loocv.imp.wilcox = do.call(rbind, lapply(oars.rf.loocv.imp.df$feature, function(x){
  data.subset = oars.all.omics.baseline[,colnames(oars.all.omics.baseline) %in% c(x, "response")]
  colnames(data.subset)[2] = "variable"
  # run wilcox
  wilcox.p = wilcox.test(subset(data.subset, response == "high")[,2],
                         subset(data.subset, response == "low")[,2])
  
  # calculate lfc
  if(x %in% c("stool_water_perc", "shannon", "fd")){
    coef = mean(subset(data.subset, response == "high")[,2]) / mean(subset(data.subset, response == "low")[,2])
  }else{
    coef = mean(log2(subset(data.subset, response == "high")[,2])) - mean(log2(subset(data.subset, response == "low")[,2]))
  }
  data.frame(feature = x,
             pval = wilcox.p$p.value,
             coef = coef)
}))
oars.rf.loocv.imp.df = merge(oars.rf.loocv.imp.wilcox,
                                       oars.rf.loocv.imp.df, by="feature")
oars.rf.loocv.imp.df$sig = ifelse(oars.rf.loocv.imp.df$pval < 0.05, "*", "")
oars.rf.loocv.imp.df = oars.rf.loocv.imp.df %>% arrange(mean.imp)

# clean up COGs
oars.rf.loocv.imp.df$clean.feature = gsub(paste(c(" and related.*", "or related.*", ", .*"), collapse="|"),"", oars.rf.loocv.imp.df$feature)

oars.rf.loocv.imp.plot = ggplot(oars.rf.loocv.imp.df,
                                             aes(x=mean.imp, y=reorder(clean.feature, mean.imp)))+
  geom_segment(aes(x=imp.low, xend=imp.high, y=clean.feature, yend=clean.feature), color="black")+
  geom_point(shape=21, aes(fill=(coef)), size=3.5)+
  scale_fill_gradient2(low="blue", high="red")+
  scale_color_manual(values=c("white", "black"))+
  theme_classic()+theme(legend.position="right",
    axis.title.y = element_blank(),
    #axis.text.y = element_blank(),
   # axis.ticks.y=element_blank(),
    plot.title = element_text(hjust = 0.5, size=12),
    strip.text = element_text(size=10),
    strip.background = element_rect(
      color="black"))+
  guides(color=FALSE)+
  labs(x="Mean Decrease in Accuracy", fill="Log2FC")+
  facet_wrap(~"Feature Importance")
oars.rf.loocv.imp.plot


# :: RF ttest plots -------------------------------------------------------


# loop through top features and plot; extract microbiome values
oars.rf.loocv.imp.ttest = do.call(rbind, lapply(slice_max(oars.rf.loocv.imp.df, mean.imp, n=10)$feature, function(x){
  data.subset = oars.all.omics.baseline[,c(x, "response")]
  rownames(data.subset) = paste(oars.all.omics.baseline$HM, oars.all.omics.baseline$phase, sep="_")
  data.subset$feature = colnames(data.subset)[1]
  data.subset$sample = rownames(data.subset)
  colnames(data.subset)[1] = "value"
  data.subset
}))
# add stats
oars.rf.loocv.imp.ttest = merge(oars.rf.loocv.imp.ttest,
                                oars.rf.loocv.imp.df, by="feature")
oars.rf.loocv.imp.ttest.p = oars.rf.loocv.imp.ttest[,c("feature", "pval")] %>% distinct()

# clean taxa names (so they fit)
oars.rf.loocv.imp.ttest$clean.feature = gsub(paste(c(" and related.*", "or related.*", ", .*"), collapse="|"),"", oars.rf.loocv.imp.ttest$feature)
oars.rf.loocv.imp.ttest.p$clean.feature = gsub(paste(c(" and related.*", "or related.*", ", .*"), collapse="|"),"", oars.rf.loocv.imp.ttest.p$feature)

# now clip names to 30 char
oars.rf.loocv.imp.ttest$short.feature = ifelse(nchar(oars.rf.loocv.imp.ttest$clean.feature)>30, paste(substr(oars.rf.loocv.imp.ttest$clean.feature, start = 1, stop = 30), "...", sep=""), oars.rf.loocv.imp.ttest$clean.feature)
oars.rf.loocv.imp.ttest.p$short.feature = ifelse(nchar(oars.rf.loocv.imp.ttest.p$clean.feature)>30, paste(substr(oars.rf.loocv.imp.ttest.p$clean.feature, start = 1, stop = 30), "...", sep=""), oars.rf.loocv.imp.ttest.p$clean.feature)
oars.rf.loocv.imp.df$short.feature = ifelse(nchar(oars.rf.loocv.imp.df$clean.feature)>30, paste(substr(oars.rf.loocv.imp.df$clean.feature, start = 1, stop = 30), "...", sep=""), oars.rf.loocv.imp.df$clean.feature)

oars.rf.loocv.imp.ttest.plot = ggplot(oars.rf.loocv.imp.ttest %>%
                                        # clean names
                                        mutate(Response = ifelse(response == "low", "Weak\nResponse", "Strong\nResponse")) %>%
                                        # add indicator of not present
                                        group_by(short.feature) %>%
                                        mutate(pseudo = ifelse(value == min(value), "pseudo", "real")) %>%
                                        # reorder taxa based on coef
                                        mutate(short.feature = factor(short.feature, levels=arrange(distinct(oars.rf.loocv.imp.df[,c("short.feature", "coef")]),-coef)$short.feature)),
                                                    aes(x=Response, y=value))+
  geom_boxplot(width=0.3, outlier.shape=NA)+
  #scale_y_log10(labels = scales::label_number(accuracy = 0.01))+
  scale_y_log10()+
  ggbeeswarm::geom_beeswarm(shape=21, aes(fill=response, alpha=pseudo), size=2)+
  scale_alpha_manual(values=c(0.2, 1))+
  theme_classic()+theme(legend.position="none",
                        axis.title.x = element_blank(),
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=8),
                        strip.background = element_rect(color="black"))+
  geom_text(data=oars.rf.loocv.imp.ttest.p%>%
              # clean names
              mutate(short.feature = factor(short.feature, levels=arrange(distinct(oars.rf.loocv.imp.df[,c("short.feature", "mean.imp")]),-mean.imp)$short.feature)),
              x=1.5, y=Inf, vjust=1.2, 
            aes(label = paste("p =", round(pval, digits=3))),
            size=2.5)+
  labs(x="", y="Feature Abundance")+
  facet_wrap(~short.feature, ncol=2 ,scales="free_y")
oars.rf.loocv.imp.ttest.plot


# plots so far

oars.loocv.models.ml.omics.rf.plot+
  oars.rf.loocv.roc.plot+
  oars.rf.loocv.imp.plot+
  oars.rf.loocv.imp.ttest.plot

# :: RF Interactions ------------------------------------------------------

# check interactions of most important features
# (otherwise, intractable to compute, even on compute canada)

oars.loocv.models.interactions = rf.function(data.types = "All",
                                             iters = 15,
                                             output = "interactions")

# analyze
ints.df = 
  # watch this
  rbind(oars.loocv.models.interactions %>% mutate(var3 = var1,
                                                  var1 = var2,
                                                  var2 = var3) %>% dplyr::select(-var3),
        oars.loocv.models.interactions) %>% data.frame() %>%
  # take mean
  group_by(features) %>%
  mutate(value = mean(Difference))%>%
  dplyr::select(`var1`, `var2`, value)

ints.df = reshape2::acast(ints.df, var1~var2, value.var="value", fun.aggregate=mean)%>%
  as.matrix() %>% reshape2::melt()

ints.df$value %>% range()

# clean feature names // defunct, but keep anyway for clarity
ints.df$clean.feature.1 = ints.df$Var1
ints.df$clean.feature.2 = ints.df$Var2
# shorten
ints.df$clean.feature.1 = gsub(paste(c(" and related.*", "or related.*", ", .*"), collapse="|"),"", ints.df$clean.feature.1)
ints.df$clean.feature.2 = gsub(paste(c(" and related.*", "or related.*", ", .*"), collapse="|"),"", ints.df$clean.feature.2)
# clip
ints.df$short.feature.1 = ifelse(nchar(ints.df$clean.feature.1)>30, paste(substr(ints.df$clean.feature.1, start = 1, stop = 30), "...", sep=""), ints.df$clean.feature.1)
ints.df$short.feature.2 = ifelse(nchar(ints.df$clean.feature.2)>30, paste(substr(ints.df$clean.feature.2, start = 1, stop = 30), "...", sep=""), ints.df$clean.feature.2)


# plot
oars.rf.loocv.interactions.plot = ints.df %>%
  #mutate(ratio = log10(Paired - Additive)) %>%
  group_by(Var2) %>%
  mutate(sum.1 = max(na.omit(value)))%>%
  ungroup()%>%
  subset(Var1 %in% oars.rf.loocv.imp.df$feature) %>% distinct()%>%
  #slice_max(order_by=na.omit(value), n=10) %>%
  subset(value > 0) %>%
  dplyr::select(clean.feature.1, short.feature.2, sum.1, value) %>% 
  rbind(data.frame(clean.feature.1 = oars.rf.loocv.imp.df$clean.feature,
                   short.feature.2 = "Marginal Importance",
                   sum.1 = 11,
                   value = oars.rf.loocv.imp.df$mean.imp)) %>%
  mutate(clean.feature.1 = factor(clean.feature.1, levels=(oars.rf.loocv.imp.df$clean.feature))) %>%
  mutate(short.feature.2 = factor(short.feature.2, levels=rev(c(oars.rf.loocv.imp.df$short.feature,"Marginal Importance")))) %>%
  ggplot(
  aes(y=clean.feature.1, x=short.feature.2))+
  geom_tile(aes(fill=scale(value)), color="white", size=2)+
  geom_vline(xintercept=1.5, color="black")+
  scale_fill_gradient2(low="blue", high="red")+
  facet_wrap(~"Interactions")+
  theme_classic()+theme(axis.text.x=element_text(angle=45,hjust=1),
                        strip.text=element_text(size=10),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank())+
  labs(fill="Scaled\nImportance")
oars.rf.loocv.interactions.plot
# fixed!
# Difference = Conditional Importance - sum(Marginal Importance)
# (negative Difference values = anti-conditional interactions)

oars.rf.loocv.interactions.plot+
oars.rf.loocv.imp.ttest.plot

# :: Predictors - PCoA ---------------------------------------------------------

# ASV
oars.asv.baseline.map = oars.asv.baseline[,c("HM", "phase", "response", "diagnosis", "adj.fiber")]
rownames(oars.asv.baseline.map) = paste(oars.asv.baseline.map$HM, oars.asv.baseline.map$phase, sep="_")
rownames(oars.asv.baseline) = rownames(oars.asv.baseline.map)
# calculate Bray-Curtis dissimilarities
oars.asv.baseline.bray = vegan::vegdist(oars.asv.baseline[,!colnames(oars.asv.baseline) %in% c("HM", "phase", "response", "diagnosis", "adj.fiber")], method="bray") 
# perform PCoA
oars.asv.baseline.bray.pcoa = ape::pcoa(oars.asv.baseline.bray)
# extract data from pcoa
oars.asv.baseline.bray.pcoa.df = data.frame(oars.asv.baseline.bray.pcoa$vectors[,c(1:2)])
# add metadata
oars.asv.baseline.bray.pcoa.df = merge(oars.asv.baseline.bray.pcoa.df,
                         oars.asv.baseline.map, by="row.names")
# extract variance explained
oars.asv.baseline.bray.pcoa.var = oars.asv.baseline.bray.pcoa$values[c(1:2),2]
oars.asv.baseline.bray.pcoa.df$var1 = round(oars.asv.baseline.bray.pcoa.var[1]*100, digits=2)
oars.asv.baseline.bray.pcoa.df$var2 = round(oars.asv.baseline.bray.pcoa.var[2]*100, digits=2)
rownames(oars.asv.baseline.bray.pcoa.df) = oars.asv.baseline.bray.pcoa.df$Row.names
# clean up "oars.on.rs" variable

set.seed(25)
t1 = Sys.time()
oars.asv.baseline.bray.permanova = vegan::adonis2(oars.asv.baseline.bray ~ response,
                                    oars.asv.baseline.bray.pcoa.df,
                                    #strata = oars.asv.baseline.bray.pcoa.df$HM,
                                    by="margin")
t2 = Sys.time()
t2 - t1
oars.asv.baseline.bray.pcoa.df # not sig

oars.asv.baseline.bray.plot <- ggplot(
  data=oars.asv.baseline.bray.pcoa.df %>% group_by(HM),
  aes(x=Axis.1, y=Axis.2))+
  geom_path(aes(group=HM), color="black", linetype=2, alpha=0.5, linewidth=0.3) + 
  geom_point(aes(fill=response), shape=21, size=2)+
  stat_ellipse(aes(group=response, color=response), alpha=0.5)+
  annotate(geom="text", x=Inf, y=Inf, hjust=1.2, vjust=1.2, size=3,
           label=paste(paste("R²: ", round(data.frame(oars.asv.baseline.bray.permanova)[1,3], 3)*100, "%",
                             "  p: ", round(data.frame(oars.asv.baseline.bray.permanova)[1,5], 3), sep=""), sep=""))+
  facet_wrap(~"ASV Taxa")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  labs(x=paste("Axis 1: ", round(unique(oars.asv.baseline.bray.pcoa.df$var1), digits=2), "%", sep=""), 
       y=paste("Axis 2: ", round(unique(oars.asv.baseline.bray.pcoa.df$var2), digits=2), "%", sep=""))
oars.asv.baseline.bray.plot


# MGX
oars.mgx.baseline.map = oars.mgx.baseline[,c("HM", "phase", "response", "diagnosis", "adj.fiber")]
rownames(oars.mgx.baseline.map) = paste(oars.mgx.baseline.map$HM, oars.mgx.baseline.map$phase, sep="_")
rownames(oars.mgx.baseline) = rownames(oars.mgx.baseline.map)
# calculate Bray-Curtis dissimilarities
oars.mgx.baseline.bray = vegan::vegdist(oars.mgx.baseline[,!colnames(oars.mgx.baseline) %in% c("HM", "phase", "response", "diagnosis", "adj.fiber")], method="bray") 
# perform PCoA
oars.mgx.baseline.bray.pcoa = ape::pcoa(oars.mgx.baseline.bray)
# extract data from pcoa
oars.mgx.baseline.bray.pcoa.df = data.frame(oars.mgx.baseline.bray.pcoa$vectors[,c(1:2)])
# add metadata
oars.mgx.baseline.bray.pcoa.df = merge(oars.mgx.baseline.bray.pcoa.df,
                                       oars.mgx.baseline.map, by="row.names")
# extract variance explained
oars.mgx.baseline.bray.pcoa.var = oars.mgx.baseline.bray.pcoa$values[c(1:2),2]
oars.mgx.baseline.bray.pcoa.df$var1 = round(oars.mgx.baseline.bray.pcoa.var[1]*100, digits=2)
oars.mgx.baseline.bray.pcoa.df$var2 = round(oars.mgx.baseline.bray.pcoa.var[2]*100, digits=2)
rownames(oars.mgx.baseline.bray.pcoa.df) = oars.mgx.baseline.bray.pcoa.df$Row.names
# clean up "oars.on.rs" variable

set.seed(25)
t1 = Sys.time()
oars.mgx.baseline.bray.permanova = vegan::adonis2(oars.mgx.baseline.bray ~ response,
                                                  oars.mgx.baseline.bray.pcoa.df,
                                                  #strata = oars.mgx.baseline.bray.pcoa.df$HM,
                                                  by="margin")
t2 = Sys.time()
t2 - t1
oars.mgx.baseline.bray.pcoa.df # not sig

oars.mgx.baseline.bray.plot <- ggplot(
  data=oars.mgx.baseline.bray.pcoa.df %>% group_by(HM),
  aes(x=Axis.1, y=Axis.2))+
  geom_path(aes(group=HM), color="black", linetype=2, alpha=0.5, linewidth=0.3) + 
  geom_point(aes(fill=response), shape=21, size=2)+
  stat_ellipse(aes(group=response, color=response), alpha=0.5)+
  annotate(geom="text", x=Inf, y=Inf, hjust=1.2, vjust=1.2,size=3,
           label=paste(paste("R²: ", round(data.frame(oars.mgx.baseline.bray.permanova)[1,3], 3)*100, "%",
                             "  p: ", round(data.frame(oars.mgx.baseline.bray.permanova)[1,5], 3), sep=""), sep=""))+
  facet_wrap(~"Species")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  labs(x=paste("Axis 1: ", round(unique(oars.mgx.baseline.bray.pcoa.df$var1), digits=2), "%", sep=""), 
       y=paste("Axis 2: ", round(unique(oars.mgx.baseline.bray.pcoa.df$var2), digits=2), "%", sep=""))
       
oars.mgx.baseline.bray.plot



# MPX KEGG
oars.mpx.kegg.baseline.map = oars.mpx.kegg.baseline[,c("HM", "phase", "response", "diagnosis", "adj.fiber")]
rownames(oars.mpx.kegg.baseline.map) = paste(oars.mpx.kegg.baseline.map$HM, oars.mpx.kegg.baseline.map$phase, sep="_")
rownames(oars.mpx.kegg.baseline) = rownames(oars.mpx.kegg.baseline.map)
# Perform PCA (not PCoA)
oars.mpx.kegg.baseline.pca = prcomp(log2(oars.mpx.kegg.baseline[,!colnames(oars.mpx.kegg.baseline) %in% c("HM", "phase", "response", "diagnosis", "adj.fiber")]))

# extract data from pcoa
oars.mpx.kegg.baseline.pca.df = data.frame(oars.mpx.kegg.baseline.pca$x[,c(1:2)])
rownames(oars.mpx.kegg.baseline.pca.df) = rownames(oars.mpx.kegg.baseline)
# add metadata
oars.mpx.kegg.baseline.pca.df = merge(oars.mpx.kegg.baseline.pca.df,
                                      oars.mpx.kegg.baseline.map, by="row.names")
# extract variance explained
oars.mpx.kegg.baseline.pca.varexp = (oars.mpx.kegg.baseline.pca$sdev)^2 / sum(oars.mpx.kegg.baseline.pca$sdev^2) * 100

rownames(oars.mpx.kegg.baseline.pca.df) = oars.mpx.kegg.baseline.pca.df$Row.names
# clean up "oars.on.rs" variable

set.seed(25)
t1 = Sys.time()
oars.mpx.kegg.baseline.pca.permanova = vegan::adonis2(dist(oars.mpx.kegg.baseline.pca$x) ~ response,
                                                 oars.mpx.kegg.baseline.pca.df,
                                                 #strata = oars.mpx.kegg.baseline.bray.pcoa.df$HM,
                                                 by="margin")
t2 = Sys.time()
t2 - t1
oars.mpx.kegg.baseline.pca.permanova # not sig

oars.mpx.kegg.baseline.pca.plot <- ggplot(
  data=oars.mpx.kegg.baseline.pca.df %>% group_by(HM),
  aes(x=PC1, y=PC2))+
  geom_path(aes(group=HM), color="black", linetype=2, alpha=0.5, linewidth=0.3) + 
  geom_point(aes(fill=response), shape=21, size=2)+
  stat_ellipse(aes(group=response, color=response), alpha=0.5)+
  annotate(geom="text", x=Inf, y=Inf, hjust=1.2, vjust=1.2,size=3,
           label=paste(paste("R²: ", round(data.frame(oars.mpx.kegg.baseline.pca.permanova)[1,3], 3)*100, "%",
                             "  p: ", round(data.frame(oars.mpx.kegg.baseline.pca.permanova)[1,5], 3), sep=""), sep=""))+
  facet_wrap(~"Pathway")+
  theme_classic()+ theme(legend.position="none",
                         plot.title = element_text(hjust = 0.5, size=12),
                         strip.text = element_text(size=10),
                         strip.background = element_rect(
                           color="black"))+
  labs(x=paste("PC1: ", round(unique(oars.mpx.kegg.baseline.pca.varexp[1]), digits=2), "%", sep=""), 
       y=paste("PC2: ", round(unique(oars.mpx.kegg.baseline.pca.varexp[2]), digits=2), "%", sep=""))
oars.mpx.kegg.baseline.pca.plot



# MPX
oars.mpx.cog.baseline.map = oars.mpx.cog.baseline[,c("HM", "phase", "response", "diagnosis", "adj.fiber")]
rownames(oars.mpx.cog.baseline.map) = paste(oars.mpx.cog.baseline.map$HM, oars.mpx.cog.baseline.map$phase, sep="_")
rownames(oars.mpx.cog.baseline) = rownames(oars.mpx.cog.baseline.map)
# Perform PCA (not PCoA)
oars.mpx.cog.baseline.pca = prcomp(log2(oars.mpx.cog.baseline[,!colnames(oars.mpx.cog.baseline) %in% c("HM", "phase", "response", "diagnosis", "adj.fiber")]))

# extract data from pcoa
oars.mpx.cog.baseline.pca.df = data.frame(oars.mpx.cog.baseline.pca$x[,c(1:2)])
rownames(oars.mpx.cog.baseline.pca.df) = rownames(oars.mpx.cog.baseline)
# add metadata
oars.mpx.cog.baseline.pca.df = merge(oars.mpx.cog.baseline.pca.df,
                                     oars.mpx.cog.baseline.map, by="row.names")
# extract variance explained
oars.mpx.cog.baseline.pca.varexp = (oars.mpx.cog.baseline.pca$sdev)^2 / sum(oars.mpx.cog.baseline.pca$sdev^2) * 100

rownames(oars.mpx.cog.baseline.pca.df) = oars.mpx.cog.baseline.pca.df$Row.names
# clean up "oars.on.rs" variable

set.seed(25)
t1 = Sys.time()
oars.mpx.cog.baseline.pca.permanova = vegan::adonis2(dist(oars.mpx.cog.baseline.pca$x) ~ response,
                                                 oars.mpx.cog.baseline.pca.df,
                                                 #strata = oars.mpx.cog.baseline.bray.pcoa.df$HM,
                                                 by="margin")
t2 = Sys.time()
t2 - t1
oars.mpx.cog.baseline.pca.permanova # not sig

oars.mpx.cog.baseline.pca.plot <- ggplot(
  data=oars.mpx.cog.baseline.pca.df %>% group_by(HM),
  aes(x=PC1, y=PC2))+
  geom_path(aes(group=HM), color="black", linetype=2, alpha=0.5, linewidth=0.3) + 
  geom_point(aes(fill=response), shape=21, size=2)+
  stat_ellipse(aes(group=response, color=response), alpha=0.5)+
  annotate(geom="text", x=Inf, y=Inf, hjust=1.2, vjust=1.2,size=3,
           label=paste(paste("R²: ", round(data.frame(oars.mpx.cog.baseline.pca.permanova)[1,3], 3)*100, "%",
                             "  p: ", round(data.frame(oars.mpx.cog.baseline.pca.permanova)[1,5], 3), sep=""), sep=""))+
  facet_wrap(~"COG")+
  theme_classic()+ theme(legend.position="none",
                         plot.title = element_text(hjust = 0.5, size=12),
                         strip.text = element_text(size=10),
                         strip.background = element_rect(
                           color="black"))+
  labs(x=paste("PC1: ", round(unique(oars.mpx.cog.baseline.pca.varexp[1]), digits=2), "%", sep=""), 
       y=paste("PC2: ", round(unique(oars.mpx.cog.baseline.pca.varexp[2]), digits=2), "%", sep=""))
oars.mpx.cog.baseline.pca.plot


 
# MPX CAZy
oars.mpx.cazy.baseline.map = oars.mpx.cazy.baseline[,c("HM", "phase", "response", "diagnosis", "adj.fiber")]
rownames(oars.mpx.cazy.baseline.map) = paste(oars.mpx.cazy.baseline.map$HM, oars.mpx.cazy.baseline.map$phase, sep="_")
rownames(oars.mpx.cazy.baseline) = rownames(oars.mpx.cazy.baseline.map)
# Perform PCA (not PCoA)
oars.mpx.cazy.baseline.pca = prcomp(log2(oars.mpx.cazy.baseline[,!colnames(oars.mpx.cazy.baseline) %in% c("HM", "phase", "response", "diagnosis", "adj.fiber")]))

# extract data from pcoa
oars.mpx.cazy.baseline.pca.df = data.frame(oars.mpx.cazy.baseline.pca$x[,c(1:2)])
rownames(oars.mpx.cazy.baseline.pca.df) = rownames(oars.mpx.cazy.baseline)
# add metadata
oars.mpx.cazy.baseline.pca.df = merge(oars.mpx.cazy.baseline.pca.df,
                                      oars.mpx.cazy.baseline.map, by="row.names")
# extract variance explained
oars.mpx.cazy.baseline.pca.varexp = (oars.mpx.cazy.baseline.pca$sdev)^2 / sum(oars.mpx.cazy.baseline.pca$sdev^2) * 100

rownames(oars.mpx.cazy.baseline.pca.df) = oars.mpx.cazy.baseline.pca.df$Row.names
# clean up "oars.on.rs" variable

set.seed(25)
t1 = Sys.time()
oars.mpx.cazy.baseline.pca.permanova = vegan::adonis2(dist(oars.mpx.cazy.baseline.pca$x) ~ response,
                                                 oars.mpx.cazy.baseline.pca.df,
                                                 #strata = oars.mpx.cazy.baseline.bray.pcoa.df$HM,
                                                 by="margin")
t2 = Sys.time()
t2 - t1
oars.mpx.cazy.baseline.pca.permanova # not sig

oars.mpx.cazy.baseline.pca.plot <- ggplot(
  data=oars.mpx.cazy.baseline.pca.df %>% group_by(HM),
  aes(x=PC1, y=PC2))+
  geom_path(aes(group=HM), color="black", linetype=2, alpha=0.5, linewidth=0.3) + 
  geom_point(aes(fill=response), shape=21, size=2)+
  stat_ellipse(aes(group=response, color=response), alpha=0.5)+
  annotate(geom="text", x=Inf, y=Inf, hjust=1.2, vjust=1.2,size=3,
           label=paste(paste("R²: ", round(data.frame(oars.mpx.cazy.baseline.pca.permanova)[1,3], 3)*100, "%",
                             "  p: ", round(data.frame(oars.mpx.cazy.baseline.pca.permanova)[1,5], 3), sep=""), sep=""))+
  facet_wrap(~"CAZy")+
  theme_classic()+ theme(legend.position="none",
                         plot.title = element_text(hjust = 0.5, size=12),
                         strip.text = element_text(size=10),
                         strip.background = element_rect(
                           color="black"))+
  labs(x=paste("PC1: ", round(unique(oars.mpx.cazy.baseline.pca.varexp[1]), digits=2), "%", sep=""), 
       y=paste("PC2: ", round(unique(oars.mpx.cazy.baseline.pca.varexp[2]), digits=2), "%", sep=""))
oars.mpx.cazy.baseline.pca.plot


# MBX
oars.mbx.baseline.map = oars.mbx.baseline[,c("HM", "phase", "response", "diagnosis", "adj.fiber")]
rownames(oars.mbx.baseline.map) = paste(oars.mbx.baseline.map$HM, oars.mbx.baseline.map$phase, sep="_")
rownames(oars.mbx.baseline) = rownames(oars.mbx.baseline.map)
# Perform PCA (not PCoA)
oars.mbx.baseline.pca = prcomp(log2(oars.mbx.baseline[,!colnames(oars.mbx.baseline) %in% c("HM", "phase", "response", "diagnosis", "adj.fiber")]))

# extract data from pcoa
oars.mbx.baseline.pca.df = data.frame(oars.mbx.baseline.pca$x[,c(1:2)])
rownames(oars.mbx.baseline.pca.df) = rownames(oars.mbx.baseline)
# add metadata
oars.mbx.baseline.pca.df = merge(oars.mbx.baseline.pca.df,
                                      oars.mbx.baseline.map, by="row.names")
# extract variance explained
oars.mbx.baseline.pca.varexp = (oars.mbx.baseline.pca$sdev)^2 / sum(oars.mbx.baseline.pca$sdev^2) * 100

rownames(oars.mbx.baseline.pca.df) = oars.mbx.baseline.pca.df$Row.names
# clean up "oars.on.rs" variable

set.seed(25)
t1 = Sys.time()
oars.mbx.baseline.pca.permanova = vegan::adonis2(dist(oars.mbx.baseline.pca$x) ~ response,
                                                      oars.mbx.baseline.pca.df,
                                                      #strata = oars.mbx.baseline.bray.pcoa.df$HM,
                                                      by="margin")
t2 = Sys.time()
t2 - t1
oars.mbx.baseline.pca.permanova # not sig

oars.mbx.baseline.pca.plot <- ggplot(
  data=oars.mbx.baseline.pca.df %>% group_by(HM),
  aes(x=PC1, y=PC2))+
  geom_path(aes(group=HM), color="black", linetype=2, alpha=0.5, linewidth=0.3) + 
  geom_point(aes(fill=response), shape=21, size=2)+
  stat_ellipse(aes(group=response, color=response), alpha=0.5)+
  annotate(geom="text", x=Inf, y=Inf, hjust=1.2, vjust=1.2,size=3,
           label=paste(paste("R²: ", round(data.frame(oars.mbx.baseline.pca.permanova)[1,3], 3)*100, "%",
                             "  p: ", round(data.frame(oars.mbx.baseline.pca.permanova)[1,5], 3), sep=""), sep=""))+
  facet_wrap(~"Metabolite")+
  theme_classic()+ theme(legend.position="none",
                         plot.title = element_text(hjust = 0.5, size=12),
                         strip.text = element_text(size=10),
                         strip.background = element_rect(
                           color="black"))+
  labs(x=paste("PC1: ", round(unique(oars.mbx.baseline.pca.varexp[1]), digits=2), "%", sep=""), 
       y=paste("PC2: ", round(unique(oars.mbx.baseline.pca.varexp[2]), digits=2), "%", sep=""))
oars.mbx.baseline.pca.plot



p1 = ((oars.asv.baseline.bray.plot+
    oars.mgx.baseline.bray.plot+
     # oars.mpx.protein.baseline.pca.plot+
      oars.mpx.kegg.baseline.pca.plot+
      oars.mpx.cog.baseline.pca.plot+
      oars.mpx.cazy.baseline.pca.plot+
      oars.mbx.baseline.pca.plot)+patchwork::plot_layout(nrow=3))

p2 = (patchwork::free(oars.loocv.models.ml.omics.rf.plot, type = "label")+
                        patchwork::free(oars.rf.loocv.imp.plot, type = "label")+
                                          patchwork::free(oars.rf.loocv.roc.plot, type = "label")+
                                                            patchwork::free(oars.rf.loocv.interactions.plot, type = "label")+
        patchwork::plot_layout(widths=c(1,1)))

p3 = oars.rf.loocv.imp.ttest.plot

(patchwork::free(p1, type = "label")| 
    patchwork::free(p2, type = "label") | 
    patchwork::free(p3, type = "label")) + patchwork::plot_layout(nrow=1, widths=c(1,2,1))

# Might need to present ASV as main plot
# and include multi-omic as supplemental to discuss other important features


# >>> SUPP ANALYSES -----------------------------------------------


# :: General correlations -------------------------------------------------

# interpet some important features using correlations
subset(all.features.cor.all, grepl("Flagellar motor switch/", Var1) & grepl("faecis", Var2))
subset(all.features.cor.all, grepl("Flagellar motor switch/", Var1) & Var2=="fcal")
# sig associated with Roseburia faecis (R = 0.67, FDR < 0.001)
subset(all.features.cor.all, grepl("Flagellar biosynthesis protein FliR", Var1) & grepl("Bitt", Var2)) 
subset(all.features.cor.all, grepl("Flagellar biosynthesis protein FliR", Var1) & grepl("fcal", Var2)) 
# most sig associated with Bittarella massiliensis (R = 0.75, FDR < 0.001)
subset(all.features.cor.all, grepl("IMP dehydro", Var1) & grepl("adol", Var2)) %>% subset(padj < 0.05)# %>% head(n=25)
# R = 0.65, FDR 0.00056
subset(all.features.cor.all, grepl("IMP dehydro", Var1) & grepl("fcal", Var2)) %>% subset(padj < 0.20)# %>% head(n=25)

subset(all.features.cor.all, grepl("UDP-N-acetylmuramyl pentapeptide phosphotransferase", Var1) & Var2 == "fcal") %>% subset(padj < 0.20)# %>% head(n=25)

subset(all.features.cor.all, grepl("Holdemania_filiformis", Var1) & Var2 == "fcal") %>% subset(padj < 0.20)# %>% head(n=25)
subset(all.features.cor.all, grepl("Holdemania_filiformis", Var1)) %>% subset(padj < 0.05)# %>% head(n=25)
subset(all.features.cor.all, grepl("adole", Var1) & grepl("fcal", Var2)) %>% subset(padj < 0.05)# %>% head(n=25)
subset(all.features.cor.all, grepl("adole", Var1) & grepl("Magn", Var2)) %>% subset(padj < 0.20)# %>% head(n=25)


subset(all.features.cor.all, grepl("Hexaethyleneglycol mono-n-tetradecyl", Var1)) %>% subset(padj < 0.05) %>% arrange(-value) 
# interestingly, strongly negatively correlated with B. adolescentis (R=-0.52, FDR = 0.019) and positively correlated with Ruthenibacterium_lactatiformans (R=0.60, FDR = 0.003)
subset(all.features.cor.all, grepl("Octaeth", Var1)) %>% subset(padj < 0.05) %>% arrange(value) 
# most sig negatively correlated with Enterocloster citroniae
subset(all.features.cor.all, grepl("13-Docosenamide | 11.69min", Var1)) %>% subset(padj < 0.05) %>% arrange(value) 
# most sig correlated with Anaeromassilibacillus_sp_An250 (R = 0.68, FDR < 0.001), but several others, too
# including neg correlated with R. gnavus (R = -0.48, FDR < 0.01)

subset(all.features.cor.all, grepl("adolescentis", Var1) & Var2 == "fcal") %>% subset(padj < 0.20) %>% arrange(-value) 
subset(all.features.cor.all, grepl("adolescentis", Var1) & grepl("Magnesium", Var2)) %>% subset(padj < 0.20) %>% arrange(-value) 

subset(all.features.cor.all, grepl("adolescentis", Var1)) %>% subset(padj < 0.05) %>% arrange(-value) 
# only positively correlated with B. bifidum
# but negatively correlated with many (Eisenbergiella, Collinsella, C. leptum, Holdemania, C. butyricum, etc)


subset(all.features.cor.all, grepl("C5-Branch", Var1)) %>% subset(padj < 0.05) %>% arrange(-value)
# extensively correlated with other proteins (n=412 with FDR < 0.05)


# :: Vancomycin == fcal/AIEC --------------------------------------------------


# stats
lmerTest::lmer(scale(log2(`Vancomycin resistance`)) ~ scale(log10(fcal)) + diagnosis + adj.fiber + (1|HM),
                                                     subset(stats.tree.features.kegg.data, timing %in% c("0M", "3M", "6M")& 
                                                              HM %in% c("HM0819", "HM0883", "HM0902","HM0899", "HM0906", "HM0924") &
                                                              compliant == TRUE)) %>%
  summary() %>% coef()

# no association with fcal

all.features.cor.all %>%
  subset(Var1=="fcal" & grepl("Vanco", Var2)) # no association with fecal

all.features.cor.all %>%
  subset(Var1=="fcal" & grepl("", Var2)) # no association with fecal


all.features.cor.all %>%
  subset(Var1=="Vancomycin resistance" & grepl("Escherichia", Var2))
# strong association with Biofilm formation - E. coli
# R = 0.62, FDR = 0.001

all.features.cor.all %>%
  subset(Var1=="Vancomycin resistance" & grepl("Kleb", Var2))
# no association with Klebsiella

all.features.cor.all %>%
  subset(Var1=="Vancomycin resistance" & grepl("Stear", Var2))
# no association with Klebsiella

## plot correlation
oars.vanco.biofilm.plot = ggplot(oars.fecalcal.omics.lfc,
       aes(x=`Vancomycin resistance`,
           y=(`Biofilm formation - Escherichia coli`)))+
  geom_point(shape=21, fill="white", size=2.5)+
  ggpubr::stat_cor(method="spearman")+
  geom_smooth(method="lm", color="black")+
  theme_classic()+
  labs(x="Vancomycin resistance Log2FC",
       y="Biofilm formation - E. coli Log2FC")
oars.vanco.biofilm.plot

# what about the fatty acids
all.features.cor.all %>%
  subset(grepl("Hexaethyleneglycol", Var1)) %>% head()
# Strongly associated with Ruthenibacterium lactatiformans (depleted in IBD, Ning,Hong); R = 0.60, FDR < 0.01

all.features.cor.all %>%
  subset(grepl("Octaeth", Var1)) %>% arrange(-value) %>% head()
# Strongly associated with Outer membrane lipoprotein SlyB; R = 0.51, FDR = 0.052


# :: Energy intake == Bacterial load? -------------------------------------

# Nutrients
redcap.data.nutrients = readRDS("./2025_08_16_oars_nutrient_table.Rds")

redcap.data.nutrients.load =  merge(metadata.oars.stool.asv,
                                    as.data.frame(redcap.data.nutrients %>%
                                                    mutate(standard.name = rownames(.)) %>%
                                                    dplyr::select(standard.name, magnesium, starch, energy)), by="standard.name")
  

lmerTest::lmer(scale(starch.x) ~ scale(starch.y) + (1|HM), 
               redcap.data.nutrients.load) %>% summary()
# no association between energy intake and predicted microbial load


# :: Starch Intake == Starch CAZy ? ---------------------------------------

colnames(redcap.data.nutrients.load)[colnames(redcap.data.nutrients.load) == "starch.x"] = "starch.cazy"
colnames(redcap.data.nutrients.load)[colnames(redcap.data.nutrients.load) == "starch.y"] = "starch.intake"

lmerTest::lmer(scale(log10(starch.cazy)) ~ scale(starch.intake) + (1|HM), 
               redcap.data.nutrients.load) %>% summary()
lmerTest::lmer(scale(log10(starch.cazy)) ~ oars.on.rs + (1|HM), 
               redcap.data.nutrients.load) %>% summary()
lmerTest::lmer(scale(log10(starch.cazy)) ~ scale(starch.intake) + oars.on.rs + (1|HM), 
               redcap.data.nutrients.load) %>% summary()
# no association overall, or when controlling for RS supplement

# :: Mg==B.adolescentis? -------------------------------------------------------------------


mg.badol = oars.mgx.taxa.filt.1 %>% as.data.frame() %>%
  mutate(standard.name = rownames(.)) %>%
  dplyr::select(standard.name, Bifidobacterium_adolescentis) %>%
  merge(as.data.frame(redcap.data.nutrients %>%
                        mutate(standard.name = rownames(.)) %>%
                        dplyr::select(standard.name, magnesium, starch, energy)), by="standard.name") %>%
  merge(metadata.oars.stool.asv.double[,c("standard.name", "response", "HM")] %>% distinct(), by="standard.name") %>%
  mutate(mg.adj = magnesium / energy * 4.184 * 1000) # mg / kcal

mg.badol.plot = ggplot(mg.badol %>% subset(Bifidobacterium_adolescentis != 0), 
                       aes(x=log10(Bifidobacterium_adolescentis), y=magnesium))+
  geom_line(aes(group = HM), linewidth=0.2, linetype=2, alpha=0.7)+
  geom_point( fill="white", shape=21, size=3)+
  geom_smooth(method="lm", color="black")+
  ggpubr::stat_cor(method="spearman")+
  theme_classic() + theme(legend.position="none")+
  labs(x = "Log2(B. adolescentis)",
       y = "Magnesium intake (mg / d)")
mg.badol.plot
# no correlation; must control for individual identity with random effect

lmerTest::lmer(scale(magnesium) ~ scale(log10(Bifidobacterium_adolescentis)) + (1|HM), 
               mg.badol %>% subset(Bifidobacterium_adolescentis != 0)) %>% summary()
# significant
lmerTest::lmer(scale(log10(Bifidobacterium_adolescentis)) ~ scale(magnesium) + (1|HM), 
               mg.badol %>% subset(Bifidobacterium_adolescentis != 0)) %>% summary()

lmerTest::lmer(scale(magnesium) ~ response + (1|HM), 
               mg.badol %>% subset(Bifidobacterium_adolescentis != 0)) %>% summary()
# no association with Mg and response

# note: association is lost when controlling for energy intake
lmerTest::lmer(scale(log10(Bifidobacterium_adolescentis)) ~ scale(magnesium) + log10(energy) + (1|HM), 
               mg.badol %>% subset(Bifidobacterium_adolescentis != 0)) %>% summary()


# :: B vitamins -----------------------------------------------------------


b.vitamins = oars.mpx.kegg.mat.filt.1 %>% as.data.frame() %>%
  mutate(standard.name = rownames(.)) %>%
  dplyr::select(standard.name, `Nicotinate and nicotinamide metabolism`, `Pantothenate and CoA biosynthesis`) %>%
  merge(as.data.frame(redcap.data.nutrients %>%
                        mutate(standard.name = rownames(.)) %>%
                        dplyr::select(standard.name, thiamin, riboflavin, niacin, niacin_equivalents)), by="standard.name") %>%
  merge(metadata.oars.stool.asv.double[,c("standard.name", "response", "HM")] %>% distinct(), by="standard.name")

Hmisc::rcorr(b.vitamins[,!colnames(b.vitamins) %in% c("response" ,"HM", "standard.name")] %>% as.matrix())$P %>%
  reshape2::melt() %>%
  subset(value < 0.05)

b.vitamins.plot = ggplot(b.vitamins,
                       aes(x=log10(`Pantothenate and CoA biosynthesis`), y=thiamin))+
  geom_line(aes(group = HM), linewidth=0.2, linetype=2, alpha=0.7)+
  geom_point( fill="white", shape=21, size=3)+
  geom_smooth(method="lm", color="black")+
  #ggpubr::stat_cor(method="spearman")+
  theme_classic() + theme(legend.position="none")+
  labs(x = "Log2(Pantothenate and CoA biosynthesis)",
       y = "Thiamin intake (mg / d)")
b.vitamins.plot
# nope; but, FFQ is not validated to estimate B vitamins (it is for Magnesium)

lmerTest::lmer(scale(log10(`Pantothenate and CoA biosynthesis`)) ~ scale(thiamin) + (1|HM), 
               b.vitamins) %>% summary()
# significant; lower Thiamin predicts Pantothante biosynthesis


lmerTest::lmer(scale(thiamin) ~ response + (1|HM), 
               b.vitamins) %>% summary()



# :: Logistic regression --------------------------------------------------

# Goal: evaluate predictive power of Starch-active CAZy, B. adolescentis / R. bromii, and Butyrogens

oars.stool.baseline.predictors = 
  metadata.oars.stool.asv.double

oars.stool.baseline.predictors = merge(oars.stool.baseline.predictors,
                                       oars.mgx.prep.ml[,c("Ruminococcus_bromii", "Bifidobacterium_adolescentis",
                                                           "HM", "phase", "reltiming")],
                                       by =c( "HM", "phase", "reltiming") )
oars.stool.baseline.predictors = oars.stool.baseline.predictors %>%
  mutate(bin.response = (as.numeric(as.factor(response))-2) * -1,
         primary.d = log2(Bifidobacterium_adolescentis+Ruminococcus_bromii)) %>%
  mutate(Bifidobacterium_adolescentis = log2(Bifidobacterium_adolescentis),
         Ruminococcus_bromii = log2(Ruminococcus_bromii))


# Starch-active CAZy
oars.logreg.others = do.call(rbind, lapply(c("richness", "shannon", "but.i", "but.ii", "fd", 
                                             "adj.fiber", "fcal", "stool_water_perc", "load.asv", 
                                             "starch", "mucin", "starch.mucin", "primary.d",
                                             "Ruminococcus_bromii", "Bifidobacterium_adolescentis"), function(x){
  
  data.subset = oars.stool.baseline.predictors %>%
    subset(reltiming == "pre") %>%
    mutate(bin.response = (as.numeric(as.factor(response))-2) * -1) %>%
    dplyr::select(bin.response, x, response)
  
  colnames(data.subset)[2] = "variable"
  
  # needs log2 transforming
  if(x %in% c("richness", "but.i", "but.ii", "fcal", "starch", "mucin")){
    oars.logreg.starch = glm(bin.response ~ scale(log2(variable)),
                             family=binomial(link='logit'),
                             data=data.subset) %>%
                               broom::tidy(conf.int = T) %>%
                               filter(term != "(Intercept)") %>%
                               mutate(OR = 2^(estimate),
                                      conf.low = 2^(conf.low),
                                      conf.high = 2^(conf.high)) %>%
                               mutate(term = x)
  } 
  # already log2 scale
  if(x %in% c("shannon", "starch.mucin", "primary.d", "Ruminococcus_bromii", "Bifidobacterium_adolescentis")){
    oars.logreg.starch = glm(bin.response ~ scale(variable),
                             family=binomial(link='logit'),
                             data=data.subset) %>%
      broom::tidy(conf.int = T) %>%
      filter(term != "(Intercept)") %>%
      mutate(OR = 2^(estimate),
             conf.low = 2^(conf.low),
             conf.high = 2^(conf.high)) %>%
      mutate(term = x)
  }
  # already log10 scale, convert to log2
  if(x %in% c("load.asv", "load.mgx")){
    oars.logreg.starch = glm(bin.response ~ scale(log2(10^(variable))),
                             family=binomial(link='logit'),
                             data=data.subset) %>%
      broom::tidy(conf.int = T) %>%
      filter(term != "(Intercept)") %>%
      mutate(OR = 2^(estimate),
             conf.low = 2^(conf.low),
             conf.high = 2^(conf.high)) %>%
      mutate(term = x)
  } # not to be log-transformed at all
    else{
    oars.logreg.starch = glm(bin.response ~ scale(variable),
                             family=binomial(link = "logit"),
                             data=data.subset) %>%
                               broom::tidy(conf.int = T) %>%
                               filter(term != "(Intercept)") %>%
                               mutate(OR = exp(estimate),
                                      conf.low = exp(conf.low),
                                      conf.high = exp(conf.high)) %>%
                               mutate(term = x)
                             
    }
  # Also, do wilcox
  wilcox.results = wilcox.test(subset(data.subset, response == "high")$variable,
                               subset(data.subset, response == "low")$variable)
  
  oars.logreg.starch$wilcox.p = wilcox.results$p.value
  
  return(oars.logreg.starch)
  
}))
# rename variables
oars.logreg.others$term = c("Richness", "Shannon", "Butyrogens", "Butyrogens (Kircher)",
                            "Functional redundancy", "Adjusted fiber intake", "Fecal calprotectin",
                            "Stool moisture", "Microbial load", "Starch-active proteins", "Mucin-active proteins",
                            "Starch:Mucin ratio", "Primary degraders", "R. bromii", "B. adolescentis")


# merge
oars.logreg.forest.plot = oars.logreg.others %>%
  data.frame() %>%
  mutate(term = factor(term, levels=rev(c("Richness", "Shannon", "Butyrogens", "Butyrogens (Kircher)",
                                      "Functional redundancy", "Adjusted fiber intake", "Fecal calprotectin",
                                      "Stool moisture", "Microbial load", "Starch-active proteins", "Mucin-active proteins",
                                      "Starch:Mucin ratio", "Primary degraders", "B. adolescentis", "R. bromii")))) %>%
  ggplot(aes(y = term, x = OR)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  geom_text(aes(label = paste("(", round(p.value, digits=2), ")", sep=""), x = conf.high * 1.2), hjust = 0, size = 3) +
  scale_x_continuous(trans = "log10", breaks = c(0.1, 0.5, 1, 2, 5, 10, 20), 
                     limits = c(0.1, 20)) +
  labs(x = "Odds Ratio of Strong Response", y = "")+
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
oars.logreg.forest.plot

# how do OR p values compare to p values from Wilcox (requires flipping X and Y)
ggplot(oars.logreg.others,
       aes(x=p.value, y=wilcox.p))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_classic()
cor.test(oars.logreg.others$p.value,
         oars.logreg.others$wilcox.p)
# sig correlate

# :: Zn-dependent Peptidase -----------------------------------------------

# are Zn-dependent peptidase and carboxypeptidase anti-correlated?
cor.test(oars.mpx.cog.baseline$`Zn-dependent peptidase, M16 (insulinase) family`,
         oars.mpx.cog.baseline$`Zn-dependent carboxypeptidase, M32 family`,
         method="spearman")
# ρ = -0.70, p = 0.001181

# :: Predictors Change? ---------------------------------------------------

# goal: plot volcano plot of predictive features as they change over intervention

oars.predictors.lmer = rbind(oars.asv.data.glom.lmer %>% mutate(feature = taxa)%>% dplyr::select(feature, estimate, padj),
                             oars.mgx.lmer%>% mutate(feature = taxa)%>% dplyr::select(feature, estimate, padj),
                             oars.mpx.kegg.lmer %>% mutate(feature = taxa)%>% dplyr::select(feature, estimate, padj),
                             oars.mpx.cog.lmer %>% mutate(feature = taxa)%>% dplyr::select(feature, estimate, padj),
                             oars.mpx.cazy.lmer %>% mutate(feature = taxa)%>% dplyr::select(feature, estimate, padj),
                             oars.mbx.lmer%>% dplyr::select(feature, estimate, padj)) %>% as.data.frame() %>%
  subset(feature %in% oars.rf.loocv.imp.ttest$feature)
oars.predictors.lmer$clean.feature = ifelse(nchar(as.character(oars.predictors.lmer$feature))>30, paste(substr(oars.predictors.lmer$feature, 1, 30), "...", sep=""), as.character(oars.predictors.lmer$feature))

# volcano plot
oars.predictors.volcano = ggplot(oars.predictors.lmer,
                                         aes(x=estimate, y=padj))+
  geom_point(shape=21, aes(fill=estimate), size=2)+
  geom_hline(yintercept=(0.20), linetype=2, alpha=0.5)+
  geom_vline(xintercept=c(-1,1), linetype=2, alpha=0.5)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.05,0.1, 0.20, 0.5, 1))+
  ggrepel::geom_text_repel(aes(label=clean.feature), size=3)+
  scale_fill_gradient2(low="blue", high="red")+
  facet_wrap(~"Interaction: Predictors")+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10))+
  labs(x="Interaction Coefficient", y="FDR")
oars.predictors.volcano



# :: Predictors Network? --------------------------------------------------

# baseline correlations

# spearman correlation (takes a min)
metadata.oars.predictors.cor = Hmisc::rcorr(as.matrix(oars.all.omics.baseline[,!colnames(oars.all.omics.baseline)%in%c("response", "HM", "phase", "diagnosis", "adj.fiber")]), type="spearman")

metadata.oars.predictors.cor.df = 
  data.frame(reshape2::melt(metadata.oars.predictors.cor$r)) %>%
  mutate(pval = reshape2::melt(metadata.oars.predictors.cor$P)$value) %>%
  subset(!is.na(pval)) %>%
  arrange(pval)

metadata.oars.predictors.cor.df = subset(metadata.oars.predictors.cor.df, Var1 %in% oars.rf.loocv.imp.df$feature) %>%
  mutate(Var2 = as.character(Var2))%>%
  mutate(padj = p.adjust(pval, method="BH")) %>%
  arrange(padj)

metadata.oars.predictors.cor.df$clean.feature = ifelse(nchar(as.character(metadata.oars.predictors.cor.df$Var1))>30, paste(substr(metadata.oars.predictors.cor.df$Var1, 1, 30), "...", sep=""), as.character(metadata.oars.predictors.cor.df$Var1))
metadata.oars.predictors.cor.df$clean.var2 = ifelse(nchar(as.character(metadata.oars.predictors.cor.df$Var2))>30, paste(substr(metadata.oars.predictors.cor.df$Var2, 1, 30), "...", sep=""), as.character(metadata.oars.predictors.cor.df$Var2))

metadata.oars.predictors.cor.df[,c("clean.feature", "clean.var2", "value", "padj")] %>%
  subset(padj < 0.20)


subset(metadata.oars.predictors.cor.df, clean.feature == "Ketopantoate hydroxymethyltran...")
# associated with:
# Zn-dependent carboxypeptidase (R = 0.80, FDR < 0.20)
# Zn-dependent peptidase (R = -0.77, FDR < 0.20)
# Phosphotransacetylase (R = -0.75, FDR < 0.20)
# Na+translocating ferredoxin:NAD+ oxidoreductase (R = -0.74, FDR < 0.20)
# R = -0.60 with Panthothenate synthetase (nom sig, but did not survive FDR)

subset(metadata.oars.predictors.cor.df, clean.feature == "Heterodisulfide reductase, sub...")
# associated with:
# grapefruit (R=0.79, FDR < 0.20)
# other Heterodisulfide reductase subunit (R=77, FDR < 0.20),
# Na+-translocating ferredoxin:NAD+ oxidoreductase (R=0.77, FDR < 0.20)
# Phosphate-selective porin (R=-0.75, FDR < 0.20)
# Prephenate dehydrogenase (R=0.75, FDR < 0.20) (intermediate in phenylalanine and tyrosine metabolism)
# Tryptophan synthase (R = 0.74, FDR < 0.20)

subset(metadata.oars.predictors.cor.df, clean.feature == "Zn-dependent carboxypeptidase,...")
# associated with:
# Ketopantoate hydroxymethyltranferase (R = 0.80, FDR < 0.20)
# Angelakisella massiliensis (ASV) (R = -0.80, FDR < 0.20)
# NADH:ubiquinone oxidoreductase (R = 0.79, FDR < 0.20)
# ABC-type sugar transport system (R = -0.76, FDR < 0.20)
# Phosphate-selective porin (R = 0.75, FDR < 0.20)

# Note: I targeted the correlation between Zn peptidases, so did not run FDR (p = 0.001)
# i.e. not exploratory

subset(metadata.oars.predictors.cor.df, clean.feature == "Ribosome biogenesis protein, N...")
# associated with:
# 2-succinyl-6-hydroxy... (R = 0.88, FDR < 0.01)
# Eubacterium_rectale (MGX) (R = 0.85, FDR < 0.05)
# Agathobacter_rectale (ASV) (R = 0.82, FDR < 0.10)
# and several others

# since this as exploratory, I used FDR

# These are correlations among baseline stool values (baseline = 0M or 3M)
# Unlike the delta correlations in section 5

# :: Methanogenesis -------------------------------------------------------

oars.mpx.kegg.prep[,c("Methane metabolism", "standard.name" ,"HM", "timing", "phase", "response")]
oars.mpx.cog.prep[,c("Heterodisulfide reductase, subunit C", "standard.name" ,"HM", "timing", "phase", "response")]
oars.mgx.prep[,c("Methanobrevibacter_smithii", "standard.name" ,"HM", "timing", "phase", "response")]
# M. smithii was only detected in 1 HM

oars.methanogenesis.data = 
  merge(metadata.oars.stool.double,
        oars.mpx.kegg.prep[,c("Methane metabolism", "HM", "phase", "reltiming")], by=c("HM", "phase", "reltiming"))
oars.methanogenesis.data = 
  merge(oars.methanogenesis.data,
        oars.mpx.cog.prep[,c("Heterodisulfide reductase, subunit C", "HM", "phase", "reltiming")], by=c("HM", "phase", "reltiming"))

lmerTest::lmer(log2(`Methane metabolism`) ~ reltiming*response + phase + (1|HM),
               data=(oars.methanogenesis.data)) %>%
  summary()

lmerTest::lmer(log2(`Heterodisulfide reductase, subunit C`) ~ reltiming*response + phase + (1|HM),
               data=(oars.methanogenesis.data)) %>%
  summary()

# no evidence of a meaningful impact on methanogenesis

# how about methane metabolism as a predictor
ggplot(oars.mpx.kegg.baseline, aes(x=response, y=log2(`Methane metabolism`))) + 
  geom_point()
wilcox.test(subset(oars.mpx.kegg.baseline, response == "low")$`Methane metabolism`,
            subset(oars.mpx.kegg.baseline, response == "high")$`Methane metabolism`)
# no association at all


oars.erectale.data = 
  merge(metadata.oars.stool.double,
        oars.mgx.prep[,c("Eubacterium_rectale", "HM", "phase", "reltiming")], by=c("HM", "phase", "reltiming"))
lmerTest::lmer(log2(`Eubacterium_rectale`) ~ reltiming*response + phase + (1|HM),
               data=(oars.erectale.data)) %>%
  summary()
ggplot(oars.mgx.baseline, aes(x=response, y=log2(`Eubacterium_rectale`))) + 
  geom_point()
wilcox.test(subset(oars.mgx.baseline, response == "low")$`Eubacterium_rectale`,
            subset(oars.mgx.baseline, response == "high")$`Eubacterium_rectale`)

# therefore, don't bother interpreting Heterodisulfide reductase as meaningfully related to methanogenesis

# :: PubChemR -------------------------------------------------------------

BiocManager::install("PubChemR")


chemical.names = gsub(" \\|.*", "", colnames(oars.mbx.raw.mat)) %>% unique()

pcr.query = do.call(rbind, lapply(chemical.names, function(x){
  print(x)
  # query database
  chemicals_by_name <- PubChemR::get_aids(
    identifier = x,
    namespace = "name",
    domain = "compound"
  )
  # if success, record CID
  if(!is.null(chemicals_by_name$result[[1]])){
    if(chemicals_by_name$result[[1]]$success == T){
    cid = chemicals_by_name$result[[1]]$result$InformationList$Information[[1]][1]
    }
    } else {
    cid = NA
  }
  data.frame(compound = x,
             CID = cid)
}))

sum(!is.na(pcr.query$CID)) / nrow(pcr.query)
# 47% of annotated metabolites can be assigned CID

# not very useful?


# :: ASV RF Importances -------------------------------------------------------


oars.loocv.models.importances.asv = rf.function(data.types = "ASV",
                                            iters = 15,
                                            output = "importances")

oars.rf.loocv.imp.asv.df = oars.loocv.models.importances.asv %>%
  #subset(imp != 0) %>%
  group_by(feature) %>%
  mutate(n.present = n()) %>%
  mutate(mean.imp = mean(na.omit(imp))) %>%
  mutate(imp.low = mean(na.omit(imp)) - (sd(na.omit(imp))/sqrt(n()) * 1.96)) %>% # 95% CI
  mutate(imp.high = mean(na.omit(imp))+ (sd(na.omit(imp))/sqrt(n()) * 1.96)) %>% # 95% CI
  subset(mean.imp != 0) %>%
  dplyr::select(feature, mean.imp, imp.low, imp.high, n.present) %>% distinct() %>%
  arrange(-mean.imp) %>% data.frame()

# log2fc
oars.rf.loocv.imp.asv.df = oars.rf.loocv.imp.asv.df %>% 
  subset(!is.na(mean.imp))%>%
  # select the top 10 mean importance
  arrange(mean.imp) %>% slice_max(mean.imp, n=10) 

# GreenGenes2 is incompatible with ML project

# >>> FIGURES ------------------------------------------------------

# :: Main -----------------------------------------------------------------

# "Ex vivo fermentation predicts biochemical inflammation and microbiome response to personalized resistant starch 
# in children with medically stable inflammatory bowel disease"

# Table 1: Patient characteristics

# Figure 1 = Stool collections
metadata.oars.stool.plot %>%
  ggsave(filename="./oars_plots/oars_1_stools.pdf",
         width=6.5, height=3.5)

# Figure 2 = Fiber intake: oars_1_ffq_fiber_rs.pdf
# Note: likely combine this with Fig 1

# Figure 3A = RS algorithm (stitched from 2025_06_09)

# Figure 3B = RS selection: oars_1_rs_selections.pdf

# Figure 4 = Butyrogens + Fecal cal(group level)
(metadata.oars.but.i.plot + metadata.oars.fcal.plot) %>%
  ggsave(filename = "./oars_plots/oars_1_butyrogens_fecalcal.pdf",
         width=10, height=5, device = cairo_pdf)

# Figure 5: Other stool variables
# Note: likely combine with Fig 4
((metadata.oars.but.i.plot + metadata.oars.fcal.plot) / 
  (metadata.oars.stool.asv.richness.plot+
    metadata.oars.stool.asv.shannon.plot+
    metadata.oars.stool.asv.fd.plot+
    metadata.oars.stool.asv.beta.between.plot+
    metadata.oars.water.plot+
    metadata.oars.load.asv.plot+
    metadata.oars.stool.starch.plot+
    metadata.oars.stool.mucin.plot+
    metadata.oars.stool.starch.mucin.plot)+patchwork::plot_layout(nrow=2, heights=c(1,2))) %>%
  ggsave(filename="./oars_plots/2026_01_15_oars_2_group_level.pdf",
         width=11, height=11, device = cairo_pdf)

# Figure 6 = Select taxa LFC heatmaps
ggsave(pheatmap::pheatmap(t(oars.asv.lfc.treatment.prev),
                          color=colorRampPalette(c("blue","white", "red"))(100),
                          clustering_distance_rows = "correlation",
                          clustering_distance_cols = "correlation",
                          breaks=c(seq(min(oars.asv.lfc.treatment.prev), 0, length.out=ceiling(100/2) + 1), 
                                   seq(max(oars.asv.lfc.treatment.prev)/100, max(oars.asv.lfc.treatment.prev), length.out=floor(100/2))),
                          annotation_col=oars.asv.lfc.treatment.prev.mapping,
                          annotation_colors = oars.asv.lfc.treatment.prev.mapping.colors,
                          border_color = "white"),
       filename = "./oars_plots/2026_01_15_oars_1_asv_heatmap_1.pdf",
       width=14, height=7, device = cairo_pdf)


# Figure 7A: Omics PCoA/PCA
((((oars.asv.pcoa.plot+
    oars.mgx.pcoa.plot+
    oars.mpx.kegg.pca.plot+
    oars.mpx.cog.pca.plot+
    oars.mpx.cazy.pca.plot+
    oars.mbx.pca.plot)+
    patchwork::plot_layout(nrow=3))|
  oars.pca.permanova.plot)+patchwork::plot_layout(widths=c(2,1)))%>%
  ggsave(filename="./oars_plots/2026_01_15_oars_2_omics_pca.pdf",
         width=9, height=8)

# Figure 7B: Omics LFC 
rownames(oars.group.mgx.mpx.mbx.lfc.map.map) = gsub(" \\|.*", "", rownames(oars.group.mgx.mpx.mbx.lfc.map.map))
pheatmap::pheatmap(reshape2::acast(oars.group.mgx.mpx.mbx.lfc %>% mutate(sample2 = paste(HM, timing)) %>%
                                     mutate(variable = gsub(" \\|.*", "", variable)),
                                   variable ~ sample2, value.var="slfc"),
                   color=colorRampPalette(c("blue","white", "red"))(100),
                   border="white",
                   #angle_col = 45,
                   #clustering_distance_rows = "correlation",
                   #clustering_distance_cols = "correlation",
                   fontsize_row = 9,
                   fontsize_col = 7,
                   annotation_row = oars.group.mgx.mpx.mbx.lfc.map.map,
                   annotation_col = oars.omics.heatmap.mapping[,colnames(oars.omics.heatmap.mapping)!="Patient"],
                   annotation_colors = list(Coefficient = colorRampPalette(c("blue","white", "red"))(100),
                                            Datatype = omics.colors,
                                            Timing = c("3M" = "white", "6M" = "lightgrey"),
                                            RS = rs.colors),
                   breaks=c(seq(min(oars.group.mgx.mpx.mbx.lfc$slfc), 0, length.out=ceiling(100/2) + 1), 
                            seq(max(oars.group.mgx.mpx.mbx.lfc$slfc)/100, max(oars.group.mgx.mpx.mbx.lfc$slfc), length.out=floor(100/2))))%>%
  ggsave(filename="./oars_plots/2026_01_15_oars_2_omics_heatmap.pdf",
         width=8, height=6)

# Response Tree (old)

# Respone Tree (new)
cowplot::plot_grid(stats.tree.features.plots.top, stats.tree.features.plots,
                   nrow=2, rel_heights=c(3,1))%>%
  ggsave(filename="./oars_plots/2026_01_15_oars_2_starch_fermentation_stream.pdf",
         width=15, height=11.5, device=cairo_pdf)


# Figure 8: Responders vs Non-responders
((oars.ph.rs.boxplot+
    oars.ph.bidist.plot.2+patchwork::plot_spacer()+oars.ph.score.change.plot+theme(legend.position="none"))+
  patchwork::plot_layout(widths=c(2,1,1,1))) %>%
  ggsave(filename="./oars_plots/2026_01_15_oars_2_ph_response_groups.pdf",
         width=12, height=4, device = cairo_pdf)

# Figure 9: Butyrogens, Fecal calprotectin, others
# Note: likely combine with Fig 8
(((metadata.oars.stool.asv.double.but.i.plot+metadata.oars.fcal.double.plot)/
    ((metadata.oars.stool.asv.double.richness.plot+
        metadata.oars.stool.asv.double.shannon.plot+
        metadata.oars.stool.asv.double.fd.plot+
        metadata.oars.stool.asv.double.beta.between.plot+
        metadata.oars.water.double.plot+
        metadata.oars.load.asv.double.plot+
        metadata.oars.stool.asv.double.starch.plot+
        metadata.oars.stool.asv.double.mucin.plot+
        metadata.oars.stool.asv.double.starch.mucin.plot)+
        patchwork::plot_layout(nrow=3)))+
  patchwork::plot_layout(heights=c(1,2))) %>%
     ggsave(filename="./oars_plots/2026_01_15_oars_2_ph_responders_nonresponders.pdf",
            width=12, height=10,device = cairo_pdf)

# Figure 10A: Omic Cor w Fecalcal
(((((oars.fecalcal.omics.cor.plot|plsda_model_plot)/
  oars.fecalcal.omics.cor.plots)+patchwork::plot_layout(heights=c(2,2.5))))+
  patchwork::plot_layout(widths=c(3,2))) %>%
  ggsave(filename="./oars_plots/2026_01_15_oars_2_fecalcal_correlates.pdf",
         width=9, height=8,device = cairo_pdf)

# Figure 10B: Omic LFC Heatmaps

pheatmap::pheatmap(reshape2::acast(subset(oars.asv.mgx.mpx.mbx.lfc, feature %in% rownames(oars.asv.mgx.mpx.mbx.lfc.map.vimp)),
                                   feature ~ sample, value.var="slfc"),
                   color=colorRampPalette(c("blue","white", "red"))(100),
                   border="white",
                   #angle_col = 45,
                   clustering_distance_rows = "correlation",
                   clustering_distance_cols = "correlation",
                   fontsize_row = 7,
                   fontsize_col = 7,
                   annotation_col = oars.omics.lfc.map %>% dplyr::select(-sample),
                   annotation_colors = list(Response = c(`Strong Response` = gg_color_hue(2)[1],
                                                         `Weak Response` = gg_color_hue(2)[2]),
                                            #`Coefficient` = colorRampPalette(c("blue", "red"))(2),
                                            #Coefficient = c(`Positive` = "red",
                                            #                `Negative` = "blue"),
                                            `Fcal Correlation` = colorRampPalette(c("blue","white", "red"))(100),
                                            `VIP` = colorRampPalette(c("white", "red"))(100),
                                            `Data type` = c(ASV = "#8DD3C7",
                                                            Species = "#FFFFB3",
                                                            Pathway = "#BEBADA",
                                                            COG = "#FB8072",
                                                            CAZy = "#80B1D3",
                                                            Metabolite = "#FDB462")),
                   annotation_row = oars.asv.mgx.mpx.mbx.lfc.map.vimp,
                   breaks=c(seq(min(na.omit(oars.asv.mgx.mpx.mbx.lfc$slfc)), 0, length.out=ceiling(100/2) + 1), 
                            seq(max(na.omit(oars.asv.mgx.mpx.mbx.lfc$slfc))/100, max(na.omit(oars.asv.mgx.mpx.mbx.lfc$slfc)), length.out=floor(100/2))))%>%
  ggsave(filename="./oars_plots/2026_01_15_oars_2_responders_fcal_heatmap.pdf",
         width=10, height=5)

# ?
(oars.badol.starch.fcal.plot +
    oars.badol.simple.plot+
    patchwork::plot_layout(heights=c(4,1))) %>%
  ggsave(filename="./oars_plots/2026_01_15_oars_2_responders_starch_responses.pdf",
         width=5, height=10)

# Figure 11: Machine learning
((patchwork::free(p1, type = "label")| 
    patchwork::free(p2, type = "label") | 
    patchwork::free(p3, type = "label")) + 
    patchwork::plot_layout(nrow=1, widths=c(1,2,1))) %>%
  ggsave(filename="./oars_plots/oars_2_ph_responders_predictors.pdf",
         width=22, height=8, device = cairo_pdf)

# 8 total figures


# :: Feature Space --------------------------------------------------------


# Group-level (compliant --> 10% prev or 80% MBX per time phase pre-on-post)
dim(oars.asv.data.median.1.filt) # ASV Bray: 1165 in 50 samples
dim(oars.asv.data.glom) # Maaslin: 1554 in 66 samples
dim(oars.mgx.taxa.filt.1) # MGX Bray + Maaslin2: 120 in 48 samples
dim(oars.mpx.kegg.mat.filt.1) # KEGG PCA + Maaslin2: 180 in 48 samples
dim(oars.mpx.cog.mat.filt.1) # COG PCA + Maaslin2: 2481 in 48 samples
dim(oars.mpx.cazy.mat.filt.1) # CAZy PCA + Maaslin2: 66 in 48 samples
dim(oars.mbx.raw.mat.filt.1) # MBX PCA + Maaslin2: 203 in 46 samples (annotated)

# Response Tree (compliant --> 10% prev or 80% MBX across all)
dim(oars.asv.data.glom.filt.3) # ASV: 341 in 50 samples (compliant; filtered in Maaslin2)
dim(oars.mgx.taxa.filt.3) # MGX: 542 in 48 samples (compliant)
dim(oars.mpx.kegg.filt.3) # KEGG: 180 in 48 samples (compliant)
dim(oars.mpx.cog.filt.3) # COG: 2481 in 48 samples (compliant)
dim(oars.mpx.cazy.filt.3) # CAZy: 66 in 48 samples (compliant)
dim(oars.mbx.filt.3) # MBX: 153 in 46 samples (annotated)

# Fecal cal (all compliant --> 10% prev or 80% MBX across all)
dim(oars.asv.data.glom.filt.2) # ASV: 340 in 66 samples
dim(oars.mgx.taxa.filt.2) # MGX: 514 in 65 samples 
dim(oars.mpx.kegg.filt.2) # KEGG: 178 in 66 samples 
dim(oars.mpx.cog.filt.2) # COG: 2478 in 66 samples
dim(oars.mpx.cazy.filt.2) # CAZy: 65 in 66 samples
dim(oars.mbx.filt.2) # MBX: 150 in 70 samples (annotated)

# Machine learning (baseline compliant --> 20% prev or 80% MBX across all)
dim(oars.asv.data.glom.prep.ml) # ASV: 257 in 21 samples
dim(oars.mgx.prep.ml) # MGX: 350 in 21 samples
dim(oars.mpx.kegg.prep.ml) # KEGG: 180 in 19 samples
dim(oars.mpx.cog.prep.ml) # COG: 2352 in 19 samples
dim(oars.mpx.cazy.prep.ml) # CAZy: 72 in 19 samples
dim(oars.mbx.prep.ml) # MBX: 173 in 21 samples
dim(redcap.data.ffq.prep.ml) # FFQ: 311 in 21 samples


# :: Supplementals --------------------------------------------------------


# SCFA correlations: oars_1b_rapidaim_validation_correlations.pdf

# SCFA Breakpoints: oars_1b_rapidaim_ph_breakpoints.pdf

# MLP Validation: oars_1_mlp_asv_mgx.pdf + oars_1_mlp_retrain.pdf

# But ~ pH cors: 


# Butyrogen validation
(metadata.oars.but.ii.plot+
    metadata.oars.stool.asv.plot+
   kircher.butyrogens.plot)%>%
  ggsave(filename="./oars_plots/oars_1b_butyrogen_validation.pdf",
         width=12, height=3.5, device = cairo_pdf)

oars.hysteresis.plot%>%
  ggsave(filename="./oars_plots/oars_supp_hysteresis.pdf",
         width=8, height=3.5, device = cairo_pdf)



# Volcano plots - Group
(((oars.asv.maaslin.036.volcano/
  oars.asv.maaslin.6912.volcano)|
  (oars.mgx.maaslin.036.volcano/
  oars.mgx.maaslin.6912.volcano)|
  (oars.mpx.kegg.maaslin.036.volcano/
  oars.mpx.kegg.maaslin.6912.volcano)|
    (oars.mpx.cog.maaslin.036.volcano/
       oars.mpx.cog.maaslin.6912.volcano)|
  (oars.mpx.cazy.maaslin.036.volcano/
  oars.mpx.cazy.maaslin.6912.volcano)|
    (oars.mbx.maaslin.036.volcano/
       oars.mbx.maaslin.6912.volcano))+
    patchwork::plot_layout(nrow=2)) %>%
  ggsave(filename="./oars_plots/2026_01_15_oars_supp_omics_volcanos.pdf",
         width=18, height=18, device = cairo_pdf)
  
# Fcal correlates network
metadata.oars.stool.omics.cor.network.pcoa %>%
  ggsave(filename="./oars_plots/2026_01_15_oars_supp_omics_fcal_network.pdf",
         width=18, height=12, device = cairo_pdf)

# Mg == B. adolescentis
mg.badol.plot %>%
  ggsave(filename="./oars_plots/oars_supp_mg_badol.pdf",
         width=4, height=3, device = cairo_pdf)

# Volcano plots - Responders
(oars.asv.data.glom.lmer.volcano+
  oars.mgx.lmer.volcano+
  oars.mpx.cog.lmer.volcano+
  oars.mpx.kegg.lmer.volcano+
  oars.mpx.cazy.lmer.volcano+
  oars.mbx.lmer.volcano+patchwork::plot_layout(nrow=2))%>%
  ggsave(filename="./oars_plots/2026_01_15_oars_supp_omics_resp_volcanos.pdf",
         width=18, height=12, device = cairo_pdf)


# Load vs Water
oars.stool.water.load.correlation.plot%>%
  ggsave(filename="./oars_plots/oars_2b_ph_load_vs_water.pdf",
         width=4, height=3, device = cairo_pdf)

oars.logreg.forest.plot%>%
  ggsave(filename="./oars_plots/oars_supp_predictors_forest_plot.pdf",
         width=4, height=4, device = cairo_pdf)

oars.predictors.volcano%>%
  ggsave(filename="./oars_plots/oars_supp_predictors_volcano.pdf",
         width=4, height=4, device = cairo_pdf)


# :: Graphical Abstract Figs --------------------------------------

# re-opt at ~ 78 days + 14 days + extra time before they started; say 100
subset(metadata.oars.stool, timing == "3M")$oars.days %>% mean()
oars.nutrient.data = read.csv("./2025_07_24_oars_fiber.csv")

# PART 1 A: RS intake
oars.1.ga.1 = ggplot(oars.nutrient.data,
                     aes(x=oars.days/32, 
                         y=total_rs))+
  annotate("rect", xmin=0, xmax=6, ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_vline(xintercept=0, linetype=2, linewidth=0.5, color="grey")+
  geom_vline(xintercept=3, linetype=2, linewidth=0.5, color="grey")+
  scale_x_continuous(breaks=c(-3, 0, 3, 6, 9, 12))+
  geom_smooth(color="black",  se=T, linewidth=1)+
  scale_fill_manual(values=c("grey",2,"grey"))+
  scale_alpha_manual(values=c(1,0.2))+
  scale_shape_manual(values=c(23,21,21))+
  facet_wrap(~"Resistant starch")+
  theme_classic()+theme(legend.position="none",
                        strip.text=element_text(size=12))+
  labs(x="Months since starting RS", y="Resistant starch (g/day)")
oars.1.ga.1

# PART 1 B: Butyrogens

oars.1.ga.2 = ggplot(subset(metadata.oars.stool.asv, compliant == TRUE),
                     aes(x=oars.days/32, 
                         y=but.i*100))+
  annotate("rect", xmin=0, xmax=6, ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_vline(xintercept=0, linetype=2, linewidth=0.5, color="grey")+
  geom_vline(xintercept=3, linetype=2, linewidth=0.5, color="grey")+
  scale_x_continuous(breaks=c(-3, 0, 3, 6, 9, 12))+
  geom_smooth(color="black")+
  # scale_y_log10()+
  #geom_text(aes(label=standard.name), size=3)+
  facet_wrap(~"Butyrogens")+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(hjust = 0.5, size=12))+
  labs(x="Months since starting RS", y="Butyrogens (%)")
oars.1.ga.2

# PART 1 C: Fecal cal
oars.1.ga.3 = ggplot(subset(metadata.oars.stool.asv, compliant == TRUE),
                     aes(x=oars.days/32, 
                         y=fcal))+
  annotate("rect", xmin=0, xmax=6, ymin=0, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_vline(xintercept=0, linetype=2, linewidth=0.5, color="grey")+
  geom_vline(xintercept=3, linetype=2, linewidth=0.5, color="grey")+
  scale_x_continuous(breaks=c(-3, 0, 3, 6, 9, 12))+
  geom_smooth(color="black")+
  scale_y_log10()+
  #geom_text(aes(label=standard.name), size=3)+
  facet_wrap(~"Calprotectin")+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(hjust = 0.5, size=12))+
  labs(x="Months since starting RS", y="Fecal calprotectin (μg/g)")
oars.1.ga.3


# PART 1 C: CAZymes
oars.1.ga.4 = ggplot(metadata.oars.stool.asv,
                     aes(x=oars.days/32, 
                         y=starch))+
  annotate("rect", xmin=0, xmax=6, ymin=-Inf, ymax=Inf, alpha=0.2, fill="salmon") +
  geom_vline(xintercept=0, linetype=2, linewidth=0.5, color="grey")+
  geom_vline(xintercept=3, linetype=2, linewidth=0.5, color="grey")+
  scale_x_continuous(breaks=c(-3, 0, 3, 6, 9, 12))+
  geom_smooth(color="black",  se=T, linewidth=1)+
  scale_fill_manual(values=c("grey",2,"grey"))+
  scale_alpha_manual(values=c(1,0.2))+
  scale_shape_manual(values=c(23,21,21))+
  facet_wrap(~"Starch CAZy")+
  theme_classic()+theme(legend.position="none",
                        strip.text=element_text(size=12))+
  labs(x="Months since starting RS", y="Σ Protein Intensity")
oars.1.ga.4

(oars.1.ga.1+oars.1.ga.4+oars.1.ga.2+oars.1.ga.3 + patchwork::plot_layout(nrow=2)) %>%
  ggsave(filename="./oars_plots/oars_1_graphical_abstract_plots.pdf",
         width=5, height=3.5, device = cairo_pdf)

# defunct // PART 1 D: Taxa LFC
set.seed(25)
pheatmap::pheatmap(t(oars.asv.lfc.treatment.prev[,sample(1:ncol(oars.asv.lfc.treatment.prev), size=10)]),
                   color=colorRampPalette(c("blue","white", "red"))(100),
                   breaks=c(seq(min(oars.asv.lfc.treatment.prev), 0, length.out=ceiling(100/2) + 1), 
                            seq(max(oars.asv.lfc.treatment.prev)/100, max(oars.asv.lfc.treatment.prev), length.out=floor(100/2))),
                   annotation_col=oars.asv.lfc.treatment.prev.mapping[,colnames(oars.asv.lfc.treatment.prev.mapping)!="Timing"],
                   annotation_colors = oars.asv.lfc.treatment.prev.mapping.colors,
                   show_rownames = F,show_colnames = F, legend=F,annotation_legend = F)


# PART 2 A: Responders Fcal
(ggplot(subset(metadata.oars.stool.double, !hm_phase %in% c("HM0924_rs2", "HM0932_rs2", "HM0759_rs2")) %>%
         mutate(reltiming = factor(ifelse(reltiming == "pre", "Pre RS", "Post RS"), levels=c("Pre RS", "Post RS"))) %>%
         mutate(Response = ifelse(response == "low", "Weak Response", "Strong Response")),
       aes(x=reltiming, y=fcal))+
  scale_y_log10()+
  geom_boxplot(width=0.5)+
  geom_line(aes(group=hm_phase), linetype=2, alpha=0.3, linewidth=0.5)+
  geom_point(shape=21, aes(fill=response), size=3)+
  geom_hline(yintercept=250, linetype=1, alpha=1, color="red")+
  #ggrepel::geom_label_repel(aes(label=RS_Name, fill=RS_Name),size=4)+
  geom_text(data = stats.fcal.ph, x=1.5, y=log10(max(na.omit(metadata.oars.stool.double$fcal))), 
            aes(label=ifelse(pval < 0.05, "p < 0.05", "n.s.")), vjust=1.2, size=3.5)+
  facet_wrap(~Response, nrow=1)+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=12),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="",
       y="Fecal Calprotectin (μg/g)")) %>%
  ggsave(filename="./oars_plots/oars_2_graphical_abstract_fcal.pdf",
         width=3.8, height=2.5, device = cairo_pdf)


# PART 2 B: ROC
(ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  # Individual ROC curves for each iteration
  #geom_line(data = oars.rf.loocv.roc %>% arrange(sens),
  #          aes(x=1-spec, y=sens, group=iter))+
  # Error ribbon (95% CI)
  geom_ribbon(data = roc_summary, 
              aes(x = fpr, ymin = lower, ymax = upper), 
              #fill = RColorBrewer::brewer.pal(n=5, "Set3")[1], alpha = 0.2) +
              fill = "black", alpha=0.2)+
  # Mean ROC curve
  geom_path(data = roc_summary, 
            aes(x = fpr, y = mean_sens), 
            #color = RColorBrewer::brewer.pal(n=5, "Set3")[1], size = 1) +
            color="black")+
  # add label
  annotate(geom="text", x=0.75, y=0.25,
           label=paste("AUC\n", round(oars.loocv.models.best.stats$mean.auc, digits=2)),
           size=4)+
  theme_classic()+theme(strip.text=element_text(size=10))+
  facet_wrap(~"Multi-Omic RandomForest")+
  xlim(0,1)+
  ylim(0,1)+
  labs(x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)")) %>%
  ggsave(filename="./oars_plots/oars_2_graphical_abstract_roc.pdf",
         width=2.5, height=2.2, device = cairo_pdf)




# >>> TABLES ---------------------------------------------------------------

# note: update figure ## once finalized

file.name = "2025_10_15_oars_supplemental_tables.xlsx"

openxlsx::write.xlsx(
  
  list(
    
    # group-level maaslin2 (on and off RS)
    "0X = maaslin asv 0 6" = oars.asv.maaslin.036,
    "0X = maaslin asv 6 12" = oars.asv.maaslin.6912,
    "0X = maaslin mgx 0 6" = oars.mgx.maaslin.036,
    "0X = maaslin mgx 6 12" = oars.mgx.maaslin.6912,
    "0X = maaslin kegg 0 6" = oars.mpx.kegg.maaslin.036,
    "0X = maaslin kegg 6 12" = oars.mpx.kegg.maaslin.6912,
    "0X = maaslin cog 0 6" = oars.mpx.cog.maaslin.036,
    "0X = maaslin cog 6 12" = oars.mpx.cog.maaslin.6912,
    "0X = maaslin cazy 0 6" = oars.mpx.cazy.maaslin.036,
    "0X = maaslin cazy 6 12" = oars.mpx.cazy.maaslin.6912,
    "0X = maaslin mbx 0 6" = oars.mbx.maaslin.036,
    "0X = maaslin mbx 6 12" = oars.mbx.maaslin.6912,
    
    # group-level interaction (strong vs weak)
    "0X = interactions asv" = oars.asv.data.glom.lmer,
    "0X = interactions mgx" = oars.mgx.lmer,
    "0X = interactions kegg" = oars.mpx.kegg.lmer,
    "0X = interactions cog" = oars.mpx.cog.lmer,
    "0X = interactions cazy" = oars.mpx.cazy.lmer,
    "0X = interactions mbx" = oars.mbx.lmer,
    
    # response tree
    "0X = conserved variable" = oars.double.response.tree.lm,
    
    # correlations (lfc)
    "0X = logfc correlations" = all.features.cor.all,
    
    # LOOCV all models
    "0X = ml aucs" = oars.loocv.models.ml.omics,
    
    # feature importances
    "0X = importances" = oars.loocv.models.importances,
    
    # top feature wilcox
    "0X = important wilcox" = oars.rf.loocv.imp.ttest
  ),
  file = file.name
)

# large file > 200 MB, takes ~2 min to save

# :: ----------------------------------------------------------------------


# :: -----------------------------------------------------------------------


# :: Weak Fermenters ------------------------------------------------------

# Do weak fermenters "strongly respond" to ANY of the RS?


# process
oars.rapidaim.scores <- readRDS("./2025_06_09_oars_scores.Rds") %>%
  # calculate LFC to PBS
  group_by(HM, timing) %>%
   # also calculate delta.pH
  mutate(delta.ph = med.ph - med.ph[RS_Name == "PBS"]) %>%
   # subset to 0M and 3M
  subset(timing %in% c("0M", "3M")) %>%
  # indicate response
  mutate(response = ifelse(delta.ph < -1.27, "strong", "weak"))
# good

# fix 618
# subset to just the selected RS
oars.rapidaim.scores = oars.rapidaim.scores %>%
  mutate(selected = ifelse(HM == "HM0618" & timing %in% c("0M", "3M"), "ActistarRT", selected))
  
subset(oars.rapidaim.scores, HM == "HM0874" & timing == "0M")

oars.rapidaim.scores %>%
  group_by(HM, timing) %>%
  filter(delta.ph == min(delta.ph)) %>%
  dplyr::select(HM, timing, RS_Name, delta.ph, response) %>%
  subset(response == "weak")
# 6 remain weak
# 22 are strong

# compared to
oars.rapidaim.scores %>%
  group_by(HM, timing) %>%
  filter(selected == RS_Name) %>%
  dplyr::select(HM, timing, RS_Name, delta.ph, response) %>%
  subset(response == "weak")  
# 16 weak
# 12 strong

# :: ----------------------------------------------------------------------



# >> PREDICT VARIANCE ------------------------------------------------------------------

# Question: Can baseline microbiome composition predict variance in response?

# Goal: Use baseline microbiome features to predict change in omics space (PC1 and PC2)

# Approach:
# 1. LFC data --> PCA (extract PC1 and PC2)
# 2. Baseline data dimension reduction
# 3. Procrustes Baseline vs LFC
# 4. Repeat for each combination (omic type, PC)


# :: LFC PCA ------------------------------------------------------------------

# first, mapping for samples
lfc.pca.mapping = data.frame(standard.name = rownames(oars.fecalcal.omics.lfc.saved)) %>%
  arrange(standard.name) %>%
  mutate(HM = substr(standard.name, 1, 6)) %>%
  group_by(HM) %>%
  mutate(lead.sample = lead(standard.name)) %>%
  subset(!is.na(lead.sample))
lfc.pca.mapping = merge(lfc.pca.mapping,
                        metadata.oars.stool.asv.double[,c("standard.name", "rs.selected", "timing",
                                                          "richness", "shannon", "fd")], by="standard.name") %>% distinct()
lfc.pca.mapping


oars.pc.extracter = function(data.type = x){
  if(data.type == "ASV"){
  data.subset = oars.fecalcal.omics.lfc.saved[,colnames(oars.fecalcal.omics.lfc.saved) %in%
                          colnames(oars.asv.data.glom)]
  }
  if(data.type == "Species"){
    data.subset = oars.fecalcal.omics.lfc.saved[,colnames(oars.fecalcal.omics.lfc.saved) %in%
                                                  colnames(oars.mgx.taxa)]
  }
  if(data.type == "Pathway"){
    data.subset = oars.fecalcal.omics.lfc.saved[,colnames(oars.fecalcal.omics.lfc.saved) %in%
                                                  colnames(oars.mpx.kegg.mat)]
  }
  if(data.type == "COG"){
    data.subset = oars.fecalcal.omics.lfc.saved[,colnames(oars.fecalcal.omics.lfc.saved) %in%
                                                  colnames(oars.mpx.cog.mat)]
  }
  if(data.type == "CAZy"){
    data.subset = oars.fecalcal.omics.lfc.saved[,colnames(oars.fecalcal.omics.lfc.saved) %in%
                                                  colnames(oars.mpx.cazy.mat)]
  }
  if(data.type == "Metabolite"){
    data.subset = oars.fecalcal.omics.lfc.saved[,colnames(oars.fecalcal.omics.lfc.saved) %in%
                                                  colnames(oars.mbx.raw.mat)]
  }
  # subset to intervention samples
  data.subset = data.subset[rownames(data.subset) %in% subset(lfc.pca.mapping, timing %in% c("0M", "3M"))$lead.sample,]
  # replace NA with 0
  data.subset[is.na(data.subset)] = 0
  # run PCA
  pca.results = prcomp(data.subset, center=T, scale=T)
  # extract components
  pca.components = pca.results$x[,c(1,2)] %>% data.frame() %>%
    rownames_to_column("lead.sample")
  # extract variance explained
  pca.results.var <- (pca.results$sdev^2 / sum(pca.results$sdev^2)) * 100
  pca.results.var = pca.results.var[1:2]
  # add to make one dataframe
  pca.components$var1 = pca.results.var[1]
  pca.components$var2 = pca.results.var[2]
  # add datatype
  pca.components$data.type = data.type
  return(pca.components)
}

oars.lfc.pc.data = do.call(rbind, lapply(c("ASV", "Species", "Pathway", "COG", "CAZy", "Metabolite"), function(x){
  oars.pc.extracter(data.type = x)
}))

# replace standard.name with baseline name
oars.lfc.pc.data = merge(oars.lfc.pc.data, 
                         lfc.pca.mapping, by="lead.sample")

# check correlations with diversity
oars.lfc.diversity = do.call(rbind, lapply(c("ASV", "Species", "Pathway", "COG", "CAZy", "Metabolite"), function(x){
  data.subset = subset(oars.lfc.pc.data, data.type == x)
  rich.data = cor.test(data.subset$PC2,
           data.subset$richness, method="spearman")
  shan.data = cor.test(data.subset$PC2,
                       data.subset$shannon, method="spearman")
  fd.data = cor.test(data.subset$PC2,
                       data.subset$fd, method="spearman")
  data.frame(data.type = x,
             cor.type = c("richness", "shannon", "fd"),
             pval = c(rich.data$p.value, shan.data$p.value, fd.data$p.value),
             cor = c(rich.data$estimate, shan.data$estimate, fd.data$estimate))
}))

# PC1: ASV diversity predicts ASV variance, but not others
# PC2: ASV functional redundancy predicts Metabolite variance

# :: PC Procrustes ----------------------------------------------------------------

oars.pc.grid = expand.grid(c("ASV", "Species", "Pathway", "COG", "CAZy", "Metabolite"),
            c("ASV", "Species", "Pathway", "COG", "CAZy", "Metabolite"))

oars.baseline.pca.data = rbind(oars.asv.pcoa.df%>%mutate(PC1 = Axis.1, PC2 = Axis.2, data.type = "ASV") %>% dplyr::select(standard.name, PC1, PC2, data.type),
                               oars.mgx.pcoa.df%>%mutate(PC1 = Axis.1, PC2 = Axis.2, data.type = "Species") %>% dplyr::select(standard.name, PC1, PC2, data.type),
                               oars.mpx.kegg.oars.pca.df %>%mutate(data.type = "Pathway")%>% dplyr::select(standard.name, PC1, PC2, data.type),
                               oars.mpx.cog.oars.pca.df %>%mutate(data.type = "COG")%>% dplyr::select(standard.name, PC1, PC2, data.type),
                               oars.mpx.cazy.pca.df%>%mutate(data.type = "CAZy") %>% dplyr::select(standard.name, PC1, PC2, data.type),
                               oars.mbx.pca.df%>%mutate(data.type = "Metabolite") %>% dplyr::select(standard.name, PC1, PC2, data.type))
oars.baseline.pca.data = subset(oars.baseline.pca.data, standard.name %in% oars.lfc.pc.data$standard.name)
  
# loop through datasets and perform procrustes
oars.baseline.lfc.protest = do.call(rbind, lapply(1:nrow(oars.pc.grid), function(x){
  print(x)
  grid.select = oars.pc.grid[x,]
  # subset baseline
  baseline.subset = subset(oars.baseline.pca.data, data.type == grid.select$Var1)
  rownames(baseline.subset) = baseline.subset$standard.name
  # subset LFC
  lfc.subset = subset(oars.lfc.pc.data, data.type == grid.select$Var2)
  rownames(lfc.subset) = lfc.subset$standard.name
  # procrustes
  protest.output = vegan::protest(dist(baseline.subset[,c("PC1", "PC2")]),
                 dist(lfc.subset[,c("PC1", "PC2")]))
  data.frame(Var1 = grid.select$Var1,
             Var2 = grid.select$Var2,
             cor = protest.output$scale,
             pval = protest.output$signif)
}))

ggplot(oars.baseline.lfc.protest,
       aes(x=Var1, y=Var2))+
  geom_tile(aes(fill=(cor)), color="white")+
  geom_text(aes(label = ifelse(pval < 0.05, "*", "")), color="white", vjust=0.7)+
  scale_fill_gradient2(low="blue", mid="white", high="red")+
  theme_classic()+theme(axis.text.x = element_text(angle=45, hjust=1))+
  labs(x="", y="", fill="Procrustes R")

# Baseline proteins and metabolites can predict change in proteins and metabolits, respectively

# :: ----------------------------------------------------------------------


# >>> OLD ------------------------------------------------------------------


# :: GAM Omics // defunct------------------------------------------------------------

# run GAM on omics
oars.asv.gam = omic.gam(oars.asv.data.glom[,c(1:10)], metadata.oars.stool.asv, 
                        transform="log2") 
oars.asv.gam.df = oars.asv.gam %>%
  #subset(variable == "days") %>%
  mutate(padj = p.adjust(pval, method="BH")) %>%
  arrange(padj)

# note: mgx must be forced to be positive after log for tweedie to work
oars.mgx.gam = omic.gam(oars.mgx.taxa*10000, metadata.oars.stool.asv, transform="log2")
oars.mgx.gam.df = oars.mgx.gam %>%
  subset(variable == "days") %>%
  mutate(padj = p.adjust(pval, method="BH")) %>%
  arrange(padj)

oars.kegg.gam = omic.gam(oars.mpx.kegg.mat, metadata.oars.stool.asv, transform="log2")
oars.kegg.gam.df = oars.kegg.gam %>%
  subset(variable == "days") %>%
  mutate(padj = p.adjust(pval, method="BH")) %>%
  arrange(padj)

oars.cog.gam = omic.gam(oars.mpx.cog.mat, metadata.oars.stool.asv, transform="log2")
oars.cog.gam.df  = oars.cog.gam %>%
  subset(variable == "days") %>%
  mutate(padj = p.adjust(pval, method="BH")) %>%
  arrange(padj)

oars.cazy.gam = omic.gam(oars.mpx.cazy.mat, metadata.oars.stool.asv, transform="log2")
oars.cazy.gam.df = oars.cazy.gam %>%
  subset(variable == "days") %>%
  mutate(padj = p.adjust(pval, method="BH")) %>%
  arrange(padj)

oars.mbx.gam = omic.gam(oars.mbx.raw.mat.filt.1, metadata.oars.stool.asv, transform="sqrt")
oars.mbx.gam.df = oars.mbx.gam %>%
  subset(variable == "days") %>%
  mutate(padj = p.adjust(pval, method="BH")) %>%
  arrange(padj)

# review
oars.asv.gam.df$data.type = "ASV"
oars.mgx.gam.df$data.type = "Species"
oars.kegg.gam.df$data.type = "Pathway"
oars.cog.gam.df$data.type = "COG"
oars.cazy.gam.df$data.type = "CAZy"
oars.mbx.gam.df$data.type = "Metabolites"

# merge
oars.omics.gam = rbind(oars.asv.gam.df,
                       oars.mgx.gam.df,
                       oars.kegg.gam.df,
                       oars.cog.gam.df,
                       oars.cazy.gam.df,
                       oars.mbx.gam.df) %>% as.data.frame()
# subset to sig
oars.omics.gam = oars.omics.gam %>%
  #mutate(padj = p.adjust(pval, method="BH")) %>%
  subset(padj < 0.20)

subset(oars.omics.gam, variable == "days")$data.type %>% table()
subset(oars.omics.gam, variable == "HM")$data.type %>% table()

# find strongest associations and plot them
subset(oars.omics.gam, variable == "days" & data.type == "ASV") %>% head()
subset(oars.omics.gam, variable == "days" & data.type == "Species") %>% head()
subset(oars.omics.gam, variable == "days" & data.type == "Pathway") %>% head()
subset(oars.omics.gam, variable == "days" & data.type == "COG") %>% head()
subset(oars.omics.gam, variable == "days" & data.type == "CAZy") %>% head()
subset(oars.omics.gam, variable == "days" & data.type == "Metabolite") %>% head()

# select to plot
omic.gam.plot(subset(oars.omics.gam, variable == "days" & data.type == "ASV")$feature[1])
omic.gam.plot(subset(oars.omics.gam, variable == "days" & data.type == "ASV")$feature[2])
omic.gam.plot(subset(oars.omics.gam, variable == "days" & data.type == "ASV")$feature[3])
omic.gam.plot(subset(oars.omics.gam, variable == "days" & data.type == "ASV")$feature[4])
omic.gam.plot(subset(oars.omics.gam, variable == "days" & data.type == "ASV")$feature[5])

omic.gam.plot(subset(oars.omics.gam, variable == "days" & data.type == "Species")$feature[1])
omic.gam.plot(subset(oars.omics.gam, variable == "days" & data.type == "Species")$feature[2])
omic.gam.plot(subset(oars.omics.gam, variable == "days" & data.type == "Species")$feature[3])
omic.gam.plot(subset(oars.omics.gam, variable == "days" & data.type == "Species")$feature[4])
omic.gam.plot(subset(oars.omics.gam, variable == "days" & data.type == "Species")$feature[5]) # *

omic.gam.plot(subset(oars.omics.gam, variable == "days" & data.type == "Pathway")$feature[1]) # **
omic.gam.plot(subset(oars.omics.gam, variable == "days" & data.type == "Pathway")$feature[2]) # **
omic.gam.plot(subset(oars.omics.gam, variable == "days" & data.type == "Pathway")$feature[3]) # **
omic.gam.plot(subset(oars.omics.gam, variable == "days" & data.type == "Pathway")$feature[4]) # **
omic.gam.plot(subset(oars.omics.gam, variable == "days" & data.type == "Pathway")$feature[5])

omic.gam.plot(subset(oars.omics.gam, variable == "days" & data.type == "COG")$feature[1])
omic.gam.plot(subset(oars.omics.gam, variable == "days" & data.type == "COG")$feature[2])
omic.gam.plot(subset(oars.omics.gam, variable == "days" & data.type == "COG")$feature[3])
omic.gam.plot(subset(oars.omics.gam, variable == "days" & data.type == "COG")$feature[4])
omic.gam.plot(subset(oars.omics.gam, variable == "days" & data.type == "COG")$feature[5]) # **

omic.gam.plot(subset(oars.omics.gam, variable == "days" & data.type == "CAZy")$feature[1])# **
omic.gam.plot(subset(oars.omics.gam, variable == "days" & data.type == "CAZy")$feature[2])# **
omic.gam.plot(subset(oars.omics.gam, variable == "days" & data.type == "CAZy")$feature[3])# *
omic.gam.plot(subset(oars.omics.gam, variable == "days" & data.type == "CAZy")$feature[4])# **
omic.gam.plot(subset(oars.omics.gam, variable == "days" & data.type == "CAZy")$feature[5])# **

omic.gam.plot(subset(oars.omics.gam, variable == "days" & data.type == "Metabolite")$feature[1])
omic.gam.plot(subset(oars.omics.gam, variable == "days" & data.type == "Metabolite")$feature[2])
omic.gam.plot(subset(oars.omics.gam, variable == "days" & data.type == "Metabolite")$feature[3])# **
omic.gam.plot(subset(oars.omics.gam, variable == "days" & data.type == "Metabolite")$feature[4])# **
omic.gam.plot(subset(oars.omics.gam, variable == "days" & data.type == "Metabolite")$feature[5])

# these results are messy



# :: MPX CAZyme taxa contributions ----------------------------------------

# Which taxa contribute to these CAZymes

oars.mpx.cazy.tax.mat = readRDS("./metaproteomics/2025_06_28_oars_mpx_cazy_tax.Rds")
# subset to sig CAZymes

sig.cazy = c(
  subset(oars.mpx.maaslin.036, padj < 0.2)$feature,
  subset(oars.mpx.maaslin.6912, padj < 0.2)$feature) %>% unique()

oars.mpx.cazy.tax.mat = oars.mpx.cazy.tax.mat[,grepl(paste(c(sig.cazy), collapse="|"), colnames(oars.mpx.cazy.tax.mat))]
# visualize contributions

oars.cazy.tax.contributions = reshape2::melt(oars.mpx.cazy.tax.mat) %>%
  tidyr::separate(col=Var2, into=c("CAZy", "Taxonomy"), sep="_", remove=F) %>%
  group_by(Var2) %>%
  mutate(sum.intensity = (sum(value))) %>%
  dplyr::select(CAZy, Taxonomy, sum.intensity) %>% distinct() %>% data.frame()
# sankey plot

library(ggalluvial)
oars.cazy.tax.contributions.plot = ggplot(oars.cazy.tax.contributions %>%
                                            # collapse low abundant to "Other"
                                            group_by(Taxonomy) %>%
                                            mutate(Taxonomy = ifelse(sum(sum.intensity) < 2e+9, "Other.", Taxonomy)) %>%
                                            group_by(CAZy) %>%
                                            mutate(CAZy = ifelse(sum(sum.intensity) < 2e+9, "Other", CAZy)),
                                          aes(y = sum.intensity, 
                                              axis1 = reorder(Taxonomy, -sum.intensity), 
                                              axis2 = reorder(CAZy, -sum.intensity))) +
  ggalluvial::geom_alluvium(aes(fill = CAZy), width = 1/12) +
  scale_fill_manual(values=RColorBrewer::brewer.pal(n=6, name= "Set2"))+
  ggalluvial::geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = gsub("\\.", "", after_stat(stratum))),
             size=3) +
  theme_void()+theme(legend.position="none")
oars.cazy.tax.contributions.plot


# :: MPX Butyrate Maaslin2---------------------------------------------------------


oars.mpx.but.mat = readRDS("./metaproteomics/2025_06_28_oars_mpx_but.Rds")

oars.mpx.but.mat = oars.mpx.but.mat[,grepl(paste(c("thl", "hbd", "crt", "bcd", "but", "buk"), collapse="|"), colnames(oars.mpx.but.mat))]

rownames(oars.mpx.but.mat) = gsub("_", "-", rownames(oars.mpx.but.mat))

subset(metadata.oars.stool.mpx, timing %in% c("0M", "3M", "6M") & compliant==TRUE) %>% rownames() %in% rownames(oars.mpx.oars.only.protein.mat)

oars.mpx.but.maaslin.036 = Maaslin2::Maaslin2(input_data = (oars.mpx.but.mat),
                                              input_metadata = subset(metadata.oars.stool.mpx, timing %in% c("0M", "3M", "6M") & compliant==TRUE),
                                              output = "~/Downloads",
                                              fixed_effects = c("oars.days", "diagnosis", "adj.fiber"),  # Example fixed effects
                                              random_effects = c("HM"),       # Example random effects
                                              normalization = "NONE",                       # Total Sum Scaling normalization
                                              transform = "LOG",                           # Log transformation
                                              analysis_method = "LM",                      # Linear model
                                              plot_scatter = FALSE,                        # Disable scatterplot generation
                                              plot_heatmap = FALSE,                         # Keep heatmap generation (optional)
                                              max_significance = 0.05,                     # Significance threshold for q-values
                                              standardize = TRUE                           # Disable standardization (optional)
)

oars.mpx.but.maaslin.036 = oars.mpx.but.maaslin.036$results %>% data.frame() %>% arrange(pval)


# recalculate padj minus diagnosis fixed effect
oars.mpx.but.maaslin.036 = subset(oars.mpx.but.maaslin.036, metadata == "oars.days") %>%
  mutate(padj = p.adjust(pval, method="BH"))

oars.mpx.but.maaslin.036.volcano = ggplot(subset(oars.mpx.but.maaslin.036, value == "oars.days"),
                                          aes(x=coef, y=-log(pval)))+
  geom_point(shape=21, aes(fill=coef))+
  geom_hline(yintercept=-log10(0.05), linetype=2, alpha=0.5)+
  scale_fill_gradient2(low="blue", high="red")+
  ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, feature, NA)),
                           size=2.5)+
  theme_minimal()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12))+
  
  labs(x="Adjusted Coefficient",
       title="Treatment: MPX (Butyrate)")
oars.mpx.but.maaslin.036.volcano

oars.mpx.but.maaslin.6912 = Maaslin2::Maaslin2(input_data = oars.mpx.but.mat,
                                               input_metadata = subset(metadata.oars.stool.mpx, timing %in% c("6M", "9M", "12M")& compliant==TRUE),
                                               output = "~/Downloads",
                                               fixed_effects = c("oars.days", "diagnosis", "adj.fiber"),  # Example fixed effects
                                               random_effects = c("HM"),       # Example random effects
                                               normalization = "NONE",                       # Total Sum Scaling normalization
                                               transform = "LOG",                           # Log transformation
                                               analysis_method = "LM",                      # Linear model
                                               plot_scatter = FALSE,                        # Disable scatterplot generation
                                               plot_heatmap = FALSE,                         # Keep heatmap generation (optional)
                                               max_significance = 0.05,                     # Significance threshold for q-values
                                               standardize = TRUE                           # Disable standardization (optional)
)

oars.mpx.but.maaslin.6912 = oars.mpx.but.maaslin.6912$results %>% data.frame() %>% arrange(pval)

# recalculate padj minus diagnosis fixed effect
oars.mpx.but.maaslin.6912 = subset(oars.mpx.but.maaslin.6912, metadata == "oars.days") %>%
  mutate(padj = p.adjust(pval, method="BH"))


oars.mpx.but.maaslin.6912.volcano = ggplot(oars.mpx.but.maaslin.6912,
                                           aes(x=coef, y=-log(pval)))+
  geom_point(shape=21, aes(fill=coef))+
  geom_hline(yintercept=-log10(0.05), linetype=2, alpha=0.5)+
  scale_fill_gradient2(low="blue", high="red")+
  ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, gsub("\\..*", "", feature), NA)),
                           size=2.5)+
  theme_minimal()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12))+
  
  labs(x="Adjusted Coefficient",
       title="Washout: MPX (Butyrate)")
oars.mpx.but.maaslin.6912.volcano
# 


# :: MPX Butyrate taxa contributions -----------------------------------------------------


oars.mpx.but.tax.mat = readRDS("./metaproteomics/2025_06_28_oars_mpx_but.Rds")

oars.mpx.but.tax.mat = oars.mpx.but.tax.mat[,grepl(paste(c("thl", "hbd", "crt", "bcd", "but", "buk"), collapse="|"), colnames(oars.mpx.but.tax.mat))]

oars.mpx.but.contributions = reshape2::melt(oars.mpx.but.tax.mat) %>%
  tidyr::separate(col=Var2, into=c("gene", "Taxonomy"), sep="_", remove=F) %>%
  group_by(Var2) %>%
  mutate(sum.intensity = (sum(value))) %>%
  dplyr::select(gene, Taxonomy, sum.intensity) %>% distinct() %>% data.frame()
# sankey plot

library(ggalluvial)
oars.mpx.but.contributions.plot = ggplot(oars.mpx.but.contributions %>%
                                           # collapse low abundant to "Other"
                                           group_by(Taxonomy) %>%
                                           mutate(Taxonomy = ifelse(sum(sum.intensity) < 2e+9, "Other.", Taxonomy)) %>%
                                           group_by(gene) %>%
                                           mutate(gene = ifelse(sum(sum.intensity) < 2e+9, "Other", gene)),
                                         aes(y = sum.intensity, 
                                             axis1 = reorder(Taxonomy, -sum.intensity), 
                                             axis2 = reorder(gene, -sum.intensity))) +
  ggalluvial::geom_alluvium(aes(fill = gene), width = 1/12) +
  scale_fill_manual(values=RColorBrewer::brewer.pal(n=6, name= "Set2"))+
  ggalluvial::geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = gsub("\\.", "", after_stat(stratum))),
             size=3) +
  theme_void()+theme(legend.position="none")
oars.mpx.but.contributions.plot

# it's almost all Clostridia



# :: Networks Prep -------------------------------------------------------------

# >> glasso: on ASV, MGX, MPX etc >>(combined?) at each timepoint 
# or 
# granger: use previous timepoint to predict future timepoint


#colnames(oars.multi) = make.names(colnames(oars.multi))

# MGX subset to nominally sig taxa
oars.mgx.glasso.data = oars.mgx.taxa[,colnames(oars.mgx.taxa) %in% 
                                       c(subset(oars.mgx.maaslin.036, pval < 0.05)$feature,
                                         subset(oars.mgx.maaslin.6912, pval < 0.05)$feature)]
dim(oars.mgx.glasso.data) # 147 features
# trace out 0 counts
oars.mgx.glasso.data.pa = oars.mgx.glasso.data
oars.mgx.glasso.data.pa[oars.mgx.glasso.data.pa !=0] = 1
# clr
oars.mgx.glasso.data = compositions::clr(oars.mgx.glasso.data)
oars.mgx.clr.pseudo = min(oars.mgx.glasso.data)-1
oars.mgx.glasso.data[oars.mgx.glasso.data==0] = oars.mgx.clr.pseudo

# MPX subset to nominally sig taxa
oars.mpx.glasso.data = oars.mpx.oars.only.protein.mat[,colnames(oars.mpx.oars.only.protein.mat) %in% 
                                                        c(slice_min(subset(oars.mpx.maaslin.036, pval < 0.05), pval, n=50)$feature,
                                                          slice_min(subset(oars.mpx.maaslin.6912, pval < 0.05), pval, n=50)$feature)]
oars.mpx.glasso.data[is.na(oars.mpx.glasso.data)] = 0
oars.mpx.glasso.data = oars.mpx.glasso.data+(min(oars.mpx.glasso.data[oars.mpx.glasso.data!=0])/2)
oars.mpx.glasso.data = log2(oars.mpx.glasso.data)

# merge into single 
oars.multi = merge(oars.mgx.glasso.data,
                   oars.mpx.glasso.data, by="row.names")
rownames(oars.multi) = oars.multi[,1]
oars.multi[,1] = NULL

# merged and processed data
oars.multi

# make map
oars.multi.map = rbind(data.frame(name = colnames(oars.mpx.glasso.data),
                                  type = "protein"),
                       data.frame(name = colnames(oars.mgx.glasso.data),
                                  type = "taxa")) %>% data.frame()


# ::  Networks --------------------------------------------------------

# first, check glasso CV curve
omic.cor(oars.multi, type="glasso", optimize = T)
# does not hit minimum

# so, use spearman

# TREATMENT
omic.cor.plot(
  # run correlation calculator within function
  omic.cor(oars.multi[rownames(oars.multi) %in% subset(metadata.oars.stool, 
                                                       timing %in% c("0M", "3M", "6M"))$standard.name,], 
           type="spearman", optimize = F),
  # resume plot function
  mapping = oars.multi.map,
  threshold = 0.5,
  label=T)

# WASHOUT
omic.cor.plot(
  # run correlation calculator within function
  omic.cor(oars.multi[rownames(oars.multi) %in% subset(metadata.oars.stool, 
                                                       timing %in% c("6M", "9M", "12M"))$standard.name,], 
           type="spearman", optimize = F),
  # resume plot function
  mapping = oars.multi.map,
  threshold = 0.5,
  label=T)

# COMBINED
omic.cor.plot(
  # run correlation calculator within function
  omic.cor(oars.multi[rownames(oars.multi) %in% metadata.oars.stool$standard.name,], 
           type="spearman", optimize = F),
  # resume plot function
  mapping = oars.multi.map,
  threshold = 0.4,
  label=T)


# :: Baseline cluster? ----------------------------------------------------

# 2025_07_08  James noticed that there are 2 clusters of Stool Water and Richness at baseline

ggplot(subset(metadata.oars.stool.asv, timing == "0M" & compliant == TRUE),
       aes(x=fd, y=stool_water_perc))+
  geom_smooth(method="lm")+
  ggpubr::stat_cor(method="spearman")+
  geom_point(shape=21)+
  theme_minimal()

pheatmap::pheatmap(as.data.frame(Hmisc::rcorr(as.matrix(metadata.oars.stool.asv[,c("richness", "shannon" ,"fd", "stool_water_perc", "but.i", "but.ii")]))$r))

# nothing really notable


# :: ASV Phylotree // defunct --------------------------------------------------------

# construct phylogenetic tree, and show LFC and ICC per HM; 1 for RS and 1 for washout
# on most prevalent taxa (rather than maaslin2)

# taxa to keep
oars.asv.selected.taxa

# calculate ICC and LFC of taxa
oars.asv.prev.icc = do.call(rbind, lapply(oars.asv.selected.taxa, function(taxa){
  # first, divide by Treatment and Washout periods
  #taxa = "f__Peptostreptococcaceae_256921"
  data.taxa = oars.asv.data.glom[,taxa] %>% data.frame()
  colnames(data.taxa)[1] = "taxa"
  data.taxa$standard.name = rownames(data.taxa)
  do.call(rbind, lapply(c("treatment", "washout"), function(phase){
    print(paste(taxa, phase))
    if(phase == "treatment"){
      data.meta = subset(metadata.oars.stool.asv, timing %in% c("0M", "3M", "6M") & compliant==TRUE)
    }
    if(phase == "washout"){
      data.meta = subset(metadata.oars.stool.asv, timing %in% c("6M", "9M", "12M") & compliant==TRUE)
    }
    # merge
    data.taxa = merge(data.taxa,
                      data.meta, by="standard.name")
    # process data
    data.taxa$taxa = log2((data.taxa$taxa+0.5)/50000)
    if(sd(data.taxa$taxa) == 0){
      lm.results.df = data.frame(
        feature = taxa,
        coef = NA,
        pval = NA,
        icc = NA,
        phase = phase)
    }else{
      # build model
      lm.results = lmerTest::lmer(scale(taxa) ~ scale(oars.days) + diagnosis + (adj.fiber) + (1|HM), data.taxa)
      lm.results.df = lm.results %>% summary() %>% coef() %>% data.frame()
      # extract ICC
      var_comp <- as.data.frame(lme4::VarCorr(lm.results))
      var_group <- var_comp$vcov[var_comp$grp == "HM"]
      var_resid <- var_comp$vcov[var_comp$grp == "Residual"]
      icc_adj <- var_group / (var_group + var_resid)
      
      #
      lm.results.df = data.frame(
        feature = taxa,
        coef = lm.results.df$Estimate[2],
        pval = lm.results.df$`Pr...t..`[2],
        icc = icc_adj,
        phase = phase)
    }
    
  }))
}))

oars.asv.prev.icc

# load tree
oars.tree = readRDS("./2025_06_29_oars_tree.Rds")


# clip down to sig features
pruned_tree <- ape::drop.tip(oars.tree, setdiff(oars.tree$tip.label, 
                                                oars.asv.selected.taxa))
# here is the final tree
pruned_tree = phyloseq::phy_tree(pruned_tree)
phyloseq::taxa_names(pruned_tree)
tree_meta = readRDS("./2025_06_29_oars_tree_meta.Rds")

tree_meta$feature = (tree_meta$LCA)

# now plot

# add group-level ICC
tree_meta = merge(tree_meta,
                  # ICC data
                  oars.asv.prev.icc,
                  by="feature")

# binarize icc_rs
tree_meta$icc_rs_bin = ifelse(tree_meta$icc >= 0.5, "high",
                              ifelse(tree_meta$icc < 0.5 & tree_meta$icc > 0, "low", "NA"))

# Add phylum data
tree_meta$phylum = as.factor(gsub("_.*", "", gsub("p__", "", tree_meta$Phylum)))
unique(tree_meta$phylum)
unique(tree_meta$feature)

# fix feature names
pruned_tree$tip.label = gsub(paste(c("p__", "o__", "c__", "f__", "g__", "s__"), collapse="|"), "", pruned_tree$tip.label)
pruned_tree$tip.label= gsub("_", " ", pruned_tree$tip.label)
pruned_tree$tip.label %>% unique()

tree_meta$feature = gsub(paste(c("p__", "o__", "c__", "f__", "g__", "s__"), collapse="|"), "", tree_meta$feature)
tree_meta$feature = gsub("_", " ", tree_meta$feature)
tree_meta$feature %>% unique()

# make extra object for tree object
tree_meta_phylum = tree_meta
tree_meta_phylum$feature_label = tree_meta_phylum$feature
pruned_tree$feature = pruned_tree$tip.label

# plot LFC heatmap with ICC bar

# PLOT
library("ggtree")
library("ggtreeExtra")
# Base circular tree with top left quarter removed

p0 <- ggtree(pruned_tree, layout = "fan", open.angle = 45, size=0.5) %<+% tree_meta_phylum[,c("feature","feature_label", "phylum")] +
  #coord_polar(theta = "y", start = pi/4, direction = -1)+
  #aes(color = phylum) +  # Move coloring to base layercoord_polar(theta = "y", start = pi/2, direction = -1) +
  #scale_color_manual(values = c("Bacillota" = "limegreen", "Bacteroidota" = "salmon" ),
  #                   na.value = "black") +  # Define colors
  guides(color = guide_legend(title = "Phylum"))+
  geom_tiplab(
    aes(label = paste(" ", feature_label, " ", sep=""), color=phylum),  # Use feature as label
    offset = 0.1,         # Adjust to place labels beyond bars (0.05 + 0.15 + extra)
    size = 3,              # Match size of previous attempts
    color = "black",       # Consistent color
    align = TRUE,          # Align labels at equal radius
    linetype = 3,          # No connecting lines
    linewidth=5
  )
p0
p0 <- rotate_tree(p0, angle=85)

# Heatmap tiles
p1 <- p0 + geom_fruit(
  data = tree_meta,
  geom = geom_tile,
  mapping = aes(y = feature, x = phase, fill = (icc)),
  color = "white",
  offset = -0.001,
  pwidth = -0.25,
  axis.params = list(axis = "none")
) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.5,
                       breaks=c(0,0.25, 0.5, 0.75, 1),
                       limits=c(0,1)) +
  #scale_color_manual(values=c("blue", "red", "white"))+
  labs(fill = "Intra-Class Correlation")
p1

margins = c(0,0,0,0)

p1 + theme(plot.margin = margin(margins, "cm"))+ expand_limits(x = c(-.01, 1))


# // defunct --------------------------------------------------------------

# :: MPX PROTEIN // defunct ----------------------------------------------------------

# load protein data
oars.mpx.oars.only.protein.mat = readRDS("./metaproteomics/2025_06_28_oars_mpx_protein.Rds")

metadata.oars.stool.mpx = subset(metadata.oars.stool,
                                 standard.name %in% rownames(oars.mpx.oars.only.protein.mat))
# create map for protein names
oars.mpx.protein.names = data.frame(good = colnames(oars.mpx.oars.only.protein.mat),
                                    feature = make.names(colnames(oars.mpx.oars.only.protein.mat)))
length(unique(oars.mpx.protein.names$good))
# 1476 proteins

# :: MPX PROTEIN PCA -------------------------------------------------------------

oars.mpx.protein.oars.pca = oars.mpx.oars.only.protein.mat[rownames(oars.mpx.oars.only.protein.mat) %in% metadata.oars.stool$standard.name,]
oars.mpx.protein.oars.pca[is.na(oars.mpx.protein.oars.pca)] = 0
# log transform
oars.mpx.protein.oars.pca = log2(oars.mpx.protein.oars.pca+(min(oars.mpx.protein.oars.pca[oars.mpx.protein.oars.pca!=0])/2))
oars.mpx.protein.oars.pca = prcomp((oars.mpx.protein.oars.pca), scale=T)
oars.mpx.protein.oars.pca.df = oars.mpx.protein.oars.pca$x[,c(1,2)] %>% data.frame() %>%
  rownames_to_column("standard.name")
oars.mpx.protein.oars.pca.df = merge(oars.mpx.protein.oars.pca.df,
                                     metadata.oars.stool.mpx, by="standard.name")
oars.mpx.protein.oars.pca.var <- (oars.mpx.protein.oars.pca$sdev^2 / sum(oars.mpx.protein.oars.pca$sdev^2)) * 100

rownames(oars.mpx.protein.oars.pca.df) = oars.mpx.protein.oars.pca.df$standard.name
# permanova
set.seed(25)
t1 = Sys.time()
oars.mpx.protein.oars.pca.permanova = vegan::adonis2(dist(oars.mpx.protein.oars.pca$x) ~ on.rs + adj.fiber,
                                                     oars.mpx.protein.oars.pca.df %>% mutate(on.rs = as.factor(ifelse(oars.on.rs=="onRS", "onRS", "offRS"))),
                                                     strata = oars.mpx.protein.oars.pca.df$HM,
                                                     by="margin")
t2 = Sys.time()
t2 - t1

oars.mpx.protein.pca.plot <- ggplot(
  data=oars.mpx.protein.oars.pca.df %>% group_by(HM) %>% arrange(stool_date_rec_v2), 
  aes(x=PC1, y=PC2))+
  geom_path(aes(group=HM), color="black", linetype=2, alpha=0.5, linewidth=0.3) + 
  geom_point(aes(fill=rs.col, shape=baseline), size=2)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values=c("grey", labelcolors$cols[c(1,1,4,5,5,8,8,9)]))+
  annotate(geom="text",
           x=-15,
           y=-20,
           label = paste(paste("Treatment\nR²: ", round(data.frame(oars.mpx.protein.oars.pca.permanova)[1,3], 3)*100, "%",
                               "\n p: ", round(data.frame(oars.mpx.protein.oars.pca.permanova)[1,5], 3), sep="")),
           size=3.3)+
  facet_wrap(~"Protein PCA")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text=element_text(size=10))+
  labs(x=paste("PC1: ", round((oars.mpx.protein.oars.pca.var)[1], digits=2), "%", sep=""), 
       y=paste("PC2: ", round((oars.mpx.protein.oars.pca.var)[2], digits=2), "%", sep=""))
oars.mpx.protein.pca.plot

# very sig


# :: MPX PROTEIN Maaslin2 ----------------------------------------------------


subset(metadata.oars.stool.mpx, timing %in% c("0M", "3M", "6M") & compliant==TRUE) %>% rownames() %in% rownames(oars.mpx.oars.only.protein.mat)

oars.mpx.protein.maaslin.036 = Maaslin2::Maaslin2(input_data = (oars.mpx.oars.only.protein.mat),
                                                  input_metadata = subset(metadata.oars.stool.mpx, timing %in% c("0M", "3M", "6M") & compliant==TRUE),
                                                  output = "~/Downloads",
                                                  fixed_effects = c("oars.days"),  # Example fixed effects
                                                  random_effects = c("HM"),       # Example random effects
                                                  normalization = "NONE",                       # Total Sum Scaling normalization
                                                  transform = "LOG",                           # Log transformation
                                                  analysis_method = "LM",                      # Linear model
                                                  plot_scatter = FALSE,                        # Disable scatterplot generation
                                                  plot_heatmap = FALSE,                         # Keep heatmap generation (optional)
                                                  max_significance = 0.05,                     # Significance threshold for q-values
                                                  standardize = TRUE                           # Disable standardization (optional)
)

oars.mpx.protein.maaslin.036 = oars.mpx.protein.maaslin.036$results %>% data.frame() %>% arrange(pval)


# recalculate padj minus diagnosis fixed effect
oars.mpx.protein.maaslin.036 = subset(oars.mpx.protein.maaslin.036, metadata == "oars.days") %>%
  mutate(padj = p.adjust(pval, method="BH"))

# Enriched: Nitrogen reg protein, Isocitrate-isopropylmalate dehydrogenase, tRNA methylthiotransferase, Histidinol dehydrogenase, Pentose-5-phosphate-3-epimerase

# fix feature names
oars.mpx.protein.maaslin.036$feature = oars.mpx.protein.names$good[match(oars.mpx.protein.maaslin.036$feature, oars.mpx.protein.names$feature)]

oars.mpx.protein.maaslin.036.volcano = ggplot(oars.mpx.protein.maaslin.036,
                                              aes(x=coef, y=(padj)))+
  geom_point(shape=21, aes(fill=coef))+
  geom_hline(yintercept=(0.2), linetype=2, alpha=0.5)+
  geom_vline(xintercept=c(-1,1), linetype=2, alpha=0.5)+
  scale_fill_gradient2(low="blue", high="red")+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.05,0.1, 0.20, 0.5, 1))+
  ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, gsub("\\..*", "", feature), NA)),
                           size=2.5)+  
  facet_wrap(~"Treatment: Protein")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text=element_text(size=10))+
  labs(x="Adjusted Coefficient",
       y="FDR")
oars.mpx.protein.maaslin.036.volcano

oars.mpx.protein.maaslin.6912 = Maaslin2::Maaslin2(input_data = oars.mpx.oars.only.protein.mat,
                                                   input_metadata = subset(metadata.oars.stool.mpx, timing %in% c("6M", "9M", "12M")& compliant==TRUE),
                                                   output = "~/Downloads",
                                                   fixed_effects = c("oars.days"),  # Example fixed effects
                                                   random_effects = c("HM"),       # Example random effects
                                                   normalization = "NONE",                       # Total Sum Scaling normalization
                                                   transform = "LOG",                           # Log transformation
                                                   analysis_method = "LM",                      # Linear model
                                                   plot_scatter = FALSE,                        # Disable scatterplot generation
                                                   plot_heatmap = FALSE,                         # Keep heatmap generation (optional)
                                                   max_significance = 0.05,                     # Significance threshold for q-values
                                                   standardize = TRUE                           # Disable standardization (optional)
)

oars.mpx.protein.maaslin.6912 = oars.mpx.protein.maaslin.6912$results %>% data.frame() %>% arrange(pval)

# recalculate padj minus diagnosis fixed effect
oars.mpx.protein.maaslin.6912 = subset(oars.mpx.protein.maaslin.6912, metadata == "oars.days") %>%
  mutate(padj = p.adjust(pval, method="BH"))

# fix feature names
oars.mpx.protein.maaslin.6912$feature = oars.mpx.protein.names$good[match(oars.mpx.protein.maaslin.6912$feature, oars.mpx.protein.names$feature)]

oars.mpx.protein.maaslin.6912.volcano = ggplot(oars.mpx.protein.maaslin.6912,
                                               aes(x=coef, y=(padj)))+
  geom_point(shape=21, aes(fill=coef))+
  geom_hline(yintercept=(0.2), linetype=2, alpha=0.5)+
  geom_vline(xintercept=c(-1,1), linetype=2, alpha=0.5)+
  scale_fill_gradient2(low="blue", high="red")+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.05,0.1, 0.20, 0.5, 1))+
  ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, gsub("\\..*", "", feature), NA)),
                           size=2.5)+  
  facet_wrap(~"Washout: Protein")+
  theme_classic()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text=element_text(size=10))+
  labs(x="Adjusted Coefficient",
       y="FDR")
oars.mpx.protein.maaslin.6912.volcano
# 


# :: CAZyme per RS // defunct --------------------------------------------------------

# Do specific RS impact CAZyme profiles differently?

oars.mpx.oars.pca.rs = oars.mpx.cazy.mat[rownames(oars.mpx.cazy.mat) %in% subset(metadata.oars.stool, !is.na(RS_Name) & compliant == TRUE)$standard.name,]
oars.mpx.oars.pca.rs[is.na(oars.mpx.oars.pca.rs)] = 0
# log transform
oars.mpx.oars.pca.rs = log2(oars.mpx.oars.pca.rs+(min(oars.mpx.oars.pca.rs[oars.mpx.oars.pca.rs!=0])/2))
oars.mpx.oars.pca.rs = prcomp((oars.mpx.oars.pca.rs), scale=T)
oars.mpx.oars.pca.rs.df = oars.mpx.oars.pca.rs$x[,c(1,2)] %>% data.frame() %>%
  rownames_to_column("standard.name")
oars.mpx.oars.pca.rs.df = merge(oars.mpx.oars.pca.rs.df,
                                metadata.oars.stool.mpx, by="standard.name")
oars.mpx.oars.pca.rs.var <- (oars.mpx.oars.pca.rs$sdev)^2 / sum(oars.mpx.oars.pca.rs$sdev^2) * 100

rownames(oars.mpx.oars.pca.rs.df) = oars.mpx.oars.pca.rs.df$standard.name
# permanova
set.seed(25)
t1 = Sys.time()
oars.mpx.oars.pca.rs.permanova = vegan::adonis2(dist(oars.mpx.oars.pca.rs$x) ~ RS_Name + adj.fiber,
                                                oars.mpx.oars.pca.rs.df,
                                                strata = oars.mpx.oars.pca.rs.df$HM,
                                                by="margin")
t2 = Sys.time()
t2 - t1

oars.mpx.oars.pca.rs.plot <- ggplot(
  data=oars.mpx.oars.pca.rs.df %>% group_by(HM) %>% arrange(stool_date_rec_v2), 
  aes(x=PC1, y=PC2))+
  geom_path(aes(group=HM), color="black", linetype=2, alpha=0.5, linewidth=0.3) + 
  geom_point(aes(fill=RS_Name, shape=baseline), size=2)+
  scale_shape_manual(values=c(21))+
  scale_fill_manual(values=labelcolors$cols[c(1,1,4,5,5,8,8,9)])+
  #geom_text(aes(label=RS_Name))+
  theme_minimal()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12))+
  labs(x=paste("PC1: ", round((oars.mpx.oars.pca.rs.var)[1], digits=2), "%", sep=""), 
       y=paste("PC2: ", round((oars.mpx.oars.pca.rs.var)[2], digits=2), "%", sep=""),
       title=paste(paste("Treatment R²: ", round(data.frame(oars.mpx.oars.pca.rs.permanova)[1,3], 3)*100, "%",
                         "  p: ", round(data.frame(oars.mpx.oars.pca.rs.permanova)[1,5], 3), sep="")), sep="")
oars.mpx.oars.pca.rs.plot

oars.mpx.rs.prep = delta.omic.prepare(oars.mpx.cazy.mat, 
                                      normalize=F,
                                      min.abun = 1,
                                      prev=0.2)
oars.mpx.rs.prep = merge(oars.mpx.rs.prep,
                         oars.mpx.oars.pca.rs.df[,c("standard.name", "RS_Name")],
                         by="standard.name", all.x=T)
# loop through and assess interaction with RS
oars.mpx.rs.lm = do.call(rbind, lapply(colnames(oars.mpx.rs.prep)[3:(ncol(oars.mpx.rs.prep)-8)], function(feature){
  # feature = colnames(oars.mpx.rs.prep)[3:(ncol(oars.mpx.rs.prep)-8)][2]
  data.subset = oars.mpx.rs.prep[,c("HM", "phase", "reltiming", "response", "diagnosis","RS_Name", "adj.fiber")]
  # process
  data.subset$feature = oars.mpx.rs.prep[,colnames(oars.mpx.rs.prep) == feature]
  data.subset$feature = log2(data.subset$feature)
  # calculate lfc
  data.subset = data.subset %>%
    subset(!HM %in% c("HM0618"))%>%
    group_by(HM, phase) %>%
    mutate(lfc = feature - feature[reltiming == "pre"]) %>%
    subset(reltiming != "pre")
  # run lm
  lm.output = lmerTest::lmer(lfc ~ RS_Name + (1|HM), data.subset) %>% summary() %>% coef() %>% data.frame()
  # save output
  data.frame(
    feature = feature,
    RS_Name = gsub("RS_Name", "", rownames(lm.output)[-8]),
    coef = lm.output[-8,1],
    pval = lm.output[-8,5])[-1,]
}))

prcomp(t(reshape2::acast(oars.mpx.rs.lm,
                         feature ~ RS_Name, value.var="coef")), scale=T)$x %>%
  data.frame() %>%
  rownames_to_column("RS_Name") %>%
  ggplot(aes(x=PC1, y=PC2))+
  geom_point()+
  ggrepel::geom_text_repel(aes(label=RS_Name))+
  theme_classic()

oars.mpx.rs.lm %>% 
  mutate(padj = p.adjust(pval, method="BH")) %>%
  subset(pval < 0.05) %>%
  arrange(pval) %>%
  mutate(RS_Name = factor(RS_Name, levels=rs.names)) %>%
  group_by(feature) %>%
  #mutate(coef = scale(coef)) %>%
  # head() %>%
  reshape2::acast(feature ~ RS_Name, value.var="coef") %>%
  as.matrix() %>%
  reshape2::melt() %>%
  mutate(value = ifelse(is.na(value), 0, value))%>%
  mutate(value = ifelse(value > 5, 5,
                        ifelse(value < -5, -5, value)))%>%
  reshape2::acast(Var1 ~ Var2, value.var="value") %>%
  
  pheatmap::pheatmap(color=colorRampPalette(c("blue","white", "red"))(100))

# plot
ggplot(aes(x=coef, y=-log10(pval)))+
  geom_point(aes(fill=scale(coef)), shape=21, size=3)+
  scale_fill_gradient2(low="blue", high="red")+
  ggrepel::geom_text_repel(aes(label=feature), size=3)+
  theme_minimal()+theme(legend.position="none")+
  facet_wrap(~RS_Name)



# :: MPX PROTEIN lmer // defunct-------------------------------------------------------------

dim(oars.mpx.oars.only.protein.mat)

# for CAZy, remove "-"
rownames(oars.mpx.oars.only.protein.mat) = gsub("_", "-", rownames(oars.mpx.oars.only.protein.mat))

oars.mpx.protein.prep = delta.omic.prepare(oars.mpx.oars.only.protein.mat, 
                                           normalize=F,
                                           min.abun = 1,
                                           prev=0.2)
# create map for protein names
oars.mpx.protein.protein.names = data.frame(good = colnames(oars.mpx.oars.only.protein.mat),
                                            feature = make.names(colnames(oars.mpx.oars.only.protein.mat))) %>%
  subset(feature %in% colnames(oars.mpx.protein.prep))

## LMER
oars.mpx.protein.lmer = delta.omic.lmer(
  split = FALSE,
  data = oars.mpx.protein.prep
)
oars.mpx.protein.lmer %>% arrange(pval)
oars.mpx.protein.lmer$feature = oars.mpx.protein.protein.names$good[match(oars.mpx.protein.lmer$taxa, oars.mpx.protein.protein.names$feature)]

# volcano plot
oars.mpx.protein.lmer.volcano = ggplot(oars.mpx.protein.lmer,
                                       aes(x=estimate, y=padj))+
  geom_point(shape=21, aes(fill=estimate), size=2)+
  geom_hline(yintercept=(0.20), linetype=2, alpha=0.5)+
  geom_vline(xintercept=c(-1,1), linetype=2, alpha=0.5)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.05,0.1, 0.20, 0.5, 1))+
  ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.2, taxa, NA)), size=3)+
  scale_fill_gradient2(low="blue", high="red")+
  facet_wrap(~"Interaction: Proteins")+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=10))+
  labs(x="Interaction Coefficient", y="FDR")

oars.mpx.protein.lmer.heatmap = ggplot(subset(oars.mpx.protein.lmer,padj<0.20),
                                       aes(x=1, y=reorder(feature, estimate)))+
  geom_tile(aes(fill=estimate), color="black")+
  scale_fill_gradient2(low="blue", high="red")+
  theme_minimal()+theme(axis.text.x = element_blank(),
                        plot.title = element_text(hjust = 0.5, size=12))+#theme(legend.position="none")+
  labs(x="", y="", title="MPX", fill="adjusted\ninteraction\ncoefficient")
oars.mpx.protein.lmer.heatmap


## Repeat with split
oars.mpx.protein.lmer.split = delta.omic.lmer(
  split = TRUE,
  data = oars.mpx.protein.prep#[,c(oars.asv.data.glom.lasso$feature, covars)]
)
oars.mpx.protein.lmer.split$feature = oars.mpx.protein.protein.names$good[match(oars.mpx.protein.lmer.split$taxa, oars.mpx.protein.protein.names$feature)]


oars.mpx.protein.lmer.split.volcano = ggplot(oars.mpx.protein.lmer.split  %>%
                                               mutate(Response = ifelse(response == "high", "Strong", "Weak")),
                                             aes(x=estimate, y=padj))+
  geom_point(shape=21, aes(fill=estimate), size=2)+
  geom_hline(yintercept=0.20, linetype=2, alpha=0.5)+
  geom_vline(xintercept=c(-1,1), linetype=2, alpha=0.5)+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.05,0.1, 0.20, 0.5, 1))+
  ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, taxa, NA)), size=3)+
  scale_fill_gradient2(low="blue", high="red")+
  facet_wrap(~Response)+
  theme_minimal()+theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="Adjusted Coefficient", y="FDR", title="Proteins")
oars.mpx.protein.lmer.split.volcano


# :: TEMPTED // defunct --------------------------------------------------------------

# Goal: Use TEMPTED to collapse patients to a point in Tensor Decomposition space
# then, identify features that most strongly contribute to variation in this space

# then, plot these trajectories to showcase variation

# first, process data (make % (not necessary) > remove not present > log+pseudo)
oars.tempted.data.asv = oars.asv.data.glom / 50000
oars.tempted.data.asv = oars.tempted.data.asv[,colSums(oars.tempted.data.asv)!=0]
oars.tempted.data.asv = log2(oars.tempted.data.asv+(min(oars.tempted.data.asv[oars.tempted.data.asv!=0])/2))

oars.tempted.data.mgx = oars.mgx.taxa / 100
oars.tempted.data.mgx = oars.tempted.data.mgx[,colSums(oars.tempted.data.mgx)!=0]
oars.tempted.data.mgx = log2(oars.tempted.data.mgx+(min(oars.tempted.data.mgx[oars.tempted.data.mgx!=0])/2))

oars.tempted.data.kegg = oars.mpx.kegg.mat
oars.tempted.data.kegg[is.na(oars.tempted.data.kegg)] = 0
oars.tempted.data.kegg = oars.tempted.data.kegg[,colSums(oars.tempted.data.kegg)!=0]
oars.tempted.data.kegg = log2(oars.tempted.data.kegg+(min(oars.tempted.data.kegg[oars.tempted.data.kegg!=0])/2))

oars.tempted.data.cog = oars.mpx.cog.mat
oars.tempted.data.cog[is.na(oars.tempted.data.cog)] = 0
oars.tempted.data.cog = oars.tempted.data.cog[,colSums(oars.tempted.data.cog)!=0]
oars.tempted.data.cog = log2(oars.tempted.data.cog+(min(oars.tempted.data.cog[oars.tempted.data.cog!=0])/2))

oars.tempted.data.cazy = oars.mpx.cazy.mat
oars.tempted.data.cazy[is.na(oars.tempted.data.cazy)] = 0
oars.tempted.data.cazy = oars.tempted.data.cazy[,colSums(oars.tempted.data.cazy)!=0]
oars.tempted.data.cazy = log2(oars.tempted.data.cazy+(min(oars.tempted.data.cazy[oars.tempted.data.cazy!=0])/2))

oars.tempted.data.mbx = oars.mbx.raw.mat.filt.1
oars.tempted.data.mbx = oars.tempted.data.mbx[,colSums(oars.tempted.data.mbx)!=0]
oars.tempted.data.mbx = log2(oars.tempted.data.mbx+(min(oars.tempted.data.mbx[oars.tempted.data.mbx!=0])/2))

# combine
oars.tempted.data = merge(oars.tempted.data.asv,oars.tempted.data.mgx,by=0)
rownames(oars.tempted.data) = oars.tempted.data$Row.names
oars.tempted.data$Row.names = NULL
oars.tempted.data = merge(oars.tempted.data,oars.tempted.data.kegg,by=0)
rownames(oars.tempted.data) = oars.tempted.data$Row.names
oars.tempted.data$Row.names = NULL
oars.tempted.data = merge(oars.tempted.data,oars.tempted.data.cog,by=0)
rownames(oars.tempted.data) = oars.tempted.data$Row.names
oars.tempted.data$Row.names = NULL
oars.tempted.data = merge(oars.tempted.data,oars.tempted.data.cazy,by=0)
rownames(oars.tempted.data) = oars.tempted.data$Row.names
oars.tempted.data$Row.names = NULL
oars.tempted.data = merge(oars.tempted.data,oars.tempted.data.mbx,by=0)
rownames(oars.tempted.data) = oars.tempted.data$Row.names
oars.tempted.data$Row.names = NULL
# 

# metadata
metadata.oars.stool.asv

# remove non.compliant (and sample with only 1 matching time point)
oars.tempted.data = oars.tempted.data[rownames(oars.tempted.data) %in% subset(metadata.oars.stool.asv, HM != "HM0618" & 
                                                                                oars.on.rs %in% c("preRS", "onRS") &
                                                                                compliant == T)$standard.name,]
# reorder
oars.tempted.data = oars.tempted.data[subset(metadata.oars.stool.asv, HM != "HM0618" & 
                                               oars.on.rs %in% c("preRS", "onRS") &
                                               compliant == T & 
                                               standard.name %in% rownames(oars.tempted.data))$standard.name,]
# n = 42 samples

# remove 0 var variables
oars.tempted.data = oars.tempted.data[,apply(oars.tempted.data, 2, sd)>0]
# scale
oars.tempted.data = scale(oars.tempted.data)

# check dim
dim(oars.tempted.data)
# check range
range(oars.tempted.data)

# make metadata
oars.tempted.meta = subset(metadata.oars.stool.asv, HM != "HM0618" & 
                             compliant == T & 
                             oars.on.rs %in% c("preRS", "onRS") &
                             standard.name %in% rownames(oars.tempted.data))

# TEMPTED
set.seed(25)
res_count <- tempted::tempted_all(oars.tempted.data,
                                  as.numeric(oars.tempted.meta$oars.days),
                                  oars.tempted.meta$HM,
                                  threshold=0.95,
                                  transform="none",
                                  #pseudo=0.5,
                                  r=2,
                                  smooth=1e-5,
                                  pct_ratio=0.1,
                                  pct_aggregate=1)
# plot
res_count_pca = as.data.frame(res_count$A_hat) %>% mutate(HM = rownames(.))

res_count_pca.plot = ggplot(data=res_count_pca, 
                            aes(x=PC1, y=PC2, fill=HM)) + 
  geom_point(shape=21, size=3) +
  ggnetwork::geom_nodetext_repel(aes(label=HM))+
  labs(x='Component 1', y='Component 2')+
  facet_wrap(~"Participant Loading")+
  theme_classic()+theme(legend.position="none",
                        strip.text=element_text(size=10))+
  labs(x=paste0("Component 1 (", round(res_count$r_square[1], digits=3)*100, "%)", sep=""),
       y=paste0("Component 2 (", round(res_count$r_square[2], digits=3)*100, "%)", sep=""))
res_count_pca.plot
# from baseline through treatment period
# "individualized effects of RS treatment"

tempted::plot_time_loading(res_count, r=2) + 
  geom_line(size=1.5) + 
  labs( x='Days on RS', y= "Loading")+
  facet_wrap(~"Temporal Loading")+
  theme_classic()+theme(#legend.position="none",
    strip.text=element_text(size=10))

# plot feature loadings 
tempted_feature_loadings = as.data.frame(res_count$B_hat) %>%
  mutate(feature = rownames(.)) %>%
  mutate(data.type = ifelse(feature %in% colnames(oars.tempted.data.asv), "ASV",
                            ifelse(feature %in% colnames(oars.tempted.data.mgx), "Species",
                                   ifelse(feature %in% colnames(oars.tempted.data.kegg), "Pathway",
                                          ifelse(feature %in% colnames(oars.tempted.data.cog), "COG", 
                                                 ifelse(feature %in% colnames(oars.tempted.data.cazy), "CAZy", 
                                                        ifelse(feature %in% colnames(oars.tempted.data.mbx), "Metabolite", "other"))))))) %>%
  mutate(data.type = factor(data.type, levels=c("ASV", "Species", "Pathway", "COG", "CAZy", "Metabolite")))


tempted_feature_loadings.plot = ggplot(tempted_feature_loadings, 
                                       aes(x=PC1, y=PC2)) + 
  geom_point(aes(shape = data.type, fill=data.type)) + 
  scale_fill_manual(values = omics.colors)+
  scale_shape_manual(values = omics.shapes)+
  geom_hline(yintercept=0, linetype=2)+
  geom_vline(xintercept=0, linetype=2)+
  #facet_wrap(~"Feature Loading")+
  facet_wrap(~data.type)+
  theme_classic()+theme(legend.position="none",
                        strip.text=element_text(size=10))+
  labs(x=paste0("Component 1 (", round(res_count$r_square[1], digits=3)*100, "%)", sep=""),
       y=paste0("Component 2 (", round(res_count$r_square[2], digits=3)*100, "%)", sep=""))

res_count_pca.plot+tempted_feature_loadings.plot

#. check features
res_count$toppct_ratio %>% as.data.frame() %>% 
  mutate(feature = rownames(.)) %>% subset(PC1 == TRUE) %>% head(n=20)

# check features manually
tempted_feature_loadings %>% slice_max(order_by = PC1, n=20)
tempted_feature_loadings %>% slice_max(order_by = PC2, n=20)


# I don't know how to establish signal from noise at individual level
# (I also don't find individualized analyses interesting)


# Stop, and check if dendrograms agree with Tempted plot


# >> 2.5. // alternative ----------------------------------------------------------
# semi-defunct

# Goal: use clinical outcome (fecal calprotectin) to classify responders
oars.fcal.vs.ph.cor.plot = metadata.oars.stool.double %>%
  #subset(HM != "HM0844" | timing != "6M") %>%
  subset(reltiming == "post") %>%
  ggplot(aes(x=lfc.fcal, y=delta.ph))+
  geom_point(shape=21, aes(fill=rs.col), size=2.5)+
  scale_fill_manual(values=c(labelcolors$cols[c(1,1,4,5,5,7,7,9)]))+
  geom_smooth(method="lm", color="black")+
  #geom_text(aes(label=paste(HM, timing)), size=3)+
  xlim(c(-6.1, 6))+
  ggpubr::stat_cor(method="spearman", size=3)+
  theme_classic()+theme(legend.position="none")+
  labs(x="Log2FC Fecal Calprotectin",
       y="Δ pH")
oars.fcal.vs.ph.cor.plot

oars.fcal.vs.ph.dist.plot = metadata.oars.stool.double %>%
  subset(reltiming == "post") %>%
  ggplot(aes(x=lfc.fcal, y=1))+
  ggridges::geom_density_ridges2(color="white")+
  xlim(c(-6.1, 6))+
  theme_classic()+theme(legend.position="none",
                        axis.title.x=element_blank())+
  labs(x="Log2FC Fecal Calprotectin",
       y="Density")
oars.fcal.vs.ph.dist.plot

oars.fcal.vs.ph.ttest.plot = metadata.oars.stool.double %>%
  #subset(HM != "HM0883" | timing != "6M") %>%
  subset(reltiming == "post") %>%
  mutate(fcal.sign = ifelse(lfc.fcal < 0, "Decrease", "Increase")) %>%
  ggplot(aes(x=fcal.sign, y=delta.ph))+
  geom_boxplot(width=0.5)+
  ggbeeswarm::geom_beeswarm(shape=21, aes(fill=rs.col), size=2.5, cex=5)+
  #geom_text(aes(label=paste(HM, timing, sep="_")))+
  scale_fill_manual(values=labelcolors$cols[c(1,1,4,5,5,7,7,9)])+
  ggpubr::stat_compare_means(method = "wilcox.test", size=3, label.x.npc = "right", hjust=1.15)+
  theme_classic()+theme(legend.position="none")+
  labs(x="Fecal Calprotectin Response",
       y="Δ pH")
oars.fcal.vs.ph.ttest.plot

(oars.fcal.vs.ph.dist.plot+
    oars.fcal.vs.ph.ttest.plot+
    oars.fcal.vs.ph.cor.plot+
    patchwork::plot_layout(nrow=2)) %>%
  ggsave(filename="./oars_plots/oars_fcal_responders.pdf",
         width=8, height=6, device = cairo_pdf)


# linear model
lmerTest::lmer(scale(delta.ph) ~ lfc.fcal + phase + diagnosis + adj.fiber + (1|HM),
               metadata.oars.stool.double %>%
                 subset(reltiming == "post")%>%
                 mutate(fcal.sign = ifelse(lfc.fcal < 0, "Decrease", "Increase")) )%>%
  summary() %>% coef()
lmerTest::lmer((delta.ph) ~ fcal.sign + phase + diagnosis + adj.fiber + (1|HM),
               metadata.oars.stool.double %>%
                 subset(reltiming == "post")%>%
                 mutate(fcal.sign = ifelse(lfc.fcal < 0, "Decrease", "Increase")) )%>%
  summary() %>% coef()


# // defunct --------------------------------------------------------------


adj_matrix <- precision_matrix
diag(adj_matrix) <- 0  # Remove self-loops
adj_matrix.df = reshape2::melt(adj_matrix)
adj_matrix.df = subset(adj_matrix.df, value != 0)
colnames(adj_matrix.df)[3] = "lasso"
# add original Spearman correlation
adj_matrix.df = merge(adj_matrix.df[,c("Var1", "Var2", "lasso")],
                      features.cor.fcal.df[,c("Var1", "Var2", "value")],
                      by=c("Var1", "Var2"))


# Convert to adjacency matrix
adj_matrix <- precision_matrix
adj_matrix[abs(adj_matrix) < 2.5e-3] <- 0  # Threshold small values
#adj_matrix <- abs(adj_matrix)  # Use absolute values for edges
diag(adj_matrix) <- 0  # Remove self-loops
pheatmap::pheatmap(adj_matrix,
                   color=colorRampPalette(c("blue","white", "red"))(100),
                   breaks=c(seq(min(adj_matrix), 0, length.out=ceiling(100/2) + 1), 
                            seq(max(adj_matrix)/100, max(adj_matrix), length.out=floor(100/2))))
# melt into df
adj_matrix.df = reshape2::melt(adj_matrix)
adj_matrix.df = subset(adj_matrix.df, value != 0)

# add original Spearman correlation
adj_matrix.df = merge(adj_matrix.df[,c("Var1", "Var2")],
                      metadata.oars.stool.omics.cor.sig[,c("Var1", "Var2", "value")],
                      by=c("Var1", "Var2"))

# Step 4: Create igraph object
# Only include nodes with significant edges
nodes <- unique(c(adj_matrix.df$Var1, adj_matrix.df$Var2))
g <- igraph::graph_from_data_frame(adj_matrix.df[, c("Var1", "Var2", "value")], 
                                   directed = FALSE, 
                                   vertices = nodes)

# Step 5: Visualize with ggnetwork
set.seed(123)  # For reproducible layout
n <- ggnetwork::ggnetwork(g, )  # Fruchterman-Reingold layout

# Plot network
library("ggnetwork")
# clean names + add omics signifier
n
n$name = ifelse(n$name == "fcal", "Fecal calprotectin", n$name)
n = n %>%
  mutate(data.type = ifelse(name %in% colnames(oars.mgx.taxa), "Species",
                            ifelse(name %in% colnames(oars.mpx.kegg.mat), "Pathway", 
                                   ifelse(name %in% colnames(oars.mpx.cog.mat), "COG", 
                                          ifelse(name %in% colnames(oars.mpx.cazy.mat), "CAZy", 
                                                 ifelse(name %in% colnames(oars.mbx.raw.mat.filt), "Metabolite",
                                                        ifelse(name %in% colnames(oars.asv.data.glom), "ASV", "Fecal calprotectin")))))))
n$omic.type = ifelse(grepl(paste(c("COG", "Pathway", "CAZy"), collapse="|"), n$data.type), "MPX", n$data.type)

n$data.type = factor(n$data.type, levels=c("Fecal calprotectin", "ASV", "Species", "Pathway", "COG", "CAZy", "Metabolite")) 
n$omic.type = factor(n$omic.type, levels=c("Fecal calprotectin", "ASV", "Species", "MPX", "Metabolite"))

library("ggnetwork")
metadata.oars.stool.omics.cor.network = ggplot(n, aes(x, y, xend = xend, yend = yend)) +
  geom_edges(aes(color = (value), size = abs(value)), alpha = 0.2) +  # Edge color by correlation, size by |corr|
  geom_nodes(size = 4, aes(shape=omic.type), fill="white", color = "white", alpha=1) +
  geom_nodes(size = 4, aes(shape=omic.type, fill=data.type), color = "black", alpha=0.9) +
  scale_shape_manual(values=c(21,22,23,24))+
  guides(shape = FALSE, fill=FALSE, size = FALSE)+
  scale_fill_manual(values= c("white", RColorBrewer::brewer.pal(n = 5, name = "Set3")[c(1,2,3,4,5)]))+
  geom_nodetext_repel(aes(label = gsub("\\/.*", "", gsub("\\,.*", "", name))), size = 2.5) +  # Node labels
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                        breaks = c(-0.5, 0, 0.5)) +  # Color for pos/neg correlations
  scale_size(range = c(0.5, 2)) +  # Edge thickness
  theme_blank() +  # ggnetwork's blank theme
  theme(legend.position="right")+
  labs(color = "Correlation",
       fill="Data type")
metadata.oars.stool.omics.cor.network


# build with UMAP instead (regularized spearman)
set.seed(25)
adj_matrix.umap = #metadata.oars.stool.omics.cor.sig %>% # or regularized with adj_matrix.df
  adj_matrix.df %>%
  # subset to strong
  #subset(abs(value) > 0.5) %>%
  # cycle through melt/cast to impute NA as 0
  reshape2::acast(Var1 ~ Var2, value.var="value") %>%
  reshape2::melt() %>%
  mutate(value = ifelse(is.na(value), 0, value)) %>% 
  #mutate(value = 1-value) %>% 
  reshape2::acast(Var1 ~ Var2, value.var="value") %>% as.matrix() %>%
  # make into correlation distance object
  #as.dist() %>% as.matrix() %>%
  # umap
  umap::umap()

# extract coords
adj_matrix.umap.df = adj_matrix.umap$layout       
adj_matrix.umap.df = adj_matrix.umap.df %>% as.data.frame() %>%
  mutate(feature = rownames(.))
# add edge data
adj_matrix.umap.df = do.call(rbind, lapply(1:nrow(adj_matrix.umap.df),function(x){
  data.subset = adj_matrix.umap.df[x,]
  # which features is this feature correlated to?
  cors = subset(adj_matrix.df, Var1 == data.subset$feature)
  colnames(cors) = c("feature", "Var2", "cor")
  data.subset = merge(data.subset, cors, by="feature")
  # add coords
  cors.coords = adj_matrix.umap.df
  colnames(cors.coords) = c("V1.B", "V2.B", "Var2")
  data.subset = merge(data.subset,
                      cors.coords, by="Var2")
  # good
  return(data.subset)
}))
adj_matrix.umap.df$feature = gsub("fcal", "Fecal calprotectin", adj_matrix.umap.df$feature)
adj_matrix.umap.df = adj_matrix.umap.df %>%
  mutate(data.type = ifelse(feature %in% colnames(oars.mgx.taxa), "Species",
                            ifelse(feature %in% colnames(oars.mpx.kegg.mat), "Pathway", 
                                   ifelse(feature %in% colnames(oars.mpx.cog.mat), "COG", 
                                          ifelse(feature %in% colnames(oars.mpx.cazy.mat), "CAZy", 
                                                 ifelse(feature %in% colnames(oars.mbx.raw.mat.filt), "Metabolite",
                                                        ifelse(feature %in% colnames(oars.asv.data.glom), "ASV", "Fecal calprotectin")))))))
adj_matrix.umap.df$omic.type = ifelse(grepl(paste(c("COG", "Pathway", "CAZy"), collapse="|"), adj_matrix.umap.df$data.type), "MPX", adj_matrix.umap.df$data.type)

adj_matrix.umap.df$data.type = factor(adj_matrix.umap.df$data.type, levels=c("Fecal calprotectin", "ASV", "Species", "Pathway", "COG", "CAZy", "Metabolite")) 

metadata.oars.stool.omics.cor.network.umap = ggplot(adj_matrix.umap.df, 
                                                    aes(V1, V2))+
  geom_segment(aes(x=V1, xend=V1.B, y=V2, yend=V2.B, color = (cor), size = abs(cor)), alpha = 0.2) +  # Edge color by correlation, size by |corr|
  geom_point(size = 4, aes(shape=data.type), fill="white", color = "white", alpha=1) +
  geom_point(size = 4, aes(shape=data.type, fill=data.type), color = "black", alpha=0.9) +
  scale_shape_manual(values=c("Fecal calprotectin" = 21, omics.shapes))+
  guides(shape = FALSE, fill=FALSE, size = FALSE)+
  scale_fill_manual(values= c("Fecal calprotectin" = "white", omics.colors))+
  geom_nodetext_repel(aes(label = gsub("\\/.*", "", gsub("\\,.*", "", feature))), size = 4) +  # Node labels
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                        breaks = c(-0.5, 0, 0.5)) +  # Color for pos/neg correlations
  scale_size(range = c(0.5, 2)) +  # Edge thickness
  theme_blank() +  # ggnetwork's blank theme
  theme(legend.position="right")+
  labs(color = "Correlation",
       fill="Data type")
metadata.oars.stool.omics.cor.network.umap


metadata.oars.stool.omics.cor.network


# :: RF Changes -----------------------------------------------

# Goal: Visualize if important features change with RS

# easier to do this manually
oars.feature.change.plot.1 = ggplot(oars.asv.data.glom.prep %>% 
                                      mutate(HM_phase = paste(HM, phase, sep="_")) %>%
                                      mutate(reltiming = factor(ifelse(reltiming == "pre", "Pre", "Post"), levels=c("Pre", "Post")))%>%
                                      mutate(Response = ifelse(response == "low", "Weak Response", "Strong Response")) %>%
                                      group_by(response, reltiming) %>%
                                      mutate(median.val = median(g__Escherichia_s__albertii)),
                                    aes(x=reltiming, y=median.val))+
  #geom_boxplot(width=0.5)+
  #scale_y_log10(labels = scales::label_number(accuracy = 0.01))+
  geom_line(aes(group=HM_phase, color=Response), linetype=2, alpha=0.5, linewidth=0.5)+
  geom_point(shape=21, aes(fill=Response), size=3)+
  geom_text(data = subset(oars.asv.data.glom.lmer, taxa == "g__Escherichia_s__albertii"), 
            x=1.5, y=0.1, vjust=1,
            aes(label=paste("p:", round(pval, digits=3))), size=3.5)+
  ylim(0.,0.1)+
  # scale_fill_manual(values=labelcolors$cols[c(1,1,4,5,5,7,7,9)])+
  facet_wrap(~"Interaction")+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=9),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="",
       y=" ")
oars.feature.change.plot.1


oars.feature.change.plot.2 = ggplot(oars.mpx.kegg.prep %>% 
                                      mutate(HM_phase = paste(HM, phase, sep="_")) %>%
                                      mutate(reltiming = factor(ifelse(reltiming == "pre", "Pre", "Post"), levels=c("Pre", "Post")))%>%
                                      mutate(Response = ifelse(response == "low", "Weak Response", "Strong Response")) %>%
                                      group_by(response, reltiming) %>%
                                      mutate(median.val = median(`Fluorobenzoate.degradation`))%>%
                                      dplyr::select(reltiming, HM_phase, Response, median.val) %>% distinct(),
                                    aes(x=reltiming, y=median.val))+
  #geom_boxplot(width=0.5)+
  #scale_y_log10(labels = scales::label_number(accuracy = 0.01))+
  geom_line(aes(group=HM_phase, color=Response), linetype=2, alpha=0.5, linewidth=0.5)+
  geom_point(shape=21, aes(fill=Response), size=3)+
  geom_text(data = subset(oars.mpx.kegg.lmer, taxa == "Fluorobenzoate.degradation"), 
            x=1.5, y=log10(923854), vjust=1,
            aes(label=paste("p:", round(pval, digits=3))), size=3.5)+ 
  scale_y_log10()+
  # scale_fill_manual(values=labelcolors$cols[c(1,1,4,5,5,7,7,9)])+
  facet_wrap(~"Interaction")+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=9),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="",
       y="Median Abundance")
oars.feature.change.plot.2

oars.feature.change.plot.3 = ggplot(oars.asv.data.glom.prep %>% 
                                      mutate(HM_phase = paste(HM, phase, sep="_")) %>%
                                      mutate(reltiming = factor(ifelse(reltiming == "pre", "Pre", "Post"), levels=c("Pre", "Post")))%>%
                                      mutate(Response = ifelse(response == "low", "Weak Response", "Strong Response")) %>%
                                      group_by(response, reltiming) %>%
                                      mutate(median.val = median(g__Bacteroides_H_857956_s__caccae)),
                                    aes(x=reltiming, y=median.val))+
  #geom_boxplot(width=0.5)+
  #scale_y_log10(labels = scales::label_number(accuracy = 0.01))+
  geom_line(aes(group=HM_phase, color=Response), linetype=2, alpha=0.5, linewidth=0.5)+
  geom_point(shape=21, aes(fill=Response), size=3)+
  geom_text(data = subset(oars.asv.data.glom.lmer, taxa == "g__Bacteroides_H_857956_s__caccae"), 
            x=1.5, y=1, vjust=1,
            aes(label=paste("p:", round(pval, digits=3))), size=3.5)+
  ylim(0,1)+
  # scale_fill_manual(values=labelcolors$cols[c(1,1,4,5,5,7,7,9)])+
  facet_wrap(~"Interaction")+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=9),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="",
       y="")
oars.feature.change.plot.3

# :: elastic net Omics ------------------------------------------------


# :: elastic net LOOCV ------------------------------------------------

# compare RF to elastic net (glmnet)

# Balance classes
# Reduce features via Spearman corr
# Train model using 5x CV
# Validate on LOOCV
# Repeat 15x

# for caret algo's [[takes ~2 hours]]
models = c("glmnet")

# test range of correlation thresholds (for feature reduction)
cor.thresholds = seq(from=0.3, to=1, by=0.1)

t1 = Sys.time()
oars.loocv.models.ml = 
  do.call(rbind, lapply(cor.thresholds, function(cor.thresh){
    do.call(rbind, lapply(models, function(model){
      do.call(rbind, lapply(1:15, function(iter){
        
        print(paste(model, iter, cor.thresh))
        # set data
        data.all = oars.all.omics.baseline
        # remove HM column
        data.all = data.all[,colnames(data.all) != "HM"]
        # dummy other vars (for xgbTree)
        data.all = data.all %>%
          mutate(phase = scale(as.numeric(as.factor(phase))),
                 response = as.factor(response),
                 diagnosis = scale(as.numeric(as.factor(diagnosis))))
        
        # log-scale vars (for non-scale invariant models)
        data.all[,!colnames(data.all) %in% c("phase", "response","diagnosis","adj.fiber")] = log2(data.all[,!colnames(data.all) %in% c("phase", "response","diagnosis","adj.fiber")])
        
        # run in parallel; loop through microbiomes (LOOCV)
        oars.loocv = do.call(rbind, parallel::mclapply(1:nrow(data.all), function(hm){
          print(paste(model, iter, hm))
          
          
          # reduce to train/test
          data.train = data.all[-hm, ]
          
          # balance data
          set.seed(iter)
          data.train = data.train %>%
            group_by(response) %>%
            sample_n(min(table(data.train$response)),replace=F)
          
          # remove features with no variance
          no_var_cols = apply(data.train[,!colnames(data.train) %in% c("phase", "response","diagnosis","adj.fiber")], 2, sd) %>% data.frame() %>% subset(. == 0)
          data.train = data.train[,!colnames(data.train) %in% rownames(no_var_cols)]
          
          # remove features with low variance
          low_var_cols <- caret::nearZeroVar(data.train[,!colnames(data.train) %in% c("phase", "response","diagnosis","adj.fiber")], 
                                             saveMetrics = TRUE)
          # aggressive; select features with >= 10% unique values
          data.train = data.train[,!colnames(data.train) %in% 
                                    rownames(subset(low_var_cols, percentUnique <10))]
          
          # feature reduction by correlation (cor is WAY faster than Hmisc::rcorr (no p vals))
          high_cor_cols <- caret::findCorrelation(cor(as.matrix(data.train[,!colnames(data.train) %in% c("phase", "response","diagnosis","adj.fiber")]), method="spearman"), 
                                                  cutoff = cor.thresh, names = TRUE)
          data.train = data.train[,!colnames(data.train) %in% c(high_cor_cols)]
          
          if(model %in% c("rf", "oob", "OOB", "RF")){
            # build model (OOB) - faster, but too many features --> overfits
            set.seed(iter)
            model.results = ranger::ranger(response ~ ., data.train, 
                                           importance = "impurity",
                                           probability = TRUE)
            
          }
          # or, build Caret model
          
          if(!model %in% c("rf", "oob", "OOB", "RF")){
            
            # build model
            set.seed(iter)
            model.results <- caret::train(
              x = data.train[,colnames(data.train) != "response"],
              y = data.train$response,
              method = model,
              metric = "ROC",
              verbosity = 0,
              verbose=F,
              trControl = caret::trainControl(method = "cv", number = 5,
                                              summaryFunction = caret::twoClassSummary,   # enables AUC
                                              classProbs = TRUE,    
                                              savePredictions = "final"),
              preProcess = c("center", "scale")) %>% suppressWarnings() %>% suppressMessages()
          }
          
          # apply model
          
          if(model %in% c("rf", "oob", "OOB", "RF")){
            # for rf oob
            pred = predict(model.results, 
                           data.all[hm, ])
            
            output = data.frame(pred = pred$predictions[1],
                                true = data.all[hm, ]$response,
                                index = hm,
                                iter = iter,
                                cor = cor.thresh,
                                model = "rf")
          } else {
            # for caret:
            pred = predict(model.results, 
                           data.all[hm, ],
                           type="prob")
            
            output = data.frame(pred = pred,
                                true = data.all[hm, ]$response,
                                index = hm,
                                iter = iter,
                                cor = cor.thresh,
                                model = model)
          }
          return(output)
        }))
      }))
    }))
  }))
t2 = Sys.time()
t2-t1 # ~40 min

# save this
#oars.loocv.models.df = oars.loocv.models.rf
#oars.loocv.models.df$model = rep(c("ranger", "glmnet", "xgbTree"), each = (15*19))

# which model is best?
oars.loocv.models.ml.plot = oars.loocv.models.ml %>%
  group_by(iter, model, cor) %>%
  summarize(auc = pROC::auc(true, pred.low)[1])%>%  
  ggplot(
    aes(x=as.factor(cor), y=auc))+
  geom_boxplot(width=0.5)+
  ylim(0.3, 1)+
  ggbeeswarm::geom_beeswarm(shape=21, fill="white")+
  theme_classic()+
  facet_wrap(~"Elastic Net")+
  labs(x="Correlation Threshold",
       y="LOOCV AUC")
oars.loocv.models.ml.plot


# calculate values and extract best
oars.loocv.models.ml.df = oars.loocv.models.ml %>%
  group_by(iter, model, cor) %>%
  summarize(auc = pROC::auc(true, pred.low)[1])%>%  
  group_by(model, cor) %>%
  mutate(mean.auc = mean(auc),
         median.auc = median(auc),
         auc.low = mean(auc) - (sd(auc)/sqrt(n()) * 1.96), # 95% CI
         auc.high = mean(auc) + (sd(auc)/sqrt(n()) * 1.96)) %>%
  dplyr::select(mean.auc, median.auc, auc.low, auc.high, model, cor) %>% distinct()

# extract best
oars.loocv.models.ml.best.stats = subset(oars.loocv.models.ml.df, mean.auc == max(oars.loocv.models.ml.df$mean.auc))
oars.loocv.models.ml.best = subset(oars.loocv.models.ml, cor == oars.loocv.models.ml.best.stats$cor)

# plot
oars.loocv.ml.roc = do.call(rbind, lapply(1:15, function(seed){
  data.subset = subset(oars.loocv.models.ml.best, iter == seed)
  data.frame(sens = pROC::roc(data.subset$true, data.subset$pred.low)$sensitivities,
             spec = pROC::roc(data.subset$true, data.subset$pred.low)$specificities,
             iter = seed)
}))


# Step 2: Compute mean and standard error for sensitivities across iterations
roc_summary.ml <- oars.loocv.ml.roc %>%
  group_by(spec) %>%  # Group by specificity (or alternatively by sens)
  summarise(
    mean_sens = mean(sens, na.rm = TRUE),
    se_sens = sd(sens, na.rm = TRUE) / sqrt(n()),  # Standard error
    lower = mean_sens - 1.96 * se_sens,           # 95% CI lower bound
    upper = mean_sens + 1.96 * se_sens            # 95% CI upper bound
  ) %>%
  mutate(fpr = 1 - spec) %>%  # False Positive Rate (1 - specificity)
  filter(!is.na(mean_sens) & !is.na(se_sens))  # Remove any NA values


oars.ml.loocv.ml.roc.plot = ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  # Individual ROC curves for each iteration
  #geom_line(data = oars.ml.loocv.roc %>% arrange(sens),
  #          aes(x=1-spec, y=sens, group=iter))+
  # Error ribbon (95% CI)
  geom_ribbon(data = roc_summary.ml, 
              aes(x = fpr, ymin = lower, ymax = upper), 
              fill = "black", alpha = 0.2) +
  # Mean ROC curve
  geom_path(data = roc_summary.ml, 
            aes(x = fpr, y = mean_sens), 
            color = "black", size = 1) +
  # add label
  annotate(geom="text", x=0.75, y=0.25,
           label=paste("AUC: ", round(oars.loocv.models.ml.best.stats$mean.auc, digits=2),
                       "\n(", round(oars.loocv.models.ml.best.stats$auc.low, digits=2), 
                       ", ", round(oars.loocv.models.ml.best.stats$auc.high,digits=2), ")", sep=""),
           size=4)+
  theme_classic()+theme(strip.text=element_text(size=10))+
  facet_wrap(~"Elastic Net 15x LOOCV")+
  xlim(0,1)+
  ylim(0,1)+
  labs(x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)") 
oars.ml.loocv.ml.roc.plot

# importances
oars.ml.loocv.ml.imp = do.call(rbind, lapply(1:15, function(iter){
  print(paste(iter))
  
  oars.ml.loocv = do.call(rbind, parallel::mclapply(1:nrow(oars.all.omics.baseline), function(hm){
    # reduce to train/test
    data.train = oars.all.omics.baseline[-hm, !colnames(oars.all.omics.baseline)%in%c("HM")]
    # balance data
    set.seed(iter)
    data.train = data.train %>%
      group_by(response) %>%
      sample_n(min(table(data.train$response)))
    
    
    # remove features with no variance
    no_var_cols = apply(data.train[,!colnames(data.train) %in% c("phase", "response","diagnosis","adj.fiber")], 2, sd) %>% data.frame() %>% subset(. == 0)
    data.train = data.train[,!colnames(data.train) %in% rownames(no_var_cols)]
    
    # remove features with low variance
    low_var_cols <- caret::nearZeroVar(data.train[,!colnames(data.train) %in% c("phase", "response","diagnosis","adj.fiber")], 
                                       saveMetrics = TRUE)
    # aggressive; select features with >= 10% unique values
    data.train = data.train[,!colnames(data.train) %in% 
                              rownames(subset(low_var_cols, percentUnique <10))]
    
    # feature reduction by correlation (cor is WAY faster than Hmisc::rcorr (no p vals))
    high_cor_cols <- caret::findCorrelation(cor(as.matrix(data.train[,!colnames(data.train) %in% c("phase", "response","diagnosis","adj.fiber")]), method="spearman"), 
                                            cutoff = oars.loocv.models.ml.best.stats$cor, names = TRUE)
    data.train = data.train[,!colnames(data.train) %in% c(high_cor_cols)]
    
    set.seed(iter)
    model.results <- caret::train(
      x = data.train[,colnames(data.train) != "response"],
      y = data.train$response,
      method = "glmnet",
      metric = "ROC",
      verbosity = 0,
      verbose=F,
      trControl = caret::trainControl(method = "cv", number = 5,
                                      summaryFunction = caret::twoClassSummary,   # enables AUC
                                      classProbs = TRUE,    
                                      savePredictions = "final"),
      preProcess = c("center", "scale")) %>% suppressWarnings() %>% suppressMessages()
    # extract feature importances
    imps = caret::varImp(model.results)[1] %>% data.frame()
    imps$feature = rownames(imps)
    
    # importance
    ml.imp = data.frame(imp = imps,
                        iter = iter,
                        index = hm) %>%
      mutate(feature = rownames(.))
    rownames(ml.imp) = NULL
    ml.imp
  }))
}))

oars.ml.loocv.ml.imp.df = oars.ml.loocv.ml.imp %>%
  #subset(imp != 0) %>%
  group_by(feature) %>%
  mutate(mean.imp = mean(na.omit(imp.Overall))) %>%
  mutate(imp.low = mean(imp.Overall)- (sd(imp.Overall)/sqrt(n()) * 1.96)) %>% # 95% CI
  mutate(imp.high = mean(imp.Overall)+ (sd(imp.Overall)/sqrt(n()) * 1.96)) %>% # 95% CI
  subset(mean.imp != 0) %>%
  dplyr::select(feature, mean.imp, imp.low, imp.high) %>% distinct() %>%
  arrange(-mean.imp) %>% data.frame()

# log2fc
oars.ml.loocv.ml.imp.df = oars.ml.loocv.ml.imp.df %>% 
  subset(!is.na(mean.imp))%>%
  arrange(mean.imp) %>% slice_max(mean.imp, n=15) 
# calculate wilcox
oars.ml.loocv.ml.imp.wilcox = do.call(rbind, lapply(oars.ml.loocv.ml.imp.df$feature, function(x){
  data.subset = oars.all.omics.baseline[,colnames(oars.all.omics.baseline) %in% c(x, "response")]
  # run wilcox
  wilcox.p = wilcox.test(subset(data.subset, response == "high")[,2],
                         subset(data.subset, response == "low")[,2])
  # calculate lfc
  if(x %in% c("stool_water_perc", "shannon", "fd")){
    coef = mean(subset(data.subset, response == "high")[,2]) / mean(subset(data.subset, response == "low")[,2])
  }else{
    coef = mean(log2(subset(data.subset, response == "high")[,2])) - mean(log2(subset(data.subset, response == "low")[,2]))
  }
  data.frame(feature = x,
             pval = wilcox.p$p.value,
             coef = coef)
}))
oars.ml.loocv.ml.imp.df = merge(oars.ml.loocv.ml.imp.wilcox,
                                oars.ml.loocv.ml.imp.df, by="feature")
oars.ml.loocv.ml.imp.df$sig = ifelse(oars.ml.loocv.ml.imp.df$pval < 0.05, "*", "")
oars.ml.loocv.ml.imp.df = oars.ml.loocv.ml.imp.df %>% arrange(mean.imp)

# add omic type
oars.ml.loocv.ml.imp.df = oars.ml.loocv.ml.imp.df %>%
  mutate(Datatype = ifelse(feature %in% colnames(oars.asv.baseline), "ASV",
                           ifelse(feature %in% colnames(oars.mgx.baseline), "Species",
                                  ifelse(feature %in% colnames(oars.mpx.kegg.baseline), "Pathway",
                                         ifelse(feature %in% colnames(oars.mpx.cazy.baseline), "CAZy", "other")))))

oars.ml.loocv.ml.imp.df$feature = gsub("\\.", " ", oars.ml.loocv.ml.imp.df$feature)
oars.ml.loocv.ml.imp.df$feature = gsub("   globo.*", "", oars.ml.loocv.ml.imp.df$feature)
oars.ml.loocv.ml.imp.df$feature = gsub("_857956", "", oars.ml.loocv.ml.imp.df$feature)
oars.ml.loocv.ml.imp.df$feature = gsub("Neomycin ", "Neomycin,", oars.ml.loocv.ml.imp.df$feature)


oars.ml.loocv.ml.imp.plot = ggplot(oars.ml.loocv.ml.imp.df %>%
                                     mutate(Datatype = factor(Datatype, levels=c("ASV", "MGX", "MPX (Pathway)", "MGX (CAZy)"))),
                                   aes(x=mean.imp, y=reorder(feature, mean.imp)))+
  geom_segment(aes(x=imp.low, xend=imp.high, y=feature, yend=feature), color="black")+
  geom_point(shape=21, aes(fill=(Datatype)), size=3.5)+
  #geom_point(shape=21, aes(fill=scale(imp)), size=3)+
  # geom_text(x=-Inf, aes(label = sig),  hjust=-1, size=5, vjust=0.75)+
  geom_text(data = oars.ml.loocv.ml.imp.df,
            aes(label=ifelse(mean.imp > 30 & mean.imp < 35, paste(" ", feature), NA), x=imp.high, y=reorder(feature, mean.imp)), size=2.5, hjust=0)+
  #geom_text(data = oars.ml.loocv.ml.imp.df,
  #          aes(label=ifelse(mean.imp < 35 & mean.imp >= 35, paste(feature," "), NA), x=imp.low, y=reorder(feature, mean.imp)), size=2.5, hjust=1)+
  geom_text(data = oars.ml.loocv.ml.imp.df,
            aes(label=ifelse(mean.imp >= 35, paste(feature, " "), NA), x=imp.low, y=reorder(feature, mean.imp)), size=2.5, hjust=1)+
  xlim(10,100)+
  #scale_fill_gradient2(low="blue", high="red")+
  scale_fill_manual(values=RColorBrewer::brewer.pal(n = 5, name = "Set3")[c(1,2,3,5)])+
  scale_color_manual(values=c("white", "black"))+
  theme_classic()+theme(legend.position=c(0.8,0.35),
                        axis.title.y = element_blank(),
                        axis.text.y = element_blank(),
                        axis.ticks.y=element_blank(),
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  guides(color=FALSE)+
  labs(x="Mean Scaled L12 Coefficient", fill="")+
  facet_wrap(~"Feature Importance")
oars.ml.loocv.ml.imp.plot


.# loop through top features and plot
oars.ml.loocv.ml.imp.ttest = do.call(rbind, lapply(slice_max(oars.ml.loocv.ml.imp.df, mean.imp, n=3)$feature, function(x){
  x = gsub(" ", "\\.", x)
  data.subset = oars.all.omics.baseline[,c(x, "response")]
  rownames(data.subset) = paste(oars.all.omics.baseline$HM, oars.all.omics.baseline$phase, sep="_")
  data.subset$feature = colnames(data.subset)[1]
  data.subset$sample = rownames(data.subset)
  colnames(data.subset)[1] = "value"
  data.subset$feature = gsub("\\.", " ", data.subset$feature)
  
  data.subset
}))
# add stats
oars.ml.loocv.ml.imp.ttest = merge(oars.ml.loocv.ml.imp.ttest,
                                   oars.ml.loocv.ml.imp.df, by="feature")
oars.ml.loocv.ml.imp.ttest.p = oars.ml.loocv.ml.imp.ttest[,c("feature", "pval")] %>% distinct()

# clean taxa names (so they fit)
oars.ml.loocv.ml.imp.ttest.p$feature = gsub("\\.", " ", gsub("s__", "", (gsub("857956_", "", oars.ml.loocv.ml.imp.ttest.p$feature))))
oars.ml.loocv.ml.imp.ttest$feature = gsub("\\.", " ", gsub("s__", "", (gsub("857956_", "", oars.ml.loocv.ml.imp.ttest$feature))))

oars.ml.loocv.ml.imp.ttest.plot = ggplot(oars.ml.loocv.ml.imp.ttest %>%
                                           # clean names
                                           mutate(Response = ifelse(response == "low", "Weak\nResponse", "Strong\nResponse")) %>%
                                           # add indicator of not present
                                           group_by(feature) %>%
                                           mutate(pseudo = ifelse(value == min(value), "pseudo", "real")) %>%
                                           # reorder taxa based on imp
                                           mutate(feature = factor(feature, levels=arrange(distinct(oars.ml.loocv.ml.imp.ttest[,c("feature", "mean.imp")]),-mean.imp)$feature)),
                                         aes(x=Response, y=value))+
  geom_boxplot(width=0.5, outlier.shape=NA)+
  #scale_y_log10(labels = scales::label_number(accuracy = 0.01))+
  #scale_y_log10()+
  ggbeeswarm::geom_beeswarm(shape=21, aes(fill=response, alpha=pseudo), size=3)+
  scale_alpha_manual(values=c(0.2, 1))+
  theme_classic()+theme(legend.position="none",
                        axis.title.x = element_blank(),
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=8),
                        strip.background = element_rect(
                          color="black"))+
  geom_text(data=oars.ml.loocv.ml.imp.ttest.p %>%
              # clean names
              mutate(feature = factor(feature, levels=arrange(distinct(oars.ml.loocv.ml.imp.ttest[,c("feature", "mean.imp")]),-mean.imp)$feature)),
            x=1.5, y=Inf, vjust=1.2, 
            aes(label = paste("p =", round(pval, digits=3))),
            size=3.5)+
  labs(x="", y="Feature Abundance")+
  facet_wrap(~feature, nrow=3, scales="free")
oars.ml.loocv.ml.imp.ttest.plot


# :: elastic net Changes -----------------------------------------------

# easier to do this manually
oars.feature.change.ml.plot.1 = ggplot(oars.mgx.prep %>% 
                                         mutate(HM_phase = paste(HM, phase, sep="_")) %>%
                                         mutate(reltiming = factor(ifelse(reltiming == "pre", "Pre", "Post"), levels=c("Pre", "Post")))%>%
                                         mutate(Response = ifelse(response == "low", "Weak Response", "Strong Response")) %>%
                                         group_by(response, reltiming) %>%
                                         mutate(median.val = median(Actinomyces_oris)),
                                       aes(x=reltiming, y=median.val))+
  #geom_boxplot(width=0.5)+
  #scale_y_log10(labels = scales::label_number(accuracy = 0.01))+
  geom_line(aes(group=HM_phase, color=Response), linetype=2, alpha=0.5, linewidth=0.5)+
  geom_point(shape=21, aes(fill=Response), size=3)+
  geom_text(data = subset(oars.mgx.lmer, taxa == "Actinomyces_oris"), 
            x=1.5, y=Inf, vjust=1.2,
            aes(label=paste("p:", round(pval, digits=3))), size=3.5)+
  ylim(0.,0.01)+
  # scale_fill_manual(values=labelcolors$cols[c(1,1,4,5,5,7,7,9)])+
  #facet_wrap(~"Interaction")+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=9),
                        axis.title.y=element_blank(),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="",
       y="")
oars.feature.change.ml.plot.1


oars.feature.change.ml.plot.2 = ggplot(oars.asv.data.glom.prep %>% 
                                         mutate(HM_phase = paste(HM, phase, sep="_")) %>%
                                         mutate(reltiming = factor(ifelse(reltiming == "pre", "Pre", "Post"), levels=c("Pre", "Post")))%>%
                                         mutate(Response = ifelse(response == "low", "Weak Response", "Strong Response")) %>%
                                         group_by(response, reltiming) %>%
                                         mutate(median.val = median(`g__Bacteroides_H_857956_s__caccae`))%>%
                                         dplyr::select(reltiming, HM_phase, Response, median.val) %>% distinct(),
                                       aes(x=reltiming, y=median.val))+
  #geom_boxplot(width=0.5)+
  #scale_y_log10(labels = scales::label_number(accuracy = 0.01))+
  geom_line(aes(group=HM_phase, color=Response), linetype=2, alpha=0.5, linewidth=0.5)+
  geom_point(shape=21, aes(fill=Response), size=3)+
  geom_text(data = subset(oars.asv.data.glom.lmer, taxa == "g__Bacteroides_H_857956_s__caccae"), 
            x=1.5, y=Inf, vjust=1.2,
            aes(label=paste("p:", round(pval, digits=3))), size=3.5)+ 
  ylim(0,0.85)+
  # scale_fill_manual(values=labelcolors$cols[c(1,1,4,5,5,7,7,9)])+
  #facet_wrap(~"Interaction")+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=9),
                        #axis.title.y=element_blank(),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="",
       y="Median Abundance")
oars.feature.change.ml.plot.2

oars.feature.change.ml.plot.3 = ggplot(oars.mpx.kegg.prep %>% 
                                         mutate(HM_phase = paste(HM, phase, sep="_")) %>%
                                         mutate(reltiming = factor(ifelse(reltiming == "pre", "Pre", "Post"), levels=c("Pre", "Post")))%>%
                                         mutate(Response = ifelse(response == "low", "Weak Response", "Strong Response")) %>%
                                         group_by(response, reltiming) %>%
                                         mutate(median.val = median(`Biosynthesis.of.unsaturated.fatty.acids`)),
                                       aes(x=reltiming, y=median.val))+
  #geom_boxplot(width=0.5)+
  #scale_y_log10(labels = scales::label_number(accuracy = 0.01))+
  geom_line(aes(group=HM_phase, color=Response), linetype=2, alpha=0.5, linewidth=0.5)+
  geom_point(shape=21, aes(fill=Response), size=3)+
  geom_text(data = subset(oars.mpx.kegg.lmer, taxa == "Biosynthesis.of.unsaturated.fatty.acids"), 
            x=1.5, y=Inf, vjust=1.2,
            aes(label=paste("p:", round(pval, digits=3))), size=3.5)+
  #geom_text(aes(label=median.val))+
  ylim(50000000, 84000000)+
  # scale_fill_manual(values=labelcolors$cols[c(1,1,4,5,5,7,7,9)])+
  #facet_wrap(~"Interaction")+
  theme_classic()+theme(legend.position="none",
                        strip.text = element_text(size=9),
                        axis.title.y=element_blank(),
                        strip.background = element_rect(
                          color="black"))+
  labs(x="",
       y="")
oars.feature.change.ml.plot.3



oars.ml.loocv.ml.imp.ttest.plot|
  (oars.feature.change.ml.plot.1/oars.feature.change.ml.plot.2/oars.feature.change.ml.plot.3)


# :: // defunct RF response tree ------------------------------------------


## RF does not work for this
set.seed(25)
rf_model <- randomForestSRC::rfsrc(
  starch.lfc ~., 
  oars.response.tree.data[rownames(oars.response.tree.data)%in%subset(metadata.oars.stool.asv, phase=="treatment")$standard.name,],
  proximity=T
)

# Proximity
rf_prox = ape::pcoa(dist(1-rf_model$proximity))$vectors[,c(1:2)] %>% as.data.frame()
rf_prox$standard.name <- rownames(oars.response.tree.data[rownames(oars.response.tree.data)%in%subset(metadata.oars.stool.asv, phase=="treatment")$standard.name,])
rf_prox$starch.lfc <- oars.response.tree.data[rownames(oars.response.tree.data)%in%subset(metadata.oars.stool.asv, phase=="treatment")$standard.name,]$starch.lfc

# add metadata
rf_prox = merge(rf_prox,
                metadata.oars.stool.asv, by="standard.name")

# plots
ggplot(rf_prox, aes(x = Axis.1, y = Axis.2, fill = rs.col)) +
  geom_point(shape=21, size=3) +
  #ggrepel::geom_text_repel(aes(label = standard.name),size = 2.5) +
  geom_text(aes(label = standard.name),size = 2.5) +
  
  #scale_fill_gradient2(low="blue", high="red")+
  scale_fill_manual(values=rs.colors, na.value="lightgrey")+
  theme_classic() +theme(legend.position="none")+
  labs(
    x = "Axis.1",
    y = "Axis.2"
  )
# correlate with starch
cor.test(scores$starch.lfc, scores$p1, method="spearman") # strong correlation
cor.test(scores$starch.lfc, scores$ortho, method="spearman") # no correlation
# Which features are associated with Axis 1 and 2?


# :: //defunct Response Tree oPLS --------------------------------------------------------

# note: uses data from Fcal (LFC), so position here

# Premise: RS must first be degraded (A) before it can be fermented (B).
# This (A) requires CAZymes.
# So, to control for the fact that some patients might not
# degrade the RS (A), it's unlikely they could ferment it (B).
# So, let's control for these non-digesters and ask:
# "which features correlate with starch CAZy"

# Side note: Ex vivo fermentation may not predict this, because
# fermentation requires both A AND B. A is not sufficient.
# (i.e. LetsDoOrganic may have been degraded, but not fermented)

# there are several approaches we could take:
# 1. oPLS --> control for starch-cazy-lfc, and see which variables
#             vary that are non-correlated to starch-cazy-lfc
# 2. Correlation --> features with LFC that correlate with starch-cazy-lfc
#                    are probably downstream of degradation (A)
#                    i.e. fermentation (B)

# Note: RF Proximity positions "off RS", and "strong RS response"
# along same position in Axis.2, which negates this analysis. Don't use.

# Starch CAZy data
oars.response.tree.meta = metadata.oars.stool.asv[,c("standard.name","HM","phase","timing", "compliant", "starch")] %>%
  subset(!is.na(HM) & compliant == T & !is.na(starch))  %>%
  mutate(starch.lfc = log2(starch / lag(starch))) %>%
  subset(phase == "treatment" & timing !="0M") %>%
  dplyr::select(starch.lfc)
nrow(oars.response.tree.meta)

# Omics LFC data
oars.response.tree.data = oars.fecalcal.omics.lfc.saved#[,colnames(oars.fecalcal.omics.lfc.saved)!="fcal"]
oars.response.tree.data = as.data.frame(oars.response.tree.data)
oars.response.tree.data[is.na(oars.response.tree.data)] = 0
# and subset to compliant patients, because now we want to 
# find associations between Starch CAZy and other features

# merge with starch
oars.response.tree.data = merge(oars.response.tree.data,
                                oars.response.tree.meta%>%
                                  dplyr::select(starch.lfc),
                                by="row.names")
rownames(oars.response.tree.data) = oars.response.tree.data$Row.names
oars.response.tree.data$Row.names = NULL

# sanity check
cor.test(oars.response.tree.data$GH13,
         oars.response.tree.data$CBM48)
cor.test(oars.response.tree.data$GH13,
         oars.response.tree.data$starch.lfc)
# GOOD

# orthogonal PLS (control for Starch CAZy Log2FC, 
# and find variables that explain variance across patients)

opls_model <- ropls::opls(
  x = oars.response.tree.data[,colnames(oars.response.tree.data) != "starch.lfc"],
  y = oars.response.tree.data$starch.lfc,
  predI = 1,         # Number of predictive components (start with 1)
  orthoI = 1,        # Number of orthogonal components (adjust as needed)
  crossvalI = 5,     # 5-fold cross-validation
  scaleC = "none" # Scale X to unit variance
)
opls_model

# Scores (patient samples)
scores <- as.data.frame(ropls::getScoreMN(opls_model))
scores$ortho = opls_model@orthoScoreMN
scores$standard.name <- rownames(oars.response.tree.data)
scores$starch.lfc <- oars.response.tree.data$starch.lfc
# add metadata
scores = merge(scores,
               metadata.oars.stool.asv, by="standard.name")

# Loadings (features)
loadings <- as.data.frame(ropls::getLoadingMN(opls_model))
loadings$ortho = opls_model@orthoLoadingMN
loadings$feature <- colnames(oars.response.tree.data[,colnames(oars.response.tree.data) != "starch.lfc"])
# find strongest variables
loadings

# plots
oars.response.tree.plot = ggplot(scores, aes(x = p1, y = ortho, fill = rs.col)) +
  geom_point(shape=21, size=3) +
  ggrepel::geom_text_repel(aes(label = paste(HM, timing, sep=" ")),size = 2.5) +
  #scale_fill_gradient2(low="blue", high="red")+
  scale_fill_manual(values=rs.colors, na.value="lightgrey")+
  theme_classic() +theme(legend.position="none")+
  labs(
    x = "Predictive Score",
    y = "Orthogonal Score"
  )
oars.response.tree.plot
# correlate with starch
cor.test(scores$starch.lfc, scores$p1, method="spearman") # strong correlation
cor.test(scores$starch.lfc, scores$ortho, method="spearman") # no correlation
# Which features are associated with Axis 1 and 2?

# plot correlation
oars.response.tree.cor.plot = ggplot(scores, aes(x=p1, y=starch.lfc))+
  geom_point(shape=21, aes(fill=rs.col), size=3)+
  geom_smooth(method="lm", color="black")+
  scale_fill_manual(values=rs.colors, na.value="lightgrey")+
  ggpubr::stat_cor(method="spearman")+
  geom_hline(yintercept=0, linetype=2, alpha=0.5)+
  theme_classic()+
  theme(legend.position="none",
        axis.title.x=element_blank())+
  labs(x="Predictive Score",
       y="Starch CAZy Log2FC")
oars.response.tree.cor.plot

# identify convex hull points
hull_indices <- chull(loadings$p1, loadings$o1)
loadings$hull = ifelse(loadings$feature%in%loadings$feature[hull_indices],"hull","not_hull")

# Compute 2D kernel density
kde <- MASS::kde2d(loadings$p1, loadings$ortho, n = 100)
contour_level <- quantile(kde$z, probs = 0.95)  # 95% density contour

# Extract contour points
contour_lines <- contourLines(kde$x, kde$y, kde$z, 
                              levels = contour_level)
boundary_kde <- as.data.frame(contour_lines[[1]][c("x", "y")])
colnames(boundary_kde) = c("p1","ortho")
# Convert boundary to a matrix for point-in-polygon test
boundary_matrix <- as.matrix(boundary_kde)
# Test which points are inside the contour
inside <- sp::point.in.polygon(loadings$p1, loadings$ortho, boundary_matrix[, 1], boundary_matrix[, 2])
loadings$ahull = inside
table(inside)

ggplot(loadings, aes(x = p1, y = ortho)) +
  geom_point(alpha = 0.7, shape=21, size=2) +
  ggpubr::stat_chull(data=boundary_kde)+
  geom_nodetext_repel(aes(label = ifelse(ahull==0, feature, NA)), vjust = 1.5, size = 3) +
  #scale_color_manual()+
  theme_classic() +theme(legend.position="none")+
  labs(
    x = "Predictive Score",
    y = "Orthogonal Score"
  )

# Let's run Correlations instead

# features to test:
oars.response.tree.features = colnames(oars.response.tree.data.pls[,1:(ncol(oars.response.tree.data.pls)-3)])
# reduce to 95% ahull
oars.response.tree.features = subset(loadings, ahull==0)$feature

# link sample LFC to scores
oars.response.tree.data.pls = oars.response.tree.data
oars.response.tree.data.pls$p1 = scores$p1
oars.response.tree.data.pls$ortho = scores$ortho

oars.response.tree.pls.cor.1 = do.call(rbind, lapply(oars.response.tree.features, function(x){
  data.subset = oars.response.tree.data.pls[,colnames(oars.response.tree.data.pls) %in% c(x, "p1")]
  cor.results = cor.test(data.subset[,1], data.subset[,2], method="spearman")
  cor.results = data.frame(cor = cor.results$estimate,
                           pval = cor.results$p.value,
                           feature = x)
})) %>%
  mutate(padj = p.adjust(pval, method="BH")) %>%
  arrange(padj)
oars.response.tree.pls.cor.1 %>% head(n=10)
# Highest covariance with Starch CAZy; no surprise, mostly Starch-degrading proteins

oars.response.tree.pls.cor.2 = do.call(rbind, lapply(oars.response.tree.features, function(x){
  print(x)
  data.subset = oars.response.tree.data.pls[,colnames(oars.response.tree.data.pls) %in% c(x, "ortho")]
  cor.results = cor.test(data.subset[,1], data.subset[,2], method="spearman")
  cor.results = data.frame(cor = cor.results$estimate,
                           pval = cor.results$p.value,
                           feature = x)
})) %>%
  mutate(padj = p.adjust(pval, method="BH")) %>%
  arrange(padj)
oars.response.tree.pls.cor.2 %>% head(n=30)
# Highest variation across individuals

# add data.type (from fcal)
oars.response.tree.lfc.map = data.frame(variable = colnames(oars.response.tree.data)) %>% distinct() %>%
  subset(variable != "starch.lfc") %>%
  mutate(datatype = ifelse(variable %in% colnames(oars.mgx.taxa), "MGX",
                           ifelse(variable %in% colnames(oars.mpx.kegg.mat), "Pathway", 
                                  ifelse(variable %in% colnames(oars.mpx.cog.mat), "COG", 
                                         ifelse(variable %in% colnames(oars.mpx.cazy.mat), "CAZy",
                                                ifelse(variable %in% colnames(oars.mbx.raw.mat.filt), "MBX",
                                                       ifelse(variable %in% colnames(oars.asv.data.glom.prep), "ASV", "Fecal calprotectin"))))))) %>% data.frame()
oars.response.tree.lfc.map$feature = oars.response.tree.lfc.map$variable
oars.response.tree.lfc.map$variable = ifelse(nchar(as.character(oars.response.tree.lfc.map$variable))>40, paste(substr(oars.response.tree.lfc.map$variable, 1, 40), "...", sep=""), as.character(oars.response.tree.lfc.map$variable))

# shorten
oars.response.tree.pls.cor.1$feature = ifelse(nchar(as.character(oars.response.tree.pls.cor.1$feature))>40, paste(substr(oars.response.tree.pls.cor.1$feature, 1, 40), "...", sep=""), as.character(oars.response.tree.pls.cor.1$feature))
oars.response.tree.pls.cor.1 = merge(oars.response.tree.pls.cor.1, oars.response.tree.lfc.map,  by="feature")

oars.response.tree.pls.cor.2$feature = ifelse(nchar(as.character(oars.response.tree.pls.cor.2$feature))>40, paste(substr(oars.response.tree.pls.cor.2$feature, 1, 40), "...", sep=""), as.character(oars.response.tree.pls.cor.2$feature))
oars.response.tree.pls.cor.2 = merge(oars.response.tree.pls.cor.2, oars.response.tree.lfc.map,  by="feature")

# plot
oars.response.tree.pls.cor.1.plot = ggplot(oars.response.tree.pls.cor.1 %>%
                                             mutate(datatype = factor(datatype, levels=c(names(omics.colors), "Fecal calprotectin")))%>%
                                             mutate(feature = ifelse(feature == "fcal", "Fecal calprotectin", feature))%>%
                                             #subset(padj < 0.20)%>%
                                             group_by(datatype) %>%
                                             slice_max(order_by=abs(cor), n=5),
                                           aes(x=cor, y=reorder(feature, cor)))+
  geom_point(aes(shape=datatype, fill=datatype), size=3)+
  scale_fill_manual(values=c(omics.colors, "Fecal calprotectin" = "#FFFFFF"))+
  scale_shape_manual(values=c(omics.shapes, "Fecal calprotectin" = 21))+
  facet_wrap(~"Correlated with Scores")+
  theme_classic()+theme(strip.text=element_text(size=10))+
  labs(x="Spearman ρ",
       y="", fill="Data type", shape="Data type")
oars.response.tree.pls.cor.1.plot

oars.response.tree.pls.cor.2.plot = ggplot(oars.response.tree.pls.cor.2 %>%
                                             mutate(datatype = factor(datatype, levels=c(names(omics.colors), "Fecal calprotectin")))%>%
                                             mutate(feature = ifelse(feature == "fcal", "Fecal calprotectin", feature))%>%
                                             #subset(padj < 0.20)%>%
                                             group_by(datatype) %>%
                                             slice_max(order_by=abs(cor), n=5),
                                           aes(x=cor, y=reorder(feature, cor)))+
  geom_point(aes(shape=datatype, fill=datatype), size=3)+
  scale_fill_manual(values=c(omics.colors, "Fecal calprotectin" = "#FFFFFF"))+
  scale_shape_manual(values=c(omics.shapes, "Fecal calprotectin" = 21))+
  facet_wrap(~"Correlated with Ortho-Scores")+
  theme_classic()+theme(strip.text=element_text(size=10))+
  labs(x="Spearman ρ",
       y="", fill="Data type", shape="Data type")
oars.response.tree.pls.cor.2.plot

# shorten feature names
loadings$feature = ifelse(nchar(as.character(loadings$feature))>40, paste(substr(loadings$feature, 1, 40), "...", sep=""), as.character(loadings$feature))
loadings = merge(loadings, oars.response.tree.lfc.map,  by="feature")

oars.response.tree.loadings.plot = ggplot(loadings, aes(x = p1, y = ortho)) +
  geom_point(aes(shape = datatype, fill=datatype), size=2) +
  scale_fill_manual(values=c(omics.colors, "Fecal calprotectin" = "#FFFFFF"))+
  scale_shape_manual(values=c(omics.shapes, "Fecal calprotectin" = 21))+
  
  geom_nodetext_repel(aes(label = ifelse(ahull==0, feature, NA)), 
                      size = 2) +
  #scale_color_manual()+
  theme_classic() +theme(legend.position="none")+
  labs(
    x = "Predictive Score",
    y = "Orthogonal Score"
  )
oars.response.tree.loadings.plot

oars.rt.p1 = (oars.response.tree.plot/
                oars.response.tree.cor.plot)

(oars.rt.p1|
    oars.response.tree.loadings.plot)+
  patchwork::plot_layout(widths=c(1,2))

# probably supplementals
oars.response.tree.pls.cor.1.plot
oars.response.tree.pls.cor.2.plot

# conclusion:
# B. adolescentis is variably associated with Starch
# R. gnavus and E. coli are more consistently negatively associated with Starch

# confirmed:
oars.response.tree.data.pls.cor.plot.1 = ggplot(oars.response.tree.data.pls %>% merge(metadata.oars.stool.asv, by="row.names"),
                                                aes(x=starch.lfc, y=Escherichia_coli))+
  geom_smooth(method="lm", color="black")+
  geom_point(shape=21, aes(fill=rs.col),size=2)+
  scale_fill_manual(values=rs.colors)+
  ggpubr::stat_cor(method="spearman")+
  theme_classic()+theme(legend.position="none",
                        text = element_text(size = 8))+
  labs(x="Starch CAZy Log2FC")

oars.response.tree.data.pls.cor.plot.2 = ggplot(oars.response.tree.data.pls %>% merge(metadata.oars.stool.asv, by="row.names"),
                                                aes(x=starch.lfc, y=g__Ruminococcus_B_s__gnavus))+
  geom_smooth(method="lm", color="black")+
  geom_point(shape=21, aes(fill=rs.col),size=2)+
  scale_fill_manual(values=rs.colors)+
  ggpubr::stat_cor(method="spearman")+
  theme_classic()+theme(legend.position="none",
                        text = element_text(size = 8))+
  labs(x="Starch CAZy Log2FC")

oars.response.tree.data.pls.cor.plot.3 = ggplot(oars.response.tree.data.pls %>% merge(metadata.oars.stool.asv, by="row.names"),
                                                aes(x=starch.lfc, y=`fcal.x`))+
  geom_smooth(method="lm", color="black")+
  geom_point(shape=21, aes(fill=rs.col),size=2)+
  scale_fill_manual(values=rs.colors)+
  ggpubr::stat_cor(method="spearman")+
  theme_classic()+theme(legend.position="none",
                        text = element_text(size = 8))+
  labs(x="Starch CAZy Log2FC",
       y="Fecal calprotectin Log2FC")

oars.response.tree.data.pls.cor.plot.1+
  oars.response.tree.data.pls.cor.plot.2+
  oars.response.tree.data.pls.cor.plot.3


ggplot(oars.response.tree.data.pls %>% merge(metadata.oars.stool.asv, by="row.names"),
       aes(x=starch.lfc, y=Waltera_intestinalis))+
  geom_smooth(method="lm", color="black")+
  geom_point(shape=21, aes(fill=rs.col),size=2)+
  scale_fill_manual(values=rs.colors)+
  ggpubr::stat_cor(method="spearman")+
  theme_classic()+theme(legend.position="none",
                        text = element_text(size = 8))+
  labs(x="Starch CAZy Log2FC")
# not sig
ggplot(oars.response.tree.data.pls %>% merge(metadata.oars.stool.asv, by="row.names"),
       aes(x=starch.lfc, y=GH77))+
  geom_smooth(method="lm", color="black")+
  geom_point(shape=21, aes(fill=rs.col),size=2)+
  scale_fill_manual(values=rs.colors)+
  ggpubr::stat_cor(method="spearman")+
  theme_classic()+theme(legend.position="none",
                        text = element_text(size = 8))+
  labs(x="Starch CAZy Log2FC")

# Heatmap?
# subset to 95% features
# add starch.lfc as an annotation
subset(loadings, ahull==0)$feature

oars.response.tree.pls.cor.1.features =oars.response.tree.pls.cor.1 %>%
  mutate(datatype = factor(datatype, levels=c(names(omics.colors), "Fecal calprotectin")))%>%
  mutate(feature = ifelse(feature == "fcal", "Fecal calprotectin", feature))%>%
  #subset(padj < 0.20)%>%
  group_by(datatype) %>%
  slice_max(order_by=abs(cor), n=5)
oars.response.tree.pls.cor.2.features =oars.response.tree.pls.cor.2 %>%
  mutate(datatype = factor(datatype, levels=c(names(omics.colors), "Fecal calprotectin")))%>%
  mutate(feature = ifelse(feature == "fcal", "Fecal calprotectin", feature))%>%
  #subset(padj < 0.20)%>%
  group_by(datatype) %>%
  slice_max(order_by=abs(cor), n=5)

oars.response.tree.heatmap.data = oars.response.tree.data
colnames(oars.response.tree.heatmap.data) = gsub("fcal", "Fecal calprotectin", colnames(oars.response.tree.heatmap.data))

pheatmap::pheatmap(t(oars.response.tree.heatmap.data[,unique(c(oars.response.tree.pls.cor.2.features$feature))]),
                   color=colorRampPalette(c("blue","white", "red"))(100),
                   clustering_distance_rows = "correlation",
                   clustering_distance_cols = "correlation",
                   breaks=c(seq(min(oars.response.tree.data[,subset(loadings, ahull==0)$feature]), 0, length.out=ceiling(100/2) + 1), 
                            seq(max(oars.response.tree.data[,subset(loadings, ahull==0)$feature])/100, max(oars.response.tree.data[,subset(loadings, ahull==0)$feature]), length.out=floor(100/2))),
                   annotation_col=oars.response.tree.data%>%select(starch.lfc),
                   #annotation_colors = oars.asv.lfc.treatment.prev.mapping.colors,
                   border_color = "white")

# :: // defunct Response Tree 2 (cor) --------------------------------------------------

# Use simple correlations
# More correlated > deeper in tree
# Less correlated > branches

# add starch.lfc
oars.fecalcal.omics.lfc.rt = merge(oars.fecalcal.omics.lfc,
                                   na.omit(oars.response.tree.meta), by="row.names")
dim(oars.fecalcal.omics.lfc.rt)


# spearman correlation (takes a min)
metadata.oars.response.tree = Hmisc::rcorr(as.matrix(oars.fecalcal.omics.lfc.rt[,colnames(oars.fecalcal.omics.lfc.rt)!="Row.names"]), type="spearman")

metadata.oars.response.tree.df = 
  data.frame(reshape2::melt(metadata.oars.response.tree$r)) %>%
  mutate(pval = reshape2::melt(metadata.oars.response.tree$P)$value) %>%
  subset(!is.na(pval)) %>%
  arrange(pval)

metadata.oars.response.tree.bh = subset(metadata.oars.response.tree.df, Var1 == "starch.lfc") %>%
  mutate(feature = as.character(Var2))%>%
  mutate(padj = p.adjust(pval, method="BH")) %>%
  arrange(padj)

oars.fecalcal.response.tree.plot = ggplot(metadata.oars.response.tree.bh,
                                          aes(x=value, y=(padj)))+
  geom_point(shape=21, size=2.5, aes(fill = ifelse(padj < 0.20, "sig", "notsig")))+
  scale_y_continuous(transform=neg_log10_trans,
                     breaks=c(0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, gsub("\\..*", "", gsub("\\(.*", "", gsub("\\/.*", "", gsub("\\,.*", "", feature)))), NA)),
                           size=3)+  
  scale_fill_manual(values=c("black", "white"))+
  geom_hline(yintercept=0.20, linetype=2, alpha=0.5)+
  facet_wrap(~"Starch CAZy Log2FC Correlation")+
  theme_classic()+theme(legend.position="none",
                        strip.text=element_text(size=10))+
  labs(x="Spearman ρ",
       y="FDR")
oars.fecalcal.response.tree.plot


oars.fecalcal.response.tree.plot = ggplot(metadata.oars.response.tree.bh,
                                          aes(x=rank(1-abs(value)), y=(1-abs(value))))+
  geom_point(shape=21, size=2.5, aes(fill = ifelse(padj < 0.20, "sig", "notsig")))+
  #scale_y_continuous(transform=neg_log10_trans,
  #                   breaks=c(0.01, 0.05,0.1, 0.20, 0.5, 1))+
  ggrepel::geom_text_repel(aes(label=ifelse(padj < 0.20, gsub("\\..*", "", gsub("\\(.*", "", gsub("\\/.*", "", gsub("\\,.*", "", feature)))), NA)),
                           size=3)+  
  scale_fill_manual(values=c("black", "white"))+
  geom_hline(yintercept=0.20, linetype=2, alpha=0.5)+
  facet_wrap(~"Starch CAZy Log2FC Correlation")+
  theme_classic()+theme(legend.position="none",
                        strip.text=element_text(size=10))+
  labs(x="Spearman ρ",
       y="FDR")
oars.fecalcal.response.tree.plot


oars.fecalcal.omics.lfc.rt.df = reshape2::melt(oars.fecalcal.omics.lfc.rt[,colnames(oars.fecalcal.omics.lfc.rt)%in%c("starch.lfc", subset(metadata.oars.response.tree.bh, pval < 0.05)$feature)],
                                               id.vars="fcal")

oars.fecalcal.omics.lfc.rt.df$clean.feature = ifelse(nchar(as.character(oars.fecalcal.omics.lfc.rt.df$variable))>30, paste(substr(oars.fecalcal.omics.lfc.rt.df$variable, 1, 30), "...", sep=""), as.character(oars.fecalcal.omics.lfc.rt.df$variable))

# add omics shape

oars.fecalcal.omics.lfc.rt.df = oars.fecalcal.omics.lfc.rt.df %>%
  mutate(data.type = ifelse(variable %in% colnames(oars.mgx.taxa), "MGX",
                            ifelse(variable %in% colnames(oars.mpx.kegg.mat), "Pathway", 
                                   ifelse(variable %in% colnames(oars.mpx.cog.mat), "COG", 
                                          ifelse(variable %in% colnames(oars.mpx.cazy.mat), "CAZy", 
                                                 ifelse(variable %in% colnames(oars.mbx.raw.mat.filt), "MBX", 
                                                        ifelse(variable %in% colnames(oars.asv.data.glom), "ASV", "Fecal calprotectin")))))))
oars.fecalcal.omics.lfc.rt.df$omic.type = ifelse(grepl(paste(c("Pathway", "COG", "CAZy"),collapse="|"), oars.fecalcal.omics.lfc.rt.df$data.type), "MPX", oars.fecalcal.omics.lfc.rt.df$data.type)

oars.fecalcal.omics.lfc.rt.df$data.type = factor(oars.fecalcal.omics.lfc.rt.df$data.type, levels=c("ASV", "MGX", "COG", "Pathway", "CAZy", "MBX", "Fecal calprotectin"))
oars.fecalcal.omics.lfc.rt.df$omic.type = factor(oars.fecalcal.omics.lfc.rt.df$omic.type, levels=c("ASV", "MGX", "MPX", "MBX", "Fecal calprotectin"))

oars.fecalcal.omics.lfc.rt.data = oars.fecalcal.omics.lfc.rt.df %>%
  subset(variable %in% slice_min(metadata.oars.response.tree.bh, order_by=abs(pval), n=8)$feature) %>%
  arrange(order_by=value)  %>%
  subset(!is.na(value)) %>%
  mutate(Var2 = variable)

oars.fecalcal.omics.lfc.rt.plots = ggplot(oars.fecalcal.omics.lfc.rt.data,
                                          aes(x=value,
                                              y=starch.lfc))+
  geom_point(size = 2, aes(shape=data.type, fill=data.type), color="black", alpha=0.9) +
  scale_shape_manual(values=omics.shapes)+
  guides(fill=FALSE, shape=FALSE, size = FALSE)+
  # colors: 
  scale_fill_manual(values=omics.colors)+
  geom_smooth(method="lm", color="black")+
  ggpubr::stat_cor(aes(label = ..r.label..), method="spearman", size=3, 
                   label.y=5.5, #vjust=1.5,
                   label.x.npc=0.5, hjust=0.5)+
  #scale_y_log10()+
  #scale_x_log10()+
  theme_classic()+theme(strip.text=element_text(size=8))+
  facet_wrap(~clean.feature, scales="free", nrow=2)+
  labs(y="Calprotectin Log2FC",
       x="Feature Log2FC")
oars.fecalcal.omics.lfc.rt.plots


# :: 4 ALT mlr // defunct ------------------------------------------------------------------

oars.individualized.lm = do.call(rbind, lapply(colnames(oars.individualized.lm.data)[!colnames(oars.individualized.lm.data)%in%c("standard.name", "oars.days", "HM", "starch", "compliant", "fcal")], function(x){
  print(x)
  data.subset = oars.individualized.lm.data[colnames(oars.individualized.lm.data)%in%c("standard.name", "oars.days", "HM", "starch", "compliant")]
  data.subset$feature = oars.individualized.lm.data[,x]
  # process (scale so that features are on the same scale)
  data.subset$time = scale(data.subset$oars.days)
  # skip if feature is not unique in < 3 samples
  if(sum(data.subset$feature!=min(data.subset$feature))>=3){
    # process (scale so that features are on the same scale)
    data.subset$feature = scale(log2(data.subset$feature))
    # build model (predict starch with feature)
    lm.output = lmerTest::lmer(feature ~ (time) + (time | HM), data.subset)
    #
    data.frame(slope = lme4::ranef(lm.output)$HM[,2],
               HM = rownames(lme4::ranef(lm.output)$HM)) %>%
      mutate(feature = x) %>%
      mutate(coef = coef(summary(lm.output))[2,1]) %>%
      mutate(pval = coef(summary(lm.output))[2,5])
  } else{
    data.frame(slope = NA,
               HM = NA,
               feature = x,
               coef = NA,
               pval = NA)
  }
}))

oars.individualized.lm = subset(oars.individualized.lm,
                                !is.na(HM))
# add omic data type
oars.individualized.lm = oars.individualized.lm %>%
  subset(feature != "fcal") %>%
  mutate(data.type = ifelse(feature %in% colnames(oars.mgx.taxa.filt.4), "Species",
                            ifelse(feature %in% colnames(oars.mpx.kegg.filt.4), "Pathway", 
                                   ifelse(feature %in% colnames(oars.mpx.cog.filt.4), "COG", 
                                          ifelse(feature %in% colnames(oars.mpx.cazy.filt.4), "CAZy", 
                                                 ifelse(feature %in% colnames(oars.mbx.filt.4), "Metabolite", 
                                                        ifelse(feature %in% colnames(oars.asv.data.glom.filt.4), "ASV", "other")))))))%>%
  mutate(data.type = factor(data.type, levels=rev(c("ASV", "Species", "Pathway", "COG", "CAZy", "Metabolite"))))

# :: first, we can reduce the dimensionality
# (PCA on feature slopes)
oars.individualized.lm.mat = reshape2::acast(oars.individualized.lm,
                                             HM ~ feature, value.var="slope")
dim(oars.individualized.lm.mat) # 3716 features
range(oars.individualized.lm.mat)

# remove 0 var
oars.individualized.lm.mat = oars.individualized.lm.mat[,apply(oars.individualized.lm.mat, 2, sd)>0]
dim(oars.individualized.lm.mat) # 3715 features remain (1 removed)

# PCA does not work (0 SD per feature)
# use Spearman correlation on samples
oars.individualized.lm.pca = prcomp(oars.individualized.lm.mat, scale=T, center=T)
oars.individualized.lm.pca.df = oars.individualized.lm.pca$x[,c(1,2)] %>%
  as.data.frame() %>% mutate(HM = rownames(.))
oars.individualized.lm.pca.var = (oars.individualized.lm.pca$sdev^2 / sum(oars.individualized.lm.pca$sdev^2)) * 100
# plot
oars.individualized.lm.pca.plot = ggplot(oars.individualized.lm.pca.df,
                                         aes(x=PC1, y=PC2))+
  geom_point(aes(fill=HM), shape=21, size=3)+
  ggnetwork::geom_nodetext_repel(aes(label=HM))+
  theme_classic()+
  facet_wrap(~"Random Slope PCA")+
  theme(legend.position="none",
        strip.text=element_text(size=10))+
  labs(x=paste0("PC1 (", round(oars.individualized.lm.pca.var[1], digits=1), "%)", sep=""),
       y=paste0("PC2 (", round(oars.individualized.lm.pca.var[2], digits=1), "%)", sep=""))
oars.individualized.lm.pca.plot



# :: second, we can identify features with the greatest variability with slopes
# 1. calculate sd of slopes
oars.individualized.features.slopes = oars.individualized.lm %>%
  #select(feature, HM, slope) %>%
  #subset(feature %in% oars.response.features$feature) %>%
  group_by(feature) %>%
  mutate(slope.ratio = abs(table(sign(slope))[1]/table(sign(slope))[2]))%>%
  mutate(slope.sd = sd(slope),
         slope.min = range(slope)[1],
         slope.max = range(slope)[2],
         slope.range = range(slope)[2] - range(slope)[1],
         slope.mean = mean(slope)) %>%
  select(feature, pval, slope.sd,slope.min, slope.max, slope.range, slope.mean, slope.ratio, coef, data.type) %>% distinct() %>%
  mutate(padj = p.adjust(pval, method="BH"))%>%
  arrange(pval)
# visualize distributions across omics data types
ggplot(oars.individualized.features.slopes %>% subset(feature != "fcal"),
       aes(x=(slope.sd), y=data.type))+
  scale_x_sqrt()+
  ggridges::geom_density_ridges2(aes(fill=data.type),
                                 quantile_lines = T, quantiles = 2)+
  scale_fill_manual(values=omics.colors)+
  theme_classic()+theme(legend.position="none")+
  labs(x="Deviation Score (Random Slope)", y="")
# metabolites have lower variance


# :: third, we can plot the most representative features per individual
# (select the largest magnitude slopes per omic per HM)
oars.individualized.lm.ind = oars.individualized.lm %>%
  group_by(feature) %>%
  mutate(slope.sd = sd(slope)) %>%
  group_by(data.type) %>%
  slice_max(order_by=slope.sd, n=11*10)

table(oars.individualized.lm.ind$data.type)
# rough plot
ggplot(oars.individualized.lm.ind %>% 
         mutate(data.type = factor(data.type, levels=c("ASV", "Species", "Pathway", "COG", "CAZy", "Metabolite"))) %>%
         mutate(feature = factor(feature, levels=unique(oars.individualized.lm.ind$feature))),
       aes(x=(slope), y=feature))+
  geom_vline(xintercept=0, linetype=2)+
  geom_point(aes(shape=data.type, fill=data.type),size=2.5)+
  #ggnetwork::geom_nodetext(aes(label=HM), size=2)+
  scale_fill_manual(values=omics.colors)+
  scale_shape_manual(values=omics.shapes)+
  facet_wrap(~data.type, scales="free_y", ncol=1)+
  theme_classic()+
  labs(x = "Deviation Score (Random Slope)")

# note: participants with only 0M and 3M were removed because they skewed results (over emphasized variation)

# :: now let's plot some of the trends
pheatmap::pheatmap(reshape2::acast(oars.individualized.lm.ind, feature ~ HM, value.var="slope"),
                   scale="column")


# >>> 4 Omic Trajectories -------------------------------------------------

# calculate random slopes for all features
# among HMs with matching omics (and 2+ timepoints)

# now run stats
oars.individualized.relmer.results = do.call(rbind, lapply(colnames(oars.individualized.lm.data)[2:(ncol(oars.individualized.lm.data)-6)], function(x){
  print(as.character(x))
  # x = "GT25"
  data.subset = oars.individualized.lm.data[,colnames(oars.individualized.lm.data) %in% c(x, "HM", "oars.days", "oars.on.rs", "standard.name")]
  colnames(data.subset)[2] = "feature"
  # log2 scale
  data.subset$feature = log2(data.subset$feature)
  # apply prev filter
  if(length(unique(data.subset$feature)) <= 3){
    data.frame(slope = NA,
               HM = NA,
               feature = x)
  } 
  if(length(unique(data.subset$feature)) > 3){
    #x = unique(oars.individualized.lfc$variable)[1]
    # subset to selected feature
    # lmer (RS vs washout)
    lm.results = lmerTest::lmer(scale(feature) ~ oars.on.rs + (oars.on.rs|HM), data.subset)
    data.frame(slope = lme4::ranef(lm.results)$HM[,2],
               HM = rownames(lme4::ranef(lm.results)$HM),
               feature = x)
  }
}))
oars.individualized.relmer.results = oars.individualized.relmer.results  %>%
  mutate(`Data type` = ifelse(feature %in% colnames(oars.mgx.taxa.filt.4), "Species",
                              ifelse(feature %in% colnames(oars.mpx.kegg.filt.4), "Pathway", 
                                     ifelse(feature %in% colnames(oars.mpx.cog.filt.4), "COG", 
                                            ifelse(feature %in% colnames(oars.mpx.cazy.filt.4), "CAZy", 
                                                   ifelse(feature %in% colnames(oars.mbx.filt.4), "Metabolite", 
                                                          ifelse(feature %in% colnames(oars.asv.data.glom.filt.4), "ASV", "Fecal calprotectin"))))))) %>% as.data.frame()

# calculate sd of each feature's random slopes
oars.individualized.relmer.results.features = oars.individualized.relmer.results %>%
  group_by(feature) %>%
  mutate(feature.sd = sd(slope)) %>%
  mutate(feature.max = max(abs((slope)))) %>%
  mutate(feature.delta = max(slope) - min(slope))%>%
  dplyr::select(feature, feature.sd, feature.max, feature.delta, `Data type`) %>% distinct() %>%
  # select top X based on parameter
  arrange(-feature.sd) %>% # feature slope standard deviation
  arrange(-abs(feature.max)) %>%
  arrange(-feature.delta)

ggplot(oars.individualized.relmer.results.features %>%
         mutate(`Data type` = factor(`Data type`, levels=names(omics.colors))),
       aes(x=`Data type`, y=(feature.sd)))+
  geom_violin(aes(fill=`Data type`), 
              draw_quantiles = 0.5)+
  geom_hline(yintercept = median(oars.individualized.relmer.results.features$feature.sd), linetype=2, linewidth=0.2)+
  scale_fill_manual(values=omics.colors)+
  theme_classic()+theme(legend.position="none",
                        axis.text.x = element_text(angle=45, hjust=1))+
  labs(x="", y="Random Slope SD")

# dunns test
oars.individualized.relmer.results.features.lm = lm(feature.sd ~ `Data type`, oars.individualized.relmer.results.features) %>%
  summary() %>% coef() %>% as.data.frame()


oars.individualized.relmer.results.features.chisq = do.call(rbind, lapply(unique(oars.individualized.relmer.results.features$`Data type`), function(x){
  data.subset = oars.individualized.relmer.results.features %>%
    mutate(datatype = ifelse(`Data type` == x, "selected", "other"))
  test.result = data.subset %>%
    mutate(above1 = ifelse(feature.sd > 1, "yes", "no")) %>%
    as.data.frame() %>%
    dplyr::select(`datatype`, above1) %>% 
    table() %>% 
    chisq.test()
  data.frame(`Data type` = x,
             statistic = test.result$statistic,
             pval = test.result$p.value)
}))
cor.test(oars.individualized.relmer.results.features.chisq$statistic,
         oars.individualized.relmer.results.features.chisq$pval)

oars.individualized.relmer.results.features.lm
oars.individualized.relmer.results.features.chisq

# >>> 4 Individualized Trajectories -----------------------------------------------

# With individualized analyzes, we cannot discriminate noise from signal
# But, we can maximize the likelihood of a change = signal
# Approach: is there a change from baseline (with RS) that returns back to baseline (during washout)

# Log2FC from baseline for each timepoint (0M-6M, 0M-9M, 0M-12M)
# Calculate slope of 6M Log2FC to 9M and 12M
# Identify strongest and most unique for each participant
# Log2FC accounts for baseline, and slope accounts for likelihood of noise 
# (i.e. an up/downregulation followed by return to baseline, rather than isolated incident of change)
# If you think of a GAM, we're looking for the biggest triangle shape per feature and individual

rownames(oars.individualized.lm.data)
# subset to HMs with all matching omics for 0-12M samples (only n=5 have all timepoints)
oars.individualized.lfc = subset(oars.individualized.lm.data, HM %in% c("HM0819", "HM0844", "HM0874", "HM0883", "HM0899"))
# melt and calculate Log2FC from baseline
oars.individualized.lfc = reshape2::melt(oars.individualized.lfc, id.vars=c("standard.name", "HM", "compliant", "oars.days", "oars.on.rs")) %>% as.data.frame() %>%
  group_by(HM, variable) %>%
  mutate(lfc = log2(value) - log2(value[oars.on.rs == "preRS"])) %>%
  subset(oars.days > 0)

# now run stats
oars.individualized.lfc.results = do.call(rbind, lapply(unique(oars.individualized.lfc$variable), function(x){
  print(as.character(x))
  #x = "GT25"
  data.subset = subset(oars.individualized.lfc, variable == x)
  # apply prev filter
  if(length(unique(data.subset$lfc)) <= 3){
    data.frame(slope = NA,
               HM = NA,
               feature = x)
  } 
  if(length(unique(data.subset$lfc)) > 3){
    #x = unique(oars.individualized.lfc$variable)[1]
    # subset to selected feature
    # lmer (RS vs washout)
    lm.results = lmerTest::lmer(scale(lfc) ~ oars.on.rs + (oars.on.rs|HM), data.subset)
    data.frame(slope = lme4::ranef(lm.results)$HM[,2],
               HM = rownames(lme4::ranef(lm.results)$HM),
               feature = x)
  }
}))

# calculate sd of each feature's random slopes
oars.individualized.lfc.results.features = oars.individualized.lfc.results %>%
  group_by(feature) %>%
  mutate(feature.sd = sd(slope)) %>%
  mutate(feature.max = max(abs((slope)))) %>%
  mutate(feature.delta = max(slope) - min(slope))%>%
  dplyr::select(feature, feature.sd, feature.max, feature.delta) %>% distinct() %>%
  # select top X based on parameter
  arrange(-feature.sd) %>% # feature slope standard deviation
  arrange(-abs(feature.max)) %>%
  arrange(-feature.delta) %>%
  head(n=20)

cor.test(oars.individualized.lfc.results.features$feature.sd,
         oars.individualized.lfc.results.features$feature.delta) # cor = 0.77

# prepare heatmap maps
# samples
oars.individualized.lfc.hm.map = subset(metadata.oars.stool, standard.name %in% oars.individualized.lfc$standard.name) %>%
  dplyr::select(HM, RS_Name, timing, standard.name)
rownames(oars.individualized.lfc.hm.map) = oars.individualized.lfc.hm.map$standard.name
oars.individualized.lfc.hm.map$standard.name = NULL
colnames(oars.individualized.lfc.hm.map) = c("HM", "RS_Name", "Timing")
# features
oars.individualized.lfc.feature.map = oars.individualized.lfc.results.features %>% dplyr::select(feature, feature.delta) %>%
  mutate(`Data type` = ifelse(feature %in% colnames(oars.mgx.taxa.filt.4), "Species",
                              ifelse(feature %in% colnames(oars.mpx.kegg.filt.4), "Pathway", 
                                     ifelse(feature %in% colnames(oars.mpx.cog.filt.4), "COG", 
                                            ifelse(feature %in% colnames(oars.mpx.cazy.filt.4), "CAZy", 
                                                   ifelse(feature %in% colnames(oars.mbx.filt.4), "Metabolite", 
                                                          ifelse(feature %in% colnames(oars.asv.data.glom.filt.4), "ASV", "Fecal calprotectin"))))))) %>% as.data.frame()
rownames(oars.individualized.lfc.feature.map) = oars.individualized.lfc.feature.map$feature
oars.individualized.lfc.feature.map$feature = NULL
colnames(oars.individualized.lfc.feature.map)[1] = "Slope Delta"

oars.individualized.lfc.heatmap = subset(oars.individualized.lfc, variable %in% oars.individualized.lfc.results.features$feature) %>%
  reshape2::acast(variable ~ standard.name, value.var="lfc")


# :: heatmap --------------------------------------------------------------


## Sept 19, 2025

(oars.spear.dist.asv.dendro.plot|
   oars.spear.dist.mgx.dendro.plot|
   oars.spear.dist.pathway.dendro.plot|
   oars.spear.dist.cog.dendro.plot|
   oars.spear.dist.cazy.dendro.plot|
   oars.spear.dist.mbx.dendro.plot) %>%
  ggsave(filename="./oars_plots/TEMP_oars_dendros.pdf",
         width=18, height=6, device = cairo_pdf)

oars.individualized.lfc.trajectory.plot%>%
  ggsave(filename="./oars_plots/TEMP_oars_ind_trajectory.pdf",
         width=12, height=10, device = cairo_pdf)


pheatmap::pheatmap(oars.individualized.lfc.heatmap,
                   color=colorRampPalette(c("blue","white", "red"))(100),
                   clustering_distance_rows = "correlation",
                   #clustering_distance_cols = "correlation",
                   cluster_cols = F,
                   fontsize_row = 8,
                   fontsize_col = 8,
                   annotation_col = oars.individualized.lfc.hm.map %>% dplyr::select(-Timing),
                   annotation_row = oars.individualized.lfc.feature.map,
                   annotation_colors = list(Phase = c(`onRS` = "black", `postRS` = "white"),
                                            HM = c(`HM0819` = viridis::inferno(5)[1],
                                                   `HM0844` = viridis::inferno(5)[2],
                                                   `HM0874` = viridis::inferno(5)[3],
                                                   `HM0883` = viridis::inferno(5)[4],
                                                   `HM0899` = viridis::inferno(5)[5]),
                                            #`Timing` = colorRampPalette(c("white","darkgreen"))(4),
                                            `RS_Name` = rs.colors,
                                            `Data type` = c(ASV = "#8DD3C7",
                                                            Species = "#FFFFB3",
                                                            Pathway = "#BEBADA",
                                                            COG = "#FB8072",
                                                            CAZy = "#80B1D3",
                                                            Metabolite = "#FDB462")),
                   cutree_rows = 5,
                   show_colnames=F,
                   breaks=c(seq(min(na.omit(oars.individualized.lfc.heatmap)), 0, length.out=ceiling(100/2) + 1), 
                            seq(max(na.omit(oars.individualized.lfc.heatmap))/100, max(na.omit(oars.individualized.lfc.heatmap)), length.out=floor(100/2))))%>%
  ggsave(filename="./oars_plots/TEMP_oars_ind_heatmap.pdf",
         width=9, height=7, device = cairo_pdf)

# among participants with complete datapoints (every time point for every data set)
# identify the most extreme differences in slopes (trajectories) of log2fc values relative to baseline samples
# Top 50 features; shows unique, individualized signatures

# e.g. HM0819: Signal transduction histidine kinase --> falls, rises
# e.g. HM0819: Ruminococcus_E --> rises, falls
# e.g. HM0844: ABC-type Zn uptake and Cu/Ag efflux pump--> rises, falls
# e.g. HM0844: Klebsiella granulomatis --> rises, falls
# e.g. HM0874: GGB3534 --> falls, rises
# e.g. HM0874: Klebsiella granulomatis --> rises, falls
# e.g. HM0883: Thiamine transporter --> falls, rises
# e.g. HM0899: Thiamine transporter --> rises, returns
# e.g. HM0899: Heme binding protein --> falls, rises

# individualized responses across metal uptake and pathobionts

oars.individualized.lfc.trajectory = subset(oars.individualized.lfc,  variable %in% rownames(oars.individualized.lfc.heatmap) &
                                              grepl(paste(c("Signal transduction", "Ruminococcus_E", "ABC-type Zn", "granulomatis", "GGB3534", "Thiamine transporter", "heme binding"), collapse="|"), variable))
# shorten names
oars.individualized.lfc.trajectory$clean.feature = ifelse(nchar(as.character(oars.individualized.lfc.trajectory$variable))>30, paste(substr(oars.individualized.lfc.trajectory$variable, 1, 30), "...", sep=""), as.character(oars.individualized.lfc.trajectory$variable))

# :: trajectories --------------------------------------------------------------

oars.individualized.lfc.trajectory.plot = ggplot(oars.individualized.lfc.trajectory,
                                                 aes(x=as.numeric(as.factor(oars.on.rs)), y=lfc))+
  geom_hline(yintercept=0, linetype=2, linewidth=0.2)+
  geom_segment(data=oars.individualized.lfc.trajectory %>% subset(oars.on.rs == "onRS") %>%
                 dplyr::group_by(HM, variable) %>% mutate(mean.on.rs = mean(lfc)), aes(x=1, xend=2, y=0, yend=(mean.on.rs)))+
  geom_smooth(method="lm", color="black")+
  geom_point(shape=21, aes(fill=oars.on.rs), size=3)+
  scale_fill_manual(values=c("onRS" = "red", "postRS" = "white"))+
  scale_x_continuous(limits=c(1,3))+
  theme_classic()+theme(axis.ticks.x=element_blank(),
                        axis.text.x=element_blank(),
                        legend.position="none")+
  facet_grid(HM~clean.feature, scales="free")+
  labs(x="On RS versus Washout", y="Log2FC to Baseline")



## Sept 19, 2025

oars.spear.dist.asv.dendro.plot|
  oars.spear.dist.mgx.dendro.plot|
  oars.spear.dist.pathway.dendro.plot|
  oars.spear.dist.cog.dendro.plot|
  oars.spear.dist.cazy.dendro.plot|
  oars.spear.dist.mbx.dendro.plot

oars.individualized.lfc.trajectory.plot



# >>> // 6.B MACHINE LEARNING -----------------------------------------------------------

# goal: see if features predict fecal cal response
# note: it can't (AUC < 0.65)

# :: process ---------------------------------------------------------

# :: ML Omics -------------------------------------------------

# Goal: compare AUC of omics
data.types = c("ASV", "Species", "Pathway", "COG", "CAZy","Metabolite", "FFQ", "All")
iters = c(1:15)

oars.loocv.models.ml.omics.fcal = rf.function(data.types = c("ASV", "Species", "Pathway", "COG", "CAZy","Metabolite", "FFQ", "All"),
                                              iters = 15,
                                              target = "fcal",
                                              output = "AUC")

# summarize
oars.loocv.models.ml.omics.fcal.df = oars.loocv.models.ml.omics.fcal %>%
  group_by(data.type) %>%
  dplyr::select(data.type, auc) %>% distinct() %>%
  mutate(mean.auc = mean(auc),
         median.auc = median(auc),
         auc.low = mean(auc) - (sd(auc)/sqrt(n()) * 1.96), # 95% CI
         auc.high = mean(auc) + (sd(auc)/sqrt(n()) * 1.96)) %>%
  dplyr::select(mean.auc, median.auc, auc.low, auc.high, data.type) %>% distinct() %>%
  arrange(-mean.auc)
oars.loocv.models.ml.omics.fcal.df
# Nice

oars.loocv.models.ml.omics.fcal.rf.plot = ggplot(oars.loocv.models.ml.omics.fcal.df %>%
                                                   mutate(data.type = ifelse(data.type == "All", "Multi-Omic", data.type))%>%
                                                   mutate(data.type = factor(data.type, levels= c("FFQ", "ASV", "Species", "COG","Pathway", "CAZy","Metabolite", "Multi-Omic"))),
                                                 aes(y=reorder(data.type, mean.auc), x=mean.auc))+
  geom_bar(stat="identity", width=0.75,
           aes(fill=data.type),color="black")+
  geom_segment(aes(x=auc.low, xend=auc.high))+
  scale_x_continuous(breaks=seq(0,1, by=0.1),
                     limits=c(0,1))+
  geom_text(data=oars.loocv.models.ml.omics.fcal.df %>%
              mutate(data.type = ifelse(data.type == "All", "Multi-Omic", data.type))%>%
              mutate(data.type = factor(data.type, levels= c("FFQ", "ASV", "Species", "COG","Pathway", "CAZy","Metabolite", "Multi-Omic"))), 
            aes(x=auc.high, y=data.type, 
                label=paste(" ", round(mean.auc, digits=3))), hjust=0, size=3)+
  scale_fill_manual(values=c(omics.colors, "Multi-Omic" = "black", "FFQ" = "white"))+
  #geom_segment()+
  facet_wrap(~"RandomForest")+
  theme_classic()+theme(legend.position="none",
                        strip.text=element_text(size=10))+
  labs(x="Mean AUC", y="")
oars.loocv.models.ml.omics.fcal.rf.plot

# :: RF LOOCV ----------------------------------------------------------------

# Now that we know COG has highest AUC, let's plot it's ROC curve

oars.loocv.models.rf.fcal = rf.function(data.types = "COG",
                                        iters = 15,
                                        target = "fcal",
                                        output = "ROC")

# calculate values and extract best
oars.loocv.models.rf.fcal.df = oars.loocv.models.rf.fcal %>%
  group_by(iter) %>%
  summarize(auc = pROC::auc(high, pred, 
                            levels=c("high", "low"),  # define case = "high" delta pH
                            direction="<")[1])%>%  
  mutate(mean.auc = mean(auc),
         median.auc = median(auc),
         auc.low = mean(auc) - (sd(auc)/sqrt(n()) * 1.96), # 95% CI
         auc.high = mean(auc) + (sd(auc)/sqrt(n()) * 1.96)) %>%
  dplyr::select(mean.auc, median.auc, auc.low, auc.high) %>% distinct()

# extract best
oars.loocv.models.fcal.best.stats = subset(oars.loocv.models.rf.fcal.df, mean.auc == max(oars.loocv.models.rf.fcal.df$mean.auc))

# plot
oars.loocv.fcal.roc = do.call(rbind, lapply(1:15, function(seed){
  data.subset = subset(oars.loocv.models.rf.fcal, iter == seed)
  # calibrate pred; not necessary
  #data.subset$pred = 1 / (1 + exp(-data.subset$pred))
  data.frame(sens = pROC::roc(data.subset$high, data.subset$pred)$sensitivities,
             spec = pROC::roc(data.subset$high, data.subset$pred)$specificities,
             iter = seed)
}))


# Step 2: Compute mean and standard error for sensitivities across iterations
roc_summary_fcal <- oars.loocv.fcal.roc %>%
  group_by(spec) %>%  # Group by specificity (or alternatively by sens)
  summarise(
    mean_sens = mean(sens, na.rm = TRUE),
    se_sens = sd(sens, na.rm = TRUE) / sqrt(n()),  # Standard error
    lower = mean_sens - 1.96 * se_sens,           # 95% CI lower bound
    upper = mean_sens + 1.96 * se_sens            # 95% CI upper bound
  ) %>%
  mutate(fpr = 1 - spec) %>%  # False Positive Rate (1 - specificity)
  filter(!is.na(mean_sens) & !is.na(se_sens))  # Remove any NA values


oars.loocv.fcal.roc.plot = ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  # Individual ROC curves for each iteration
  #geom_line(data = oars.rf.loocv.roc %>% arrange(sens),
  #          aes(x=1-spec, y=sens, group=iter))+
  # Error ribbon (95% CI)
  geom_ribbon(data = roc_summary_fcal, 
              aes(x = fpr, ymin = lower, ymax = upper), 
              #fill = RColorBrewer::brewer.pal(n=5, "Set3")[1], alpha = 0.2) +
              fill = "black", alpha=0.2)+
  # Mean ROC curve
  geom_path(data = roc_summary_fcal, 
            aes(x = fpr, y = mean_sens), 
            #color = RColorBrewer::brewer.pal(n=5, "Set3")[1], size = 1) +
            color="black")+
  # add label
  annotate(geom="text", x=0.75, y=0.25,
           label=paste("15x LOOCV\n",
                       "AUC: ", round(oars.loocv.models.fcal.best.stats$mean.auc, digits=2),
                       "\n(", round(oars.loocv.models.fcal.best.stats$auc.low, digits=2), 
                       ", ", round(oars.loocv.models.fcal.best.stats$auc.high,digits=2), ")", sep=""),
           size=4)+
  theme_classic()+theme(strip.text=element_text(size=10))+
  facet_wrap(~"COG RandomForest")+
  xlim(0,1)+
  ylim(0,1)+
  labs(x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)") 
oars.loocv.fcal.roc.plot



# :: RF Importances -------------------------------------------------------


oars.loocv.models.importances.fcal = rf.function(data.types = "COG",
                                                 iters = 15,
                                                 target = "fcal",
                                                 output = "importances")

oars.rf.loocv.imp.fcal.df = oars.loocv.models.importances.fcal %>%
  #subset(imp != 0) %>%
  group_by(feature) %>%
  mutate(n.present = n()) %>%
  mutate(mean.imp = mean(na.omit(imp))) %>%
  mutate(imp.low = mean(na.omit(imp)) - (sd(na.omit(imp))/sqrt(n()) * 1.96)) %>% # 95% CI
  mutate(imp.high = mean(na.omit(imp))+ (sd(na.omit(imp))/sqrt(n()) * 1.96)) %>% # 95% CI
  subset(mean.imp != 0) %>%
  dplyr::select(feature, mean.imp, imp.low, imp.high, n.present) %>% distinct() %>%
  arrange(-mean.imp) %>% data.frame()

# Note: need to be careful selecting which features to highlight,
# as their frequency of being selected impacts results considerably
ggplot(oars.rf.loocv.imp.fcal.df,
       aes(x=n.present, y=mean.imp))+
  scale_x_log10()+geom_smooth()+
  geom_point()+theme_classic()


# log2fc
oars.rf.loocv.imp.fcal.df = oars.rf.loocv.imp.fcal.df %>% 
  subset(!is.na(mean.imp))%>%
  # select the top 10 mean importance
  arrange(mean.imp) %>% slice_max(mean.imp, n=10) 

# calculate wilcox test pvalue of these features between response groups
oars.rf.loocv.imp.wilcox.fcal = do.call(rbind, lapply(oars.rf.loocv.imp.fcal.df$feature, function(x){
  data.subset = oars.all.omics.baseline[,colnames(oars.all.omics.baseline) %in% c(x, "response")]
  colnames(data.subset)[2] = "variable"
  # run wilcox
  wilcox.p = wilcox.test(subset(data.subset, response == "high")[,2],
                         subset(data.subset, response == "low")[,2])
  
  # calculate lfc
  if(x %in% c("stool_water_perc", "shannon", "fd")){
    coef = mean(subset(data.subset, response == "high")[,2]) / mean(subset(data.subset, response == "low")[,2])
  }else{
    coef = mean(log2(subset(data.subset, response == "high")[,2])) - mean(log2(subset(data.subset, response == "low")[,2]))
  }
  data.frame(feature = x,
             pval = wilcox.p$p.value,
             coef = coef)
}))
oars.rf.loocv.imp.fcal.df = merge(oars.rf.loocv.imp.wilcox.fcal,
                                  oars.rf.loocv.imp.fcal.df, by="feature")
oars.rf.loocv.imp.fcal.df$sig = ifelse(oars.rf.loocv.imp.fcal.df$pval < 0.05, "*", "")
oars.rf.loocv.imp.fcal.df = oars.rf.loocv.imp.fcal.df %>% arrange(mean.imp)

# clean up COGs
oars.rf.loocv.imp.fcal.df$clean.feature = gsub(paste(c(" and related.*", "or related.*", ", .*"), collapse="|"),"", oars.rf.loocv.imp.fcal.df$feature)

oars.rf.loocv.imp.fcal.plot = ggplot(oars.rf.loocv.imp.fcal.df,
                                     aes(x=mean.imp, y=reorder(clean.feature, mean.imp)))+
  geom_segment(aes(x=imp.low, xend=imp.high, y=clean.feature, yend=clean.feature), color="black")+
  geom_point(shape=21, aes(fill=(coef)), size=3.5)+
  scale_fill_gradient2(low="blue", high="red")+
  scale_color_manual(values=c("white", "black"))+
  theme_classic()+theme(legend.position="right",
                        axis.title.y = element_blank(),
                        #axis.text.y = element_blank(),
                        # axis.ticks.y=element_blank(),
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=10),
                        strip.background = element_rect(
                          color="black"))+
  guides(color=FALSE)+
  labs(x="Mean Decrease in Accuracy", fill="Log2FC")+
  facet_wrap(~"Feature Importance")
oars.rf.loocv.imp.fcal.plot


# :: RF ttest plots -------------------------------------------------------


# loop through top features and plot; extract microbiome values
oars.rf.loocv.imp.ttest.fcal = do.call(rbind, lapply(slice_max(oars.rf.loocv.imp.fcal.df, mean.imp, n=10)$feature, function(x){
  data.subset = oars.all.omics.baseline[,c(x, "response")]
  rownames(data.subset) = paste(oars.all.omics.baseline$HM, oars.all.omics.baseline$phase, sep="_")
  data.subset$feature = colnames(data.subset)[1]
  data.subset$sample = rownames(data.subset)
  colnames(data.subset)[1] = "value"
  data.subset
}))
# add stats
oars.rf.loocv.imp.ttest.fcal = merge(oars.rf.loocv.imp.ttest.fcal,
                                     oars.rf.loocv.imp.fcal.df, by="feature")
oars.rf.loocv.imp.ttest.fcal.p = oars.rf.loocv.imp.ttest.fcal[,c("feature", "pval")] %>% distinct()

# clean taxa names (so they fit)
oars.rf.loocv.imp.ttest.fcal$clean.feature = gsub(paste(c(" and related.*", "or related.*", ", .*"), collapse="|"),"", oars.rf.loocv.imp.ttest.fcal$feature)
oars.rf.loocv.imp.ttest.fcal.p$clean.feature = gsub(paste(c(" and related.*", "or related.*", ", .*"), collapse="|"),"", oars.rf.loocv.imp.ttest.fcal.p$feature)

# now clip names to 30 char
oars.rf.loocv.imp.ttest.fcal$short.feature = ifelse(nchar(oars.rf.loocv.imp.ttest.fcal$clean.feature)>30, paste(substr(oars.rf.loocv.imp.ttest.fcal$clean.feature, start = 1, stop = 30), "...", sep=""), oars.rf.loocv.imp.ttest.fcal$clean.feature)
oars.rf.loocv.imp.ttest.fcal.p$short.feature = ifelse(nchar(oars.rf.loocv.imp.ttest.fcal.p$clean.feature)>30, paste(substr(oars.rf.loocv.imp.ttest.fcal.p$clean.feature, start = 1, stop = 30), "...", sep=""), oars.rf.loocv.imp.ttest.fcal.p$clean.feature)
oars.rf.loocv.imp.fcal.df$short.feature = ifelse(nchar(oars.rf.loocv.imp.fcal.df$clean.feature)>30, paste(substr(oars.rf.loocv.imp.fcal.df$clean.feature, start = 1, stop = 30), "...", sep=""), oars.rf.loocv.imp.fcal.df$clean.feature)

oars.rf.loocv.imp.ttest.fcal.plot = ggplot(oars.rf.loocv.imp.ttest.fcal %>%
                                             # clean names
                                             mutate(Response = ifelse(response == "low", "Weak\nResponse", "Strong\nResponse")) %>%
                                             # add indicator of not present
                                             group_by(short.feature) %>%
                                             mutate(pseudo = ifelse(value == min(value), "pseudo", "real")) %>%
                                             # reorder taxa based on coef
                                             mutate(short.feature = factor(short.feature, levels=arrange(distinct(oars.rf.loocv.imp.fcal.df[,c("short.feature", "coef")]),-coef)$short.feature)),
                                           aes(x=Response, y=value))+
  geom_boxplot(width=0.3, outlier.shape=NA)+
  #scale_y_log10(labels = scales::label_number(accuracy = 0.01))+
  scale_y_log10()+
  ggbeeswarm::geom_beeswarm(shape=21, aes(fill=response, alpha=pseudo), size=2)+
  scale_alpha_manual(values=c(0.2, 1))+
  theme_classic()+theme(legend.position="none",
                        axis.title.x = element_blank(),
                        plot.title = element_text(hjust = 0.5, size=12),
                        strip.text = element_text(size=8),
                        strip.background = element_rect(color="black"))+
  geom_text(data=oars.rf.loocv.imp.ttest.fcal.p%>%
              # clean names
              mutate(short.feature = factor(short.feature, levels=arrange(distinct(oars.rf.loocv.imp.fcal.df[,c("short.feature", "mean.imp")]),-mean.imp)$short.feature)),
            x=1.5, y=Inf, vjust=1.2, 
            aes(label = paste("p =", round(pval, digits=3))),
            size=2.5)+
  labs(x="", y="Feature Abundance")+
  facet_wrap(~short.feature, ncol=2 ,scales="free_y")
oars.rf.loocv.imp.ttest.fcal.plot

# These results suggest that fecal cal response cannot be predicted
# by the microbiome
# It also lends more support to predicting strong response being "real"

# >> Extra ----------------------------------------------------------------

# create mapping file with responses as 2 columns
oars.data.clustering.meta = metadata.oars.stool.double %>% 
  #subset(compliant==T)%>%
  #subset(reltiming == "pre") %>% 
  dplyr::select(standard.name, response) %>% distinct() %>%
  arrange(standard.name)
rownames(oars.data.clustering.meta) = oars.data.clustering.meta$standard.name

# :: Spearman Cluster ASV -----------------------------------------------------

Samples_Used_in_Metabo = read.csv("~/Documents/PhD/For others/Samples_Used_in_Metabo.csv")

# convert to CLR
oars.asv.data.glom.clr = (oars.asv.data.glom)[rownames(oars.asv.data.glom) %in% Samples_Used_in_Metabo[,1],]
#oars.asv.data.glom.clr = (oars.asv.data.glom.clr)[rownames(oars.asv.data.glom.clr) %in% oars.data.clustering.meta$standard.name,]

oars.asv.data.glom.clr.pa = oars.asv.data.glom.clr
oars.asv.data.glom.clr.pa[oars.asv.data.glom.clr.pa > 0] = 1
oars.asv.data.glom.clr = compositions::clr(oars.asv.data.glom.clr) %>% as.data.frame()
oars.asv.data.glom.clr[oars.asv.data.glom.clr.pa == 0] = min(oars.asv.data.glom.clr[oars.asv.data.glom.clr!=0])-1

oars.asv.data.glom.clr.dist = as.dist(1-Hmisc::rcorr(t(oars.asv.data.glom.clr), type="spearman")$r)
oars.asv.data.glom.clr.dist.tree <- ape::as.phylo(hclust(oars.asv.data.glom.clr.dist, method = "complete"))
oars.asv.data.glom.clr.dist.tree$tip.label <- rownames(as.matrix(oars.asv.data.glom.clr.dist))
library("ggtree")
oars.asv.data.glom.clr.dist.tree.plot <- ggtree(oars.asv.data.glom.clr.dist.tree) %<+% oars.data.clustering.meta +
  geom_tiplab(size = 3, nudge_x=0.01) +
  #coord_flip()+
  geom_tippoint(aes(fill = response), shape=21, size = 3, stroke = 0.3) + 
  theme_tree()+
  xlim(0, max(tree$edge.length) *2.4)+  # widen x-axis beyond tree height
  facet_wrap(~"OARS ASV")+
  theme(strip.text = element_text(size=10),
        strip.background = element_rect(
          color="black", fill="white"))
oars.asv.data.glom.clr.dist.tree.plot


# :: Spearman Cluster MGX -----------------------------------------------------

# convert to CLR
oars.mgx.data.clr = (oars.mgx.taxa)[rownames(oars.mgx.taxa) %in% Samples_Used_in_Metabo[,1],]
#oars.mgx.data.clr = (oars.mgx.data.clr)[rownames(oars.mgx.data.clr) %in% oars.data.clustering.meta$standard.name,]
oars.mgx.data.clr.pa = oars.mgx.data.clr
oars.mgx.data.clr.pa[oars.mgx.data.clr.pa > 0] = 1
oars.mgx.data.clr = compositions::clr(oars.mgx.data.clr) %>% as.data.frame()
oars.mgx.data.clr[oars.mgx.data.clr.pa == 0] = min(oars.mgx.data.clr[oars.mgx.data.clr!=0])-1

oars.mgx.data.clr.dist = as.dist(1-Hmisc::rcorr(t(oars.mgx.data.clr), type="spearman")$r)
oars.mgx.data.clr.dist.tree <- ape::as.phylo(hclust(oars.mgx.data.clr.dist, method = "complete"))
oars.mgx.data.clr.dist.tree$tip.label <- rownames(as.matrix(oars.mgx.data.clr.dist))
library("ggtree")
oars.mgx.data.clr.dist.tree.plot <- ggtree(oars.mgx.data.clr.dist.tree) %<+% oars.data.clustering.meta +
  geom_tiplab(size = 3, nudge_x=0.01) +
  #coord_flip()+
  geom_tippoint(aes(fill = response), shape=21, size = 3, stroke = 0.3) + 
  theme_tree()+
  xlim(0, max(tree$edge.length) *3.2)+  # widen x-axis beyond tree height
  facet_wrap(~"OARS MGX")+
  theme(strip.text = element_text(size=10),
        strip.background = element_rect(
          color="black", fill="white"))
oars.mgx.data.clr.dist.tree.plot


# :: Spearman Cluster MPX -----------------------------------------------------

#oars.mpx.oars.only.protein.mat = readRDS("./metaproteomics/2025_06_28_oars_mpx_protein.Rds")
oars.mpx.oars.only.protein.mat = readRDS("./metaproteomics/2025_06_28_oars_mpx_cog.Rds")
#oars.mpx.oars.only.protein.mat = readRDS("./metaproteomics/2025_06_28_oars_mpx_cazy.Rds")

# convert to log
oars.mpx.data.log = (oars.mpx.oars.only.protein.mat)[rownames(oars.mpx.oars.only.protein.mat) %in% Samples_Used_in_Metabo[,1],]
oars.mpx.data.log = (oars.mpx.data.log)[rownames(oars.mpx.data.log) %in% oars.data.clustering.meta$standard.name,]
oars.mpx.data.log[is.na(oars.mpx.data.log)] = 0
oars.mpx.data.log.pseudo = min(oars.mpx.data.log[oars.mpx.data.log!=0])/2
oars.mpx.data.log = log2(oars.mpx.data.log+oars.mpx.data.log.pseudo)
oars.mpx.data.log[oars.mpx.data.log == 0] = min(oars.mpx.data.log[oars.mpx.data.log!=0])-1

oars.mpx.data.log.dist = as.dist(1-Hmisc::rcorr(t(oars.mpx.data.log), type="spearman")$r)
oars.mpx.data.log.dist.tree <- ape::as.phylo(hclust(oars.mpx.data.log.dist, method = "complete"))
oars.mpx.data.log.dist.tree$tip.label <- rownames(as.matrix(oars.mpx.data.log.dist))
library("ggtree")
oars.mpx.data.log.dist.tree.plot <- ggtree(oars.mpx.data.log.dist.tree) %<+% oars.data.clustering.meta +
  geom_tiplab(size = 3, nudge_x=0.005) +
  #coord_flip()+
  geom_tippoint(aes(fill = response), shape=21, size = 3, stroke = 0.3) + 
  theme_tree()+
  xlim(0, max(tree$edge.length) *0.9)+  # widen x-axis beyond tree height
  facet_wrap(~"OARS MPX")+
  theme(strip.text = element_text(size=10),
        strip.background = element_rect(
          color="black", fill="white"))
oars.mpx.data.log.dist.tree.plot


(oars.asv.data.glom.clr.dist.tree.plot+
    oars.mgx.data.clr.dist.tree.plot+
    oars.mpx.data.log.dist.tree.plot+
    patchwork::plot_layout(guides = "collect"))


(lsarp.cd.asv.data.glom.clr.plot+
    lsarp.cd.asv.mgx.clr.plot+
    lsarp.cd.asv.mpx.log.plot+
    patchwork::plot_layout(guides = "collect"))

