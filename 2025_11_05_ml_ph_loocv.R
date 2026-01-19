### 2025_11_05  ML Redux 2 for Compute Canada
# 2025_11_05_ml_ph_loocv.R

library("ggplot2"); library("ggridges"); library("dplyr"); library("tidyverse")

# upload code
# scp ~/Documents/PhD/git_ml_archfolder/ml_git_main/2025_11_05_ml_ph_loocv.R 'pdobrano@nibi.alliancecan.ca:/home/pdobrano/scratch/rs_ml/ml_git_main/'


# upload files
# scp ~/Documents/PhD/git_ml_archfolder/ml_git_data/2025_11_07_meta_rs_features.Rds 'pdobrano@nibi.alliancecan.ca:/home/pdobrano/scratch/rs_ml/ml_git_data/'

# scp ~/Documents/PhD/git_ml_archfolder/ml_git_data/2025_11_07_ml_ph_target.Rds 'pdobrano@nibi.alliancecan.ca:/home/pdobrano/scratch/rs_ml/ml_git_data/'

# scp ~/Documents/PhD/git_ml_archfolder/ml_git_data/2025_11_07_n117_s25_asv_input_data.rds 'pdobrano@nibi.alliancecan.ca:/home/pdobrano/scratch/rs_ml/ml_git_data/'



home.dir = "."
print(getwd())
# load data
ml_rs_features = readRDS("./ml_git_data/2025_11_07_meta_rs_features.Rds")
ph.data = readRDS("./ml_git_data/2025_11_07_ml_ph_target.Rds")
ml.input.data = readRDS("./ml_git_data/2025_11_07_n117_s25_asv_input_data.rds")

# establish variables
iter.length = 15
iters = 1:iter.length
rs.names.pbs <- c("PBS", "Authentic", "BobsRedMill", "MSPrebiotic", "LetsDoOrganic", "HiMaize260", "Novelose330", "ActistarRT", "FibersymRW", "Versafibe1490")
hms.ml <- rownames(ml.input.data)
hms.length = length(hms.ml) # n = 117


# reduce to algorithms capable of both classification + regression:
algos = c("glmnet", "knn", "ranger", "gbm", "rpart", "nnet",        # FAST algos
          "pls", "svmLinear", "svmRadial", "xgbTree", "stack") # SLOW algos
# note: using n=16 cores + 4G memory

# FAST: run < 3 hours: 16 core 4gb ram: 1 - 90
# glmnet 1-15, knn 16-30, ranger 31-45, gbm 46-60, rpart 61-75, nnet 76-90


# SLOW: need 12 hrs and fewer cores: 4 core 8gb ram, : 91 - 165
# pls 91-105, svmL 106-120, svmR 121-135, xgbTree 136-150, stack 150-165

# prepare grid of conditions 
array.options = data.frame(
  array_no = seq(1:((length(iters)*length(algos)))),
  iter = rep(iters, times=length(algos)),
  algos = rep(algos, each=length(iters))
)
# table(array.options$type) # 165

# select parallel condition
args <- commandArgs(trailingOnly = TRUE)
print(paste0("args = ", args))
array_no <- as.numeric(args)[as.numeric(args) %in% c(1:nrow(array.options))]
print(paste0("array_no = ", array_no))
algo = array.options$algos[array_no]
iter = array.options$iter[array_no]
type = array.options$type[array_no]


print(paste("iter ", iter, ". algo = ", algo, ".", sep=""))


t1 <- Sys.time()

if(algo %in% c("rpart", "pls", "svmLinear", "svmRadial", 
               "xgbTree", "stack")){ core.n = 1 } else { core.n = 48 }
# Note: these algos don't parallelize well


# log transform (+ min/2) predictors here (though, not necessary for tree-based methods)
ml.input.data = ml.input.data + min(ml.input.data[ml.input.data!=0])/2
ml.input.data = log10(ml.input.data)

hms.ml = rownames(ml.input.data)

conditions = "delta_pH"

feature_select = T

#n.features = 15

# Run code chunk for one algo + one iter 
ml.loocv.regression.results  <-
  do.call(rbind, lapply(conditions, function(z){ # conditions
    # subset to one target
    # z = "delta_butyrogen"
    target.subset = ph.data %>% mutate(value = med.ph)
    # subset to RS
    do.call(rbind, lapply(rs.names.pbs, function(x){ # rs.names
      do.call(rbind, parallel::mclapply(hms.ml, function(hm){ # parallel::mc # hms.ml
        # subset to LOOCV patient and RS
        # x=rs.names[1]
        # iter = 1
        # hm = hms.ml[1]
        
        # fix names
        ml.input.data = data.frame(ml.input.data)
        
        # apply feature selection
        # conditionally apply feature selection
        if (feature_select == TRUE & x != "PBS"){
          ml.input.data = ml.input.data[,colnames(ml.input.data) %in% ml_rs_features]
        }
        
        train.data = ml.input.data[rownames(ml.input.data) != hm,] %>% data.frame()
        train.response = subset(target.subset, HM != hm & RS_Name == x)
        test.data = ml.input.data[rownames(ml.input.data) == hm,] %>% data.frame() 
        test.response = subset(target.subset, HM == hm & RS_Name == x)
        # replace with scaled value
        test.response$value = subset(target.subset, RS_Name == x) %>% 
          #mutate(value = scale(value)) %>% 
          subset(HM == hm) %>% dplyr::select(med.ph) %>% as.numeric()
        # create shadow data
        shadow.data <- do.call(cbind, lapply(1:ncol(train.data), function (w) {
          shadow.data.column <- train.data[,w]
          shadow.data.column <- data.frame(shadow.data.column[sample(1:length(shadow.data.column), size=1)])
          colnames(shadow.data.column) <- colnames(train.data)[w]
          shadow.data.column
        }))
        # good
        
        
        # prepare data (merge, etc)
        rownames(train.response) <- train.response$HM
        train.data$HM <- rownames(train.data)
        full.train.data = merge(train.data, train.response[,c("HM", "value")], by="HM")
        rownames(full.train.data) <- full.train.data$HM
        full.train.data$Row.names <- NULL
        full.train.data$HM <- NULL
        
        set.seed(iter)
        if(!algo %in% c("stack", "nnet")){ # for all base-models except nnet
          ml.model <-  caret::train((value) ~.,
                                    data = full.train.data, 
                                    method = algo, 
                                    metric = "RMSE",
                                    preProc = c("center", "scale", "zv"),
                                    trControl = caret::trainControl(method="cv", number=5)
          ) %>% suppressWarnings() %>% suppressMessages()
          t2 <- Sys.time()
          
          # apply base model to test
          ml.test.results = data.frame(predicted_response = predict(ml.model, newdata = test.data))
          ml.shadow.results = data.frame(shadow_response = predict(ml.model, newdata = shadow.data))
        } # end base model
        
        if(algo == "nnet"){ # for nnet
          ml.model <-  caret::train((value) ~.,
                                    data = full.train.data, 
                                    method = algo, 
                                    metric = "RMSE",
                                    preProc = c("center", "scale", "zv"),
                                    linout = T, # added uniquely for nnet
                                    trControl = caret::trainControl(method="cv", number=5)
          ) %>% suppressWarnings() %>% suppressMessages()
          t2 <- Sys.time()
          
          # apply base model to test
          ml.test.results = data.frame(predicted_response = predict(ml.model, newdata = test.data))
          ml.shadow.results = data.frame(shadow_response = predict(ml.model, newdata = shadow.data))
        } # end base model
        
        # optionally: build stack
        if(algo == "stack"){ # for stacked ensemble
          set.seed(iter)
          model_list <- caretEnsemble::caretList(
            (value) ~., 
            data = full.train.data,
            metric = "RMSE", 
            methodList = c("ranger", "glmnet", "svmRadial"),
            #trace=FALSE,
            preProc = c("center", "scale", "zv"),
            trControl = caret::trainControl(method="cv", number=5) 
          ) %>% suppressWarnings()
          # note: ranger returns numbers as.character(); must fix with as.numeric()
          model_list$ranger$pred$pred <- as.numeric(model_list$ranger$pred$pred)
          
          set.seed(iter)
          # build meta-model
          ml.model <- caretEnsemble::caretStack(
            model_list,
            method = "glm", 
            metric = "RMSE",
            trControl = caret::trainControl(method="cv", number=5))
          # apply ensemble model to test
          ml.test.results = data.frame(predict(ml.model, newdata = test.data))
          colnames(ml.test.results) = "predicted_response"
          ml.shadow.results = data.frame(predict(ml.model, newdata = shadow.data))
          colnames(ml.shadow.results) = "shadow_response"
          
        } # end ensemble model
        
        # combine and melt results
        ml.results = data.frame(prediction = c(ml.test.results$predicted_response, ml.shadow.results$shadow_response))
        
        # fill out data
        ml.results$true.response = test.response$value
        ml.results$data.type = c("real", "shadow")
        ml.results$RS_Name = x
        ml.results$HM = hm
        ml.results$algo = algo
        ml.results$iter = iter
        ml.results$target = z
        
        t2 <- Sys.time()
        
        ml.results$time = t2-t1
        
        print(paste(x, " ", iter, " ", z, " ", hm, " ", algo, " ", t2-t1))
        
        ml.results
      #}))}))})) 
  }, mc.cores=core.n))}))}))

t2 = Sys.time()

print(paste("Completed:  algo =", algo, "  Iter =", iter, round(t2-t1, digits=3), "  FS =", feature_select))
ml.loocv.regression.results$elapsed = t2 - t1
ml.loocv.regression.results$feature_selected = ifelse(feature_select == TRUE, TRUE, FALSE)

saveRDS(ml.loocv.regression.results, paste0(home.dir, "/ml_git_data/ml_loocv_data/2025_11_07_ml_loocv_regression_results", 
                                            ifelse(feature_select == TRUE, "_fs_", "_"), algo, "_", iter, ".rds", sep=""))
print("Complete LOOCV run.")
