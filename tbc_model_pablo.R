source("scripts/init_pablo.R")

#### AUXILIARY FUNCTIONS ##########################################################################
# Coordinates transformation
coordTrans <- function(coords, from.epsg, to.epsg){
  coordinates(coords) <- c(1,2)
  proj4string(coords) <- CRS(paste("+init=epsg:", from.epsg, sep=""))
  df.trans <- data.frame(spTransform(coords, CRS(paste("+init=epsg:", to.epsg, sep=""))))
  names(df.trans) <- c("x", "y")
  return(df.trans)
}

# Transform tifs to ascs for Maxent
transformPredictors <- function(predictors){
  for(predictor in predictors){
    f.asc <- paste0("data_processed/predictors/", predictor, ".asc")
    if(file.exists(f.asc))
      next
    f.tif<- paste0(data.processed.dir, predictor, ".tif")
    r <- raster(f.tif)
    writeRaster(r, f.asc)
  }    
}

getModellingSetFileName <- function(group, model.name, percentile, speciesOrBackground, testOrTrain){
  if(speciesOrBackground == 'background'){
    f <- paste0(paste0(output.dir, "maxent/pbg_sets/"), group, "_", speciesOrBackground, "_", testOrTrain, ".csv")  
  } else {
    f <- paste0(paste0(output.dir, "maxent/pbg_sets/"), group, "_", model.name, "_eqOrGreaterThanPerc", percentile, "_", 
                speciesOrBackground, "_", testOrTrain, ".csv")
  }
  return(f)
}

#### DATA PREPARATION #############################################################################
#### CLIMATE AND BIO PREDICTORS ----------------------------------------------------------------- #
generatePredictorLayers <- function(){
  # Reference raster
  cat("Preparing predictor layers in TIF and ASC formats ...\n")
  r.ref <- raster(f.raster.ref)
  for(i in 1:dim(predictors)[1]){
    r <- raster(predictors[i, "original.file"])
    if(predictors[i, "type"]=='clim'){
      r <- resample(r, r.ref, method='ngb')
    } else {
      r <- projectRaster(r, r.ref)
    }
    dir.create(paste0(data.processed.dir, "predictors/"), recursive = T, showWarnings = F)
    f.out <- paste0(data.processed.dir, "predictors/", predictors[i, "processed.file"], ".tif")
    writeRaster(r, f.out, overwrite=T)
    cat("  ==>", f.out, "\n")
    f.out <- paste0(data.processed.dir, "predictors/", predictors[i, "processed.file"], ".asc")
    writeRaster(r, f.out, overwrite=T)
    cat("  ==>", f.out, "\n")
  }
}
#### GROWTH MODEL ------------------------------------------------------------------------------- #

# Regression formula provided by MR on 2017.07.05
generateGrowthModel <- function(){
  cat("Generating growth model (Inc.AB <- -240.409 + 0.453 * PP_Spring + 11.021 * AutTmax) ...\n")
  cat("  calculating spring precipitation ...\n")
  MarP <- raster(paste0(dcaip.dir, "prec3.tif"))/10
  AprP <- raster(paste0(dcaip.dir, "prec4.tif"))/10
  MayP <- raster(paste0(dcaip.dir, "prec5.tif"))/10
  PP_Spring <- MarP + AprP + MayP
  
  cat("  calculating mean fall temperature ...\n")
  SepT <- raster(paste0(dcaip.dir, "tmax9.tif"))/10
  OctT <- raster(paste0(dcaip.dir, "tmax10.tif"))/10
  NovT <- raster(paste0(dcaip.dir, "tmax11.tif"))/10
  AutTmax <- mean(SepT, OctT, NovT)
  
  cat("  calculating and writing model ...\n")
  Inc.AB <- -240.409 + 0.453 * PP_Spring + 11.021 * AutTmax
  f.out <- paste0(data.processed.dir, "MeAGD.tif")
  writeRaster(Inc.AB, f.out, overwrite=T)
  cat("    ==>", f.out, "\n")
  
  df.MeAGD <- data.frame(as(Inc.AB, "SpatialPixelsDataFrame"))
  colnames(df.MeAGD) <- c("Growth_AB", "x", "y")
  f.out <- paste0(data.processed.dir, "MeAGD.csv")
  write.csv(df.MeAGD, f.out, row.names = F)
  cat("    ==>", f.out, "\n")
}
#### MODELLING DATASET -------------------------------------------------------------------------- #
generateModellingDataset <- function(){
  cat("Preparing occurrence modelling dataset ...\n")
  cat("  reading occurrences raw data ...\n")
  occ <- bind_rows( # CREAF occurrences
    read_excel(f.tbc.occ.creaf, sheet="TAXBAC") %>% 
      dplyr::select(long=LONG, lat=LAT, location=TOPONIMO, res=RES_INIC),
    # University of Santiago occurrences
    read_excel(f.tbc.occ.usc) %>% 
      mutate(res="puntual") %>% 
      dplyr::select(long=`Arbol-CoordX`, lat=`Arbol-CoordY`, location=Población, res))
  occ <- occ %>% filter(!is.na(long) & !is.na(lat) & res %in% c("puntual","1km"))
  occ <- bind_cols(occ, coordTrans(occ[, c("long", "lat")], "4326", "3035")) %>% 
    dplyr::select(long, lat, x, y, location) %>% 
    # Filter out Balearic Islands
    filter(!(x > 3600000 & y < 2000000))
  cat("Total number of occurrences for the Iberian Peninsula: ", length(occ$x))
  cat("  separating occurrences into adaptive groups ...\n")
  ### Adaptive group separation based on Winter Temperature
  MeWiT <- raster(paste0(data.processed.dir, "predictors/MeWiT.tif"))
  MeWiT_value <- raster::extract(MeWiT, occ[,c("x","y")])
  occ <- cbind(occ[,c("x", "y")], MeWiT_value)
  
  # Resolve T values for pixels that fall outside MeWiT raster
  occ.na <- occ[is.na(occ$MeWiT_value),]
  occ <- occ[!is.na(occ$MeWiT_value),]
  occ.na$MeWiT_value <- raster::extract(MeWiT, occ.na[,c("x","y")], buffer = 5000, fun=mean, na.rm=TRUE)
  occ <- rbind(occ, occ.na)
  
  # Group populations by T distance to experimental site
  T_threshold <- 3.5
  
  occ$group <- ifelse(occ$MeWiT_value <= T_threshold, "Continental", "Mild")
  occ <- occ[, c("x", "y", "group")]
  
  cat("  filtering occurrences to one per grid cell ...\n")
  # Keep one occurrence per grid cell
  occ$x<-(floor(occ$x/1000)*1000)+500
  occ$y<-(floor(occ$y/1000)*1000)+500
  occ <- unique(occ[,c("x","y", "group")])
  cat("Number of occurrences once filtered (1/km2) ", length(occ$x), 
      "\nContinental: ", length(occ[which(occ$group == "Continental"), "x"]),
      "\nMild: ", length(occ[which(occ$group == "Mild"), "x"]))
  
  cat("  adding predictor and estimated growth values to occurrences ...\n")
  tifs <- c(paste0(data.processed.dir, "predictors/", predictors$processed.file, ".tif"),
            paste0(data.processed.dir, "MeAGD.tif"))
  st <- stack(tifs)
  
  layers <- c(predictors$processed.file, "MeAGD")
  occ.pred <- raster::extract(st, occ[, c("x", "y")])
  occ.pred <- cbind(occ, occ.pred)
  for(i in 1:length(layers)){
    cat("    processing", layers[i], "\n")
    for(j in 1:dim(occ.pred)[1]){
      pop <- occ.pred[j,]
      if(is.na(occ.pred[j, layers[i]])){
        v <- raster::extract(st[[layers[i]]], occ.pred[j, c("x", "y")], buffer=4500, na.rm=T, fun = mean)
        if(is.na(v)){
          cat("  WARNING:", paste(occ.pred[j, c("x", "y")]), "still is NA, maybe the extracting buffer needs to be increased\n")
        }
        occ.pred[j, layers[i]] <- v
      }
    }
  }
  
  cat("  readying environmental stack as data frame ...\n")  
  f.out <- paste0(data.processed.dir, "climate_spdf_ip.csv")
  spdf.pred <- data.frame(as(st, "SpatialPixelsDataFrame"))
  
  cat("  generating final modelling dataset as csv file ...\n")
  occ.tbc <- occ.pred
  occ.tbc$p <- 1
  occ.tbc <- occ.tbc[, c("p", "x", "y", "group", layers)]
  set.seed(123)
  bg.tbc <- sample_n(spdf.pred, 10000)
  bg.tbc$p <- 0
  # Prepare dataset for modelling
  df.tbc <- bind_rows(
    occ.tbc %>% dplyr::select(p, x, y, group, layers) %>% mutate(group='Sp'),
    bg.tbc %>% mutate(group='Sp') %>% dplyr::select(p, x, y, group, layers),
    occ.tbc %>% filter(group == 'Continental') %>% dplyr::select(p, x, y, group, layers),
    bg.tbc %>% mutate(group='Continental') %>% dplyr::select(p, x, y, group, layers),
    occ.tbc %>% filter(group == 'Mild') %>% dplyr::select(p, x, y, group, layers),
    bg.tbc %>% mutate(group='Mild') %>% dplyr::select(p, x, y, group, layers)
  ) %>% rename(growth_AB=MeAGD)
  
  write.csv(df.tbc, f.modelling.data, row.names = F)
  cat("  ==>", f.modelling.data, "\n")
}
#### CHECK CORRELATIONS AMONG PREDICTORS -------------------------------------------------------- #
checkPredictorCorrelations <- function(){
  cat("Checking predictors correlations ...\n")
  occ.env <- read.csv(f.modelling.data) %>% filter(group == 'Sp' & p == 1) %>% 
    dplyr::select(predictors$processed.file)
  cor_Vsp <- cor(occ.env)
  f.out <- paste0(output.dir, "cor_sp_current.csv")
  write.csv(cor_Vsp, f.out, row.names=F, quote=F)
  cat("  ==>", f.out, "\n")
}

checkPredictorCorrelationsBg <- function(){
  cat("Checking predictors correlations including also background points ...\n")
  occ.env <- read.csv(f.modelling.data) %>% filter(group == 'Sp' & p == 1) %>% 
    dplyr::select(predictors$processed.file)
  cor_Vsp <- cor(occ.env)
  f.out <- paste0(output.dir, "cor_sp_current_bg.csv")
  write.csv(cor_Vsp, f.out, row.names=F, quote=F)
  cat("  ==>", f.out, "\n")
}

#### MODEL FITTING AND EVALUATION #################################################################
fitMaxentModels <- function(groups, models, n.percentiles, maxent.args, force.run=F){
  cat("\nFitting Maxent models ...")
  dir.create(paste0(output.dir, "maxent/pbg_sets"), recursive=T, showWarnings=F)
  df <- read.csv(f.modelling.data)
  for(g in groups){
    for(m in 1:length(models)){
      cat("  \n----------------------------------------------------------------------------------------------------------------\nmodels of", g, names(models[m]), "\n")
      # Presences per quartile
      df.pres <- df %>% filter(group == g & p == 1) %>% mutate(percentile = ntile(growth_AB, n.percentiles))
      # Background
      df.bg <- df %>% filter(group == g & p == 0)
      
      # Save background as csv files for later referene
      cat("    preparing and writing training and test csv files for background ...\n")
      set.seed(1234)
      df.bg.train <- df.bg
      # We take number of background points for testing in the same proportion as the number of training points to 10000 background points
      # However, I tested results at different nr of background points (10000, 50000, 10000) for testing and Test AUC does not vary.
      indx.bg.test <- sample(1:dim(df.bg)[1], round((1-0.7) * dim(df.bg)[1] / 0.7, 0))
      df.bg.test <- df.bg[indx.bg.test,]
      
      f.bg.train <- getModellingSetFileName(g, names(models[m]), '-', 'background', 'train')
      f.bg.test <- getModellingSetFileName(g, names(models[m]), '-', 'background', 'test')
      write.csv(df.bg.train, f.bg.train, quote=F, row.names=F)
      cat("      --> training bg:", f.bg.train, "\n")
      write.csv(df.bg.test, f.bg.test, quote=F, row.names=F)
      cat("      --> testing bg:", f.bg.test, "\n\n")
      
      # If model already exists do next model, unless force.run=T
      for(perc in 1:n.percentiles){
        f.in <- paste0(output.dir, "maxent/", names(models[m]), "/", g, "/eqOrGreaterThanPerc", perc, "/maxentResults.csv")
        f.inPred <- paste0(output.dir, "maxentPred/", names(models[m]), "/", g, "/eqOrGreaterThanPerc", perc, "/maxentResults.csv")
        if(file.exists(f.in) & file.size(f.in) != 0 & file.exists(f.inPred) & file.size(f.inPred) != 0 & !force.run){
          cat(paste("    --> !!", g, names(models[m]), perc, "exists, skipping ...\n"))
          next
        }
        
        cat(paste("    percentile >=", perc, "\n"))
        df.pres.perc <- df.pres[df.pres$percentile >= perc, names(df.pres) != "percentile"]
        # Separate training, test and background and save files for further use
        cat("      preparing and writing training and test csv files for presences ...\n")
        set.seed(1234)
        indx.sp.train <- sample(1:dim(df.pres.perc)[1], round(dim(df.pres.perc)[1] * 0.7, 0))
        df.sp.train <- df.pres.perc[indx.sp.train,]
        df.sp.test <- df.pres.perc[-indx.sp.train,]
        
        f.sp.train <- getModellingSetFileName(g, names(models[m]), perc, 'species', 'train')
        f.sp.test <- getModellingSetFileName(g, names(models[m]), perc, 'species', 'test')
        write.csv(df.sp.train, f.sp.train, quote=F, row.names=F)
        cat("        --> training presences:", f.sp.train, "\n")
        write.csv(df.sp.test, f.sp.test, quote=F, row.names=F)
        cat("        --> training presences:", f.sp.test, "\n")
        
        # Vector of presences and background
        pbg <- c(rep(1, nrow(df.sp.train)), rep(0, nrow(df.bg.train)))
        
        # Data frame of environmental data for presences and background
        df.mod <- rbind(df.sp.train, df.bg.train)
        
        # Where maxent should store its outputs
        maxent.results.dir <- paste0(output.dir, "maxent/", names(models[m]), "/", g, "/eqOrGreaterThanPerc", perc)
        maxentPred.results.dir <- paste0(output.dir, "maxentPred/", names(models[m]), "/", g, "/eqOrGreaterThanPerc", perc)
        dir.create(maxent.results.dir, recursive=T, showWarnings = F)
        dir.create(maxentPred.results.dir, recursive=T, showWarnings = F)
        
        # Maxent model with 30% data held as test data
        cat("      running Maxent ...\n")
        maxent.model <- dismo::maxent(x=df.mod[, models[[m]]], p=pbg, path=maxent.results.dir, maxent.args)
        
        f.out <- paste0(output.dir, "maxent/maxent_model_", names(models[m]), "_", g, "_eqOrGreaterThanPerc", perc, ".rds")
        saveRDS(maxent.model, f.out)
        cat(paste("      ==>", f.out, "\n\n"))
        
        cat("      running Maxent using all presences ...\n")
        head(df.bg)
        df.mod <- rbind(df.pres.perc, df.bg)
        
        maxent.pred.model <- dismo::maxent(x=df.mod[, models[[m]]], p=df.mod$p, path=maxentPred.results.dir, maxent.args)
        
        f.out <- paste0(output.dir, "maxentPred/maxent_prediction_model_", names(models[m]), "_", g, "_eqOrGreaterThanPerc", perc, ".rds")
        saveRDS(maxent.pred.model, f.out)
        cat(paste("      ==>", f.out, "\n\n"))
      }
    }
  }
}
# Generate species composites from adaptive groups
generateSpeciesComposites <- function(){
  cat("\nGenerating species composites from adaptive groups ...")
  for(model in names(models)){
    for(perc in 1:n.percentiles){
      cat("\n  ", model, perc)
      fs <- paste0(output.dir, "maxentPred/", model, "/", groups[groups != 'Sp'], "/", 
                   paste0("eqOrGreaterThanPerc", perc), "/species_predictors.asc")
      if(paste0(file.exists(fs), collapse = "") != "TRUETRUE"){
        cat("    ! one of these files does not exist: ", fs)
      }
      st <- stack(fs)
      r <- max(st)
      dir.out <- paste0(output.dir, "maxentPred/", model, "/SpC/eqOrGreaterThanPerc", perc, "/")
      dir.create(dir.out, recursive = T, showWarnings = F)
      f.out <- paste0(dir.out, "species_predictors.asc")
      writeRaster(r, f.out, overwrite = T)
      cat("\n   ==>", f.out, "\n")
    }
  }
}

evaluateMaxentModels <- function(groups, models, n.percentiles, force.run=F){
  cat("\nEvaluating models with test files ...\n")
  for(g in groups){
    f.bg.test <- getModellingSetFileName(g, names(models[m]), '-', 'background', 'test')
    for(m in 1:length(models)){
      for(perc in 1:n.percentiles){
        f.model <-  paste0(output.dir, "maxent/maxent_model_", names(models[m]), "_", g, "_eqOrGreaterThanPerc", perc, ".rds")
        f.sp.test <- getModellingSetFileName(g, names(models[m]), perc, 'species', 'test')
        if(!file.exists(f.model)){
          cat("  --> !! file", f.model, "does not exist, skipping ...\n")
          next()
        }
        if(!file.exists(f.sp.test)){
          cat("  --> !! file", f.sp.test, "does not exist, skipping ...\n")
          next()
        }
        if(!file.exists(f.bg.test)){
          cat("  --> !! file", f.bg.test, "does not exist, skipping ...\n")
          next()
        }
        df.sp.test <- read.csv(f.sp.test)
        df.bg.test <- read.csv(f.bg.test)
        mx.model <- readRDS(f.model)
        df.sp.test.predictions <- predict(mx.model, df.sp.test)
        df.bg.test.predictions <- predict(mx.model, df.bg.test)
        eval <- evaluate(p=df.sp.test.predictions, a=df.bg.test.predictions)
        f.out <- paste0(output.dir, "maxent/maxent_model_evaluation_", names(models[m]), "_", g, "_eqOrGreaterThanPerc", perc, ".rds")
        saveRDS(eval, f.out)
        cat("  ==>", f.out, "\n")
      }
    }
  }
  # Evaluate the SpC composites with the same test and background files used for the Species
  for(m in 1:length(models)){
    for(perc in 1:n.percentiles){
      f.r <- paste0(output.dir, "maxent/", names(models[m]), "/SpC/eqOrGreaterThanPerc", perc, "/species_predictors.asc")
      r <- raster(f.r)
      occ.test <- read.csv(getModellingSetFileName("Sp", names(models[m]), perc, 'species', 'test')) %>% 
        dplyr::select(x, y)
      bg.test <- read.csv(getModellingSetFileName("Sp", names(models[m]), perc, 'background', 'test')) %>% 
        dplyr::select(x, y)
      v.occ <- raster::extract(r, occ.test)
      v.bg <- raster::extract(r, bg.test)
      eval <- evaluate(p=v.occ, a=v.bg)
      f.out <- paste0(output.dir, "maxent/maxent_model_evaluation_", names(models[m]), "_SpC_eqOrGreaterThanPerc", perc, ".rds")
      saveRDS(eval, f.out)
      cat("  ==>", f.out, "\n")
    }
  }
}

prepareMaxentModelsResults <- function(groups, models, n.percentiles, force.run = T){
  cat("\nGathering results from all models ...\n")
  if(file.exists(f.maxent.models.results) & !force.run){
    return(read.csv(f.maxent.models.results))
  }
  else {
  df.results <- data.frame()
  for(g in c("Continental", "Mild", "Sp", "SpC")){
    for(m in 1:length(models)){
      for(perc in 1:n.percentiles){
        cat("  -->", g, m, perc, "\n")
        if(!isTRUE(g == "SpC")){
          f.in <- paste0(output.dir, "maxent/", names(models[m]), "/", g, "/eqOrGreaterThanPerc", perc, "/maxentResults.csv") 
        }
        f.ev.in <- paste0(output.dir, "maxent/maxent_model_evaluation_", names(models[m]), "_", g, "_eqOrGreaterThanPerc", perc, ".rds")
        model <- names(models[m])
        predictors.permutation.importance <- NA
        training.auc <- NA
        training.occ.nr <- NA
        training.bg.nr <- NA
        test.auc <- NA
        test.occ.nr <- NA
        test.bg.nr <- NA
        test.cor <- NA
        if(!isTRUE(g == "SpC")){
          if(!file.exists(f.in) | file.size(f.in) == 0){
            cat("      !", f.in, "does not exist. Skipping ...\n")
            next
          } else {
            maxent.results <- read.csv(f.in)
          }
        }
        if(!file.exists(f.ev.in) | file.size(f.ev.in) == 0){
          cat("      !", f.ev.in, "does not exist. Skipping ...\n")
          next
        }
        evaluation.results <- readRDS(f.ev.in)
        if(!isTRUE(g == "SpC")){
        predictors.permutation.importance <- paste(models[[m]], maxent.results[,paste0(models[[m]], ".permutation.importance")], collapse=", ")
        training.auc <- round(maxent.results$Training.AUC, 3)
        training.occ.nr <- maxent.results$X.Training.samples
        training.bg.nr <- maxent.results$X.Background.points
        }
        test.auc <- round(evaluation.results@auc, 3)
        test.occ.nr <- evaluation.results@np
        test.bg.nr <- evaluation.results@na
        
        results <- cbind("model"=names(models[m]),
                         "group"=g,
                         "percentile"=perc,
                         "predictors"= paste(models[[m]], collapse=", "),
                         "predictors.permutation.importance"=predictors.permutation.importance,
                         "training.auc" = training.auc,
                         "training.occ.nr" = training.occ.nr,
                         "training.bg.nr" = training.bg.nr,
                         "test.auc" = test.auc,
                         "test.occ.nr" = test.occ.nr,
                         "test.bg.nr"= test.bg.nr)
        df.results <- rbind(df.results, results)
      }
    }
  }
  write.csv(df.results, f.maxent.models.results, quote = T, row.names = F)
  cat("  ==>", f.maxent.models.results, "\n")
  }
}

prepareMaxentPredPermImport <- function(groups, models, n.percentiles, force.run = T){
  cat("\nGathering results from all models ...\n")
  if(file.exists(f.maxentPred.models.results) & !force.run){
    return(read.csv(f.maxentPred.models.results))
  }
  else {
    df.results <- data.frame()
    for(g in c("Continental", "Mild", "Sp", "SpC")){
      for(m in 1:length(models)){
        for(perc in 1:n.percentiles){
          cat("  -->", g, m, perc, "\n")
          if(!isTRUE(g == "SpC")){
            f.in <- paste0(output.dir, "maxentPred/", names(models[m]), "/", g, "/eqOrGreaterThanPerc", perc, "/maxentResults.csv") 
          }
          model <- names(models[m])
          predictors.permutation.importance <- NA
          test.cor <- NA
          if(!isTRUE(g == "SpC")){
            if(!file.exists(f.in) | file.size(f.in) == 0){
              cat("      !", f.in, "does not exist. Skipping ...\n")
              next
            } else {
              maxent.results <- read.csv(f.in)
            }
          }
          if(!isTRUE(g == "SpC")){
            predictors.permutation.importance <- paste(models[[m]], maxent.results[,paste0(models[[m]], ".permutation.importance")], collapse=", ")
          }
          
          results <- cbind("model"=names(models[m]),
                           "group"=g,
                           "percentile"=perc,
                           "predictors"= paste(models[[m]], collapse=", "),
                           "predictors.permutation.importance"=predictors.permutation.importance)
          df.results <- rbind(df.results, results)
        }
      }
    }
    write.csv(df.results, f.maxentPred.models.results, quote = T, row.names = F)
    cat("  ==>", f.maxentPred.models.results, "\n")
  }
}


getObservedGrowthData <- function(){
  df.growth <- read_excel(f.growth) %>% 
    dplyr::select(long=Long, lat=Lat, growth=Inc.AB5Anya)
  xy <- coordTrans(df.growth[, c("long", "lat")], "4326", "3035")
  names(xy) <- c("x", "y")
  df.growth <- cbind(df.growth, xy) %>% dplyr::select(x, y, growth)
  return(df.growth)    
}
# Returns points with predicted growth data at the given test dataset
# Here is where we decide which points use for the correlation. For now, I am using all occurrences per group (e.g., 1817 for species)
getPredictedGrowthData <- function(group, model.name, samePerc = T, percentile, onlyTest = F){
  if(group == 'SpC'){
    group <- 'Sp'
  }
  if(samePerc == T){
    percentile <- 1
  }
  if(onlyTest == T){
    f.test <- getModellingSetFileName(group, model.name, percentile, "species", "test")
    df <- read.csv(f.test) %>% dplyr::select(x,y,growth=growth_AB)
  } else {
    f.test <- getModellingSetFileName(group, model.name, percentile, "species", "test")
    f.train <- getModellingSetFileName(group, model.name, percentile, "species", "train")
    df <- read.csv(f.test) %>% dplyr::select(x,y,growth=growth_AB)
    df2 <- read.csv(f.train) %>% dplyr::select(x,y,growth=growth_AB)
    df <- rbind(df, df2)
  }
  return(df)
}
# Returns a dataset with columns: group, model, above.perc, suitability, growth.ds, growth
getGrowthSuitabilityDataSet <- function(groups, models, n.percentiles, force.run = T){
  cat("\nGetting dataset for analysing relation growth-suitability ...\n")
  
  f.out <- paste0(output.dir, "maxentPred/suitability-growth_data.csv")
  if(file.exists(f.out) & !force.run)
    return(read.csv(f.out))
  
  df.growth.observed <- getObservedGrowthData()
  
  df.res <- data.frame()
  for(g in c(groups, "SpC")){
    for(m in names(models)){
      for(perc in 1:n.percentiles){
        cat("  processing growth data for model", g, names(models[m]), perc, "...\n")
        f <- paste0(output.dir, "maxentPred/", names(models[m]), "/", g, "/eqOrGreaterThanPerc", perc, "/species_predictors.asc")
        if(!file.exists(f)){
          cat("    --> !!! Model does not exist.\n")
          next
        }
        # Growth-suitability relation with 25 growth observations
        r.suit <- raster(f)          
        v <- raster::extract(r.suit, df.growth.observed[, c("x", "y")])
        df.go <- data.frame(group=rep(g, length(v)), model=rep(names(models[m]), length(v)), above.perc=rep(perc, length(v)), 
                            x=df.growth.observed$x, y=df.growth.observed$y, 
                            suitability=v, growth.dataset=rep('observed', length(v)), growth=df.growth.observed$growth)
        
        # Growth-suitability relation with test data points
        df.growth.predicted <- getPredictedGrowthData(g, names(models[m]), percentile = perc)
        v <- raster::extract(r.suit, df.growth.predicted[, c("x", "y")])
        df.pr <- data.frame(group=rep(g, length(v)), model=rep(names(models[m]), length(v)), above.perc=rep(perc, length(v)), 
                            x=df.growth.predicted$x, y=df.growth.predicted$y, 
                            suitability=v, growth.dataset=rep('predicted', length(v)), growth=df.growth.predicted$growth)
        
        df.res <- rbind(df.res, df.go, df.pr)
      }
    }
  }
  write.csv(df.res, f.out, quote = F, row.names = F)
  cat("  ==>", f.out, "\n")
  return(df.res)
}
# Returns a data frame with correlation test data between growth and suitability
generateGrowthSuitabilityCorr <- function(groups, models, n.percentiles, force.run = T){
  cat("\nTesting growth ~ suitability Spearman correlation for all models ...\n")
  if(file.exists(f.suit.growth.corr) & !force.run)
    return(read.csv(f.suit.growth.corr))
  
  df.growth.suit <- getGrowthSuitabilityDataSet(groups, models, n.percentiles, force.run)
  
  df.res <- data.frame()
  for(g in c(groups, "SpC")){
    for(m in names(models)){
      for(perc in 1:n.percentiles){
        for(gds in c('observed', 'predicted')){
          df <- df.growth.suit %>% filter(group == g & model == names(models[m]) 
                                          & above.perc == perc & growth.dataset == gds)
          df <- na.omit(df)
          if(dim(df)[1]==0)
            next
          cor <- cor.test(df$growth, df$suitability, method='spearman')
          cor.df <- data.frame(cbind(group=g, model=m, above.perc=perc, n=nrow(df), growth.dataset=gds, 
                                     statistic=round(cor$statistic,4), 
                                     p.value=round(cor$p.value, 4), 
                                     rho.spearman=round(cor$estimate, 4)))
          df.res <- rbind(df.res, cor.df)
        }
      }
    }
  }
  
  write.csv(df.res, f.suit.growth.corr, quote = F, row.names = F)
  cat("  ==>", f.suit.growth.corr, "\n")
}

# Generate growth correlation using backgraound random points covering the whole study area for Sp and SpC
generateGrowthSuitabilityCorrRandPoint <- function(){
  df.res <- data.frame()
  growth.predicted <- raster("data_processed/MeAGD.tif")
  for(g in c("Sp", "SpC")){
    for(m in names(models)){
      for(perc in 1:n.percentiles){
        gr <- g
        if(gr == 'SpC'){
          gr <- 'Sp'
        }
        occ <- read.csv(f.modelling.data, header = T) %>% filter(p == 0, group == gr)
        occ <- sample_n(occ, 1000)
        v.growthPred <- cbind(occ[, c("x", "y")], "growth" = raster::extract(growth.predicted, occ[, c("x", "y")]))
        f <- paste0(output.dir, "maxent/", names(models[m]), "/", g, "/eqOrGreaterThanPerc", perc, "/species_predictors.asc")
        #Predicted
        r.suit <- raster(f)          
        df <- cbind(v.growthPred, "suitability" = raster::extract(r.suit, occ[, c("x", "y")]))
        # Density (to check normality)
        ggplot(df) + geom_density(aes(growth))
        ggplot(df) + geom_density(aes(suitability))
        # Correlations
        cor <- cor.test(df$growth, df$suitability, method='spearman')
        cor.df <- data.frame(cbind(group=g, model=m, above.perc=perc, n=nrow(df), growth.dataset = "predicted",
                                   statistic=round(cor$statistic,4), 
                                   p.value=round(cor$p.value, 4), 
                                   rho.spearman=round(cor$estimate, 4)))
        df.res <- rbind(df.res, cor.df)
      }
    }
  }
  write.csv(df.res, "outputs/maxent/suitability-growth_corr_randPoints.csv", row.names = F)
  cat("  ==> outputs/maxent/suitability-growth_corr_RandPoints.csv\n")
}

# Calculates range sizes for all models
generateRangeSizes <- function(force.run=F){
  cat("\nCalculating range sizes for all models\n")
  if(file.exists(f.range.size) & !force.run)
    return(read.csv(f.range.size))
  
  df.res <- data.frame()
  for(g in c(groups, "SpC")){
    for(m in names(models)){
      for(perc in 1:n.percentiles){
        cat("  ", g, names(models[m]), perc, "...\n")
        f <- paste0(output.dir, "maxentPred/", names(models[m]), "/", g, "/eqOrGreaterThanPerc", perc, "/species_predictors.asc")
        if(!file.exists(f)){
          cat("     --> !!! Model does not exist.\n")
          next
        }
        r.suit <- raster(f)
        v <- na.omit(getValues(r.suit))
        area <- round(sum(v),2)
        df.res <- rbind(df.res, cbind(group=g, model=names(models[m]), above.perc=perc, area=area))
      }
    }
  }
  df.res$group <- as.character(df.res$group)
  df.res$model <- as.character(df.res$model)
  df.res$above.perc <- as.character(df.res$above.perc)
  df.res$area <- as.numeric(as.character(df.res$area))
  write.csv(df.res, f.range.size, quote = F, row.names = F)
  cat("  ==>", f.range.size, "\n")
}

tic()
# --- Data preparation spends about 3' in execution time --- #
generatePredictorLayers()
generateGrowthModel()
generateModellingDataset()
checkPredictorCorrelations()
checkPredictorCorrelationsBg()

# --- Model fitting spends about 1h 45' in execution time --- #
fitMaxentModels(groups, models, n.percentiles, maxent.args, F)
generateSpeciesComposites()
evaluateMaxentModels(groups, models, n.percentiles, T)

variablesImpMaxentPredModels(groups, models, n.percentiles, T, testOrPred.dir = paste0(output.dir, testOrPred.dir, "/maxent_prediction_model_", names(models[m]), "_", g, "_eqOrGreaterThanPerc", perc, ".rds"))

# --- Results preparation spends about 1' in execution time --- #
prepareMaxentModelsResults(groups, models, n.percentiles, T)

prepareMaxentPredPermImport(groups, models, n.percentiles, T)

generateGrowthSuitabilityCorr(groups, models, n.percentiles, T)

generateRangeSizes(T)
toc()
