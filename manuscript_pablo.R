source("scripts/init_pablo.R")

# Auxiliar functions

coordTrans <- function(coords, from.epsg, to.epsg){
  coordinates(coords) <- c(1,2)
  proj4string(coords) <- CRS(paste("+init=epsg:", from.epsg, sep=""))
  df.trans <- data.frame(spTransform(coords, CRS(paste("+init=epsg:", to.epsg, sep=""))))
  names(df.trans) <- c("x", "y")
  return(df.trans)
}

getObservedGrowthData <- function(){
  df.growth <- read_excel(f.growth) %>% 
    dplyr::select(long=Long, lat=Lat, growth=Inc.AB5Anya)
  xy <- coordTrans(df.growth[, c("long", "lat")], "4326", "3035")
  names(xy) <- c("x", "y")
  df.growth <- cbind(df.growth, xy) %>% dplyr::select(x, y, growth)
  return(df.growth)    
}

getOccurrencesByGroupAndPercentile <- function(){
  df <- bind_rows(
    read.csv(f.modelling.data) %>% filter(p==1 & group == 'Continental') %>% 
      dplyr::select(x, y, group, growth=growth_AB) %>% 
      mutate(percentile = ntile(growth, n.percentiles)),
    read.csv(f.modelling.data) %>% filter(p==1 & group == 'Mild') %>% 
      dplyr::select(x, y, group, growth=growth_AB) %>% 
      mutate(percentile = ntile(growth, n.percentiles)),
    read.csv(f.modelling.data) %>% filter(p==1 & group == 'Sp') %>% 
      dplyr::select(x, y, group, growth=growth_AB) %>% 
      mutate(percentile = ntile(growth, n.percentiles))) %>% 
    mutate(percentile=as.character(percentile))
  return(df)
}

getSuitabilityModels <- function(perc = 1){
  df.res <- data.frame()
  for(m in names(models)){
    for(g in c("Continental", "Mild", "Sp", "SpC")){
      f <- paste0(output.dir, "maxentPred/", m, "/", g, "/eqOrGreaterThanPerc", perc, "/species_predictors.asc")
      cat("Processing", f, "\n")
      df <- data.frame(as(raster(f), "SpatialPixelsDataFrame")) %>% setNames(., c("v", "x", "y")) %>% 
        mutate(model = m, group=g)
      df.res <- df.res %>% bind_rows(df)
    }
  }
  return(df.res)
}

getSuitabilityBestModels <- function(bestMdls){
  df.res <- data.frame()
  bestMdls <- bestMdls %>% mutate(group =  recode(group,"SP" = "Sp", "C" = "Continental", "M" = "Mild", "M+C" = "SpC"))
  for(i in 1:length(bestMdls$model)){
    m <- bestMdls[i, "model"]
    g <- bestMdls[i, "group"]
    perc <- bestMdls[i, "percentile"]
    f <- paste0(output.dir, "maxentPred/", m, "/", g, "/eqOrGreaterThanPerc", perc, "/species_predictors.asc")
    cat("Processing", f, "\n")
    df <- data.frame(as(raster(f), "SpatialPixelsDataFrame")) %>% setNames(., c("v", "x", "y")) %>% 
        mutate(model = m, group=g)
    df.res <- df.res %>% bind_rows(df)
  }
  return(df.res)
}

numextract <- function(string){ 
  str_extract(string, "\\-*\\d+\\.*\\d*")
} 

## TABLES

# Models results
bestMdlFun <- function(){
  bestMdls <- data.frame()
  df <- read.csv(f.maxent.models.results)
  df$Occurrences <- df$percentile
  df[df$group == "SpC", ]
  df <- df %>% mutate(Occurrences = recode(Occurrences, "1" = "100%", "2" = "90%", "3" = "80%", "4" = "70%", "5" = "60%", "6" = "50%", "7" = "40%", "8" = "30%", "9" = "20%", "10" = "10%"),
                      group =  recode(group,"Sp" = "SP", "Continental" = "C", "Mild" = "M", "SpC" = "M+C"))
  # All models results
  allMdls <- df[, c("model", "group", "Occurrences", "training.auc", "training.occ.nr", "training.bg.nr", "test.auc", "test.occ.nr")]
  write.csv(allMdls, paste0(manuscript.dir, "tables/supplementary/allModels.csv"), row.names = F)
  cat("all models results saved on ===> ", manuscript.dir, "tables/supplementary/allModels.csv \n")
  # All models variable permutation importance
  df.var <- read.csv(f.maxentPred.models.results)
  allPerm <- df.var[, c("model", "group", "percentile", "predictors.permutation.importance")]
  write.csv(allPerm, paste0(manuscript.dir, "tables/supplementary/allModels_permutationImportance.csv"), row.names = F)
  cat("all models permutation importance results saved on ===> ", manuscript.dir, "tables/supplementary/allModels.csv \n")
  
  df.growth <- read.csv(f.suit.growth.corr)
  df.growth <- df.growth[which(df.growth$growth.dataset == "predicted"),]
  df.compl <- cbind(df, "growth_cor" = df.growth$rho.spearman, "gr_cor_pvalue" = df.growth$p.value)
  df.compl$test.auc <- allMdls$test.auc
  df.compl$mdl_name <- paste0(df$model, "_", df$group)
  for(mdl in unique(df.compl$mdl_name)){
    mdlGroup <- df.compl[which(df.compl$mdl_name == mdl), ]
    mdlGroup <- mdlGroup[which(mdlGroup$growth_cor > 0), ]
    if(any(mdlGroup$gr_cor_pvalue < 0.05)){
    mdlGroup <- mdlGroup[which(mdlGroup$gr_cor_pvalue < 0.05), ]
    } else {
      mdlGroup <- mdlGroup[which(mdlGroup$gr_cor_pvalue == min(mdlGroup$gr_cor_pvalue)), ] 
    }
    bestMdl <- mdlGroup[which(mdlGroup$growth_cor == max(mdlGroup$growth_cor)), ]
    if(bestMdl$test.auc < 0.9){
    bestMdl <- mdlGroup[which(mdlGroup$test.auc == max(mdlGroup$test.auc)), ] 
    }
    bestMdls <- rbind(bestMdls, bestMdl)
    print(paste0(bestMdl$mdl_name, "_",bestMdl$percentile,". test.AUC: ", bestMdl$test.auc, ". Growth_cor:", bestMdl$growth_cor, ". p_value:", bestMdl$gr_cor_pvalue))
  }
  # Table manuscript
  # Non filtered and best models comparation
  df.perc1 <- df.compl[which(df.compl$percentile == 1), ]
  df.perc1 <- df.perc1[, c("model", "group", "Occurrences", "training.auc", "growth_cor", "gr_cor_pvalue")]
  bestMdlsT <- bestMdls[, c("model", "group", "Occurrences", "training.auc", "growth_cor", "gr_cor_pvalue")]
  table <- rbind(df.perc1, bestMdlsT)
  names(table) <- c("Model", "Group", "Occurrences", "AUC", "Correlation", "p-value")
  write.csv(table, paste0(manuscript.dir, "tables/results_table.csv"), row.names = F)
  cat("results table comparison models saved on ===> ", manuscript.dir, "tables/results_table.csv \n")
  return(bestMdls)
}

# Permutation importance table
permutation_importanceTableFun <- function(){
  df <- read.csv(f.maxentPred.models.results, header = T)
  bMdls <- bestMdlFun()
  bMdls <- bMdls %>% mutate(group =  recode(group,"SP" = "Sp", "C" = "Continental", "M" = "Mild", "M+C" = "SpC"))
  clim.res <- data.frame()
  climbio.res <- data.frame()
  for(m in names(models)){
    pi.df <- as.data.frame(df[which(df$model == m), ])
    for(g in groups){
      pi_gr.df <- as.data.frame(pi.df[which(pi.df$group == g), ])
      bestPerc <- bMdls %>% filter(group == g, model == m)
      bestPerc <- bestPerc$percentile
      for(perc in c(1, bestPerc)){
        pi.vec <- as.vector(pi_gr.df[which(pi_gr.df$percentile == perc), "predictors.permutation.importance"])
        pi <- unlist(str_split(pi.vec, ", "))
        if(m == "Clim"){
          clim.df <- data.frame("model" = m, "group" = g, "percentile" = perc, "AP" = numextract(pi[1]), "MeWiT" = numextract(pi[2]), "PS" = numextract(pi[3]), "SuP" = numextract(pi[4]), "TS" = numextract(pi[5]))
          clim.res <- rbind(clim.res, clim.df)
        } else {
          climbio.df <- data.frame("model" = m, "group" = g, "percentile" = perc, "AP" = numextract(pi[1]), "MeWiT" = numextract(pi[2]), "PS" = numextract(pi[3]), "SuP" = numextract(pi[4]), "TS" = numextract(pi[5]),
                                   "que_hum_curr" = numextract(pi[6]), "que_pyr_curr" = numextract(pi[7]), "pin_pina_curr" = numextract(pi[8]), "pin_nig_curr" = numextract(pi[9]),  "fag_syl_curr" = numextract(pi[10]))
          climbio.res <- rbind(climbio.res, climbio.df)
        }
      }
    }
  }
  clim.res <- clim.res %>% mutate(percentile = recode(percentile, "1" = "100%", "2" = "90%", "3" = "80%", "4" = "70%", "5" = "60%", "6" = "50%", "7" = "40%", "8" = "30%", "9" = "20%", "10" = "10%"),
                      group =  recode(group,"Sp" = "SP", "Continental" = "C", "Mild" = "M", "SpC" = "M+C"))
  write.csv(clim.res, paste0(manuscript.dir, "tables/permutation_importance_clim.csv"), row.names = T)
  print(paste0("====>", manuscript.dir, "tables/permutation_importance_clim.csv"))
  climbio.res <- climbio.res %>% mutate(percentile = recode(percentile, "1" = "100%", "2" = "90%", "3" = "80%", "4" = "70%", "5" = "60%", "6" = "50%", "7" = "40%", "8" = "30%", "9" = "20%", "10" = "10%"),
                                     group =  recode(group,"Sp" = "SP", "Continental" = "C", "Mild" = "M", "SpC" = "M+C"))
  write.csv(climbio.res, paste0(manuscript.dir, "tables/permutation_importance_climbio.csv"), row.names = T)
  print(paste0("====>", manuscript.dir, "tables/permutation_importance_climbio.csv"))
}

# All models growth
allModelsGrowth <- function(){
  df <- read.csv("outputs/maxentPred/suitability-growth_corr.csv") %>% mutate(sig = as.character(ifelse(p.value < 0.05, 1, 0))) %>%
    mutate(sig = recode(sig, "1" = "p.value <= 0.05", "0" = "p.value > 0.05"),group =  recode(group,"Sp" = "SP", "Continental" = "C", "Mild" = "M", "SpC" = "M+C"))
  df <- df %>% filter(growth.dataset == "predicted")
  df$Occurrences <- df$above.perc
  df <- df %>% mutate(Occurrences = recode(Occurrences, "1" = "100%", "2" = "90%", "3" = "80%", "4" = "70%", "5" = "60%", "6" = "50%", "7" = "40%", "8" = "30%", "9" = "20%", "10" = "10%"))
  allCor <- df[, c("group", "model", "Occurrences", "n", "rho.spearman", "p.value")]
  names(allCor) <- c("Group", "Model", "Occurrences", "N", "Spearman_correlation_coefficient", "p-value")
  write.csv(allCor, paste0(manuscript.dir, "tables/supplementary/allModels_cor.csv"), row.names = F)
  cat("===> tables/supplementary/allModels_cor.csv \n")
}

## FIGURES 

# Occurrence distribution 
plotOccurrencesDistribution <- function(){
  f.out <- paste0(manuscript.dir, "figures/OccurrencesByGroup/occurrences_by_group.png")
  df <- read.csv(f.modelling.data) %>% filter(p==1) %>% 
    filter(group != 'Sp') %>% 
    dplyr::select(group, x, y)
  png(f.out, width=1900, height=1482, units = 'px') # For 8 cm wide at 600 dpi
  p <- ggplot() + geom_point(data=df, aes(x=x, y=y, col = group), shape = 20, size = 8)
  obsPoints <- getObservedGrowthData()
  obsGr.df <- read.csv("outputs/maxent/pbg_sets/Sp_ClimBio_eqOrGreaterThanPerc1_species_test.csv", header = T)
  p <- p + scale_color_manual(values=c("#00B1A0", "#DE3D3D"))
  p <- p + geom_path(data=ip.df, aes(x=x, y=y)) + coord_equal()
  p <- p +geom_point(x = obsPoints$x, y = obsPoints$y, data = obsPoints, shape = 17, size = 10, alpha = 0.9)
  p <- p + theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
                 axis.text = element_blank(), panel.grid = element_blank(), legend.position = 'none')
  print(p)
  dev.off()
  cat("==>", f.out, "\n")
}

plotOccurrencesByGrowthPercentileDistribution <- function(){
  for(perc in 1:10){
    df <- getOccurrencesByGroupAndPercentile()
    df <- df %>% filter(group != 'Sp', as.numeric(percentile) >= perc)
    count.df <- df %>% mutate(group = recode(group, "Continental" = "C", "Mild" = "M")) %>% group_by(group) %>% summarize(count = n())
    count.df <- rbind(count.df, cbind("group" = "SP", "count"  = sum(count.df$count)))
    write.csv(count.df, paste0(manuscript.dir, "figures/OccurrencesByGrowth/numberOccurrencesByGroup_",perc, ".csv"))
    f.out <- paste0(manuscript.dir, "figures/OccurrencesByGrowth/occurrrences_by_growth_perc_", perc,".png")
    png(f.out, height = 6.5, width = 6.5,res = 1000, units = "cm") # For 8 cm wide at 600 dpi
    p <- ggplot(df) + geom_point(aes(x=x, y=y, colour=group), shape=20, size=1.5, alpha = 0.7)
    p <- p + scale_color_manual(values=c("#00B1A0", "#DE3D3D"))
    p <- p + geom_path(data=ip.df, aes(x=x, y=y)) + coord_equal()
    p <- p + theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
                   axis.text = element_blank(), panel.grid = element_blank(), legend.position = 'none',
                   strip.text.x = element_text(size = 22))
    # p <- p + facet_wrap(~ group)
    print(p, res = 1000)
    dev.off()
    cat("==>", f.out, "\n") 
  }
}


# Percentiles per group
plotNumberOfOccurrencesPerGroupAndPercentile <- function(){
  df2 <- getOccurrencesByGroupAndPercentile() %>% group_by(group, percentile) %>% summarise(n=n()) %>% mutate(percentile=as.integer(percentile))
  df.res <- data.frame()
  for(g in groups){
    for(perc in 1:n.percentiles){
      n <- df2 %>% filter(group==g & percentile >= perc) %>% pull(n) %>% sum()
      df.res <- rbind(df.res, cbind(group=g, perc=perc, n=n))
    }
  }
  df.res <- df.res %>% mutate(n=as.integer(as.character(n)), group=as.character(group), perc = as.integer(perc))
  f.out <- paste0(manuscript.dir, "figures/number_occurrrences_per_growth_and_percentile.png")
  png(f.out, width=1890, height=1228, units = 'px') # For 8 cm wide at 600 dpi
  p <- ggplot(df.res) + geom_point(aes(x=perc, y=n, group=group, colour=group), size=7)
  p <- p + scale_color_manual(values=c("#00B1A0", "#DE3D3D", "grey30"))
  p <- p + geom_text(aes(x=perc, y=n, label=n), hjust=-0.3, vjust=0.5, size=10)
  p <- p + xlab("\nPercentage of occurrences considered according to growth values")
  p <- p + ylab("\nNumber of occurrences\n")
  p <- p + theme(panel.grid.minor = element_blank(),
                 axis.title = element_text(size=24),
                 axis.text.x = element_text(size=24, vjust=-0.5),
                 axis.text.y = element_text(size=24),
                 axis.ticks.length = unit(0.3, 'cm'),
                 legend.title = element_blank(), 
                 legend.text = element_text(size = 26),
                 legend.position = c(0.92,0.9), legend.background = element_rect(fill = "white"),
                 legend.key.size = unit(1.5, 'cm'),
                 legend.key = element_rect(colour = NA, fill = NA))
  print(p)
  dev.off()
  cat("==>", f.out, "\n")
}


# To plot single maps showing best model projection

# Plot maps suitability

plotSuitabilityModelsByGroup <- function(){
  df <- getSuitabilityModels(perc = 1)
  df$occurrences <- "100%"
  bestMdls <- bestMdlFun()
  df2 <- getSuitabilityBestModels(bestMdls)
  df2$occurrences <- "Best models"
  dfAll <- rbind(df, df2)
  dfAll <- dfAll %>% mutate(group = recode(group, "Continental" = "C", "Mild" = "M", "Sp" = "SP", "SpC" = "M+C"))
  f.out <- paste0(manuscript.dir, "figures/SuitabilityMaps/suitability_models_by_model_percentile_and_group.png")
  png(f.out, width=3780/2, height=2279/2, units = 'px') # For 8 cm wide at 600 dpi
  pdf(f.out)
  p <- ggplot(dfAll) + geom_raster(aes(x=x, y=y, fill=v))
  p <- p + scale_fill_gradientn("", colours = c("snow1", "yellow3", "green", "darkgreen"), limits=c(0,1))
  p <- p + coord_equal()
  p <- p + facet_grid(occurrences + model ~  group)
  p <- p + theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 axis.title = element_blank(),
                 axis.text = element_blank(),
                 axis.ticks = element_blank(),
                 strip.text = element_text(colour="black", size=20),
                 # strip.background = element_rect(fill= "white"),
                 panel.background = element_rect(fill= "#ebebebff"),
                 legend.position = 'bottom',
                 legend.text = element_text(size = 30))
  p <- p + guides(fill = guide_colourbar(barwidth = 60, barheight = 2))
  print(p)
  dev.off()
  cat("==>", f.out, "\n")
  # pdf
  f.out <- paste0(manuscript.dir, "figures/SuitabilityMaps/suitability_models_by_model_percentile_and_group.pdf")
  pdf(f.out, width=20, height=20)
  print(p)
  dev.off()
  cat("==>", f.out, "\n")
  return(bestMdls)
}

# Plot range sizes
plotRangeSizeByModelAndGroup <- function(){
  for(perc in c(1, 7)){
  f.out <- paste0(manuscript.dir, "figures/RangeSize/range_size_by_model_and_group_", perc, ".png")
  png(f.out, width=1890, height=1228, units = 'px') # For 8 cm wide at 600 dpi
  df <- read.csv(f.range.size) %>% filter(above.perc == perc) %>% 
    mutate(area = round(area, 0), group = recode(group, "Continental" = "C", "Mild" = "M", "Sp" = "SP", "SpC" = "C+M"))
  df <- left_join(data.frame(group = c("C", "C", "M", "M", "C+M", "C+M", "SP", "SP")),    # Reorder data frame
                         df,
                         by = "group")
  df$group <- factor(df$group, levels = unique(df$group))
  p <- ggplot(df, aes(x = group, y = area, fill = model)) 
  p <- p + geom_bar(stat = 'identity', position = position_dodge())
  p <- p + geom_text(aes(x = group, y = area, label = format(area, big.mark=" ")), 
                     stat = 'identity', position = position_dodge(width=1), vjust=-0.50, hjust=0.5, size=12)
  # p <- p + xlab("\nGroup") + ylab("Range size\n")
  p <- p + scale_fill_manual(values=c("#af8dc3", "#7fbf7b"))
  p <- p + theme(panel.grid.minor = element_blank(),
                 legend.position = "none",
                 legend.title = element_blank(),
                 axis.ticks.length = unit(0.3, 'cm'),
                 axis.title.x = element_blank(),
                 axis.title.y = element_blank(),
                 axis.text.x = element_text(size=40, vjust=-0.5),
                 axis.text.y = element_text(size=40))
  p <- p + ylim(0, 100000)
  p <- p + scale_color_npg()
  print(p)
  dev.off()
  cat("==>", f.out, "\n")
  }
}

plotRangeSizeBestModels <- function(){
  f.out <- paste0(manuscript.dir, "figures/RangeSize/range_size_by_model_and_group_bestModels.png")
  png(f.out, width=1890, height=1228, units = 'px') # For 8 cm wide at 600 dpi
  bMdls <- bestMdlFun()
  bMdls$area <- numeric(length(bMdls$model))
  bMdls <- bMdls %>% mutate(group = recode(group, "C" = "Continental", "M" = "Mild", "M+C" = "SpC", "SP" = "Sp"))
  for(mdl in bMdls$mdl_name){
    range.df <- read.csv(f.range.size) %>% filter(group == bMdls[which(bMdls$mdl_name == mdl), "group"], model == bMdls[which(bMdls$mdl_name == mdl), "model"], above.perc == bMdls[which(bMdls$mdl_name == mdl), "percentile"])
    bMdls[which(bMdls$mdl_name == mdl), "area"] <- round(range.df$area, 0)
  }
  df <- bMdls %>% mutate(group = recode(group, "Continental" = "C", "Mild" = "M", "Sp" = "SP", "SpC" = "C+M"))
  df <- left_join(data.frame(group = c("C", "M", "C+M", "SP")),    # Reorder data frame
                  df,
                  by = "group")
  df$group <- factor(df$group, levels = unique(df$group))
    p <- ggplot(df, aes(x = group, y = area, fill = model)) 
    p <- p + geom_bar(stat = 'identity', position = position_dodge())
    p <- p + geom_text(aes(x = group, y = area, label = format(area, big.mark=" ")), 
                       stat = 'identity', position = position_dodge(width=1), vjust=-0.50, hjust=0.5, size=12)
    # p <- p + xlab("\nGroup") + ylab("Range size\n")
    p <- p + scale_fill_manual(values=c("#af8dc3", "#7fbf7b"))
    p <- p + theme(panel.grid.minor = element_blank(),
                   legend.position = "none",
                   legend.title = element_blank(),
                   axis.ticks.length = unit(0.3, 'cm'),
                   axis.title.x = element_blank(),
                   axis.title.y = element_blank(),
                   axis.text.x = element_text(size=40, vjust=-0.5),
                   axis.text.y = element_text(size=40))
    p <- p + scale_color_npg()
    p <- p + ylim(0, 100000)
    print(p)
    dev.off()
    cat("==>", f.out, "\n")
}


# Growth suitability correlation
plotGrowthSuitabilityCorrelation <- function(){
  df <-  read.csv("outputs/maxentPred/suitability-growth_corr.csv") %>% mutate(sig = as.character(ifelse(p.value < 0.05, 1, 0))) %>% 
    mutate(sig = recode(sig, "1" = "p.value <= 0.05", "0" = "p.value > 0.05"), group = recode(group, "Continental" = "C", "Mild" = "M", "Sp" = "SP", "SpC" = "C+M"))
  df <- left_join(data.frame(group = c("C", "M", "C+M", "SP")),    # Reorder data frame
                  df,
                  by = "group")
  df$group <- factor(df$group, levels = unique(df$group))
  df$sig <- as.factor(df$sig)
  df$model <- as.factor(df$model)
  df <- df %>% filter(growth.dataset == "predicted")
  f.out <- paste0(manuscript.dir, "figures/GrowthSuitabilityCorrelation/growth_suitability_correlation.png")
  png(f.out, width=1800/3, height=1800/3, units = 'px') # For 8 cm wide at 600 dpi
  p <- ggplot(df) + geom_point(aes(x=above.perc, y=rho.spearman, shape=sig, colour=model), size=4, alpha = 0.7)  + scale_shape_manual(values = c(17, 16))
  p <- p + scale_color_manual(values=c("#af8dc3", "#7fbf7b"))
  p <- p + scale_y_continuous(limits=c(-0.3, 1))
  p <- p + scale_x_discrete(limits = c("100", "90", "80", "70", "60", "50", "40", "30", "20", "10"))
  # p <- p + xlab("\nPercentage according to growth values")
  # p <- p + ylab("Spaearman correlation coefficient\n")
  p <- p + theme(legend.position = 'none',
                 panel.grid.minor = element_blank(),
                 legend.title = element_blank(),
                 axis.title.x = element_blank(),
                 axis.title.y = element_blank(),
                 axis.text.x = element_text(size=12, vjust=-0.5),
                 axis.text.y = element_text(size=12),
                 strip.text = element_text(colour="black", size=14),
                 # strip.background = element_rect(fill= "grey40"),
                 panel.background = element_rect(fill= "#ebebebff")

  )
  p <- p + facet_wrap(~group, ncol = 2) # ~ growth.dataset
  print(p)
  dev.off()
  cat("==>", f.out, "\n")
  
  f.out <- paste0(manuscript.dir, "figures/GrowthSuitabilityCorrelation/growth_suitability_correlation.pdf")
  pdf(f.out)
  print(p)
  dev.off()
  cat("==>", f.out, "\n")
}


# AUC plot
plotTestAUCByGroupAndModel <- function(){
  df <- read.csv(f.maxent.models.results) %>% mutate(group = recode(group, "Continental" = "C", "Mild" = "M", "Sp" = "SP", "SpC" = "M+C"))
  df$percentile <- df$percentile * 10
  f.out <- paste0(manuscript.dir, "figures/AUCs/models_test_auc.png")
  png(f.out, width=2400, height=1320, units = 'px') # For 8 cm wide at 600 dpi
  p <- ggplot(df)
  p <- p + geom_point(aes(x=percentile, y=test.auc, colour=group), shape=20, size=10)
  # p <- p + stat_smooth(aes(x=percentile, y=test.auc, colour = group), method = "gam")
  p <- p + ylab("Test AUC\n")
  p <- p + xlab("\nPercentage of occurrences considered according to growth values")
  p <- p + theme(panel.grid.minor = element_blank(),
                 axis.title = element_text(size=24),
                 axis.text.x = element_text(size=24, vjust=-0.5),
                 axis.text.y = element_text(size=24),
                 axis.ticks.length = unit(0.3, 'cm'),
                 strip.text = element_text(size=36),
                 strip.background = element_rect(fill= "grey80"),
                 legend.title = element_blank(), 
                 legend.text = element_text(size = 26),
                 legend.position = "bottom", legend.background = element_rect(fill = "white"),
                 legend.key.size = unit(1.5, 'cm'),
                 legend.key = element_rect(colour = NA, fill = NA))
  p <- p + facet_wrap(~ model)
  print(p)
  dev.off()
  cat("==>", f.out, "\n")
}


# Growth projeciton plot
growthPlotFun <- function(){
  gr.df <- read.csv(paste0(data.processed.dir, "MeAGD.csv"), header = T)
  obsPoints <- getObservedGrowthData()
  # obsGr.df <- read.csv("outputs/maxent/pbg_sets/Sp_ClimBio_eqOrGreaterThanPerc1_species_test.csv", header = T)
  f.out <- paste0(manuscript.dir, "figures/GrowthProjection/growthProjection.png")
  png(f.out, width=1900/2, height=1482/2, units = 'px') # For 8 cm wide at 600 dpi
  p <- ggplot(gr.df) + geom_raster(aes(x=x, y=y, fill=Growth_AB))
  p <- p + scale_fill_continuous(type = "viridis")
  # p <- p + scale_fill_gradient("Relative value", low = "Moccasin", high = "Dark Green", na.value = "grey50")
  # p <- p +geom_point(x = obsPoints$x, y = obsPoints$y, data = obsPoints, shape = 20, size = 8)
  p <- p + coord_equal()
  p <- p + theme(panel.grid = element_blank(),
                 axis.title = element_blank(),
                 axis.text = element_blank(),
                 axis.ticks = element_blank(),
                 legend.title = element_blank())
  print(p)
  dev.off()
  print(paste("==>", f.out))
}

# VARIABLE RESPONSE
respData.fun <- function(plotPredictors = c("AP", "MeWiT")){
  whole.df <- data.frame()
  bMdls <- bestMdlFun() %>% filter(group != "SpC")
  bMdls <- bMdls %>% mutate(group = recode(group, "C" = "Continental", "M" = "Mild", "SP" = "Sp", "M+C" = "SpC"))
  for(model in names(models)){
    predictors <- models[[model]]
    for(group in groups){
      for(predictor in predictors){
        # Perc1
        perc <- 1
        file.data <- paste0("outputs/maxentPred/", model, "/", group, "/eqOrGreaterThanPerc", perc,"/plots/species_", predictor, ".dat")
        predictor.df <- read.csv(file.data, header = T)
        predictor.df$model <- as.factor(model)
        predictor.df$group <- as.factor(group)
        predictor.df$percentile <- as.factor(perc)
        # Perc2
        perc2 <- bMdls[which(bMdls$model == model), ]
        perc2 <- perc2[which(perc2$group == group), "percentile"]
        file.data <- paste0("outputs/maxentPred/", model, "/", group, "/eqOrGreaterThanPerc", perc2,"/plots/species_", predictor, ".dat")
        predictorPerc2.df <- read.csv(file.data, header = T)
        predictorPerc2.df$model <- as.factor(model)
        predictorPerc2.df$group <- as.factor(group)
        predictorPerc2.df$percentile <- as.factor(perc2)
        whole.df <- bind_rows(whole.df, predictor.df, predictorPerc2.df)
      } # predictors
    } # group
  } # model
  whole.df$model <- as.factor(whole.df$model)
  whole.df$variable <- as.factor(whole.df$variable)
  whole.df$x <- as.numeric(whole.df$x)
  whole.df$y <- as.numeric(whole.df$y)
  # PLot
  for(predictor in plotPredictors){
    for(model in names(models)){
      for(group in groups){
        varPlot.df <- whole.df[whole.df$variable == predictor, ]
        varPlot.df <- varPlot.df[varPlot.df$model == model, ]
        varPlot.df <- varPlot.df[varPlot.df$group == group, ]
        f.out <- paste0("manuscript/figures/ClimVariableResponse/", model, "_", group, "_", predictor, "_response.pdf")
        pdf(f.out, width=2, height=3)
        p <- ggplot(data = varPlot.df, aes(x = x, y = y, col = percentile))

          if(predictor == "MeWiT" && group == "Continental"){
            p <- p + xlim(min(varPlot.df$x), 3.5)
          }
        
        if(predictor == "MeWiT" && group == "Mild"){
          p <- p + xlim(3.5, max(varPlot.df$x))
        }
        p <- p + geom_line(lty = 1, lwd = 2)
        p <- p + scale_color_grey()
        p <- p + theme(axis.title = element_blank(), legend.position = 'none')
        p <- p + facet_wrap(~ group)
        print(p, res = 1000)
        dev.off()
        print(paste("==>", f.out))
      }
    }
  }
}

respDataBio.fun <- function(plotPredictors = c("fag_syl_curr", "pin_pina_curr", "que_hum_curr")){
  whole.df <- data.frame()
  bMdls <- bestMdlFun() %>% filter(group != "SpC")
  bMdls <- bMdls %>% mutate(group = recode(group, "C" = "Continental", "M" = "Mild", "SP" = "Sp", "M+C" = "SpC"))
  model <- "ClimBio"
  predictors <- models[[model]]
    for(group in groups){
      for(predictor in predictors){
        # Perc1
        perc <- 1
        file.data <- paste0("outputs/maxentPred/", model, "/", group, "/eqOrGreaterThanPerc", perc,"/plots/species_", predictor, ".dat")
        predictor.df <- read.csv(file.data, header = T)
        predictor.df$model <- as.factor(model)
        predictor.df$group <- as.factor(group)
        predictor.df$percentile <- as.factor(perc)
        # Perc2
        perc2 <- bMdls[which(bMdls$model == model), ]
        perc2 <- perc2[which(perc2$group == group), "percentile"]
        file.data <- paste0("outputs/maxentPred/", model, "/", group, "/eqOrGreaterThanPerc", perc2,"/plots/species_", predictor, ".dat")
        predictorPerc2.df <- read.csv(file.data, header = T)
        predictorPerc2.df$model <- as.factor(model)
        predictorPerc2.df$group <- as.factor(group)
        predictorPerc2.df$percentile <- as.factor(perc2)
        whole.df <- bind_rows(whole.df, predictor.df, predictorPerc2.df)
      } # predictors
    } # group
  whole.df$model <- as.factor(whole.df$model)
  whole.df$variable <- as.factor(whole.df$variable)
  whole.df$x <- as.numeric(whole.df$x)
  whole.df$y <- as.numeric(whole.df$y)
  # PLot
  for(predictor in plotPredictors){
      varPlot.df <- whole.df[whole.df$variable == predictor, ]
      f.out <- paste0("manuscript/figures/ClimBioVariableResponse/ClimBio", "_", predictor, "_response.pdf")
      pdf(f.out, width=6, height=3)
      p <- ggplot(data = varPlot.df, aes(x = x, y = y, col = percentile))
      p <- p + geom_line(lty = 1, lwd = 2)
      p <- p + scale_color_grey()
      p <- p + theme(axis.title = element_blank(), legend.position = 'none')
      p <- p + facet_wrap(~ group)
      print(p, res = 1000)
      dev.off()
      print(paste("==>", f.out))
  }
}


# Best models table
bMdls <- bestMdlFun()

#### FIGURES ------------------------------------------------------------------------------------ ####

### FIGURE 1. Occurrence distribution

## a) Group
plotOccurrencesDistribution()

## b) Group and growth (also Fig. S2)
plotOccurrencesByGrowthPercentileDistribution()


### FIGURE 4. Suitability maps

plotSuitabilityModelsByGroup()


### FIGURE 5. Range size

# 100%
plotRangeSizeByModelAndGroup()

# Best models
plotRangeSizeBestModels()


### FIGURE 6. Growth suitability correlation

# Same test dataset (Fig. 6)
plotGrowthSuitabilityCorrelation()

# Random points (Fig. S4)
# plotGrowthSuitabilityCorrelationRandPoint()


### FIGURE S1. Number of occurences

growthPlotFun()


### FIGURE 4. AUCs

plotTestAUCByGroupAndModel()


### FIGURE S5. Biotic variable response

respDataBio.fun(plotPredictors = c("fag_syl_curr", "pin_pina_curr", "que_hum_curr"))


### FIGURE S6. Climate variable response

respData.fun(c("AP", "MeWiT"))


### Others

# plotNumberOfOccurrencesPerGroupAndPercentile()

#### TABLES --------------------------------------------------------------------------------------------------------------------------####

### TABLE 2.

permutation_importanceTableFun()

### TABLE 3.

## Best models results

bMdls <- bestMdlFun()

### Table S3.

allModelsGrowth()


#### OUTPUTS TO FURTHER USE -------------------------------------------------------------------------------------------------------------------------- ####

df.res <- data.frame()

for(i in 1:10){
  df <- getSuitabilityModels(perc = i)
  df$growthPercentile <- i
  df.res <- rbind(df.res, df)
  }

tail(df.res)
