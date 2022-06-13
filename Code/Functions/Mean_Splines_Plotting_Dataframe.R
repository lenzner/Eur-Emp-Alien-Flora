#=========================================================================#
# Naturalized alien floras still carry the legacy of European colonialism #
# Lenzner et al.                                                          #
# Functions                                                               #
# Developed by Bernd Lenzner                                              #
# Contact: bernd.lenzner@univie.ac.at                                     #
#=========================================================================#

msgdm_spline <- function(model, zeta.order, predictor, values){
  # model = MSGDM model
  # zeta.order = zeta order of interest (= original zeta order - 1)
  # predictor = predictor of interest
  # values = either extract Ispline values ("Ispline") or rescaled predictor values ("env")
  model[[zeta.order]] %>% # Select model output for one zeta oder where i indicates the zeta order (i ranges form 1-4, representing zeta orders 2-5)
    map(~ .x[[values]]) %>% # Get $env object from each list (each list is one repetition; there are currently 6 repetitions)
    map(~ .x[predictor]) %>% # Get only the first column from each Ispline object
    bind_cols(.name_repair = "unique") %>% # unlist all Ispline objects into one dataframe and create unique column names
    rowMeans() # calculate row means
  }


mean_splines <- function(model, zeta.order, predictor, emp.name){
  
  emp = emp.name #gsub("mod.ispline.","",paste(print(substitute(model)))) # Extract empire name from the model object name
  zeta = paste0("Zeta ", zeta.order+1)
  
  tibble(
    empire = emp,
    zeta = zeta,
    predictor = colnames(model[[zeta.order]][[1]]$env[predictor]), # Extract predictor name
    spline = msgdm_spline(model = model, zeta.order = zeta.order, predictor = predictor, values = "Ispline"), # extract mean Ispline values
    pred.range = msgdm_spline(model = model, zeta.order = zeta.order, predictor = predictor, values = "env") # extract rescaled predictor values
    )
}


Plot_mean_isplines <- function(models, num.zeta, pred.num){
  
  models <- models
  names(models) <- c("mod.ispline.GBR", "mod.ispline.ESP", "mod.ispline.PRT", "mod.ispline.NED")
  
  num.zeta = num.zeta
  pred.num = pred.num
  mod_name <- c("GBR", "ESP", "PRT", "NED")
  df_m <- list()
  
    for(m in 1:length(models)){
      
      df_z <- list()
      emp <- mod_name[m] # set name of empire 
      
        for(z in 1:num.zeta){
        
        df_l <- list() # Build empty list. Each predictor will be stored in one list object
        
          for(p in 1:pred.num){ # loop through nuber of predictors
            
            df_l[[p]] <- mean_splines(model = models[[m]], zeta.order = z, predictor = p, emp.name = emp) # calculate mean splines per predictor per zeta order and store in list object
            
            }
          
        df_z[[z]] <- bind_rows(df_l) # unlist into long data tibble
        
        }
      
      df_zeta <- bind_rows(df_z)
      
      
      df_m[[m]] <- df_zeta # store in list object for the respective empire
      
      
      
    }
  
  df_final <- bind_rows(df_m)
  
  # Renaming and final edits to dataframe (e.g., variable names, empire names, ...)
  df_final$empire <- as.factor(df_final$empire)
  df_final$empire <- sub("GBR", "British Empire", df_final$empire)
  df_final$empire <- sub("ESP", "Spanish Empire", df_final$empire)
  df_final$empire <- sub("PRT", "Portuguese Empire", df_final$empire)
  df_final$empire <- sub("NED", "Dutch Empire", df_final$empire)
  df_final$empire <- factor(df_final$empire, levels = c("British Empire", "Spanish Empire", "Portuguese Empire", "Dutch Empire"))
  
  df_final$predictor <- sub("area", "Area", df_final$predictor)
  df_final$predictor <- sub("^aridity.*", "Aridity index", df_final$predictor)
  df_final$predictor <- sub("^crop.*", "Agricultural land", df_final$predictor)
  df_final$predictor <- sub("^dem.*", "Habitat heterogeneity", df_final$predictor)
  df_final$predictor <- sub("distance", "Geographic distance", df_final$predictor)
  df_final$predictor <- sub("^gdp.*", "GDPpc", df_final$predictor)
  df_final$predictor <- sub("^pop.*", "Population density", df_final$predictor)
  df_final$predictor <- sub("^temp.*", "Mean annual temperature", df_final$predictor)
  df_final$predictor <- sub("^time.*", "Occupation time", df_final$predictor)
  
  df_final <- unique(df_final)
  
  return(df_final)
}




