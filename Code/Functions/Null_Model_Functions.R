#=========================================================================#
# Naturalized alien floras still carry the legacy of European colonialism #
# Lenzner et al.                                                          #
# Functions                                                               #
# Developed by Bernd Lenzner                                              #
# Contact: bernd.lenzner@univie.ac.at                                     #
#=========================================================================#


# Randomization - Null model

## Perfroming a randomized draw of regions for a specific empire.
## The total number of colonial regions per empire remains constant
## The number of colonial regions per UN georegion remains constant
## The number of mainland and island colonial regions remains constant

# Datasets ----

## Load datasets from "analysis" Script


# Null model 1 ----
# The null model selects random regions globally. The regions selection is constrained by the number of regions in a geospatial region in the original empire and by the number of mainland and island regions in the geospatial region. Besically the null model builds a random empire with the same spatial structure (i.e. number of mainland and island regions per geospatial region as in the observed empire)

random.empire <- function(colonizer = c("Great_Britain", "Portugal", "Spain", "France", "Netherlands"), iteration = iteration){
  
  COL <-  c("Great_Britain", "Portugal", "Spain", "France", "Netherlands")
  if(!colonizer %in% COL) stop("empire not found")
  
  
  emp <- empires[which(empires$colonizer == colonizer),]    # select empire
  emp <- unique(emp[,c(1,3,4,5,8)])                      # remove duplicates based on multiple colonial affiliation periods
  n.emp <- dim(emp)[1]                                   # number of colonial regions
  n.emp.geospat <- table(emp$UN_geospat)                 # number of colonial regions per geospatial region
  n.emp.reg <- table(emp$MainIsl)                        # number of colonial regions that are islands and mainlands
  n.emp.reg.geospat <- table(emp[,c("UN_geospat", "MainIsl")])
  
  un <- rownames(n.emp.reg.geospat)                                                                  # UN geospatial regions covered by the empire
  mi <- colnames(n.emp.reg.geospat)                                                                  # islands - mainland identifier
  
  
  
  result <- list()
  
  for(k in 1:iteration){
    
    res <- data.frame(OBJIDsic = character(),
                      region_id = character(),
                      code = character(),
                      name = character(),
                      island = character(),
                      IDregion = character(),
                      regionF = character(),
                      LAT = character(),
                      LON = character(),
                      GeodAREA = character(),
                      UN_geospat = character())
    
    for(i in 1:length(un)){
      
      reg <- regions[which(regions$UN_geospat == un[i]),]                                                       # subset regions to a specific UN geospatial region
      reg_ran <- list()
      
      for(j in 1:length(mi)){
        
        reg_ran[[j]] <- sample_n(reg[which(reg$isl == mi[j]),], n.emp.reg.geospat[un[i],mi[j]], replace = F)    # random sample of island/mainland regions per UN gespatial region based on the observed numbers
        
      }
      
      dat <- do.call(rbind,reg_ran)                                                                            # rbind random sample of islands and mainlands for the specific geospatial region
      
      res <- rbind(res, dat)
      
    }
    
    result[[k]] <- res
    
  }
  
  return(result)
  
}


# Null model 2 ----

# The null model selects random regions globally. 
### The regions selection is constrained by the number of regions in a geospatial region in the original empire and by the number of mainland and island regions in the geospatial region.
### The first region per geospatial region is selected randomly, all subsequent regions will be selected sequentially based on the shortest distance to the other regions. The idea is that regions are colonized sequentially 
### Besically the null model builds a random empire with the same spatial structure (i.e. number of mainland and island regions per geospatial region as in the observed empire) but the selection of regions is not totally random as in the first null model

random.empire.seq <- function(colonizer = c("Great_Britain", "Portugal", "Spain", "France", "Netherlands"), iteration = iteration){
  
  # Check if colonizer is defined in the function
  COL <-  c("Great_Britain", "Portugal", "Spain", "France", "Netherlands")
  # Error message if provided "colonizer"-name differs from input options
  if(!colonizer %in% COL) stop("empire not found")
  
  
  
  # subset empire
  emp <- empires[which(empires$colonizer == colonizer),]    
  # remove duplicates based on multiple colonial affiliation periods
  emp <- unique(emp[,c(1,3,4,5,8)])
  # extract number of colonial regions
  n.emp <- dim(emp)[1]
  # extract number of colonial regions per geospatial region
  n.emp.geospat <- table(emp$UN_geospat)
  # extract number of clonial regions that are islands and mainlands
  n.emp.reg <- table(emp$MainIsl)
  # extract number of islands and mainland regions per UN geospatial region
  n.emp.reg.geospat <- table(emp[,c("UN_geospat", "MainIsl")])
  # extract UN geospatial regions covered by the empire
  un <- rownames(n.emp.reg.geospat)
  # create an islands - mainland identifier
  mi <- colnames(n.emp.reg.geospat)
  
  
  
  result <- list()
  
  # for loop to run through multiple iterations
  for(k in 1:iteration){
    
    # Create empty output dataframe with all relevant information
    res <- data.frame(OBJIDsic = character(),
                      region_id = character(),
                      code = character(),
                      name = character(),
                      island = character(),
                      IDregion = character(),
                      regionF = character(),
                      LAT = character(),
                      LON = character(),
                      GeodAREA = character(),
                      UN_geospat = character())
    
    
    # for loop to run through all UN geospatial regions to select random regions based on the observed number of mainland and island regions in that UN geospatial region
    for(i in 1:length(un)){
      
      # subset ALL regions to a specific UN geospatial region
      reg <- regions[which(regions$UN_geospat == un[i]),]
      # subset MAINLAND regions to a specific UN geospatial region
      reg.main <- reg[which(reg$island == "MAIN"),]
      # subset ISLAND regions to a specific UN geospatial region
      reg.isl <- reg[which(reg$island == "ISL"),]
      
      # Subset distance matrix for MAINLAND regions
      reg.dist.main <- m[which(rownames(m) %in% reg.main$OBJIDsic),which(colnames(m) %in% reg.main$OBJIDsic)]
      # Subset distance matrix for ISLAND regions
      reg.dist.isl <- m[which(rownames(m) %in% reg.isl$OBJIDsic),which(colnames(m) %in% reg.isl$OBJIDsic)]  # Subset distance matrix for island regions per UN geospatial region
      
      # Create list-object to store the random draws for mainland or island regions
      reg_ran <- list()
      
      
      # for loop to run calculations individually for mainland and island regions within a UN geospatial region
      for(j in 1:length(mi)){
        
        # if - CONDITION: if mi defines for ISLAND regions run the code below
        if(mi[j] == "ISL"){
          
          # if - CONDITION: if there is no ISLAND region in the observed empire, write an empty dataframe in the list
          if(n.emp.reg.geospat[un[i],mi[j]] == 0){
            
            d <- data.frame(OBJIDsic = character(),
                            region_id = character(),
                            code = character(),
                            name = character(),
                            island = character(),
                            IDregion = character(),
                            regionF = character(),
                            LAT = character(),
                            LON = character(),
                            GeodAREA = character(),
                            UN_geospat = character())
            
            reg_ran[[j]] <- d
            
            
          }
          
          else{
            
            # for loop to run through the observed number of mainland/island regions within a UN geospatial region
            for(p in 1:n.emp.reg.geospat[un[i],mi[j]]){
              
              # if - CONDITION: if this is the draw of the first region then draw one completely at random
              if(p == 1){
                
                # select first region as a random sample of island/mainland regions within the respective UN geospatial region
                reg_ran[[j]] <- sample_n(reg[which(reg$isl == mi[j]),], 1, replace = F)
                
              } 
              
              # else - CONDITION: if this is NOT the first draw then follow the procedure afterwards
              else {
                
                # remove region from first draw from distance matrix
                reg.dist.isl.2 <- reg.dist.isl[,-which(colnames(reg.dist.isl) %in% reg_ran[[j]]$OBJIDsic)]
                # OBJIDsic of last region colonized
                last.col <- tail(reg_ran[[j]],1)$OBJIDsic
                
                # if CONDITION: if all regions are colonized in the observed empire, for the last region the distance dataframe becomes a vector (the if then else statement accounts for that)
                if(is.vector(reg.dist.isl.2) == FALSE){
                  # find minimun distance for last region
                  min.dist <- min(reg.dist.isl.2[rownames(reg.dist.isl.2) == last.col,], na.rm = T)
                  # extract column number of the minimum distance value
                  new.col <- names(which(reg.dist.isl.2[rownames(reg.dist.isl.2) == last.col,] == min.dist, arr.ind = TRUE))
                }
                else{
                  # find minimun distance for last region
                  min.dist <- min(reg.dist.isl.2, na.rm = T)
                  # extract column number of the minimum distance value
                  new.col <- names(which(reg.dist.isl.2 == min.dist, arr.ind = TRUE))
                  
                }
                
                # extract record for the next (subsequently) colonized island
                next.reg <- reg[which(reg$OBJIDsic == new.col),]
                # add ne island to the new dataframe
                reg_ran[[j]] <- rbind(reg_ran[[j]], next.reg)
                
              }
            } 
          }
        }
        
        # else - CONDITION: if mi defines for MAINLAND regions run code below
        else{
          
          # if - CONDITION: if there is no MAINLAND region in the observed empire, write an empty dataframe in the list
          if(n.emp.reg.geospat[un[i],mi[j]] == 0){
            
            d <- data.frame(OBJIDsic = character(),
                            region_id = character(),
                            code = character(),
                            name = character(),
                            island = character(),
                            IDregion = character(),
                            regionF = character(),
                            LAT = character(),
                            LON = character(),
                            GeodAREA = character(),
                            UN_geospat = character())
            
            reg_ran[[j]] <- d
            
          }
          
          else{
            
            # for loop to run through the observed number of mainland/island regions within a UN geospatial region
            for(p in 1:n.emp.reg.geospat[un[i],mi[j]]){
              
              # if - CONDITION: if this is the draw of the first region then draw one completely at random
              if(p == 1){
                
                # select first region as a random sample of island/mainland regions within the respective UN geospatial region
                reg_ran[[j]] <- sample_n(reg[which(reg$isl == mi[j]),], 1, replace = F)
                
              } 
              
              # else - CONDITION: if this is NOT the first draw then follow the procedure afterwards
              else {
                
                # remove region from first draw from distance matrix
                reg.dist.main.2 <- reg.dist.main[,-which(colnames(reg.dist.main) %in% reg_ran[[j]]$OBJIDsic)]
                # OBJIDsic of last region colonized
                last.col <- tail(reg_ran[[j]],1)$OBJIDsic
                
                # if CONDITION: if all regions are colonized in the observed empire, for the last region the distance dataframe becomes a vector (the if then else statement accounts for that)
                if(is.vector(reg.dist.main.2) == FALSE){
                  # find minimun distance for last region
                  min.dist <- min(reg.dist.main.2[rownames(reg.dist.main.2) == last.col,], na.rm = T)
                  # extract column number of the minimum distance value
                  new.col <- names(which(reg.dist.main.2[rownames(reg.dist.main.2) == last.col,] == min.dist, arr.ind = TRUE))
                }
                else{
                  # find minimun distance for last region
                  min.dist <- min(reg.dist.main.2, na.rm = T)
                  # extract column number of the minimum distance value
                  new.col <- names(which(reg.dist.main.2 == min.dist, arr.ind = TRUE))
                  
                }
                
                
                # extract record for the next (subsequently) colonized island
                next.reg <- reg[which(reg$OBJIDsic == new.col),]
                # add ne island to the new dataframe
                reg_ran[[j]] <- rbind(reg_ran[[j]], next.reg)
                
              }
              
            }
          }
          
          
        }
        
        
      }
      
      # rbind random sample of islands and mainlands for the specific geospatial region
      dat <- do.call(rbind,reg_ran)
      # write results from UN geospatial region in final results dataframe
      res <- rbind(res, dat)
      
      result[[k]] <- res
    }
    
    
    
    
  }
  
  return(result)
  
}

