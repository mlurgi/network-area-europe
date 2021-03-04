
##############################################################################################
## european-network-area-relationships.r
## Constructing Network-Area relationships for the European vertebrate food web
##
## This script implements different map cells aggregation algorithms to simulate
## increases in spatial scales across Europe and calculates food web properties
## on the food web determined by the species present in each given area. To do
## this it manipulates species distribution maps and a meta-network (metaweb) of 
## trophic interactions constructed at the European regional scale.
##
## For this example, an artificial dataset was created to illustrate the scaling
## procedure.
##
## Authors: Dr Miguel Lurgi
## Lecturer in Biosciences (Computational Ecology)
## Computational Ecology Lab - Department of Biosciences
## Swansea University, UK
## 
## and
##
## Centre for Biodiversity Theory and Modelling
## Theoretical and Experimental Ecology Station, CNRS, France
##
## Date Created: June-2016
##
## Copyright (c) Miguel Lurgi, 2016-2021
## Email: miguel.lurgi@swansea.ac.uk
##
## AND
##
## Dr Nuria Galiana
##
## Centre for Biodiversity Theory and Modelling
## Theoretical and Experimental Ecology Station, CNRS, France
##
## Notes:
##
## This script is provided as supplementary material for the paper:
##
## 'The spatial scaling of food web structure across European biogeographical regions'
## By Galiana N, Barros C, Braga J, Ficetola GF, Maiorano L, Thuiller W, Montoya JM & Lurgi M
## Ecography (2021) https://doi.org/10.1111/ecog.05229
##
## which contains additional information about the data manipulated in this script and their 
## sources (including other supplementary material)
##
##############################################################################################

################################ HELPER FUNCTION (By Joao Braga) ################################

###### This function is used throughout the code to perform different actions

# Function to transform spp distribution (it can be network properties per pixel) into raster
# The code for this function was provided by Joao Braga
fun.dbf2raster <- function(SPPPA, mask.dir = NULL){
  # SPPPA must be in this format - first colmun with CELL ID (PAGENAME) and second column with the value (to plot)
  if(is.null(mask.dir)) stop("Must specify the mask file directory")
  require(raster)
  require(foreign)
  
  maskID <- read.dbf(list.files(path = mask.dir, full.names = TRUE, pattern = ".img.vat.dbf$"))
  maskk <- raster(x = list.files(path = mask.dir, full.names = TRUE, pattern = ".img$"))
  
  spp <- maskID
  spp$val <- NA
  spp$PageName <- as.character(spp$PageName)
  row.names(spp) <- spp$PageName
  
  SPPPA$PAGENAME <- as.character(SPPPA$PAGENAME)
  SPPPA[,2] <- as.numeric(as.character(SPPPA[,2]))
  row.names(SPPPA) <- SPPPA$PAGENAME
  
  cellID <- as.character(SPPPA$PAGENAME)
  if( nrow(spp[cellID,]) != nrow(SPPPA[cellID,])) stop("Cell IDs do not match")
  spp <- spp[cellID,] 
  spp$val <- SPPPA[,2]
  
  xx <- raster::values(maskk)
  
  if( length(xx[!is.na(xx)]) != nrow(spp)) stop("Mask size inadequate")
  xx[!is.na(xx)] <- spp$val
  
  raster::values(maskk) <- xx
  return(maskk)
}

############################################################################################

#### Loading required libraries
require(raster)
require(rgdal)
require(igraph)
require(betapart)
require(ggplot2)

##### this is for building the network-area by bioregion

##### To run this script we need to load two databases that hold the interaction matrix (i.e., food web)
##### and the geographical distribution of species in that network respectively

load('test-data.rda')

##### we also need a bunch of utility functions to calculate food web properties which are implemented in
##### the source file 'utils.r'
source('utils.r')

######### even though not a function per se, this piece of code is necessary to build
######### data structures used by the spiral aggregation algorithm to go through the landscape
latitudes <- c()
longitudes <- c()

for(pn in master$PAGENAME){
  code_lat <- gsub("[0-9]", "", pn);
  code_lon <- gsub("[A-Z]", "", pn);
  
  if(! code_lat %in% latitudes)latitudes <- append(latitudes, code_lat);
  if(! code_lon %in% longitudes)longitudes <- append(longitudes, code_lon);
}

longitudes <- sort(as.numeric(longitudes))

#the extent of the raster
rows <- 589
cols <- 671


##### to obtain the mapping information we read the european grid and the maps of the bioregions obtained
##### from the European Environment Agency

europeRaster <- raster(x="./mask10k/reference_grid_10km.img")
cells_info <- read.dbf('./mask10k/reference_grid_10km.img.vat.dbf')

shape <- readOGR(dsn = "./BiogeoRegions2016_shapefile", layer = "BiogeoRegions2016")

region <- data.frame(PAGENAME = all_cells, SPP = 0)
regions_to_remove <- c('outside')


colours <- rainbow(length(shape@data$PK_UID))
idx <- 1
regions_properties <- NULL
for(i in shape@data$PK_UID){
  
  if(shape@data[which(shape@data$PK_UID == i),]$short_name %in% regions_to_remove) next;
  
  cur_pol <- shape[which(shape@data$PK_UID == i),]
  region_name <- as.character(shape@data[which(shape@data$PK_UID == i),]$code)
  cur_ids <- extract(europeRaster, cur_pol)[[1]]
  
  cur_codes <- as.character(cells_info[which(cells_info$Value %in% cur_ids),]$PageName)
  
  region$SPP <- 0
  region[which(region$PAGENAME %in% cur_codes),]$SPP <- 1
  
  # bioregion_raster <- fun.dbf2raster(SPPPA = region, mask.dir = "./mask10k/")
  # 
  # plot(bioregion_raster, main=region_name)
  
  cur_comm <- master[which(master$PAGENAME %in% cur_codes),]
  
  if(dim(cur_comm)[1] == 0) next
  
  cur_comm[is.na(cur_comm)] <- 0
  
  row.names(cur_comm) <- cur_comm$PAGENAME
  cur_comm <- cur_comm[,-which(names(cur_comm) == 'PAGENAME')]
  
  if(length(which(colSums(cur_comm) == 0)) > 0){
    cur_comm <- cur_comm[,-which(colSums(cur_comm) == 0)]  
  }
  
  beta_div <- beta.multi(cur_comm)$beta.SOR
  cur_comm <- colSums(cur_comm)
  cur_comm[cur_comm > 1] <- 1
  
  sum(cur_comm)
  # Defining the species pool
  sppSTRING <- names(cur_comm)[which(cur_comm == 1)]
  
  ##### Uncomment the following line if using the real (i.e., complete) dataset
  # sppSTRING <- setdiff(sppSTRING, absent_species)
  
  # Sampling the metaweb
  M <- BARMdiet.binary[sppSTRING,sppSTRING]
  diag(M) <- 0
  net_bio <- graph_from_adjacency_matrix(M)
  
  # pdf('network.pdf')
  # plot(net_bio, main=region_name, layout=layout.kamada.kawai, vertex.label=NA, vertex.size=.5, edge.arrow.size=.2, edge.width=.2)
  # dev.off()
  
  if(idx > 1) lines(degree.distribution(net_bio, cumulative=T), col=colours[idx])
  else plot(degree.distribution(net_bio, cumulative=T), type='l', col=colours[idx])
  
  idx <- idx + 1
  
  print(region_name)
  S <- vcount(net_bio)
  L <- ecount(net_bio)
  L.S <- ecount(net_bio)/vcount(net_bio)
  
  C <- sum(M)/((dim(M)[1])*(dim(M)[1]-1)) #Connectance if we prevent cannibalism
  
  mfcl <- tryCatch({
    MeanFoodChainLength(M)
  }, warning = function(w) {
    'NA'
  }, error = function(e) {
    'NA'
  }, finally = {
  })
  
  mfcl <- 0
  
  indegree <- Generality(M)
  
  outdegree <- Vulnerability(M)
  maxsim <- MaximumSimilarity(M)
  
  basal <- FractionOfBasal(M)
  top <- FractionOfTop(M)
  intermediate <- FractionOfIntermediate(M)
  sd_gen <- SDGenerality(M)
  sd_vul <- SDVulnerability(M)
  
  overlap <- CalculatePredatorOverlap(M)
  
  omnivory <- Omnivory(M)
  S_top <- NumberOfTop(M)
  S_intermediate <- NumberOfIntermediate(M)
  S_basal <- NumberOfBasal(M)
  
  cur_props <- data.frame(region_name, area=length(cur_codes), S, L, L.S, C, mfcl, indegree, outdegree, maxsim, basal, top, intermediate, sd_gen, sd_vul, overlap, omnivory, S_top, S_basal, S_intermediate, beta_div)
  
  if(is.null(regions_properties)) regions_properties <- cur_props
  else regions_properties <- rbind(regions_properties, cur_props)
  
}

regions_order <- c('Macaronesia', 'Mediterranean', 'Anatolian', 'BlackSea', 'Continental', 'Pannonian', 'Alpine', 'Atlantic', 'Steppic', 'Boreal', 'Arctic')

ggplot(regions_properties, aes(region_name, S)) + geom_point()


##### If several cpu's are available this code can run in parallel to make replicates at the same time
##### uncomment this code to enable paralellisation
# library(foreach)
# library(doParallel)
# cl <- makeCluster(15)     # this function receives as an argument the number of CPUs to use
# registerDoParallel(cl)
# start_time <- proc.time()

##### Here we declare the spatial aggregation method (see the paper's methods for specificaitons on how this work)
agg_types <- c('random', 'spiral', 'linear')
output_bioregions_types <- NULL

#### this is relevant for linear aggreagation only, whether increase area from south to north
south_north <- FALSE

for(agg_type in agg_types){
  print(agg_type)
  ##### Change here the number of replicates
  if(agg_type == 'linear') replicates <- 1
  else replicates <- 1      #change this to more replicates
  
  cur_dat <- NULL
  
  #### let's make it a bit fun and parallelise! 
  #### if parallel, uncomment the following line and comment the following one
  #cur_dat <- foreach(r = 1:replicates, .combine=rbind) %dopar%{
  for(r in 1:replicates){
    print(r)
    require(igraph)
    require(cheddar)
    require(raster)
    output <- NULL
    #### and here we extract the cell values
    for(i in shape@data$PK_UID){
      if(shape@data[which(shape@data$PK_UID == i),]$short_name %in% regions_to_remove) next;
      region_name <- as.character(shape@data[which(shape@data$PK_UID == i),]$code)
      
      #### comment the following line if you are using the complete dataset
      if(region_name != 'Steppic'){
        next
      }
      
      cur_pol <- shape[which(shape@data$PK_UID == i),]
      cur_ids <- extract(europeRaster, cur_pol)[[1]]
      
      cur_codes <- as.character(cells_info[which(cells_info$Value %in% cur_ids),]$PageName)
      
      # region$SPP <- 0
      # region[which(region$PAGENAME %in% cur_codes),]$SPP <- 1
      # 
      # bioregion_raster <- fun.dbf2raster(SPPPA = region, mask.dir = "./mask10k/")
      # plot(bioregion_raster, main=region_name)
      
      ##### this part of the code was added on July 18th 2018 when we realised we needed to
      ##### obtain the metaweb per bioregion beforehand to incorporate the information of 
      ##### potential indegree of the species to the output results
      
      ##### here we obtain the metaweb for the bioregion using the information of the cells
      meta_comm <- master[which(master$PAGENAME %in% cur_codes),]
      
      if(length(which(is.na(rowSums(meta_comm[,-which(names(cur_comm) == 'PAGENAME')])))) > 0){
        meta_comm <- meta_comm[-which(is.na(rowSums(meta_comm[,-which(names(cur_comm) == 'PAGENAME')]))),]  
      }
      
      meta_comm <- colSums(meta_comm[,-which(names(meta_comm) == 'PAGENAME')])
      meta_comm[meta_comm > 1] <- 1
      sum(meta_comm)
      # Defining the species pool at the bioregion level
      species_pool <- names(meta_comm)[which(meta_comm == 1)]
      
      ##### Uncomment the following line if using the real (i.e., complete) dataset
      # species_pool <- setdiff(species_pool, absent_species)
      
      M <- BARMdiet.binary[species_pool,species_pool] # Web adjacency matrix (prey x predator)
      metaweb <- graph_from_adjacency_matrix(M)
      ######################### here ends the metaweb ##############################
      
      
      if(agg_type == 'linear' & south_north){
        cur_codes <- rev(cur_codes)  
      }
      
      if(agg_type == 'random') cur_codes <- sample(cur_codes)
      if(agg_type == 'random' || agg_type == 'linear') j <- 1
      if(agg_type == 'spiral'){
        ### choose a random cell from the bioregion
        cur_cellid <- sample(cur_codes, 1)
        
        y_coord <- match(gsub("[0-9]", "", cur_cellid), latitudes);
        x_coord <- match(gsub("[A-Z]", "", cur_cellid), longitudes);
        
        x_index <- x_coord - 1;
        y_index <- y_coord - 1;
        x_count <- 1
        y_count <- 0
        x_next <- 1
        y_next <- 1
        x_offset <- -1
        y_offset <- 1
      }
      
      # region$SPP <- 0
      # region[which(region$PAGENAME %in% cur_codes),]$SPP <- 1
      # 
      # bioregion_raster <- fun.dbf2raster(SPPPA = region, mask.dir = "./mask10k/")
      # plot(bioregion_raster, main=region_name)
      
      cells_to_traverse <- length(cur_codes) #rows*cols
      traversed_cells <- 1
      visited_cells <- c()
      traversed_set <- c()
      
      ###### These are vectors for storing the values of the different network properties analysed
      species <- c()
      connectances <- c()
      links <- c()
      links_per_sp <- c()
      indegrees <- c()
      outdegrees <- c()
      maxsims <- c()
      basals <- c()
      tops <- c()
      intermediates <- c()
      sd_gen <- c()
      sd_vul <- c()
      overlaps <- c()
      omnivorys <- c()
      S_tops <- c()
      S_intermediates <- c()
      S_basals <- c()
      modularities <- c()
      normalised_indegree <- c()
      normalised_outdegree <- c()
      
      consumers <- c()
      resources <- c()
      
      potential_indegree <- c()
      
      areas <- c()
      
      ##### current region considered and number of cells to visit
      print(region_name)
      print(paste0('cells to traverse = ', cells_to_traverse))
      
      sppSTRING_prev <- NULL
      
      #### colors and cc are only used if degree distributions are plotted
      colors <- rev(rainbow(6, alpha=0.7))
      cc <- 1
      while(length(visited_cells) < cells_to_traverse){
        # if(traversed_cells %% 1000 == 0){
        #   print(traversed_cells)
        #   region$SPP <- 0
        #   region[which(region$PAGENAME %in% visited_cells),]$SPP <- 1
        # 
        #   bioregion_raster <- fun.dbf2raster(SPPPA = region, mask.dir = "./mask10k/")
        #   plot(bioregion_raster, main=region_name)
        # }
        if(agg_type == 'random' || agg_type == 'linear'){
          cur_cellid <- cur_codes[j]
          j <- j + 1
        }else if(agg_type == 'spiral') {
          cur_cellid <- paste0(latitudes[y_coord],longitudes[x_coord])
        }
        
        #print(cur_cellid)
        
        traversed_set <- append(traversed_set, cur_cellid)
        
        skip <- FALSE
        
        if((! cur_cellid %in% cur_codes) | (cur_cellid %in% visited_cells)) skip <- TRUE
        if(! skip){
          if(is.na(sum(master[which(master$PAGENAME == cur_cellid),-which(names(master) == 'PAGENAME')]))){
            skip <- TRUE
            cells_to_traverse <- cells_to_traverse - 1
          }
        }
        
        if(! skip){
          visited_cells <- append(visited_cells, cur_cellid)
          
          if(TRUE){
            cur_comm <- master[which(master$PAGENAME %in% visited_cells),]
            
            cur_comm <- colSums(cur_comm[-which(names(cur_comm) == 'PAGENAME')])
            cur_comm[cur_comm > 1] <- 1
            print(sum(cur_comm))
            # Defining the species pool
            sppSTRING <- names(cur_comm)[which(cur_comm == 1)]
            
            ##### UNCOMMENT IF COMPLETE DATASET IS USED
            # sppSTRING <- setdiff(sppSTRING, absent_species)
            
            if(length(sppSTRING) <= 1){
              S <- length(sppSTRING)
              ls <- NA
              C <- NA
              l.s <- NA
              gen <- NA
              vul <- NA
              max_sim <- NA
              f_basal <- NA
              f_top <- NA
              f_int <- NA
              sd_g <- NA
              sd_v <- NA
              pred_overlap <- NA
              omni <- NA
              s_top <- NA
              s_int <- NA
              s_basal <- NA
              modularity <- NA
              normal_id <- NA
              normal_od <- NA
              
              cons <- NA
              resrcs <- NA
              
              pot_indeg <- NA
              
              sppSTRING_prev <- sppSTRING
            }else if(is.null(sppSTRING_prev) | length(setdiff(sppSTRING, sppSTRING_prev)) != 0){
              # Sampling the metaweb
              M <- BARMdiet.binary[sppSTRING,sppSTRING] # Web adjacency matrix (predator x prey)
              S <- length(sppSTRING)
              ls <- sum(M)
              C <- (2*ls/((S-1)*S))
              l.s <- ls/S
              gen <- Generality(M)
              vul <- Vulnerability(M)
              max_sim <- 0 #MaximumSimilarity(M)
              f_basal <- FractionOfBasal(M)
              f_top <- FractionOfTop(M)
              f_int <- FractionOfIntermediate(M)
              sd_g <- SDGenerality(M)/l.s
              sd_v <- SDVulnerability(M)/l.s
              pred_overlap <- CalculatePredatorOverlap(M)
              omni <- 0 #Omnivory(M)
              s_top <- NumberOfTop(M)
              s_int <- NumberOfIntermediate(M)
              s_basal <- NumberOfBasal(M)
              
              
              cur_net <- graph_from_adjacency_matrix(M)
              modularity <- tryCatch({
                max(walktrap.community(cur_net)$modularity)
              }, warning = function(w) {
                NA
              }, error = function(e) {
                NA
              }, finally = {
              })
              
              normal_id <- mean(degree(cur_net, V(cur_net), mode='in')[-which(degree(cur_net, V(cur_net), mode='in') == 0)]/l.s)
              normal_od <- mean(degree(cur_net, V(cur_net), mode='out')[-which(degree(cur_net, V(cur_net), mode='out') == 0)]/l.s)
              
              cons <- length(which(degree(cur_net, V(cur_net), mode='in') != 0))
              resrcs <- length(which(degree(cur_net, V(cur_net), mode='out') != 0))
              
              pot_indeg <- mean(degree(metaweb, names(which(degree(cur_net, mode='in') != 0)), mode='in'))
              
              sppSTRING_prev <- sppSTRING
            }
            
            if(is.null(sppSTRING_prev))
              sppSTRING_prev <- sppSTRING
            
            species <- append(species, S)
            links <- append(links, ls)
            connectances <- append(connectances, C)
            links_per_sp <- append(links_per_sp, l.s)
            indegrees <- append(indegrees, gen)
            outdegrees <- append(outdegrees, vul)
            maxsims <- append(maxsims, max_sim)
            basals <- append(basals, f_basal)
            tops <- append(tops, f_top)
            intermediates <- append(intermediates, f_int)
            sd_gen <- append(sd_gen, sd_g)
            sd_vul <- append(sd_vul, sd_v)
            overlaps <- append(overlaps, pred_overlap)
            omnivorys <- append(omnivorys, omni)
            S_tops <- append(S_tops, s_top)
            S_intermediates <- append(S_intermediates, s_int)
            S_basals <- append(S_basals, s_basal)
            modularities <- append(modularities, modularity)
            normalised_indegree <- append(normalised_indegree, normal_id)
            normalised_outdegree <- append(normalised_outdegree, normal_od)
            
            consumers <- append(consumers, cons)
            resources <- append(resources, resrcs)
            
            potential_indegree <- append(potential_indegree, pot_indeg)
            
            areas <- append(areas, j)
            
            
            local_net <- graph_from_adjacency_matrix(M)
            
            degs <- igraph::degree(local_net, mode='all')
            
            ############## TO DO DEGREE DISTRIBUTIONS UNCOMMENT THE FOLLOWING CODE ########################
            norm_term <- ecount(local_net)/vcount(local_net)
            
            occur = as.vector(table(degs/norm_term))
            occur = occur/sum(occur)
            p = occur/sum(occur)
            y = rev(cumsum(rev(p)))
            x = as.numeric(names(table(degs/norm_term)))
            
            if(j == 1){
              plot(x,y, log='xy', ylim=c(10^-3,1), xlim=c(10^-2, 10^1), pch=16, cex=1.5, cex.lab=1.5, cex.axis=1.2, xlab='', ylab='log Pc(k)', col=colors[cc], tck=0, xaxt='n', yaxt='n')
              title(xlab='log k', line=2, cex.lab=1.5)
              
              box(lwd=1.5)
              x1 <- floor(log10(range(x)))
              pow <- seq(x1[1], x1[2])
              
              ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
              axis(1, 10^pow, tck=0.02, lwd=1.5, lwd.ticks=1.5, cex.axis=1.5, labels=c(expression(10^{-2}), expression(10^{-1}), expression(10^{0}), expression(10^{1})))
              # axis(1, 10^pow, tck=0.02, lwd=1.5, lwd.ticks=1.5, labels=c(expression(10^{0}), expression(10^{1}), expression(10^{2})))
              axis(1, ticksat, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1, tck=0.01)
              
              
              y1 <- floor(log10(range(y)))
              pow <- seq(y1[1], y1[2]+1)
              ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
              axis(2, 10^pow, tck=0.02, lwd=1.5, lwd.ticks=1.5, cex.axis=1.5, labels=c(expression(10^{-3}), expression(10^{-2}), expression(10^{-1}), expression(10^{0}), expression(10^{1})), las=1, mgp=c(1,.5,0))
              axis(2, ticksat, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1, tck=0.01)
              
            }else{
              points(x,y,pch=16, cex.lab=1.5, cex.axis=1.5, col=colors[cc], cex=1.5)
            }
            
            cc <- cc + 1
            ############## HERE END THE DEGREE DISTRIBUTIONS ########################
            
          }
        }
        
        if(agg_type == 'spiral'){
          if(x_count > 0){
            x_index <- (x_index+x_offset) %% cols
            x_coord <- x_index + 1;
            x_count <- x_count - 1
            
            if(x_count == 0){
              x_next <- x_next + 1
              if(x_offset == 1) x_offset <- -1
              else if(x_offset == -1) x_offset <- 1
              y_count <- y_next
              
            }
            
          }else if(y_count > 0){
            y_index <- (y_index + y_offset) %% rows
            y_coord <- y_index + 1;
            y_count <- y_count - 1
            
            if(y_count == 0){
              y_next <- y_next + 1
              if(y_offset == 1) y_offset <- -1
              else if(y_offset == -1) y_offset <- 1
              x_count <- x_next
            }
          }
        }
        
        
        traversed_cells <- traversed_cells + 1
        
        ###### this if is for plotting the current progress in map exploration
        if((length(visited_cells) %% 1000) == 0){
          region$SPP <- 0
          region[which(region$PAGENAME %in% traversed_set),]$SPP <- 1
          
          bioregion_raster <- fun.dbf2raster(SPPPA = region, mask.dir = "./mask10k/")
          
          plot(bioregion_raster, main=region_name)  
        }
        
      }
      
      # par(mar=c(4,4,2,2))
      # plot(links, xlab='area', main=region_name)
      
      cur_out <- data.frame(region=rep(region_name, length(species)), area=1:length(species), species, links, links_per_sp, connectances,
                            indegrees, outdegrees, maxsims, basals, tops, intermediates, sd_gen, sd_vul, overlaps,
                            omnivorys, S_tops, S_intermediates, S_basals, modularities, normalised_indegree, normalised_outdegree, 
                            resources, consumers, potential_indegree)
      
      if(is.null(output)) output <- cur_out
      else output <- rbind(output, cur_out)
    }
    
    ###### UNCOMMENT THIS IF PARALLEL ######
    #cbind(data.frame(replicate=rep(r,dim(output)[1])),output)
    
    
    ####### COMMENT THIS IF PARALLEL #########
    if(is.null(cur_dat))
      cur_dat <- cbind(data.frame(replicate=rep(r,dim(output)[1])),output)
    else
      cur_dat <- rbind(cur_dat, cbind(data.frame(replicate=rep(r,dim(output)[1])),output))
    ##########################################
  
  }
  
  cur_dat <- cbind(data.frame(agg_type = rep(agg_type,dim(cur_dat)[1])), cur_dat)
  
  if(is.null(output_bioregions_types)) output_bioregions_types <- cur_dat
  else output_bioregions_types <- rbind(output_bioregions_types, cur_dat)
}


# stop_time <- proc.time()
# print(stop_time - start_time)
# 
# stopCluster(cl)

###### To save the output into a table file
# write.csv(output_bioregions_types, file = 'test-output.csv')

require(ggplot2)

###### To plot some results to check out --- change the response variable for any other in output_bioregions_types
ggplot(subset(output_bioregions_types, replicate==1), aes(area, connectances, colour=region, shape=agg_type)) +
  geom_point()


