##############################################################################################
##
## null_model.r
## Constructing null models for Network-Area relationships of the European terrestrial food webs across bioregions
##
## This script implements two null models for the scaling of network properties of the European food web
## across bioregions.
##
## For this example, an subset of the entire dataset is used to illustrate the scaling
## procedure.
##
## Author: Dr Miguel Lurgi
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


##### To run this script we need to load two databases that hold the interaction matrix (i.e., food web)
##### and the geographical distribution of species in that network respectively

load('test-data.rda')

## First we load the necessary libraries
require(raster)
require(rgdal)
require(igraph)
require(cheddar)
require(AICcmodavg)

## and the spatial data for the European bioregions

europeRaster <- raster(x="./mask10k/reference_grid_10km.img")
cells_info <- read.dbf('./mask10k/reference_grid_10km.img.vat.dbf')

shape <- readOGR(dsn = "./BiogeoRegions2016_shapefile", layer = "BiogeoRegions2016")

region <- data.frame(PAGENAME = master$PAGENAME, SPP = 0)
regions_to_remove <- c('outside', 'macaronesia')

############################### NULL MODEL ###########################################
## One of the main criticisms from one of the reviewers of the paper for
## European bioregions was on the fact that we need to prove that the patterns
## observed are solely driven by species richness in a more robust way.
## For this, he / she suggested a null model based on choosing a set of species randomly
## from the metaweb based on species richness at each spatial scale and then build
## NARs from that. In the following code I implement such a null model.

## We complemented this with the implementation of another null model based on the generation of a random
## graph using the same number of species and links from the corresponding empirical
## network at each spatial scale.

## To take the number of species we use the outputs from the original runs
output_runs <- read.csv('./partial-output.csv', stringsAsFactors = FALSE)

replicates <- unique(output_runs$replicate)
output_null_model <- NULL
colors <- c('firebrick', 'navyblue', 'orange')

degree_dist_params <- NULL

#### and here we extract the cell values
for(i in as.numeric(as.character(shape@data$PK_UID))){
  if(shape@data[which(shape@data$PK_UID == i),]$short_name %in% regions_to_remove) next;
  
  cur_pol <- shape[which(shape@data$PK_UID == i),]
  region_name <- as.character(shape@data[which(shape@data$PK_UID == i),]$code)
  
  print(region_name)
  
  cur_ids <- extract(europeRaster, cur_pol)[[1]]
  
  cur_codes <- as.character(cells_info[which(cells_info$Value %in% cur_ids),]$PageName)
  
  ##### this part of the code was added on July 18th 2018 when we realised we needed to
  ##### obtain the metaweb per bioregion beforehand to incorporate the information of 
  ##### potential indegree of the species to the output results
  
  ##### here we obtain the metaweb for the bioregion using the information of the cells
  meta_comm <- master[which(master$PAGENAME %in% cur_codes),]
  
  if(length(which(is.na(rowSums(meta_comm[-1])))) > 0){
    meta_comm <- meta_comm[-which(is.na(rowSums(meta_comm[-1]))),]  
  }
  
  meta_comm <- colSums(meta_comm[-1])
  meta_comm[meta_comm > 1] <- 1
  # sum(meta_comm)
  # Defining the species pool at the bioregion level
  species_pool <- names(meta_comm)[which(meta_comm == 1)]
  species_pool <- setdiff(species_pool, absent_species)
  
  M <- BARMdiet.binary[species_pool,species_pool] # Web adjacency matrix (predator x prey)
  
  metaweb <- graph_from_adjacency_matrix(M)
  
  ######################### here ends the metaweb ##############################
  
  for(r in replicates){
    print(r)
    cur_out_runs <- output_runs[which(output_runs$dataset == region_name & output_runs$replicate == r),] 
    
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
    prev_sps <- 0
    
    for(a in cur_out_runs$areas){
      # print(a)
      cur_sps <- cur_out_runs[which(cur_out_runs$area == a),]$species
      cur_links <- cur_out_runs[which(cur_out_runs$area == a),]$links
      
      if(prev_sps != cur_sps){
        # Sampling the metaweb
        sppSTRING <- sample(V(metaweb)$name, cur_sps)
        
        ######### These are the lines for the null model
        #### Null model 1
        # cur_net <- induced_subgraph(metaweb, sppSTRING)
        
        #### Null model 2
        cur_net <- sample_gnm(cur_sps, cur_links, directed=TRUE)
        ##########################
        
#         if(a == 1) {
#           FW <- list()
#           FW$M <- get.adjacency(cur_net, sparse = FALSE)
#           FW <- TrophicLevels(FW)
#           FW$TrLevels <- FW$TrLevels-1
#           tr_layout <- GetRankLayout(FW)
#           
#           pdf(paste0('sample-net-', region_name,'.pdf'))
#           par(mar=c(0,0,1,0))
#           plot(cur_net, vertex.size=5, vertex.color=colors[FW$TrLevels+1], vertex.label=NA, edge.arrow.size=.3, layout=tr_layout, main=region_name, edge.width=0.5)
#           dev.off()
#         }
#         
#       }
#     }
#   }
# }
        M <- as_adjacency_matrix(cur_net, sparse = FALSE)
        S <- length(sppSTRING)
        ls <- sum(M)
        C <- (2*ls/((S-1)*S))
        l.s <- ls/S
        gen <- Generality(M)
        vul <- Vulnerability(M)
        f_basal <- FractionOfBasal(M)
        f_top <- FractionOfTop(M)
        f_int <- FractionOfIntermediate(M)
        sd_g <- SDGenerality(M)/l.s
        sd_v <- SDVulnerability(M)/l.s
        pred_overlap <- CalculatePredatorOverlap(M)
        s_top <- NumberOfTop(M)
        s_int <- NumberOfIntermediate(M)
        s_basal <- NumberOfBasal(M)

        modularity <- tryCatch({
          max(walktrap.community(cur_net)$modularity)
        }, warning = function(w) {
          NA
        }, error = function(e) {
          NA
        }, finally = {
        })

        # normal_id <- mean(igraph::degree(cur_net, V(cur_net), mode='in')[-which(igraph::degree(cur_net, V(cur_net), mode='in') == 0)]/l.s)
        # normal_od <- mean(igraph::degree(cur_net, V(cur_net), mode='out')[-which(igraph::degree(cur_net, V(cur_net), mode='out') == 0)]/l.s)
        # 
        cons <- length(which(igraph::degree(cur_net, V(cur_net), mode='in') != 0))
        resrcs <- length(which(igraph::degree(cur_net, V(cur_net), mode='out') != 0))
        
        # pot_indeg <- mean(igraph::degree(metaweb, names(which(igraph::degree(cur_net, mode='in') != 0)), mode='in'))
        
        prev_sps <- cur_sps
        
        
        degs <- igraph::degree(cur_net, mode='all')
        if(length(which(degs == 0)) != 0){
          degs <- degs[-which(degs == 0)]
        }
        
        ###################### this is for the degree distributions
        norm_term <- 1 #ecount(local_net)/vcount(local_net)
        
        occur = as.vector(table(degs/norm_term))
        occur = occur/sum(occur)
        p = occur/sum(occur)
        y = rev(cumsum(rev(p)))
        x = as.numeric(names(table(degs/norm_term)))
        temp <- data.frame(x,y)
        
        failed <- tryCatch({
          mod1 <- tryCatch({
            nls(y ~ ( (x^-a) *(exp(-x/b))), data = temp, start = list(a = .1, b = 20), control=nls.control(maxiter = 1e3))
          }, error = function(e) {
            NA
          }, finally = {
          })
          mod2 <- tryCatch({
            nls(y ~ (exp(-x/b)), data = temp, start = list(b = 35), control=nls.control(maxiter = 1000))
          }, error = function(e) {
            NA
          }, finally = {
          })
          mod3 <- tryCatch({
            nls(y ~ (x^-a), data = temp, start = list(a = .1), control=nls.control(maxiter = 1e3))
          }, error = function(e) {
            NA
          }, finally = {
          })
          mod4 <- tryCatch({
            nls(y ~ ( (1/ (x * b * sqrt(2*pi) )) * exp(- ( ((log(x) - a)^2) / (2*(b^2)) )) ), data = temp, start = list(a = 1.5, b = .3), control=nls.control(maxiter = 1000))
          }, error = function(e) {
            NA
          }, finally = {
          })
          FALSE
        }, warning = function(w) {
          FALSE
        }, error = function(e) {
          TRUE
        }, finally = {
          
        })
        
        if(! failed){
          model_list <- list(mod1,mod2,mod3,mod4)
          names_list <- c('mod1','mod2','mod3','mod4')
          if(length(which(is.na(model_list))) != 0){
            names_list <- names_list[-which(is.na(model_list))]
            model_list <- model_list[-which(is.na(model_list))]
          }
          names(model_list) <- names_list
          
          if(length(model_list) != 0){
            
            if(length(model_list) == 1){
              model_name <- names_list[1]
              model <- summary(eval(as.symbol(model_name)))
            }else{
              aic_comp <- aictab(model_list)
              model_name <- tolower(as.character(aic_comp$Modnames[1]))
              model <- summary(eval(as.symbol(model_name)))
            }
            
            if(model_name == 'mod1' | model_name == 'mod4'){
              cur_out <- data.frame(replicate=r+70, area=a, region=region_name,
                                    model=model_name, a=model$coefficients[1,1],
                                    a.std.err=model$coefficients[1,2], a.tval=model$coefficients[1,3],
                                    a.pval=model$coefficients[1,4], b=model$coefficients[2,1],
                                    b.std.err=model$coefficients[2,2], b.tval=model$coefficients[2,3],
                                    b.pval=model$coefficients[2,4])
            }else if(model_name == 'mod2'){
              cur_out <- data.frame(replicate=r+70, area=a, region=region_name,
                                    model=model_name, a=NA, a.std.err=NA, a.tval=NA, a.pval=NA,
                                    b=model$coefficients[1,1], b.std.err=model$coefficients[1,2],
                                    b.tval=model$coefficients[1,3], b.pval=model$coefficients[1,4])
            }else{
              cur_out <- data.frame(replicate=r+70, area=a, region=region_name,
                                    model=model_name, a=model$coefficients[1,1],
                                    a.std.err=model$coefficients[1,2], a.tval=model$coefficients[1,3],
                                    a.pval=model$coefficients[1,4], b=NA, b.std.err=NA, b.tval=NA,
                                    b.pval=NA)
            }
            
            if(is.null(degree_dist_params)){
              degree_dist_params <- cur_out
            }else{
              degree_dist_params <- rbind(degree_dist_params, cur_out)
            }
          }
        }
      }
      
      areas <- append(areas, a)
      species <- append(species, S)
      links <- append(links, ls)
      connectances <- append(connectances, C)
      links_per_sp <- append(links_per_sp, l.s)
      indegrees <- append(indegrees, gen)
      outdegrees <- append(outdegrees, vul)
      basals <- append(basals, f_basal)
      tops <- append(tops, f_top)
      intermediates <- append(intermediates, f_int)
      sd_gen <- append(sd_gen, sd_g)
      sd_vul <- append(sd_vul, sd_v)
      overlaps <- append(overlaps, pred_overlap)
      S_tops <- append(S_tops, s_top)
      S_intermediates <- append(S_intermediates, s_int)
      S_basals <- append(S_basals, s_basal)
      modularities <- append(modularities, modularity)
      # normalised_indegree <- append(normalised_indegree, normal_id)
      # normalised_outdegree <- append(normalised_outdegree, normal_od)
      
      consumers <- append(consumers, cons)
      resources <- append(resources, resrcs)
      
      # potential_indegree <- append(potential_indegree, pot_indeg)
      
    }
    
    cur_out <- data.frame(region=region_name, replicate=r, area=areas, species, links, links_per_sp, connectances,
                          indegrees, outdegrees, basals, tops, intermediates, sd_gen, sd_vul, overlaps,
                          S_tops, S_intermediates, S_basals, modularities,
                          resources, consumers)
    
    
    if(is.null(output_null_model)) output_null_model <- cur_out
    else output_null_model <- rbind(output_null_model, cur_out)
    
  }
}

write.csv(output_null_model, file='output-europe-null-model-2.csv')
write.csv(degree_dist_params, file='fits-degree-dists-europe-null-model-2.csv')



