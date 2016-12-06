# Libraries
library(raster)

#################################################
# Metaweb with diet categories
load(file = "BARMdiet_bin.RData")

#### This is the metaweb ######

#### Row and column names are species codes according to the species names dictionary
#### which is kept in 'SppID.txt'

head(BARMdiet.binary)   # Format: rows are predators ;  columns are preys 

#################################################
# Species codes and species names
# Function to Identify a spp by the code
whois <- function(SPPCODE = NULL, SPPNAME = NULL) {
  # Function to Identify a spp by the code
  
  if(is.null(SPPCODE) & is.null(SPPNAME)) stop("Must specify a species code or name(Genus_species)")
  if(!is.null(SPPCODE) & !is.null(SPPNAME)) stop("Must specify a species code or name(Genus_species)")
  
  SppID <- read.table(file = "SppID.txt", header = TRUE, stringsAsFactors = F)
  
  if(length(SPPCODE) > 1){
    SPPCODE <- paste0(SPPCODE, "$", collapse = "|")
  }
  if(length(SPPNAME) > 1){
    SPPNAME <- paste0(SPPNAME, "$", collapse = "|")
  }
  
  if(!is.null(SPPCODE))    who <- SppID[which(SppID$ID==SPPCODE),]$SPPname
  if(!is.null(SPPNAME))    who <- SppID[grepl(pattern = SPPNAME,x = SppID$SPPname),]
  
  return(who)
}



# Function to transform spp distribution (it can be network properties per pixel) into raster
fun.dbf2raster <- function(SPPPA, mask.dir = NULL){
  # SPPPA must be in this format - first colmun with CELL ID (PAGENALE) and second column with the value (to plot)
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
  
  xx <- values(maskk)
  
  if( length(xx[!is.na(xx)]) != nrow(spp)) stop("Mask size inadequate")
  xx[!is.na(xx)] <- spp$val    #[xx[!is.na(xx)]]
  
  values(maskk) <- xx
  return(maskk)
}

whois(SPPCODE = "B122")
#      ID          SPPname
# 27 B122 Falco_rusticolus

whois(SPPNAME = "Falco")

#################################################
# Species presence/absence per pixel 
# e.g. at @10km
load(file = "MASTER.bin10000_allhab_tresh0.Rdata")
head(master) # rows are individual pixels; 1st columm (PAGENAME) is the unique name of the pixel, other columns are species

# Convert a dbf table to raster
# E.g. getting the spp "B122" distribution at @10k
sppDIST <- data.frame(PAGENAME = master$PAGENAME, SPP = master$B258) # 1st Column: PAGENAME, necessary to match the pixels; 2nd can be any value, in here is presence and absence, but it can be network properties
sppDIST$PAGENAME <- as.character(sppDIST$PAGENAME) # important: PAGENAME must be as character for this to work correctly
str(sppDIST)


sppRaster <- fun.dbf2raster(SPPPA = sppDIST, mask.dir = "./mask10k/")
sppRaster

plot(sppRaster, main= whois(SPPCODE = "B258"))

#################################################
# Local webs (pixel level): sample of the metaweb for a given pixel species pool
pixel <- master[52000,]

head(pixel)

sum(pixel[-1])              # Note: first value is the pixel name
# 241 species in this pixel 

# Defining the species pool
sppSTRING <- names(pixel)[which( pixel==1, arr.ind=FALSE)]

head(sppSTRING)

length(sppSTRING)
# 241

# Sampling the metaweb
localWEB <- BARMdiet.binary[sppSTRING,sppSTRING] # Web adjacency matrix (predator x prey)
dim(localWEB)
# 241 241

# end



####### THIS IS WHERE THE NETWORK-AREA CODE BEGINS
####### This code was added to Joao's on the 12/10/2016


#######   this is a test to see what the letters in the pagename mean

letters_in_codes <- c()

for(pn in sppDIST$PAGENAME){
  code <- gsub("[A-Z]", "", pn);
  
  if(! code %in% letters_in_codes){
    letters_in_codes <- append(letters_in_codes, code);
  }
}

letters_in_codes <- as.character(sort(as.numeric(letters_in_codes)))

for(pn in 1:(dim(sppDIST)[1])){
  pname <- sppDIST[pn,]$PAGENAME
  
  code <- gsub("[A-Z]", "", pname);
  
  sppDIST[pn,]$SPP <- which(letters_in_codes == code)
}

sppRaster <- fun.dbf2raster(SPPPA = sppDIST, mask.dir = "./mask10k/")
sppRaster

plot(sppRaster)


latitudes <- c()
longitudes <- c()

for(pn in master$PAGENAME){
  code_lat <- gsub("[0-9]", "", pn);
  code_lon <- gsub("[A-Z]", "", pn);
  
  if(! code_lat %in% latitudes)latitudes <- append(latitudes, code_lat);
  if(! code_lon %in% longitudes)longitudes <- append(longitudes, code_lon);
}

longitudes <- sort(as.numeric(longitudes))


###### the actual code for the runs starts here!!!!

#the extent of the raster
rows <- 589
cols <- 671

longitudes <- 1:cols
latitudes <- LETTERS

cur_lat <- length(latitudes)
done <- FALSE

for(f in LETTERS){
  for(s in LETTERS){
    cur_name <- paste(f,s,sep='')
    latitudes <- append(latitudes, cur_name)
    
    cur_lat <- cur_lat+1
    print(cur_lat)
    if(cur_lat == rows) done <- TRUE;
    
    if(done) break;
  }
  if(done) break;
}

absent_species <- c("B511", "B512", "R250")

x_coord <- sample(1:cols, 1)
y_coord <- sample(1:rows, 1)
cur_cellid <- paste0(latitudes[y_coord],longitudes[x_coord])

while(! cur_cellid %in% master$PAGENAME){
  x_coord <- sample(1:cols, 1)
  y_coord <- sample(1:rows, 1)
  cur_cellid <- paste0(latitudes[y_coord],longitudes[x_coord])
}

x_index <- x_coord - 1
y_index <- y_coord - 1

x_count <- 1
y_count <- 0

x_next <- 1
y_next <- 1

x_offset <- -1
y_offset <- 1

cells_to_traverse <- 10000 #rows*cols
traversed_cells <- 1

visited_cells <- c()

species <- c()
connectances <- c()
links <- c()
links_per_sp <- c()


total_cells <- length(master$PAGENAME)

while(length(visited_cells) < total_cells){
  print(traversed_cells)
 # print(paste('x = ',x_coord, ' --- y = ', y_coord,sep=''))
  
  cur_cellid <- paste0(latitudes[y_coord],longitudes[x_coord])
  #print(cur_cellid)
  
  skip <- FALSE
  
  if(! cur_cellid %in% master$PAGENAME) skip <- TRUE
  if(! skip){
    if(is.na(sum(master[which(master$PAGENAME == cur_cellid),-1]))){
      skip <- TRUE
      total_cells <- total_cells - 1
    }
  }
  
  if(! skip){
    visited_cells <- append(visited_cells, cur_cellid)
    
    cur_comm <- master[which(master$PAGENAME %in% visited_cells),]
     
    cur_comm <- colSums(cur_comm[-1])
    cur_comm[cur_comm > 1] <- 1
    sum(cur_comm)
    # Defining the species pool
    sppSTRING <- names(cur_comm)[which(cur_comm == 1)]
    
    sppSTRING <- setdiff(sppSTRING, absent_species)
    
    # Sampling the metaweb
    cur_web <- BARMdiet.binary[sppSTRING,sppSTRING] # Web adjacency matrix (predator x prey)
    
    S <- length(sppSTRING)
    species <- append(species, S)
    ls <- sum(cur_web)
    links <- append(links, ls)
    connectances <- append(connectances, (2*ls/((S-1)*S)))
    links_per_sp <- append(links_per_sp, ls/S)
  }
  
  if(x_count > 0){
    x_index <- (x_index+x_offset)%%cols
    x_coord <- x_index + 1
    x_count <- x_count - 1
    
    if(x_count == 0){
      x_next <- x_next + 1
      if(x_offset == 1) x_offset = -1
      else if(x_offset == -1) x_offset = 1
      
      y_count <- y_next
      
    }
    
  }else if(y_count > 0){
    y_index <- (y_index + y_offset) %% rows
    y_coord <- y_index + 1
    y_count <- y_count - 1
    
    if(y_count == 0){
      y_next <- y_next + 1
      if(y_offset == 1) y_offset <- -1
      else if(y_offset == -1) y_offset <- 1
      x_count <- x_next
    }
  }
  
  traversed_cells <- traversed_cells + 1
  
  
  if(length(visited_cells) %% 1000 == 0){
    print(paste0('I have visited ',length(visited_cells), ' cells'))
    # sppDIST$SPP <- 0
    # to_draw <- visited_cells[1:42000]
    # sppDIST[which(sppDIST$PAGENAME %in% to_draw),]$SPP <- 1
    # 
    # sppRaster <- fun.dbf2raster(SPPPA = sppDIST, mask.dir = "./mask10k/")
    # plot(sppRaster, main=42000)
  }
  
}

#### draw selected cells on the map
sppDIST$SPP <- 0

sppDIST[which(sppDIST$PAGENAME %in% visited_cells),]$SPP <- 1

sppRaster <- fun.dbf2raster(SPPPA = sppDIST, mask.dir = "./mask10k/")
sppRaster

plot(sppRaster)

par(mar=c(4,4,2,2))
plot(links, xlab='area', pch=20, cex=.3, col='red', log='xy')





################# this is for building the network-area by bioregion

require(raster)
require(rgdal)

europeRaster <- raster(x="./mask10k/reference_grid_10km.img")
cells_info <- read.dbf('./mask10k/reference_grid_10km.img.vat.dbf')

shape <- readOGR(dsn = "./BiogeoRegions2016_shapefile", layer = "BiogeoRegions2016")

region <- data.frame(PAGENAME = master$PAGENAME, SPP = 0)
regions_to_remove <- c('outside')

output <- NULL

#### and here we extract the cell values
for(i in 5:5){
  
  if(shape@data[which(shape@data$PK_UID == i),]$short_name %in% regions_to_remove) next;
  
  cur_pol <- shape[i,]
  region_name <- as.character(shape@data[which(shape@data$PK_UID == i),]$code)
  cur_ids <- extract(europeRaster, cur_pol)[[1]]
  
  
  cur_codes <- as.character(cells_info[which(cells_info$Value %in% cur_ids),]$PageName)
  
  region$SPP <- 0
  region[which(region$PAGENAME %in% cur_codes),]$SPP <- 1
  
  bioregion_raster <- fun.dbf2raster(SPPPA = region, mask.dir = "./mask10k/")
  
  plot(bioregion_raster, main=region_name)
  
  ### choose a random cell from the bioregion
  cur_cellid <- sample(cur_codes, 1)

  y_coord <- match(gsub("[0-9]", "", cur_cellid), latitudes);
  x_coord <- match(gsub("[A-Z]", "", cur_cellid), longitudes);
  
  # y_coord <- floor(rows/2)
  # x_coord <- floor(cols/2)
  # 
  x_index <- x_coord - 1;
  y_index <- y_coord - 1;
  
  x_count <- 1
  y_count <- 0
  
  x_next <- 1
  y_next <- 1
  
  x_offset <- -1
  y_offset <- 1
  
  cells_to_traverse <- length(cur_codes) #rows*cols
  traversed_cells <- 1
  
  visited_cells <- c()
  
  traversed_set <- c()
  
  species <- c()
  connectances <- c()
  links <- c()
  links_per_sp <- c()
  
  
  while(length(visited_cells) < cells_to_traverse){
    print(traversed_cells)
    # print(paste('x = ',x_coord, ' --- y = ', y_coord,sep=''))
    
    cur_cellid <- paste0(latitudes[y_coord],longitudes[x_coord])
    #print(cur_cellid)
    
    traversed_set <- append(traversed_set, cur_cellid)
    
    skip <- FALSE
    
    if((! cur_cellid %in% cur_codes) | (cur_cellid %in% visited_cells)) skip <- TRUE
    if(! skip){
      if(is.na(sum(master[which(master$PAGENAME == cur_cellid),-1]))){
        skip <- TRUE
        cells_to_traverse <- cells_to_traverse - 1
      }
    }
    
    if(! skip){
      visited_cells <- append(visited_cells, cur_cellid)
      
      cur_comm <- master[which(master$PAGENAME %in% visited_cells),]
      
      cur_comm <- colSums(cur_comm[-1])
      cur_comm[cur_comm > 1] <- 1
      sum(cur_comm)
      # Defining the species pool
      sppSTRING <- names(cur_comm)[which(cur_comm == 1)]
      
      sppSTRING <- setdiff(sppSTRING, absent_species)
      
      # Sampling the metaweb
      cur_web <- BARMdiet.binary[sppSTRING,sppSTRING] # Web adjacency matrix (predator x prey)
      
      S <- length(sppSTRING)
      species <- append(species, S)
      ls <- sum(cur_web)
      links <- append(links, ls)
      connectances <- append(connectances, (2*ls/((S-1)*S)))
      links_per_sp <- append(links_per_sp, ls/S)
    }
    
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
    traversed_cells <- traversed_cells + 1
    
  }
  
  par(mar=c(4,4,2,2))
  plot(links, xlab='area', main=region_name)
  
  cur_out <- data.frame(region=rep(region_name, length(species)), area=1:length(species), species, links, links_per_sp, connectances)
  
  if(is.null(output)) output <- cur_out
  else output <- rbind(output, cur_out)
  
  region$SPP <- 0
  region[which(region$PAGENAME %in% traversed_set),]$SPP <- 1
  
  bioregion_raster <- fun.dbf2raster(SPPPA = region, mask.dir = "./mask10k/")
  
  plot(bioregion_raster, main=region_name)
  
}


####### here we try a different way of aggregating cells, first within each bio-region and
####### then put together the regions
output_region_aggregation <- NULL
species <- c()
connectances <- c()
links <- c()
links_per_sp <- c()

visited_cells <- c()

for(i in unique(shape@data$PK_UID)){
  
  if(shape@data[which(shape@data$PK_UID == i),]$short_name %in% regions_to_remove) next;
  
  cur_pol <- shape[i,]
  region_name <- as.character(shape@data[which(shape@data$PK_UID == i),]$code)
  cur_ids <- extract(europeRaster, cur_pol)[[1]]
  
  cur_codes <- as.character(cells_info[which(cells_info$Value %in% cur_ids),]$PageName)
  
  region$SPP <- 0
  region[which(region$PAGENAME %in% cur_codes),]$SPP <- 1
  
  bioregion_raster <- fun.dbf2raster(SPPPA = region, mask.dir = "./mask10k/")
  
  plot(bioregion_raster, main=region_name)
  
  ### choose a random cell from the bioregion
  cur_cellid <- sample(cur_codes, 1)
  
  y_coord <- match(gsub("[0-9]", "", cur_cellid), latitudes);
  x_coord <- match(gsub("[A-Z]", "", cur_cellid), longitudes);

  cells_to_traverse <- length(cur_codes) #rows*cols
  traversed_cells <- 1
  
  
  traversed_set <- c()
  
  while(length(visited_cells) < cells_to_traverse){
    print(traversed_cells)
    cur_cellid <- paste0(latitudes[y_coord],longitudes[x_coord])
    traversed_set <- append(traversed_set, cur_cellid)
    
    skip <- FALSE
    
    if((! cur_cellid %in% cur_codes) | (cur_cellid %in% visited_cells)) skip <- TRUE
    if(! skip){
      if(is.na(sum(master[which(master$PAGENAME == cur_cellid),-1]))){
        skip <- TRUE
        cells_to_traverse <- cells_to_traverse - 1
      }
    }
    
    if(! skip){
      visited_cells <- append(visited_cells, cur_cellid)
      
      cur_comm <- master[which(master$PAGENAME %in% visited_cells),]
      
      cur_comm <- colSums(cur_comm[-1])
      cur_comm[cur_comm > 1] <- 1
      sum(cur_comm)
      # Defining the species pool
      sppSTRING <- names(cur_comm)[which(cur_comm == 1)]
      
      sppSTRING <- setdiff(sppSTRING, absent_species)
      
      # Sampling the metaweb
      cur_web <- BARMdiet.binary[sppSTRING,sppSTRING] # Web adjacency matrix (predator x prey)
      
      S <- length(sppSTRING)
      species <- append(species, S)
      ls <- sum(cur_web)
      links <- append(links, ls)
      connectances <- append(connectances, (2*ls/((S-1)*S)))
      links_per_sp <- append(links_per_sp, ls/S)
    }
    
     
  }
  
  cur_out <- data.frame(region=rep(region_name, length(species)), area=1:length(species), species, links, links_per_sp, connectances)
  
  if(is.null(output)) output <- cur_out
  else output <- rbind(output, cur_out)
  
}




###### this bit of code is for generating the plots
###### first maps and then properties

for(r in unique(output$region)){
  temp <- output[which(output$region == r),]
  for(p in names(temp)){
    if(p == 'region'){
      region <- data.frame(PAGENAME = master$PAGENAME, SPP = 0)
      cur_poly <- shape[which(shape$short_name == tolower(r)),]
      cur_ids <- extract(europeRaster, cur_poly)[[1]]
      cur_codes <- as.character(cells_info[which(cells_info$Value %in% cur_ids),]$PageName)
      
      region$SPP <- 0
      region[which(region$PAGENAME %in% cur_codes),]$SPP <- 1
      
      bioregion_raster <- fun.dbf2raster(SPPPA = region, mask.dir = "./mask10k/")
      
      pdf(file=paste0(r,'_map.pdf'), width = 7, height = 5)
      par(mar=c(4,4,2,2))
      plot(bioregion_raster, main=r)
      dev.off()
    }else if((p != 'X') & (p != 'area')){
      pdf(file=paste0(r,'_',p,'.pdf'), width = 7, height = 5)
      par(mar=c(4,4,2,2))
      plot(temp[,which(names(temp) == p)], xlab='area', ylab=p, main='')
      dev.off()
    }
    
  }
  
}

output_all_europe <- read.csv('output_all_europe.csv')
output_all_europe <- output_all_europe[,-1]

region <- data.frame(PAGENAME = master$PAGENAME, SPP = 1)
region$SPP[1] <- 0
all_europe_raster <- fun.dbf2raster(SPPPA = region, mask.dir = "./mask10k/")

pdf(file='europe_map.pdf', width = 7, height = 5)
par(mar=c(4,4,2,2))
plot(all_europe_raster, main='Europe')
dev.off()

for(p in names(output_all_europe)){
  pdf(file=paste0('all-europe_',p,'.pdf'), width = 7, height = 5)
  par(mar=c(4,4,2,2))
  plot(output_all_europe[,which(names(output_all_europe) == p)], xlab='area', ylab=p, main='')
  dev.off()
}


#### this is for testing the maps!


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

cells_to_traverse <- length(cur_codes) #rows*cols
traversed_cells <- 0

visited_cells <- c()
traversed_set <- c()

while(length(visited_cells) < cells_to_traverse){
 
  cur_cellid <- paste0(latitudes[y_coord],longitudes[x_coord])
  traversed_set <- append(traversed_set, cur_cellid)
  
  skip <- FALSE
  
  if((! cur_cellid %in% cur_codes) | (cur_cellid %in% visited_cells)) skip <- TRUE
  if(! skip){
    if(is.na(sum(master[which(master$PAGENAME == cur_cellid),-1]))) skip <- TRUE
  }
  
  if(! skip){
    visited_cells <- append(visited_cells, cur_cellid)
  }
  
  if(x_count > 0){
    x_index = (x_index+x_offset) %% cols
    x_coord <- x_index + 1;
    x_count <- x_count - 1
    
    if(x_count == 0){
      x_next <- x_next + 1
      if(x_offset == 1) x_offset = -1
      else if(x_offset == -1) x_offset = 1
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
  
  traversed_cells <- traversed_cells + 1
  
  print(traversed_cells)
  
  # if(traversed_cells %% 1000 == 0){
  #   region$SPP <- 0
  #   region[which(region$PAGENAME %in% traversed_set),]$SPP <- 1
  #   
  #   bioregion_raster <- fun.dbf2raster(SPPPA = region, mask.dir = "./mask10k/")
  #   
  #   plot(bioregion_raster, main=region_name)
  # }
  
}


###### here we start with network analyses 

##### first over the whole network

int_matrix <- BARMdiet.binary

rownames(int_matrix) <- unlist((apply(as.matrix(rownames(int_matrix)), 1, function(x){ n <- whois(SPPCODE=x); if(length(n) == 0) {return (x)} else{return (n)} }  ) ))
colnames(int_matrix) <- unlist((apply(as.matrix(colnames(int_matrix)), 1, function(x){ n <- whois(SPPCODE=x); if(length(n) == 0) {return (x)} else{return (n)} }  ) ))

### as an adjacency matrix
int_matrix <- t(int_matrix)

L <- sum(int_matrix)
S <- dim(int_matrix)[1]
C <- L/(S**2)
L.S <- L/S

require(igraph)

### as an igraph network
network <- graph.adjacency(int_matrix)
igraph::walktrap.community(network)
igraph::cluster_spinglass(network)
igraph::cluster_fast_greedy(as.undirected(network))
igraph::cluster_infomap(network)
igraph::cluster_edge_betweenness(network)


#### compute modularity using the netcarto algorithm (Guimerá and Amaral)
require(rnetcarto)

##### we need to create a new adjaceny matrix because in the original one col and row names
##### are not in the same order. We do this because netcarto requires both col and row names
##### to be in the same order
new_int_matrix <- as.matrix(igraph::as_adj(network))
mod_mat <- netcarto(new_int_matrix)



#### with adjacency list...
adj_list <- as_edgelist(network)
mod_adj <- netcarto(list(adj_list[,1], adj_list[,2]))


### netcarto gets as input an adjacency list (above)
### or an upper triangular network. Below we create the undirected representation
### of the network to create an upper triangular matrix. Then use this as an input.
upper_tr_net <- as.matrix(as_adj(as.undirected(network, mode='each'), type='upper'))
mod_undirected <- netcarto(upper_tr_net)


###### this part of the code calculates the division in groups based on
###### Allesina's algorithm with and without Solow & Beet's group number selection
###### heuristic

##### to be implemented!!!!!!!

write.table(new_int_matrix, file='adjacency-matrix-europe', sep=" ", row.names = F, col.names = F)

#### for some reason we obtain different values of modularity (using netcarto) 
#### when calculated over the adjacency matrix vs when calculated using edgelist

##### to dig deeper into this issue I use the netcarto library from python and see if
##### there are differences

library(rPython)

python.load('../network-area/src/py_script.py')

python.assign('seed', as.integer(floor(runif(1, 1, 100000001))))
python.assign('iter_factor', 1.0)
python.assign('cooling_factor', 0.995)
python.assign('randoms', 0)
python.assign('file', 'net-edge-list')
python.exec('mod = modularity_rgraph(file, seed, iter_factor, cooling_factor, randoms)')
mod_adjlist_python <- python.get('mod')




###### let's repeat the modularity analyses without the all-inclusive nodes

pan_nodes <- c("Mushrooms", "MossesLichens", "Algae", "Detritus", "SeedsNutsGrains", 
               "Fruits", "OtherPlantParts", "Invertebrates", "Fish", "DomesticAnimals",
               "Carrion", "Coprofagous")

remove_rows <- which(rownames(int_matrix) %in% pan_nodes)
remove_cols <- which(colnames(int_matrix) %in% pan_nodes)

int_matrix_no_hubs <- int_matrix[-remove_rows, -remove_cols]

L2 <- sum(int_matrix_no_hubs)
S2 <- dim(int_matrix_no_hubs)[1]
C2 <- L2/(S2**2)
L.S2 <- L2/S2

### as an igraph network
network_no_hubs <- graph.adjacency(int_matrix_no_hubs)
igraph::walktrap.community(network_no_hubs)
igraph::cluster_spinglass(network_no_hubs)
igraph::cluster_fast_greedy(as.undirected(network_no_hubs))
igraph::cluster_infomap(network_no_hubs)


new_int_matrix <- as.matrix(igraph::as_adj(network_no_hubs))
mod_mat_no_hubs <- netcarto(new_int_matrix)

#### with adjacency list...
adj_list <- as_edgelist(network_no_hubs)
mod_adj_no_hubs <- netcarto(list(adj_list[,1], adj_list[,2]))


### netcarto gets as input an adjacency list (above)
### or an upper triangular network. Below we create the undirected representation
### of the network to create an upper triangular matrix. Then use this as an input.
upper_tr_net <- as.matrix(as_adj(as.undirected(network_no_hubs, mode='each'), type='upper'))
mod_undirected_no_hubs <- netcarto(upper_tr_net)






################################################################
######## Modules analysis in relation to bioregions ########
################################################################

##### let's find out what bioregions are species located in
species_codes <- names(master)[-1]
region <- data.frame(PAGENAME = master$PAGENAME, SPP = 0)
regions_to_remove <- c('outside')

absent_species <- c("B511", "B512", "R250")

species_regions <- NULL
# species_regions_unique <- data.frame(species_codes, region=NA, cell_count=0)

species_regions_unique <- data.frame(species_codes, primary_region=NA, fraction_region=0)

for(i in shape@data$PK_UID){
  
  if(shape@data[which(shape@data$PK_UID == i),]$short_name %in% regions_to_remove) next;
  
  cur_pol <- shape[i,]
  region_name <- as.character(shape@data[which(shape@data$PK_UID == i),]$code)
  cur_ids <- extract(europeRaster, cur_pol)[[1]]
  cells_codes_region <- as.character(cells_info[which(cells_info$Value %in% cur_ids),]$PageName)
  
  region$SPP <- 0
  region[which(region$PAGENAME %in% cells_codes_region),]$SPP <- 1
  
  bioregion_raster <- fun.dbf2raster(SPPPA = region, mask.dir = "./mask10k/")
  
  plot(bioregion_raster, main=region_name)
  cells_in_region <- length(cells_codes_region)
  for(sp_code in species_codes){
    if(!sp_code %in% absent_species){
      cells_codes_species <- as.character(master[which(master[sp_code] == 1),]$PAGENAME)
      
      # sppDIST <- data.frame(PAGENAME = master$PAGENAME, SPP = master[sp_code])
      
      # if(length(intersect(cells_codes_region, cells_codes_species)) > 0){
      #   cur_out <- data.frame(species_code=sp_code, species_name=as.character(whois(SPPCODE = sp_code)), region=region_name)
      # 
      #   if(is.null(species_regions)) species_regions <- cur_out
      #   else species_regions <- rbind(species_regions, cur_out)
      # 
      # }
      
      species_cells <- length(cells_codes_species)
      
      common_cells <- length(intersect(cells_codes_region, cells_codes_species))
      # if(common_cells > species_regions_unique[which(species_regions_unique$species_codes == sp_code),]$cell_count){
      #   species_regions_unique[which(species_regions_unique$species_codes == sp_code),]$region <- region_name
      #   species_regions_unique[which(species_regions_unique$species_codes == sp_code),]$cell_count <- common_cells
      # }
      
      fraction_cells <- (species_cells/cells_in_region)
      if(common_cells > 0 & species_regions_unique[which(species_regions_unique$species_codes == sp_code),]$fraction_region < fraction_cells){
        species_regions_unique[which(species_regions_unique$species_codes == sp_code),]$primary_region <- region_name
        species_regions_unique[which(species_regions_unique$species_codes == sp_code),]$fraction_region <- fraction_cells
      }
      
      
      
      
    }
    # sppDIST$PAGENAME <- as.character(sppDIST$PAGENAME) # important: PAGENAME must be as character for this to work correctly
    # 
    # sppRaster <- fun.dbf2raster(SPPPA = sppDIST, mask.dir = "./mask10k/")
    # 
    # plot(sppRaster, main= whois(SPPCODE = sp_code)$SPPname)
    
    # cur_pol <- shape[i,]
    # region_name <- as.character(shape@data[which(shape@data$PK_UID == i),]$code)
    # cur_ids <- extract(europeRaster, cur_pol)[[1]]
    # 
    # 
    # cur_codes <- as.character(cells_info[which(cells_info$Value %in% cur_ids),]$PageName)
    # 
    # region$SPP <- 0
    # region[which(region$PAGENAME %in% cur_codes),]$SPP <- 1
    # 
    # bioregion_raster <- fun.dbf2raster(SPPPA = region, mask.dir = "./mask10k/")
    # 
    # plot(bioregion_raster, main=region_name)
    
  }
}

species_regions_unique$species_name <- NA
for(idx in 1:(dim(species_regions_unique)[1])){
  code <- as.character(species_regions_unique[idx,]$species_codes)
  if(!code %in% absent_species)
    species_regions_unique[idx,]$species_name <- whois(SPPCODE = code)
}


par(mar=c(4,4,4,4))
bins <- c(0,1,2,3,4,5,6,7,8,9,10,11)
hist(table(species_regions$species_name), main='Distribution of regions per species', xlab='Number of regions per species', breaks=bins, xlim=c(0,12))

output_netcarto <- mod_undirected[[1]]

output_netcarto$regions <- 0

for(sp in output_netcarto$name){
  output_netcarto[which(output_netcarto$name == sp),]$regions <- length(which(species_regions$species_name == sp))
}

###### these are the species that for some reason do not have a distribution

setdiff(output_netcarto$name, as.character(species_regions$species_name))

require(ggplot2)

ggplot(output_netcarto, aes(as.factor(regions), degree)) + geom_boxplot() + theme_bw()

ggplot(output_netcarto, aes(role, (degree))) + geom_boxplot() + theme_bw()


output_netcarto$degree <- 0
output_netcarto$degree <- degree(network, output_netcarto$name, mode='all')

ggplot(output_netcarto[-which(output_netcarto$regions == 0),], aes(degree, regions)) + geom_point() + theme_bw() + geom_smooth(method='lm')



#### to explore whether modules accumulate species from the same habitats
#### we do some pie charts

#### we iterate over modules
par(mar=c(0,0,2,0))
for(m in unique(output_netcarto$module)){
  
  sps_in_m <- as.character(output_netcarto[which(output_netcarto$module == m),]$name)
  
  regions_in_m <- as.character(species_regions_unique[which(species_regions_unique$species_name %in% sps_in_m),]$primary_region)
  
  counts <- sort(table(regions_in_m))
  
  # slices <- c(10, 12, 4, 16, 8) 
  lbls <- names(counts)
  pct <- round(counts/sum(counts)*100)
  lbls <- paste(lbls, pct) # add percents to labels 
  lbls <- paste(lbls,"%",sep="") # ad % to labels 
  pie(counts, labels = lbls, col=rainbow(length(lbls)),
      main=paste("Composition of regions in module", m+1))
}





















