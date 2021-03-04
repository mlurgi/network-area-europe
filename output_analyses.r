##############################################################################################
##
## output_analyses.r
## Analyzing network area relationships across European bioregions
##
## This script performs fits of power functions to the network area relationships obtained from the
## empirical data of the spatial scaling of food web structure across European bioregions.
## It can be used to process outputs from both empirical data analysis ('european-network-area-relationships.r')
## and the correspoding null models (null_model.r).
## It also plots the data as in the figures presented in the corresponding paper (see below).
##
## For this example, an subset of the original dataset was created to illustrate the analyses.
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

# Load required libraries
require(ggplot2)

# We define the power function 
power_function <- function(x,c,z){
  return(c*x^z)
}

# Read the output dataframe
output <- read.table("./partial-output.csv", sep=",", header=TRUE, stringsAsFactors = F)

# We define the properties we want to fit with the power function
properties <- c('species','links','links_per_sp','indegrees','outdegrees','basals','tops','intermediates','overlaps', 'sd_gens','sd_vuls')

# Empty dataframe to store the resutls of the statistical analyses
power_fits <- NULL

# Loop to perform the statistical analyses across all bioregions
for(d in unique(output$region)){
  cur_data <- na.omit(output[output$region==d,])
# And across the selected properties  
  for(p in properties){
    if(!p %in% names(cur_data)) next
    if(sum(cur_data[p], na.rm = T) == 0) next
    
    cur_model <- tryCatch({
      eval(parse(text=paste0('nls(',p,'~power_function(area, c, z), data=cur_data, start=list(c=10, z=.3))')))
    }, warning = function(w) {
      NA
    }, error = function(e) {
      NA
    }, finally = {
      
    })
    
    if(is.na(cur_model)) next
    
    # Selecting the parameters we want to keep from the analyses
    cur_model <- cbind(dataset=d,property=p,as.data.frame(summary(cur_model)$coefficients))
    
    # Storing the data
    if(is.null(power_fits)){
      power_fits <- cur_model
    }else{
      power_fits <- rbind(power_fits, cur_model)
    }
  }
}

# Uncomment the line below in case you want to write the dataframe generated as a .csv in your computer
write.csv(power_fits, "output_power_fits.csv")

# Plot of the scaling of a given network property with area for all biogeographical regions
ggplot(output, aes(area, species, color=as.factor(region))) +
  geom_smooth(se=T) + #stat_summary(fun.data=mean_se, geom="ribbon", alpha=0.25)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"), axis.text=element_text(size=22),
        axis.title=element_text(size= 22 ,face="bold"),
        legend.title = element_blank(),
        legend.box.background = element_rect(colour = "black", size = 1.5),
        #legend.title = element_text(size=12, face="bold"),
        #legend.text = element_text(size = 10, face = "bold"),
        panel.background = element_rect(fill = "white", color = 'black', size = 1.5))+ 
  labs(x = "Area", y = "species", title= "") + 
  guides(color=guide_legend(override.aes=list(fill=NA)))
