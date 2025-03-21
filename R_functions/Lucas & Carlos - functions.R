

#########################################################################################################
#Function to filter selected columns
#########################################################################################################


# Define a function to filter selected columns from a dataset
filter_columns = function(data, selected_columns) {
  # Ensure only existing columns are selected
  existing_columns = intersect(selected_columns, colnames(data))
  filtered_data = data[, existing_columns, drop = FALSE]
  
  # Return the filtered dataset
  return(filtered_data)
}



#########################################################################################################
#########################################################################################################
#########################################################################################################



#########################################################################################################
#Function to search sites
#########################################################################################################

search.site = function(data, site_coding=F){
  #checking input:
  var_acept = c("Species", "MaxT", "MinT")
  if (sum(colnames(data) %in% var_acept)<3){
    stop("Please assign data column names in the format: \"Species\" \"Status\", \"MaxT\", \"MinT\"")
  }
  
  sp_col = which(colnames(data)=="Species")
  info_site = data[,-sp_col]
  sites = unique(do.call(paste, c(info_site, sep="")))
  matchs = as.data.frame(do.call(paste, c(info_site, sep="")))
  site_number = apply(matchs, 1, match, table=sites)
  
  data$Site = site_number
  if (site_coding == F){
    return(data)
  }else{
    coding = cbind(data, site_info, site_number)
    result = list(data = data, site_coding = coding)
    return(result)
  }
}


#########################################################################################################
#########################################################################################################
#########################################################################################################


#########################################################################################################
# Function to calculate resolution and plot distribution without filtering any data
#########################################################################################################

# This function will mark the occurences that have a temporal resolution bigger than the threshold 
# informed. At the end, you will have two more columns in the dataset: Resolution and Resolution_flag
# The Resolution column is the MaxT - MinT
# The Resolution_flag will informe if the occurence could be problematic or not, based on your treshold
# and the user can create a new dataset latter with only the occurences without flags

# ARGUMENTS:
# data = dataset containing occurrence data
# threshold = temporal value defining the maximum accepted resolution

filter_resolution_and_plot = function(data, threshold) {
  library(ggplot2)
  library(dplyr)
  library(stringr)
  
  # Ensure MaxT and MinT are numeric before calculations
  data = data %>%
    mutate(MaxT = as.numeric(MaxT),
           MinT = as.numeric(MinT))
  
  # Calculate lifespan (MaxT - MinT) for each occurrence
  data = data %>%
    mutate(Resolution = MaxT - MinT)
  
  # Add a flag "Resolution_flag": TRUE if Resolution > threshold, FALSE otherwise
  data = data %>%
    mutate(Resolution_flag = Resolution > threshold)
  
  # Define x-axis limit to ensure consistency across plots
  x_limit = range(data$Resolution, na.rm = TRUE)
  
  # Generate histogram of all dataset resolution
  hist_obj = hist(data$Resolution, 
                  xlab="Resolution (Myr)", 
                  main="Resolution distribution", 
                  col=rgb(1, 0, 0, 0.5), 
                  border="black",
                  breaks=20, 
                  xlim=x_limit)
  
  # Overlay occurrences that exceed the threshold in blue
  hist(data$Resolution[data$Resolution > threshold], 
       col="blue",
       border="black",
       breaks=hist_obj$breaks, 
       add=TRUE)
  
  # Add legend
  legend("topright", legend=c("Below Threshold", "Above Threshold"), 
         fill=c(rgb(1, 0, 0, 0.5), "blue"), 
         border="black")
  
  # Return the full dataset with additional columns
  return(data)
}


#########################################################################################################
#########################################################################################################


#########################################################################################################
# Function to compute excess resolution, plot histogram, and return counts as a table
#########################################################################################################

plot_excess_resolution_and_counts <- function(data, threshold) {
  
  # Ensure the dataset has the necessary columns
  if (!("Resolution" %in% colnames(data)) || !("Resolution_flag" %in% colnames(data))) {
    stop("The dataset must contain 'Resolution' and 'Resolution_flag' columns. Run the previous processing steps first.")
  }
  
  # Create the Resolution_excess column and round values
  data <- data %>%
    mutate(Resolution_excess = ifelse(Resolution_flag, round(Resolution - threshold), NA))
  
  # Filter only occurrences where Resolution_flag is TRUE
  data_excess <- data %>% filter(Resolution_flag == TRUE)
  
  # Generate histogram for excess resolution values
  plot <- ggplot(data_excess, aes(x = Resolution_excess)) +
    geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
    labs(title = "Excess Resolution Beyond Threshold",
         x = "Myr Exceeding Threshold",
         y = "Frequency") +
    theme_minimal()
  
  # Print the plot
  print(plot)
  
  # Create a table where columns represent Resolution_excess values and the row is the count
  excess_counts <- table(data_excess$Resolution_excess) %>% as.data.frame()
  excess_counts <- spread(excess_counts, key = Var1, value = Freq, fill = 0)
  
  # Return the table
  return(excess_counts)
}



#########################################################################################################
#########################################################################################################



#########################################################################################################
# Function to get species list grouped by Resolution_excess intervals
#########################################################################################################

get_species_vectors_by_excess_interval <- function(data, threshold) {
  
  # Ensure necessary columns exist
  if (!all(c("Resolution", "Resolution_flag", "Species") %in% colnames(data))) {
    stop("The dataset must contain 'Resolution', 'Resolution_flag', and 'Species' columns.")
  }
  
  # Create Resolution_excess column (rounded)
  data <- data %>%
    dplyr::mutate(Resolution_excess = ifelse(Resolution_flag, round(Resolution - threshold), NA))
  
  # Check if Resolution_excess was created properly
  if (!"Resolution_excess" %in% colnames(data)) {
    stop("Failed to create 'Resolution_excess' column.")
  }
  
  # Filter only TRUE values for Resolution_flag
  data <- data %>% filter(Resolution_flag == TRUE)
  
  # Create a named list where each Resolution_excess interval has a vector of species
  species_by_excess <- data %>%
    dplyr::group_by(Resolution_excess) %>%
    dplyr::summarise(Species_Vector = list(unique(Species)), .groups = "drop") %>%
    dplyr::arrange(Resolution_excess) %>%
    dplyr::pull(Species_Vector, name = Resolution_excess)  # Convert to named list
  
  return(species_by_excess)
}



#########################################################################################################
#########################################################################################################


#########################################################################################################
# Function to find pseudoreplicas
#########################################################################################################


find.pseudoreplicas = function(dtst1, dtst2=NULL, lat_threshold=Inf, long_threshold=Inf, distance_threshold=Inf, time_threshold=Inf, JW_threshold=1, print_progress=T, print_i=T){
  
  if(is.null(dtst2)){
    dtst2=dtst1
  }
  
  if( sum(colnames(dtst1) %in% c("Species", "MaxT", "MinT", "lng", "lat")) < 5 ){
    stop("The columns in dtst1 must be named in the format \"Species\", \"MaxT\", \"MinT\", \"lng\", \"lat\"")
  }
  
  if( sum(colnames(dtst2) %in% c("Species", "MaxT", "MinT", "lng", "lat")) < 5 ){
    stop("The columns in dtst2 must be named in the format \"Species\", \"MaxT\", \"MinT\", \"lng\", \"lat\"")
  }
  
  if(JW_threshold <0 | JW_threshold >1){
    stop("JW_threshold must be between 0 and 1")
  }
  
  #adjusting data
  
  corr_ids=c(which(colnames(dtst1)=="Species"),
             which(colnames(dtst1)=="MinT"),
             which(colnames(dtst1)=="MaxT"),
             which(colnames(dtst1)=="lng"),
             which(colnames(dtst1)=="lat"))
  dtst1=dtst1[,corr_ids]
  
  corr_ids=c(which(colnames(dtst2)=="Species"),
             which(colnames(dtst2)=="MinT"),
             which(colnames(dtst2)=="MaxT"),
             which(colnames(dtst2)=="lng"),
             which(colnames(dtst2)=="lat"))
  dtst2=dtst2[,corr_ids]
  
  
  results=list()
  ids=vector()
  
  if(print_progress){pb <- txtProgressBar(min = 1, max=nrow(dtst1))}
  for(i in 1:nrow(dtst1)){
    JW=RecordLinkage::jarowinkler(str1 = dtst1$Species[i], str2= dtst2$Species)
    
    if(sum(JW>=JW_threshold)>0){ #if there is more than one occurrence from that lineage, given the JW threshold:
      
      #i=55
      #dtst1=now_for_OUR
      #dtst2=pbdb_for_OUR
      #time_threshold=1
      #JW_threshold=1
      #long_threshold=.1
      #lat_threshold=.1
      
      
      #separate which occs in dtst2 are within that threshold
      lines=which(JW>=JW_threshold)
      aux=dtst2[lines,]
      
      #for test:
      #aux=data.frame(Species="Eotragus_artenensis", MaxT=16.8, MinT=15.97, lng=1.856, lat=49, dtst="now")
      
      #difference between the occs in dtst2 and the occ in dtst1:
      aux_2=data.frame(diff_MaxT=abs(aux$MaxT-dtst1$MaxT[i]),
                       diff_MinT=abs(aux$MinT-dtst1$MinT[i]),
                       diff_lng=abs(aux$lng-dtst1$lng[i]),
                       diff_lat=abs(aux$lat-dtst1$lat[i]))
      
      
      aux_log=cbind(aux_2[,1:2]<=time_threshold,
                    aux_2$diff_lng<long_threshold,
                    aux_2$diff_lat<lat_threshold,
                    apply(aux[,4:5], 1, geosphere::distGeo, p2=c(dtst2$lng[i], dtst2$lat[i])) < distance_threshold)
      aux_log=as.data.frame(aux_log)
      
      #now we do the logic test to see if the occurrences may be pseudoreplicas
      test=(apply(aux_log[,1:2], 1, sum)==2) & #time
        (apply(aux_log[,3:4], 1, sum)==2) & #space
        aux_log[,5] #distance
      
      #for test:
      #aux_2[test,]
      
      if(nrow(aux_2[test,])>0){
        ids=c(ids, i)
        res=rbind.fill(dtst1[i,], aux[test,])
        res$dataset=c(1, rep(2, times=nrow(res)-1))
        res$row=c(i, lines[test])
        results[[i]]=res
      }
    }
    
    if(print_progress==T){
      setTxtProgressBar(pb, i)
    }
  }
  
  
  
  message(paste0(length(results)-sum(unlist(lapply(results, is.null))),
                 " (",100-100*round(sum(unlist(lapply(results, is.null)))/length(results), digits = 2) ,
                 "%) possible pseudoreplicas found in dtst1" ))
  return(plyr::compact(results))
}




#########################################################################################################
#########################################################################################################


#########################################################################################################
# Function to removed the indicated pseudoreplicas from one database
#########################################################################################################

remove.pseudoreplicas = function(pseudorep, dtst1){
  id=vector()
  for(i in 1:length(pseudorep)){
    id=c(id, pseudorep[[i]]$row[1])
  }
  dtst1=dtst1[-id,]
  return(dtst1)
}

  
#########################################################################################################
######################################################################################################### 
