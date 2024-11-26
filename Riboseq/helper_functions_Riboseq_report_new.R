if (!requireNamespace("types", quietly = TRUE)) {
    install.packages("types", repos = "https://cloud.r-project.org/")
  }

if (!requireNamespace("UpSetR", quietly = TRUE)) {
    install.packages("UpSetR", repos = "https://cloud.r-project.org/")
  }

if (!requireNamespace("lemon", quietly = TRUE)) {
    install.packages("lemon", repos = "https://cloud.r-project.org/")
  }

if (!requireNamespace("cmapR", quietly = TRUE)) {
    if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org/")

    BiocManager::install("cmapR")
  }

library(glue)
library(types)
library(stringr)
library(dplyr)
library(UpSetR)
library(knitr)
library(lemon)
library(cmapR)


calculate_background_threshold <- function(unique_region_type = ? character, input_path){
  #get directory list of subdirectories of the current working directory
  directories<-list.dirs(path = input_path, full.names = TRUE, recursive = TRUE)
  
  #get random intersect counts of the region type of interest
  randomfiles=list()
  for(i in directories){
    if(!(identical(list.files(i,pattern=paste0("*", unique_region_type, "_random_intersect_counts_relative_sorted.bed")), character(0)))){
      randomfiles=c(randomfiles,paste0(i,"/",list.files(i,pattern=paste0("*", unique_region_type, "_random_intersect_counts_relative_sorted.bed"))))
    }
  }
  print(randomfiles)
  
  #read in intersections with random alignments and calculate a threshold for each
  randomdataframes <- lapply(randomfiles, read.csv, header = FALSE, sep="\t")
  background=list()
  for(randomframe in randomdataframes){
    colnames(randomframe)=c("ID","start","stop","read_count", "relative_count")
    robust_z_scores_background <- robust_zscore(randomframe$relative_count)
    #print(robust_z_scores_background)
    filter_zscore_mask <- robust_z_scores_background < 3
    #print(length(filter_zscore_mask))
    #print(length(randomframe$relative_count))
    #print(typeof(randomframe$relative_count))
    #print(randomframe$relative_count[filter_zscore_mask])
    background_filtered <- randomframe$relative_count[filter_zscore_mask]
    #print(background_filtered)
    backgorund_95_percentile <- quantile(background_filtered, probs = c(0.95))
    print(backgorund_95_percentile)
    
    background = c(background,backgorund_95_percentile)
  }
  return(background)
}


print_thresholds <- function(unique_region_type = ? character, background, input_path){
  #get directory list of subdirectories of the current working directory
  directories<-list.dirs(path = input_path, full.names = TRUE, recursive = TRUE)
  #obtain names of the datasets for which the thresholds are calculated
  randomSetnames=list()
  for(i in directories){
    if(!(identical(list.files(i,pattern=paste0("*", unique_region_type, "_random_intersect_counts_relative_sorted.bed")), character(0)))){
      randomSetnames=c(randomSetnames,list.files(i,pattern=paste0("*", unique_region_type, "_random_intersect_counts_relative_sorted.bed")))
    }
  }
  randomnames=stringr::str_replace(randomSetnames, pattern = "_random_intersect_counts_relative_sorted.bed", replacement = "")
  
  printbackground=data.frame(x=background)
  printbackground=as.data.frame(t(printbackground))
  rownames(printbackground)=randomnames
  colnames(printbackground)="Threshold"
  print(kable(printbackground, caption="Thresholds for each dataset"))
}


get_intersect_files <- function(unique_region_type = ? character, input_path){
  directories=list.dirs(path = input_path, full.names = TRUE, recursive = TRUE)
  files=list()
  setnames=list()
  for(i in directories){
    if(!(identical(list.files(i,pattern=paste0("*", unique_region_type, "_intersect_counts_relative_sorted.bed")), character(0)))){
      files=c(files,paste0(i,"/",list.files(i,pattern=paste0("*", unique_region_type, "_intersect_counts_relative_sorted.bed"))))
      setnames=c(setnames,list.files(i,pattern=paste0("*", unique_region_type, "_intersect_counts_relative_sorted.bed")))
    }
  }
  
  dataframes <- lapply(unlist(files), read.csv, header = FALSE, sep="\t")
  datalist=list()
  for(a in dataframes){
    datalist=c(datalist,list(a))
  }
  dataframes=datalist
  names(dataframes) <- stringr::str_replace(setnames, pattern = "_intersect_counts_relative_sorted.bed", replacement = "")
  
  return(dataframes)
}


count_unique_regions_above_threhsold <- function(dataframes, threshold){
  
  relevantregionscount=list()
  relevantregions=list()
  unique_names <- list()
  unique_names_per_sample <- list()
  i=1
  for(frame in dataframes){
    colnames(frame)=c("ID", "start", "stop", "ORF", "read_count", "relative_count")
    temp=0
    c=1
    unique_names_per_sample[[as.character(i)]] <- c()
    for(count in frame$relative_count){
      if(count >= threshold[[i]]){
        if(frame$read_count[[c]] > 2){#read count at least 3 required
          temp = temp + 1
          unique_names <- c(unique_names, paste(frame$ID[[c]], frame$start[[c]], frame$stop[[c]], sep ="-"))
          unique_names_per_sample[[as.character(i)]] <- c(unique_names_per_sample[[as.character(i)]], paste(frame$ID[[c]], frame$start[[c]], frame$stop[[c]], sep ="-"))
        }
      }
      c = c + 1
    }
    i = i + 1
    relevantregionscount = c(relevantregionscount,temp)
  }
  
  
  
  upsetlist=list()
  j=1
  for(f in dataframes){
    bed=f
    bed = filter(bed, V6 >= threshold[[j]])
    bed = filter(bed, V5 > 2)#read count at least 3 required
    unique_regions<- rep(NA, length(bed[,1]))
    for(i in seq_along(1:length(bed[,1]))){
      unique_regions[i]=paste0(bed[i,1],"-",bed[i,2],"-",bed[i,3])
    }
    upsetlist=c(upsetlist,list(unique_regions))
    
    j = j + 1
  }
  
  print(unique_names[!(unique_names %in% unlist(upsetlist))])
  print(length(unique_names) == length(unlist(upsetlist)))
  
  counts_and_names_list <- list()
  counts_and_names_list$counts_above_t_relative_length<- relevantregionscount
  counts_and_names_list$unique_names <- unique_names_per_sample
  return(counts_and_names_list)
}


print_unique_region_above_t <- function(relevantregionscount, dataframes){
  printframe=data.frame(x=relevantregionscount)
  printframe=as.data.frame(t(printframe))
  rownames(printframe)=names(dataframes)
  colnames(printframe)=paste0("Number of unique regions with relative count >= threshold")
  print(kable(printframe, caption="Regions above the threshold", escape = TRUE))
}


get_top_5_unique_regions <- function(dataframes_list, threshold){
  upsetlist=list()
  j=1
  names <- names(dataframes_list)
  for(f in dataframes_list){
    bed=f
    bed = filter(bed, V6 >= threshold[[j]])
    bed = filter(bed, V5 > 2)
    if (dim(bed)[1] > 0){
      unique_regions<- rep(NA, length(bed[,1]))
      for(i in seq_along(1:length(bed[,1]))){
        unique_regions[i]=paste(bed[i,1],bed[i,2],bed[i,3], bed[i,4],sep =":")
      }
      upsetlist=c(upsetlist,list(unique_regions))
      colnames(bed)=c("ID","start","stop","ORF","read_count",  "relative_count")
      
      bed$ID <- paste0(str_split_fixed(bed$ID,  "\\|", 2)[,1], "-", str_split_fixed(bed$ID,  "\\|", 2)[,2])
      print(kable(bed[1:5,], caption=names(dataframes_list)[j], row.names = FALSE, escape = TRUE))
    }
    else {
      names <- names[! names %in% c(names[j])]
    }
    j = j + 1
  }
    names(upsetlist) <- names
    return(upsetlist)
}


count_unique_regions_get_count <- function(dataframes, unique_names_per_sample, outdir){
  relevantregionscount=list()
  frame_number = 1
  for(frame in dataframes){
    colnames(frame)=c("ID","start","stop","ORF","read_count",  "relative_count")
    frame$unique_name <- paste0(frame$ID,"-",frame$start,"-",frame$stop)
    frame$signficant <- ifelse(frame$unique_name %in% unique_names_per_sample[[as.character(frame_number)]], 1, 0)
    write.csv(frame, file.path(outdir, paste(names(dataframes)[frame_number],  "unique_regions.csv", sep = "_")))
    above_5_counter=0
    c=1
    for(count in frame$relative_count){
        if(frame$read_count[[c]] >= 4){
          above_5_counter = above_5_counter + 1
      }
      c = c + 1
    }
    relevantregionscount = c(relevantregionscount, above_5_counter)
    frame_number = frame_number + 1
  }
  
  print = t(as.data.frame(relevantregionscount))
  rownames(print)=names(dataframes)
  colnames(print)="Count get 4"
  print(kable(print, caption="Nr of unique regions with count get 4"))
}

Split_ORFs_validation <- function(dataframes, unique_names_per_sample){
  #create empty dataframe to concatenate the dfs
  df <- data.frame(matrix(ncol = 6, nrow = 0))
  # Assign column names
  colnames(df) <- c("ID", "start", "stop", "ORF", "read_count", "relative_count")
  
  frame_number = 1
  for (frame in dataframes){
    colnames(frame)=c("ID", "start", "stop", "ORF", "read_count", "relative_count")
    frame$unique_name <- paste0(frame$ID,"-",frame$start,"-",frame$stop)
    frame$significant <- ifelse(frame$unique_name %in% unique_names_per_sample[[as.character(frame_number)]], 1, 0)
    #print(frame)
    frame <- frame %>%
      filter(significant == 1)
    
     SO_2_uniques_validated <- frame %>%
      group_by(ID) %>%
      filter(n() >= 2) %>%
      ungroup() %>%
       arrange(ID)
    
    print(SO_2_uniques_validated[, c("ID", "start", "stop", "ORF", "read_count")], max = nrow(SO_2_uniques_validated))
    frame_number <- frame_number + 1
    df <- rbind(df, frame)
    
  }
  df <- df %>%
    group_by(ID) %>%
    filter(n() >= 2) %>%
    ungroup() %>%
    arrange(ID)
  
  
  print(df[, c("ID", "start", "stop", "ORF", "read_count")], max = nrow(df))
}