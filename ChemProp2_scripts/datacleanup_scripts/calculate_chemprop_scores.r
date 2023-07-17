prepare_meta_chemprop <- function(subset_md, new_md, md_sequence) {
    
    # the user has to give 2 metadata names subset_md, new_md, the column number in the corresponding md for meta_ChemProp. 
    #If they don't have subset_md, the function will not choose that
    
  if(exists("subset_md")==T){
    my_data <- subset_md
  } else {
    my_data <- new_md
  }
  Meta_ChemProp <- my_data %>% select(contains(colnames(my_data)[md_sequence]))
  
  # Removing any characters in the Meta_ChemProp column and converting it to numeric type
  Meta_ChemProp[,1] <- as.numeric(gsub("\\D+", "", Meta_ChemProp[,1]))
  rownames(Meta_ChemProp) <- my_data$filename

  View(head(Meta_ChemProp,3))
  message("The number of rows and columns in Meta_ChemProp column is: ", paste(dim(Meta_ChemProp), collapse = " x "), "\n")
    
  #seeing the levels of the attribute for ChemProp calculation  
  message("The levels in the Meta_ChemProp column is:", paste(levels(as.factor(Meta_ChemProp[,1])), collapse = ", ")) 

  return(Meta_ChemProp)
}

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
prepare_ChemProp2_feature_table <- function(subset_ft, new_ft) {
  if(exists("subset_ft") == TRUE){
    Feature_file <- subset_ft
  } else {
    Feature_file <- new_ft
  }
  
  # Check for existence of column 'row ID', 'rowID', or 'RowID'
  row_id_colname <- intersect(tolower(colnames(Feature_file)), c("row id", "rowid", "row.id"))
    
  if(length(row_id_colname) > 0) {
    # If the row IDs have an 'X' prefix, remove it
    if(all(grepl("^X", Feature_file[[row_id_colname]]))) {
      rownames(Feature_file) <- gsub("^X", "", Feature_file[[row_id_colname]])
    } else {
      # Otherwise, use the row IDs directly
      rownames(Feature_file) <- Feature_file[[row_id_colname]]
    }
    # Remove the row ID column
    Feature_file <- Feature_file[, !(tolower(colnames(Feature_file)) %in% c("row id", "rowid", "row.id"))]
    message("Your Feature file had a 'row ID' column which was taken as its row names.")
  } else if (grepl("_", rownames(Feature_file)[1], fixed = TRUE)) {
      # Replace row names with the ID part (before the first underscore)
      rownames(Feature_file) <- sub("_.*", "", rownames(Feature_file))
      # Remove 'X' prefix from both row and column names
      colnames(Feature_file) <- gsub("^X", "", colnames(Feature_file))
      rownames(Feature_file) <- gsub("^X", "", rownames(Feature_file))
      message("Your Feature file had combined row names (ID_mz_RT_annotation). Only the 'ID' info was taken as row names.")
  }
  
  message("The number of rows and columns in the Feature_table is: ", paste(dim(Feature_file), collapse = " x "))
  
  View(head(Feature_file))
  return(Feature_file)
}












#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
calculate_chemprop2_scores <- function(Feature_file, Meta_ChemProp, nw){
  ChemProp2 <- c()
  ChemProp_spearman <-c()
  ChemProp_log <- c()
  ChemProp_sqrt <- c()

  for (i in 1:nrow(nw)) {
      
      # rownames(Feature_file) is the feature ID or cluster ID. 
      #The subset command gets the 'Feature ID 1' from the first column of Nw_edge. Then picks the row from the Feature_file corresponding to the 'Feature ID 1'
      
      x <- subset(Feature_file, rownames(Feature_file) == nw[i,1]) 
      x <- rbind(x,subset(Feature_file, rownames(Feature_file) == nw[i,2]))
      
      # x is the subset data which has the Feature ID 1 and 2 specified according to Nw_edge file.
    
      #Match gives the position in which B (the column names of Meta data) is present in A (subset data) and store the position info in reorder_id 
      reorder_id <- match(rownames(Meta_ChemProp),colnames(x))
      
      #Rearranging x (subset data) with respect to the new positions and transposing it
      reordered_x <- data.frame(t(x[reorder_id]))
      
      #Combining the metadata column (here, timepoint) with reordered_x
      reordered_x <- cbind(Meta_ChemProp[,1],reordered_x)
      reordered_x <- as.matrix(reordered_x)
    
      #Thus, the resulting reordered_x contains 3 columns, such as: 'Metadata info(eg., Timepoint)', 'Feature ID 1', 'Feature ID 2'
      
      # Performing Pearson correlation
      corr_result <- cor(reordered_x, method = "pearson")
      ChemProp_score <- (corr_result[1,3] - corr_result[1,2]) / 2 # ChemProp2 score is obtained by: (Pearson(Feature ID 2) - Pearson(Feature ID 1)) / 2

      # Performing Spearman correlation
      corr_2 <- cor(reordered_x, method = "spearman")
      Score_spearman <- (corr_2[1,3] - corr_2[1,2]) / 2
  
      # Performing natural log transformations on Feature IDs 1 and 2
      log_reorderedX <- cbind(reordered_x[,1],log(reordered_x[,2:3]+1))
      corr_3 <- cor(log_reorderedX)  # performing (pearson) correlation on the log transformed data
      Score_log <-(corr_3[1,3] - corr_3[1,2]) / 2
  
      # Taking square roots of Feature IDs 1 and 2
      sqrt_reorderedX <- cbind(reordered_x[,1],sqrt(reordered_x[,2:3]))
      corr_4 <- cor(sqrt_reorderedX) # performing (pearson) correlation on the square roots
      Score_sqrt <- (corr_4[1,3] - corr_4[1,2])/2
      
      ChemProp2 <- rbind(ChemProp2, ChemProp_score, deparse.level = 0) 
      # deparse.level = 0 constructs no labels; if not given, the resultant matrix has row names (for all rows) created from the input arguments such as 'ChemProp_score' here.
      
      ChemProp_spearman <- rbind(ChemProp_spearman,Score_spearman,  deparse.level = 0)
      ChemProp_log <- rbind(ChemProp_log,Score_log,  deparse.level = 0)
      ChemProp_sqrt <- rbind(ChemProp_sqrt, Score_sqrt, deparse.level = 0)
  }
    
    nw_new <- cbind (nw, ChemProp2,ChemProp_spearman,ChemProp_log,ChemProp_sqrt)
    rownames(nw_new) <- NULL

    Abs_values <- abs(nw_new[,which( colnames(nw_new)=="ChemProp2" ):length(nw_new)])
    colnames(Abs_values) <- paste("abs", colnames(Abs_values), sep = "_")

    Sign_ChemProp2 <- sign(nw_new$ChemProp2) #getting only the sign of ChemProp2 as 1 or -1
         
    ChemProp2_file <- cbind(nw_new,Abs_values,Sign_ChemProp2)
  
    return(ChemProp2_file)
}