#Goal:  we bring feature table and metadata in the correct format:
# Make the rownames of md and column names of the ft to be the same. 
# Filenames and the order of files need to correspond in both tables, as we will match md attributes to the ft. In that way, both md and ft, can easily be filtered.

#storing the files under different names to preserve the original files
new_ft <- ft 
new_md <- md

colnames(new_ft) <- gsub(' Peak area','',colnames(new_ft)) #removing Peak area extensions from the column names of ft
new_md$filename <- gsub(' Peak area','',new_md$filename) #removing Peak area extensions from the filenames in md

new_ft <- new_ft[order(new_ft$`row ID`),,drop=F] #arranging the rows of ft file by  by ascending order of row ID
new_ft <- new_ft[,colSums(is.na(new_ft))<nrow(new_ft)] #removing if any NA columns present in the ft file
new_md <- new_md[,colSums(is.na(new_md))<nrow(new_md)] #removing if any NA columns present in the md file

#remove the (front & tail) spaces, if any present, from all columns of md
new_md <- data.frame(lapply(new_md, trimws, which = "both"))

#should return TRUE if you have annotation file
if(exists("merge_table")){identical(merge_table$`row ID`,new_ft$`row ID`)} 

#Changing the row names of the files into the combined name as "XID_mz_RT":
rownames(new_ft) <- paste(paste0("X",new_ft$`row ID`),
                          round(new_ft$`row m/z`,digits = 3),
                          round(new_ft$`row retention time`,digits = 3),
                          if(exists("merge_table")){merge_table$Compound_Name}, 
                          sep = '_') 

rownames(new_ft) <- sub("_$", "", rownames(new_ft)) #to remove the trailing underscore at rownames

message('Your new_ft table:')
show(head(new_ft,2))

message('Your new_md table:')
show(head(new_md,2))

# Check if all column names end with 'mzML' or 'mzXML'
mzML_count <- sum(grepl('\\.mzML$', colnames(new_ft)))
mzXML_count <- sum(grepl('\\.mzXML$', colnames(new_ft)))

# If there's a mix of 'mzML' and 'mzXML', stop and print an error message
if (mzML_count > 0 && mzXML_count > 0) {
  stop("Error: The column names contain a mix of 'mzML' and 'mzXML'. Please correct this.")
} else if (mzML_count > 0) { # If all column names end with 'mzML'
  print("You have all mzML files.")
  new_ft <- new_ft[, grepl('\\.mzML$', colnames(new_ft))]
} else if (mzXML_count > 0) { # If all column names end with 'mzXML'
  print("You have all mzXML files.")
  new_ft <- new_ft[, grepl('\\.mzXML$', colnames(new_ft))]
} else {
  stop("Error: No columns with '.mzML' or '.mzXML' found.")
}

# Order the feature table by its column names
new_ft <- new_ft[, order(colnames(new_ft)), drop = FALSE]

# Order the metadata by the filename column
new_md <- new_md[order(new_md$filename), , drop = FALSE]

# Check how many files in the metadata are also present in the feature table
overlap <- sum(new_md$filename %in% colnames(new_ft))
message("Number of files in the metadata also present in the feature table: ", overlap)

# Check if all filenames in the metadata are in the feature table
if (!identical(new_md$filename, colnames(new_ft))) {
  message("The filenames in the metadata and feature table do not match.")
  
  # Identify which file names in the metadata are not in the feature table
  diff <- setdiff(new_md$filename, colnames(new_ft))
  message("The following filenames are in the metadata but not in the feature table:\n", paste(diff, collapse = ", "), "\nPlease check for spelling mistakes, case-sensitive errors, and re-upload the files with correct metadata.")
} else {
  message("All filenames in the metadata match those in the feature table.")
}

# Display the dimensions of the original and new feature tables and new metadata
message("The number of rows and columns in our original feature table is: ", paste(dim(ft), collapse = " x "), "\n")
message("The number of rows and columns in our new feature table is: ", paste(dim(new_ft), collapse = " x "), "\n")
message("The number of rows and columns in our new metadata is: ", paste(dim(new_md), collapse = " x "))

#THE END
