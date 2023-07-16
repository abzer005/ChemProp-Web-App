   
#The first column in feature table: row ID (the unique ID for each detected feature) is given in different columns in different files:
#          In DB result and DB analog result files, the row ID is given in the column: #Scan#('Compound_Name' column has the annotation information)
#          For SIRIUS and CANOPUS summary files, the row ID of the feature table is given in the column 'id'. A typical feature 'id' would be: "3_ProjectName_MZmine3_SIRIUS_1_16", where the last number 16 representing the row ID.
    


# Define the list of tables
tables <- data.frame(ID=1:5,Tables=c("GNPS annotation", "GNPS analog annotation", "SIRIUS", "CANOPUS", "Escape"))
View(tables)

# Print the options for the user
cat("Please enter the numbers corresponding to the tables you want to merge with the feature table, separated by commas:\n")

flush.console()

# Get the user's choices
choices <- readline(prompt = "Your choices: ")

# Split the choices by comma and convert to numeric
choices <- as.numeric(strsplit(choices, ",")[[1]])

# Subset the tables based on the user's choices
chosen_tables <- tables[choices,]
View(chosen_tables)

# Initialize the merged table as the feature table
merge_table <- ft

# Check the user's choices and perform the corresponding actions
if (1 %in% choices) {
  an <- gnps # Use the already loaded gnps
  
  if (identical(class(merge_table$`row ID`), class(an$`#Scan#`))) {
      merge_table <- merge(merge_table, an, by.x="row ID", by.y="#Scan#", all.x= TRUE)
      cat("This is the merged table 'merge_table'")
      View(head(merge_table, 2))
  } else {
      print("The classes of the merging columns do not match.")
  }
}

if (2 %in% choices) {
  an <- gnps_an # Use the already loaded gnps_an
    
  if (identical(class(merge_table$`row ID`), class(an$`#Scan#`))) {
      merge_table <- merge(merge_table, an, by.x="row ID", by.y="#Scan#", all.x= TRUE)
      cat("This is the merged table 'merge_table'")
      View(head(merge_table, 2))
  } else {
      print("The classes of the merging columns do not match.")
  }
}

if (3 %in% choices || 4 %in% choices) {
  # Load the SIRIUS or CANOPUS summary table as summary

  # Split the 'id' column at the underscore and keep only the last part
  summary$`id` <- sapply(strsplit(summary$`id`, "_"), tail, 1)
  
  if (identical(class(merge_table$`row ID`), class(summary$`id`))) {
      merge_table <- merge(merge_table, summary, by.x="row ID", by.y="id", all.x= TRUE)
      cat("This is the merged table 'merge_table'")
      View(head(merge_table, 2))
  } else {
      print("The classes of the merging columns do not match.")
  }
}
 

# Now merged_table contains the feature table merged with all the chosen tables