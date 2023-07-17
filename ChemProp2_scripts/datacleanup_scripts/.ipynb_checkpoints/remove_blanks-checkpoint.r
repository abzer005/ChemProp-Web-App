remove_blanks <- function(ft_t, new_md) {
    
    sample_type_col <- readline('Enter the column name for Sample Type in metadata: ') 
    levels <- unique(new_md[[sample_type_col]])  # Get unique levels in the sample_type_col
    levels_df <- data.frame(ID=1:length(levels), levels=levels)
    
    # Display the levels to the user
    message("Please enter the ID numbers corresponding to the labels for blanks and samples, separated by commas:\n")
    View(levels_df)
    
    # Get the user's choices for blanks and samples
    blank_labels_indices <- readline('Enter the numbers for blanks: ')
    sample_labels_indices <- readline('Enter the numbers for samples: ')
    
    # Convert the choices to numeric and get the corresponding labels
    blank_labels <- levels[as.numeric(strsplit(blank_labels_indices, ",")[[1]])]
    sample_labels <- levels[as.numeric(strsplit(sample_labels_indices, ",")[[1]])]
    
    # Getting the blanks and samples based on the metadata
    md_Blank <- new_md %>% filter((!!sym(sample_type_col)) %in% blank_labels)
    Blank <- ft_t[which(rownames(ft_t) %in% (md_Blank$`filename`)),,drop=F]
    
    md_Samples <- new_md %>% filter((!!sym(sample_type_col)) %in% sample_labels)
    Samples <- ft_t[which(rownames(ft_t) %in% (md_Samples$`filename`)),,drop=F]
    
    # When cutoff is low, more noise (or background) detected; With higher cutoff, less background detected, thus more features observed
    Cutoff <- as.numeric(readline('Enter Cutoff value between 0.1 & 1: ')) 
    # (i.e. 10% - 100%). Ideal cutoff range: 0.1-0.3

    # Getting mean for every feature in blank and Samples in a data frame named 'Avg_ft'
    Avg_ft <- data.frame(Avg_blank=colMeans(Blank, na.rm= F)) # set na.rm = F to check if there are NA values. When set as T, NA values are changed to 0
    Avg_ft$`Avg_samples` <- colMeans(Samples, na.rm= F) # adding another column 'Avg_samples' for feature means of samples

    # Getting the ratio of blank vs Sample
    Avg_ft$`Ratio_blank_Sample` <- (Avg_ft$`Avg_blank`+1)/(Avg_ft$`Avg_samples`+1)

    # Creating a bin with 1s when the ratio>Cutoff, else put 0s
    Avg_ft$`Bg_bin` <- ifelse(Avg_ft$`Ratio_blank_Sample` > Cutoff, 1, 0 )

    # Calculating the number of background features and features present
    print(paste("Total no.of features:",nrow(Avg_ft)))
    print(paste("No.of Background or noise features:",sum(Avg_ft$`Bg_bin` ==1,na.rm = T)))
    print(paste("No.of features after excluding noise:",(ncol(Samples) - sum(Avg_ft$`Bg_bin` ==1,na.rm = T))))
                                    
    blk_rem <- merge(as.data.frame(t(Samples)), Avg_ft, by=0) %>%
                        filter(Bg_bin == 0) %>% #picking only the features
                        select(-c(Avg_blank,Avg_samples,Ratio_blank_Sample,Bg_bin)) %>% #removing the last 4 columns
                        column_to_rownames(var="Row.names")

    write.csv(blk_rem, paste0(Sys.Date(),'_Blanks_Removed_with_cutoff_',Cutoff,'.csv'),row.names =TRUE)
    
    message("Your Blank removed table: head(blk_rem,3)")
    View(head(blk_rem,3))
    
     message("head(md_Samples,3)")
    View(head(md_Samples,3))
    
    # return a list containing the two objects
    return(list(blk_rem = blk_rem, md_Samples = md_Samples, Cutoff = Cutoff))
}

impute_table <- function(blk_rem){
    
    #creating bins from -1 to 10^10 using sequence function seq()
    bins <- c(-1,0,(1 * 10^(seq(0,10,1)))) 

    #cut function cuts the give table into its appropriate bins
    scores_gapfilled <- cut(as.matrix(blk_rem),bins, labels = c('0','1',paste("1E",1:10,sep="")))

    #transform function convert the tables into a column format: easy for visualization
    FreqTable <- transform(table(scores_gapfilled)) #contains 2 columns: "scores_x1", "Freq"
    FreqTable$Log_Freq <- log(FreqTable$Freq+1) #Log scaling the frequency values
    colnames(FreqTable)[1] <- 'Range_Bins' #changing the 1st colname to 'Range Bins'

    ## GGPLOT2
    BarPlot <- ggplot(FreqTable, aes(x=Range_Bins, y=Log_Freq)) + 
    geom_bar(stat = "identity", position = "dodge", width=0.3) + 
    ggtitle(label = "Frequency plot - Gap Filled") +
    xlab("Range") + 
    ylab("(Log)Frequency") + 
    theme(plot.title = element_text(hjust = 0.5))
    
    
    Cutoff_LOD <- round(min(blk_rem[blk_rem > 0]))
    message(paste0("The limit of detection (LOD) is: ",Cutoff_LOD)) 
    
    message('Creating random variables (between 0 and LOD) to replace all zeros...')

    set.seed(141222) # by setting a seed, we generate the same set of random number all the time
    ran_val <-round(runif(length(blk_rem),0,Cutoff_LOD),digits=1)
    message("The first 10 of the generated random values: ", paste(ran_val[1:10], collapse = ", "))
    
    imp <- blk_rem  %>% mutate(across(everything(),
                                      ~replace(., . == 0 , # here '.' represents the dataframe 'blk_rem'
                                           sample(ran_val, size=1)))) # sample() picks 1 random sample from 'ran_val'
    
    message("Resulting imputed table: head(imp,3)")
    View(head(imp, 3))
    
    # Return the plot and imp table in a list
    return(list(Plot = BarPlot, Table = imp))
    }