#Function: InsideLevels
InsideLevels <- function(metatable){
    LEVELS <- c() #creating empty vector to store information
    typ <-c()
    COUNT <- c()
    for(i in 1:ncol(metatable)){ # for each metadata column
      temp <- as.data.frame(table(metatable[,i])) #table function gives the category in each column and the count of each category
      x <- temp$Var1 #getting the name of each category in every column
      if(is.double(metatable[,i])==T){ # for numeric columns in metadata table, round the category values
          x=round(as.double(x),2)
      } 
      LEVELS <- rbind(LEVELS,toString(x)) # adding all the category values in a row
      COUNT <- rbind(COUNT,toString(temp$Freq)) # getting the frequency of each level in every column
      typ <- rbind(typ,class(metatable[,i])) # getting the class of each column
    }
    out <- data.frame(INDEX=c(1:ncol(metatable)), #creating an output dataframe with 1st column as INDEX
                      ATTRIBUTES=colnames(metatable), #2nd column ATTRIBUTES will be the column name of metadata table
                      LEVELS, #3rd column LEVELS will give the different categories in each ATTRIBUTE
                      COUNT, #4th column COUNT will give the number of files present with each category
                      'ATTRIBUTE_CLASS'=typ,row.names=NULL) #Final column indicating the Class or datatype of each ATTTRIBUTE
    return(out)
  }

#Function: SubsetLevels  
SubsetLevels <- function(metatable){
    IRdisplay::display(InsideLevels(metatable)) # use show() when working in RStudio
    Condition <- as.double(unlist(strsplit(readline("Enter the IDs of interested attributes to subset (separaed by commas if more than one attribute):"), split=',')))
    
    for( i in 1:length(Condition)){
      list_final <- c() 
      #Shows the different levels within each selected condition:
      Levels_Cdtn <- levels(droplevels(as.factor(metatable[,Condition[i]])))
      subset_meta <- data.frame(Index=1:length(Levels_Cdtn),Levels_Cdtn)
      colnames(subset_meta)[2] <- paste("Levels_",colnames(metatable[Condition[i]]))
      IRdisplay::display(subset_meta) # use show() when working in RStudio
      
      
      #Among the shown levels of an attribute, select the ones to keep or exclude:
      Read_cdtn <- readline("Do you want to keep or exclude few conditions? K/E: ")
      if( Read_cdtn=="K"){
        temp <- as.double(unlist(strsplit(readline("Enter the index numbers of condition(s) you want to KEEP (separated by commas):"), split=',')))
        ty <- class(metatable[,Condition[i]])
        list_keep <-Levels_Cdtn[temp, drop=F]
        list_exc <-Levels_Cdtn[-temp, drop=F]
        cat("The condition(s) you want to exclude in ",colnames(metatable)[Condition[i]]," :",list_exc,"\n")
        cat("The condition(s) you want to keep in ",colnames(metatable)[Condition[i]]," :",list_keep,"\n")
        
        }else if(Read_cdtn=="E"){
          temp <- as.double(unlist(strsplit(readline("Enter the index numbers of condition(s) you want to EXCLUDE (separated by commas):"), split=',')))
          ty <- class(metatable[,Condition[i]])
          list_exc <-Levels_Cdtn[temp, drop=F]
          list_keep <-Levels_Cdtn[-temp, drop=F]
          cat("The condition(s) you want to exclude in ",colnames(metatable)[Condition[i]]," :",list_exc,"\n")
          cat("The condition(s) you want to keep in ",colnames(metatable)[Condition[i]]," :",list_keep,"\n")
        
          }else{
            print("Sorry! You have given a wrong input!! Please enter either K or E")
            break
      }
      
      #In order to keep the original datatype of the columns, we define the following conditions, else it would all become 'characters or factors'
      if(ty=="integer"){list_keep <- as.integer(list_keep);list_exc <- as.integer(list_exc)
      }else if(ty=="double"){list_keep <- as.double(list_keep);list_exc <- as.double(list_exc)
      }else if(ty=="numeric"){list_keep <- as.numeric(list_keep);list_exc <- as.numeric(list_exc)}
      
      #Gets all the elements in list_keep into list_final 
      for(j in 1:length(list_keep)){
        sub_list <- metatable[(metatable[,Condition[i]] == list_keep[j]),]
        list_final <- rbind(list_final,sub_list)
      }
      metatable <- list_final #list_final again called as metatable in order to keep it in the for-loop for further subsetting
    }
    return(metatable)
  }
