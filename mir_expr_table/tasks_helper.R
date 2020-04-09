# Purpose: 3 task functions

setwd("/Users/helending/Documents/BCB330/mir_expr_table")


# ============= TASK 1 ============
give_tissues_return_mir <- function(tissues, confidence_level){
  # tissues: character string of tissues
  # confidence: string of characters, comma deliminated. From either "Very High", "High", "Medium", "Low"
  # Return the miR that are in these tissues and also their numbers according to the specified confidence level.

  
  if (confidence_level == "Very High"){
    confidence <- read.table("/Users/helending/Documents/BCB330/mir_expr_table/data/combined_all_expr_tables_0.99.csv")
  } else if (confidence_level == "High"){
    confidence <- read.table("/Users/helending/Documents/BCB330/mir_expr_table/data/combined_all_expr_tables_0.95.csv")
  } else if (confidence_level == "Medium"){
    confidence <- read.table("/Users/helending/Documents/BCB330/mir_expr_table/data/combined_all_expr_tables_2_3.csv")
  } else if (confidence_level == "Low"){
    confidence <- read.table("/Users/helending/Documents/BCB330/mir_expr_table/data/combined_all_expr_tables_1_3.csv")
  }
  confidence <- data.frame(confidence)
  
  tissues_copy <- tissues
  if (regexpr(",", tissues_copy, perl=TRUE) == -1){
    tissues_individual <- tissues_copy
  }else{
    tissues_individual <- unlist(strsplit(tissues_copy, ","))
  }

  
  return_table <- data.frame(matrix("NA", nrow=length(tissues_individual), ncol=3), stringsAsFactors = FALSE)
  colnames(return_table) <- c("Tissues", "miR", "num_miR")
  return_table$Tissues <- tissues_individual
  #return_table$miR <- return_miR
  
  for (y in 1:length(tissues_individual)){
    col_ind <- which(colnames(confidence) == tissues_individual[y])
    return_table$miR[y] <- paste(levels(confidence$miR[which(confidence[,col_ind] == col_ind)]), collapse = ", ")
    return_table$num_miR[y] <- as.character(length(which(confidence[,col_ind] == 1)))
  }
  return_table
}

# ============= TASK 2 ============

give_mir_return_tissues <- function(mirs, confidence_level){
  # mirs: character string of user mirs
  # Returns all tissues that these mirs are in (1), a table

  
  if (confidence_level == "Very High"){
    confidence <- read.table("/Users/helending/Documents/BCB330/mir_expr_table/data/combined_all_expr_tables_0.99.csv")
  } else if (confidence_level == "High"){
    confidence <- read.table("/Users/helending/Documents/BCB330/mir_expr_table/data/combined_all_expr_tables_0.95.csv")
  } else if (confidence_level == "Medium"){
    confidence <- read.table("/Users/helending/Documents/BCB330/mir_expr_table/data/combined_all_expr_tables_2_3.csv")
  } else if (confidence_level == "Low"){
    confidence <- read.table("/Users/helending/Documents/BCB330/mir_expr_table/data/combined_all_expr_tables_1_3.csv")
  }
  
  confidence <- data.frame(confidence)
  mirs_copy <- mirs
  
  if (regexpr(",", mirs_copy, perl = TRUE) == -1){
    tissues <- as.character()
    mirs_individual <- mirs_copy
    return_table <- confidence[which(confidence$miR == mirs_individual[1]),]
    tissues <- paste(colnames(return_table)[which(return_table[1,] == 1)], tissues, collapse=",") 
    return_table$All_tissues <- tissues
    return_table_copy <- return_table[,which(return_table[1,] == 1)]
    return_table_copy <- cbind(return_table$miR, return_table_copy)
    return_table_copy <- cbind(return_table_copy, "Tissues"=return_table$All_tissues)
    return_table_copy
  } else{
    mirs_individual <- unlist(strsplit(mirs_copy, ","))
    ind_contain_mir <- which(confidence$miR == mirs_individual[1])
    ind <- ind_contain_mir
    for (i in 2:length(mirs_individual)){
      ind_contain_mir <- which(confidence$miR == mirs_individual[i])
      ind <- c(ind, ind_contain_mir)
      return_table <- confidence[ind,]
    }
    
    return_table$All_tissues <- "NA"

    for (k in 1:nrow(return_table)){
      tissues <- as.character()
      tissues <- paste(colnames(return_table)[which(return_table[k,] == 1)], tissues, collapse=",") 
      return_table$All_tissues[k] <- tissues
    }
    rownames(return_table) <- 1:nrow(return_table)
    return_table_copy <- return_table[,which(return_table[1,] == 1)]
    return_table_copy <- cbind(return_table$miR, return_table_copy)
    return_table_copy <- cbind(return_table_copy, "Tissues"=return_table$All_tissues)
    colnames(return_table_copy)[1] <- "miR"
    return_table_copy
  }
  
}


# ============= TASK 3 ============

give_tissue_mir_return_info <- function(tissues, mirs, confidence_level){
  # tissues: character string of tissues
  # mirs: character string of mirs
  # Returns a table of information (a subset of unique_mirdip_all_one_copy)
  
  
  if (confidence_level == "Very High"){
    confidence <- read.table("/Users/helending/Documents/BCB330/mir_expr_table/data/combined_all_expr_tables_0.99.csv")
  } else if (confidence_level == "High"){
    confidence <- read.table("/Users/helending/Documents/BCB330/mir_expr_table/data/combined_all_expr_tables_0.95.csv")
  } else if (confidence_level == "Medium"){
    confidence <- read.table("/Users/helending/Documents/BCB330/mir_expr_table/data/combined_all_expr_tables_2_3.csv")
  } else if (confidence_level == "Low"){
    confidence <- read.table("/Users/helending/Documents/BCB330/mir_expr_table/data/combined_all_expr_tables_1_3.csv")
  }
  confidence <- data.frame(confidence)
  
  
  tissues_copy <- tissues
  if (regexpr(",", tissues_copy, perl=TRUE) == -1){
    tissues_individual <- tissues_copy
  }else{
    tissues_individual <- unlist(strsplit(tissues_copy, ","))
  }
  
  mirs_copy <- mirs
  if (regexpr(",", mirs_copy, perl=TRUE) == -1){
    mirs_individual <- mirs_copy
  }else{
    mirs_individual <- unlist(strsplit(mirs_copy, ","))
  }
  
  # select mir 
  ind <- which(confidence$miR == mirs_individual[1])
  ind_all <- ind
  for (i in 2:length(mirs_individual)){
    ind <- which(confidence$miR == mirs_individual[i])
    ind_all <- c(ind_all, ind)
  }
  
  selected_mir <- confidence[ind_all,]
  rownames(selected_mir) <- 1:nrow(selected_mir)
  
  # select tissue
  ind_t <- which(colnames(confidence) == tissues_individual[1])
  ind_t_all <- ind_t
  for (i in 2:length(tissues_individual)){
    ind <- which(colnames(confidence)== tissues_individual[i])
    ind_t_all <- c(ind_t_all, ind)
  }
  
  selected_both <- selected_mir[,c(1, ind_t_all)]
  selected_both
}
