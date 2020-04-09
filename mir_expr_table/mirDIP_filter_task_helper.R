library(httr)
library(miRNAmeConverter)
library(dplyr)
library(ggplot2)
library(shiny)

define_expressivity <- function(percent, table){
  # n: percentage of expression
  # table: table of miRNA expression
  # Return: table which is a modified table on top of table, added a column: SUM_USER
  
  table_copy <- table
  table_copy <- read.table(table_copy)
  table_copy <- data.frame(table_copy)
  table_copy <- cbind("SUM_USER"=1:nrow(table_copy), table_copy)
  num_sample <- table_copy[1,ncol(table_copy)]
  user_selected_thres <- num_sample * percent
  for (i in 1:nrow(table_copy)){
    if (table_copy$SUM[i] > user_selected_thres){
      table_copy$SUM_USER[i] <- 1
    } else {
      table_copy$SUM_USER[i] <- 0
    }
  }
  table_copy[,c(1,2)] <- table_copy[,c(2, 1)]
  colnames(table_copy)[c(1, 2)] <- colnames(table_copy)[c(2, 1)]
  colnames(table_copy)[2] <- colnames(table_copy)[3]
  
  table_copy[,c(1,2)]
}



call_expressivity <- function(tables, percentage){
  # tables: a string of tables
  # percentage: user defined threshold
  # Return filtered tables for all input tables
  tables_copy <- tables

  if (length(tables_copy) == 1){
    new_table <- define_expressivity(percent=percentage, table=tables_copy[[1]])
    output <- list(new_table)
    
  }else{
    new_table_first <- define_expressivity(percent=percentage, table=tables_copy[[1]])
    output <- list(new_table_first)
    for (i in 2: length(tables_copy)){
      new_table <- define_expressivity(percent=percentage, table=tables_copy[[i]])
      output[i] <- list(new_table)
    }
  }
  
  return(output) 
}


combine_all_expr_tables <- function(new_tables, datalist){
  # new_table: the tables after user has defined them
  # datalist: 
  # Updates the expression_all_tissues table, which starts from an empty table
  
  
  for (i in 1:length(list_of_files)){
    new <- c(new, as.character(new_tables[[i]][,1])) # 18478
  }
  unique_together <- data.frame("unique_mir" = unique(new)) # 2704
  
  expression_all_tissue <- data.frame(matrix("NA", nrow=length(unique_together$unique_mir), 
                                             ncol=ncol(gene_interaction_dup) + 1), 
                                      stringsAsFactors = FALSE)
  
  colnames(expression_all_tissue)[2:ncol(expression_all_tissue)] <- colnames(gene_interaction_dup)
  colnames(expression_all_tissue)[1] <- "miR"
  expression_all_tissue$miR <- as.character(unique_together$unique_mir)
  expression_all_tissue <- expression_all_tissue[,-2]
  
  for(k in 1:length(new_tables)){
    new_table_df <- data.frame(new_tables[[k]])
    row_names <- as.character(colnames(new_table_df)[2])
    
    for (i in 1:nrow(new_table_df)){
      for (j in 1:nrow(expression_all_tissue)){
        if (identical(as.character(new_table_df$miR[i]), expression_all_tissue$miR[j])){
          expression_all_tissue[j, row_names] <- new_table_df[i, row_names]
        }
      }
    }
    # the rest: mark as "?"
    expression_all_tissue[which(expression_all_tissue[, row_names] == "NA"), row_names] <- "?"
  }
  expression_all_tissue
}



expr_filter_mirdip <- function(confidence, expression_all_tissue){
  # confidence: string in "Very High", "High", "Medium", "Low"
  # expression_all_tissue: combined expression data from my tissues
  # Returns tables that are filtered according to different confidence levels
  
  # my tissue expr
  expression_all_tissue_copy <- expression_all_tissue
  #[,c(1, 5, 6, 9, 11, 16, 19, 21)]
  
  # mirdip expr
  confidence_copy <- data.frame(matrix("NA", ncol=ncol(gene_interaction_dup)+1, nrow=nrow(confidence)), stringsAsFactors = FALSE)
  colnames(confidence_copy)[2:ncol(confidence_copy)] <- colnames(gene_interaction_dup)
  colnames(confidence_copy)[1] <- "miR"
  confidence_copy$miR <- confidence$miR
  confidence_copy$miR <- as.character(confidence_copy$miR)
  confidence_copy <- confidence_copy[,-2]
  confidence_copy <- confidence_copy
  
  for (i in 1:nrow(expression_all_tissue_copy)){
    for(j in 1:nrow(confidence_copy)){
      if (identical(confidence_copy$miR[j], expression_all_tissue_copy$miR[i])){
        confidence_copy[j,] <- expression_all_tissue_copy[i,]
      }
    }
  }
  
  for (p in 2:ncol(confidence_copy)){
    confidence_copy[which(confidence_copy[,p] == "NA"), p] <- "not on expr"
  }
  
  confidence_copy
}
