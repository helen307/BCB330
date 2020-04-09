library(limma)
library(ggplot2)
library(miRNAmeConverter)


read_data_by_limma <- function(norm_ind){
  # norm_ind: integer
  # Returns data that is read by limma
  
  read_data <- limma::read.maimages(files=SDRF[norm_ind, "Array Data File"], 
                                    source="agilent", 
                                    green.only=TRUE, 
                                    other.columns = "gIsWellAboveBG")
  read_data
}


plot_negctrl_vs_rest <- function(read_data){
  # read_data: Agilent data read by limma, type: dataframe
  # ind: the index that the miR column is 
  # Return a plot that result by only picking negative control versus the rest.
  ind <- 2
  neg_ctrl_ind <- which(read_data$SystematicName == "NegativeControl")
  read_data <- cbind("LABEL" = 1:nrow(read_data), read_data)
  read_data$LABEL[neg_ctrl_ind] <- "negControl"
  neg_ctrl <- read_data[neg_ctrl_ind, ]
  rows <- nrow(neg_ctrl) # 288
  
  # the rest
  rest <- read_data[-neg_ctrl_ind,]
  read_data$LABEL[which(read_data$LABEL != "negControl")] <- "rest"
  
  # Build plot
  plot <- boxplot(read_data[, ind + 1] ~ read_data$LABEL, ylab="Expression", xlab="labels")
  for (i in (ind + 2):ncol(read_data)){
    plot <- boxplot(read_data[, i] ~ read_data$LABEL, ylab="Expression", xlab="labels", add=TRUE)
  }
  plot
}

plot_negctrl_plus_nonhuman_vs_rest <- function(read_data){
  # read_data: Agilent data read by limma, type: dataframe
  # ind: ind for neg_vs_human graph
  # Return a plot that result by only picking negative control versus the rest.
  ind <- 2
  # exclude human
  human <- grep(pattern="hsa", read_data$SystematicName)
  length(human) # 11569
  read_data <- cbind("LABEL_OTHERS"=1:nrow(read_data), read_data)
  read_data$LABEL_OTHERS[human] <- "human"
  
  # exclude positive controls
  pos_control <- which(read_data$ControlType == 1)
  length(pos_control) # 281
  
  # only non-human and neg control
  neg_ctrl_with_other_spec <- read_data[-c(human, pos_control), ]
  read_data$LABEL_OTHERS[which(read_data$LABEL_OTHERS != "human")] <- "rest (negControl+non-human)"
  
  plot2 <- boxplot(read_data[, ind + 2] ~ read_data$LABEL_OTHERS, col=c("red", "blue"), ylab="Expression", xlab="labels")
  for (i in (ind + 3):ncol(read_data)){
    plot2 <- boxplot(read_data[, i] ~ read_data$LABEL_OTHERS, col=c("red", "blue"), ylab="Expression", xlab="labels", add=TRUE)
  }
  plot2
}

data_cleaning <- function(read_data){
  # read_data: Agilent data read by limma, type: dataframe
  # Return a report (# human miR, corner removal report, picking negative control report)
  
  # find human miRNA
  human <- grep(pattern="hsa", read_data$SystematicName)
  human_miR <- read_data[human,]
  cat("number of hsa miRNAs: ", nrow(read_data), "\n")
  
  # Remove corners
  remove <- grep(pattern="Corner", read_data$SystematicName)
  read_data <- read_data[-remove, ]
  cat("number of corners removed: ", length(remove), "\n")
  cat("number of miR left after removal: ", nrow(read_data), "\n")
  
  # Pick negative control
  y_neg_ctrl <- read_data[which(read_data$SystematicName == "NegativeControl"),]
  y_neg_ctrl <- data.frame(y_neg_ctrl)
  cat("number of negative control selected: ", nrow(y_neg_ctrl))
  y_neg_ctrl
}


nlargest <- function(m, n) {
  # Return n largest values and position for matrix m
  
  copy_m <- m
  copy_m <- as.vector(copy_m)
  res <- order(copy_m, decreasing=TRUE)[1:n]
  pos <- arrayInd(res, dim(m), useNames = TRUE)
  list(values = m[res],position = pos)
}

define_expression_threshold <- function(y_neg_ctrl, norm_ind_len){
  # top 1% 
  rows <- nrow(y_neg_ctrl)
  cols <- ncol(y_neg_ctrl)
  n <- floor(rows * (cols-1) * 0.01) # 9
  cat("number of negative control to select: ", n, "\n")
  
  y_neg_ctrl <- as.matrix(y_neg_ctrl[,2:cols])
  cat("number of samples consistent: ", (ncol(y_neg_ctrl) == norm_ind_len), "\n")
  
  # inspect the top 1% values
  greatest_0.01 <- nlargest(y_neg_ctrl, n)
  
  # sd
  sd(greatest_0.01$position[,2]) # 1.922094
  cat("SD is", sd(greatest_0.01$position[,2]), "in", ncol(y_neg_ctrl), "samples", "\n")
  
  # thres
  thres <- min(greatest_0.01$values) # 54
  cat("Threshold: ", thres, "\n")
  
  greatest_0.01
}


define_miRNA_expression <- function(human_miR_table, thres){
  
  # read in human mirna expression data
  copy_y <- human_miR_table[, 2:ncol(human_miR_table)]
  
  # compare each expression value to thres
  for(i in 1:nrow(copy_y)){
    for (j in 1:ncol(copy_y)){
      if (copy_y[i, j] >= thres){
        copy_y[i, j] <- 1
      } else{
        copy_y[i, j] <- 0
      }
    }
  }
  
  # sum the rows
  copy_y <- cbind("SUM"=rowSums(copy_y), copy_y)
  copy_y <- cbind("miR"=human_miR_table$SystematicName, copy_y)
  rownames(copy_y) <- 1:nrow(copy_y)
  copy_y
}

if_ctrl_all_zero <- function(human_miR_table, copy_y){
  # human_miR_table: human miR 
  # Returns a report.
  
  copy_y <- data.frame(copy_y)
  human_miR_table <- data.frame(human_miR_table)
  
  # Select negative control
  neg_control <- copy_y$SUM[which(copy_y$miR == "NegativeControl")]
  cat("expression value of negative control are: ", length(neg_control), "\n")
  
  # Select positive control
  pos_control <- copy_y$SUM[which(copy_y$ControlType == 1)]
  cat("expression value of positive control are: ", length(pos_control), "\n")
  
  copy_y
}

process_human_miR <- function(copy_y, norm_ind){
  # Select only human miR
  copy_y <- data.frame(copy_y)
  human_mir <- copy_y[grep(pattern="hsa", copy_y$miR),]
  rownames(human_mir) <- 1:nrow(human_mir)
  
  cat("number of unique human miR: ", length(unique(human_mir$miR)), 
      ", but total number of human miR: ", length(human_mir$miR), "\n")
  
  unique <- data.frame("miR"= human_mir$miR, 
                       "SUM" = human_mir$SUM)
  
  unique <- with(unique, aggregate(list(SUM = SUM, num_dup=rep(1,nrow(unique))), list(miR = tolower(miR)), sum))
  unique <- cbind(unique, "avg" = floor(unique$SUM/unique$num_dup))
  unique_avg <- unique[, c(1, 4)]
  unique_avg <- cbind(unique_avg, "num_samples"=length(norm_ind))
  unique_avg
}

save_known_and_novel <- function(unique_all_table, tissue_name){
  unique_all_table <- data.frame(unique_all_table)
  # check name conversion
  nc = MiRNANameConverter()
  known_miR <- checkMiRNAName(nc, as.character(unique_all_table$miR), verbose = FALSE)
  novel_len <- length(unique_all_table$miR) - length(known_miR)
  cat("Novel miRs: ", novel_len, "\n")
  
  filename_novel <- "/Users/helending/Documents/BCB330/saved_expr_data/expr_data_novel/final_novel_"
  filename_novel <- paste(filename_novel, tissue_name, sep="")
  filename_novel <- paste(filename_novel, ".csv", sep="")
  
  filename_known <- "/Users/helending/Documents/BCB330/saved_expr_data/expr_data_known/final_known_"
  filename_known <- paste(filename_known, tissue_name, sep="")
  filename_known <- paste(filename_known, ".csv", sep="")
  
  # save novel
  if (novel_len > 0){
    `%notin%` <- Negate(`%in%`)
    ind_novel <- which(unique_all_table$miR %notin% known_miR)
    novel <- unique_all_table[ind_novel, ]
    rownames(novel) <- 1:nrow(novel)
    # save
    colnames(novel)[2] <- tissue_name
    write.table(novel, file=filename_novel)
    
    # save known
    known <- unique_all_table[-which(unique_all_table$miR %in% novel$miR),]
    colnames(known)[2] <- tissue_name
    write.table(known, file=filename_known)
    cat("number of known miR: ", nrow(known))
    known
  } else if (novel_len == 0){
    # save known
    colnames(unique_all_table)[2] <- tissue_name
    write.table(unique_all_table, file=filename_known)
    cat("number of known miR: ", nrow(unique_all_table))
    unique_all_table
  }
}
