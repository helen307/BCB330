library(httr)

# const values
url <- "http://ophid.utoronto.ca/mirDIP"

mapScore <- list("0", "1", "2", "3");
names(mapScore) <- c("Very High", "High", "Medium", "Low")


unidirectionalSearchOnGenes <- function(geneSymbols, minimumScore) {
  
  parameters <- list(
    genesymbol = geneSymbols,
    microrna = "",
    scoreClass = mapScore[minimumScore]
  )
  
  # ... send http POST
  res <- POST(paste(url, "/Http_U", sep = ""), body = parameters, encode = "form", verbose())
}

unidirectionalSearchOnMicroRNAs <- function(microRNAs, minimumScore) {
  
  parameters <- list(
    genesymbol = "",
    microrna = microRNAs,
    scoreClass = mapScore[minimumScore]
  )
  
  # ... send http POST
  res <- POST(paste(url, "/Http_U", sep = ""), body = parameters, encode = "form", verbose())
}

# make results-map as keyword - value
makeMap <- function(res) {
  
  ENTRY_DEL = "\001"
  KEY_DEL = "\002"
  
  response = content(res, "text")
  
  arr = unlist(strsplit(response, ENTRY_DEL, fixed = TRUE))
  
  list_map <- list("")
  vec_map_names <- c("");
  
  for (str in arr) {
    arrKeyValue = unlist(strsplit(str, KEY_DEL, fixed = TRUE));
    
    if (length(arrKeyValue) > 1) {
      list_map[length(list_map) + 1] <- arrKeyValue[2]
      vec_map_names[length(vec_map_names) + 1] <- arrKeyValue[1]
    }
  }
  
  names(list_map) <- vec_map_names
  
  list_map
}

get_mirDIP_microrna <- function(gene_string, userScore){
  # gene_string: string of gene names, comma delimitaed
  # userScore: one of 'Very High', 'High', 'Medium', 'Low'
  # Return: miR associated with the genes
  
  geneSymbols = gene_string
  minimumScore = userScore
  res <- unidirectionalSearchOnGenes(geneSymbols, minimumScore)
  responseCode = status_code(res)
  
  if (responseCode != 200) {
    cat("Error: Response Code : ", responseCode, "\r\n")
  } else {
    list_map <- makeMap(res)
    write(unlist(list_map$results), file="test3.txt", sep="\r\n")
    x <- read.delim("test3.txt")
    miR_obtained <- as.character(x$MicroRNA)
    miR_obtained
  }
}


query_for_mir_from_mirdip <- function(confidence){
  # confidence: string, only from: "Very High", "High", "Medium", "Low", mind spelling!
  # Update df: fills in the miR column
  
  last_ind = 0
  for (j in 1:nrow(gene_interaction)){
    a <- get_mirDIP_microrna(as.character(gene_interaction$SYMBOL[j]), userScore=confidence)
    
    # bind new miR to the first column
    if (length(a) != 0){
      df$miR[last_ind:(last_ind+length(a))] <- a
      
      # find which tissues each interaction has a 1, 
      # and fill that df cell accordingly, only 0/1
      for (p in 2:ncol(gene_interaction)){
        if (gene_interaction[j, p] == 1){
          df[last_ind:(last_ind+length(a)),p] <- 1
        }
        else if(gene_interaction[j, p] == 0){ # == 0
          df[last_ind:(last_ind+length(a)),p] <- 0
        } else { # == ?
          df[last_ind:(last_ind+length(a)),p] <- "?"
        }
      }
      last_ind = last_ind + length(a) + 1
    } else {
      j = j + 1
    }
  }
}


query_mir_all_confidence <- function(confidence){
  # confidence: one of "Very High", "High", "Medium", "Low"
  # Write to file in: "mirdip_very_high.csv", "mirdip_high.csv", "mirdip_medium.csv", "mirdip_low.csv"
  
  confidence_level <- confidence
  # how many miR will I finally get?
  all_genes <- paste0(gene_interaction_dup_approved$SYMBOLS, collapse = ', ')
  result <- get_mirDIP_microrna(all_genes, confidence) 
  
  # save df in various file names
  if (length(grep(" ", confidence_level)) == 0){
    name <- confidence_level
  }else{
    name <- unlist(strsplit(confidence_level, " "))
    name <- paste(name[1], name[2], sep="_")
  }
  filename_mod <- paste("mirdip", name, sep="_")
  filename_mod <- paste(filename_mod, ".csv", sep="")
  result <- result[order(result)]
  write.table(result, filename_mod)
}

get_mirDIP_gene <-function(microrna_string, userScore){
  
  microRNAs = microrna_string

  minimumScore = userScore
  
  res <- unidirectionalSearchOnMicroRNAs(microRNAs, minimumScore)
  responseCode = status_code(res)
  if (responseCode != 200) {
    cat("Error: Response Code : ", responseCode, "\r\n")
  } else {
    list_map <- makeMap(res)
    write(unlist(list_map$results), file="test4.txt", sep="\r\n")
    x <- read.delim("test4.txt")
  }
}

