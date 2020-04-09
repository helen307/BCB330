# Purpose: get data from GEO directly for "diseases"
# Date: 2020-02-19

#  ================== PACKAGES ==================
library(GEOmetadb)


if(!file.exists('GEOmetadb.sqlite')){
  getSQLiteFile()
}

# ================== QUERY ==================
file.info('GEOmetadb.sqlite')
# connect
con <- dbConnect(SQLite(),'GEOmetadb.sqlite')
# query
sql <- paste("SELECT DISTINCT gse.title,gse.gse, gpl.title,",
             " gse.submission_date,",
             " gse.supplementary_file",
             "FROM",
             "  gse JOIN gse_gpl ON gse_gpl.gse=gse.gse",
             "  JOIN gpl ON gse_gpl.gpl=gpl.gpl",
             "WHERE",
             "  gse.submission_date > '2010-01-01' AND",
             "  gse.title LIKE '%colorectal cancer%' AND",
             "  gse.title LIKE '%mirna%' AND",
             "  gpl.organism LIKE '%Homo sapiens%' AND",
             "  gpl.title LIKE '%RNA%' ",
             "  ORDER BY gse.submission_date DESC",sep=" ")

rs <- dbGetQuery(con,sql)
rs$title
rs$gse[]
# break the file names up and just get the actual file name
unlist(lapply(rs$supplementary_file,
              FUN = function(x){x <- unlist(strsplit(x,";")) ;
              x <- x[grep(x,pattern="txt",ignore.case = TRUE)];
              tail(unlist(strsplit(x,"/")),n=1)})) [1:10]

counts_files <- rs$supplementary_file[grep(rs$supplementary_file,
                                           pattern = "count",ignore.case = TRUE)]

rs$title
rs$gse[]
