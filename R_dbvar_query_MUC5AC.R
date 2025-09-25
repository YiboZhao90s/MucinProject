#install.packages("rentrez")
library(rentrez)
search_res <- entrez_search(db = "dbvar",
                            term = "MUC5AC",
                            retmax = 1000,
                             use_history = T)
batch_size <- 100   # fetch 100 at a time
all_records <- list()

for(start in seq(0, search_res$count-1, by=batch_size)) {
  message("Fetching records ", start+1, " to ", min(start+batch_size, search_res$count))
  
  records <- entrez_fetch(
    db = "dbvar",
    web_history = search_res$web_history,
    rettype = "docsum",
    retmode = "xml",
    retstart = start,
    retmax = batch_size
  )
  all_records[[length(all_records)+1]] <- records
}
all_rows <- list()  # to store rows from all batches

for(records in all_records) {
  xml_file <- read_xml(records)
  docs <- xml_find_all(xml_file, ".//DocumentSummary")
  
  # Extract info for each DocumentSummary
  rows <- lapply(docs, function(doc) {
    uid <- xml_attr(doc, "uid")
    ST <- xml_text(xml_find_first(doc, "./ST"))
    SV <- xml_text(xml_find_first(doc, "./SV"))
    Variant_type <- xml_text(xml_find_first(doc, ".//dbVarVariantTypeList/string"))
    Genes <- xml_text(xml_find_all(doc, ".//dbVarGeneList/dbVarGene/name")) %>% paste(collapse=",")
    
    # Extract placements
    placements <- xml_find_all(doc, ".//dbVarPlacement")
    if(length(placements) == 0) {
      data.frame(uid, ST, SV, Variant_type, Genes,
                 Chr=NA, Start=NA, End=NA, Assembly=NA, stringsAsFactors = FALSE)
    } else {
      do.call(rbind, lapply(placements, function(p) {
        Chr <- xml_text(xml_find_first(p, "./Chr"))
        Start <- xml_text(xml_find_first(p, "./Chr_start"))
        End <- xml_text(xml_find_first(p, "./Chr_end"))
        Assembly <- xml_text(xml_find_first(p, "./Assembly"))
        data.frame(uid, ST, SV, Variant_type, Genes, Chr, Start, End, Assembly, stringsAsFactors = FALSE)
      }))
    }
  })
  
  all_rows <- c(all_rows, rows)
}


# Combine all batches into a single dataframe
df <- bind_rows(all_rows)

head(df)

keep <- ifelse(substr(df$Assembly,1,6)=="GRCh38",1,0)
table(keep)
df <- df[keep==1,]
df <- df[!duplicated(df$uid),]
df <- df[df$Genes=="MUC5AC",]
str(df)
df$chr <- as.numeric(df$Chr)
df$Start <- as.numeric(df$Start)
df$End <- as.numeric(df$End)
df$length <- df$End-df$Start+1
df <- df[!is.na(df$Chr)&!is.na(df$Start)&!is.na(df$length),]
table(df$Variant_type)
entrez_db_summary("dbvar")

summary(df$length>=50)
library(ggplot2)
ggplot(df, aes(x = length)) +
  geom_histogram(binwidth = 2000, fill = "grey", color = "black")+
  labs(
    title = "Distribution of Variant Lengths",
    x = "Variant Length (bp)",
    y = "Count"
  ) +
  theme_minimal()
write.csv(df, "SV_MUC5AC_dbvar.csv", row.names = F, quote = F)
