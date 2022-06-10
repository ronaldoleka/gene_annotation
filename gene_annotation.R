# Ronaldo Leka
# 2022-06-10
# Ensembl Gene ID to HGNC Gene Symbol mapping
# All 68016 Ensembl Gene IDs were obtained from Ensembl Biomart.
# 62279 of them could be mapped to HGNC Gene Symbol using the ensembldb package and EnsDb.Hsapiens.v86 data.
# The remaining 5737 were annotated using web scrapping of the Ensembl.org webpage for each transcript.
# The information given at the top of the page ("XYZ" part of the "Gene: XYZ") was taken as Gene Symbol.

library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(htmltools)
library(textreadr)
library(rvest)


# Read Ensembl Gene IDs and delete duplicates.
ensembl_ids <- read.table("mart_export.txt", fill = TRUE, header = TRUE)
ensembl_ids <- ensembl_ids[, 1]
ensembl_ids <- unique(ensembl_ids)


# Query database using the genes() function from ensembldb and the EnsDb.Hsapiens.v86 object.
ids <- GeneIdFilter(ensembl_ids)
names <- genes(EnsDb.Hsapiens.v86, filter = ids, columns = "gene_name")
names <- cbind.data.frame(names$gene_id, names$gene_name)

colnames(names) <- c("Ensembl Gene ID", "Gene Symbol")


# Get the unmapped Ensembl Gene IDs.
unmapped_ids <- setdiff(ensembl_ids, names$`names$gene_id`)

# Web scraping of Ensembl.org for each unmapped ID.
# Warning: Takes a long time.
scraped_gene_symbols <- c()

for (i in 1:length(unmapped_ids)) {
  
  info <- read_html(paste0("http://useast.ensembl.org/Homo_sapiens/Gene/Summary?g=", unmapped_ids[i]))
  
  infonodes <- info %>% html_nodes(xpath = "//h1[@class='summary-heading']")
  gene <- strsplit(strsplit(toString(infonodes), "<")[[1]][2], " ")[[1]][length(strsplit(strsplit(toString(infonodes), "<")[[1]][2], " ")[[1]])]
  scraped_gene_symbols <- c(scraped_gene_symbols, gene)
  
  # Optional: Will print the progress of the scraping job.
  print(i)
  # print(gene) # Prints current scraped symbol to make sure the progress is correct.
}

# Bind all mapped IDs into a single data frame.
scraped_gene_symbols.df <- cbind.data.frame(unmapped_ids, scraped_gene_symbols)
colnames(scraped_gene_symbols.df) <- c("Ensembl Gene ID", "Gene Symbol")


mapped_gene_ids <- rbind.data.frame(names, scraped_gene_symbols.df)

mapped_gene_ids <- mapped_gene_ids[order(mapped_gene_ids$`Ensembl Gene ID`), ]

rownames(mapped_gene_ids) <- NULL

# Write out a csv file.
write.csv(mapped_gene_ids, "mapped_gene_ids.csv", row.names = FALSE)

# End of script