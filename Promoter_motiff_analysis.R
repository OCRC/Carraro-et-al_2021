# Analyzing FASTA sequences obtained from Nhlh2 targets promoters

setwd("~/Documents/Bioinfo/Carraro/test_genes.fasta_SPLIT FILES")

list_of_files <- list.files(path = ".", recursive = TRUE,
                            pattern = "\\.txt$", 
                            full.names = TRUE)

#install.packages("data.table", repos = "https://cran.rstudio.com")
library(data.table)

list_of_transcripts <- list()

for( i in (1:length(list_of_files))){
  list_of_transcripts[i] <- fread(list_of_files[[i]], header = F)
}

promoter_lists <- list()
transcript_names <- list()
for( i in (1:length(list_of_transcripts))){
  transcript <- unlist(list_of_transcripts[i])
  transcript_name <- gsub(" .*","", transcript[1])
  transcript_promoter_seq <- paste(unlist(transcript[-1]), collapse="")
  promoter_lists[i] <- transcript_promoter_seq
  transcript_names[i] <- transcript_name
}

transcript_names_cleaned <- gsub(".*_","",unlist(transcript_names))
promoter_lists_vector <- unlist(promoter_lists)

names(promoter_lists) <- transcript_names_cleaned
names(promoter_lists_vector) <- transcript_names_cleaned

# Now that we've handled this messy data

motiff1 <- 'CAGCTG'
motiff2 <- 'CACGTG'
motiff3 <- 'CAAATG'
motiff4 <- 'CATGTG'
motiff5 <- 'CACATG'

has_motiff_CAGCTG <- grepl(motiff1, promoter_lists_vector, fixed = T)
has_motiff_CACGTG <- grepl(motiff2, promoter_lists_vector, fixed = T)
has_motiff_CAAATG <- grepl(motiff3, promoter_lists_vector, fixed = T)
has_motiff_CATGTG <- grepl(motiff4, promoter_lists_vector, fixed = T)
has_motiff_CACATG <- grepl(motiff5, promoter_lists_vector, fixed = T)


setwd("~/Documents/Bioinfo/Carraro")
convertion_ids <- read.table('id_conversion.txt', row.names = 1)

names(promoter_lists_vector) <- convertion_ids$V4
names(has_motiff_CAGCTG) <- convertion_ids$V4
names(has_motiff_CACGTG) <- convertion_ids$V4
names(has_motiff_CAAATG) <- convertion_ids$V4
names(has_motiff_CATGTG) <- convertion_ids$V4
names(has_motiff_CACATG) <- convertion_ids$V4

has_motiff_CAGCTG_int <- as.integer(has_motiff_CAGCTG)
has_motiff_CACGTG_int <- as.integer(has_motiff_CACGTG)
has_motiff_CAAATG_int <- as.integer(has_motiff_CAAATG)
has_motiff_CATGTG_int <- as.integer(has_motiff_CATGTG)
has_motiff_CACATG_int <- as.integer(has_motiff_CACATG)

names(has_motiff_CAGCTG_int) <- convertion_ids$V4
names(has_motiff_CACGTG_int) <- convertion_ids$V4
names(has_motiff_CAAATG_int) <- convertion_ids$V4
names(has_motiff_CATGTG_int) <- convertion_ids$V4
names(has_motiff_CACATG_int) <- convertion_ids$V4

genes <- unique(convertion_ids$V4)

motiff_matrix_int <- cbind(has_motiff_CAGCTG_int,
                           has_motiff_CACGTG_int,
                           has_motiff_CAAATG_int,
                           has_motiff_CATGTG_int,
                           has_motiff_CACATG_int)

has_motiff <- rowSums(motiff_matrix_int)

has_any_motiff <- has_motiff >= 1
has_any_motiff_int <- as.integer(has_any_motiff)
names(has_any_motiff_int) <- convertion_ids$V4

write.csv(motiff_matrix_int, file = 'Has_Motiff_Matrix.csv', row.names = T)
write.csv(has_any_motiff_int, file = 'Has_Any_Motiff_Matrix.csv', row.names = convertion_ids$V4)

names(has_any_motiff_int) == 'Pomc'

sum(has_any_motiff_int[names(has_any_motiff_int) == 'Pomc']) > 0


targets <- read.table('Nhlh2_targets_1.txt')
target_genes <- targets$V1
has_promoter_motiff <- list()
how_many_promoters <- list()
for(i in 1:length(target_genes)){
  if ( sum(has_any_motiff_int[names(has_any_motiff_int) == target_genes[i]]) > 0 ){
    has_promoter_motiff[i] <- T
  } else {
    has_promoter_motiff[i] <- F
  }
  how_many_promoters[i] <- sum(has_any_motiff_int[names(has_any_motiff_int) == target_genes[i]])
}

has_promoter_motiff <- unlist(has_promoter_motiff)
names(has_promoter_motiff) <- target_genes
write.csv(as.integer(has_promoter_motiff), file = 'Gene_Has_Any_Promoter_Motiff_Matrix.csv', row.names = target_genes)
write.csv(has_promoter_motiff, file = 'Gene_Has_Any_Promoter_Motiff_Matrix_Boolean.csv', row.names = target_genes)

how_many_promoters <- unlist(how_many_promoters)
names(how_many_promoters) <- target_genes
write.csv(how_many_promoters, file = 'How_Many_Transcripts_Have_Promoter_Motiffs.csv', row.names = target_genes)


