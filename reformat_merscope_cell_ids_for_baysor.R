library(tidyverse)

file <- "detected_transcripts.csv"
transcripts <- read_csv(file, col_types = list(cell_id=col_character()) )

#for baysor, transcripts without a cell ID assigned have to be indicated with a  zero
transcripts$cell_id [ transcripts$cell_id == "-1"] <- "0"

write.csv(transcripts, file="detected_transcripts_for_baysor.csv", row.names = FALSE)
