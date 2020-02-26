authors: Tobias Guldberg Frøslev & Susana Santos

#### Initial steps on server - using scripts from other repositories on this github
#Demultiplexing
#    dada2_demultiplexing.sh

# DADA2 processing
#    R < dada2_sickle.r --no-save &>dada2log.txt

#Empty samples gave some problems. These were deleted manually from the filtered/matched folders, and the last part of the scripts were run alone

# ITSx extraction
#    itsx_fungi.parallel.sh

# Clustering
#    Cluster_fungi_vsearch_985.sh

# Remove abundance information from fasta and otutable
#    sed -i 's/;size.*//g' Fungi.centroids
#    sed -i 's/;size.*;\t/\t/g' Fungi.otutable

#make matchlist for LULU curation and get fungal IDs by matching with UNITE.
# bash ~/data/Biowide/scripts/BW_batch_matchlists.sh 80 84
# bash get_id_vsearch.sh Fungi.centroids

#move filed to laptop, and do the rest in R


library(devtools)
library(lulu)
library(readr)
library(stringr)
library(dplyr)
library(here)

####################################################################################################################################
##########################################               ######################################################
########################################## LULU curation ######################################################
##########################################               ######################################################
####################################################################################################################################

# Removes erroneous molecular operational taxonomic units (OTUs) by combining sequence similarity and co-occurrence patterns

tabL <- read.csv("Fungi.otutable",sep='\t',header=T,as.is=TRUE,row.names = 1)
mlL <- read.csv("Fungi.centroids.matchlist",sep='\t',header=F,as.is=TRUE)
proc_min <- lulu(tabL,mlL, minimum_match = 84, minimum_relative_cooccurence = 1)
#Number of OTUs retained:
proc_min$curated_count
##  7882
#Number of OTUs discarded
proc_min$discarded_count
##  914
saveRDS(proc_min,"LULU_fungiRDS") 

####################################################################################################################################
########################################## attach taxonomic IDs    ######################################################
########################################## adjust for match rate   ######################################################
########################################## Produce final OTU table ######################################################
####################################################################################################################################

edna_tax <- NULL
edna_tax <- read_tsv("Fungi.centroids.hits",col_names=F) # vsearch version
#taxonomy$X1 <- gsub(";size=.*","",taxonomy$X1)

names(edna_tax) <- c("OTUid","pident","SH_header")
SH_split <- as.data.frame(str_split_fixed(edna_tax$SH_header, "\\|", 5))
names(SH_split) <- c("Taxon","Accession","SH","dataset","taxonstring")

#sum(is.na(SH_split))

SH_split$taxonstring <- gsub(".__","",SH_split$taxonstring)
Taxa <- as.data.frame(str_split_fixed(SH_split$taxonstring, ";", 7))
names(Taxa) <- c("kingdom","phylum","class","order","family","genus","species")

sum(is.na(Taxa))

edna_tax <- cbind(edna_tax[,1:2],SH_split[,1:4])
edna_tax <- cbind(edna_tax,Taxa)

sum(is.na(edna_tax))

edna_tax[,-2] <- lapply(edna_tax[,-2], as.character)
sum(is.na(edna_tax))

#If you want to throw out non-fungal sequences. Use with care...
#edna_tax <- edna_tax[!is.na(edna_tax$kingdom) & edna_tax$kingdom == "Fungi",]
#edna_tax <- edna_tax[edna_tax$pident > 70,]

cutoff_index <- which(edna_tax$pident < 75 & edna_tax$phylum != "unidentified")
reclassify_index <- c("class","order","family","genus")
edna_tax[cutoff_index,reclassify_index] <- "unidentified"
edna_tax$species[cutoff_index] <- paste0(edna_tax$phylum[cutoff_index], "_sp")

cutoff_index <- which(edna_tax$pident < 80 & edna_tax$pident >= 75 & edna_tax$class != "unidentified")
reclassify_index <- c("order","family","genus")
edna_tax[cutoff_index,reclassify_index] <- "unidentified"
edna_tax$species[cutoff_index] <- paste0(edna_tax$class[cutoff_index], "_sp")

cutoff_index <- which(edna_tax$pident < 85 & edna_tax$pident >= 80 & edna_tax$order != "unidentified")
reclassify_index <- c("family","genus")
edna_tax[cutoff_index,reclassify_index] <- "unidentified"
edna_tax$species[cutoff_index] <- paste0(edna_tax$order[cutoff_index],"_sp")

cutoff_index <- which(edna_tax$pident < 90& edna_tax$pident >= 85 & edna_tax$family != "unidentified")
reclassify_index <- c("genus")
edna_tax$genus[cutoff_index] <- "unidentified"
edna_tax$species[cutoff_index] <- paste0(edna_tax$family[cutoff_index], "_sp")

cutoff_index <- which(edna_tax$pident < 98 & edna_tax$pident >= 90 & edna_tax$genus != "unidentified")

edna_tax$species[cutoff_index] <- paste0(edna_tax$genus[cutoff_index], "_sp")

#adjust classification to match rate: 98, 90, 85, 80, and 75% sequence identity as a criterion for assigning OTUs with names of a species, genus, family, order, or class, respectively.
for(i in 1:nrow(edna_tax)){
  if (edna_tax[i,"genus"] == "unidentified"){
    if (edna_tax[i,"family"] != "unidentified"){
      edna_tax[i,"species"] <- paste0(edna_tax[i,"family"],"_sp")
    } else if (edna_tax[i,"order"] != "unidentified"){
      edna_tax[i,"species"] <- paste0(edna_tax[i,"order"],"_sp")
    } else if (edna_tax[i,"class"] != "unidentified"){
      edna_tax[i,"species"] <- paste0(edna_tax[i,"class"],"_sp")
    } else if (edna_tax[i,"phylum"] != "unidentified"){
      edna_tax[i,"species"] <- paste0(edna_tax[i,"phylum"],"_sp")
    } else if (edna_tax[i,"kingdom"] != "unidentified"){
      edna_tax[i,"species"] <- paste0(edna_tax[i,"kingdom"],"_sp")
    }
  }
}


#JOIN Taxonomic info and OTU table
proc_min <- readRDS(here::here("LULU_fungiRDS")) 
tabL$OTUid <- row.names(tabL)
MiseqRun2_fungal_table <- left_join(tabL,edna_tax, by = "OTUid")
write_tsv(MiseqRun2_fungal_table,here::here("MiseqRun2_fungal_table.txt"))

#Join sequences to table
allcentroids <- read.csv("Fungi.centroids", sep="\t", header=F, as.is = TRUE)

otusID <- seq(1,length(allcentroids$V1),2)
seqsID <- seq(2,length(allcentroids$V1),2)
otus <- allcentroids[otusID,]
seqs <- allcentroids[seqsID,]
otus <- gsub(">","",otus)
centroid_df <- data.frame(centroid = otus, sequence = seqs)

MiseqRun2_fungal_table2 <- left_join(MiseqRun2_fungal_table,centroid_df, by = c("OTUid" = "centroid"))

#Add curation "advice" from LULU
MiseqRun2_fungal_table2$lulu_result <- "real_otu"
MiseqRun2_fungal_table2[MiseqRun2_fungal_table2$OTUid %in% proc_min$discarded_otus,]$lulu_result <- "error"
write_tsv(MiseqRun2_fungal_table2,here::here("MiseqRun2_fungal_table.txt"))

