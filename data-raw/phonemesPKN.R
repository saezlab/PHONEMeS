library(OmnipathR)
library(magrittr)
library(dplyr)
### Construct Protein-protein interaction network
omniR <- import_omnipath_interactions()

omnipath_sd <- omniR %>% dplyr::filter(consensus_direction == 1 &
                                         (consensus_stimulation == 1 |
                                            consensus_inhibition == 1
                                         ))

# changing 0/1 criteria in consensus_stimulation/inhibition to -1/1
omnipath_sd$consensus_stimulation[which(omnipath_sd$consensus_stimulation == 0)] <- -1
omnipath_sd$consensus_inhibition[which(omnipath_sd$consensus_inhibition == 1)] <- -1
omnipath_sd$consensus_inhibition[which(omnipath_sd$consensus_inhibition == 0)] <- 1

# check consistency on consensus sign and select only those in a SIF format
sif <- omnipath_sd[, c("source_genesymbol", "consensus_stimulation", "consensus_inhibition", "target_genesymbol")] %>%
  dplyr::filter(consensus_stimulation == consensus_inhibition) %>%
  unique.data.frame()

sif$consensus_stimulation <- NULL
colnames(sif) <- c("source", "interaction", "target")

# remove special characters (complexes are now separated with __, protein names with a minus are replaces with ___)
sif$source <- gsub("_", "__", sif$source)
sif$source <- gsub("-", "___", sif$source)
sif$target <- gsub("_", "__", sif$target)
sif$target <- gsub("-", "___", sif$target)

### Construct kinase-substrate interaction network
omnipath_ptm <- get_signed_ptms()
omnipath_ptm <- omnipath_ptm[omnipath_ptm$modification %in% c("dephosphorylation", "phosphorylation"), ]

# Filter out ProtMapper
omnipath_ptm_filtered <- omnipath_ptm %>%
  dplyr::filter(!(stringr::str_detect(omnipath_ptm$source, "ProtMapper") & n_resources == 1))

# select target (substrate_genesymbol) and source (enzyme_genesymbol)
KSN <- omnipath_ptm_filtered[, c(4, 3)]

# add phosphorylation site to target
KSN$substrate_genesymbol <- paste(KSN$substrate_genesymbol, omnipath_ptm_filtered$residue_type, sep = "_")
KSN$substrate_genesymbol <- paste(KSN$substrate_genesymbol, omnipath_ptm_filtered$residue_offset, sep = "")

# set direction and likelihood of interaction
KSN$mor <- ifelse(omnipath_ptm_filtered$modification == "phosphorylation", 1, -1)
KSN$likelihood <- 1

# we remove ambiguous modes of regulations
KSN$id <- paste(KSN$substrate_genesymbol, KSN$enzyme_genesymbol, sep = "")
KSN <- KSN[!(duplicated(KSN$id) | duplicated(KSN$id, fromLast = TRUE)), ]
KSN <- KSN[, -5]

# rename KSN to fit decoupler format
names(KSN)[1:3] <- c("target", "source", "interaction")
KSN <- KSN[c("source", "interaction", "target")]

phonemesPKN <- rbind(sif, KSN)
rm(KSN, omnipath_ptm, omnipath_ptm_filtered, omnipath_sd, omniR, sif)
