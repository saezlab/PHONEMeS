library(OmnipathR)
library(magrittr)
library(dplyr)

### Construct kinase-substrate interaction network
omnipath_ptm <- get_signed_ptms()
omnipath_ptm <- omnipath_ptm[omnipath_ptm$modification %in% c("dephosphorylation", "phosphorylation"), ]

# select target (substrate_genesymbol) and source (enzyme_genesymbol)
KSN <- omnipath_ptm[, c(4, 3)]

# add phosphorylation site to target
KSN$substrate_genesymbol <- paste(KSN$substrate_genesymbol, omnipath_ptm$residue_type, sep = "_")
KSN$substrate_genesymbol <- paste(KSN$substrate_genesymbol, omnipath_ptm$residue_offset, sep = "")

# set direction and likelihood of interaction
KSN$mor <- 1
KSN$likelihood <- 1

# we remove duplicated edges
KSN$id <- paste(KSN$substrate_genesymbol, KSN$enzyme_genesymbol, sep = "")
KSN <- KSN[!(duplicated(KSN$id)), ]
KSN <- KSN[, -5]

# rename KSN to fit decoupler format
names(KSN)[1:3] <- c("target", "source", "interaction")
KSN <- KSN[c("source", "interaction", "target")]

# re-attach phosphosites
phospho_prots <- data.frame(KSN$target)
names(phospho_prots) <- "source"
phospho_prots$interaction <- 1
phospho_prots$target <- gsub("_.*","",phospho_prots$source)


phonemesKSN <- rbind(KSN, phospho_prots)
rm(KSN, phospho_prots, omnipath_ptm)
