# implementation of the (un)wheighted z-test method for kinase activity inference
# pubmed IDs: 28200105, 23532336


# load R packages
library(tidyverse)
library(parallel)


# source match pipeline
# contains the functions required for the weighted z-test
source("./src/utils/matchTool.R")



z_test <- function(all_substr, kinase_substr){
  
  muT <- mean(all_substr, na.rm = T)
  sdT <- sd(all_substr, na.rm = T)
  
  Sm <- mean(kinase_substr, na.rm = T)
  m <- length(kinase_substr)
  
  z_score <- ((Sm - muT) * sqrt(m))/sdT
  
  p_value <- 2*pnorm(-abs(z_score))
  if(p_value == 0) p_value <- 0.0001
  
  if (z_score < 0) activity <- log10(p_value)
  if (z_score > 0) activity <- -log10(p_value)
  #if (Sm < 0) activity <- log10(p_value)
  #if (Sm > 0) activity <- -log10(p_value)
  
  results <- tibble(n=m, log10P=activity)
  
  return(results)
}


infer_by_kinase <- function(substrates, pwm=NULL, phospho){
  
  # get the substrates of a given kinase that were quantified
  # on a given sample
  substrates <- substrates %>%
    inner_join(phospho, by = c("substrate" = "gene", "position", "residue"))
  
  # remove the kinase substrates from the all phosphosites
  # quantified
  phospho <- phospho %>%
    anti_join(substrates, by = c("gene" = "substrate", "position", "residue"))
  
  # if no substrates were quantified for this kinase return NA
  if(nrow(substrates) == 0){
    p <- tibble(n=0, log10P=NA)
  }else{
    # if there is a PWM for a given kinase use it to weigh the substrates fold-changes
    # by affinity to this kinase
    if(!is.null(pwm)){
      substrates <- substrates %>%
        mutate(score = map(.x=kmer, .f = matchScore, pwm = pwm, look_central_aa = T)) %>%
        mutate(score = map_dbl(.x=score, .f = ~ .x$mss)) %>%
        #mutate(score = replace_na(score, 0)) %>%
        filter(!is.na(score))
      if(nrow(substrates) == 0){
        p <- tibble(n=0, log10P=NA)
      }else{
        substrates <- substrates %>%
          mutate(log2fc = log2fc*score)
        p <- z_test(phospho$log2fc, substrates$log2fc)
      }
    }else{
      #if there is not use the kinase substrates fold-changes "as is"
      p <- z_test(phospho$log2fc, substrates$log2fc)
    }
  }
  return(p)
}


infer_by_sample <- function(phospho, ks_list, weights=F){
  
  # remove NAs from phospho data for a given sample
  phospho <- phospho %>%
    filter(!is.na(log2fc))
  
  # infer the activation of each kinase on a given sample
  
  # if to compute without weights just feed the phosphorylation data
  # plus substrate list to infer_by_kinase
  
  # if to compute with weights feed the phosphorylation data,
  # substrate list and PWM of the kinase to infer_by_kinase
  
  if(!weights){
    ka <- ks_list %>%
      mutate(data2 = map(.x = data, .f=infer_by_kinase, phospho=phospho)) %>%
      select(-data) %>%
      unnest(cols = data2)
  } else{
    ka <- ks_list %>%
      mutate(data2 = map2(.x = data, .y=pwm, .f=infer_by_kinase, phospho=phospho)) %>%
      select(-data, -pwm) %>%
      unnest(cols = data2)
  }
  
  return(ka)
}


infer_withoutW <- function(phospho, ks_list, multi_core=F){
  
  # infer kinase activation across samples
  if(!multi_core){
    inference <- phospho %>%
      mutate(data2 = map(.x = data, .f = infer_by_sample, ks_list = ks_list)) %>%
      select(-data) %>%
      unnest(cols = data2)
  }else{
    cores = detectCores()-1
    #data2 <- mclapply(
    #  X = phospho$data,
    #  FUN = infer_by_sample,
    #  ks_list = ks_list,
    #  mc.cores = cores)
    cl <- makeCluster(cores)
    clusterExport(cl, c("infer_by_kinase", "z_test"))
    clusterEvalQ(cl, library(tidyverse))
    data2 <- parLapply(
      cl = cl,
      X = phospho$data,
      fun = infer_by_sample,
      ks_list = ks_list)
    stopCluster(cl)
    inference <- phospho %>%
      select(-data) %>%
      mutate(data2 = data2) %>%
      unnest(cols = data2)
  }
  return(inference)
}


infer_withW <- function(phospho, ks_list, prot_seqs, kinase_pwm, multi_core=F){
  
  # extract k and kmer length (k*2+1) used to build the kinase PWMs
  k = attr(kinase_pwm, "k")
  kmer_l = attr(kinase_pwm, "kmer_l")
  
  # extract kmers for each phosphosite
  # remove phosphosites with kmer length smaller than kmer_l
  phospho <- phospho %>%
    unnest(cols = data) %>%
    pivot_wider(names_from = "sample", values_from = "log2fc") %>%
    inner_join(prot_seqs[, c("gene_name", "seq")], by = c("gene" = "gene_name")) %>%
    select(gene, position, residue, seq, everything()) %>%
    mutate(kmer = map2_chr(.x=seq, .y=position, .f=getKmers, k=k)) %>%
    select(-seq) %>%
    select(gene, position, residue, kmer, everything()) %>%
    filter(map_dbl(.x=kmer, nchar) == kmer_l) %>%
    pivot_longer(names_to = "sample", values_to = "log2fc", -c(gene, position, residue, kmer)) %>%
    select(sample, gene, everything()) %>%
    group_by(sample) %>%
    nest() %>%
    ungroup()
  
  # ensure all kinases to predict activity are in the table with the kinase PWMs
  # get the PWM for each kinase
  ks_list <- ks_list %>%
    inner_join(kinase_pwm, by = "kinase")
  
  
  # infer kinase activation across samples
  if(!multi_core){
    inference <- phospho %>%
      mutate(data2 = map(.x = data, .f = infer_by_sample, ks_list = ks_list, weights = T)) %>%
      select(-data) %>%
      unnest(cols = data2)
  }else{
    cores = detectCores()-1
    #data2 <- mclapply(
    #  X = phospho$data,
    #  FUN = infer_by_sample,
    #  ks_list = ks_list,
    #  weights = T,
    #  mc.cores = cores)
    cl <- makeCluster(cores)
    clusterExport(cl, c("infer_by_kinase", "z_test", "matchScore", "scoreArray", "revComp", "strReverse", "bestSequence", "worstSequence"))
    clusterEvalQ(cl, library(tidyverse))
    data2 <- parLapply(
      cl = cl,
      X = phospho$data,
      fun = infer_by_sample,
      ks_list = ks_list,
      weights = T)
    stopCluster(cl)
    inference <- phospho %>%
      select(-data) %>%
      mutate(data2 = data2) %>%
      unnest(cols = data2)
  }
  return(inference)
}


infer_kinase_activity <- function(ks_list, phospho, with_weights=F, prot_seqs=NULL, kinase_pwm=NULL, multi_core=F){
  
  if(with_weights & c(is.null(prot_seqs) | is.null(kinase_pwm)))
    stop("with weights you have to supply a data frame/tibble with the protein sequences/kinase PWMs")
  else if(with_weights & c(!is.data.frame(prot_seqs) | !is.data.frame(kinase_pwm)))
    stop("with weights you have to supply a data frame/tibble with the protein sequences/kinase PWMs")
  
  # for each kinase build a tibble with all the substrates, positions and residues
  # store them in a list-column
  ks_list <- ks_list %>%
    nest(data = c(substrate, position, residue))
  
  # for each sample build a tibble with all substrates, positions, residues and log2FC
  # store them in a list-column
  phospho <- phospho %>%
    pivot_longer(names_to = "sample", values_to = "log2fc", -c(position, gene, residue)) %>%
    select(sample, gene, everything()) %>%
    group_by(sample) %>%
    nest() %>%
    ungroup()

  
  # infer kinase activation
  if(!with_weights){
    inferKA <- infer_withoutW(phospho, ks_list, multi_core=multi_core)}
  else{
    
    inferKA <- infer_withW(phospho, ks_list, prot_seqs, kinase_pwm, multi_core=multi_core)}
  
  return(inferKA)
}



# wrapper function for infer_kinase_activity
inferKA_W <- function(ks_lists, phospho, with_weights=F, prot_seqs=NULL, kinase_pwm=NULL, multi_core=F){
  
  inferKA <- ks_lists %>%
    group_by(source_type) %>%
    nest() %>%
    ungroup() %>%
    mutate(
      inference = map(
        .x=data,
        .f=infer_kinase_activity,
        phospho=phospho,
        with_weights=with_weights,
        prot_seqs=prot_seqs,
        kinase_pwm=kinase_pwm,
        multi_core=multi_core)) %>%
    select(-data) %>%
    unnest(cols = inference)
  
  return(inferKA)
}
