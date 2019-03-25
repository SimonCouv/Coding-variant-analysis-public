map_peptides <- function(protein_coding_tx, protein_peptides, get_coding_tx=TRUE, ens_mart_db = NA, ensdb, max_tsl=2){
  #protein_coding_tx: columns 
  # 'ens_gene': ens gene IDs for genes to retrieve
  # (optional) 'tx_id': coding transcripts to analyse, not used if get_coding_tx==TRUE
  
  # see ensembldb vignettes https://bioconductor.org/packages/release/bioc/html/ensembldb.html
  # https://bioconductor.org/packages/release/bioc/vignettes/ensembldb/inst/doc/proteins.html
  # https://bioconductor.org/packages/release/bioc/vignettes/ensembldb/inst/doc/coordinate-mapping.html
  
  
  # ANNOTATION AND FILTERING ---------------------------------------------------------------------------------------------------------------------
  
  # get Ensembl annotations (always) -------------------------------------------------------------
  if (is.na(ens_mart_db)){
    warning("No ens_mart_db provided to get coding transcripts, defaulting to dataset hsapiens_gene_ensembl of Mart ENSEMBL_MART_ENSEMBL at http://jan2019.archive.ensembl.org")
    ens_mart_db <- useDataset("hsapiens_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL", host = "http://jan2019.archive.ensembl.org"))
  }
  protein_coding_tx_anno <- getBM(mart = ens_mart_db,
                                  filters = c("ensembl_gene_id", "chromosome_name"),
                                  values = list(ensembl_gene_id=unique(protein_coding_tx$ens_gene),
                                                chromosome_name=c(1:22, "X", "Y", "MT")),
                                  attributes = c("ensembl_gene_id", "ensembl_transcript_id", "transcript_biotype", "transcript_tsl", "ensembl_peptide_id")) %>%
    mutate(tsl = as.numeric(str_replace(transcript_tsl, pattern = "tsl", replacement = ""))) %>%
    dplyr::filter((transcript_biotype == "protein_coding") & (tsl <= max_tsl)) %>%
    dplyr::select(
      ens_gene = ensembl_gene_id,
      tx_id = ensembl_transcript_id,
      p_id = ensembl_peptide_id,
      tsl
      )
  
  if (nrow(protein_coding_tx_anno) == 0)
    stop("No coding transcipts returned.")
  
  genes_not_returned <- setdiff(protein_coding_tx$ens_gene, protein_coding_tx_anno$ens_gene)
  
  if (length(genes_not_returned) > 0){
    warning(sprintf("No coding transcripts returned for gene IDs %s", paste0(genes_not_returned, collapse = ",")))
  }
  
  
  # build tx_overview  -------------------------------------------------------------
  # get coding_tx
  if (get_coding_tx){
    cat("getting all coding transcripts from Ensembl\n")
    
    if("tx_id" %in% names(protein_coding_tx)) {
      warning("tx_id column ignored. Getting all coding tx because get_coding_tx == TRUE")
    }
    
    tx_overview <- protein_coding_tx_anno
    
  } else { #work with provided tx but filter for coding biotype and tsl > tsl_cutoff
    tx_overview <- inner_join(protein_coding_tx, protein_coding_tx_anno, by=c("ens_gene", "tx_id"))
  }
  
  tx_overview <- distinct(tx_overview)
  
  
  # browser()
  
  # MAPPING ---------------------------------------------------------------------------------------------------------------------
  
  pept_prot_pos <- lost_pepts <- lost_tx <- list()
  for (ens_gene in unique(tx_overview$ens_gene)){
    
    #debug
    # browser()
    # ens_gene <- "ENSG00000267467"
    
    cat(sprintf("\nmapping peptides for %s\n",ens_gene))
    
    pepts <- protein_peptides$Peptide_clean[protein_peptides$ens_gene==ens_gene]
    
    # select all transcipts on which select rs IDs occur
    rs_txs <- unique(tx_overview$tx_id[tx_overview$ens_gene==ens_gene])
    
    prots <- ensembldb::proteins(ensdb, filter = ~ gene_id == ens_gene) %>%
      as_data_frame()
    
    pept_prot_pos[[ens_gene]] <- list()
    
    for (i in 1:length(rs_txs)){  # loop over transcripts
      tx <- rs_txs[i]
      cat(sprintf(" - %s\n", tx))
      
      #debug
      # tx <- "ENST00000592954"
      
      pept_prot_pos[[ens_gene]][[tx]] <- list()
      tx_prots <- prots[prots$tx_id==tx,]
      
      
      if (nrow(tx_prots)==0){   # transcripts for which no protein sequence in ensDB, i.e. non-coding/untranslated transcripts?
        x <- list(ens_gene=ens_gene, tx_id = tx)
        lost_tx <- c(lost_tx, list(x))
        
      } else{
        for (j in 1:nrow(tx_prots)){  # loop over proteins (conservative design, in unexpected case multiple proteins would be listed for a certain transcript)
          prot_iso_id <- tx_prots$protein_id[j]
          cat(sprintf(" -- %s\n", prot_iso_id))
          prot_iso_seq <- tx_prots$protein_sequence[j]
          
          # browser()
          
          starts <- matched_pepts <- ends <- c()
          for (k in 1:length(pepts)){  # loop over peptides
            pept <- pepts[k]
            pos <- str_locate_all(prot_iso_seq, pattern = pept)[[1]]
            
            
            if ((nrow(pos) > 1) & (ens_gene != "ENSG00000198670")){   #apo(a) [ENSG00000198670] is only expected case
              warning(glue("\n{pept} matches multiple positions of the protein translation of transcript {tx}\n"))
              browser()
            }
            
            # some peptides will not be retrieved in some proteins,
            # eg in truncated proteins or for skipped exons
            # These are not included in the IRanges
            
            if (!all(is.na(pos))){ 
              starts <- c(starts, pos[,1])
              ends <- c(ends, pos[,2])
              matched_pepts <- c(matched_pepts, rep(pept, nrow(pos)))
            } else {
              # browser()
              x <- list(ens_gene=ens_gene, tx_id=tx, prot=prot_iso_id, pept=pept)
              lost_pepts <- c(lost_pepts, list(x))
            }
          }
          
          if (length(starts)>0){
            prot_r <- IRanges(
                start=starts, end=ends, 
                names=rep(prot_iso_id,length(starts))
            )
            
            prot_rd <- RangedData(
              prot_r,
              peptide = matched_pepts,
              space = "protein"
            )
            
            
            # browser()
            #genomic positions in GRCh37
            genome_r <- proteinToGenome(prot_r, ensdb)  #respects splicing
            names(genome_r) <- matched_pepts
            
            
            
            
            pept_prot_pos[[ens_gene]][[tx]][[prot_iso_id]] <- 
              list(protein_level= prot_rd,
                   genome_level= genome_r)
          }
          
          
          
          
        }
      }
    }
  }
  
  
  # browser()
  lost_pepts <- lost_pepts %>%
    bind_rows() %>%
    mutate(mapped_success = FALSE)
  
  # Overview of which peptides were sucessfully mapped to which protein isoforms
  mapping_success_overview <- tx_overview %>%
    left_join(dplyr::select(protein_peptides, ens_gene, Peptide_clean), 
              by=c("ens_gene"="ens_gene")) %>%
    left_join(dplyr::select(lost_pepts, ens_gene, tx_id=tx_id, Peptide_clean=pept, mapped_success),
              by = c("ens_gene", "tx_id", "Peptide_clean")) %>%
    mutate(mapped_success = is.na(mapped_success))
  
  
  return(list(
    mapping_success_overview=mapping_success_overview,
    pept_prot_pos=pept_prot_pos, 
    lost_pepts=lost_pepts,
    lost_tx=bind_rows(lost_tx)))
}


