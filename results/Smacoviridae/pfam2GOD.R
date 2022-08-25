#BiocManager::install("edgeR")

library("taxize")
library("usethis")
library("myTAI")
library("tidyverse")

#PosArgs <- "Geminiviridae_all_homologues_19062021_0435"
PosArgs <- as.character(commandArgs(trailingOnly = TRUE))
#PosArgs <- "Anelloviridae_all_homologues_19062021_0409"
GETHOMS.outdir <- PosArgs[1] 

files_in_GETHOMS <- list.files(path = GETHOMS.outdir)
taxids_in_GETHOMS <- unique(gsub("\\..*","",gsub("^.*taxid_","",files_in_GETHOMS[grep("taxid_",files_in_GETHOMS)]))) 

#GET INFO DEPARTING FROM TAXID
taxids_in_GETHOMS.df <- data.frame(taxid=taxids_in_GETHOMS)

#fill with id2names info
#taxize::use_entrez()

for(i in 1:length(taxids_in_GETHOMS.df$taxid)){
  start_time <- Sys.time()
  id2name_out <- id2name(as.numeric(taxids_in_GETHOMS.df$taxid[i]),db="ncbi")
  name <- id2name_out[[1]][2]
  rank <- id2name_out[[1]][3]
  status <- id2name_out[[1]][4]
  taxids_in_GETHOMS.df$name[i] <- name
  taxids_in_GETHOMS.df$rank[i] <- rank
  taxids_in_GETHOMS.df$status[i] <- status
  end_time <- Sys.time()
  time_needed <- end_time - start_time 
  taxids_in_GETHOMS.df$time[i]<-time_needed
  if(i>=8){
    last_8_steps <- sum(taxids_in_GETHOMS.df$time[(i-8):i])
    if(last_8_steps<1){
      Sys.sleep(1-last_8_steps)
    }
  }
}

taxids_in_GETHOMS.df$time <- NULL

# Up to this point I got species name info of contained files within a
# pangenomic input folder

# Ahora quiero hacer un listado de los clusters de homologos por output 
# pangenomico

# Cluster list of algorithm COG for taxon
cluster_list_algCOG.file <- files_in_GETHOMS[grep("algCOG",files_in_GETHOMS)][grep(".cluster_list",files_in_GETHOMS[grep("algCOG",files_in_GETHOMS)])]
cluster_list_algOMCL.file <- files_in_GETHOMS[grep("algOMCL",files_in_GETHOMS)][grep(".cluster_list",files_in_GETHOMS[grep("algOMCL",files_in_GETHOMS)])]
Command_tmp_tsv <- " | sed s/^.*taxid_/:/g | sed s/.gbk//g | sed s/cluster\\\ /cluster/g | sed s/=//g | sed s/file:\\\ /file/g | sed -e ':a' -e 'N' -e '$!ba' -e 's/\\n:/_/g' | sed -e ':a' -e 'N' -e '$!ba' -e 's/\\\ /\\t/g' | sed -e ':a' -e 'N' -e '$!ba' -e 's/fna_/fna\\t/g'"
Command_COG <- paste("cat ",GETHOMS.outdir,"/",cluster_list_algCOG.file," ",Command_tmp_tsv," > temporary_COG.tsv", sep="")
Command_OMCL <- paste("cat ",GETHOMS.outdir,"/",cluster_list_algOMCL.file," ",Command_tmp_tsv," > temporary_OMCL.tsv", sep="")

system(Command_COG)
system(Command_OMCL)

Correct_COG <- "sed -i '/Pfam/!s/file/\\tfile/1' temporary_COG.tsv "
Correct_OMCL <- "sed -i '/Pfam/!s/file/\\tfile/1' temporary_OMCL.tsv "

system(Correct_COG)
system(Correct_OMCL)

COG.tsv <- read.table(file = 'temporary_COG.tsv', sep = '\t', header = FALSE, fill = TRUE)
OMCL.tsv <- read.table(file = 'temporary_OMCL.tsv', sep = '\t', header = FALSE, fill = TRUE)

# Give proper format to data frames of clusters
COG.df <- data.frame(cluster_name = as.character(gsub("cluster","",COG.tsv$V1)),
                     size = as.numeric(gsub("size","",COG.tsv$V2)),
                     taxa = as.numeric(gsub("taxa","",COG.tsv$V3)),
                     Pfam_IDs = as.character(gsub(",","_",gsub("Pfam","",COG.tsv$V4))),
                     aa_file = as.character(gsub("file","",COG.tsv$V5)),
                     na_file = as.character(gsub("dnafile","",COG.tsv$V6)),
                     taxids = as.character(gsub("","",COG.tsv$V7)))

OMCL.df <- data.frame(cluster_name = as.character(gsub("cluster","",OMCL.tsv$V1)),
                     size = as.numeric(gsub("size","",OMCL.tsv$V2)),
                     taxa = as.numeric(gsub("taxa","",OMCL.tsv$V3)),
                     Pfam_IDs = as.character(gsub(",","_",gsub("Pfam","",OMCL.tsv$V4))),
                     aa_file = as.character(gsub("file","",OMCL.tsv$V5)),
                     na_file = as.character(gsub("dnafile","",OMCL.tsv$V6)),
                     taxids = as.character(gsub("","",OMCL.tsv$V7)))

# Add percentage of prevalence for each cluster
total_taxa <- length(taxids_in_GETHOMS)
COG.df$Percentage_prevalence <- (COG.df$taxa / total_taxa)*100
OMCL.df$Percentage_prevalence <- (OMCL.df$taxa / total_taxa)*100

# Build a presence matrix for taxids and cluster
COG_Taxids_clusters.mx <- matrix(nrow = length(taxids_in_GETHOMS), ncol = length(COG.df$cluster_name)) 
OMCL_Taxids_clusters.mx <- matrix(nrow = length(taxids_in_GETHOMS), ncol = length(OMCL.df$cluster_name)) 
rownames(COG_Taxids_clusters.mx) <- taxids_in_GETHOMS
rownames(OMCL_Taxids_clusters.mx) <- taxids_in_GETHOMS
colnames(COG_Taxids_clusters.mx) <- COG.df$cluster_name
colnames(OMCL_Taxids_clusters.mx) <- OMCL.df$cluster_name

# function to check if taxid has a gethoms cluster
IsClusterInTaxid_COG <- function(t,c){ #for COG clusters
  ID <- t
  CL <- c
  CL_IDs <- COG.df$taxids[which(COG.df$cluster_name == c)]
  if(grepl(ID,CL_IDs)){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

IsClusterInTaxid_OMCL <- function(t,c){ #for OMCL clusters
  ID <- t
  CL <- c
  CL_IDs <- OMCL.df$taxids[which(OMCL.df$cluster_name == c)]
  if(grepl(ID,CL_IDs)){
    return(TRUE)
  } else {
    return(FALSE)
  }
}


for(t in 1:nrow(COG_Taxids_clusters.mx)){
  for(c in 1:ncol(COG_Taxids_clusters.mx)){
    COG_Taxids_clusters.mx[t,c]<-IsClusterInTaxid_COG(rownames(COG_Taxids_clusters.mx)[t],colnames(COG_Taxids_clusters.mx)[c])
  }
}

for(t in 1:nrow(OMCL_Taxids_clusters.mx)){
  for(c in 1:ncol(OMCL_Taxids_clusters.mx)){
    OMCL_Taxids_clusters.mx[t,c]<-IsClusterInTaxid_OMCL(rownames(OMCL_Taxids_clusters.mx)[t],colnames(OMCL_Taxids_clusters.mx)[c])
  }
}

# Presence matrix done

# Ahora necesito buscar las funciones atribuidas a los Pfams

# Read GODomaninMiner
# molecular function
GOD_mf <- read.table(file = "../../../Pfam_associations/pfam_gomf_most_specific.txt", sep = ';', header = TRUE, fill = TRUE)
# biological process
GOD_bp <- read.table(file = "../../../Pfam_associations/pfam_gobp_most_specific.txt", sep = ';', header = TRUE, fill = TRUE)
# celular component
GOD_cc <- read.table(file = "../../../Pfam_associations/pfam_gocc_most_specific.txt", sep = ';', header = TRUE, fill = TRUE)

# Read pfam2go
raw_pfam2go <- readLines("../../../Pfam_associations/pfam2go")
skip_com_pfam2go <- raw_pfam2go[-grep("^!",raw_pfam2go)]
pfam2goheader <- c("Pfam_ID", "Name", "GO_function","GO_ID")
#skip_com_pfam2go[1]
pfam2go_pfams <- str_match(skip_com_pfam2go, "Pfam:(.*?)\\s")[,2]
pfam2go_names <- str_match(skip_com_pfam2go, "\\s(.*?)\\s")[,2]
pfam2go_gofunction <- str_match(skip_com_pfam2go, ">\\sGO:(.*?)\\s;")[,2]
pfam2go_goID <- gsub(" ; GO:","",str_match(skip_com_pfam2go, "\\s;\\sGO:(.*?).*")[,1])

pfam2go.df <- data.frame(Pfam_ID = pfam2go_pfams,
                         Pfam_name = pfam2go_names,
                         GO_function = pfam2go_gofunction,
                         GO_ID = pfam2go_goID)

# Now that I have many sources to associate with pfam its time to build an 
# inclusive database and plot results

# I have to create matrices of presence of GO_IDs similar to 
# OMCL_Taxids_clusters.mx

###############################################################################
#Function: cluster_hits_in_links()
#
#ARGUMENTS:
#
#   COG_or_OMCL: character string, either "COG" or "OMCL"
#
#   cluster: character string, the name of a gethoms cluster, (g.e. "12_MP")
#
#   link_db: character string, one of: 
#         "pfam2go" pfam2go maintained by Interpro team
#         "GOD_bp"  Gene Ontology Domain Miner biological process
#         "GOD_cc"  Gene Ontology Domain Miner celular component
#         "GOD_mf"  Gene Ontology Domain Miner molecular function
#
#   which_out: character string, one of:
#         "names"   function returns names of GO hits (only if link_db = "pfam2go")
#         "functions"   function returns list of functions (only if link_db = "pfam2go")
#         "GO_ids" function returns list of GO IDs 
#         "full_hits" function returns full rows of hits in selected link_db
#         "pfam" function returns only pfam IDs list
cluster_hits_in_links <- function(COG_or_OMCL,cluster,link_db,which_out){
  if(COG_or_OMCL == "COG"){
    indice_cluster <- grep(cluster,COG.df$cluster_name)
    pfams_in_clust <- COG.df$Pfam_IDs[indice_cluster]
  }
  if(COG_or_OMCL == "OMCL"){
    indice_cluster <- grep(cluster,OMCL.df$cluster_name)
    pfams_in_clust <- OMCL.df$Pfam_IDs[indice_cluster]
  }
  pfams_vector <- str_split(pfams_in_clust,"_")[[1]]
  pfams2check <- pfams_vector[pfams_vector != ""]
  if(length(pfams2check) < 1){
    warning("Cluster could not be annotated with Pfam")
    return(NA)
  }
  if(which_out == "pfam"){
    return(pfams2check)
  }
  #pfam2go
  if(link_db == "pfam2go"){
    pfam2go_indexohits <- which(as.character(pfam2go.df$Pfam_ID)%in%pfams2check)
    pfam_names <- unique(pfam2go.df$Pfam_name[pfam2go_indexohits])
    pfam_functions <- unique(pfam2go.df$GO_function[pfam2go_indexohits])
    pfam_GOids <- unique(pfam2go.df$GO_ID[pfam2go_indexohits])
    pfam_fullhits <- pfam2go.df[pfam2go_indexohits,]
    if(which_out == "names"){
      return(pfam_names)
    }
    if(which_out == "functions"){
      return(pfam_functions)
    }
    if(which_out == "GO_ids"){
      return(pfam_GOids)
    }
    if(which_out == "full_hits"){
      return(pfam_fullhits)
    }
  }
  #GOD_hits
  if(link_db == "GOD_bp"){
    GOD_bp_indexohits <- which(as.character(GOD_bp$PFAM)%in%pfams2check)
    GOD_bp_GOids <- gsub("GO:","",unique(GOD_bp$GO[GOD_bp_indexohits]))
    GOD_bp_fullhits <- GOD_bp[GOD_bp_indexohits,]
    if(which_out == "GO_ids"){
      return(GOD_bp_GOids)
    }else{
      return(GOD_bp_fullhits)
    }
  }
  if(link_db == "GOD_cc"){
    GOD_cc_indexohits <- which(as.character(GOD_cc$PFAM)%in%pfams2check)
    GOD_cc_GOids <- gsub("GO:","",unique(GOD_cc$GO[GOD_cc_indexohits]))
    GOD_cc_fullhits <- GOD_cc[GOD_cc_indexohits,]
    if(which_out == "GO_ids"){
      return(GOD_cc_GOids)
    }else{
      return(GOD_cc_fullhits)
    }
  }
  if(link_db == "GOD_mf"){
    GOD_mf_indexohits <- which(as.character(GOD_mf$PFAM)%in%pfams2check)
    GOD_mf_GOids <- gsub("GO:","",unique(GOD_mf$GO[GOD_mf_indexohits]))
    GOD_mf_fullhits <- GOD_mf[GOD_mf_indexohits,]
    if(which_out == "GO_ids"){
      return(GOD_mf_GOids)
    }else{
      return(GOD_mf_fullhits)
    }
  }
}
###########################################################################
# Prueba de la funcion
#cluster_hits_in_links(COG_or_OMCL = "OMCL",cluster = OMCL_clusts[1],link_db = "pfam2go",which_out = "pfam")

#List of cluster names of each algorithm


OMCL_clusts <- OMCL.df$cluster_name[grep("PF",OMCL.df$Pfam_IDs)]
COG_clusts <- COG.df$cluster_name[grep("PF",COG.df$Pfam_IDs)]

#FILL OMCL MATRICES
#OMCL_names_pfam2go.mx
OMCL_names_pfam2go.mx <- as.data.frame(matrix(nrow = length(OMCL_clusts),ncol = length(unique(pfam2go.df$Pfam_name))))
colnames(OMCL_names_pfam2go.mx) <- unique(pfam2go.df$Pfam_name)
rownames(OMCL_names_pfam2go.mx) <- OMCL_clusts
for(f in 1:length(OMCL_clusts)){
  for(c in 1:length(unique(pfam2go.df$Pfam_name))){
    res_myfun <- cluster_hits_in_links(COG_or_OMCL = "OMCL",cluster = OMCL_clusts[f],link_db = "pfam2go",which_out = "names")
    if(unique(pfam2go.df$Pfam_name)[c]%in%res_myfun){
      OMCL_names_pfam2go.mx[f,c]<-TRUE
    }
  }
}

tochos <- as.numeric(nrow(OMCL_names_pfam2go.mx)*ncol(OMCL_names_pfam2go.mx))
nases <- as.numeric(table(is.na(OMCL_names_pfam2go.mx))["TRUE"] %>% unname())

if(tochos != nases){
  OMCL_names_pfam2go.mx.sub <- tryCatch({OMCL_names_pfam2go.mx[,unname(which(apply(OMCL_names_pfam2go.mx, 2, function(x) any(grepl(TRUE, x)))))]})
  OMCL_names_pfam2go.mx.sub <- as.data.frame(OMCL_names_pfam2go.mx.sub)
  OMCL_names_pfam2go.mx.sub.2 <- tryCatch({OMCL_names_pfam2go.mx.sub[unname(which(apply(OMCL_names_pfam2go.mx.sub, 1, function(x) any(grepl(TRUE, x))))),]})
  OMCL_names_pfam2go.mx.sub.2 <- as.data.frame(OMCL_names_pfam2go.mx.sub.2)
}

#OMCL_GO_function_pfam2go.mx
OMCL_GO_function_pfam2go.mx <- as.data.frame(matrix(nrow = length(OMCL_clusts),ncol = length(unique(pfam2go.df$GO_function))))
colnames(OMCL_GO_function_pfam2go.mx) <- unique(pfam2go.df$GO_function)
rownames(OMCL_GO_function_pfam2go.mx) <- OMCL_clusts
for(f in 1:length(OMCL_clusts)){
  for(c in 1:length(unique(pfam2go.df$GO_function))){
    res_myfun <- cluster_hits_in_links(COG_or_OMCL = "OMCL",cluster = OMCL_clusts[f],link_db = "pfam2go",which_out = "functions")
    if(unique(pfam2go.df$GO_function)[c]%in%res_myfun){
      OMCL_GO_function_pfam2go.mx[f,c]<-TRUE
    }
  }
}

tochos <- as.numeric(nrow(OMCL_GO_function_pfam2go.mx)*ncol(OMCL_GO_function_pfam2go.mx))
nases <- as.numeric(table(is.na(OMCL_GO_function_pfam2go.mx))["TRUE"] %>% unname())

if(tochos != nases){
  OMCL_GO_function_pfam2go.mx.sub <- tryCatch({OMCL_GO_function_pfam2go.mx[,unname(which(apply(OMCL_GO_function_pfam2go.mx, 2, function(x) any(grepl(TRUE, x)))))]})
  OMCL_GO_function_pfam2go.mx.sub <- as.data.frame(OMCL_GO_function_pfam2go.mx.sub)
  OMCL_GO_function_pfam2go.mx.sub.2 <- tryCatch({OMCL_GO_function_pfam2go.mx.sub[unname(which(apply(OMCL_GO_function_pfam2go.mx.sub, 1, function(x) any(grepl(TRUE, x))))),]})
  OMCL_GO_function_pfam2go.mx.sub.2 <- as.data.frame(OMCL_GO_function_pfam2go.mx.sub.2)
}


#OMCL_GO_ID_pfam2go.mx
OMCL_GO_ID_pfam2go.mx <- as.data.frame(matrix(nrow = length(OMCL_clusts),ncol = length(unique(pfam2go.df$GO_ID))))
colnames(OMCL_GO_ID_pfam2go.mx) <- unique(pfam2go.df$GO_ID)
rownames(OMCL_GO_ID_pfam2go.mx) <- OMCL_clusts
for(f in 1:length(OMCL_clusts)){
  for(c in 1:length(unique(pfam2go.df$GO_ID))){
    res_myfun <- cluster_hits_in_links(COG_or_OMCL = "OMCL",cluster = OMCL_clusts[f],link_db = "pfam2go",which_out = "GO_ids")
    if(unique(pfam2go.df$GO_ID)[c]%in%res_myfun){
      OMCL_GO_ID_pfam2go.mx[f,c]<-TRUE
    }
  }
}

tochos <- as.numeric(nrow(OMCL_GO_ID_pfam2go.mx)*ncol(OMCL_GO_ID_pfam2go.mx))
nases <- as.numeric(table(is.na(OMCL_GO_ID_pfam2go.mx))["TRUE"] %>% unname())

if(tochos != nases){
  OMCL_GO_ID_pfam2go.mx.sub <- tryCatch({OMCL_GO_ID_pfam2go.mx[,unname(which(apply(OMCL_GO_ID_pfam2go.mx, 2, function(x) any(grepl(TRUE, x)))))]})
  OMCL_GO_ID_pfam2go.mx.sub <- as.data.frame(OMCL_GO_ID_pfam2go.mx.sub)
  OMCL_GO_ID_pfam2go.mx.sub.2 <- tryCatch({OMCL_GO_ID_pfam2go.mx.sub[unname(which(apply(OMCL_GO_ID_pfam2go.mx.sub, 1, function(x) any(grepl(TRUE, x))))),]})
  OMCL_GO_ID_pfam2go.mx.sub.2 <- as.data.frame(OMCL_GO_ID_pfam2go.mx.sub.2)
}

#OMCL_GO_ID_GOD_bp.mx
OMCL_GO_ID_GOD_bp.mx <- as.data.frame(matrix(nrow = length(OMCL_clusts),ncol = length(unique(GOD_bp$GO))))
colnames(OMCL_GO_ID_GOD_bp.mx) <- unique(unique(GOD_bp$GO))
rownames(OMCL_GO_ID_GOD_bp.mx) <- OMCL_clusts
for(f in 1:length(OMCL_clusts)){
  res_myfun <- cluster_hits_in_links(COG_or_OMCL = "OMCL",cluster = OMCL_clusts[f],link_db = "GOD_bp",which_out = "GO_ids")
  for(c in 1:length(unique(GOD_bp$GO))){
    if(gsub("GO:","",unique(GOD_bp$GO))[c]%in%res_myfun){
      OMCL_GO_ID_GOD_bp.mx[f,c]<-TRUE
    }
  }
}

tochos <- as.numeric(nrow(OMCL_GO_ID_GOD_bp.mx)*ncol(OMCL_GO_ID_GOD_bp.mx))
nases <- as.numeric(table(is.na(OMCL_GO_ID_GOD_bp.mx))["TRUE"] %>% unname())

if(tochos != nases){
  OMCL_GO_ID_GOD_bp.mx.sub <- tryCatch({OMCL_GO_ID_GOD_bp.mx[,unname(which(apply(OMCL_GO_ID_GOD_bp.mx, 2, function(x) any(grepl(TRUE, x)))))]})
  OMCL_GO_ID_GOD_bp.mx.sub <- as.data.frame(OMCL_GO_ID_GOD_bp.mx.sub)
  OMCL_GO_ID_GOD_bp.mx.sub.2 <- tryCatch({OMCL_GO_ID_GOD_bp.mx.sub[unname(which(apply(OMCL_GO_ID_GOD_bp.mx.sub, 1, function(x) any(grepl(TRUE, x))))),]})
  OMCL_GO_ID_GOD_bp.mx.sub.2 <- as.data.frame(OMCL_GO_ID_GOD_bp.mx.sub.2)
}

#OMCL_GO_ID_GOD_cc.mx
OMCL_GO_ID_GOD_cc.mx <- as.data.frame(matrix(nrow = length(OMCL_clusts),ncol = length(unique(GOD_cc$GO))))
colnames(OMCL_GO_ID_GOD_cc.mx) <- unique(unique(GOD_cc$GO))
rownames(OMCL_GO_ID_GOD_cc.mx) <- OMCL_clusts
for(f in 1:length(OMCL_clusts)){
  res_myfun <- cluster_hits_in_links(COG_or_OMCL = "OMCL",cluster = OMCL_clusts[f],link_db = "GOD_cc",which_out = "GO_ids")
  for(c in 1:length(unique(GOD_cc$GO))){
    if(gsub("GO:","",unique(GOD_cc$GO))[c]%in%res_myfun){
      OMCL_GO_ID_GOD_cc.mx[f,c]<-TRUE
    }
  }
}

tochos <- as.numeric(nrow(OMCL_GO_ID_GOD_cc.mx)*ncol(OMCL_GO_ID_GOD_cc.mx))
nases <- as.numeric(table(is.na(OMCL_GO_ID_GOD_cc.mx))["TRUE"] %>% unname())

if(tochos != nases){
  OMCL_GO_ID_GOD_cc.mx.sub <- tryCatch({OMCL_GO_ID_GOD_cc.mx[,unname(which(apply(OMCL_GO_ID_GOD_cc.mx, 2, function(x) any(grepl(TRUE, x)))))]})
  OMCL_GO_ID_GOD_cc.mx.sub <- as.data.frame(OMCL_GO_ID_GOD_cc.mx.sub)
  OMCL_GO_ID_GOD_cc.mx.sub.2 <- tryCatch({OMCL_GO_ID_GOD_cc.mx.sub[unname(which(apply(OMCL_GO_ID_GOD_cc.mx.sub, 1, function(x) any(grepl(TRUE, x))))),]})
  OMCL_GO_ID_GOD_cc.mx.sub.2 <- as.data.frame(OMCL_GO_ID_GOD_cc.mx.sub.2)
}

#OMCL_GO_ID_GOD_mf.mx
OMCL_GO_ID_GOD_mf.mx <- as.data.frame(matrix(nrow = length(OMCL_clusts),ncol = length(unique(GOD_mf$GO))))
colnames(OMCL_GO_ID_GOD_mf.mx) <- unique(unique(GOD_mf$GO))
rownames(OMCL_GO_ID_GOD_mf.mx) <- OMCL_clusts
for(f in 1:length(OMCL_clusts)){
  res_myfun <- cluster_hits_in_links(COG_or_OMCL = "OMCL",cluster = OMCL_clusts[f],link_db = "GOD_mf",which_out = "GO_ids")
  for(c in 1:length(unique(GOD_mf$GO))){
    if(gsub("GO:","",unique(GOD_mf$GO))[c]%in%res_myfun){
      OMCL_GO_ID_GOD_mf.mx[f,c]<-TRUE
    }
  }
}



tochos <- as.numeric(nrow(OMCL_GO_ID_GOD_mf.mx)*ncol(OMCL_GO_ID_GOD_mf.mx))
nases <- as.numeric(table(is.na(OMCL_GO_ID_GOD_mf.mx))["TRUE"] %>% unname())


if(tochos != nases){
  OMCL_GO_ID_GOD_mf.mx.sub <- tryCatch({OMCL_GO_ID_GOD_mf.mx[,unname(which(apply(OMCL_GO_ID_GOD_mf.mx, 2, function(x) any(grepl(TRUE, x)))))]})
  OMCL_GO_ID_GOD_mf.mx.sub <- as.data.frame(OMCL_GO_ID_GOD_mf.mx.sub)
  OMCL_GO_ID_GOD_mf.mx.sub.2 <- tryCatch({OMCL_GO_ID_GOD_mf.mx.sub[unname(which(apply(OMCL_GO_ID_GOD_mf.mx.sub, 1, function(x) any(grepl(TRUE, x))))),]})
  OMCL_GO_ID_GOD_mf.mx.sub.2 <- as.data.frame(OMCL_GO_ID_GOD_mf.mx.sub.2)
}




#FILL COG MATRICES
#COG_names_pfam2go.mx
COG_names_pfam2go.mx <- as.data.frame(matrix(nrow = length(COG_clusts),ncol = length(unique(pfam2go.df$Pfam_name))))
colnames(COG_names_pfam2go.mx) <- unique(pfam2go.df$Pfam_name)
rownames(COG_names_pfam2go.mx) <- COG_clusts
for(f in 1:length(COG_clusts)){
  for(c in 1:length(unique(pfam2go.df$Pfam_name))){
    res_myfun <- cluster_hits_in_links(COG_or_OMCL = "COG",cluster = COG_clusts[f],link_db = "pfam2go",which_out = "names")
    if(unique(pfam2go.df$Pfam_name)[c]%in%res_myfun){
      COG_names_pfam2go.mx[f,c]<-TRUE
    }
  }
}

tochos <- as.numeric(nrow(COG_names_pfam2go.mx)*ncol(COG_names_pfam2go.mx))
nases <- as.numeric(table(is.na(COG_names_pfam2go.mx))["TRUE"] %>% unname())


if(tochos != nases){
  COG_names_pfam2go.mx.sub <- tryCatch({COG_names_pfam2go.mx[,unname(which(apply(COG_names_pfam2go.mx, 2, function(x) any(grepl(TRUE, x)))))]})
  COG_names_pfam2go.mx.sub <- as.data.frame(COG_names_pfam2go.mx.sub)
  COG_names_pfam2go.mx.sub.2 <- tryCatch({COG_names_pfam2go.mx.sub[unname(which(apply(COG_names_pfam2go.mx.sub, 1, function(x) any(grepl(TRUE, x))))),]})
  COG_names_pfam2go.mx.sub.2 <- as.data.frame(COG_names_pfam2go.mx.sub.2)
}


#COG_GO_function_pfam2go.mx
COG_GO_function_pfam2go.mx <- as.data.frame(matrix(nrow = length(COG_clusts),ncol = length(unique(pfam2go.df$GO_function))))
colnames(COG_GO_function_pfam2go.mx) <- unique(pfam2go.df$GO_function)
rownames(COG_GO_function_pfam2go.mx) <- COG_clusts
for(f in 1:length(COG_clusts)){
  for(c in 1:length(unique(pfam2go.df$GO_function))){
    res_myfun <- cluster_hits_in_links(COG_or_OMCL = "COG",cluster = COG_clusts[f],link_db = "pfam2go",which_out = "functions")
    if(unique(pfam2go.df$GO_function)[c]%in%res_myfun){
      COG_GO_function_pfam2go.mx[f,c]<-TRUE
    }
  }
}

tochos <- as.numeric(nrow(COG_GO_function_pfam2go.mx)*ncol(COG_GO_function_pfam2go.mx))
nases <- as.numeric(table(is.na(COG_GO_function_pfam2go.mx))["TRUE"] %>% unname())


if(tochos != nases){
  COG_GO_function_pfam2go.mx.sub <- tryCatch({COG_GO_function_pfam2go.mx[,unname(which(apply(COG_GO_function_pfam2go.mx, 2, function(x) any(grepl(TRUE, x)))))]})
  COG_GO_function_pfam2go.mx.sub <- as.data.frame(COG_GO_function_pfam2go.mx.sub)
  COG_GO_function_pfam2go.mx.sub.2 <- tryCatch({COG_GO_function_pfam2go.mx.sub[unname(which(apply(COG_GO_function_pfam2go.mx.sub, 1, function(x) any(grepl(TRUE, x))))),]})
  COG_GO_function_pfam2go.mx.sub.2 <- as.data.frame(COG_GO_function_pfam2go.mx.sub.2)
}





#COG_GO_ID_pfam2go.mx
COG_GO_ID_pfam2go.mx <- as.data.frame(matrix(nrow = length(COG_clusts),ncol = length(unique(pfam2go.df$GO_ID))))
colnames(COG_GO_ID_pfam2go.mx) <- unique(pfam2go.df$GO_ID)
rownames(COG_GO_ID_pfam2go.mx) <- COG_clusts
for(f in 1:length(COG_clusts)){
  for(c in 1:length(unique(pfam2go.df$GO_ID))){
    res_myfun <- cluster_hits_in_links(COG_or_OMCL = "COG",cluster = COG_clusts[f],link_db = "pfam2go",which_out = "GO_ids")
    if(unique(pfam2go.df$GO_ID)[c]%in%res_myfun){
      COG_GO_ID_pfam2go.mx[f,c]<-TRUE
    }
  }
}


tochos <- as.numeric(nrow(COG_GO_ID_pfam2go.mx)*ncol(COG_GO_ID_pfam2go.mx))
nases <- as.numeric(table(is.na(COG_GO_ID_pfam2go.mx))["TRUE"] %>% unname())


if(tochos != nases){
  COG_GO_ID_pfam2go.mx.sub <- tryCatch({COG_GO_ID_pfam2go.mx[,unname(which(apply(COG_GO_ID_pfam2go.mx, 2, function(x) any(grepl(TRUE, x)))))]})
  COG_GO_ID_pfam2go.mx.sub <- as.data.frame(COG_GO_ID_pfam2go.mx.sub)
  COG_GO_ID_pfam2go.mx.sub.2 <- tryCatch({COG_GO_ID_pfam2go.mx.sub[unname(which(apply(COG_GO_ID_pfam2go.mx.sub, 1, function(x) any(grepl(TRUE, x))))),]})
  COG_GO_ID_pfam2go.mx.sub.2 <- as.data.frame(COG_GO_ID_pfam2go.mx.sub.2)
}


#COG_GO_ID_GOD_bp.mx
COG_GO_ID_GOD_bp.mx <- as.data.frame(matrix(nrow = length(COG_clusts),ncol = length(unique(GOD_bp$GO))))
colnames(COG_GO_ID_GOD_bp.mx) <- unique(unique(GOD_bp$GO))
rownames(COG_GO_ID_GOD_bp.mx) <- COG_clusts
for(f in 1:length(COG_clusts)){
  res_myfun <- cluster_hits_in_links(COG_or_OMCL = "COG",cluster = COG_clusts[f],link_db = "GOD_bp",which_out = "GO_ids")
  for(c in 1:length(unique(GOD_bp$GO))){
    if(gsub("GO:","",unique(GOD_bp$GO))[c]%in%res_myfun){
      COG_GO_ID_GOD_bp.mx[f,c]<-TRUE
    }
  }
}

tochos <- as.numeric(nrow(COG_GO_ID_GOD_bp.mx)*ncol(COG_GO_ID_GOD_bp.mx))
nases <- as.numeric(table(is.na(COG_GO_ID_GOD_bp.mx))["TRUE"] %>% unname())


if(tochos != nases){
  COG_GO_ID_GOD_bp.mx.sub <- tryCatch({COG_GO_ID_GOD_bp.mx[,unname(which(apply(COG_GO_ID_GOD_bp.mx, 2, function(x) any(grepl(TRUE, x)))))]})
  COG_GO_ID_GOD_bp.mx.sub <- as.data.frame(COG_GO_ID_GOD_bp.mx.sub)
  COG_GO_ID_GOD_bp.mx.sub.2 <-tryCatch({COG_GO_ID_GOD_bp.mx.sub[unname(which(apply(COG_GO_ID_GOD_bp.mx.sub, 1, function(x) any(grepl(TRUE, x))))),]})
  COG_GO_ID_GOD_bp.mx.sub.2 <- as.data.frame(COG_GO_ID_GOD_bp.mx.sub.2)
}




#COG_GO_ID_GOD_cc.mx
COG_GO_ID_GOD_cc.mx <- as.data.frame(matrix(nrow = length(COG_clusts),ncol = length(unique(GOD_cc$GO))))
colnames(COG_GO_ID_GOD_cc.mx) <- unique(unique(GOD_cc$GO))
rownames(COG_GO_ID_GOD_cc.mx) <- COG_clusts
for(f in 1:length(COG_clusts)){
  res_myfun <- cluster_hits_in_links(COG_or_OMCL = "COG",cluster = COG_clusts[f],link_db = "GOD_cc",which_out = "GO_ids")
  for(c in 1:length(unique(GOD_cc$GO))){
    if(gsub("GO:","",unique(GOD_cc$GO))[c]%in%res_myfun){
      COG_GO_ID_GOD_cc.mx[f,c]<-TRUE
    }
  }
}

tochos <- as.numeric(nrow(COG_GO_ID_GOD_cc.mx)*ncol(COG_GO_ID_GOD_cc.mx))
nases <- as.numeric(table(is.na(COG_GO_ID_GOD_cc.mx))["TRUE"] %>% unname())


if(tochos != nases){
  COG_GO_ID_GOD_cc.mx.sub <- tryCatch({COG_GO_ID_GOD_cc.mx[,unname(which(apply(COG_GO_ID_GOD_cc.mx, 2, function(x) any(grepl(TRUE, x)))))]})
  COG_GO_ID_GOD_cc.mx.sub <- as.data.frame(COG_GO_ID_GOD_cc.mx.sub)
  COG_GO_ID_GOD_cc.mx.sub.2 <-tryCatch({COG_GO_ID_GOD_cc.mx.sub[unname(which(apply(COG_GO_ID_GOD_cc.mx.sub, 1, function(x) any(grepl(TRUE, x))))),]})
  COG_GO_ID_GOD_cc.mx.sub.2 <- as.data.frame(COG_GO_ID_GOD_cc.mx.sub.2)
}

#COG_GO_ID_GOD_mf.mx
COG_GO_ID_GOD_mf.mx <- as.data.frame(matrix(nrow = length(COG_clusts),ncol = length(unique(GOD_mf$GO))))
colnames(COG_GO_ID_GOD_mf.mx) <- unique(unique(GOD_mf$GO))
rownames(COG_GO_ID_GOD_mf.mx) <- COG_clusts
for(f in 1:length(COG_clusts)){
  res_myfun <- cluster_hits_in_links(COG_or_OMCL = "COG",cluster = COG_clusts[f],link_db = "GOD_mf",which_out = "GO_ids")
  for(c in 1:length(unique(GOD_mf$GO))){
    if(gsub("GO:","",unique(GOD_mf$GO))[c]%in%res_myfun){
      COG_GO_ID_GOD_mf.mx[f,c]<-TRUE
    }
  }
}

tochos <- as.numeric(nrow(COG_GO_ID_GOD_mf.mx)*ncol(COG_GO_ID_GOD_mf.mx))
nases <- as.numeric(table(is.na(COG_GO_ID_GOD_mf.mx))["TRUE"] %>% unname())


if(tochos != nases){
  COG_GO_ID_GOD_mf.mx.sub <- tryCatch({COG_GO_ID_GOD_mf.mx[,unname(which(apply(COG_GO_ID_GOD_mf.mx, 2, function(x) any(grepl(TRUE, x)))))]})
  COG_GO_ID_GOD_mf.mx.sub <- as.data.frame(COG_GO_ID_GOD_mf.mx.sub)
  COG_GO_ID_GOD_mf.mx.sub.2 <- tryCatch({COG_GO_ID_GOD_mf.mx.sub[unname(which(apply(COG_GO_ID_GOD_mf.mx.sub, 1, function(x) any(grepl(TRUE, x))))),]})
  COG_GO_ID_GOD_mf.mx.sub.2 <- as.data.frame(COG_GO_ID_GOD_mf.mx.sub.2)
}


# Create folder to generate output
OUTPOST.dir <- paste(GETHOMS.outdir,"_postgethoms",sep="")
system(paste("mkdir -p",OUTPOST.dir))

# Write 16 objects
sink(paste(OUTPOST.dir,"/pfam2GOD.log",sep=""))
cat("Error log (unfinished files):\n\n")

tryCatch(
  {write.csv(COG.df,paste(OUTPOST.dir,"/COG.df.csv",sep=""))}
  ,error = cat("COG.df.csv"))

tryCatch(
  {write.csv(as.data.frame(COG_Taxids_clusters.mx),paste(OUTPOST.dir,"/COG_Taxids_clust_presence.mx.csv",sep=""))}
  ,error = cat("COG_Taxids_clust_presence.mx.csv"))

tryCatch(
  {write.csv(COG_GO_function_pfam2go.mx.sub.2,paste(OUTPOST.dir,"/COG_GO_function_pfam2go.csv",sep=""))}
  ,error = cat("COG_GO_function_pfam2go.csv"))

tryCatch(
  {write.csv(COG_GO_ID_GOD_bp.mx.sub.2,paste(OUTPOST.dir,"/COG_GO_ID_GOD_bp.csv",sep=""))},
  error = cat("COG_GO_ID_GOD_bp.csv"))

tryCatch(
  {write.csv(COG_GO_ID_GOD_cc.mx.sub.2,paste(OUTPOST.dir,"/COG_GO_ID_GOD_cc.csv",sep=""))},
  error = cat("COG_GO_ID_GOD_cc.csv"))

tryCatch(
  {write.csv(COG_GO_ID_GOD_mf.mx.sub.2,paste(OUTPOST.dir,"/COG_GO_ID_GOD_mf.csv",sep=""))}
  ,error = cat("COG_GO_ID_GOD_mf.csv"))

tryCatch(
  {write.csv(COG_GO_ID_pfam2go.mx.sub.2,paste(OUTPOST.dir,"/COG_GO_ID_pfam2go.csv",sep=""))}
  ,error = cat("COG_GO_ID_pfam2go.csv"))

tryCatch(
  {write.csv(COG_names_pfam2go.mx.sub.2,paste(OUTPOST.dir,"/COG_names_pfam2go.csv",sep=""))}
  ,error = cat("COG_names_pfam2go.csv"))

tryCatch(
  {write.csv(OMCL.df,paste(OUTPOST.dir,"/OMCL.df.csv",sep=""))}
  ,error = cat("OMCL.df.csv"))

tryCatch(
  {write.csv(as.data.frame(OMCL_Taxids_clusters.mx),paste(OUTPOST.dir,"/OMCL_Taxids_clust_presence.mx.csv",sep=""))}
  ,error = cat("OMCL_Taxids_clust_presence.mx.csv"))

tryCatch(
  {write.csv(OMCL_GO_function_pfam2go.mx.sub.2,paste(OUTPOST.dir,"/OMCL_GO_function_pfam2go.csv",sep=""))}
  ,error = cat("OMCL_GO_function_pfam2go.csv"))

tryCatch(
  {write.csv(OMCL_GO_ID_GOD_bp.mx.sub.2,paste(OUTPOST.dir,"/OMCL_GO_ID_GOD_bp.csv",sep=""))}
  ,error = cat("OMCL_GO_ID_GOD_bp.csv"))

tryCatch(
  {write.csv(OMCL_GO_ID_GOD_cc.mx.sub.2,paste(OUTPOST.dir,"/OMCL_GO_ID_GOD_cc.csv",sep=""))}
  ,error = cat("OMCL_GO_ID_GOD_cc.csv"))

tryCatch(
  {write.csv(OMCL_GO_ID_GOD_mf.mx.sub.2,paste(OUTPOST.dir,"/OMCL_GO_ID_GOD_mf.csv",sep=""))}
  ,error = cat("OMCL_GO_ID_GOD_mf.csv"))

tryCatch(
  {write.csv(OMCL_GO_ID_pfam2go.mx.sub.2,paste(OUTPOST.dir,"/OMCL_GO_ID_pfam2go.csv",sep=""))}
  ,error = cat("OMCL_GO_ID_pfam2go.csv"))

tryCatch(
  {write.csv(OMCL_names_pfam2go.mx.sub.2,paste(OUTPOST.dir,"/OMCL_names_pfam2go.csv",sep=""))}
  ,error = cat("OMCL_names_pfam2go.csv"))

sink()


