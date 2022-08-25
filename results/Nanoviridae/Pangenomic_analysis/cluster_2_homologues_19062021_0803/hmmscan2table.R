library("taxize")
library("usethis")
library("myTAI")
library("tidyverse")

#PosArgs <- "Geminiviridae_all_homologues_19062021_0435"
#PosArgs <- "Nanoviridae_all_homologues_19062021_0803"
PosArgs <- as.character(commandArgs(trailingOnly = TRUE))
GETHOMS.outdir <- PosArgs[1] 

gbk.fasta.pfam <- system(paste("cat ",GETHOMS.outdir,"/*.gbk.fasta.pfam",sep = ""), intern = TRUE)

total_hits = length(which(grepl("^>> PF",gbk.fasta.pfam)))

Domain_annotation <- as.data.frame(matrix(nrow = total_hits,ncol = 17))
colnames(Domain_annotation) <- c("Taxid","Num_q_global","Num_q_local","aa_length","Pfam_id","Pfam_name","HMMER_score","biasq","c_value","i_value","hmmfrom","hmmto","alifrom","alito","envfrom","envto","acc")

locount = 0
hitcount = 0
for (l in gbk.fasta.pfam){
  if (l == gbk.fasta.pfam[2]){
    locount = 0
  }
  locount = locount + 1
  if (any(grepl("# raw hmmscan command:",l))){
    t = str_match(l, "taxid_(.*?).gbk")[2] # if line contains "# raw hmmscan command:", update taxid
    qcount = 0
  }
  if (any(grepl("^Query:",l))){
    qcount = qcount + 1
    num_query_total = as.numeric(strsplit(gsub("[^[:alnum:] ]", "", l), " +")[[1]][2])
    aa_length = as.numeric(gsub("L","",strsplit(gsub("[^[:alnum:] ]", "", l), " +")[[1]][3]))
  }
  if (any(grepl("^>> PF",l))){
    PFid = strsplit(l, " +")[[1]][2]
    PF = strsplit(l, " +")[[1]][2]
    endline = length(strsplit(l, " +")[[1]])
    PFindex = grep(PF,strsplit(l, " +")[[1]])
    PFname = paste(strsplit(l, " +")[[1]][(PFindex + 1):endline], collapse = ' ')
    PFname = gsub(" ","_",PFname)
    #print(paste(t,num_query_total,qcount,aa_length,PFid,PFname))
    parametros <- strsplit(gsub("\\]","",gsub("\\[","",gbk.fasta.pfam[locount+4]))," +")
    parametros2 <- parametros[[1]][-which(gsub("\\.","",parametros[[1]])=="")]
    scoreq = parametros2[3]
    biasq = parametros2[4]
    c_value = parametros2[5]
    i_value = parametros2[6]
    hmmfrom = parametros2[7]
    hmmto = parametros2[8]
    alifrom = parametros2[9]
    alito = parametros2[10]
    envfrom = parametros2[11]
    envto = parametros2[12]
    acc = parametros2[13]
    hitcount = hitcount + 1
    print(paste(t,num_query_total,qcount,aa_length,PFid,PFname,scoreq,biasq,c_value,i_value,hmmfrom,hmmto,alifrom,alito,envfrom,envto,acc))
    Domain_annotation$Taxid[hitcount] <- t
    Domain_annotation$Num_q_global[hitcount] <- num_query_total
    Domain_annotation$Num_q_local[hitcount] <- qcount
    Domain_annotation$aa_length[hitcount] <- aa_length
    Domain_annotation$Pfam_id[hitcount] <- PFid
    Domain_annotation$Pfam_name[hitcount] <- PFname
    Domain_annotation$HMMER_score[hitcount] <- scoreq
    Domain_annotation$biasq[hitcount] <- biasq
    Domain_annotation$c_value[hitcount] <- c_value
    Domain_annotation$i_value[hitcount] <- i_value
    Domain_annotation$hmmfrom[hitcount] <- hmmfrom
    Domain_annotation$hmmto[hitcount] <- hmmto
    Domain_annotation$alifrom[hitcount] <- alifrom
    Domain_annotation$alito[hitcount] <- alito
    Domain_annotation$envfrom[hitcount] <- envfrom
    Domain_annotation$envto[hitcount] <- envto
    Domain_annotation$acc[hitcount] <- acc
  }
}

# LEER TODOS LOS .gbk.faa

gbk.faa <- system(paste("cat ",GETHOMS.outdir,"/*.gbk.faa",sep = ""), intern = TRUE)
pc_count <- length(grep("^>",gbk.faa))
oddsv <- seq(1, pc_count*2, 2)
evensv <- seq(2, pc_count*2, 2)
COG_GETHOMS.outdir <- system(paste("ls",GETHOMS.outdir, "| grep \"algCOG\" | grep -v \"cluster_list\""), intern = TRUE)
OMCL_GETHOMS.outdir <- system(paste("ls",GETHOMS.outdir, "| grep \"algOMCL\" | grep -v \"cluster_list\""), intern = TRUE)

for (i in 1:length(Domain_annotation$Num_q_local)){
  Domain_annotation$t_faa_index[i] <- oddsv[Domain_annotation$Num_q_local[i]]
  seqindex <- evensv[Domain_annotation$Num_q_local[i]]
  t = Domain_annotation$Taxid[i]
  t_headers <- grep(t,gbk.faa)
  t_sequences <- t_headers + 1
  t_faa <- gbk.faa[sort(c(t_headers,t_sequences))]
  Domain_annotation$Query_header[i] <- t_faa[Domain_annotation$t_faa_index[i]]
  Domain_annotation$Query_ID[i] <- str_match(Domain_annotation$Query_header[i], ">ID:(.*?)\\s")[2]
  Domain_annotation$Other_ID[i] <- gsub("\\^","",str_match(Domain_annotation$Query_header[i], "\\^(.*?)\\^")[1])
  Domain_annotation$Prot_name[i] <- str_match(Domain_annotation$Query_header[i], "\\.gbk\\|(.*?)\\|")[2]
  Domain_annotation$Full_sequence[i] <- t_faa[seqindex]
  toclust_COG_command <- paste("grep -l \"", Domain_annotation$Query_ID[i], "\" ",GETHOMS.outdir, "/", COG_GETHOMS.outdir,"/*.faa",sep = "")
  toclust_OMCL_command <- paste("grep -l \"", Domain_annotation$Query_ID[i], "\" ",GETHOMS.outdir, "/", OMCL_GETHOMS.outdir,"/*.faa",sep = "")
  toclust_COG_out <- gsub("/","",gsub(COG_GETHOMS.outdir,"",gsub(GETHOMS.outdir,"",system(toclust_COG_command, intern = TRUE))))
  toclust_OMCL_out <- gsub("/","",gsub(OMCL_GETHOMS.outdir,"",gsub(GETHOMS.outdir,"",system(toclust_OMCL_command, intern = TRUE))))
  Domain_annotation$COG_cluster[i] <- toclust_COG_out
  Domain_annotation$OMCL_cluster[i] <- toclust_OMCL_out
}


export1 <- Domain_annotation
export1name <- paste(GETHOMS.outdir,"_dannot.csv",sep="")

write.csv(export1,export1name)

##### EJEMPLO NANOVIRIDAE

#Functions <- c("Replication recombination repair","General function prediction","General function prediction","Unknown","Unknown","Capsid","Viral pathogenesis","Viral pathogenesis")
#Pg_clust <- c("Core","Cloud","Shell","Core","Cloud","Core","Core","Shell")
#Prev <- c(100,75,25,33,67,100,90,10)
#length(Prev)
#Nanoviridae <- data.frame(Funciones = as.factor(Functions),Grupo_Pangenomico = as.factor(Pg_clust),Porcentaje = Prev)
#library("ggplot2")
#x11()
#ggplot(Nanoviridae, aes(x=Funciones,y=Porcentaje,fill = Grupo_Pangenomico))+geom_bar(aes(width=.65),stat = "identity")

