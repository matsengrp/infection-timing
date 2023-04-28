args =  commandArgs(TRUE)
filepath=args[1]
outdir=args[2]
if (file.exists(filepath) == FALSE || length(args) < 2){
  writeLines ("Usage:\nRscript AAVF_3rdcodon_function.R filepath outdir\n");
  quit()
}

library(stringr)
library(readr)
library(stringr)
library(magrittr)
library(tidyverse)
library(Biostrings)
library(seqinr)

if(str_detect(basename(filepath),"F1") == TRUE){
	frag="F1"
}
if(str_detect(basename(filepath),"F2") == TRUE){
        frag="F2"
}
if(str_detect(basename(filepath),"F3") == TRUE){
        frag="F3"
}

codon_variant_frac_min=0.01
coverage_min=5000
fragment=frag

  my_output_vector=vector()
  AAVF_metadata=read_lines(filepath, n_max=10)
  AAVF_raw=read.table(filepath, sep="\t", header=FALSE)
  AAVF_names=unlist(strsplit(AAVF_metadata[10], "\\t"))
  AAVF_names[1]="CHROM"
  
  
  #pull out the name of the sample from line 4
  sample_name=unlist(strsplit((unlist(strsplit(AAVF_metadata[4], "[.]"))[1]),"="))[2]
  my_output_vector[1]=sample_name
  
  
  colnames(AAVF_raw)=AAVF_names
  AAVF_raw$POS=as.numeric(AAVF_raw$POS)
  AAVF_raw$COVERAGE=as.numeric(AAVF_raw$COVERAGE)
  
  #For each AA, this lists the number of codons found. There are up to 6 codons per AA (leucine, arginine) so we'll split this column into 6. Unused columns are filled with NA to the right.
  AAVF_raw%<>%separate(INFO, c("AC", "ACC", "ACF"), sep=";")
  AAVF_raw=separate(AAVF_raw, AC, c("AC1", "AC2", "AC3", "AC4", "AC5", "AC6"), sep=",", fill="right")
  Variants=separate(AAVF_raw, ACC, c("ACC1", "ACC2", "ACC3", "ACC4", "ACC5", "ACC6"), sep="," , fill="right")
  Variants=separate(Variants, ACF, c("ACF1", "ACF2", "ACF3", "ACF4", "ACF5", "ACF6"), sep="," , fill="right")
  #need to remove strings "AC=", "ACC=" "ACF=" from the first set of columns
  Variants$AC1=str_replace(Variants$AC1, "AC=", "")
  Variants$ACC1=str_replace(Variants$ACC1, "ACC=", "")
  Variants$ACF1=str_replace(Variants$ACF1, "ACF=", "")
  
  
  Variants$ACC1%<>%as.numeric()
  Variants$ACC2%<>%as.numeric()
  Variants$ACC3%<>%as.numeric()
  Variants$ACC4%<>%as.numeric()
  Variants$ACC5%<>%as.numeric()
  Variants$ACC6%<>%as.numeric()
  
  Variants%<>% mutate(Syn_Coverage=rowSums(Variants[,c("ACC1", "ACC2", "ACC3", "ACC4", "ACC5", "ACC6")],  na.rm=T))
  
  
  Variants%<>%unite(Codon_1, c(AC1, ACC1, ACF1), sep = "_", remove = TRUE)
  Variants%<>%unite(Codon_2, c(AC2, ACC2, ACF2), sep = "_", remove = TRUE)
  Variants%<>%unite(Codon_3, c(AC3, ACC3, ACF3), sep = "_", remove = TRUE)
  Variants%<>%unite(Codon_4, c(AC4, ACC4, ACF4), sep = "_", remove = TRUE)
  Variants%<>%unite(Codon_5, c(AC5, ACC5, ACF5), sep = "_", remove = TRUE)
  Variants%<>%unite(Codon_6, c(AC6, ACC6, ACF6), sep = "_", remove = TRUE)
  
  Variants %<>% gather(key=codon_number, value=codon_info, Codon_1, Codon_2, Codon_3, Codon_4, Codon_5, Codon_6)
  
  Variants=separate(Variants, codon_info, c("Codon", "Codon_coverage", "Codon_frequency"), sep="_")
  Variants=Variants[which(Variants$Codon!="NA"),]
  
  #Codon_frequency=how often a codon was seen out of the total number of reads that covered that position
  
  Variants$Codon_coverage=as.numeric(Variants$Codon_coverage)
  Variants%<>% group_by(POS, ALT) %>% mutate(dominant_allele_coverage=max(Codon_coverage,na.rm=T))
  #here I want to get a nucleotide consensus sequence--will use this for making a blast database to check for cross-contamination
  #I'll use the ref AA and then the most common codon for that AA
  Variants$REF=as.character(Variants$REF)
  Variants$ALT=as.character(Variants$ALT)
  get_nt_consensus=Variants[which(Variants$REF==Variants$ALT),]
  get_nt_consensus=get_nt_consensus[which(get_nt_consensus$Codon_coverage==get_nt_consensus$dominant_allele_coverage),]
  #make sure it's in order of position
  get_nt_consensus=get_nt_consensus[order(get_nt_consensus$POS),]
  #squish them all into one string
  nt_consensus=paste(get_nt_consensus$Codon, collapse="")
  
  #now to look at third codon positions only--separate our codons into nucleotides
  Variants=separate(Variants, Codon, c("First", "Second", "Third"), sep=c(1,2))
  Variants$Codon_coverage=as.numeric(Variants$Codon_coverage)
  Variants$Codon_frequency=as.numeric(Variants$Codon_frequency)

  #looking at the third codon in a position, no matter what the amino acid is; 
  #for example, reads of aaa and cta will be counted together
  Variants%<>% group_by(POS, Third) %>% mutate(third_nuc_sum=sum(Codon_coverage, na.rm=T))
  
  Variants%<>%mutate(freq_in_syn=Codon_coverage/Syn_Coverage)
  Variants%<>%mutate(freq_syn_dom=dominant_allele_coverage/Syn_Coverage)
  
  #Select only those that pass the filter
  Pass_vars=Variants[which(Variants$FILTER=="PASS"),]
  Pass_vars=Variants[which(Variants$COVERAGE>=coverage_min),]
  Pass_vars%<>% group_by(POS) %>% mutate(Pass_denom=sum(Codon_coverage, na.rm=T))
  
  #just want one row for each nuc at the third position
  #note this part just ignores the amino acids--it's as if we never translated it 
  #and only knew the nucleotides
  third_pass_vars=select(Pass_vars, POS,COVERAGE,Third, third_nuc_sum, Pass_denom)
  third_pass_vars=unique(third_pass_vars)
  third_pass_vars %<>% mutate(third_freq=third_nuc_sum/Pass_denom)
  third_pass_vars %<>% group_by(POS) %>% mutate(freq_third_dom=max(third_freq))
  #use this to ignore sites where diversity doesn't meet our minimum variant calling threshhold
  third_pass_vars%<>%mutate(to_use=ifelse(1-freq_third_dom-codon_variant_frac_min > 0, 1, 0 )) 
  third_pass_vars%<>%mutate(contribution=third_freq*(1-third_freq)*to_use)
  #Now only those that are synonymous
  #need to know number of sites to divde by
  #note this number is the number of amino acids, so equal to the number of 3rd codon positions
  covered_sites=length(unique(third_pass_vars$POS))
  APD=sum(third_pass_vars$contribution)/covered_sites
  
  #my_output_vector is what the function will return when you run it. 
  #Any info we'll want in the end needs to go in there.
  
  my_output_vector[2]=APD
  

  #the coverage column lists overall sequencing depth at each position
  #We want to just check quickly that that is reasonable--but it's repeated if there are different AAs at the same pos, so we just want unique position/coverage values
  POS_COVERAGE=unique(AAVF_raw[,c("POS", "COVERAGE", "REF")])
  POS_COVERAGE%<>%arrange(POS)
  POS_COVERAGE$REF=as.character(POS_COVERAGE$REF)
  #plot coverage by positionls
  
  covplot=ggplot(POS_COVERAGE, aes(y=COVERAGE, x=POS)) + 
    geom_bar(position="dodge", stat="identity") + 
    ggtitle(sample_name) +
    geom_hline(yintercept=coverage_min, color="red", linetype="dashed") +
    theme_bw()
  covplot
  consensus=paste(POS_COVERAGE$REF, collapse="")

  good_consensus=paste(POS_COVERAGE$REF[(which(POS_COVERAGE$COVERAGE>coverage_min))], collapse="")
  
  #make sure our sequence of amino acids with significant coverage is continuous, and extract 1st and last codons with good covereage, relative to 
  #sample specific reference
  check_continuous=function (POS_COVERAGE, coverage_min) {
    goodAAs=POS_COVERAGE$POS[(which(POS_COVERAGE$COVERAGE>coverage_min))]
    start_pos=min(goodAAs)
    print (paste0("starting position ", start_pos))
    end_pos=max(goodAAs)
    print (paste0("ending position ", end_pos))
    print(paste("There are",  length(goodAAs), "amino acids with >5000x coverage for this sample"))
    is_cont=ifelse(end_pos-start_pos+1==length(goodAAs), "Continuous sequence", "Not continuous sequence" )
    seq_ok=ifelse(end_pos-start_pos+1==length(goodAAs), 1, 0)
    print(is_cont)
    return_vector=c(start_pos, end_pos, length(goodAAs), seq_ok)
  }
  
  
  check_vector=check_continuous(POS_COVERAGE, coverage_min)

  
  #translate those start and end codons to the corresponding HXB2 nucleotide start and stop positions
  #this is important for plugging them into Neher's clock, found at:
  #https://hiv.biozentrum.unibas.ch/ETI/
  
  #now to figure out where our sample specific consensus lives relative to HXB2. This was tricky, I suspect I'm missing something easy.
  #but, I aligned them
  #found where the aligned portion of HXB2 sits in the HXB2 translated pol gene 
  #(use str split, size of the first split is how many AAs skipped before alignment starts.)
  #take that length +1 to get first aa, multiply by 3 to get nt position, then add nt position of beginning of gene, subtracting 1..
  
  
  
  get_start_HXB2_AA_coord=function(consensus, frag){
    if (frag=="F1"){
      #this name might not be right!
      myRef=read.fasta(file="gag_HXB2_aa.fasta")
    } else {
      myRef=read.fasta(file="HXB2_pol_aa.fasta")
    }
    consens_string=toupper(c2s(consensus[[1]]))
    myref_string=toupper(c2s(myRef[[1]]))
    myalign=pairwiseAlignment(consens_string, myref_string, type="local")
    myalign
    search_string=as.character(subject(myalign))
    HXB2aastart=nchar(str_split(myref_string, search_string)[[1]][1])
    return(HXB2aastart)
  }
  
  HXB2aastart=get_start_HXB2_AA_coord(good_consensus, fragment)
  
  #get hxb2 gene start site depending on fragment
  gene_start=ifelse(fragment=="F1", 790, 2085)
  
  HXB_nt_start=(HXB2aastart+1)*3+gene_start-1
  HXB_nt_end=HXB_nt_start + check_vector[3]*3
  print(HXB_nt_start)
  print(HXB_nt_end)
  print(APD)
  my_output_vector=c(my_output_vector, HXB2aastart, HXB2aastart+check_vector[3])
  my_output_vector=c(my_output_vector, check_vector[3:4])
  my_output_vector=c(my_output_vector, HXB_nt_start, HXB_nt_end, nt_consensus, fragment)
  write.table(t(as.data.frame(my_output_vector)),file=paste0(outdir,sample_name,"_APD.txt"),sep="\t",row.names=FALSE)
  #plot coverage by positionls
  covplot=ggplot(POS_COVERAGE, aes(y=COVERAGE, x=POS)) + 
    geom_bar(position="dodge", stat="identity") +  
    theme_bw() + geom_hline(yintercept = 5000, color="red", linetype="dashed")
  ggsave(paste0(outdir,sample_name,"_cov.png"), height=8, width=12, units="in")
  


