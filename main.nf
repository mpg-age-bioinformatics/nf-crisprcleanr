#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process get_images {
  stageInMode 'symlink'
  stageOutMode 'move'

  script:
    """

    if [[ "${params.containers}" == "singularity" ]] ; 

      then

        cd ${params.image_folder}

        if [[ ! -f crisprcleanr-3.0.0.sif ]] ;
          then
            singularity pull crisprcleanr-3.0.0.sif docker://index.docker.io/mpgagebioinformatics/crisprcleanr:3.0.0
        fi

    fi


    if [[ "${params.containers}" == "docker" ]] ; 

      then
        docker pull mpgagebioinformatics/crisprcleanr:3.0.0
    fi

    """

}

process subsetting_counts_file{
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val label
    val control
    val treatment

  output:
    val label

  when:
    ( ! file("${params.project_folder}/input/${label}.counts.txt").exists() ) 

  script:
  """
#!/usr/local/bin/python
import pandas as pd
import os

counts=pd.read_csv("${params.mageck_counts_file}", sep="\\t")
counts=counts.rename(columns={"Gene":"gene"})

control="${control}"
control_cols=control.split(",")

treatment="${treatment}"
treartment_cols=treatment.split(",")

counts_sub=counts[["sgRNA","gene"]+control_cols+treartment_cols]
counts_sub.to_csv("${params.project_folder}/input/${label}.counts.txt", sep="\\t", index=None, header=True)

  """

}

process cleanR_pipe {
  stageInMode 'symlink'
  stageOutMode 'move'

  input:
    val label
    val control
    val treatment

  when:
    ( ! file("${params.crisprcleanr_output}/${label}.gene_signatures.tsv").exists() )

  script:
  """
#!/usr/bin/env Rscript
library(CRISPRcleanR)

setwd("${params.project_folder}")

library_file <- "${params.lib_file_cleanR}"
fn <- "${params.counts_file_cleanR}"

#############################################
## Pipeline
#############################################

#"FUN": ["ccr.getLibrary"]
#"desc": ["Load Library Annotation"]
lib <- ccr.getLibrary(library_file = library_file, library_builtin = NULL)

#"FUN": ["ccr.getCounts"]
#"desc": ["Load Count data"]
counts <- ccr.getCounts(file_counts = fn, 
                    verbose = TRUE, 
                    files_FASTQ_controls=NULL, 
                    files_FASTQ_samples=NULL,
                    files_BAM_controls=NULL,
                    files_BAM_samples=NULL)

write.table(counts, file = "${params.crisprcleanr_output}/${label}.raw_counts.tsv", sep="\\t", quote = FALSE, row.names = FALSE)

#"FUN": ["ccr.checkCounts"]
#"desc": ["Check Library/Count data"]
####ccr.checkCounts(counts = counts, libraryAnnotation = lib, ncontrols = 2, min_reads = 30)
### add if statement here that if false do not proceed

#"FUN": ["ccr.NormfoldChanges"]
#"desc": ["Run Count normalization"]
normCountsAndFCs <- ccr.NormfoldChanges(filename = fn, 
                            saveToFig=TRUE,
                            display = TRUE,
                            method = 'ScalingByTotalReads', ## One of (ScalingByTotalReads, MedRatios)
                            ncontrols = ${params.ncontrols},
                            libraryAnnotation = lib,
                            EXPname = '${label}',
                            min_reads = 30,
                            outdir = "${params.crisprcleanr_output}")

write.table(normCountsAndFCs\$norm_counts, file="${params.crisprcleanr_output}/${label}.count_norm.tsv", sep="\\t", quote = FALSE, row.names = FALSE)
write.table(normCountsAndFCs\$logFCs, file="${params.crisprcleanr_output}/${label}.logFCs.tsv", sep="\\t", quote = FALSE, row.names = FALSE)


#### CN correction

#"FUN": ["ccr.logFCs2chromPos"]
#"desc": ["Run sgRNA Sorting"]
gwSortedFCs <- ccr.logFCs2chromPos(foldchanges = normCountsAndFCs\$logFCs, libraryAnnotation = lib)

# "FUN": ["ccr.GWclean"]
#"desc": ["Run sgRNA Sorting"]
correctedFCs <- ccr.GWclean(gwSortedFCs = gwSortedFCs ,label='${label}',display=TRUE,
                        saveTO="${params.crisprcleanr_output}",
                        ignoredGenes=NULL,
                        min.ngenes=3,
                        alpha = 0.01,
                        nperm = 10000,
                        p.method ="hybrid", ## one of (perm, hybrid)
                        min.width=2,
                        kmax=25,
                        nmin=200, 
                        eta=0.05,
                        trim = 0.025,
                        undo.splits = "none",
                        undo.prune=0.05, 
                        undo.SD=3,
                        return.segments.unadj = TRUE,
                        return.segments.adj = TRUE)

write.table(correctedFCs\$segments, file="${params.crisprcleanr_output}/${label}.segments.tsv", sep="\\t", quote = FALSE, row.names = FALSE)
write.table(correctedFCs\$segments_adj, file="${params.crisprcleanr_output}/${label}.segments_adj.tsv", sep="\\t", quote = FALSE, row.names = FALSE)
write.table(correctedFCs\$corrected_logFCs, "${params.crisprcleanr_output}/${label}.logFCs_adj.tsv", sep="\\t", quote = FALSE, row.names = FALSE)


#"FUN": ["ccr.correctCounts"]
#"desc": ["Correct counts"]
correctedCounts <- ccr.correctCounts(CL = '${label}',
                                normalised_counts = normCountsAndFCs\$norm_counts,
                                correctedFCs_and_segments = correctedFCs,
                                libraryAnnotation = lib,
                                minTargetedGenes=3,
                                OutDir="${params.crisprcleanr_output}",
                                ncontrols=${params.ncontrols})

write.table(correctedCounts, file="${params.crisprcleanr_output}/${label}.counts_corrected.tsv", sep="\\t", quote = FALSE, row.names = FALSE)



#"FUN": ["ccr.sgRNAmeanFCs"]
#"desc": ["Get corrected FCs by probe"]
sgRNAmeanFCs <- ccr.sgRNAmeanFCs(foldchanges = correctedFCs\$corrected_logFCs) 

#"FUN": ["ccr.geneMeanFCs"]
#"desc": ["Get corrected FCs by gene"]
geneMeanFCs <- ccr.geneMeanFCs(sgRNA_FCprofile = sgRNAmeanFCs, libraryAnnotation = lib)

#"FUN": ["ccr.geneSummary"]
#"desc": ["Get corrected FCs by gene"]
geneSummary <- ccr.geneSummary(sgRNA_FCprofile = sgRNAmeanFCs,
                        libraryAnnotation = lib,
                        FDRth=0.05)
write.table(geneSummary, file="${params.crisprcleanr_output}/${label}.gene_summary.tsv", sep="\\t", quote = FALSE, row.names = FALSE)

data("BAGEL_essential")
data("BAGEL_nonEssential")

BAGEL_essential_sgRNAs <- ccr.genes2sgRNAs(libraryAnnotation = lib, BAGEL_essential)
BAGEL_nonEssential_sgRNAs <- ccr.genes2sgRNAs(libraryAnnotation = lib,BAGEL_nonEssential)

#"FUN": ["ccr.ROC_Curve"]
#"desc": ["ROC curve at sgRNA level"]
pdf(paste0("${params.crisprcleanr_output}/${label}.ROC_by_sgRNA.pdf"), width = 10, height = 10)
sgRNA_level_ROC <- ccr.ROC_Curve(FCsprofile = sgRNAmeanFCs, 
                            positives = BAGEL_essential_sgRNAs,
                            negatives = BAGEL_nonEssential_sgRNAs,
                            display = TRUE,
                            FDRth = 0.05,
                            expName = '${label}')
dev.off()

#"FUN": ["ccr.ROC_Curve"]
#"desc": ["ROC curve at gene level"]
pdf(paste0("${params.crisprcleanr_output}/${label}.ROC_by_gene.pdf"), width = 10, height = 10)
gene_level_ROC <- ccr.ROC_Curve(FCsprofile = geneMeanFCs,
                        positives = BAGEL_essential,
                        negatives = BAGEL_nonEssential,
                        display = TRUE,
                        FDRth = 0.05,
                        expName = '${label}')
dev.off()                 

#"FUN": ["ccr.PrRc_Curve"]
#"desc": ["PrRc curve at sgRNA level"]
pdf(paste0("${params.crisprcleanr_output}/${label}.PrRc_by_sgRNA.pdf"), width = 10, height = 10)
sgRNA_level_PrRc <- ccr.PrRc_Curve(FCsprofile = sgRNAmeanFCs,
                            positives = BAGEL_essential_sgRNAs,
                            negatives = BAGEL_nonEssential_sgRNAs,
                            display = TRUE,
                            FDRth = 0.05,
                            expName = "${label}")
dev.off()

#"FUN": ["ccr.PrRc_Curve"]
#"desc": ["PrRc curve at gene level"]
pdf(paste0("${params.crisprcleanr_output}/${label}.PrRc_by_gene.pdf"), width = 10, height = 10)
gene_level_PrRc<-ccr.PrRc_Curve(FCsprofile = geneMeanFCs,
                            positives = BAGEL_essential,
                            negatives = BAGEL_nonEssential,
                            display = TRUE,
                            FDRth = 0.05,
                            expName = "${label}")
dev.off()


#"FUN": ["ccr.VisDepAndSig"]
#"desc": ["Visialize dependecy by gene signatures"]

data(EssGenes.ribosomalProteins)
data(EssGenes.DNA_REPLICATION_cons)
data(EssGenes.KEGG_rna_polymerase)
data(EssGenes.PROTEASOME_cons)
data(EssGenes.SPLICEOSOME_cons)

SIGNATURES <- list(Ribosomal_Proteins=EssGenes.ribosomalProteins,
            DNA_Replication = EssGenes.DNA_REPLICATION_cons,
            RNA_polymerase = EssGenes.KEGG_rna_polymerase,
            Proteasome = EssGenes.PROTEASOME_cons,
            Spliceosome = EssGenes.SPLICEOSOME_cons,
            CFE = BAGEL_essential,
            non_essential = BAGEL_nonEssential)

pdf(paste0("${params.crisprcleanr_output}/${label}.gene_signatures.pdf"), width = 10, height = 10)
Recall_scores <- ccr.VisDepAndSig(FCsprofile = geneMeanFCs,
                            SIGNATURES = SIGNATURES,
                            TITLE = '${params.label}',
                            pIs = 6,
                            nIs = 7,
                            th=0.05,
                            plotFCprofile=TRUE)
dev.off()
write.table(as.data.frame(Recall_scores), file="${params.crisprcleanr_output}/${label}.gene_signatures.tsv", sep="\\t", quote = FALSE, row.names = TRUE)



  """
}

workflow images {
  main:
    get_images()
}


workflow cleanR_workflow {
  if ( ! file("${params.crisprcleanr_output}").isDirectory() ) {
    file("${params.crisprcleanr_output}").mkdirs()
  }

  rows=Channel.fromPath("${params.samples_tsv}", checkIfExists:true).splitCsv(sep:';')
  rows=rows.filter{ ! file( "${params.crisprcleanr_output}//${it[0]}.gene_signatures.tsv" ).exists() }
  
  label=rows.flatMap { n -> n[0] }
  // paired=rows.flatMap { n -> n[1] }
  control=rows.flatMap { n -> n[2] }
  control=control.map{ "$it".replace(".fastq.gz","") }
  treatment=rows.flatMap { n -> n[3] }
  treatment=treatment.map{ "$it".replace(".fastq.gz","") }
  
  // subsetting_counts_file( label, control, treatment )
  // cleanR_pipe( subsetting_counts_file.out.collect() , control, treatment )

  cleanR_pipe( label , control, treatment )

}
