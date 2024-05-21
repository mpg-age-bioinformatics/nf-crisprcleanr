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

process cleanR_pipe {
  stageInMode 'symlink'
  stageOutMode 'move'

//   output:
//     val "/workdir/deseq2_output/deseq2.part1.Rdata", emit: rdata

  when:
    ( ! file("${params.project_folder}/CRISPRcleanR_correctedCounts.RData").exists() ) 

  script:
  """
#!/usr/bin/Rscript
library(CRISPRcleanR)
library(xlsx)

setwd("/nexus/posix0/MAGE-flaski/service/posit/home/aiqbal/count/")

library_file <- "${params.lib_file_cleanR}"
###fn <- file.path(getwd(), "JEKO.txt", fsep="/")
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

#"FUN": ["ccr.checkCounts"]
#"desc": ["Check Library/Count data"]
ccr.checkCounts(counts = counts, libraryAnnotation = lib, ncontrols = 2, min_reads = 30)

### add if statement here that if false do not proceed

#"FUN": ["ccr.NormfoldChanges"]
#"desc": ["Run Count normalization"]
normCountsAndFCs <- ccr.NormfoldChanges(filename = fn, 
                            display = TRUE,  ##for pipe set saveToFig=TRUE
                            method = 'ScalingByTotalReads', ## One of (ScalingByTotalReads, MedRatios)
                            ncontrols = 1,
                            libraryAnnotation = lib,
                            EXPname = 'CRISPRcleanR',
                            min_reads = 30,
                            outdir = "${params.crisprcleanr_output}")

#### CN correction

#"FUN": ["ccr.logFCs2chromPos"]
#"desc": ["Run sgRNA Sorting"]
gwSortedFCs <- ccr.logFCs2chromPos(foldchanges = normCountsAndFCs\$logFCs, libraryAnnotation = lib)

# "FUN": ["ccr.GWclean"]
#"desc": ["Run sgRNA Sorting"]
correctedFCs <- ccr.GWclean(gwSortedFCs = gwSortedFCs ,label='CRISPRcleanR',display=TRUE,
                        saveTO="${params.crisprcleanr_output}", ### set path in pipe
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
                        undo.SD=3)

#"FUN": ["ccr.correctCounts"]
#"desc": ["Correct counts"]
correctedCounts <- ccr.correctCounts(CL = 'CRISPRcleanR',
                                normalised_counts = normCountsAndFCs\$norm_counts,
                                correctedFCs_and_segments = correctedFCs,
                                libraryAnnotation = lib,
                                minTargetedGenes=3,
                                OutDir="${params.crisprcleanr_output}",
                                ncontrols=1)

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

data("BAGEL_essential")
data("BAGEL_nonEssential")

BAGEL_essential_sgRNAs <- ccr.genes2sgRNAs(libraryAnnotation = lib, BAGEL_essential)
BAGEL_nonEssential_sgRNAs <- ccr.genes2sgRNAs(libraryAnnotation = lib,BAGEL_nonEssential)

#"FUN": ["ccr.ROC_Curve"]
#"desc": ["ROC curve at sgRNA level"]
sgRNA_level_ROC <- ccr.ROC_Curve(FCsprofile = sgRNAmeanFCs, 
                            positives = BAGEL_essential_sgRNAs,
                            negatives = BAGEL_nonEssential_sgRNAs,
                            display = TRUE,
                            FDRth = 0.05,
                            expName = 'CRISPRcleanR')


#"FUN": ["ccr.ROC_Curve"]
#"desc": ["ROC curve at gene level"]
gene_level_ROC <- ccr.ROC_Curve(FCsprofile = geneMeanFCs,
                        positives = BAGEL_essential,
                        negatives = BAGEL_nonEssential,
                        display = TRUE,
                        FDRth = 0.05,
                        expName = 'CRISPRcleanR')
                        

#"FUN": ["ccr.PrRc_Curve"]
#"desc": ["PrRc curve at sgRNA level"]
sgRNA_level_PrRc <- ccr.PrRc_Curve(FCsprofile = sgRNAmeanFCs,
                            positives = BAGEL_essential_sgRNAs,
                            negatives = BAGEL_nonEssential_sgRNAs,
                            display = TRUE,
                            FDRth = 0.05,
                            expName = "CRISPRcleanR")

#"FUN": ["ccr.PrRc_Curve"]
#"desc": ["PrRc curve at gene level"]
gene_level_PrRc<-ccr.PrRc_Curve(FCsprofile = geneMeanFCs,
                            positives = BAGEL_essential,
                            negatives = BAGEL_nonEssential,
                            display = TRUE,
                            FDRth = 0.05,
                            expName = "CRISPRcleanR")

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

Recall_scores <- ccr.VisDepAndSig(FCsprofile = geneMeanFCs,
                            SIGNATURES = SIGNATURES,
                            TITLE = 'CRISPRcleanR',
                            pIs = 6,
                            nIs = 7,
                            th=0.05,
                            plotFCprofile=TRUE)


  """
}

workflow cleanR_workflow {
    cleanR_pipe()
}