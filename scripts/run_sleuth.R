args <- commandArgs(trailingOnly = TRUE) #read command line arguments from snakemake
output_file <- args[1] #first argument is output path where sleuth results go

#make sure directory for output file exists, otherwise make it
if (!dir.exists(dirname(output_file))) {
    dir.create(dirname(output_file), recursive = TRUE)
}

library(sleuth) #load sleuth

#create sample/condition dataframe for experimental data
s2c <- data.frame(
sample=c(
"SRR5660030",
"SRR5660044",
"SRR5660033",
"SRR5660045"
),
condition=c(
"dpi2",
"dpi2",
"dpi6",
"dpi6"
),

#paths to kallisto directories
path=c(
"output/kallisto/SRR5660030",
"output/kallisto/SRR5660044",
"output/kallisto/SRR5660033",
"output/kallisto/SRR5660045"
),

#prevent R from converting strings to factors
stringsAsFactors=FALSE
)

#convert condition to factor
s2c$condition <- factor(s2c$condition)

#need to run in tryCatch, otherwise snakemake pipeline crashes if anything fails
tryCatch({
    #create sleuth object
    so <- sleuth_prep(
        s2c,
        ~condition,
        extra_bootstrap_summary = TRUE
        )
    #fit linear model to expression data
    so <- sleuth_fit(so)
    #wald test comparing 6 and 2 dpi
    so <- sleuth_wt(so,"conditiondpi6")
    #get differential expression results with stats for each transcript
    res <- sleuth_results(so,"conditiondpi6")

    #compute test stat bc it doesnt already exist in res??
    res$test_stat <- res$b / res$se_b

    print(colnames(res))

    #grab columns we want for output
    sig <- res[,c("target_id", "test_stat", "pval", "qval")]
    #filter for only significant results
    sig_results <- sig[!is.na(sig$qval) & sig$qval < 0.05, ]

    #write significant results to output
    write.table(
        sig_results,
        output_file,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE
        )
    }, 
    #handle errors
    error = function(e) {
    #if it fails put error message in log file
    writeLines(as.character(e), "sleuth_error_log.txt")
    #stop running so we know it failed
    stop(e)
})

