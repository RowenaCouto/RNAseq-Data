###############################################
## DESeq2 Pipeline — Multi‑Tissue RNA‑seq
## Author: Rowena
## Description:
##   - Reads featureCounts outputs
##   - Builds count matrix
##   - Creates colData automatically
##   - Runs DESeq2
##   - Saves LOW vs CTRL, HIGH vs CTRL, HIGH vs LOW
###############################################

library(DESeq2)

message("=== Starting DESeq2 pipeline ===")

###############################################
## 1. Detect files in current directory
###############################################

files <- list.files(pattern = "txt$", full.names = TRUE)

if (length(files) != 9) {
  stop("ERROR: Expected 9 count files (CTRL1‑3, HIGH1‑3, LOW1‑3). Found: ", length(files))
}

message("Found 9 count files.")

###############################################
## 2. Extract sample names
###############################################

sampleNames <- gsub("_Star_Aligned.out.bam_all.txt", "", basename(files))
message("Sample names detected:")
print(sampleNames)

###############################################
## 3. Read each file and extract Geneid + Count
###############################################

count_list <- lapply(seq_along(files), function(i) {
  df <- read.table(files[i], header = FALSE, comment.char = "#")
  df <- df[-1, ]  # remove header line
  colnames(df) <- c("Geneid","Chr","Start","End","Strand","Length","Count")
  df <- df[, c("Geneid","Count")]
  colnames(df)[2] <- sampleNames[i]
  df
})

###############################################
## 4. Merge all count tables
###############################################

counts <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE), count_list)

rownames(counts) <- counts$Geneid
counts$Geneid <- NULL

###############################################
## 5. Convert counts to numeric
###############################################

counts[] <- lapply(counts, function(x) as.numeric(as.character(x)))

message("Count matrix dimensions:")
print(dim(counts))

###############################################
## 6. Build colData automatically
###############################################

condition <- c(
  "CTRL","CTRL","CTRL",
  "HIGH","HIGH","HIGH",
  "LOW","LOW","LOW"
)

colData <- data.frame(
  row.names = colnames(counts),
  condition = factor(condition)
)

message("colData:")
print(colData)

###############################################
## 7. Run DESeq2
###############################################

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = colData,
  design = ~ condition
)

dds <- DESeq(dds)

###############################################
## 8. Extract contrasts
###############################################

res_low      <- results(dds, contrast = c("condition","LOW","CTRL"))
res_high     <- results(dds, contrast = c("condition","HIGH","CTRL"))
res_high_low <- results(dds, contrast = c("condition","HIGH","LOW"))

###############################################
## 9. Save results
###############################################

tissue <- basename(getwd())

write.csv(as.data.frame(res_low),
          paste0("DESeq2_", tissue, "_LOW_vs_CTRL.csv"),
          row.names = TRUE)

write.csv(as.data.frame(res_high),
          paste0("DESeq2_", tissue, "_HIGH_vs_CTRL.csv"),
          row.names = TRUE)

write.csv(as.data.frame(res_high_low),
          paste0("DESeq2_", tissue, "_HIGH_vs_LOW.csv"),
          row.names = TRUE)

message("=== DESeq2 pipeline completed successfully for tissue: ", tissue, " ===")

