install.packages("BiocManager")

library(org.Mm.eg.db)
library(AnnotationDbi)

BiocManager::install("DESeq2")
library(DESeq2)
## 1) Listar arquivos de contagem
## Seus arquivos terminam com: _Star_Aligned.out.bam_all.txt
files <- list.files(".", pattern="Star_Aligned.out.bam_all.txt")

## Verificar se encontrou os 9 arquivos
print(files)
print(length(files))

## 4) Ler cada arquivo e extrair apenas a coluna de counts
## featureCounts gera 7 colunas:
## Geneid | Chr | Start | End | Strand | Length | counts
## A última coluna é sempre o count da amostra
cts <- lapply(files, function(f) {
  x <- read.table(f, header=TRUE, row.names=1)
  x[, ncol(x), drop=FALSE]   # mantém só a coluna de counts
})
####ver colunas
for (i in seq_along(cts)) {
  cat(files[i], ": ", ncol(cts[[i]]), " colunas\n", sep = "")
}
###ver nome das col;unas
for (i in seq_along(cts)) {
  cat(files[i], ": ", colnames(cts[[i]]), "\n", sep = "")
}

## 5) Combinar todas as colunas de counts em uma matriz
cts <- do.call(cbind, cts)

## 6) Nomear colunas com o nome da amostra (removendo o sufixo do arquivo)
colnames(cts) <- gsub("_Star_Aligned.out.bam_all.txt", "", files)

## Conferir estrutura final
colnames(cts)
head(cts)

###############################################
## 7) Criar colData automaticamente pelo nome ##
###############################################

sample_names <- colnames(cts)

## condition

condition <- ifelse(grepl("CTRL", sample_names, ignore.case=TRUE), "CTRL",
         ifelse(grepl("LOW",  sample_names, ignore.case=TRUE), "LOW",
         ifelse(grepl("HIGH", sample_names, ignore.case=TRUE), "HIGH", NA)))

condition
table(condition)

## tissue
tissue <- ifelse(grepl("Adipose", sample_names), "Adipose",
          ifelse(grepl("Duodenum", sample_names), "Duodenum",
          ifelse(grepl("Kidney", sample_names), "Kidney",
          ifelse(grepl("Liver", sample_names), "Liver",
          ifelse(grepl("Muscle", sample_names), "Muscle", NA)))))

tissue
table(tissue)

#type
type <- rep("paired-end", length(sample_names))
type <- factor(type)

#colData
colData <- data.frame(
  row.names = sample_names,
  condition = factor(condition, levels = c("CTRL","LOW","HIGH")),
  tissue    = factor(tissue),
  type = type
)

colData$condition <- factor(colData$condition)
colData$type <- factor(colData$type)


print(colData)
all(rownames(colData) == colnames(cts))
#########################################
## 8) Criar objeto DESeq2 e rodar DESeq ##
#########################################

## Aqui você NÃO define contrastes.
## Aqui você só ajusta o modelo estatístico.
## O DESeq() prepara tudo para os contrastes,
## mas não calcula log2FC ainda.

dds_global <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = colData,
  design = ~ condition
)

dds_global <- estimateSizeFactors(dds_global)

##Subset por tecido
dds_muscle <- dds_global[, dds_global$tissue == "Muscle"]

#Prefiltrar genes (≥10 counts em ≥3 amostras)
smallestGroupSize <- 3
keep <- rowSums(counts(dds_muscle) >= 10) >= smallestGroupSize
dds_muscle <- dds_muscle[keep, ]
dds_muscle

######Note on factor leves/Relevel#####
dds_muscle$condition <- relevel(dds_muscle$condition, ref = "CTRL")

#Rodar DESeq no subset
dds_muscle <- DESeq(dds_muscle)

dds_muscle
mcols(dds_muscle)

## 9) Extrair resultados dos dois contrastes do seu estudo ##
###########################################################

## Aqui SIM você define os contrastes explicitamente.

## LOW vs CTRL
res_low <- results(dds_muscle, contrast = c("condition", "LOW", "CTRL"))
res_low

## HIGH vs CTRL
res_high <- results(dds_muscle, contrast = c("condition", "HIGH", "CTRL"))
res_high

####Shrink#######################################################
install.packages("ashr")

resultsNames(dds_muscle)

res_high_lfc <- lfcShrink(
  dds_muscle,
  coef = "condition_HIGH_vs_CTRL",
  type = "ashr"
)
res_high_lfc

res_low_lfc <- lfcShrink(
  dds_muscle,
  coef = "condition_LOW_vs_CTRL",
  type = "ashr"
)
res_low_lfc

#Normalized counts
norm_counts <- counts(dds_muscle, normalized = TRUE)

norm_counts_df <- data.frame(
  gene = rownames(norm_counts),
  norm_counts
)

write.csv(
  norm_counts_df,
  file = "/group/sbms006/rrubimsilvadocouto/NovaSeq_Data/gene-counts-S637-S633/muscle/normalizedcounts_muscle.csv",
  row.names = FALSE
)


head(read.csv("/group/sbms006/rrubimsilvadocouto/NovaSeq_Data/gene-counts-S637-S633/muscle/normalizedcounts_muscle.csv"))

###############################################################
## 11) DIFERENTIAL com symbol, ensembl sem versão e User_ID com versão
###############################################################

norm <- counts(dds, normalized=TRUE)
norm_counts <- log2(norm + 1)

## Identificar genes
zero_all  <- apply(norm, 1, function(x) all(x == 0))
zero_some <- apply(norm, 1, function(x) any(x == 0))

## IDs
ensembl_with_version <- rownames(counts)
ensembl_no_version <- sub("\\..*", "", ensembl_with_version)

symbols <- mapIds(
  org.Mm.eg.db,
  keys = ensembl_no_version,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
symbols[is.na(symbols)] <- "NA"

## Parte principal
keep_main <- !zero_some

counts_main <- counts[keep_main, ]
norm_counts_main <- norm_counts[keep_main, ]
res_low_main <- res_low[keep_main, ]
res_high_main <- res_high[keep_main, ]

processed_main <- rep("", nrow(counts_main))
mean_expression_main <- rowMeans(norm_counts_main)

final_main <- data.frame(
  symbol = symbols[keep_main],
  ensembl_ID = ensembl_no_version[keep_main],
  User_ID = ensembl_with_version[keep_main],
  LOW_vs_CTRL_log2FC = res_low_main$log2FoldChange,
  LOW_vs_CTRL_adjPval = res_low_main$padj,
  HIGH_vs_CTRL_log2FC = res_high_main$log2FoldChange,
  HIGH_vs_CTRL_adjPval = res_high_main$padj,
  Processed_data = processed_main,
  norm_counts_main,
  mean_expression = mean_expression_main
)

final_main <- final_main[order(-final_main$mean_expression), ]
final_main$mean_expression <- NULL

## Parte intermediária
keep_mid <- zero_some & !zero_all

counts_mid <- counts[keep_mid, ]
norm_counts_mid <- norm_counts[keep_mid, ]
res_low_mid <- res_low[keep_mid, ]
res_high_mid <- res_high[keep_mid, ]

processed_mid <- rep("", nrow(counts_mid))

final_mid <- data.frame(
  symbol = symbols[keep_mid],
  ensembl_ID = ensembl_no_version[keep_mid],
  User_ID = ensembl_with_version[keep_mid],
  LOW_vs_CTRL_log2FC = res_low_mid$log2FoldChange,
  LOW_vs_CTRL_adjPval = res_low_mid$padj,
  HIGH_vs_CTRL_log2FC = res_high_mid$log2FoldChange,
  HIGH_vs_CTRL_adjPval = res_high_mid$padj,
  Processed_data = processed_mid,
  norm_counts_mid
)

## Parte final (zero em todas)
keep_low <- zero_all

counts_low <- counts[keep_low, ]
norm_counts_low <- norm_counts[keep_low, ]
res_low_low <- res_low[keep_low, ]
res_high_low <- res_high[keep_low, ]

processed_low <- rep("", nrow(counts_low))

final_low <- data.frame(
  symbol = symbols[keep_low],
  ensembl_ID = ensembl_no_version[keep_low],
  User_ID = ensembl_with_version[keep_low],
  LOW_vs_CTRL_log2FC = res_low_low$log2FoldChange,
  LOW_vs_CTRL_adjPval = res_low_low$padj,
  HIGH_vs_CTRL_log2FC = res_high_low$log2FoldChange,
  HIGH_vs_CTRL_adjPval = res_high_low$padj,
  Processed_data = processed_low,
  norm_counts_low
)

## Juntar tudo
final <- rbind(final_main, final_mid, final_low)

## Salvar
write.csv(
  final,
  "/group/sbms006/rrubimsilvadocouto/NovaSeq_Data/genecount_all/adipose-tissue_all/Adipose-Tissue_DIFFERENTIAL.csv",
  row.names = FALSE
)

###############################################################
## 12) Criar segunda aba "REGULATION" igual ao artigo
##     mantendo symbol, ensembl_ID e User_ID
###############################################################

## Usar o mesmo data.frame final da aba DIFFERENTIAL
df <- final

## Extrair log2FC e padj
lfc_low  <- df$LOW_vs_CTRL_log2FC
padj_low <- df$LOW_vs_CTRL_adjPval

lfc_high  <- df$HIGH_vs_CTRL_log2FC
padj_high <- df$HIGH_vs_CTRL_adjPval

## Função atualizada
regulation_call <- function(lfc, padj) {
  if (is.na(padj)) {
    return("None")
  } else if (padj >= 0.05) {
    return("None")
  } else if (abs(lfc) < 1) {   # <---- filtro novo
    return("None")
  } else if (lfc > 0) {
    return("Up")
  } else if (lfc < 0) {
    return("Down")
  } else {
    return("None")
  }
}

## Aplicar a função para cada gene
reg_low  <- mapply(regulation_call, lfc_low,  padj_low)
reg_high <- mapply(regulation_call, lfc_high, padj_high)

## Construir a aba REGULATION com as colunas adicionais
regulation <- data.frame(
  symbol = df$symbol,
  ensembl_ID = df$ensembl_ID,
  User_ID = df$User_ID,
  `20 mg CBD/kg b.w. vs. CTR` = reg_low,
  `100 mg CBD/kg b.w. vs. CTR` = reg_high
)

## Salvar a aba REGULATION
write.csv(
  regulation,
  "/group/sbms006/rrubimsilvadocouto/NovaSeq_Data/genecount_all/muscle_all/Muscle_REGULATION_log2.csv",
  row.names = FALSE
)



###############################################################
## 13) Criar terceira aba "VENN" com symbol, ensembl_ID e User_ID
###############################################################

df <- final   ## mesma ordem da aba DIFFERENTIAL

## IDs
ids_no_version <- df$ensembl_ID
ids_with_version <- df$User_ID
symbols <- df$symbol

## Recalcular classificação diretamente do final
lfc_low  <- df$LOW_vs_CTRL_log2FC
padj_low <- df$LOW_vs_CTRL_adjPval

lfc_high  <- df$HIGH_vs_CTRL_log2FC
padj_high <- df$HIGH_vs_CTRL_adjPval

## Função atualizada novamente
regulation_call <- function(lfc, padj) {
  if (is.na(padj)) {
    return("None")
  } else if (padj >= 0.05) {
    return("None")
  } else if (abs(lfc) < 1) {   # <---- filtro novo
    return("None")
  } else if (lfc > 0) {
    return("Up")
  } else if (lfc < 0) {
    return("Down")
  } else {
    return("None")
  }
}

reg_low  <- mapply(regulation_call, lfc_low,  padj_low)
reg_high <- mapply(regulation_call, lfc_high, padj_high)

## Conjuntos
LOW_DOWN  <- ids_no_version[reg_low  == "Down"]
LOW_UP    <- ids_no_version[reg_low  == "Up"]
HIGH_DOWN <- ids_no_version[reg_high == "Down"]
HIGH_UP   <- ids_no_version[reg_high == "Up"]

## Interseções
LOW_DOWN_HIGH_DOWN <- intersect(LOW_DOWN, HIGH_DOWN)
LOW_UP_HIGH_UP     <- intersect(LOW_UP, HIGH_UP)

## Função para criar bloco longo com todas as colunas
make_block <- function(name, elements) {

  if (length(elements) == 0) {
    return(data.frame(
      Names = name,
      total = 0,
      symbol = NA,
      ensembl_ID = NA,
      User_ID = NA
    ))
  }

  idx <- match(elements, ids_no_version)

  data.frame(
    Names = name,
    total = length(elements),
    symbol = symbols[idx],
    ensembl_ID = ids_no_version[idx],
    User_ID = ids_with_version[idx]
  )
}

## Criar blocos
block1 <- make_block("20mg_DOWN 100mg_DOWN", LOW_DOWN_HIGH_DOWN)
block2 <- make_block("20mg_UP 100mg_UP", LOW_UP_HIGH_UP)
block3 <- make_block("20mg_DOWN", LOW_DOWN)
block4 <- make_block("20mg_UP", LOW_UP)
block5 <- make_block("100mg_DOWN", HIGH_DOWN)
block6 <- make_block("100mg_UP", HIGH_UP)

## Juntar tudo
venn_df <- rbind(block1, block2, block3, block4, block5, block6)

## Salvar
write.csv(
  venn_df,
  "/group/sbms006/rrubimsilvadocouto/NovaSeq_Data/genecount_all/muscle_all/Muscle_VENN_log2.csv",
  row.names = FALSE
)

########################## VENN DIAGRAM #####################

install.packages("VennDiagram")
library(VennDiagram)

## Conjuntos
A <- LOW_DOWN      # 20mg_DOWN
B <- LOW_UP        # 20mg_UP
C <- HIGH_DOWN     # 100mg_DOWN
D <- HIGH_UP       # 100mg_UP

## Criar arquivo PNG
png(
  "/group/sbms006/rrubimsilvadocouto/NovaSeq_Data/genecount_all/muscle_all/Muscle_VENN_plot_log2.png",
  width = 6000,
  height = 4000,
  res = 300
)

draw.quad.venn(
  area1 = length(A),
  area2 = length(B),
  area3 = length(C),
  area4 = length(D),

  ## interseções duplas
  n12 = length(intersect(A, B)),
  n13 = length(intersect(A, C)),
  n14 = length(intersect(A, D)),
  n23 = length(intersect(B, C)),
  n24 = length(intersect(B, D)),
  n34 = length(intersect(C, D)),

  ## interseções triplas (todas zero)
  n123 = 0,
  n124 = 0,
  n134 = 0,
  n234 = 0,

  ## interseção quádrupla (zero)
  n1234 = 0,

  category = c("20mg_DOWN", "20mg_UP", "100mg_DOWN", "100mg_UP"),

  fill = c("dodgerblue", "red", "forestgreen", "gold"),
  alpha = 0.5,
  lwd = 3,
  cex = 2.5,
  cat.cex = 2.5,
  cat.col = c("dodgerblue", "red", "forestgreen", "gold")
)

dev.off()

###############################################
############## BAR PLOT#######################
###############################################

library(ggplot2)

df <- final   # sua tabela DIFFERENTIAL já carregada

## Extrair vetores
lfc_low  <- df$LOW_vs_CTRL_log2FC
padj_low <- df$LOW_vs_CTRL_adjPval

lfc_high  <- df$HIGH_vs_CTRL_log2FC
padj_high <- df$HIGH_vs_CTRL_adjPval

## Função atualizada novamente
regulation_call <- function(lfc, padj) {
  if (is.na(padj)) {
    return("None")
  } else if (padj >= 0.05) {
    return("None")
  } else if (abs(lfc) < 1) {   # <---- filtro novo
    return("None")
  } else if (lfc > 0) {
    return("Up")
  } else if (lfc < 0) {
    return("Down")
  } else {
    return("None")
  }
}

reg_low  <- mapply(regulation_call, lfc_low,  padj_low)
reg_high <- mapply(regulation_call, lfc_high, padj_high)

## Contagens
LOW_UP    <- sum(reg_low  == "Up")
LOW_DOWN  <- sum(reg_low  == "Down")
HIGH_UP   <- sum(reg_high == "Up")
HIGH_DOWN <- sum(reg_high == "Down")

## Criar data frame para o plot
plot_df <- data.frame(
  Treatment = c("100 mg/kg UP", "100 mg/kg DOWN",
                "20 mg/kg UP",  "20 mg/kg DOWN"),
  Count = c(HIGH_UP, HIGH_DOWN, LOW_UP, LOW_DOWN),
  Type = c("UP", "DOWN", "UP", "DOWN")
)

## Cores
color_up   <- "#C2185B"   # dark pink/red
color_down <- "#5C6BC0"   # light blue/purple

plot_df$Color <- ifelse(plot_df$Type == "UP", color_up, color_down)

############# Bar plot########################

png(
  "/group/sbms006/rrubimsilvadocouto/NovaSeq_Data/genecount_all/muscle_all/Muscle_BarPlot_log2.png",
  width = 3000,
  height = 2500,
  res = 300
)

ggplot(plot_df, aes(x = Treatment, y = Count, fill = Color)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = Count), vjust = -0.5, size = 6) +
  scale_fill_identity() +
  theme_minimal(base_size = 16) +
  labs(
    title = "Number of Differentially Expressed Genes",
    x = "",
    y = "Number of Genes"
  ) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 14, angle = 20, hjust = 1)
  )

dev.off()

##################################
## PCA PLOT #####################
#################################

library(ggplot2)

## Transformação VST (segura para PCA)
vsd <- vst(dds, blind = FALSE)

## PCA
pca <- prcomp(t(assay(vsd)))

## Data frame
pca_df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  group = colData$group
)

## Renomear grupos para legenda
pca_df$group <- factor(
  pca_df$group,
  levels = c("CTRL", "LOW", "HIGH"),
  labels = c("CTRL", "20 mg/kg b.w.", "100 mg/kg b.w.")
)

## Cores novas
color_ctrl <- "grey50"      # cinza médio
color_low  <- "#4CAF50"     # verde planta
color_high <- "#D4A017"     # mostarda

## Plot PCA
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 6) +
  theme_minimal(base_size = 16) +
  labs(
    title = "Muscle PCA",
    x = paste0("PC1 (", round(summary(pca)$importance[2,1] * 100, 1), "% variance)"),
    y = paste0("PC2 (", round(summary(pca)$importance[2,2] * 100, 1), "% variance)")
  ) +
  scale_color_manual(values = c(
    "CTRL" = color_ctrl,
    "20 mg/kg b.w." = color_low,
    "100 mg/kg b.w." = color_high
  )) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 14)
  )


png(
  "/group/sbms006/rrubimsilvadocouto/NovaSeq_Data/genecount_all/muscle_all/Muscle_PCA.png",
  width = 3000,
  height = 2500,
  res = 300
)

print(p)

dev.off()


####################################################################
##### HIGH VOLCANO PLOT ##################
###############################################################
library(ggplot2)

df <- final

lfc  <- df$HIGH_vs_CTRL_log2FC
padj <- df$HIGH_vs_CTRL_adjPval

## Remover NA
keep <- !is.na(lfc) & !is.na(padj)
lfc  <- lfc[keep]
padj <- padj[keep]

## Remover outliers absurdos (caso existam)
lfc <- lfc[lfc > -50 & lfc < 50]
padj <- padj[lfc > -50 & lfc < 50]

neglog10 <- -log10(padj)

## Classificação
reg <- ifelse(padj < 0.05 & lfc > 1, "Up-regulated",
       ifelse(padj < 0.05 & lfc < -1, "Down-regulated",
              "Not significant"))

volcano_df <- data.frame(
  lfc = lfc,
  neglog10 = neglog10,
  reg = reg
)

## Cores
color_up   <- "#C2185B"
color_down <- "#5C6BC0"
color_none <- "grey60"

p <- ggplot(volcano_df, aes(x = lfc, y = neglog10, color = reg)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c(
    "Up-regulated" = color_up,
    "Down-regulated" = color_down,
    "Not significant" = color_none
  )) +
  theme_minimal(base_size = 16) +
  labs(
    title = "Muscle Volcano Plot",
    subtitle = "CBD 100 mg/kg b.w. × CTRL",
    x = "log2 Fold Change",
    y = "-log10(padj)"
  ) +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 16, hjust = 0.5),
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 14)
  )
############SALVAR A IMAGEM###########
png(
  "/group/sbms006/rrubimsilvadocouto/NovaSeq_Data/genecount_all/muscle_all/Muscle_HighVolcanoPlot.png",
  width = 3000,
  height = 2500,
  res = 300
)
print(p)
dev.off()
#########################LOW VOLCANO PLOT ########################
library(ggplot2)

## Use SEMPRE o objeto res_low do DESeq2
df <- as.data.frame(res_low)

## Remover NA corretamente
df <- df[!is.na(df$log2FoldChange) & !is.na(df$padj), ]

## Criar colunas
df$lfc      <- df$log2FoldChange
df$padj     <- df$padj
df$neglog10 <- -log10(df$padj)

## Classificação
df$reg <- ifelse(
  df$padj < 0.05 & df$lfc > 0, "Up-regulated",
  ifelse(df$padj < 0.05 & df$lfc < 0, "Down-regulated",
         "Not significant")
)

## FORÇAR categorias a existirem na legenda
df$reg <- factor(df$reg,
                 levels = c("Up-regulated", "Down-regulated", "Not significant"))

## Cores
colors <- c(
  "Up-regulated"   = "#C2185B",
  "Down-regulated" = "#5C6BC0",
  "Not significant" = "grey60"
)

## Volcano plot
p <- ggplot(df, aes(x = lfc, y = neglog10, color = reg)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = colors, drop = FALSE) +   # <--- mantém azul na legenda
  theme_minimal(base_size = 16) +
  labs(
    title = "Volcano Plot — LOW vs CTRL",
    subtitle = "CBD 20 mg/kg b.w. × CTRL",
    x = "log2 Fold Change",
    y = "-log10(padj)"
  ) +
  coord_cartesian(
    ylim = c(0, max(df$neglog10, na.rm = TRUE))          # <--- mostra padj=1 (y=0)
  ) +
  theme(
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 16, hjust = 0.5),
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 14)
  )

## Salvar
png(
  "/group/sbms006/rrubimsilvadocouto/NovaSeq_Data/genecount_all/adipose-tissue_all/Adipose-Tissue_LowVolcanoPlot.png",
  width = 3000,
  height = 2500,
  res = 300
)
print(p)
dev.off()









