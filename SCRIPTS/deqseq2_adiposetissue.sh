## 1) Listar arquivos de contagem

## Seus arquivos terminam com: _Star_Aligned.out.bam_all.txt
files <- list.files(".", pattern="_all.txt")

## Verificar se encontrou os 9 arquivos
print(files)
print(length(files))

#### Checar número de linhas e colunas de todos os arquivos

files <- list.files(pattern = "_all.txt$")

dims <- lapply(files, function(f) {
  df <- read.table(f, header = TRUE)
  c(linhas = nrow(df), colunas = ncol(df))
})

names(dims) <- files
dims

##################################################
#######Ver a primeira linha de todos os arquivos

files <- list.files(pattern = "_all.txt$")

primeiras <- sapply(files, function(f) readLines(f, n = 1))

data.frame(
  arquivo = files,
  primeira_linha = primeiras,
  row.names = NULL
)

##################################################
####### remove a primeira linha dos arquivos ######
files <- list.files(pattern = "_all.txt$")

for (f in files) {
  linhas <- readLines(f)
  writeLines(linhas[-1], f)
}

#######Ver a primeira linha de todos os arquivos

files <- list.files(pattern = "_all.txt$")

primeiras <- sapply(files, function(f) readLines(f, n = 1))

data.frame(
  arquivo = files,
  primeira_linha = primeiras,
  row.names = NULL
)
#### Checar número de linhas de cada arquivo ###

num_linhas <- sapply(files, function(f) length(readLines(f)))

data.frame(
  arquivo = files,
  linhas = num_linhas,
  row.names = NULL
)

#### Checar número de colunas de cada arquivo
num_colunas <- sapply(files, function(f) {
  ncol(read.table(f, header = TRUE, comment.char = "#", nrows = 1))
})

data.frame(
  arquivo = files,
  colunas = num_colunas,
  row.names = NULL
)
###### Manter apenas Geneid + última coluna (counts)

for (f in files) {
  df <- read.table(f, header = FALSE)
  df2 <- df[, c(1, ncol(df))]   # coluna 1 = geneid, última coluna = counts
  write.table(df2, f, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

########Checar número de colunas novamente
num_colunas <- sapply(files, function(f) {
  ncol(read.table(f, header = FALSE))
})

data.frame(
  arquivo = files,
  colunas = num_colunas,
  row.names = NULL
)

########Checar número de linhas novamente
num_linhas <- sapply(files, function(f) {
  nrow(read.table(f, header = FALSE))
})

data.frame(
  arquivo = files,
  linhas = num_linhas,
  row.names = NULL
)

