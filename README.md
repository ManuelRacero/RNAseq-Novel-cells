RNA-seq_unassigned
Manuel Racero de la Rosa
2024-06-28
# Preparación del entorno
library(DESeq2)
## Loading required package: S4Vectors
## Loading required package: stats4
## Loading required package: BiocGenerics
## 
## Attaching package: 'BiocGenerics'
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
##     table, tapply, union, unique, unsplit, which.max, which.min
## 
## Attaching package: 'S4Vectors'
## The following object is masked from 'package:utils':
## 
##     findMatches
## The following objects are masked from 'package:base':
## 
##     expand.grid, I, unname
## Loading required package: IRanges
## 
## Attaching package: 'IRanges'
## The following object is masked from 'package:grDevices':
## 
##     windows
## Loading required package: GenomicRanges
## Loading required package: GenomeInfoDb
## Loading required package: SummarizedExperiment
## Loading required package: MatrixGenerics
## Loading required package: matrixStats
## 
## Attaching package: 'MatrixGenerics'
## The following objects are masked from 'package:matrixStats':
## 
##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
##     colWeightedMeans, colWeightedMedians, colWeightedSds,
##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
##     rowWeightedSds, rowWeightedVars
## Loading required package: Biobase
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
## 
## Attaching package: 'Biobase'
## The following object is masked from 'package:MatrixGenerics':
## 
##     rowMedians
## The following objects are masked from 'package:matrixStats':
## 
##     anyMissing, rowMedians
library(clusterProfiler)
## 
## clusterProfiler v4.10.1  For help: https://yulab-smu.top/biomedical-knowledge-mining-book/
## 
## If you use clusterProfiler in published research, please cite:
## T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu. clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. The Innovation. 2021, 2(3):100141
## 
## Attaching package: 'clusterProfiler'
## The following object is masked from 'package:IRanges':
## 
##     slice
## The following object is masked from 'package:S4Vectors':
## 
##     rename
## The following object is masked from 'package:stats':
## 
##     filter
library(org.Hs.eg.db)
## Loading required package: AnnotationDbi
## 
## Attaching package: 'AnnotationDbi'
## The following object is masked from 'package:clusterProfiler':
## 
##     select
## 
library(AnnotationDbi)

# Carga de datos
data <- read.csv("unassigned_expression_data.csv", row.names = 1)


# Transponer los datos para que los genes sean filas y las muestras sean columnas
data <- t(data)

# Verificar si hay valores no cero en los datos originales


# Asegurarse de que todos los valores sean enteros
data <- round(data)

data <- data + 1

# Filtrar genes con baja expresión
# Aquí se elimina cualquier gen que no tiene al menos 10 lecturas en al menos 2 muestras
keep <- rowSums(data >= 5) >= 2  # Al menos 10 lecturas en al menos 2 muestras
data <- data[keep, ]


# Asignar nombres a las columnas (opcional, si no están ya asignados)
colnames(data) <- paste0("sample", 1:ncol(data))
# Crear condiciones arbitrarias para las muestras
# Ajustar según el número de muestras que tienes
conditions <- factor(rep(c("control", "treated"), length.out = ncol(data)))

# Crear el DataFrame con las condiciones
colData <- DataFrame(condition = conditions)
# Verificar que el número de condiciones coincida con el número de columnas en data
if(ncol(data) != nrow(colData)) {
  stop("El número de columnas en 'data' no coincide con el número de filas en 'colData'")
}
# Crear el objeto DESeq2
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = colData,
                              design = ~ condition)
## converting counts to integer mode
# Realizar análisis de expresión diferencial
dds <- DESeq(dds)
## estimating size factors
## estimating dispersions
## gene-wise dispersion estimates
## mean-dispersion relationship
## -- note: fitType='parametric', but the dispersion trend was not well captured by the
##    function: y = a/x + b, and a local regression fit was automatically substituted.
##    specify fitType='local' or 'mean' to avoid this message next time.
## Warning in lfproc(x, y, weights = weights, cens = cens, base = base, geth =
## geth, : Estimated rdf < 1.0; not estimating variance
## final dispersion estimates
## fitting model and testing
res <- results(dds)
res
## log2 fold change (MLE): condition treated vs control 
## Wald test p-value: condition treated vs control 
## DataFrame with 691 rows and 6 columns
##          baseMean log2FoldChange     lfcSE      stat    pvalue      padj
##         <numeric>      <numeric> <numeric> <numeric> <numeric> <numeric>
## A2M       1.77080    -0.00583070 0.0582568 -0.100086  0.920276   0.99741
## ABCC9     1.26949     0.01031812 0.0694593  0.148549  0.881909   0.99741
## ACAA2     1.39813    -0.00703209 0.0654952 -0.107368  0.914497   0.99741
## ACTA2     1.13227     0.00774524 0.0717106  0.108007  0.913990   0.99741
## ACTB      4.55057    -0.01449539 0.0358058 -0.404834  0.685600   0.99741
## ...           ...            ...       ...       ...       ...       ...
## XRCC6     2.48660      0.0141373 0.0485794  0.291015  0.771040   0.99741
## YWHAB     2.07955     -0.0282949 0.0530670 -0.533193  0.593900   0.99741
## YWHAZ     2.99643      0.0122722 0.0442482  0.277348  0.781513   0.99741
## ZFP36     1.43083      0.0343230 0.0652205  0.526260  0.598707   0.99741
## ZFP36L1   1.85511     -0.0544784 0.0563914 -0.966077  0.334006   0.99741
hist(res$pvalue, breaks = 50, main = "Distribución de p-valores", xlab = "p-valores")


# Verificar el diseño experimental
head(colData)
## DataFrame with 6 rows and 1 column
##   condition
##    <factor>
## 1   control
## 2   treated
## 3   control
## 4   treated
## 5   control
## 6   treated
table(colData$condition)
## 
## control treated 
##     733     732
# Cargar las bibliotecas necesarias
library(clusterProfiler)
library(org.Hs.eg.db) # Usamos la base de datos de anotaciones para Homo sapiens
library(AnnotationDbi)

# Extraer los genes diferencialmente expresados con un ajuste de p-valor significativo
# Filtramos por ejemplo para un p-valor ajustado < 0.05
sig_genes <- subset(res, padj > 0.8)
gene_list <- rownames(sig_genes)

# Convertir los nombres de los genes a identificadores de Entrez
entrez_genes <- mapIds(org.Hs.eg.db, keys = gene_list, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
## 'select()' returned 1:1 mapping between keys and columns
# Filtrar los NA si es necesario
entrez_genes <- na.omit(entrez_genes)

# Análisis de enriquecimiento GO
go_enrich <- enrichGO(gene         = entrez_genes,
                      OrgDb        = org.Hs.eg.db,
                      keyType      = "ENTREZID",
                      ont          = "BP",  # Puedes cambiar a "MF" o "CC" para otras ontologías
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.2)

# Visualizar los resultados
dotplot(go_enrich, showCategory=20)
