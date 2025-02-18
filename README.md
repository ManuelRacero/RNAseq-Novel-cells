---
title: "RNA-seq_unassigned"
author: "Manuel Racero de la Rosa"
date: "2024-06-28"
output: html_document
---

```{r, warning=FALSE}

# Preparación del entorno
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
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


```

```{r}
# Crear condiciones arbitrarias para las muestras
# Ajustar según el número de muestras que tienes
conditions <- factor(rep(c("control", "treated"), length.out = ncol(data)))

# Crear el DataFrame con las condiciones
colData <- DataFrame(condition = conditions)
```

```{r}
# Verificar que el número de condiciones coincida con el número de columnas en data
if(ncol(data) != nrow(colData)) {
  stop("El número de columnas en 'data' no coincide con el número de filas en 'colData'")
}
```


```{r}

# Crear el objeto DESeq2
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = colData,
                              design = ~ condition)

# Realizar análisis de expresión diferencial
dds <- DESeq(dds)
res <- results(dds)
```
```{r}
res
```

```{r}
hist(res$pvalue, breaks = 50, main = "Distribución de p-valores", xlab = "p-valores")
```

```{r}
# Verificar el diseño experimental
head(colData)
table(colData$condition)


```


```{r}
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

```
