# Introducción

Esta práctica introduce **MOFA2 (Multi-Omics Factor Analysis)**, una herramienta para integrar datos multi-ómicos y descubrir **factores latentes** que explican la variabilidad en distintos tipos de datos biológicos.

Analizaremos un dataset de pacientes con **leucemia linfocítica crónica (CLL)** que contiene:

- Respuesta a fármacos ex vivo  
- Metilación del ADN  
- RNA-seq normalizado  
- Mutaciones (presencia/ausencia)  

Se explorarán asociaciones de los factores con mutaciones clave, expresión génica, respuesta a fármacos y tiempo hasta un segundo tratamiento.

---

# Paquetes y datos

```{r paquetes-datos}
library(MOFA2)
library(MOFAdata)
library(data.table)
library(ggplot2)
library(tidyverse)
library(fitdistrplus)
library(pheatmap)
library(survival)
library(survminer)
library(mygene)

# Cargamos los datos
utils::data("CLL_data")
source("src/functions.R")

# Metadatos
CLL_metadata <- fread("ftp://ftp.ebi.ac.uk/pub/databases/mofa/cll_vignette/sample_metadata.txt")
head(CLL_metadata)
```

---

# 1. Preprocesado y exploración de datos

## 1.1 Respuesta a fármacos

```{r prepro-drugs}
drug <- CLL_data$Drugs[1,]
drug <- drug[!is.na(drug)]

# Histograma y QQ
plotdist(drug, histo = TRUE, demp = TRUE)
plotdist(log2(drug), histo = TRUE, demp = TRUE)

qqnorm(drug); qqline(drug, col="red")
qqnorm(log2(drug)); qqline(log2(drug), col="red")

# Prueba de normalidad
ks.test(drug, "pnorm", mean(drug), sd(drug))
ks.test(log2(drug), "pnorm", mean(log2(drug)), sd(log2(drug))

# Aplicamos log2
CLL_data$Drugs[!is.na(CLL_data$Drugs)] <- log2(CLL_data$Drugs[!is.na(CLL_data$Drugs)])
```

## 1.2 Metilación

```{r prepro-methylation}
cpg <- CLL_data$Methylation[1,]
cpg <- cpg[!is.na(cpg)]
ks.test(cpg, "pnorm", mean(cpg), sd(cpg))
```

## 1.3 RNA-seq

```{r prepro-rna}
gene <- CLL_data$mRNA[1,]
gene <- gene[!is.na(gene)]

plotdist(gene, histo=TRUE, demp=TRUE)
plotdist(log2(gene), histo=TRUE, demp=TRUE)

qqnorm(gene); qqline(gene, col="red")
qqnorm(log2(gene)); qqline(log2(gene), col="red")

ks.test(gene, "pnorm", mean(gene), sd(gene))
ks.test(log2(gene), "pnorm", mean(log2(gene)), sd(log2(gene))

CLL_data$mRNA[!is.na(CLL_data$mRNA)] <- log2(CLL_data$mRNA[!is.na(CLL_data$mRNA)] +1)
```

## 1.4 Mutaciones

```{r prepro-mutations}
# Ya en formato presencia/ausencia → Bernoulli
CLL_data$Mutations[1:10,1:10]
```

---

# 2. Creación y entrenamiento del modelo MOFA

```{r crear-modelo}
MOFAobject <- create_mofa(CLL_data)

plot_data_overview(MOFAobject)

data_opts <- get_default_data_options(MOFAobject)
model_opts <- get_default_model_options(MOFAobject)
model_opts$likelihoods["Mutations"] <- "bernoulli"
model_opts$num_factors <- 10
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "fast"
train_opts$seed <- 123

MOFAobject <- prepare_mofa(MOFAobject,
                           data_options = data_opts,
                           model_options = model_opts,
                           training_options = train_opts)

MOFAobject <- run_mofa(MOFAobject, outfile="MOFA2_CLL.hdf5", use_basilisk = TRUE)
samples_metadata(MOFAobject) <- CLL_metadata
```

---

# 3. Exploración del modelo

## 3.1 Correlación y varianza de factores

```{r exploracion-factores}
plot_factor_cor(MOFAobject)
plot_variance_explained(MOFAobject, max_r2 = 15)
plot_variance_explained(MOFAobject, plot_total = TRUE)[[2]]
```

## 3.2 Asociación con variables clínicas

```{r asociacion-clinica}
correlate_factors_with_covariates(MOFAobject, 
  covariates = c("Gender", "died", "trisomy12", "IGHV"), 
  plot = "log_pval")
```

---

# 4. Pesos de features y interpretación biológica

```{r pesos-features}
Z <- get_expectations(MOFAobject, variable = "Z")
W.mutations <- as.data.frame(get_expectations(MOFAobject, variable = "W")$Mutations)

plot_weights(MOFAobject, view="Mutations", factor=1, nfeatures=10, scale=TRUE)
plot_top_weights(MOFAobject, view="mRNA", factor=1, nfeatures=10)
```

## 4.1 Anotación de genes

```{r anotacion-genica}
MOFAobject_symbol <- MOFAobject
updated_features_names <- features_names(MOFAobject_symbol)

genesID <- mygene::queryMany(updated_features_names[["mRNA"]],
                             scopes = "ensembl.gene",
                             fields = "symbol",
                             species = "human")

genesID <- genesID[!duplicated(genesID$query), ]
updated_features_names[["mRNA"]] <- ifelse(is.na(genesID$symbol), genesID$query, genesID$symbol)
features_names(MOFAobject_symbol) <- updated_features_names

plot_top_weights(MOFAobject_symbol, view="mRNA", factor=1, nfeatures=10)
```

## 4.2 Enriquecimiento funcional

```{r enrichment}
utils::data("MSigDB_v6.0_C2_human")
enrichment.parametric <- run_enrichment(MOFAobject,
  view = "mRNA",
  factors = 1,
  feature.sets = MSigDB_v6.0_C2_human,
  sign = "positive",
  statistical.test = "parametric")

plot_enrichment(enrichment.parametric, factor=1, max.pathways=15)
```

---

# 5. Análisis clínico tiempo-dependiente

```{r survival-analysis}
SurvObject <- Surv(MOFAobject@samples_metadata$TTT, MOFAobject@samples_metadata$treatedAfter)
Z <- get_factors(MOFAobject)[[1]]
fit <- coxph(SurvObject ~ Z)
summary(fit)

# Forest plot
s <- summary(fit)
coef <- s[["coefficients"]]
df <- data.frame(
  factor = factor(rownames(coef), levels = rev(rownames(coef))),
  p      = coef[,"Pr(>|z|)"], 
  coef   = coef[,"exp(coef)"], 
  lower  = s[["conf.int"]][,"lower .95"], 
  higher = s[["conf.int"]][,"upper .95"]
)
ggplot(df, aes(x=factor, y=coef, ymin=lower, ymax=higher)) +
  geom_pointrange(col='#619CFF') + coord_flip() + geom_hline(aes(yintercept=1), linetype="dotted") +
  labs(y="Hazard Ratio", x="") + theme_bw()

# Kaplan-Meier por Factor 1
df <- data.frame(time = SurvObject[,1], event = SurvObject[,2], Z1 = Z[,1])
cut <- surv_cutpoint(df, variables='Z1')
df$FactorCluster <- df$Z1 > cut$cutpoint$cutpoint
fit <- survfit(Surv(time, event) ~ FactorCluster, df)
ggsurvplot(fit, data=df,
           conf.int = TRUE, pval = TRUE,
           fun = function(y) y * 100,
           legend = "top", legend.labs = c("low LF1", "high LF1"),
           xlab = "Time to treatment", ylab="Survival probability (%)",
           title= "Factor 1")
```

---

# Recursos oficiales

- MOFA2 Documentation: [https://biofam.github.io/MOFA2/](https://biofam.github.io/MOFA2/)  
- Repositorio MOFA: [https://github.com/bioFAM/MOFA](https://github.com/bioFAM/MOFA)  
- Datos de ejemplo (MOFAdata): [https://github.com/bioFAM/MOFAdata](https://github.com/bioFAM/MOFAdata)  
- Artículos científicos:  
  - Argelaguet et al., 2018 – Multi‑Omics Factor Analysis ([DOI](https://doi.org/10.15252/msb.20178124))  
  - Argelaguet et al., 2020 – MOFA+ ([DOI](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02015-1))  

---

# Conclusión

El análisis con MOFA2 permite:

- Integrar múltiples capas ómicas.  
- Identificar factores latentes explicativos.  
- Relacionar factores con mutaciones, expresión génica, respuesta a fármacos y supervivencia clínica.  
- Realizar interpretaciones funcionales mediante enriquecimiento.  

Esta práctica proporciona un flujo completo de **preprocesado → modelo → interpretación biológica y clínica**.