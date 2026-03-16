


# Cargamos los paquetes necesarios para la práctica

library(MOFA2)
library(MOFAdata)
library(data.table)
library(ggplot2)
library(tidyverse)
library(fitdistrplus)
library(pheatmap)
library(survival)
library(survminer)
# Cargamos los datos 
utils::data("CLL_data")       
source("src/functions.R")
## 1. Cargar, visualizar y transformar los datos

# Echamos un vistazo al objeto,
# que está formado por una lista de matrices.

# Respuesta a drogas ex vivo
CLL_data$Drugs[1:10,1:10]

drug <- CLL_data$Drugs[1,]
drug <- drug[!is.na(drug)]

# MOFA funciona mejor con distribuciones gaussianas.
# En ocasiones, es necesario transformar los datos
# utilizando el logaritmo para ajustarlos a una distribución normal.

# Estudiamos el histograma antes y después de aplicar el logaritmo
plotdist(drug, histo = TRUE, demp = TRUE)
plotdist(log2(drug), histo = TRUE, demp = TRUE)

# Estudiamos el gráfico QQ para ver en qué caso se ajusta mejor a una normal
qqnorm(drug)
qqline(drug, col = "red")

qqnorm(log2(drug))
qqline(log2(drug), col = "red")

# Por último, podemos usar la prueba de Kolmogorov-Smirnov
# para evaluar si la distribución se aproxima a una normal.
# La hipótesis nula es que sigue una distribución normal,
# por lo tanto, valores p mayores a 0.05 indicarían que no se rechaza la hipótesis nula
# y, por lo tanto, que la distribución es normal.

ks.test(drug, "pnorm", mean(drug), sd(drug))
ks.test(log2(drug), "pnorm", mean(log2(drug)), sd(log2(drug)))

# Vemos que la transformación logarítmica mejora los resultados,
# por lo tanto, modificamos los datos.

CLL_data$Drugs[!is.na(CLL_data$Drugs)] <- log2(CLL_data$Drugs[!is.na(CLL_data$Drugs)])

# Metilación

# La metilación debe utilizarse en valores M
# para que siga una distribución lo más parecida posible a la gaussiana.

CLL_data$Methylation[1:10,1:10]

cpg <- CLL_data$Methylation[1,]
cpg <- cpg[!is.na(cpg)]

# Verificamos que no sigue una normal, pero es la distribución más parecida.
ks.test(cpg, "pnorm", mean(cpg), sd(cpg))

# RNA

# Vamos a averiguar la distribución que tiene el
# RNA-seq normalizado para decidir qué verosimilitud
# utilizar. Seguimos la misma estrategia que con las drogas.

CLL_data$mRNA[1:10,1:10]

gene <- CLL_data$mRNA[1,]
gene <- gene[!is.na(gene)]

plotdist(gene, histo = TRUE, demp = TRUE)
plotdist(log2(gene), histo = TRUE, demp = TRUE)

qqnorm(gene)
qqline(gene, col = "red")

qqnorm(log(gene))
qqline(log(gene), col = "red")

ks.test(gene, "pnorm", mean(gene), sd(gene))
ks.test(log(gene), "pnorm", mean(log(gene)), sd(log(gene)))

CLL_data$mRNA[!is.na(CLL_data$mRNA)] <- log2(CLL_data$mRNA[!is.na(CLL_data$mRNA)] +1)

# Mutaciones
CLL_data$Mutations[1:10,1:10]

# Las mutaciones están en formato presencia/ausencia,
# por lo que siguen una distribución de Bernoulli.

# Cargamos los metadatos
CLL_metadata <- fread("ftp://ftp.ebi.ac.uk/pub/databases/mofa/cll_vignette/sample_metadata.txt")

head(CLL_metadata)

# Creamos el objeto MOFA
MOFAobject <- create_mofa(CLL_data)

# Revisamos el conjunto de nuestros datos: número de vistas
# y si hay datos desconocidos.
plot_data_overview(MOFAobject)

#################################################
#################################################
# 2. Generar el modelo

data_opts <- get_default_data_options(MOFAobject)
data_opts

# Dejamos las opciones de los datos como están.

# El escalado de las vistas permite ajustar
# el hecho de que una vista aporte mucha más información que otra.
# Sin embargo, hay que tener cuidado con la interpretación,
# ya que esto fuerza a que las vistas sean comparables.

model_opts <- get_default_model_options(MOFAobject)
model_opts

# En primer lugar, cambiamos la distribución de
# las mutaciones a Bernoulli (presencia/ausencia).
model_opts$likelihoods["Mutations"] <- "bernoulli"
model_opts$num_factors <- 10

train_opts <- get_default_training_options(MOFAobject)

# Opciones de entrenamiento:
# En esta práctica hemos seleccionado "fast"
# para reducir el tiempo de ejecución,
# pero para un modelo final se recomienda "slow"
# para obtener mayor precisión.

# stochastic = TRUE puede ser útil y necesario
# para conjuntos de datos grandes, como en single-cell.

train_opts$convergence_mode <- "fast"
train_opts$seed <- 123

# Entrenamos el modelo
MOFAobject <- prepare_mofa(MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

# Generamos el modelo con las opciones escogidas
MOFAobject <- run_mofa(MOFAobject, outfile="MOFA2_CLL.hdf5", use_basilisk = TRUE)

# Añadimos los metadatos
samples_metadata(MOFAobject) <- CLL_metadata

#####################################################################################
# 3. Investigar el modelo

# Buscamos factores que expliquen variabilidades únicas
# y que no deben correlacionarse entre sí.
plot_factor_cor(MOFAobject)

# Estudio de la varianza explicada por cada factor en cada vista.

plot_variance_explained(MOFAobject, max_r2 = 15)

# El Factor 1 captura una fuente de variabilidad presente en todas las modalidades de datos.
# Por lo tanto, su etiología probablemente sea clave para la enfermedad.

# Los Factores 3 y 4 capturan variabilidad en fármacos y mRNA.

# El Factor 2 captura variabilidad en múltiples modalidades de datos,
# aunque en menor medida en la metilación del ADN.
# Esto también puede ser importante.

# Analizamos qué vistas están más explicadas por el modelo.
plot_variance_explained(MOFAobject, plot_total = TRUE)[[2]]

##############################################################################
# Asociación entre factores y variables clínicas con MOFA

correlate_factors_with_covariates(MOFAobject, 
  covariates = c("Gender", "died", "trisomy12", "IGHV"), 
  plot = "log_pval"
)

# Si no convencen los métodos utilizados por los autores,
# se pueden extraer los pesos y realizar los análisis con código propio.

Z <- get_expectations(MOFAobject, variable = "Z") 

meta <- MOFAobject@samples_metadata
meta <- cbind(meta, Z)

# Se escriben funciones para calcular la matriz de p-valores usando ANOVA.

# Se extrae la matriz de p-valores    
p.values <- matrix.pvalues(data = meta,
  variables_of_interest = c("Gender", "died", "trisomy12", "IGHV"))

# Se ajustan los p-valores por pruebas múltiples.
p_adjusted_matrix <- apply(p.values, 2, function(x) p.adjust(x, method = "fdr"))

# Se transforman los valores a -log10(p)
p_transformed_matrix <- -log10(p_adjusted_matrix)
# Convertimos la matriz a un formato largo
p_matrix_long <- as.data.frame(as.table(p_transformed_matrix)) %>%
  rename(Row = Var1, Column = Var2, p_value = Freq) %>% 
  mutate(
    Row = gsub("group1\\.", "", Row),
    Column = gsub("\\.p.value", "", Column)
  ) 

# Establecemos un umbral para los valores de p transformados
p_matrix_long$p_value[p_matrix_long$p_value < -log10(0.05)] <- 0

p_matrix_long$Row <- factor(p_matrix_long$Row, levels = c("Factor1","Factor2","Factor3",
"Factor4","Factor5","Factor6","Factor7","Factor8","Factor9","Factor10"))
# Creamos el heatmap
ggplot(p_matrix_long, aes(x = Column, y = Row, fill = p_value)) +
  geom_tile(color = "black", lwd = 0.5) +
  scale_fill_gradientn(colors = c("grey", "white", "red"), values = c(0, 0.000001, 1)) +
  theme_minimal() +
  labs(title = "Heatmap de p-valores", x = "Columnas", y = "Filas") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

######################################################################################
# Una vez que identificamos qué variables están relacionadas con qué factores,
# podemos representarlo con diferentes gráficos.

# En este caso, el Factor 1 y la mutación en IGHV muestran una clara asociación.

plot_factor(MOFAobject, 
  factors = 1, 
  color_by = "IGHV"
)

plot_factor(MOFAobject, 
  factors = 1, 
  color_by = "IGHV",
  add_violin = TRUE,
  dodge = TRUE
)

# Verificamos en la matriz de mutaciones si IGHV realmente está explicando
# variabilidad en la vista de mutaciones.

W.mutations <- as.data.frame(get_expectations(MOFAobject, variable = "W")$Mutations)

head(W.mutations[order(abs(W.mutations$Factor1), decreasing = TRUE),])

# Podemos usar las funciones de MOFA para representarlo gráficamente.

plot_weights(MOFAobject,
  view = "Mutations",
  factor = 1,
  nfeatures = 10,     # Número de características principales a resaltar
  scale = TRUE        # Escala los pesos de -1 a 1
)

plot_top_weights(MOFAobject,
  view = "Mutations",
  factor = 1,
  nfeatures = 10,
  scale = TRUE
)

# Ahora sabemos que el Factor 1 diferencia las muestras según la mutación en IGHV.
# Queremos analizar cómo se relaciona esto con las demás capas.

plot_top_weights(MOFAobject, 
  view = "mRNA", 
  factor = 1, 
  nfeatures = 10
)

# Observamos que hay un gen asociado a este factor en el polo positivo
# y otro grupo de genes en el otro sentido.
# Es decir, hay un gen más expresado en las muestras mutadas
# y otros genes más expresados en las muestras no mutadas.

## Anotamos los genes para mostrar sus símbolos y facilitar la interpretación biológica.

MOFAobject_symbol <- MOFAobject

updated_features_names <- features_names(MOFAobject_symbol)

genesID <- mygene::queryMany(updated_features_names[["mRNA"]], 
                             scopes = "ensembl.gene", 
                             fields = "symbol", 
                             species = "human")

genesID <- genesID[!duplicated(genesID$query), ]
updated_features_names[["mRNA"]] <- ifelse(is.na(genesID$symbol), genesID$query, genesID$symbol)
features_names(MOFAobject_symbol) <- updated_features_names

plot_top_weights(MOFAobject_symbol, 
  view = "mRNA", 
  factor = 1, 
  nfeatures = 10
)

# Vemos que la mutación en IGHV (cadena de inmunoglobulina pesada)
# podría estar afectando el transcriptoma, provocando:
# - Mayor expresión de ADAM29

# Ahora analizamos su impacto en la respuesta a fármacos:

plot_weights(MOFAobject, 
  view = "Drugs", 
  factor = 1, 
  nfeatures = 10
)

# Vemos que la respuesta a varios fármacos está afectada por esta relación.

# Para realizar un análisis de enriquecimiento con el método de los autores,
# cargamos los conjuntos de genes predefinidos en el paquete.
# Si quisiéramos usar otros, tendríamos que cargarlos y ajustarlos 
# al mismo formato (genes en columnas, 1 si pertenece a la ruta, 0 si no).

utils::data("MSigDB_v6.0_C2_human")

enrichment.parametric <- run_enrichment(MOFAobject,
  view = "mRNA", 
  factors = 1,
  feature.sets = MSigDB_v6.0_C2_human,
  sign = "positive",
  statistical.test = "parametric"
)

plot_enrichment(enrichment.parametric, 
  factor = 1, 
  max.pathways = 15
)

####################################################################
#Por último vamos a hacer un análisis tiempo dependiente tratando de relacionar
# Los factores con el riesgo de necesitar un segundo tratamiento 

SurvObject <- Surv(MOFAobject@samples_metadata$TTT, MOFAobject@samples_metadata$treatedAfter)
Z <- get_factors(MOFAobject)[[1]]
fit <- coxph(SurvObject ~ Z) 
fit

#Vemos varios factores asociados 
# Sin embargo, el factor 1 es el que se encuentra mas asociado donde cuando mayor sea el factor 1
# menor riesgo de tener que recibir un segundo tratamiento
#Y sabemos que la mutación en IGHV y la mayor expresión de ADAM29 se asocia respuesta de varios fármacos
# Tambien lo hace en los tratamientos de los pacientes


#Podemos ver los resultados de la regresiñon de cox en un forestplot
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
  geom_pointrange( col='#619CFF') + 
  coord_flip() +
  scale_x_discrete() + 
  labs(y="Hazard Ratio", x="") + 
  geom_hline(aes(yintercept=1), linetype="dotted") +
  theme_bw()

# Podemos dividir los pacienets en dos grupos en base al factor 1
# Y hacer un análisis de Kaplan Meier
  df <- data.frame(
  time = SurvObject[,1], 
  event = SurvObject[,2], Z1 = Z[,1]
)
cut <- surv_cutpoint(df, variables='Z1')
df$FactorCluster <- df$Z1 > cut$cutpoint$cutpoint
fit <- survfit(Surv(time, event) ~ FactorCluster, df)

ggsurvplot(fit, data = df,
  conf.int = TRUE, pval = TRUE,
  fun = function(y) y * 100,
  legend = "top", legend.labs = c(paste("low LF 1"), paste("high LF 1")),
  xlab = "Time to treatment", ylab="Survival probability (%)", title= "Factor 1"
)$plot