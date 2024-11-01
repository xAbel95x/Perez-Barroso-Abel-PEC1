paquetes <- c("readr", "Biobase", "SummarizedExperiment", "dplyr", "ggplot2", "tidyr", "grDevices") # Instalar librerias
for (paquete in paquetes) {
  if (!requireNamespace(paquete, quietly = TRUE)) {
    install.packages(paquete, dependencies = TRUE)
  }
  library(paquete, character.only = TRUE)
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("Biobase", "SummarizedExperiment"))
library(readr)
library(Biobase)
library(SummarizedExperiment)
library(dplyr)
library(grDevices)
library(ggplot2)
library(tidyr)
setwd("C:/Users/Administrador/Desktop/UOC/Repositorio Git/metaboData/Datasets/2024-Cachexia")
human_cachexia <- read_csv("human_cachexia.csv")


human_cachexia_df <- as.data.frame(human_cachexia) # Creación SummarizedExperiment
metabolite_data <- as.matrix(human_cachexia_df[, -c(1, 2)]) 
rownames(metabolite_data) <- human_cachexia_df$`Patient ID`
metabolite_data <- t(metabolite_data)
colData <- DataFrame(
  Patient_ID = human_cachexia_df$`Patient ID`,
  Muscle_loss = human_cachexia_df$`Muscle loss`)
se <- SummarizedExperiment(
  assays = list(counts = metabolite_data),
  colData = colData,
)

colSums(is.na(human_cachexia_df)) # Mirar si hay NA
str(human_cachexia_df) # Variables
summary(human_cachexia_df) # Resumen datos
datos_todo <- human_cachexia_df[, -c(1, 2)] # Quitar ID y muscle loss
datos_cachexia <- subset(human_cachexia_df, `Muscle loss` == "cachexic")[, -c(1, 2)] # Grupo caquexia
datos_control <- subset(human_cachexia_df, `Muscle loss` == "control")[, -c(1, 2)] # Grupo control
f<- function(x){ # Función para histogramas
  ifelse (is.numeric(x), 
          hist(x, breaks=5),
          barplot(table(x))
  )
}
par(mfrow=c(3,2))
apply(human_cachexia_df,2,f) # Gráficos sin normalizar para control calidad
apply(datos_todo,2,f)
apply(datos_cachexia,2,f)
apply(datos_control,2,f)
Muscle_loss <- ifelse(human_cachexia_df$`Muscle loss` == "cachexic", 0, 1) # Grupos pasan a número 0 y 1
hist(Muscle_loss, # Histograma grupos
     main = "Histograma Pacientes",  
     xlab = "Grupos",               
     ylab = "Pacientes",
     breaks = 2)
hist(datos_cachexia$Creatinine, # Histogramas nivel creatinina
     main = "Histograma de Creatinina Caquexia", 
     xlab = "Nivel de Creatinina",        
     ylab = "Frecuencia",                  
     breaks = 10)
hist(datos_control$Creatinine, 
     main = "Histograma de Creatinina Control", 
     xlab = "Nivel de Creatinina",        
     ylab = "Frecuencia",                  
     breaks = 10)
datos_todo_norm <- scale(datos_todo) # Normalizar datos
datos_todo_normalizados <- cbind(human_cachexia_df[, c(1, 2)], datos_todo_norm)
datos_cachexia_norm <- scale(datos_cachexia)
datos_control_norm <- scale(datos_control)
datos_cachexia_norm <- as.data.frame(datos_cachexia_norm)
par(mfrow=c(3,2))
apply(datos_todo_norm,2,f) # Gráficos normalizados para control calidad
apply(datos_todo_normalizados,2,f)
apply(datos_cachexia_norm,2,f)
apply(datos_control_norm,2,f)
hist(datos_cachexia_norm$Creatinine, # Histogramas normalizado nivel creatinina
     main = "Histograma de Creatinina Normalizado Caquexia", 
     xlab = "Nivel de Creatinina",        
     ylab = "Frecuencia")
datos_control_norm <- as.data.frame(datos_control_norm)
hist(datos_control_norm$Creatinine, 
     main = "Histograma de Creatinina Normalizado Control", 
     xlab = "Nivel de Creatinina",        
     ylab = "Frecuencia",                  
     breaks = 10)
paleta_cor <- colorRampPalette(c("blue", "white", "red"))(256) # Paleta para heatmap
cor_todo <- cor(datos_todo_norm) # Correlación
heatmap(cor_todo, main = "Matriz de Correlación de Metabolitos", # Gráficos heatmap
        Colv = NA, Rowv = NA, scale = "none", margins = c(5, 5),col = paleta_cor)
cor_cachexia <- cor(datos_cachexia_norm)
heatmap(cor_cachexia, main = "Matriz de Correlación de Metabolitos Cachexia", 
        Colv = NA, Rowv = NA, scale = "none", margins = c(5, 5), col = paleta_cor)
cor_control <- cor(datos_control_norm)
heatmap(cor_control, main = "Matriz de Correlación de Metabolitos Control", 
        Colv = NA, Rowv = NA, scale = "none", margins = c(5, 5), col = paleta_cor)
pca_todo <- prcomp(datos_todo_norm) # PCA
summary(pca_todo) # Resumen PCA
pca_todo$sdev
pca_todo$rotation[,1]
pca_todo$rotation[order(pca_todo$rotation[,1], decreasing = TRUE), 1] # Ordenar primera componente
pca_todo$rotation[order(pca_todo$rotation[,2], decreasing = TRUE), 1] # Ordenar segunda componente
pca_scores <- pca_todo$x
head(pca_todo$x)
varianza_acumulada <- cumsum(pca_todo$sdev^2 / sum(pca_todo$sdev^2)) * 100 # Varianza acumulada
print(varianza_acumulada)
barplot(pca_todo$sdev^2 / sum(pca_todo$sdev^2) * 100, # Gráfico barras de varianza
        main = "Varianza Componentes",
        xlab = "Componente Principal",
        ylab = "Varianza (%)",
        col = "lightblue")
var_todo <- pca_todo$sdev^2 / sum(pca_todo$sdev^2) * 100
plot(pca_scores[, 1:2], # Gráfico de PCA puntos
     col = ifelse(human_cachexia_df$`Muscle loss` == "cachexic", "red", "green"),
     pch = 19, 
     xlab = paste0("Componente Principal 1 (", round(var_todo[1], 2), "%)"), 
     ylab = paste0("Componente Principal 2 (", round(var_todo[2], 2), "%)"),
     main = "PCA: Primeros Dos Componentes Principales")
legend("topright", legend = c("Cachexia", "Control"), col = c("red", "green"), pch = 19) # Leyenda
totales_cachexia <- colSums(datos_cachexia) # Suma columnas
totales_control <- colSums(datos_control)
head(sort(totales_cachexia, decreasing = TRUE)) # Concentraciones más altas
head(sort(totales_control, decreasing = TRUE)) 
comparacion_metabolitos <- data.frame( # Data frame para comparar concentraciones
  Metabolito = names(totales_cachexia),
  Cachexia = totales_cachexia,
  Control = totales_control
)
comparacion_larga <- comparacion_metabolitos %>% # Reorganizar datos
  pivot_longer(cols = -Metabolito, 
               names_to = "Grupo", 
               values_to = "Cantidad")
ggplot(comparacion_larga, aes(x = Metabolito, y = Cantidad, fill = Grupo)) + # Gráfico concentraciones
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Comparación de Concentracion de Metabolitos entre Cachexia y Control",
       x = "Metabolito",
       y = "Concentracion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Cachexia" = "red", "Control" = "blue"))

save(se, file = "human_cachexia.Rda")
