# Codi per a l'exploració de dades - Genòmica Computacional

#Carrego les llibreries SummarizedExperiment i ggplot2 per treballar amb l'objecte se i realitzar gràfics per visualitzar les dades:
library(SummarizedExperiment)
library(ggplot2) 
#Carrego l'objecte se (arxiu:summarized_experiment.Rda) i comprovo que té les dimensions correctes:
load("summarized_experiment.Rda")
dim(se)
#Comprovo si hi ha valors nuls en el conjunt de dades:
sum(is.na(assay(se)))
#Calculo la variància de cada fosfopèptid (fila) en el conjunt de dades:
row_variances <- apply(assay(se), 1, var) 
#Determino quantes files tenen una variància zero:
zero_var_rows <- sum(row_variances == 0)
print(zero_var_rows)
#Creo un nou objecte 'filtered_se' que elimina les files amb variància zero:
filtered_se <- se[row_variances > 0, ]
#Es crea un dataframe anomenat df amb dues columnes: una amb els valors d'intensitat de senyal dels fosfopèptids (convertits a vector) i una altra amb els noms de les mostres tumorals, 
#amb la funció rep que repeteix els noms de les mostres tantes vegades com el nombre de metabolits per tal que cada metabolit tingui el nom de mostra associat:
df <- data.frame(Intensity = as.vector(assay(filtered_se)),
                 Sample = rep(colnames(filtered_se), each = nrow(filtered_se)))
#Creo un boxplot per mostrar la distribució de les intensitats de senyal dels fosfopèptids analitzats per cada mostra:
ggplot(df, aes(x = Sample, y = Intensity)) +
  geom_boxplot(fill = "lightblue") +
  theme_minimal() +
  labs(title = "Distribució dels valors d'intensitat dels metabolits") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#Es realitza el PCA sobre la matriu transposada (els components a separar, en aquest cas les mostres, han d'estar en files) de dades d'intensitat dels fosfopèptids filtrades amb la funció prcomp() i es guarden els resultats de l'anàlisi a pca_res:
pca_res <- prcomp(t(assay(filtered_se)), scale. = TRUE)
#Calculo el percentatge de la variància explicada per cada component del PCA:
explained_variance <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100
#Creo un dataframe amb dos columnes: la columna PC amb el nom dels components principals i la columna Variance amb el percentatge de la variància explicada per cada component:
explained_variance_df <- data.frame(
  PC = paste0("PC", 1:length(explained_variance)),
  Variance = explained_variance)
print(explained_variance_df)
#Creo una vector de noms pels components principals de la variable explained_variance:
pc_labels <- paste0("PC", 1:length(explained_variance))
#creo un gràfic de barres que mostra el percentatge de la variància explicada per cada component principal:
barplot(explained_variance, names.arg = pc_labels, main = "Variància explicada per cada PC",
        xlab = "Components Principals", ylab = "Variància (%)", col = "steelblue", 
        ylim = c(0, max(explained_variance) * 1.1), las = 2) 
#Calculo la matriu de correlació entre els tres primers components principals (PC1, PC2 i PC3) i la imprimeixo per determinar si hi ha correlació entre aquests components del PCA:
correlation_matrix <- cor(pca_res$x[, 1:3]) 
print(correlation_matrix)
#Creo un dataframe anomenat pca_df amb dos columnes: una que conté els resultats del PCA pels components PC1 i PC2 i, la columna Sample, que conté el SampleID de cada mostra, extret de les metadades associades a l'objecte filtered_se: 
pca_df <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2], Sample = colData(filtered_se)$SampleID)
#creo un gràfic de dispersió que mostra les mostres en un espai bidimensional reduït per PCA, utilitzant els components PC1 i PC2. Cada punt del gràfic representa una mostra, i els colors indiquen les diferents mostres segons el SampleID. 
ggplot(pca_df, aes(x = PC1, y = PC2, color = Sample)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(title = "PCA: PC1 vs PC2")
#Carrego la llibreria Plotly per realitzar gràfics avançats:
library(plotly)
#creo un dataframe anomenat df_pca amb dues columnes: una amb els resultats del PCA dels components PC1,PC2 i PC3, i els IDs de les mostres extrets de les metadades a través de col_data$SampleID:
df_pca <- data.frame(pca_res$x[, 1:3], SampleID = col_data$SampleID)
#Creo un gràfic interactiu en 3D amb la funció plot_ly() que assigna els components als eixos X, Y i Z. Cada punt representa una mostra i es distingeixen per colors segons el seu SampleID:
plot_ly(df_pca, x = ~PC1, y = ~PC2, z = ~PC3, color = ~SampleID, colors = "Paired") %>%
  add_markers() %>%
  layout(title = "PCA: PC1 vs PC2 vs PC3")
# Calculo la distància euclidiana entre les mostres en el pla de PC1 i PC2:
dist_2d <- dist(df_pca[, c("PC1", "PC2")])
#Calculo la distància euclidiana entre les mostres en un espai tridimensional format per PC1, PC2 i PC3:
dist_3d <- dist(df_pca[, c("PC1", "PC2", "PC3")])
#Es mostra un resum estadístic de les distàncies calculades:
summary(dist_2d)
summary(dist_3d)
#converteixo les coordenades de les mostres a l'espai PCA emmagatzemades a pca_res$x en un dataframe anomenat pca_scores:
pca_scores <- as.data.frame(pca_res$x)
#Afegeixo una nova columna anomenada Sample al dataframe pca_scores que conté els noms de les mostres:
pca_scores$Sample <- rownames(pca_scores)
#Especifico un valor de seed per garantir la reproduibilitat de l'anàlisi:
set.seed(123)
#Especifico el nombre de clústers
num_clusters <- 3 
#Realitzo un anàlisi de clustering k-means utilitzant els components PC1, PC2, i PC3 i guardo els resultats a kmeans_result:
kmeans_result <- kmeans(pca_scores[, 1:3], centers = num_clusters)
#Creo el dataframe "clustered_data":
clustered_data <- as.data.frame(pca_scores)
#AfegeixO una nova columna anomenada Cluster a clustered_data que indica a quin cluster ha estat assignada cada mostra segons k-means:
clustered_data$Cluster <- as.factor(kmeans_result$cluster)
#Amb la funció plot_ly, creo un gràfic en 3D amb PC1, PC2, PC3 de les mostres. Els punts de les mostres tenen un color assignat segons el cluster otorgat pel k-means:
fig <- plot_ly(clustered_data, x = ~PC1, y = ~PC2, z = ~PC3, 
               color = ~Cluster, colors = c('#1f77b4', '#ff7f0e', '#2ca02c'), 
               type = 'scatter3d', mode = 'markers', marker = list(size = 5)) %>%
  layout(title = "Clustering K-means a 3D (PC1, PC2, PC3)",
         scene = list(xaxis = list(title = 'PC1'),
                      yaxis = list(title = 'PC2'),
                      zaxis = list(title = 'PC3')))
fig
#Assigno els noms de les files (fosfopèptids) a l'objecte se basant-me en les modificacions de seqüència de cada fosfopèptid:
rownames(se) <- rowData(se)$SequenceModifications
#Elimino les files que corresponen a fosfopèptids amb variància zero:
rownames(filtered_se) <- rownames(se)[row_variances > 0]
#Creo un nou dataframe (row_data_filtered) que conté només les dades associades als metabolits amb variància > 0:
row_data_filtered <- rowData(se)[rownames(filtered_se), ]
#Creo un dataframe (loadings_df) a partir de la matriu de rotació del PCA que conté els coeficients de cada component principal per cada fosfopèptid:
loadings_df <- as.data.frame(pca_res$rotation)
#Afegeixo una columna a loadings_df que conté els noms dels metabolits seleccionats com a identificador:
loadings_df$Metabolite <- rownames(filtered_se)
#Afegeixo les dades associades a cada metabolit al dataframe loadings_df:
loadings_df <- cbind(loadings_df, row_data_filtered)
#Ordeno el dataframe loadings_df segons el valor absolut del loading de PC1, de major a menor, per identificar quins fosfopèptids tenen més influència sobre aquest component:
top_pc1 <- loadings_df[order(abs(loadings_df$PC1), decreasing = TRUE), ]
#Extrec les modificacions de seqüència dels primers 10 elements ordenats segons la seva influència sobre PC1 i els imprimeixo:
top_metabolites_pc1 <- top_pc1$SequenceModifications
print(head(top_metabolites_pc1, 10))
#Segueixo el mateix procediment per identificar els fosfopèptids més influents sobre PC2 i PC3:
top_pc2 <- loadings_df[order(abs(loadings_df$PC2), decreasing = TRUE), ]
top_metabolites_pc2 <- top_pc2$SequenceModifications
print(head(top_metabolites_pc2, 10))
top_pc3 <- loadings_df[order(abs(loadings_df$PC3), decreasing = TRUE), ]
top_metabolites_pc3 <- top_pc3$SequenceModifications
print(head(top_metabolites_pc3, 10))
