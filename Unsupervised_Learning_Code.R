#PART2: PCA and Clustering

#Read database
db<-read.csv(file="hepatitis.csv",header = T, sep = ",", na.strings = "NA", dec = "." )
#Columns with no information will be given NA
db[db == ""] <- NA

#Edit class variable
db$class[db$class=="live"]<- "1"
db$class[db$class=="die"]<- "0"
db$class<- as.numeric(db$class)

#Identify columns of class character
char_cols <- sapply(db, is.character)

#Convert them to factor
db[, char_cols] <- lapply(db[, char_cols], as.factor)
summary(db)

#Handle NAN
#Install mice package
#install.packages("mice")
library(mice)

methods <- make.method(db)
#Impute missing values with NAN choose m=3 since our database is not very big
imputed_data <- mice(db, m = 3, method = methods, maxit = 50, seed = 123)
summary(imputed_data)

#Know which methods have been used for each imputation
imputed_data$method
db <- complete(imputed_data)
summary(db)

#load Package
#install.packages('dplyr')
library(dplyr)

#PCA
#columns to transform
columns_to_transform <- c("steroid", "antivirals", "fatigue", "malaise", 
                          "anorexia", "ascites", "spiders", 
                          "liver_firm", "liver_big", "varices", "histology","spleen_palpable")

#Get numerical features in order to get dummies
db <- db %>%
  mutate(across(all_of(columns_to_transform), 
                ~ as.factor(recode_factor(.x,"False" = "0", "True" = "1"))))
summary(db)
#Not take gender into account
db_numeric <- db[,c(1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)]
db_numeric[] <- lapply(db_numeric, function(col) {
  if (is.factor(col)) {
    return(as.numeric(as.character(col)))
  }
  return(col)
})
#Turn into numeric to compute covariance matrix
library(ggplot2)
library(reshape2)

# Calculate the correlation matrix
a <- cor(db_numeric)

# Convert the correlation matrix into long format
a_long <- melt(a)

# Create the heatmap
ggplot(data = a_long, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
                       limit = c(-1, 1), space = "Lab") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Covariance Matrix")


# Filter covariances
filtered_cov <- ifelse(abs(a) > 0.3, a, NA)
print(filtered_cov)

library(ggplot2)
library(reshape2)

filtered_cov_long <- melt(filtered_cov, na.rm = TRUE)

#Create heatmap
ggplot(data = filtered_cov_long, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
                       limit = c(-1, 1), space = "Lab") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Select most correlated variables(clinical variables)
mat<- db_numeric[,c(9,10,11,12,13,14,15,16,17,19)]
y<-mat[,c(10)]
mat1<-mat[,c(1,2,3,4,5,6,7,8,9)]
#Scale variables
mat1<-scale(mat1)
mat1 <- as.data.frame(mat1)
sapply(mat, class)

##Examine normality of numerical variables
library(nortest)

shapiro.test(mat1$bilirubin) 
shapiro.test(mat1$alk_phosphate) 
shapiro.test(mat1$sgot) 
shapiro.test(mat1$albumin) 
shapiro.test(mat1$protime)
shapiro.test(mat1$spleen_palpable)
shapiro.test(mat1$spiders)
shapiro.test(mat1$ascites)
shapiro.test(mat1$varices)
#none of them is normal. so we can not do probabilistic pca

#Compute correlation matrix
cor1<- cor(mat1)

#get eigenvalues of correlation matrix
eigen(cor1)
eig <- eigen(cor1)
eigenvalues <- eig$values

#Plot to see the elbow of the curve
ggplot(data = data.frame(x = 1:length(eigenvalues), y = eigenvalues), aes(x = x, y = y)) +
  geom_point(color = "blue", size = 3) +  
  geom_line(color = "blue", linewidth = 1) +   
  labs(
    title = "Scree Plot",               
    x = "Number of Components",          
    y = "Eigenvalue"                     
  ) +
  theme_minimal()  
#The elbow takes place when we have only one eigenvalue. As said before, it is senseless so we take three.
library(psych)
#PCA with the number of components specified to check the variance explained
pca_result <- prcomp(mat1, scale. = TRUE)
#Plot the variance explained by each component
explained_variance <- summary(pca_result)$importance[2, ]  
explained_variance

# Create a data frame to explain it
variance_df <- data.frame(
  Component = 1:length(explained_variance),
  Variance_Explained = explained_variance
)
#Plot the results
ggplot(variance_df, aes(x = as.factor(Component), y = Variance_Explained)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_line(aes(group = 1, y = Variance_Explained), color = "red", linewidth = 1) +
  geom_point(color = "red", size = 2) +
  theme_minimal() +
  labs(title = "Proportion of Explained Variance",
       x = "Principal Component",
       y = "Explained Variance") +
  scale_x_discrete(labels = as.character(variance_df$Component))


#loadings
loadings<-pca_result$rotation
loadings[,1:3]
#scores
pca_rotated_scores <- pca_result$x[, 1:3]

#Save loadings in a data frame
loadings_df <- as.data.frame(loadings[,1:3])
loadings_df$Variable <- rownames(loadings_df)

#plot the results
library(tidyr)
loadings_long <- gather(loadings_df, key = "Component", value = "Loading", -Variable)
#Barplot per each component
ggplot(loadings_long, aes(x = Variable, y = Loading, fill = Component)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  theme_minimal() +
  labs(title = "PCA Loadings",
       x = "Variables",
       y = "Loadings") +
  scale_fill_brewer(palette = "Set2")

#Save each component for cluster analysis
mat1$PC1 <- pca_rotated_scores[,1]
mat1$PC2 <- pca_rotated_scores[,2]
mat1$PC3 <- pca_rotated_scores[,3]
mat1$class<-y

#plot the results
library(car)
# Set up the plotting area for three plots in one row
par(mfrow = c(1, 3), mar = c(5, 5, 3, 1))  

# Plot 1: PC2 vs PC1
plot(PC2~ PC1, data = subset(mat1, class == "1"), 
     col = "blue", pch = 16, cex = 1.5,  
     xlab = "PC1", ylab = "PC2", 
     main = "PC2 vs PC1", 
     xlim = c(min(mat1$PC1) - 1, max(mat1$PC1) + 1),  
     ylim = c(min(mat1$PC2) - 1, max(mat1$PC2) + 1),  
     cex.lab = 1.2, cex.main = 1.5, cex.axis = 1.1)  

# Add points for the second class (Death)
points(PC2 ~ PC1, data = subset(mat1, class == "0"), 
       col = "red", pch = 17, cex = 1.5)  # Red color, triangle shape, and larger points

# Add reference lines at 0
abline(h = 0, v = 0, col = "gray", lty = 2, lwd = 1)

# Customize the legend
legend("topright", legend = c("Alive", "Death"), 
       pch = c(16, 17), col = c("blue", "red"), 
       cex = 1.2, bg = "white", box.lwd = 1.5)


# Plot 2: PC3 vs PC1
plot(PC3 ~ PC1, data = subset(mat1, class == "1"), 
     col = "blue", pch = 16, cex = 1.5,  
     xlab = "PC1", ylab = "PC3", 
     main = "PC3 vs PC1", 
     xlim = c(min(mat1$PC1) - 1, max(mat1$PC1) + 1),  
     ylim = c(min(mat1$PC3) - 1, max(mat1$PC3) + 1),  
     cex.lab = 1.2, cex.main = 1.5, cex.axis = 1.1)  

# Add points for the second class (Death)
points(PC3 ~ PC1, data = subset(mat1, class == "0"), 
       col = "red", pch = 17, cex = 1.5)  

# Add reference lines at 0
abline(h = 0, v = 0, col = "gray", lty = 2, lwd = 1)

# Customize the legend
legend("topright", legend = c("Alive", "Death"), 
       pch = c(16, 17), col = c("blue", "red"), 
       cex = 1.2, bg = "white", box.lwd = 1.5)


# Plot 3: PC3 vs PC2
plot(PC3 ~ PC2, data = subset(mat1, class == "1"), 
     col = "blue", pch = 16, cex = 1.5,  #
     xlab = "PC2", ylab = "PC3", 
     main = "PC3 vs PC2", 
     xlim = c(min(mat1$PC2) - 1, max(mat1$PC2) + 1),  
     ylim = c(min(mat1$PC3) - 1, max(mat1$PC3) + 1),  
     cex.lab = 1.2, cex.main = 1.5, cex.axis = 1.1)  

# Add points for the second class (Death)
points(PC3 ~ PC2, data = subset(mat1, class == "0"), 
       col = "red", pch = 17, cex = 1.5)  

# Add reference lines at 0
abline(h = 0, v = 0, col = "gray", lty = 2, lwd = 1)

# Customize the legend
legend("topright", legend = c("Alive", "Death"), 
       pch = c(16, 17), col = c("blue", "red"), 
       cex = 1.2, bg = "white", box.lwd = 1.5)

# Reset the plot layout to default
par(mfrow = c(1, 1))


#CLUSTERING
mat2<-mat1
datuak <- mat2[,-c(1:9)]  # Remove unnecessary columns from mat2
datuak2 <- mat2[,-c(1:9)] #Create a copy for second analysis
# Get the distance matrix with Euclidean distance (squared)
distmatrix <- dist(datuak, method = "euclidean")
distmatrix2 <- distmatrix^2
distmatrix2  # Display the squared distance matrix
#set seed for reproducibility
set.seed(123) 

# Use Ward method for hierarchical clustering based on the factors F1 and F2
HClust.1 <- hclust(dist(model.matrix(~-1 + PC1 + PC2+ PC3, datuak))^2, method = "ward.D")

# Plot dendrogram
plot(HClust.1, main = "Cluster Dendrogram", xlab = "Individual number in the dataset", ylab = "", sub = "")

# Set 4 clusters from the dendrogram
# Improved hierarchical clustering plot

# Set plot size and margins
par(mar = c(5, 4, 4, 2) + 0.1) 

# Create the basic dendrogram plot
plot(HClust.1, 
     main = "Hierarchical Clustering: 4 Clusters",  # Title
     xlab ="",    # X-axis label
     ylab = "Height",                              # Y-axis label
     sub = "",                                     # Subtitle (leave empty)
     col.main = "darkblue",                        # Title color
     col.lab = "darkgreen",                        # Axis labels color
     cex.main = 1.5,                               # Title size
     cex.lab = 1.2,                                # Axis labels size
     cex.axis = 1.1,                               # Axis tick labels size
     font.lab = 2,                                 # Bold axis labels
     font.main = 2,                                # Bold title
     lwd = 2,                                      # Line width for dendrogram
     hang = -1,                                    # Make the labels aligned at the bottom
     col.axis = "black")                           # Color for axis tick labels

# Add colored rectangle around the 4 clusters
rect.hclust(HClust.1, k = 4, border = "blue")

# Reset plotting parameters (optional)
par(mfrow = c(1, 1))


# Print the number of individuals in each cluster
summary(as.factor(cutree(HClust.1, k =4 ))) 

# Calculate Cluster Centroids
by(model.matrix(~-1 + PC1 + PC2+ PC3, datuak), as.factor(cutree(HClust.1, k = 4)), colMeans) 


# Graphical representation 
# Calculate Cluster Centroids
centroids <- by(model.matrix(~-1 + PC1 + PC2 + PC3, datuak), 
                as.factor(cutree(HClust.1, k = 4)), 
                colMeans)

# Ensure centroids is a list of numeric vectors
centroids <- lapply(centroids, unlist)  # Convert each element to a vector if needed

# Combine them into a data frame
centroids_df <- do.call(rbind, centroids)

# Convert the result to a data frame if necessary
centroids_df <- as.data.frame(centroids_df)

# Assign column names
colnames(centroids_df) <- c("PC1", "PC2", "PC3")

# Add cluster identifiers (assuming there are 4 clusters)
centroids_df$cluster <- factor(1:4)  # Add cluster identifiers

# Plotly 3D scatter plot
library(plotly)

# Get the cluster assignments
clusters <- cutree(HClust.1, k = 4)

# Create the 3D scatter plot for the data points
fig <- plot_ly(
  x = mat1$PC1,  # First principal component (PC1)
  y = mat1$PC2,  # Second principal component (PC2)
  z = mat1$PC3,  # Third principal component (PC3)
  color = factor(clusters),  # Cluster colors
  colors = c("blue", "red1", "green1", "purple1"),  # Custom colors for clusters
  type = "scatter3d",  # 3D scatter plot
  mode = "markers",  # Display points as markers
  marker = list(size = 5, opacity = 0.7)  # Point size and opacity for data points
)

# Add centroids to the plot as a separate trace with larger and distinct markers
fig <- fig %>% add_trace(
  x = centroids_df$PC1,  # Centroid positions for PC1
  y = centroids_df$PC2,  # Centroid positions for PC2
  z = centroids_df$PC3,  # Centroid positions for PC3
  color = factor(centroids_df$cluster),  # Color the centroids by cluster
  colors = c("blue4", "red4", "green4", "purple4"),  # Matching centroid colors
  type = "scatter3d",
  mode = "markers+text",  # Display markers and text labels
  text = paste("Cluster", centroids_df$cluster),  # Cluster label
  marker = list(
    size = 10,  # Larger marker size for centroids
    symbol = "x",  # Use "x" symbol for centroids
    line = list(width = 3),  # Thicker border for centroid markers
    opacity = 2  # Full opacity for centroids
  ),
  showlegend = FALSE  # Hide the legend for centroids
)

# Customize layout
fig <- fig %>% layout(
  title = "3D Plot of Principal Components with Centroids",
  scene = list(
    xaxis = list(title = "PC1"),
    yaxis = list(title = "PC2"),
    zaxis = list(title = "PC3")
  )
)

# Show the plot
fig



# Add the cluster assignment to the 'datuak' dataframe
datuak$hclus.label <- cutree(HClust.1, k = 4)

# Check if the new cluster labels have been added
head(datuak)


###########
# K-means #
###########

#ELBOW METHOD

# Create the feature matrix
data_matrix <- model.matrix(~-1 + PC1 + PC2 +PC3, datuak2)

# Test k-means for different numbers of clusters
set.seed(123)  # Ensures reproducibility

total_withinss <- numeric()  

# Iterate
for (k in 2:10) {
  kmeans_result <- kmeans(data_matrix, centers = k, iter.max = 10)  # Ejecuta k-means para cada k
  total_withinss[k] <- kmeans_result$tot.withinss  # Guarda el within-cluster sum of squares
}
# Create a data frame to store k values and the corresponding total within-cluster sum of squares
elbow_data <- data.frame(
  k = 2:10,
  total_withinss = total_withinss[2:10]
)

# Create the elbow plot with ggplot2
ggplot(elbow_data, aes(x = k, y = total_withinss)) +
  geom_point(color = "darkblue", size = 3) +  # Points for each k
  geom_line(color = "blue", linewidth = 1) +  # Line connecting the points
  geom_smooth(method = "loess", color = "red", linetype = "dashed", size = 1) +  # Smoothed line
  labs(
    title = "Elbow Method for Determining Optimal k",
    x = "Number of Clusters",
    y = "Total Within-Cluster Sum of Squares"
  ) +
  theme_minimal(base_size = 14) +  
  theme(
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),  # Center title, larger font
    axis.title = element_text(size = 14),  # Axis title font size
    axis.text = element_text(size = 12),  # Axis labels font size
    plot.margin = margin(20, 20, 20, 20)  # Increase plot margins for readability
  ) +
  scale_x_continuous(breaks = 2:10)  # Ensure x-axis ticks correspond to the


#GAP STATISTIC
# Load necessary library
if (!require(cluster)) install.packages("cluster")
library(cluster)

data_matrix <- model.matrix(~-1 + PC1 + PC2+ PC3, datuak2)

# Set parameters for the gap statistic
set.seed(123)  # For reproducibility
max_clusters <- 10  # Maximum number of clusters to test

# Compute the gap statistic
gap_stat <- clusGap(data_matrix, 
                    FUN = kmeans, 
                    K.max = max_clusters, 
                    B = 50,  # Number of bootstrap samples for reference data
                    iter.max = 10)

# Display the gap statistic results
print(gap_stat)

# Improved gap statistic plot
plot(gap_stat, 
     main = "Gap Statistic for Optimal Number of Clusters",   # Title with a larger font
     xlab = "Number of Clusters",                             # X-axis label
     ylab = "Gap Statistic",                                  # Y-axis label
     col = "blue",                                           # Line color
     lwd = 2,                                                # Line width for better visibility
     pch = 16,                                               # Point type (filled circles)
     cex = 1.2,                                              # Point size
     col.main = "darkblue",                                   # Title color
     col.lab = "darkgreen",                                   # Axis label color
     cex.main = 1.5,                                          # Title size
     cex.lab = 1.2,                                           # Axis labels size
     cex.axis = 1.1,                                          # Axis tick labels size
     font.lab = 2,                                            # Bold axis labels
     font.main = 2,                                           # Bold title
     xlim = c(1, length(gap_stat$Tab[,1]))                    # Ensure the x-axis includes all clusters
)

# Add gridlines for better readability
grid(col = "gray", lty = "dotted")

#Add horizontal line at y = 0 to better highlight the gaps
abline(h = 0, col = "red", lty = 2)


#The maximum gap is observed at k=4 with a value of 0.5739

#K-MEANS IMPLEMENTATION
#Do 4 clusters
set.seed(142)
kmeans_result <- kmeans(model.matrix(~-1 +PC1 + PC2+ PC3, datuak2), centers = 4, iter.max = 1000000)

#Cluster sizes
kmeans_result$size 

# Calculate Cluster Centroids from kmeans_result$centers
centroids2_df <- as.data.frame(kmeans_result$centers)

# Assign column names based on your principal components
colnames(centroids2_df) <- c("PC1", "PC2", "PC3")

# Add cluster identifiers 
centroids2_df$cluster <- factor(1:nrow(centroids2_df))  

# View the centroids data frame
print(centroids2_df)


#inerties
kmeans_result$withinss # Within Cluster Sum of Squares
kmeans_result$tot.withinss # Total Within Sum of Squares
kmeans_result$betweenss # Between Cluster Sum of Squares

# Add cluster assignments to the original dataset
datuak2$cluster <- kmeans_result$cluster

assignCluster <- function(kmeans_result, data) {
  data$cluster <- kmeans_result$cluster
  return(data)
}
assignCluster(kmeans_result, datuak2)
# Use the function
datuak2 <- assignCluster(kmeans_result, datuak2)

# 3D Plot with Plotly 
library(plotly)

#Plot with colours each individual
fig <- plot_ly(
  data = datuak2, 
  x = ~PC1, 
  y = ~PC2, 
  z = ~PC3, 
  type = "scatter3d", 
  mode = "markers", 
  marker = list(size = 5),
  color = ~factor(cluster),  # Mapear clusters como factores
  colors = c("blue", "red", "green", "purple"),  # Definir colores para 4 clusters
  showlegend = TRUE  # Asegurar que la leyenda esté visible
)

#Add centroids
fig <- fig %>% add_trace(
  data = centroids2_df,  
  x = ~PC1, 
  y = ~PC2, 
  z = ~PC3, 
  type = "scatter3d", 
  mode = "markers+text", 
  text = paste("Centroid", centroids2_df$cluster),  
  marker = list(
    size = 10, 
    symbol = "x", 
    line = list(width = 3),
    color = ~factor(cluster), 
    colors = c("blue", "red", "green", "purple")  
  ),
  showlegend = FALSE  
)

#Legend
fig <- fig %>% layout(
  title = "3D K-Means Clustering with Centroids (4 Clusters)",
  scene = list(
    xaxis = list(title = "PC1"),
    yaxis = list(title = "PC2"),
    zaxis = list(title = "PC3")
  )
)

# Mostrar el gráfico
fig



#FINAL INTERPRETATION#
#Add the corresponding variable
datuak$class<-y
datuak2$class<-y
library(ggplot2)

#For K-MEANS
datuak2$survived=datuak2$class==1
# Calculate the proportion of survivors within each cluster
survival_proportion <- tapply(datuak2$survived, datuak2$cluster, function(x) mean(x, na.rm = TRUE))

# Plot the survival proportions using a bar plot
barplot(survival_proportion,
        names.arg = paste("Cluster", 1:length(survival_proportion)),
        col = "blue",
        main = "Proportion of Survivors in Each Cluster",
        xlab = "Cluster",
        ylab = "Proportion Survived",
        border = "black")

# Add text labels to the bars showing the proportion
text(x = seq_along(survival_proportion), 
     y = survival_proportion + 0.05,  
     labels = round(survival_proportion, 2), 
     col = "black")

#Using ggplot2 for better visualization
survival_data <- data.frame(
  Cluster = factor(1:length(survival_proportion)),  # Cluster labels
  Proportion = survival_proportion  # Survival proportions
)

# Plot using ggplot2
ggplot(survival_data, aes(x = Cluster, y = Proportion)) +
  geom_bar(stat = "identity", fill = "blue", color = "black") +
  geom_text(aes(label = round(Proportion, 2)), vjust = -0.5, color = "black") +
  labs(title = "Proportion of Survivors in Each Cluster",
       x = "Cluster",
       y = "Proportion Survived") +
  theme_minimal()


#For Ward's method
datuak$survived=datuak$class==1
# Calculate the proportion of survivors within each cluster
survival_proportion <- tapply(datuak$survived, datuak$hclus.label, function(x) mean(x, na.rm = TRUE))

# Plot the survival proportions using a bar plot
barplot(survival_proportion,
        names.arg = paste("Cluster", 1:length(survival_proportion)),
        col = "blue",
        main = "Proportion of Survivors in Each Cluster",
        xlab = "Cluster",
        ylab = "Proportion Survived",
        border = "black")

# Add text labels to the bars showing the proportion
text(x = seq_along(survival_proportion), 
     y = survival_proportion + 0.05,  # Position the text slightly above the bars
     labels = round(survival_proportion, 2), 
     col = "black")

# Using ggplot2 for better visualization
survival_data <- data.frame(
  Cluster = factor(1:length(survival_proportion)),  # Cluster labels
  Proportion = survival_proportion  # Survival proportions
)

# Plot using ggplot2
ggplot(survival_data, aes(x = Cluster, y = Proportion)) +
  geom_bar(stat = "identity", fill = "blue", color = "black") +
  geom_text(aes(label = round(Proportion, 2)), vjust = -0.5, color = "black") +
  labs(title = "Proportion of Survivors in Each Cluster",
       x = "Cluster",
       y = "Proportion Survived") +
  theme_minimal()
  












