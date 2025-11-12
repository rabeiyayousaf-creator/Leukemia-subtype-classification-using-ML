# ---- Install Libraries ----
install.packages(c("Rtsne", "umap", "mclust", "dbscan", "pheatmap", "rpart", "rpart.plot", "ipred", "randomForest", "gbm", "ggplot2", "GGally", "caret"))

# ---- Load Libraries ----
library(Rtsne)
library(umap)
library(mclust)
library(dbscan)
library(pheatmap)
library(rpart)
library(rpart.plot)
library(ipred)
library(randomForest)
library(gbm)
library(ggplot2)
library(GGally)
library(caret)

# ---- Step 1: Load dataset ----
setwd("C:/Users/Rabia/Desktop")  # Modify path if needed

# Load the dataset
dataset1 <- read.csv("Leukemia_dataset.csv")

# Ensure target variable is categorical
dataset1$type <- as.factor(dataset1$type)

# Keep only numeric gene expression columns
dataset_numeric <- dataset1[, !(names(dataset1) %in% c("samples", "type"))]

# ---- Visualization ----
ggplot(dataset1, aes(x = type)) + 
  geom_bar(fill = "skyblue") + 
  theme_minimal() + 
  labs(title = "Distribution of Leukemia Subtypes", 
       x = "Leukemia Subtypes", y = "Frequency")

# Pairwise plot of first 5 genes
ggpairs(dataset_numeric[, 1:5], title = "Pairwise Plot of First 5 Genes")

# ---- Step 2: PCA ----
pr.out <- prcomp(dataset_numeric, scale. = TRUE)
pve <- pr.out$sdev^2 / sum(pr.out$sdev^2)

# Plot PVE and cumulative PVE
plot(pve, type = "b", xlab = "Principal Component", 
     ylab = "Proportion of Variance Explained", main = "PVE per PC")

# Cumulative variance explained
plot(cumsum(pve), type = "b", xlab = "Principal Component", 
     ylab = "Cumulative PVE", main = "Cumulative Proportion of Variance Explained")

# Scatterplot of first two PCs
plot(pr.out$x[, 1:2],
     col = as.factor(dataset1$type),
     pch = 19, cex = 1.3,
     xlab = "PC1", ylab = "PC2",
     main = "PCA: First Two Components")
legend("bottomright", legend = unique(dataset1$type),
       col = 1:length(unique(dataset1$type)), pch = 19, cex = 0.7, bty = "o")

# Reduced PC data for clustering
pc_data <- pr.out$x[, 1:50]

# ---- Step 3: Train-Test Split ----
set.seed(42)
train_idx <- sample(1:nrow(dataset1), 0.7 * nrow(dataset1))
train_data <- dataset1[train_idx, ]
test_data <- dataset1[-train_idx, ]

# ---- Step 4: t-SNE ----
set.seed(42)
tsne_out <- Rtsne(pc_data[train_idx, ], dims = 2, perplexity = 10, verbose = TRUE)
plot(tsne_out$Y, col = as.factor(train_data$type), pch = 19,
     xlab = "t-SNE 1", ylab = "t-SNE 2", main = "t-SNE on First 50 PCs")
legend("topright", legend = unique(train_data$type),
       col = 1:length(unique(train_data$type)), pch = 19, cex = 0.6, bty = "o")

# ---- Step 5: UMAP ----
set.seed(42)
umap_out <- umap(pc_data[train_idx, ])
plot(umap_out$layout, col = as.factor(train_data$type), pch = 19,
     xlab = "UMAP 1", ylab = "UMAP 2", main = "UMAP on First 50 PCs")
legend("bottomright", legend = unique(train_data$type),
       col = 1:length(unique(train_data$type)), pch = 19, cex = 0.6, bty = "o")

# ---- Step 6: K-Means ----
set.seed(1)
km.out <- kmeans(pc_data[train_idx, ], centers = 2, nstart = 20)
plot(pc_data[train_idx, 1:2], col = km.out$cluster,
     main = "K-Means Clustering (K=2)", xlab = "PC1", ylab = "PC2", pch = 20, cex = 2)

# ---- Step 7: Hierarchical Clustering ----
dist_matrix <- dist(pc_data[train_idx, ])
hc.complete <- hclust(dist_matrix, method = "complete")
plot(hc.complete, main = "Complete Linkage Dendrogram", xlab = "", sub = "", cex = 0.5)

clusters_complete <- cutree(hc.complete, k = 2)
plot(pc_data[train_idx, 1:2], col = clusters_complete, pch = 20, cex = 2,
     main = "Hierarchical Clustering (Complete Linkage)", xlab = "PC1", ylab = "PC2")

# ---- Step 8: Gaussian Mixture Models (GMM) ----
gmm_model <- Mclust(pc_data[train_idx, ], G = 1:5)
summary(gmm_model)
plot(gmm_model, what = "classification")

# ---- Step 9: DBSCAN ----
dbscan_out <- dbscan(pc_data[train_idx, ], eps = 5, minPts = 3)
plot(pc_data[train_idx, 1:2], col = dbscan_out$cluster + 1, pch = 19,
     main = "DBSCAN Clustering on First 50 PCs", xlab = "PC1", ylab = "PC2")

# ---- Step 10: Heatmap of Top Genes ----
gene_variances <- apply(dataset_numeric, 2, var)
top_genes <- names(sort(gene_variances, decreasing = TRUE))[1:50]
annotation_df <- data.frame(Cluster = factor(gmm_model$classification))
rownames(annotation_df) <- rownames(dataset_numeric)

pheatmap(dataset_numeric[train_idx, top_genes],
         annotation_row = annotation_df, scale = "row",
         clustering_method = "complete", show_rownames = FALSE,
         main = "Top 50 Variable Genes by GMM Cluster")

top_gene_df <- data.frame(type = train_data$type, dataset_numeric[train_idx, top_genes])

# ---- Step 11: Classification Tree ----
class_tree <- rpart(type ~ ., data = top_gene_df, method = "class",
                    control = rpart.control(cp = 0.01, maxdepth = 5))
rpart.plot(class_tree, type = 2, extra = 106, fallen.leaves = TRUE,
           main = "Classification Tree (Top 50 Genes)")

pred_class <- predict(class_tree, test_data, type = "class")
conf_matrix <- table(Predicted = pred_class, Actual = test_data$type)
print(conf_matrix)
accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)
cat("Classification Accuracy:", round(accuracy * 100, 2), "%\n")

# ---- Step 12: Regression Tree (predict most variable gene) ----
target_gene <- top_genes[1]
predictor_genes <- setdiff(top_genes, target_gene)

reg_tree <- rpart(dataset_numeric[, target_gene] ~ .,
                  data = dataset_numeric[, predictor_genes, drop = FALSE],
                  method = "anova",
                  control = rpart.control(cp = 0.01, maxdepth = 5))
rpart.plot(reg_tree, type = 2, extra = 101, fallen.leaves = TRUE,
           main = paste("Regression Tree for", target_gene))

# ---- Step 13: Random Forest with Cross-Validation ----
set.seed(123)
rf_model <- randomForest(type ~ ., data = top_gene_df, ntree = 500, mtry = 7, importance = TRUE)

# Cross-validation using caret
train_control <- trainControl(method = "cv", number = 10)
rf_cv_model <- train(type ~ ., data = top_gene_df, method = "rf", trControl = train_control)
print(rf_cv_model)

# Prediction on test set
rf_pred <- predict(rf_model, test_data)
rf_conf <- table(Predicted = rf_pred, Actual = test_data$type)
print(rf_conf)
rf_acc <- sum(diag(rf_conf)) / sum(rf_conf)
cat("Random Forest Accuracy:", round(rf_acc * 100, 2), "%\n")

varImpPlot(rf_model, main = "Random Forest Variable Importance (Top 50 Genes)")

# ---- Step 14: Evaluate with Function ----
evaluate_model <- function(predictions, actual) {
  cm <- table(Predicted = predictions, Actual = actual)
  accuracy <- sum(diag(cm)) / sum(cm)
  cat("Accuracy:", round(accuracy * 100, 2), "%\n")
  print(cm)
}

evaluate_model(pred_class, test_data$type)
