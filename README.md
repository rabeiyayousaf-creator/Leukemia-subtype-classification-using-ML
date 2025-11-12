# Leukemia Gene Expression Analysis

This repository contains an R script for the comprehensive analysis of gene expression data related to leukemia. The project involves various machine learning techniques, including dimensionality reduction, clustering, and classification algorithms, to explore patterns within the dataset and build predictive models for leukemia subtypes.

## Project Overview

The main goals of this project are:
*   To load and preprocess a leukemia gene expression dataset.
*   To visualize the dataset characteristics and gene distributions.
*   To apply dimensionality reduction techniques (PCA, t-SNE, UMAP) for data exploration and feature engineering.
*   To perform various clustering analyses (K-Means, Hierarchical Clustering, GMM, DBSCAN) to identify natural groupings within the data.
*   To build and evaluate classification models (Decision Trees, Random Forest) for predicting leukemia subtypes.
*   To identify important genes contributing to the classification.

## Dataset

The analysis uses a `Leukemia_dataset.csv` file, which is expected to contain gene expression levels and a 'type' column indicating the leukemia subtype for each sample. The script assumes the first column is 'samples' and the 'type' column is the target variable.

## Script Structure

The R script (`leukemia_analysis.R`) is organized into the following main sections:

1.  **Install & Load Libraries**: Installs and loads all necessary R packages for the analysis.
2.  **Load Dataset**: Reads the `Leukemia_dataset.csv` file and performs initial data preparation, including converting the target variable to a factor and isolating numeric gene expression data.
3.  **Visualization**: Generates initial plots to understand the distribution of leukemia subtypes and pairwise relationships between the first few genes.
4.  **Principal Component Analysis (PCA)**: Performs PCA to reduce dimensionality and visualizes the proportion of variance explained and the first two principal components.
5.  **Train-Test Split**: Divides the dataset into training and testing sets for model development and evaluation.
6.  **t-SNE**: Applies t-distributed Stochastic Neighbor Embedding for visualizing high-dimensional data in 2D.
7.  **UMAP**: Implements Uniform Manifold Approximation and Projection for dimensionality reduction and visualization.
8.  **K-Means Clustering**: Performs K-Means clustering on the reduced data.
9.  **Hierarchical Clustering**: Conducts hierarchical clustering and visualizes the dendrogram and clusters.
10. **Gaussian Mixture Models (GMM)**: Applies GMM for probabilistic clustering.
11. **DBSCAN**: Utilizes DBSCAN for density-based spatial clustering of applications with noise.
12. **Heatmap of Top Genes**: Identifies top variable genes and visualizes their expression patterns using a heatmap, often grouped by GMM clusters.
13. **Classification Tree**: Builds and evaluates a decision tree classifier to predict leukemia subtypes using the top variable genes.
14. **Regression Tree**: Demonstrates a regression tree model to predict the expression of a single gene based on others.
15. **Random Forest with Cross-Validation**: Implements a Random Forest model for classification, including cross-validation for robust evaluation and variable importance plotting.
16. **Evaluate with Function**: Provides a reusable function to calculate and display model accuracy and confusion matrix.

## Getting Started

### Prerequisites

*   R statistical software installed on your system.
*   RStudio (recommended IDE).

### Installation

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/YourUsername/Leukemia-Gene-Expression-Analysis.git
    cd Leukemia-Gene-Expression-Analysis
    ```
2.  **Place your dataset:**
    Ensure `Leukemia_dataset.csv` is in the same directory as the R script, or update the `setwd()` path in the script.
3.  **Run the R script:**
    Open `leukemia_analysis.R` in RStudio and run the entire script. The script will install any missing packages automatically.

    ```R
    # ---- Install Libraries ----
    install.packages(c("Rtsne", "umap", "mclust", "dbscan", "pheatmap", "rpart", "rpart.plot", "ipred", "randomForest", "gbm", "ggplot2", "GGally", "caret"))

    # ... rest of your script
    ```

## Contributing

Feel free to fork this repository, open issues, or submit pull requests.

