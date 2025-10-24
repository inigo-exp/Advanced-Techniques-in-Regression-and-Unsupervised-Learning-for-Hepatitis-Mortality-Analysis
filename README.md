# Advanced-Techniques-in-Regression-and-Unsupervised-Learning-for-Hepatitis-Mortality-Analysis

This project explores statistical and machine learning approaches to predict mortality in hepatitis patients. It combines regression and unsupervised learning methods to identify the most relevant clinical and biochemical factors influencing survival outcomes.

🎯 Objectives

Predict mortality risk based on patient clinical data.

Identify the most relevant medical features associated with liver dysfunction.

Group patients with similar clinical profiles using clustering techniques.

⚙️ Methods

Supervised Learning:

Logistic Regression

LASSO and Adaptive LASSO for feature selection and regularization

Bayesian Model Averaging (BMA) for model inference and uncertainty

Unsupervised Learning:

Principal Component Analysis (PCA) for dimensionality reduction

K-Means and Hierarchical Clustering for patient grouping

📁 Repository Structure
File	Description
Data.csv	Hepatitis dataset used for analysis
Regression_Techniques.R	Regression models and model selection procedures
Unsupervised_Learning.R	PCA and clustering implementation
Report.pdf	Full report with methodology, results, and interpretation
README.md	Project documentation
🧰 Tools

Language: R

Packages: glmnet, mice, BAS, mombf, ggplot2, psych, pROC

🚀 How to Run

Clone the repository

git clone https://github.com/inigo-exp/hepatitis-analysis.git
cd hepatitis-analysis


Open R or RStudio

Run the scripts in order:

Regression_Techniques.R → regression models and variable selection

Unsupervised_Learning.R → PCA and clustering

👥 Authors

Iñigo Expósito & Oriol Gelabert
📅 December 2024
