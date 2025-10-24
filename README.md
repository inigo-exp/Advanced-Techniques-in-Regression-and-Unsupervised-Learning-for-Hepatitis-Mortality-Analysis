# Hepatitis Mortality Analysis  
**Advanced Techniques in Regression and Unsupervised Learning**

This project applies advanced statistical and machine learning techniques to predict mortality in hepatitis patients. It combines regression and unsupervised learning methods to identify the most relevant clinical and biochemical factors associated with patient survival.

---

## Objectives
- Predict mortality risk using patient clinical and biochemical data.  
- Identify the most influential medical features related to liver dysfunction.  
- Group patients with similar clinical profiles through clustering techniques.  

---

## Methods
**Supervised Learning:**  
- Logistic Regression  
- LASSO and Adaptive LASSO for feature selection and regularization  
- Bayesian Model Averaging for model inference and uncertainty  

**Unsupervised Learning:**  
- Principal Component Analysis for dimensionality reduction  
- K-Means and Hierarchical Clustering for patient grouping  

---

## Repository Structure
| File | Description |
|------|--------------|
| `Data.csv` | Hepatitis dataset used for analysis |
| `Regression_Techniques.R` | Implementation of regression models and model selection |
| `Unsupervised_Learning.R` | PCA and clustering analysis |
| `Report.pdf` | Full report with methodology, results, and interpretation |
| `README.md` | Project documentation |

---

## Tools
- **Language:** R  
- **Packages:** `glmnet`, `mice`, `BAS`, `mombf`, `ggplot2`, `psych`, `pROC`  

---

## How to Run
1. Clone the repository  
   ```bash
   git clone https://github.com/inigo-exp/hepatitis-analysis.git
   cd hepatitis-analysis

