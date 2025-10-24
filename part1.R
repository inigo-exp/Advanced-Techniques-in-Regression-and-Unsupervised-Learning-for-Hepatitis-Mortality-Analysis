#PART1: LOGISTIC REGRESSION, NON LINEAR MODELS AND LASSO
set.seed(142)
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

db$class[db$class=="live"]<- "1"
db$class[db$class=="die"]<- "0"
db$class<- as.numeric(db$class)
plot(db$class,ylab = "")

#Check for outliers in continuous variables
boxplot(db$age)
boxplot(db$bilirubin)
boxplot(db$alk_phosphate)
boxplot(db$sgot)
boxplot(db$albumin)
boxplot(db$protime)


#LASSO REGRESSION

#install.packages("glmnet")
library(glmnet)

# Use colnames()
colnames(db)

#Load neccesary libraries
library(glmnet)

# Prepare data
db$class <- as.factor(db$class) 


#Function to compute AIC, BIC and EBIC
calculate_criteria <- function(y, y_hat, coef_matrix, X, gamma = 1) {
  n <- length(y)  
  p <- ncol(X)    
  rss <- sum((y - y_hat)^2)  
  df <- sum(coef_matrix != 0)  
  log_lik <- -n/2 * (log(2 * pi) + log(rss/n) + 1)  
  
  #information criterion
  aic <- -2 * log_lik + 2 * df
  bic <- -2 * log_lik + log(n) * df
  ebic <- bic + 2 * gamma * log(p) * df
  
  return(list(AIC = aic, BIC = bic, EBIC = ebic))
}


#ordinary logictic regression
modelo1 <- glm(class ~ age + sex + steroid + antivirals + fatigue + malaise + anorexia + liver_big + 
                liver_firm + spleen_palpable + spiders + ascites + varices + bilirubin + alk_phosphate + 
                sgot + albumin + protime + histology, 
              family = binomial(link = "logit"), 
              data = db)
summary(modelo1)

#take into account second order ineractions
model_formula <- as.formula(paste(
  "class ~ (", 
  paste(names(db)[-length(names(db))], collapse = " + "), 
  ")^2"
))
#ordinary logistic regression with interactions(more parameters than data)
modelo2 <- glm(model_formula, family = binomial(link = "logit"), data = db)
summary(modelo2)

#LASSO REGRESSION
set.seed(142)
X <- model.matrix(model_formula,data=db)[,-1]
#target variable
y <- db$class

#Fit Lasso ith cross validation
cv_lasso <- cv.glmnet(X, y, family = "binomial", alpha = 1)
#Check best value for lamda
cv_lasso$lambda.min

#Get coefficients for best lambda value
coef_lasso<-coef(cv_lasso, s = "lambda.min")[,1]


#Create a data frame for coefficients
coef_df <- as.data.frame(as.matrix(coef_lasso))
coef_df$Variable <- rownames(coef_df)  

#Filter non cero coeficients
coef_df_non_zero <- coef_df[coef_df$`V1` != 0, ]

#Compute AIC,BIC and EBIC for this model
y<-as.numeric(db$class)-1
fitted_probs <- predict(cv_lasso, newx = X, s = cv_lasso$lambda.min, type = "response")
fitted_probs <- as.numeric(as.character(fitted_probs))
coef_matrix <- as.matrix(coef(cv_lasso, s = cv_lasso$lambda.min))
matrix <- as.matrix(db[, -1])  # Asumiendo que la primera columna es la variable respuesta
criteria<-calculate_criteria(y,fitted_probs,coef_matrix,matrix)
print(criteria)

#ADAPTIVE LASSO
set.seed(142)
coef_initial <- as.vector(coef(cv_lasso, s = "lambda.min"))[-1] 
weights <- 1 / abs(coef_initial)  
weights[is.infinite(weights)] <- 1e10  #How to handle null coefficients

#Fix adaptive Lasso
adaptive_lasso <- glmnet(X, y, alpha = 1, family = "binomial", penalty.factor = weights)

cv_adaptive_lasso <- cv.glmnet(X, y, alpha = 1, family = "binomial", penalty.factor = weights)

#Final coefficients with best lambda in cross validation
best_lambda <- cv_adaptive_lasso$lambda.min
coef_adaptative<-coef(cv_adaptive_lasso, s = best_lambda)
coef_adaptative_df<- as.data.frame(as.matrix(coef_adaptative))
coef_adaptative_df$Variable <- rownames(coef_adaptative_df)
coef_adaptative_df_non_zero <- coef_adaptative_df[coef_adaptative_df$`s1` != 0, ]

#Compute AIC,BIC and EBIC for this model
y<-as.numeric(db$class)-1
fitted_probs_adapt <- predict(cv_adaptive_lasso, newx = X, s = cv_adaptive_lasso$lambda.min, type = "response")
fitted_probs_adapt <- as.numeric(as.character(fitted_probs_adapt))
coef_matrix_adapt <- as.matrix(coef(cv_adaptive_lasso, s = cv_adaptive_lasso$lambda.min))
matrix_adapt <- as.matrix(db[, -1])  # Asumiendo que la primera columna es la variable respuesta
criteria_adapt<-calculate_criteria(y,fitted_probs_adapt,coef_matrix_adapt,matrix_adapt)
print(criteria_adapt)


################################################################################
################           ACCURACY FOR EACH MODEL            ##################
################################################################################

#AUC and CLASSIFICATION TABLES
#install.packages("pROC")
library(pROC)
roc_el<-roc(db$class,fitted_probs, quiet=TRUE)
roc_el$auc
plot(roc_el,print.thres=TRUE,print.thres.col="dark green")

db$probs<-ifelse (fitted_probs < 0.757, 0, 1)
.Tablew1 <- xtabs(~db$class+db$probs)
.Tablew1
prop.table(.Tablew1,1) 

roc_el2<-roc(db$class,fitted_probs_adapt, quiet=TRUE)
roc_el2$auc
plot(roc_el2,print.thres=TRUE,print.thres.col="dark green")
db$probs_adapt<-ifelse (fitted_probs_adapt < 0.665, 0, 1)
.Tablew1 <- xtabs(~db$class+db$probs_adapt)
.Tablew1
prop.table(.Tablew1,1) 

#remove previous information about probabilities
colnames(db)
db<-db[,-c(21,22,23)]


#BAYESIAN MODEL AVERAGING

#install.packages('BMA')
library(BMA)
bma_model<-bic.glm(x = db[, -c(20)],y=db$class,glm.family=binomial(link="logit"))
summary(bma_model)
formula_full <- as.formula(
  class ~ (age + sex + steroid + antivirals + fatigue + malaise + anorexia + 
             liver_big + liver_firm + spleen_palpable + spiders + ascites + 
             varices + bilirubin + alk_phosphate + sgot + albumin + protime + 
             histology)^2
)
interaction_data <- model.matrix(formula_full, data = db)

#Fit BMA model
bma_model1 <- bic.glm(interaction_data, y = db$class, glm.family = binomial(link = "logit"), data = db)

summary(bma_model1)
#with lasso
formula_full <- as.formula(
  albumin ~ age + spiders + malaise + histology + liver_big + ascites + spleen_palpable +
    bilirubin + varices + age*spiders + malaise*histology + liver_big*ascites + 
    liver_big*histology + spleen_palpable*bilirubin + ascites*varices + ascites*bilirubin
)
adaptive_lasso_df <- model.matrix(formula_full, data = db)

#Set data frame format
adaptive_lasso_df <- as.data.frame(adaptive_lasso_df)
adaptive_lasso_df<-adaptive_lasso_df[,-c(1)]

#Use the correct formula
bma_model1 <- bic.glm(x = adaptive_lasso_df, y = db$class, glm.family = binomial(link = "logit"))

#Summary of the model
summary(bma_model1)

#####
install.packages('mombf')
library(mombf)
calculateEBIC <- function(logLik, n, p, gamma = 1) {
  bic <- -2 * logLik + log(n) * p
  ebic <- bic + 2 * gamma * p * log(n)
  return(ebic)
}

bestModels <- function(X, y, priorCoef, gamma = 1) {
  n <- nrow(X)
  p <- ncol(X)
  
  # Ajustar modelo lineal
  fit <- lm(y ~ X)  # Incluye intercepto automáticamente
  
  # Log-verosimilitud usando residuos
  residuals <- y - fitted(fit)  # Predicción del modelo
  logLik <- -sum(residuals^2) / 2
  
  # Calcular EBIC
  ebic <- calculateEBIC(logLik, n, p, gamma)
  
  # Retornar modelo con EBIC calculado
  return(list(bestModel = coef(fit), ebic = ebic))
}


# Ejecutar la función
result <- bestModels(X, y, priorCoef = imomprior(tau = 0.1), gamma = 1)
print(result)

#Instalamos el paquete mombf para model selection
#dependencia necessaria

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("sparseMatrixStats")
install.packages("sparseMatrixStats")
install.packages("matrixStats")
install.packages("SparseM")

install.packages("mombf")

library(sparseMatrixStats)
library(matrixStats)
library(mombf)

#usamos el criterio BIC
formula <- as.formula(paste("y ~ (", paste((colnames(X), collapse=" + "), ")^2")))

fitbic= bestBIC(formula,data=X)

fitbic

summary(fitbic)
# Remove intercept column if it exists
X <- X[, !apply(X, 2, function(col) all(col == 1))]

#usamos el criterio EBIC
fitebic= bestEBIC(y~X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5] + X[, 6] + X[, 7] + 
               X[, 8] + X[, 9] + X[, 10] + X[, 11] + X[, 12] + X[, 13] + 
               X[, 14] + X[, 15] + X[, 16] + X[, 17] + X[, 18] + X[, 19]
)
# Fórmula con todas las interacciones de segundo orden
fitbic <- bestBIC(y ~ (X[, 1] + X[, 2] + X[, 3] + X[, 4] + X[, 5] + 
                         X[, 6] + X[, 7] + X[, 8] + X[, 9] + X[, 10] + 
                         X[, 11] + X[, 12] + X[, 13] + X[, 14] + X[, 15] + 
                         X[, 16] + X[, 17] + X[, 18] + X[, 19])^2)

fitbic
# Supongamos que X es un data.frame, lo usamos directamente
fitbic <- bestBIC(y ~ .^2 - 1, data = data.frame(X))

summary(fitbic)


#BAM????
#First model
priorCoef <- momprior(tau=0.368)  # Default MOM prior on parameters
priorDelta <- modelbbprior(1,1)   # Beta-Binomial prior for model space
fit1 <- modelSelection(y~db[, 1] + db[, 2] + db[, 3] + db[, 4] + db[, 5] + db[, 6] + db[, 7] + 
                         db[, 8] + db[, 9] + db[, 10] + db[, 11] + db[, 12] + db[, 13] + 
                         db[, 14] + db[, 15] + db[, 16] + db[, 17] + db[, 18] + db[, 19])
postProb(fit1)
coef(fit1)
#Second model
normalprior <- function(beta, mean = 0, sd = 100) {
  # Calculate the log of the normal prior density for the coefficient
  log_prior <- dnorm(beta, mean = mean, sd = sd, log = TRUE)
  return(log_prior)
}
# Fit a logistic regression model
model_logistic <- glm(y ~db[, 1] + db[, 2] + db[, 3] + db[, 4] + db[, 5] + db[, 6] + db[, 7] + 
                        db[, 8] + db[, 9] + db[, 10] + db[, 11] + db[, 12] + db[, 13] + 
                        db[, 14] + db[, 15] + db[, 16] + db[, 17] + db[, 18] + db[, 19] , family = binomial)

# Extract regression coefficients
priorCoef <- normalprior(coefficients(model_logistic),mean=0, sd=100)  # Weakly informative prior
priorDelta <- modelbbprior(a=1, b=6)  

fit2 <- modelSelection(y~db[, 1] + db[, 2] + db[, 3] + db[, 4] + db[, 5] + db[, 6] + db[, 7] + 
                         db[, 8] + db[, 9] + db[, 10] + db[, 11] + db[, 12] + db[, 13] + 
                         db[, 14] + db[, 15] + db[, 16] + db[, 17] + db[, 18] + db[, 19])
postProb(fit2)
coef(fit2)


