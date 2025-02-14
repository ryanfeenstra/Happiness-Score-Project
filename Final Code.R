g1_no_na <- g1[, !colSums(is.na(g1))> 0]
h1_no_na <- h1[,!colSums(is.na(h1))>0] 
m_data <- merge(g1,h1,by.x = "CountryName", by.y = "Country name")
m_data_no_na <- m_data[, !colSums(is.na(m_data)) > 0] 
m_data_imputed <- m_data_no_na %>%
  mutate(across(everything(), ~ifelse(is.na(.), mean(., na.rm = TRUE), .))) 
final_data_m <- subset(m_data_imputed, select = -c(CountryName,X,CountryCode)) 
g1_imputed <- g1_no_na %>%
  mutate(across(everything(), ~ifelse(is.na(.), mean(., na.rm = TRUE), .))) 
h1_imputed <- h1_no_na %>%
  mutate(across(everything(), ~ifelse(is.na(.), mean(., na.rm = TRUE), .))) 
merged_data <- merge(g1_imputed,h1_imputed, by.x = "CountryName", by.y = "Country name")
final_data <- subset(merged_data, select = -c(CountryName,X,CountryCode)) 
sample_size <- floor(0.8*nrow(final_data_m)) 
set.seed(777) 
choice = sample(seq_len(nrow(final_data_m)),size=sample_size) 
train_data = final_data_m[choice, ]
test_data= final_data_m[-choice,] 
x1 <- model.matrix(~. - 1, data = train_data[, -which(names(train_data) ==
                                                        "Happiness_score")])
y1 <- train_data$Happiness_score #Assigns the dependent variable to y1
lasso_regression <- glmnet(x1,y1, alpha = 1, lambda = 10^seq(-2,5, length.out = 50),
                           standardize = TRUE) 
plot(lasso_regression, xvar = "lambda", label = TRUE, las = 1, main =
       "Figure 1: Coefficient Path for Lasso Regression", cex.main = 0.8) 
title(sub = "Figure 1.1: Coefficient Path for the Lasso Regression,
      showing the convergence of coefficients toward zero", cex.sub = 0.6) 
lasso_cv <- cv.glmnet(x1,y1, lambda = 10^seq(-2, 5, length.out = 50), nfolds = 10) 
lasso1se <- lasso_cv$lambda.1se 
lassomin <- lasso_cv$lambda.min 
lasso_min <- glmnet(x1,y1,alpha = 1, lambda = lassomin) 
lasso_1se <- glmnet(x1,y1,alpha = 1, lambda = lasso1se) 
x_t <- model.matrix(~. - 1, data = test_data[, -which(names(test_data) ==
                                                        "Happiness_score")]) 
y_t <- test_data$Happiness_score 
lasso1se_prediction <- predict(lasso_1se, x_t) 
lasso1se_rmse <- round(sqrt(mean((y_t - lasso1se_prediction)^2)),2) 
lassomin_prediction <- predict(lasso_min, x_t) 
lassomin_rmse <- round(sqrt(mean(y_t - lassomin_prediction)^2),2) 
rmse_table1 <- data.frame(Regression_Method = c("Lasso1SE",
                                                "Lasso Min"),
                          RMSE = c(lasso1se_rmse,lassomin_rmse)) 
kable(rmse_table1, caption =
        "Lasso1SE and LassoMin RMSE Values Table") 
pred_variables <- colnames(x_t) 
pred_data <- p1[,pred_variables] 
x_p <- model.matrix(~. - 1, data = pred_data) 
lasso_prediction <- predict(lasso_min, s = lassomin, newx = x_p) 
lasso_prediction_df <- data.frame(Countries = c("Costa Rica", "Croatia", "Peru"),
                                  Predictions = round(lasso_prediction,2)) 
colnames(lasso_prediction_df) <- c("Countries","Happiness Score") 
kable(lasso_prediction_df, caption = 
        length_checker <- function(x){ 
          if (length(unique(x)) == 2){ 
            return(TRUE) }else{ return(FALSE) 
            }
        }
      constant_variance_checker <- function(x){ 
        if(length(unique(x)) == 1){ 
          return(TRUE) }else{ return(FALSE) 
          }
      }
      zero_checker <- function(x){ #
        if(any(x==0)){ 
          return(TRUE) }else{ return(FALSE) 
          }
      })
      binary_data <- sapply(final_data, length_checker) 
      b1 <- final_data[, !binary_data] 
      variance_data <- sapply(b1, constant_variance_checker) 
      b1 <- b1[, !variance_data] 
      zero_data <- sapply(b1, zero_checker)  
      b1 <- b1[,!zero_data] 
      b1_y <- b1
      b1 <- b1[, names(b1) != "Happiness_score"] 
      b1_final <- b1[1:48] 
      pca <- princomp(b1_final, cor = TRUE, scores = T) 
      permtestPCA <- function (X, nTests = 100, alpha = 0.05, center.data = TRUE,
                               scale.data = TRUE){ 
        n <- nrow(X) 
        m <- ncol(X) 
        X <- scale(X, center = center.data, scale = scale.data) 
        if (scale.data) {a <- 1/(n - 1)} else {a <- 1} 
        res.X <- prcomp(X) 
        eigs.X<- res.X$sdev^2 
        eigs.Xperm <- matrix(0, m, nTests) 
        Xperm <- matrix(0, n, m) 
        Xperm[, 1] <- X[, 1]; 
        for (i in 1:nTests){ 
          for (j in 2:m) { 
            ind <- sort(runif(n), index.return = TRUE)$ix 
            Xperm[, j] <- X[ind, j] }
          res.Xperm  <- prcomp(Xperm) 
          eigs.Xperm[, i] <- res.Xperm$sdev^2}
        perc.alpha <- matrix(0, m, 2) 
        for (s in 1:m){ 
          perc.alpha[s,] <- quantile(eigs.Xperm[s,], c(alpha/2, 1 - alpha/2) ) 
        }
        plot(1:m, eigs.X, type = "b", col = "red",
             main = " Figure 2: Permutation test PCA", xlab = "Component",
             ylab = "Eigenvalue", sub = "Figure 2: Permutation Test to
       find optimal amount of components", cex.sub = 0.8) 
        lines(1:m, perc.alpha[, 1], type = "b", col="blue") 
        lines(1:m, perc.alpha[, 2], type = "b", col="green") 
        string1 <- paste("Confidence: ",formatC(alpha/2, digits=3, width=1, format="f")) 
        string2 <- paste("Confidence: ",formatC(1-alpha/2, digits=3, width=1, format="f")) 
        legend("topright", inset=.05, c("Observed", string1, string2), 
               lty = c(1, 1, 1), col = c("red", "green", "blue"), pch = c("o", "o", "o")) 
        return(perc.alpha) 
      }
      perm_test <- permtestPCA(b1_final) 
      kaiser_bootstrap <- rep(NA, 1000) 
      for (i in 1:1000){ 
        boot_s <- sample(1:nrow(b1_final), replace = TRUE) 
        boot_d <- b1_final[boot_s,] 
        boot_pca <- princomp(boot_d,cor = TRUE, scores = TRUE) 
        kaiser_bootstrap[i] <- sum(boot_pca$sdev^2 > 1) 
      }
      actual_value <- sum(pca$sdev^2 > 1) 
      lower_ci <- quantile(kaiser_bootstrap, c(0.025,0.975))  
      hist(kaiser_bootstrap, main = "Figure 3:
     Bootstrapping for Kaiser Rule on the 3 Principal Components",
           xlab = "Components Greater than 1", ylab = "Frequency",
           sub = "Figure 3: Bootstrapping the 3 principal components chosen to
     test if the Kaiser Rule is fulfilled", cex.sub = 0.6, cex.main = 0.7) 
      abline(v = lower_ci, col = "red", lwd = 2) 
      abline(v = actual_value, col = "green", lwd = 2) 
      b1_final$happiness_score <- b1_y$Happiness_score 
      sample_size_pcr <- floor(0.8*nrow(b1_final))
      set.seed(777) 
      choice_pcr = sample(seq_len(nrow(b1_final)),size=sample_size_pcr) 
      train_data_pcr = b1_final[choice_pcr, ] 
      test_data_pcr = b1_final[-choice_pcr,] 
      pcr_model <- pcr(data = train_data_pcr, happiness_score ~ ., validation = NULL,
                       scale = TRUE) 
      pcr_pred <- predict(pcr_model, test_data_pcr, ncomp = 3) 
      y_pcr <- test_data_pcr$happiness_score 
      pcr_rmse <- round(sqrt(mean((y_pcr - pcr_pred)^2)),2) 
      rmse_table2 <- data.frame(Regression_Method = c("Lasso1SE",
                                                      "Lasso Min", "PCR Regression"),
                                RMSE = c(lasso1se_rmse,lassomin_rmse,pcr_rmse)) 
      kable(rmse_table2, caption = "RMSE Values for Lasso1se, Lasso Min and PCR") 
      pred_variables_pcr <- colnames(b1_final[1:ncol(b1_final)-1])
      pred_data_pcr <- p1[,pred_variables_pcr] 
      pcr_prediction <- predict(pcr_model, pred_data_pcr, ncomp = 3) 
      pcr_prediction_df <- data.frame(Countries = c("Costa Rica", "Croatia", "Peru"),
                                      Predictions = round(pcr_prediction,2)) 
      colnames(pcr_prediction_df) <- c("Countries","Happiness Score")
      kable(pcr_prediction_df, caption = "PCR Regression Prediction for the Happiness
      Score of Costa Rica, Croatia, and Peru") 
      lasso1se_coefficients <- round(coef(lasso_1se)[coef(lasso_1se)[,1]!=0,],2)
      lassomin_coefficients <- round(coef(lasso_min)[coef(lasso_min)[,1]!=0,],2) 
      kable(lasso1se_coefficients, caption = 
              kable(lassomin_coefficients, caption = 
                      kable(round((pca$sdev[1:5]),2), caption = 
                              pca_sum <- summary(pca) #
                            kable((round(pca_sum$loadings[, 1:3][1:10,],2)), caption = "First 10 Component
      Loadings for the three PCA Components")