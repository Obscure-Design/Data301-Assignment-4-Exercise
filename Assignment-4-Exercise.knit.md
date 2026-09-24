---
title: "Assignment-4-Exercise"
author: "Vishnun-300665889"
date: "2026-09-24"
output: pdf_document
---

Exercise 1

1.
Computing of Order Statistic Summary

(H&F, 1996) had an in-depth look into the nine sample quantile definition, and their 'positive' over other definitions and their coverage of dealing with distributions. In the end, (H&F) noted that the current variation caused confusion, and aimed for standardized definition across the package, opting for the unbiased variation of denominator (n-1), adopted by the statistical community. Wicklin's surmised that despite the many statistical definitions to calculate the quantile range, by using the same small dataset for their calculation, the resulting quantile-quantile plots generated differ very slightly. Depending on the platform used such as R, had a standardized form of Type 7 as the default form (confirmed using ?quantile)

Personal Opinion: Defintion of Q8, gives a median, unbiased estimator, while the definition defines the properties are more difficult to obtain, the resulting use case is more robust and presents a more accurate result regardless of the distribution, other definitions have some caveats in that they cannot be used outside the normal distribution or other types of distributions. Most other definitions that are distribution-free fall short in a biased estimator. If the goal is for true accuracy, we require that the order statistic definition do not fall short of the true population parameter used for estimation.

2.
Suspicion of IQR is biased, Derive and implement first-order bootstrap bias:

``` r
lambda_IQR <- function(x) {
  q <- quantile(x, probs = c(0.25, 0.75),
                type = 7) #Since type 7 is default
  log(3) / q[2] - q[1]
}
#First Order Bootstrap
lambda_hat <- function(x, B = 200) {
  n <- length(x)
  
  lam_hat <- lambda_IQR(x)
  
  boot_estimates <- numeric(B)
  for (b in 1:B) {
    x_star <- sample(x, size = n, replace = TRUE)
    boot_estimates[b] <- lambda_IQR(x_star)
  }
  bias_hat <- mean(boot_estimates) - lam_hat
  
  lam_T <- lam_hat - bias_hat
  
  return(lam_T)
  }
```

``` r
#Use Case:
set.seed(1)
x <- rexp(n = 50, rate = 2) #Sample with true lambda = 2

lambda_IQR(x) #Original (Biased) Estimate
```

```
##      75% 
## 1.597416
```

``` r
lambda_hat(x) #Bias Corrected Estimate
```

```
##     75% 
## 1.58207
```
3. Monte Carlo Experiment 

``` r
lambda_ML <- function(x) {
  1 / mean(x)
}

lambda_tilde <- function(x, B = 200) {
  n <- length(x)
  
  lam_hat <- lambda_IQR(x)
  
  boot_estimates <- numeric(B)
  for (b in 1:B) {
    x_star <- sample(x, size = n, replace = TRUE)
    boot_estimates[b] <- lambda_IQR(x_star)
  }
  bias_hat <- mean(boot_estimates) - lam_hat
  lam_hat - bias_hat

  }

one_replication <- function(lambda, n, B = 200) {
  x <- rexp(n, rate = lambda)
  
  est_ml <- lambda_ML(x)
  est_IQR <- lambda_IQR(x)
  est_tilde <- lambda_tilde(x, B = B)

  c(
    ml = (est_ml - lambda)^2,
    iqr =(est_IQR - lambda)^2,
    hat = (est_tilde - lambda)^2
  )
}
```


``` r
mc_converg <- function(lambda, n , R_max, B = 200) {
  sq_errors <- matrix(NA, nrow = R_max, ncol = 3)
  colnames(sq_errors) <- c("ML", "IQR", "Tilde")
  
  for (r in 1:R_max) {
    sq_errors[r, ] <- one_replication(lambda, n, B = B)
  }

fall_MSE <- apply(sq_errors, 2, function(col) cumsum(col) / seq_along(col))

foll_se <- matrix(NA, nrow = R_max, ncol = 3)
colnames(fall_MSE) <- c("ML", "IQR", "Tilde")

for(j in 1:3) {
  for(r in 2:R_max) {
    v_r <- var(sq_errors[1:r, j])
    foll_se[r, j] <- sqrt(v_r / r)
  }
} 

list(sq_errors = sq_errors, fall_MSE = fall_MSE, foll_se = foll_se)
}
```


``` r
set.seed(42)

result <- mc_converg(lambda = 1, n = 5, R_max = 2000, B = 200)
```
Visual of Convergence Check simulation methods

``` r
matplot(result$fall_MSE, type = "l", lty = 1, lwd = 2,
        col = c("blue", "red", "darkgreen"),
        xlab = "Number of replications (R)", ylab = "MSE estimate",
        main = "Convergence of MSE Estimate")
legend("topright", legend = c("ML", "IQR", "Tilde"),
       col = c("blue", "red", "darkgreen"), lty = 1, lwd = 2)
```

![](Assignment-4-Exercise_files/figure-latex/unnamed-chunk-6-1.pdf)<!-- --> 
Monte Carlo SE

``` r
matplot(result$foll_se,  type = "l", lty = 1, lwd = 2,
        col = c("blue", "red", "darkgreen"),
        xlab = "Number of replications (R)", ylab = "Monte Carlo SE of MSE estimate",
        main = "Monte Carlo SE vs R")
legend("topright", legend = c("ML", "IQR", "Tilde"),
       col = c("blue", "red", "darkgreen"), lty = 1, lwd = 2)
```

![](Assignment-4-Exercise_files/figure-latex/unnamed-chunk-7-1.pdf)<!-- --> 

Full Experiment

``` r
lambda_grid <- c(0.1, 1, 10)
n_grid <- c(5, 30, 100)
R_grid <- c(`5` = 5000, `30` = 2000, `100` = 1000)

results <- data.frame()

set.seed(42)

for (lambda in lambda_grid) {
 for (n in n_grid) {
   R <- R_grid[as.character(n)]
   
   sq_errors <- matrix(NA, nrow = R, ncol = 3)
   colnames(sq_errors) <- c("ML", "IQR", "Tilde")
   
   for (r in 1:R) {
     sq_errors[r, ] <- one_replication(lambda, n, B = 200)
   }
   mse <- colMeans(sq_errors)
   
   results <- rbind(results, data.frame(
     lambda = lambda,
     n = n,
     R = R,
     MSE_ML = mse["ML"],
     MSE_IQR = mse["IQR"],
     MSE_Tilde = mse["Tilde"]
   ))
 } 
}

rownames(results) <- NULL
clean_results <- results
clean_results[ ,4:6] <- round(clean_results[ ,4:6], 4)
print(clean_results)
```

```
##   lambda   n    R  MSE_ML MSE_IQR MSE_Tilde
## 1    0.1   5 5000  0.0060 30.5247   32.4515
## 2    0.1  30 2000  0.0004 11.3234   10.6321
## 3    0.1 100 1000  0.0001  9.4610    9.2450
## 4    1.0   5 5000  0.6071  1.2056   13.1227
## 5    1.0  30 2000  0.0394  0.2949    0.3197
## 6    1.0 100 1000  0.0105  0.2567    0.2633
## 7   10.0   5 5000 50.7400 63.8940  917.6861
## 8   10.0  30 2000  3.9489  6.4778    8.5829
## 9   10.0 100 1000  1.0245  5.0457    5.6132
```


Plot Convergence Report:
ML (Blue) settling the fastest and lowest, and is the 'best' estimator when the model is correctly specified, outcome has both low bias and low variance
IQR (Red) slightly higher than ML, only taking slightly longer, relatively just as smooth.
Tilde (Green) is the highest, still have spikes even with R = 500 replications. It didn't converge properly.

Monte Carlo 
Tilde (Green) with large standard error is still large, and wasn't flatten by R = 500, while ML and IQR have stablised well much sooner.

Why the Lambda^Tilde (Green) still perform badly in both iterations, the difference in quartile 3/4 and quartile 1/4 can be small by chance. This makes the $\lambda_{IQR}$ = ln 3 / $q_{3/4(X)}$ − $𝑞_{1/4(X)}$ blow up. This bootstrap correction than compounds this, resampling from the same small, already unstable sample. Son on unfavorable replications, the correction can massively overshoot, producing occasional huge squared errors as shown as sharp spikes within the plots.

Therefore, increasing the R in theory should help the SE and potentially flatten the lines within the plot. For n = 5, lambda should hold the instability.

However, despite increasing the R = 2000, no meaningful changes is observed, and no smoothing effect. At n = 5, the Monte Carlo estimate of MSE for $\lambda$ tilde did not stablise at R = 2000, observing extreme outlying values consistent with a heavy tailed sampling distribution for this estimator in small samples. 

Final Report: 
MSE grows with n growth, at n = 100, MSE_ML at $\lambda = 0,1, 1, 10$ is roughly 0.0001, 0.0105, and 1.0245, all three values matching the MLE asymptotic variance of $\lambda^{hat}$ is $\lambda^2 / n$. This match proves simulation is behaving correctly and ML performing uniformly the best.

The $\lambda^{tilde}$ vs $\lambda^{hat}$ comparison, the bootstrap bias correction did not improve on the plain IQR estimator in this experiment. At n increases, the value of MSE performance improves, smaller n, performs horribly particularly paired with large $\lambda$




2.

``` r
library(fpp2)
```

```
## Warning: package 'fpp2' was built under R version 4.5.2
```

```
## -- Attaching packages -------------------------------------------- fpp2 2.5.1 --
```

```
## v ggplot2   4.0.3     v fma       2.5  
## v forecast  9.0.1     v expsmooth 2.3
```

```
## Warning: package 'ggplot2' was built under R version 4.5.2
```

```
## Warning: package 'forecast' was built under R version 4.5.2
```

```
## 
```

``` r
data("ausbeer")
head(ausbeer)
```

```
##      Qtr1 Qtr2 Qtr3 Qtr4
## 1956  284  213  227  308
## 1957  262  228
```

1.

``` r
autoplot(ausbeer) + 
  ggtitle("Quarterly Australian Beer Production") +
  xlab("Year") + ylab("Megalitres")
```

![](Assignment-4-Exercise_files/figure-latex/unnamed-chunk-11-1.pdf)<!-- --> 
Interpretation: 
Frequency based on observation is "Quarterly" beer production, hence for our time series, Frequency observation per year is four, and each one unit "peak" along the (X-axis) are individual years that contain one full up and down cycle made roughly of 4 data points, one rise, one peak, one fall, and one lowest minimum point. Each pattern change of the individual "peak" show us that four point per year.

2.
We can note the trend of beer production increases significantly and follow an upward trend from 1960s to 1990s, followed by gradual decline from the earl 1990s to 2010s.

3.
I can comment on the seasonal changes of the individual years, this consistent peak point and low points tells us there is significant contribution from seasonal changes that affects beer production. This consistency aligns with high quarters and low quarters each year that can be explained by seasonal demands such as summer vs winter for beer production. The size of the seasonal changes is not always constant over time, the gap between peak and lowest production grows from 1960s to the 1980s, appearing around the largest around 1980s to 1990s, and shortening again around 2010s.

4.
Variability refer to the seasonal fluctuation gap of the peak and the lowest point here, the levels seen here locally grow and shrink as the time progresses hence variability is not constant also noting the growth from 1960s and shrinking somewhat later around 2010s. A pattern consistent with variance that scales with the mean. This is exactly the kind of situation a log transformation is appropriate, and its application would help stabilise the seasonal swings across the series.

5.
Based on observation, there seems to be two points in time that beer production that reaches the highest level of around 600 megalitres in 1988 and 1991 respectively, in comparison to the peaks before and after, which sits around 540 to 580 megalitres, breaking an otherwise smooth, gradual progression of peak heights that was seen in seasonal patterns. Afterward, there is a steady decline in production, and this suggest a long term trend was coincidential in timing, events that is worth flagging, perhaps reasoning may be due to anomalies in production or consumer demand.

6.

``` r
ggseasonplot(ausbeer, polar = TRUE) +
  ylab("Megalitres") +
  ggtitle("Polar Seasonal Plot of Quarterly Beer Production in Australia")
```

![](Assignment-4-Exercise_files/figure-latex/unnamed-chunk-12-1.pdf)<!-- --> 

``` r
ggseasonplot(ausbeer, year.labels = TRUE, year.labels.left = TRUE) +
  ylab("Megalitres") +
  ggtitle("Seasonal Plot of Beer Production in Australia")
```

![](Assignment-4-Exercise_files/figure-latex/unnamed-chunk-13-1.pdf)<!-- --> 
Based on observation of the two seasonal plots above, we can see a clearer picture of seasonal patterns itself, each year follow a similar basic shape, high in Q1, lowest point in Q2, rising in Q3, peaking in Q4 and these patterns repeat consistently over the years. This advantage in over the original time series plot, where seasonal fluctuation and long term trend are combined into a single line, making it harder to isolate the shape of the season. In the seasonal pot, the trend show as a vertical stack of the colored lines by year, highest line following 1980s to 1990s, while earlier stacks in 1950s to 1960s and 2010s sit lower. This matches the patterns seen in the time series plot, but this plot separates the seasonal shape rather than mixing it.  

7.

``` r
ggsubseriesplot(ausbeer) + 
  ylab("Megalitres") + 
  ggtitle("Seasonal Subseries plot: Beer Production in Australia")
```

![](Assignment-4-Exercise_files/figure-latex/unnamed-chunk-14-1.pdf)<!-- --> 
Based on observation, we can see the means for each quarter, showing an underlying pattern in seasonal changes and shown more clearly. Q1 as second highest following the previous quarter, Q2 has the lowest mean, Q3 is rising again, and Q4 has the highest mean in production. this consistency can be seen with the seasonal plot and and following the trend of the time series plot. Looking at each quarter's panel, the level is not perfectly stable over time, each panel's line rise and falls following the same long term trend seen in the over series, climbing from 1950s to 1960s, peaking around 1980s to 1990s, and declining toward 2010s, with the largest deviations from the quarter's mean occurring during that peak time. So in terms of relative ordering of the quarters Q4 being highest, Q2 being lowest is seen consistently throughout the series plotting. The level of production within each quarter is not stable, and is similar to the broader trend and changing variability identified earlier.


8.

``` r
library(ggplot2) 
library(gridExtra)

plot1 <- autoplot(ausbeer) + 
  labs(title = "Original Series", y = "Megalitres")

plot2<- autoplot(log(ausbeer)) +
  labs(title = "Log-transformed series", y = "log(Megalitres)")

grid.arrange(plot1, plot2, ncol = 1)
```

![](Assignment-4-Exercise_files/figure-latex/unnamed-chunk-15-1.pdf)<!-- --> 
Detrending is worth doing due to the evidence seen in time series plot and the seasonal plot color stacking. The time series show a clear non-stationary trend, rising sharply during 1950s to 1990s, then declining in 2010s. The season plot backs this with clearer stacking visuals, each line by year. The time series models assume a constant mean level roughly. This series doesn't have one, so detrending is needed before modelling.

Transforming is also worth doing due to our earlier evidence in variability, the seasonal spiking isn't constant, its small in 1960s, growing towards 1990s, and declining again later. This growth swinging tracks the growth in the series overall level over the same period. The sub series plot making this visible during each quarter's panel, showing the deviation from the horizontal line. This spike swinging increasing alongside the level is where a log transformation is commonly used and help with stabilisation of the variance, making it more constant regardless of the series level.

Based on the observation of the two plots above, we can see the swing in the original series, the same pattern spiking in the 600 megalitres being produced, the growth from 1960s to 1990s and slow decline to 2010s. 

However the log series swing in the 1980s to 1990s spike is relatively smaller to the original. This change in log transformation while not the most dramatic in change, were worth being investigated, and proved a slight improvement in stablisation of the variance, while not fully getting rid of the changing spread on its own.
