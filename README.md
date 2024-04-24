# BNDtest
## Barnard's exact test for large unbalanced 2 by 2 contingency tables
Barnard's exact test is an unconditional test for the association between a binary exposure and binary outcome. It is generally more powerful than Fisher's exact test when the data are collected in unconditional manners. The original Barnard's algorithm (CSM) is computationally intensive. This package is designed to deal with relatively large 2 by 2 contingency tables using the CSM method.  

To install this R package, run `install.packages("remotes"); remotes::install_github("limintao-pku/BNDtest")`.  
To load this package, run `library(BNDtest)`.  
See how to use this package by running `?Barnard_test`.  

There are other alternative algorithms to determine more extreme tables (e.g., Z-pooled, Z-unpooled, Santner & Snell, and Boschloo; see the "Exact" R package; https://CRAN.R-project.org/package=Exact). The CSM method is often the most powerful.  

We have used this package mainly for unbalanced 2 by 2 tables. The sample size of the smaller exposure group is < 200, and that of the larger group is ~ 2000. The algorithm typically finishes within seconds or minutes. This function can also be used for smaller or larger tables.  
