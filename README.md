# BNDtest
## Barnard's exact test for large unbalanced 2 by 2 contingency tables
Barnard's exact test is an unconditional test for the association between a binary exposure and binary outcome. It is generally more powerful than Fisher's exact test when the data are collected in unconditional manners. The original Barnard's algorithm (CSM) is computationally intensive. This package is designed to deal with relatively large 2 by 2 contingency tables using the CSM method.

To install this R package, run `install.packages("remotes"); remotes::install_github("limintao-pku/BNDtest")`.  
To load this package, run `library(BNDtest)`.  
See how to use this package by running `?Barnard_test`.  
