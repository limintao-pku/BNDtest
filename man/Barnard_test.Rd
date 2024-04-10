\name{Barnard_test}
\alias{Barnard_test}
\title{
Barnard's exact test for 2 by 2 contingency table
}
\description{
This function is designed to implement Barnard's original algorithm (CSM) for relatively large tables.
}

\usage{
Barnard_test(mat, alternative=c("two.sided","greater","less"), n_p=101L, trace=F,
             return_data=F, data_previous=NULL)
}

\arguments{
\item{mat}{
A 2 by 2 contingency table (matrix). \code{mat[1,1]} and \code{mat[2,1]} denote the numbers of cases among the exposed and unexposed, respectively.
}
\item{alternative}{
The alternative hypothesis. Must be one of "two.sided", "greater" or "less".
}
\item{n_p}{
The number of nuisance parameters considered. This argument is passed to \code{seq(0, 1, length.out = n_p)}.
}
\item{trace}{
A logical argument indicating whether to track the progress. Tracing will increase the time consumed.
}
\item{return_data}{
A logical argument indicating whether to return the data that can be passed to \code{data_previous}.
}
\item{data_previous}{
Data output by this function (\code{return_data = T}) that can be used to accelerate the test for another equally sized contingency table (sample sizes of both the exposure groups are the same as before).
}
}

\details{
Barnard's exact test is an unconditional test for the association between a binary exposure and binary outcome. It is generally more powerful than Fisher's exact test when the data are collected in unconditional manners. The original Barnard's algorithm (CSM) is computationally intensive. This package is designed to deal with relatively large 2 by 2 contingency tables using the CSM method. If both the sample sizes of the exposure groups are large, use \code{\link[stats]{chisq.test}} if appropriate. More details can be found in this paper (Barnard, 1947) and the help document of \code{\link[Exact]{exact.test}} in the 'Exact' R package.
}

\value{
If \code{return_data = F}, the \emph{P} value is directly returned. Otherwise, a list containing the \emph{P} value and other data will be returned.
}

\note{
We have used this function mainly for unbalanced 2 by 2 tables. The sample size of the smaller exposure group is < 200, and that of the larger group is ~ 2000. The algorithm typically finishes within seconds or minutes. This function can also be used for smaller or larger tables.
}

\references{
Barnard, G.A. (1947). Significance tests for 2x2 tables. Biometrika 34, 123-138.
}

\seealso{
\code{\link[Exact]{exact.test}}, \code{\link[stats]{chisq.test}}, and \code{\link[stats]{fisher.test}}.
}
\examples{
#A small table
mat <- rbind(c(5, 0), c(7, 8))
Barnard_test(mat, trace=F, return_data=F)
fisher.test(mat)$p.value

#A large unbalanced table
mat <- rbind(c(5, 46), c(63, 2122))
test_out1 <- Barnard_test(mat, trace=T, return_data=T)
test_out1$p.value

mat <- rbind(c(8, 43), c(65, 2120))
test_out2 <- Barnard_test(mat, trace=T, return_data=T, data_previous=test_out1$data)
test_out2$p.value

test_out3 <- Barnard_test(mat, trace=T, return_data=T, data_previous=NULL)
test_out3$p.value
}