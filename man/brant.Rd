\name{brant}
\alias{brant}

\title{
Brant Test
}
\description{
The function calculates the brant test for  ordinal logit models to test the parallel regression assumption.
}
\usage{
brant(model,by.var=F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
the polr-Object generated with polr()
}
  \item{by.var}{
OPTIONAL if set to true, the tests are made for each variable instead of each coefficient. Default: FALSE.
}
}
\details{
  The function calculates the brant test for parallel regression assumption. The brant test was published by Brant (1990). The function works with models generated with the function polr() from the package 'MASS'.
}
\value{
  The output is the brant test, which shows if the parallel assumption holds or not.
}

\author{
  Benjamin Schlegel, \email{kontakt@benjaminschlegel.ch}
}

\examples{
data = MASS::survey
data$Smoke = ordered(MASS::survey$Smoke,levels=c("Never","Occas","Regul","Heavy"))
model1 = MASS::polr(Smoke ~ Sex + Height, data=data, Hess=TRUE)
summary(model1)
brant(model1)
}

\references{
Brant, R. (1990) Assessing proportionality in the proportional odds model for ordinal logistic regression. \emph{Biometrics},  \bold{46}, 1171--1178.
}

\keyword{ models, brant, brant test, parallel regression assumption, ordinal logit }
