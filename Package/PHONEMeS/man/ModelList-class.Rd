\name{ModelList-class}
\alias{ModelList}
\alias{ModelList-class}
\docType{class}
\title{
Sampled models and scores
}
\description{
This class contains a model sampled during optimisation, and associated scores
}
\section{Objects from the Class}{ 
Objects of class \code{ModelList} can be created from \code{
new("ModelList", model=model, scores=scoresT, scoresList=scores.dataOn, intAnd=intermediates.AND.or, intgAnd=integrators.AND.or)}
}
\section{Slots}{
	\describe{	
		\item{\code{model}:}{
		data frame, subset of an \code{interactions} element from a \code{PKNlist}
		}
		
		\item{\code{scores}:}{
		numeric, total score for a network
		}
		
		\item{\code{scoresList}:}{
		list, scores by conditions and treatments
		}

		\item{\code{intAnd}:}{
		logical, for each intermediate, are the incoming edges combined with an AND or not
		}

		\item{\code{intgAnd}:}{
		logical, for each integrator, are the incoming edges combined with an AND or not
		}	
	}
}
\section{Methods}{
	An overview of methods with class-specific functionality:
	More detailed introduction can be found in help for each specific
	function.
	\describe{	
		\item{\code{run1model}}{ 
		samples and runs one model
		}
	}
}
\author{C. Terfve \email{terfve@ebi.ac.uk}}
\seealso{
  \code{\link[PHONEMeS: run1model]{run1model}}
}
\examples{
}
\keyword{classes}
