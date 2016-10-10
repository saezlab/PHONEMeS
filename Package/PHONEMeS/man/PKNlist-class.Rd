\name{PKNlist-class}
\alias{PKNlist}
\alias{PKNlist-class}
\docType{class}
\title{
K/P-S network 
}
\description{
This class contains the background of K/P-S interactions networks
}
\section{Objects from the Class}{ 
Objects of class \code{PKNlist} can be created from \code{
new("KPSbg", interactions=allD, species=unique(c(allD$K.ID, allD$S.cc)))}
}
\section{Slots}{
	\describe{	
		\item{\code{interactions}:}{
		an interactions data frame
		}
		
		\item{\code{interactionsD}:}{
		a subset of \code{interactions} that only contains experimentally supported direct interactions between drug targets and data points
		}
		
		\item{\code{species}:}{
		a character vector of nodes
		}	

		\item{\code{sinks}:}{
		a character vector of sink nodes
		}
		\item{\code{integrators}:}{
		a character vector of integrator nodes
		}
		\item{\code{intermediates}:}{
		a character vector of intermediate nodes
		}

	}
}
\section{Methods}{
	An overview of methods with class-specific functionality:
	More detailed introduction can be found in help for each specific
	function.
	\describe{	
		\item{\code{buildPKN}}{ 
		build a PKN for a specific data set, from a full background network
		}
		\item{\code{PKNlist}}{ 
		create the PKNlist from a PKN created by buildPKN
		}
	}
}
\author{C. Terfve \email{terfve@ebi.ac.uk}}
\seealso{
  \code{\link[PHONEMeS: buildPKN]{buildPKN}}
  \code{\link[PHONEMeS: PKNlist]{PKNlist}}
}
\examples{
}
\keyword{classes}
