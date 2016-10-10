\name{GMMbyCond-class}
\alias{GMMbyCond}
\alias{GMMbyCond-class}
\docType{class}
\title{
Processed MS data, organised by condition
}
\description{
This class contains the MS data to be used, after pre-processing and Gaussian Mixture Modelling, organised by condition
}
\section{Objects from the Class}{ 
Objects of class \code{GMMbyCond} can be created from \code{dataBycond}
}
\section{Slots}{
	\describe{	
		\item{\code{dataBC}:}{
		list with one element for each condition, and inside one data frame (data for all peptides under one drug) for each drug treatment that matches that condition
		}
		
		\item{\code{speciesP}:}{
		list matching to the above, containing the nodes in the data (in the same format as S.cc, i.e. "UPID.res.pos") that are perturbed under the drug treatments used. The data in both of these elements only relates to those nodes that are reachable in the full background network.
		}	
	}
}
\section{Methods}{
	An overview of methods with class-specific functionality:
	More detailed introduction can be found in help for each specific
	function.
}
\author{C. Terfve \email{terfve@ebi.ac.uk}}
\seealso{
  \code{\link[PHONEMeS:dataBycond]{dataBycond}}
}
\examples{
}
\keyword{classes}
