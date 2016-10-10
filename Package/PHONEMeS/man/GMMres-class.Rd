\name{GMMres-class}
\alias{GMMres}
\alias{GMMres-class}
\docType{class}
\title{
Processed MS data 
}
\description{
This class contains the MS data to be used, after pre-processing and Gaussian Mixture Modelling 
}
\section{Objects from the Class}{ 
Objects of class \code{GMMres} can be created from \code{
new("GMMres", res=GMM.res, IDmap=GMM.res.ID)}
}
\section{Slots}{
	\describe{	
		\item{\code{res}:}{
		result of the Gaussian mixture modelling process, i.e. a list with an element for each peptide, each element containing a matrix with: the P/C cluster assignment ("clus"), log likelihood of P/C ("Indiv"), adjusted p value of FC drug vs ctrl ("FCvCaPval") and OK/FP/FN/I status ("status")
		}
		
		\item{\code{IDmap}:}{
		a data frame  that maps the data identifiers (column "dataID", names of the elements in the \code{res} list) to the network identifiers (column "S.cc", UPID.res.pos)
		}
		
		\item{\code{resFC}:}{
		same as \code{res} but with an extra column in each matrix, containing the value of the fold changes
		}	
	}
}
\section{Methods}{
	An overview of methods with class-specific functionality:
	More detailed introduction can be found in help for each specific
	function.
	\describe{	
		\item{\code{dataBycond}}{ 
		reorganise the data to have a matrix for each drug treatment, as elements of a list with an entry for each condition (set of drug targets that are shared between different treatments)
		}
	}
}
\author{C. Terfve \email{terfve@ebi.ac.uk}}
\seealso{
  \code{\link[PHONEMeS:dataBycond]{dataBycond}}
}
\examples{
}
\keyword{classes}
