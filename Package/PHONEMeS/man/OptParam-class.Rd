\name{OptParam-class}
\alias{OptParam}
\alias{OptParam-class}
\docType{class}
\title{
PHONEMeS analysis parameters 
}
\description{
This class contains all parameters for an analysis
}
\section{Objects from the Class}{ 
Objects of class \code{OptParam} can be created from \code{
new("OptParam", sizeP=0, nG1=100, cstart=TRUE, intgAsintm=FALSE, nScripts=50, absTol=FALSE, tol=0.15, cap=20, resN=2, nG=50)}
}
\section{Slots}{
	\describe{	
		\item{\code{sizeP}:}{
		numeric, size penalty factor
		}
		
		\item{\code{nG1}:}{
		numeric, number of models sampled in each script
		}
		
		\item{\code{cstart}:}{
		logical, should the direct interactions between drug targets and data sites be added to every model
		}
		\item{\code{intgAsintm}:}{
		logical, TRUE if we want integrators to be sampled the same way intermediates are
		}
		\item{\code{nScripts}:}{
		numeric, number of parallel scripts that are run at each generation
		}
		\item{\code{absTol}:}{
		logical, should absolute tolerance be used
		}
		\item{\code{tol}:}{
		numeric, tolerance
		}
		\item{\code{cap}:}{
		numeric, maximum number of times that an edge is copied through at each generation
		}
		\item{\code{resN}:}{
		index of the optimisation (for naming files and folders)
		}
		\item{\code{nG}:}{
		number of generations that the optimisation is run for
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
}
\examples{
}
\keyword{classes}
