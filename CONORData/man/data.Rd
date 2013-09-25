\name{data}
\alias{MAQC.ABI}
\alias{MAQC.AFX}
\alias{MAQC.AFX.transfer}
\alias{MAQC.AG1}
\alias{MAQC.ILM}
\alias{MAQC.ILM.transfer}
\alias{SPERM.AFX}
\alias{SPERM.AFX.transfer}
\alias{SPERM.ILM}
\alias{SPERM.ILM.transfer}



%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Example Data Sets for Package CONOR
}
\description{
The data sets described below are included in the CONOR package.
}

\usage{
data(MAQC.ABI)
data(MAQC.AFX)
data(MAQC.AFX.transfer)
data(MAQC.AG1)
data(MAQC.ILM)
data(MAQC.ILM.transfer)
data(SPERM.AFX)
data(SPERM.AFX.transfer)
data(SPERM.ILM)
data(SPERM.ILM.transfer)
}

\arguments{
\item{MAQC.ABI}{A data.frame containing normalized Applied Biosystems (ABI) data from the Microarray Quality Control Project (MAQC).  The ABI data were downloaded already normalized from GEO using the GEOquery package from Bioconductor. These data were normalized by the MAQC authors using the Expression Array System Software suite from Applied Biosystems, which implements a platform specific normalization sequence based on the specific properties of the ABI array and the 1700 Chemiluminescent Microarray Analyzer. The steps in this normalization sequence take advantage of co-localized control probes and signal to noise ratios obtained during image processing. Details can be found in the supplemental materials of the MAQC publication and in the document entitled "User Bulletin: Applied Biosystems 1700 Chemiluminescent Microarray Analyzer" issued by ABI. The downloaded data were natural log transformed. The ABI data set was the only set not subjected to quantile normalization.}

\item{MAQC.AFX}{A data.frame containing normalized Affymetrix (AFX) data from the Microarray Quality Control Project (MAQC). Raw Affymetrix CEL files were pre-processed using the function justRMA from the affy package of Bioconductor.}

\item{MAQC.AFX.transfer}{The same data set as MAQC.AFX, but normalized together with the SPERM.AFX data set.}

\item{MAQC.AG1}{A data.frame containing normalized Agilent one color (AG1) data from the Microarray Quality Control Project (MAQC). The raw Agilent data were processed using the limma package from Bioconductor. Background correction was performed using the "normexp" and "mle" options of the backgroundCorrect function. Quantile normalization was performed by the normalizeBetweenArrays function. Data from duplicate probes were then averaged and natural log transformed.}

\item{MAQC.ILM}{A data.frame containing normalized Illumina data from the Microarray Quality Control Project (MAQC). Raw Illumina data were acquired from GEO in text format. The mean signal for each probe type was extracted and subjected to quantile normalization (provided by the normalize.quantiles function from Bioconductor's preprocessCore package) and natural log transformation.}

\item{MAQC.ILM.transfer}{The same data set as MAQC.ILM, but normalized together with the SPERM.ILM data set.}

\item{SPERM.AFX}{A data.frame containing normalized Affymetrix (AFX) data from the Platts et al. human sperm data set.  Data were downloaded from GEO and processed and normalized in the same way as MAQC.AFX.}

\item{SPERM.AFX.transfer}{The same data set as SPERM.AFX, but normalized together with the MAQC.AFX data.}

\item{SPERM.ILM}{A data.frame containing normalized Illumina (ILM) data from the Platts et al. human sperm data set.  Data were downloaded from GEO and processed and normalized in the same way as MAQC.ILM.}

\item{SPERM.ILM.transfer}{The same data set as SPERM.ILM, but normalized together with the SPERM.ILM data.}

}

\references{
Platts A, Dix D, Chemes H, Thompson K, Goodrich R, Rockett J, Rawe V, Quintana S, Diamond M, Strader L: Success and failure in human spermatogenesis as revealed by teratozoospermic RNAs. Human Molecular Genetics 2007, 16(7):763.

Shi et al. The MicroArray Quality Control (MAQC) project shows inter- and intraplatform reproducibility of gene expression measurements. Nature Biotechnology 2006, 24(9):1151-1161.

Davis S, Meltzer P: GEOquery: a bridge between the Gene Expression Omnibus (GEO) and BioConductor. Bioinformatics 2007, 23(14):1846.

Smyth G: Linear models and empirical Bayes methods for assessing differential expression in microarray experiments. Statistical Applications in Genetics and Molecular Biology 2004, 3.

Gautier L, Cope L, Bolstad B, Irizarry R: Affy: analysis of Affymetrix GeneChip data at the probe level. Bioinformatics 2004, 20(3):307.

Applied Biosystems: User Bulletin: Applied Biosystems 1700 Chemiluminescent Microarray Analyzer 2005.
}

\examples{
data(MAQC.ABI)
}

\keyword{microarray}
\keyword{gene}
\keyword{expression}
\keyword{normalization}
\keyword{cross-platform}
\keyword{inter-platform}