######################################################################
#Copyright Jason Rudy & Faramarz Valafar 2009-2010

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
######################################################################

######################################################################
#This file contains code adapted from the WebArrayDB program,
#available from http://www.webarraydb.org/webarray/index.html  Please 
#cite appropriately when using these functions.  Type 
#citation("CONOR") at the R prompt for details.


#Implements the median rank scores algorithm for cross
#platform normalization.

mrs = function(platform1.data, platform2.data, p1.names=0, p2.names=0, skip.match = FALSE){
	#This function is basically a wrapper for normalizeGQ
	
	#Match names
	input = processplatforms(list(x=platform1.data,y=platform2.data),namesvec = c(p1.names, p2.names), skip.match=skip.match)
	
	#Prepare for normalizeGQ
	combined = cbind(input$x,input$y)
	pf = c(seq(1,1,length.out=dim(input$x)[2]),seq(2,2,length.out=dim(input$y)[2]))
	
	#Call normalizeGQ
	nmrs = norm_mrs(combined,pf)
	
	#Split the results and return
	out=split(seq(pf),pf)
	out[[1]] = nmrs[,out[[1]]]
	out[[2]] = nmrs[,out[[2]]]
	names(out) <- c("x","y")
	return(out)
}



norm_mrs <- function(M, pf, ...) { # median rank  scores algorithm  - a modification of quantile for multi-platform data
	idx <- tapply(seq(pf), pf, function(x) x)
	if (length(pf)<=0) return(M)
	#if (length(idx)==1) return(normalizeQuantiles(M))
	imax <- which.max(sapply(idx, length)) # for reference
	ref_med_srt <- sort(apply(M[, idx[[imax]]], 1, function(x) median(x, na.rm=TRUE) ))
	idx[imax] <- NULL
	lapply(unlist(idx), function(i) M[,i] <<- ref_med_srt[rank(M[,i])])
	invisible(M)
	}
	
	
	
	
