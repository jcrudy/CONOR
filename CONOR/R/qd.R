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


#Implements the quantile discretization algorithm for 
#cross platform normalization.


norm_qd <- function(M, nbin, ...) { # quantile discretization algorithm
	#lapply(seq(ncol(M)), function(i) {ra <- rank(M[,i]); rlt<-as.integer(ra/(length(ra)/nbin)); M[,i] <- rlt - median(rlt)} )
	#lapply(seq(ncol(M)), function(i) {ra <- rank(M[,i]); rlt<-ceiling(ra/(length(ra)/nbin)); M[,i] <- rlt - median(rlt)} )
	ref <- c(-1*((nbin/2)-1):0, 0:((nbin/2)-1))
	ref <- ref[ceiling(seq(nrow(M))/(nrow(M)/nbin))]
	lapply(seq(ncol(M)), function(i) M[,i] <<- ref[rank(M[,i])])
	invisible(M)
}

qd = function(platform1.data, platform2.data, b=8, p1.names=0, p2.names=0, skip.match=FALSE) {
	
	#Match names
	input = processplatforms(list(x=platform1.data,y=platform2.data),namesvec = c(p1.names, p2.names), skip.match=skip.match)
	
	return(list(x=norm_qd(input$x, nbin=b),y=norm_qd(input$y,nbin=b)))
}
