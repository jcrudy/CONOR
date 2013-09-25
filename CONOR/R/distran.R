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

#Implementation of the distran cross platform
#normalization algorithm


distran = function(platform1.data, platform2.data, p1.assayclasses=NULL, L=4, cluster="pam",corr="pearson", p1.names=0, p2.names=0, skip.match=FALSE){
	#platform1.data is used as a reference set, which is why p1.sampleinds is present.  
	#p1.assayclass should be a vector with the same number of elements as there are
	#assays in platform1.data.  If not given, the same clustering procedure used by xpn
	#will be used to automatically compute assay clusters, giving an estimate of which 
	#assays come from the same sample or sample class.  Additional parameters apply to 
	#the xpn sample clustering procedure.
	
	#Match names
	input = processplatforms(list(x=platform1.data,y=platform2.data),namesvec = c(p1.names, p2.names), skip.match=skip.match)
	
	genenames <- rownames(input$x)
	colnamesx <- colnames(input$x)
	colnamesy <- colnames(input$y)

	#dimensions
	m = dim(input$x)[1]
	nx = dim(input$x)[2]
	ny = dim(input$y)[2]
	
	#Cluster if necessary
	if (is.null(p1.assayclasses)){
		xpnclusters = xpnassaycluster(input$x,input$y,L=L,cluster=cluster,corr=corr)
		p1.assayclasses = xpnclusters[1:nx]
	}

	
	#Construct the reference data
	classes = unique(p1.assayclasses)
	nclasses = length(classes)
	refdat = numeric(m)
	for (class in classes){
		idx = p1.assayclasses==class
		refdat = refdat + (1/nclasses)*rowMeans(as.matrix(input$x[,idx],m,length(idx)))
	}
	refdat = sort(refdat)
	
	#Get the ranks of the data
	input$x = colRanks(input$x)
	input$y = colRanks(input$y)
	
	#Do the substitution
	input$x = as.matrix(input$x)
	input$y = as.matrix(input$y)
	for (i in 1:(m*ny)){
		input$y[i] = refdat[input$y[i]]
	}
	for (i in 1:(m*nx)){
		input$x[i] = refdat[input$x[i]]
	}
	
	colnames(input$x) <- colnamesx
	colnames(input$y) <- colnamesy
	rownames(input$x) <- genenames
	rownames(input$y) <- genenames

	
	#All done!
	return(list(x=data.frame(input$x),y=data.frame(input$y)))
	
}

colRanks = function(aDataFrame){
	return(apply(aDataFrame,2,rank))
}
