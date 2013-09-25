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

#Implements a function for using quantile normalization for
#cross platform normalization.



qn = function(platform1.data, platform2.data, p1.names=0, p2.names=0, skip.match=FALSE){
	
	
	#Match names
	input = processplatforms(list(x=platform1.data,y=platform2.data),namesvec = c(p1.names, p2.names), skip.match=skip.match)
	
	#Get dimensions
	nx = dim(input$x)[2]
	ny = dim(input$y)[2]
	
	#combine the data
	combined = as.matrix(cbind(input$x,input$y))
	
	#Do it
	out = as.data.frame(normalize.quantiles(combined, copy=FALSE))
	
	#Return the data
	return(list(x=out[,1:nx],y=out[,(nx+1):(nx+ny)]))
	
}

