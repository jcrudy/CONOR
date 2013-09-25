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
#This file contains code adapted from the NorDi program available from 
#http://keia.i3s.unice.fr/?Logiciels:Nordi.  Please cite appropriately 
#when using these functions.  Type citation("CONOR") at the R prompt 
#for details.


nordi = function(platform1.data, platform2.data, pvalue = .01, alpha = .05, p1.names=0, p2.names=0, skip.match=FALSE){
	
	
	#Match names
	input = processplatforms(list(x=platform1.data,y=platform2.data),namesvec = c(p1.names, p2.names), skip.match=skip.match)
	
	sampnamesx <- colnames(input$x)
	sampnamesy <- colnames(input$y)
	
	genenamesx <- rownames(input$x)
	genenamesy <- rownames(input$y)
	
	#Do the nordi
	input$x = nordi1p(input$x,pvalue,alpha)
	input$y = nordi1p(input$y,pvalue,alpha)
	
	sampnamesx -> colnames(input$x)
	sampnamesy -> colnames(input$y)

	genenamesx -> rownames(input$x)
	genenamesy -> rownames(input$y)

	#Return the values
	return(list(x=as.data.frame(input$x),y=as.data.frame(input$y)))
}

getoutliers<-function(D,pvalue){
	Out=TRUE;
	Noa=TRUE;
	while (Out && Noa){
		J<-jarque.bera.test(D);
		J0<-J$statistic;
		G<-grubbs.test(D,type=10);
		DBack<-D;
		if (G$p.value<pvalue) {
			D<-rm.outlier(D, fill = FALSE, median = FALSE)
		}else{
			Out<-FALSE
		};
	
		J1<-jarque.bera.test(D);
		J1<-J1$statistic;
		if (J0-J1>=0) {
			Noa<-TRUE
		}else{
			D<-DBack; 
			Noa<-FALSE
		};
	}; 
	return(D)
};

nordi1p = function(M, pvalue = .01, alpha = .05){
	#The pvalue determines the sensitivity for detecting outliers 
	#in each column of the gene expression matrix.
	#The alpha value determines the size of the tails of the 
	#normal distributions considered to be over or under expressed.

	M<-as.matrix(M);
	C<-split(M,col(M));


	Ot<- c(0*1:length(C));
	Ut<- c(0*1:length(C));
	for (i in 1:length(C)){
		C[[i]]<-na.omit(C[[i]]);
		C[[i]]<-getoutliers(C[[i]],pvalue);
		m<-mean(C[[i]]);
		stdev<-sqrt(var(C[[i]]));
		zcutoff = -1*qnorm(alpha/2)
		Ot[i]<-m + zcutoff*stdev;
		Ut[i]<-m - zcutoff*stdev;
	}
	Mfin<-M*0;    
	for(i in 1:length(colnames(M))){
		for (j in 1:length(rownames(M))){
			ifelse(M[j,i]<Ut[i],Mfin[j,i]<-(-1),ifelse(M[j,i]>Ot[i],Mfin[j,i]<-(1),Mfin[j,i]<-0));
  		}
	}
	return(Mfin)
}