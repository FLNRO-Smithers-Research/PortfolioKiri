//Kiri Daust
//Cpp function for portfolio optimisation


#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector simGrowthCpp(DataFrame DF){
	NumericVector Returns(101);
	int nTrees = 100;
	NumericVector Growth = DF["Growth"]; //convert to vectors
	NumericVector NoMort  = DF["NoMort"];
	NumericVector MeanDead = DF["MeanDead"];
	double height, percentDead;
	int prevTrees, numDead, i;
	for(i = 0; i < 101; i++){
		height = sum(Growth[Rcpp::Range(0,i)]);
		Returns[i] = nTrees*height;
		if(Rcpp::runif(1,0,100)[0] > NoMort[i]){
			prevTrees = nTrees;
			percentDead = Rcpp::rgamma(1, 1, MeanDead[i])[0];
			numDead = (percentDead/100)*prevTrees;
			nTrees = prevTrees - numDead;
		}
	}
	return(Returns);
}
// // [[Rcpp::export]]
// NumericVector calcAssetsCpp(NumericVector Returns){
//   NumericVector assets(101);
//   assets[0] = Returns[0];
//   int z;
//   for(z = 0; z < 101; z++){
//     assets[z+1] = Returns[z+1] - Returns[z];
//   }
//   return(assets);
// }