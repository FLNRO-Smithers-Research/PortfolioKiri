//Kiri Daust 2019
//Cpp function for portfolio optimisation


#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector simGrowthCpp(DataFrame DF){
	int nTrees = 100;
	NumericVector Growth = DF["Growth"]; //convert to vectors
	NumericVector NoMort  = DF["NoMort"];
	NumericVector MeanDead = DF["MeanDead"];
	NumericVector Ruin = DF["RuinSeverity"];
	int numYears = Growth.length();
	NumericVector Returns(numYears);
	double height, percentDead;
	int prevTrees, numDead, i;
	for(i = 0; i < numYears; i++){
		height = sum(Growth[Rcpp::Range(0,i)]);
		Returns[i] = nTrees*height;
		if(Rcpp::runif(1,0,110)[0] >= 109){//drastic pest loss
		  percentDead = rnorm(1,Ruin[i],0.1)[0];
		  if(percentDead < 0){
		    percentDead = 0;
		  }
		  numDead = percentDead*prevTrees;
		  nTrees = prevTrees - numDead;
		}else if(Rcpp::runif(1,0,100)[0] > NoMort[i]){//regular environmental loss
			prevTrees = nTrees;
			percentDead = Rcpp::rgamma(1, 1.5, MeanDead[i])[0];
			//Rcout << "Percent Dead:" << percentDead << "\n";
			numDead = (percentDead/100)*prevTrees;
			nTrees = prevTrees - numDead;
		}
	}
	return(Returns);
}

// [[Rcpp::export]]
NumericVector simGrowthCBST(DataFrame DF){
  int nTrees = 100;
  NumericVector Growth = DF["Growth"]; //convert to vectors
  NumericVector NoMort  = DF["NoMort"];
  NumericVector MeanDead = DF["MeanDead"];
  int numYears = Growth.length();
  NumericVector Returns(numYears);
  double height, percentDead;
  int prevTrees, numDead, i;
  for(i = 0; i < numYears; i++){
    height = sum(Growth[Rcpp::Range(0,i)]);
    Returns[i] = nTrees*height;
    if(Rcpp::runif(1,0,100)[0] > NoMort[i]){
      prevTrees = nTrees;
      percentDead = Rcpp::rgamma(1, 1, MeanDead[i]/10)[0];
      numDead = (percentDead/100)*prevTrees;
      nTrees = prevTrees - numDead;
    }
  }
  return(Returns);
}

// [[Rcpp::export]]
NumericVector gs2gw(NumericVector x, double a, double b){
  int len = x.length();
  NumericVector out(len);
  for(int i = 0; i < len; i++){
    out[i] = a*exp(x[i]*b);
  }
  return(out);
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