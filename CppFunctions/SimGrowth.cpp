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
	int numYears = Growth.length();
	NumericVector Returns(numYears);
	double height, percentDead;
	int prevTrees, numDead, i;
	for(i = 0; i < numYears; i++){
		height = sum(Growth[Rcpp::Range(0,i)]);
		Returns[i] = nTrees*height;
		if(Rcpp::runif(1,0,100)[0] > NoMort[i]){//regular environmental loss
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
NumericVector SimGrowth_v2(DataFrame DF, double ProbPest, double cmdMin, 
                           double cmdMax, double tempMin, double tempMax, double climLoss){
  NumericVector Growth = DF["Growth"]; //convert to vectors
  NumericVector NoMort  = DF["NoMort"];
  NumericVector MeanDead = DF["MeanDead"];
  NumericVector Ruin = DF["Suit"];//think about this one
  NumericVector climCMD = DF["CMD"];
  NumericVector climMax = DF["Tmax_sm"];
  NumericVector climMin = DF["Tmin_sp"];
  
  int numYears = Growth.length();
  NumericVector Returns(numYears);
  double height, percentDead, percentRuin, climDead, climDiff, diffProp;
  int prevTrees, numDead, i;
  int nTrees = 100;
  for(i = 0; i < numYears; i++){
    height = sum(Growth[Rcpp::Range(0,i)]);
    Returns[i] = nTrees*height;
    climDead = 0;
    if(climCMD[i] > cmdMax){ // CMD max 
      climDiff = climCMD[i] - cmdMax;
      diffProp = (climDiff*100)/(cmdMax - cmdMin);
      climDiff = 1 - exp(-climLoss*diffProp);//1 - exp function
      climDead += climDiff*nTrees;
      //Rcout << "Too dry \n";
    }
    if(climCMD[i] < cmdMin){
      climDiff = cmdMin - climCMD[i];
      diffProp = (climDiff*100)/(cmdMax - cmdMin);
      climDiff = 1 - exp(-climLoss*diffProp);//1 - exp function
      climDead += climDiff*nTrees;
      //Rcout << "Too wet \n";
    }
    if(climMax[i] > tempMax){
      climDiff = climMax[i] - tempMax;
      diffProp = (climDiff*100)/(tempMax - tempMin);
      climDiff = 1 - exp(-climLoss*diffProp);//1 - exp function
      climDead += climDiff*nTrees;
      //Rcout << "Too hot \n";
    }
    if(climMin[i] < tempMin){
      climDiff = tempMin - climMin[i];
      diffProp = (climDiff*100)/(cmdMax - cmdMin);
      climDiff = 1 - exp(-climLoss*diffProp);//1 - exp function
      climDead += climDiff*nTrees;
      //Rcout << "Too cold \n";
    }
    if(climDead > nTrees){
      climDead = nTrees;
    }
    nTrees = nTrees - climDead;
    if(Rcpp::runif(1,0,1)[0] <= ProbPest){//pest outbreak
      percentRuin = rgamma(1, Ruin[i], 0.07)[0];
      if(percentRuin > 1){
        percentRuin = 1;
      }
      numDead = percentRuin*nTrees;
      nTrees = nTrees - numDead;
    }else if(Rcpp::runif(1,0,100)[0] > NoMort[i]){//regular environmental loss
      //prevTrees = nTrees;
      percentDead = Rcpp::rgamma(1, 1.5, MeanDead[i])[0];
      //Rcout << "Percent Dead:" << percentDead << "\n";
      numDead = (percentDead/100)*prevTrees;
      nTrees = nTrees - numDead;
    }
  }
  // NumericVector out = NumericVector::create(_["Volume"] = Returns[numYears-1], 
  //                                           _["NumTrees"] = nTrees);
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

// // [[Rcpp::export]]
// NumericVector SimGrowth_Volume(DataFrame DF){
//   int nTrees = 100;
//   NumericVector Growth = DF["Growth"]; //convert to vectors
//   NumericVector NoMort  = DF["NoMort"];
//   NumericVector MeanDead = DF["MeanDead"];
//   NumericVector Ruin = DF["RuinSeverity"];
//   int numYears = Growth.length();
//   NumericVector Returns(numYears);
//   double height, percentDead, percentRuin;
//   int prevTrees, numDead, i;
//   for(i = 0; i < numYears; i++){
//     height = sum(Growth[Rcpp::Range(0,i)]);
//     Returns[i] = nTrees*height;
//     if(Rcpp::runif(1,0,100)[0] > NoMort[i]){//regular environmental loss
//       prevTrees = nTrees;
//       percentDead = Rcpp::rgamma(1, 1.5, MeanDead[i])[0];
//       if(percentDead > 5){
//         percentRuin = rnorm(1,Ruin[i],0.1)[0];
//         if(percentRuin < 0){
//           percentRuin = 0;
//         }
//         numDead = percentRuin*prevTrees;
//         nTrees = prevTrees - numDead;
//       }else{
//         //Rcout << "Percent Dead:" << percentDead << "\n";
//         numDead = (percentDead/100)*prevTrees;
//         nTrees = prevTrees - numDead;
//       }
//     }
//   }
//   NumericVector out = NumericVector::create(_["Volume"] = Returns[numYears-1], 
//                                             _["NumTrees"] = nTrees);
//   return(out);
// }