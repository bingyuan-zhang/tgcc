#include <RcppArmadillo.h>
#include <tuple>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
void updateBoundaries(
    const double nodeValue,
    const double nodeWeight,
    const double lambda,
    const std::deque<double>& boundaryPoints,
    const std::deque<double>& intercepts,
    const std::deque<double>& slopes) {
  
  std::deque<double> boundaryPointsNew = boundaryPoints;
  std::deque<double> interceptsNew = intercepts;
  std::deque<double> slopesNew = slopes;
  double lb, ub;
  
  // Calculate the initial lower bound
  lb = (-lambda + nodeWeight * nodeValue -
    interceptsNew.front()) / (slopesNew.front() + nodeWeight);
  
  // Update lower bound by iterating through boundaryPoints
  while (boundaryPointsNew.front() < lb) {
    interceptsNew.pop_front();
    slopesNew.pop_front();
    boundaryPointsNew.pop_front();
    lb = (-lambda + nodeWeight * nodeValue -
      interceptsNew.front()) / (slopesNew.front() + nodeWeight);
    if (slopesNew.size() == 1) break;
  }
  
  // Calculate the initial upper bound
  if (slopesNew.size() == 1) {
    ub = (lambda + nodeWeight * nodeValue -
      interceptsNew.front()) / (slopesNew.front() + nodeWeight);
  } else {
    ub = (lambda + nodeWeight * nodeValue -
      interceptsNew.back()) / (slopesNew.back() + nodeWeight);
    while (boundaryPointsNew.back() > ub) {
      interceptsNew.pop_back();
      slopesNew.pop_back();
      boundaryPointsNew.pop_back();
      ub = (lambda + nodeWeight * nodeValue -
        interceptsNew.back()) / (slopesNew.back() + nodeWeight);
      if (slopesNew.size() == 1) break;
    }
  }
  
  // Update boundary points
  boundaryPointsNew.push_front(lb);
  boundaryPointsNew.push_back(ub);
  
  // Update intercepts
  for (int i = 0; i < interceptsNew.size(); i++) {
    interceptsNew[i] -= nodeWeight * nodeValue;
  }
  interceptsNew.push_front(-lambda);
  interceptsNew.push_back(lambda);
  
  // Update slopes
  for (int i = 0; i < slopes.size(); i++) {
    slopesNew[i] += nodeWeight;
  }
  slopesNew.push_front(0);
  slopesNew.push_back(0);
  
}

