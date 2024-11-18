#include <RcppArmadillo.h>
#include "computeTheta.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// Calculate the loss function for FISTA
double lossFuncFista(
    const arma::colvec &y,
    const arma::colvec &x,
    double lambda) {

  double term1 = 0.5 * std::pow(arma::norm(y - x, 2), 2);
  double term2 = lambda * arma::norm(x, 2);

  return term1 + term2;
}

// Perform FISTA (Fast Iterative Shrinkage-Thresholding Algorithm) descent
arma::colvec fistaDescent(
    const arma::colvec & y,
    double lambda,
    int maxIter) {

  int n = y.n_elem;
  arma::colvec xOld(n);
  arma::colvec xNew(n);
  arma::colvec vOld(n);
  arma::colvec vNew(n);
  arma::colvec uTmp(n);
  arma::colvec uu(n);
  arma::colvec u(n);
  arma::colvec loss(maxIter);

  double lossOld = arma::datum::inf;
  double t;
  double threshold;
  double lossNew;
  bool descent;
  int noDescentTol = 10; // Constant for no descent tolerance
  int noDescent = 0;

  for (int iter = 0; iter < maxIter; ++iter) {

    // FISTA step
    t = 2.0 / (iter + 2);

    uu = (1 - t) * xOld + t * vOld;
    uTmp = uu - 0.5 * (uu - y);
    threshold = std::max(0.0, 1.0 - lambda / arma::norm(uTmp,2) );
    u = uTmp * threshold;

    lossNew = lossFuncFista(y, u, lambda);
    descent = lossNew < lossOld;

    if (descent) {
      xNew = u;
      lossOld = lossNew;
      noDescent = 0;
    } else {
      xNew = xOld;
      noDescent++;
    }

    vNew = xOld + 1 / t * (u - xOld);

    if (noDescent >= noDescentTol) break;

    xOld = xNew;
    vOld = vNew;
    loss[iter] = lossOld;
  }

  return xOld;
}

// Calculate the sparse clustering loss function
double SparseClusteringLoss(
    const arma::mat& data,
    const arma::mat& U,
    double lambda,
    double gamma) {

  double term1 = 0.5 * std::pow(arma::norm(data - U, "fro"), 2);

  int numRows = data.n_rows;
  int numCols = data.n_cols;
  double term2 = 0;
  double term3 = 0;

  // term 2
  for (int i = 0; i < numCols; ++i) {
    for (int j = i + 1; j < numCols; ++j) {
      term2 += arma::sum(arma::abs(U.col(i) - U.col(j)));
    }
  }
  // term 3
  for (int i = 0; i < numRows; ++i) {
    term3 += arma::norm(U.row(i), 2);
  }

  return (term1 + lambda * term2 + gamma * term3);
}

// [[Rcpp::export]]
arma::mat updateUSC(
    const arma::mat& dataMatrix,
    const arma::mat& initialMatrix,
    double lambda, double gamma, int maxIterations,
    double precisionThreshold,
    const std::vector<double>& vertices,
    const std::vector<double>& nodeTypes,
    const std::vector<double>& parents,
    const std::vector<double>& nodeWeights,
    const std::vector<double>& edgeWeights,
    const std::vector<std::vector<double>>& childrenList) {

  int numRows = dataMatrix.n_rows;
  int numCols = dataMatrix.n_cols;

  arma::mat updatedMatrix = initialMatrix;
  arma::mat P = arma::zeros(numRows, numCols);
  arma::mat Q = arma::zeros(numRows, numCols);
  arma::mat Y = arma::zeros(numRows, numCols);

  arma::vec diffIter = arma::zeros(maxIterations);
  arma::vec lossIter = arma::zeros(maxIterations);

  double diff = arma::datum::inf;
  double previousLoss = SparseClusteringLoss(initialMatrix, initialMatrix, lambda, gamma);
  double currentLoss;

  std::vector<double> tempColumn(numRows);
  std::vector<double> tempColumnResult(numRows);

  for (int iter = 0; iter < maxIterations; ++iter) {

    // Update Y
    for (int j = 0; j < numCols; ++j) {
      Y.col(j) = fistaDescent(updatedMatrix.col(j) + P.col(j), gamma, maxIterations);

      P.col(j) = updatedMatrix.col(j) + P.col(j) - Y.col(j);

      tempColumn = arma::conv_to<std::vector<double>>::from(Y.col(j) + Q.col(j));
      tempColumnResult = computeTheta(tempColumn, lambda, vertices, nodeTypes, parents, nodeWeights, edgeWeights, childrenList);
      updatedMatrix.col(j) = arma::conv_to<arma::colvec>::from(tempColumnResult);

      Q.col(j) = Y.col(j) + Q.col(j) - updatedMatrix.col(j);
    }

    // Calculate current loss and check convergence
    currentLoss = SparseClusteringLoss(initialMatrix, updatedMatrix, lambda, gamma);

    diff = std::abs(currentLoss - previousLoss) / previousLoss;
    lossIter(iter) = currentLoss;
    diffIter(iter) = diff;
    previousLoss = currentLoss;

    if (diff < precisionThreshold) {
      break;
    }

  }

  return updatedMatrix;
}

// [[Rcpp::export]]
double BiClusteringLoss(const arma::mat& data,
                        const arma::mat& U,
                        double lambda,
                        double gamma) {

  double term1 = 0.5 * std::pow(arma::norm(data - U, "fro"), 2);
  int numRows = data.n_rows;
  int numCols = data.n_cols;
  double term2 = 0;
  double term3 = 0;

  // term 2
  for (int i = 0; i < numCols; ++i) {
    for (int j = i + 1; j < numCols; ++j) {
      term2 += arma::sum(arma::abs(U.col(i) - U.col(j)));
    }
  }
  // term 3
  for (int i = 0; i < numRows; ++i) {
    for (int j = i + 1; j < numRows; ++j) {
      term3 += arma::sum(arma::abs(U.row(i) - U.row(j)));
    }
  }

  return (term1 + lambda * term2 + gamma * term3);
}

// [[Rcpp::export]]
arma::mat updateUBC(
    const arma::mat& dataMatrix, const arma::mat& initialMatrix,
    double lambda, double gamma, int maxIterations, double precisionThreshold,
    const std::vector<double>& verticesSamples,
    const std::vector<double>& nodeTypesSamples,
    const std::vector<double>& parentsSamples,
    const std::vector<double>& nodeWeightsSamples,
    const std::vector<double>& edgeWeightsSamples,
    const std::vector<std::vector<double>>& childrenListSamples,
    const std::vector<double>& verticesFeatures,
    const std::vector<double>& nodeTypesFeatures,
    const std::vector<double>& parentsFeatures,
    const std::vector<double>& nodeWeightsFeatures,
    const std::vector<double>& edgeWeightsFeatures,
    const std::vector<std::vector<double>>& childrenListFeatures) {

  int numRows = dataMatrix.n_rows;
  int numCols = dataMatrix.n_cols;

  arma::mat updatedMatrix = initialMatrix;
  arma::mat P = arma::zeros(numRows, numCols);
  arma::mat Q = arma::zeros(numRows, numCols);
  arma::mat Y = arma::zeros(numRows, numCols);

  arma::vec diffIter = arma::zeros(maxIterations);
  arma::vec lossIter = arma::zeros(maxIterations);

  double diff = arma::datum::inf;
  double previousLoss = BiClusteringLoss(initialMatrix, initialMatrix, lambda, gamma);
  double currentLoss;

  std::vector<double> tempCol(numRows);
  std::vector<double> tempColResult(numRows);
  std::vector<double> tempRow(numCols);
  std::vector<double> tempRowResult(numCols);

  for (int iter = 0; iter < maxIterations; ++iter) {

    // Update Y
    for (int j = 0; j < numCols; ++j) {
      tempCol = arma::conv_to<std::vector<double>>::from(updatedMatrix.col(j) + P.col(j));
      tempColResult = computeTheta(tempCol, lambda, verticesSamples, nodeTypesSamples, parentsSamples, nodeWeightsSamples, edgeWeightsSamples, childrenListSamples);
      Y.col(j) = arma::conv_to<arma::colvec>::from(tempColResult);
    }
    // Update P
    P = updatedMatrix + P - Y;

    // Update U
    for (int j = 0; j < numRows; ++j) {
      tempRow = arma::conv_to<std::vector<double>>::from(Y.row(j) + Q.row(j));
      tempRowResult = computeTheta(tempRow, gamma, verticesFeatures, nodeTypesFeatures, parentsFeatures, nodeWeightsFeatures, edgeWeightsFeatures, childrenListFeatures);
      updatedMatrix.row(j) = arma::conv_to<arma::rowvec>::from(tempRowResult);
    }

    // Update Q
    Q = Y + Q - updatedMatrix;

    // Calculate current loss and check convergence
    currentLoss = BiClusteringLoss(initialMatrix, updatedMatrix, lambda, gamma);
    diff = std::abs((currentLoss - previousLoss) / previousLoss);

    lossIter(iter) = currentLoss;
    diffIter(iter) = diff;
    previousLoss = currentLoss;

    if (diff < precisionThreshold) break;
  }
  return updatedMatrix;
}
