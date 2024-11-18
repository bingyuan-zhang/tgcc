#include <RcppArmadillo.h>
#include <tuple>

#include "computeTheta.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

struct Boundary {
  std::deque<double> boundaryPoints;
  std::deque<double> intercepts;
  std::deque<double> slopes;
};

std::tuple<Boundary, double> updateRootTheta(
    const double nodeValue,
    const double nodeWeight,
    const Boundary& boundary) {

  Boundary b = boundary;

  // Calculate the initial value of theta
  double theta = (nodeWeight * nodeValue - b.intercepts.front()) / (b.slopes.front() + nodeWeight);

  // Iterate through boundaryPoints to update theta
  while (b.boundaryPoints.front() < theta) {
    b.intercepts.pop_front();
    b.slopes.pop_front();
    b.boundaryPoints.pop_front();
    // Recalculate theta
    theta = (nodeWeight * nodeValue - b.intercepts.front()) / (b.slopes.front() + nodeWeight);
    // Break the loop if only one element remains in slopes
    if (b.slopes.size() == 1) break;
  }

  return std::make_tuple(b, theta);
}

std::tuple<Boundary, double, double>  updateBoundariesForLeaf(
    const double nodeValue,
    const double nodeWeight,
    const double lambda,
    const Boundary& boundary) {

  Boundary b = boundary;

  // Calculate the lower and upper bounds
  double lb = nodeValue - lambda / nodeWeight;
  double ub = nodeValue + lambda / nodeWeight;

  // Update boundary points (bp = c(l, u))
  b.boundaryPoints.push_front(lb);
  b.boundaryPoints.push_back(ub);

  // Update intercepts (i = c(-lam, -nw*xnode, lam))
  b.intercepts.push_front(-lambda);
  b.intercepts.push_back(-nodeWeight * nodeValue);
  b.intercepts.push_back(lambda);

  // Update slopes (s = c(0, nw, 0))
  b.slopes.push_front(0);
  b.slopes.push_back(nodeWeight);
  b.slopes.push_back(0);

  return std::make_tuple(b, lb, ub);
}

std::tuple<Boundary, double, double> updateBoundaries(
    const double nodeValue,
    const double nodeWeight,
    const double lambda,
    const Boundary& boundary) {

  Boundary b = boundary;
  double lb, ub;

  // Calculate the initial lower bound
  lb = (-lambda + nodeWeight * nodeValue -
    b.intercepts.front()) / (b.slopes.front() + nodeWeight);

  // Update lower bound by iterating through boundaryPoints
  while (b.boundaryPoints.front() < lb) {
    b.intercepts.pop_front();
    b.slopes.pop_front();
    b.boundaryPoints.pop_front();
    lb = (-lambda + nodeWeight * nodeValue -
      b.intercepts.front()) / (b.slopes.front() + nodeWeight);
    if (b.slopes.size() == 1) break;
  }

  // Calculate the initial upper bound
  if (b.slopes.size() == 1) {
    ub = (lambda + nodeWeight * nodeValue -
      b.intercepts.front()) / (b.slopes.front() + nodeWeight);
  } else {
    ub = (lambda + nodeWeight * nodeValue -
      b.intercepts.back()) / (b.slopes.back() + nodeWeight);
    while (b.boundaryPoints.back() > ub) {
      b.intercepts.pop_back();
      b.slopes.pop_back();
      b.boundaryPoints.pop_back();
      ub = (lambda + nodeWeight * nodeValue -
        b.intercepts.back()) / (b.slopes.back() + nodeWeight);
      if (b.slopes.size() == 1) break;
    }
  }

  // Update boundary points
  b.boundaryPoints.push_front(lb);
  b.boundaryPoints.push_back(ub);

  // Update intercepts
  for (int i = 0; i < b.intercepts.size(); i++) {
    b.intercepts[i] -= nodeWeight * nodeValue;
  }
  b.intercepts.push_front(-lambda);
  b.intercepts.push_back(lambda);

  // Update slopes
  for (int i = 0; i < b.slopes.size(); i++) {
    b.slopes[i] += nodeWeight;
  }
  b.slopes.push_front(0);
  b.slopes.push_back(0);

  return std::make_tuple(b, lb, ub);
}

Boundary mergeQueues(
    const std::vector<double>& children,
    const std::vector<Boundary>& Boundaries) {

  int numChildren = children.size();
  std::vector<double> bpLengths(numChildren);
  std::deque<double> mergedBoundaryPoints;

  // Concatenate all boundary points from children
  for (int i = 0; i < numChildren; i++) {
    int child = children[i];
    bpLengths[i] = Boundaries[child].boundaryPoints.size();
    mergedBoundaryPoints.insert(
      mergedBoundaryPoints.end(),
      Boundaries[child].boundaryPoints.begin(),
      Boundaries[child].boundaryPoints.end());
  }

  int totalLength = mergedBoundaryPoints.size();

  // Obtain the rank of merged boundary points
  std::vector<std::pair<double, int>> bpSort(totalLength);
  std::vector<double> ranks(totalLength);
  for (int i = 0; i < totalLength; i++) {
    bpSort[i] = std::make_pair(mergedBoundaryPoints[i], i);
  }
  std::sort(bpSort.begin(), bpSort.end());
  for (int i = 0; i < totalLength; i++) {
    ranks[bpSort[i].second] = i;
  }
  std::sort(mergedBoundaryPoints.begin(), mergedBoundaryPoints.end());

  // Merge queues
  int start1 = 0;
  std::vector<std::vector<double>> interceptsToSum(numChildren, std::vector<double>(totalLength + 1));
  std::vector<std::vector<double>> slopesToSum(numChildren, std::vector<double>(totalLength + 1));
  for (int i = 0; i < numChildren; i++) {
    int child = children[i];
    int start2 = 0;
    for (int j = 0; j < bpLengths[i]; j++) {
      std::fill(interceptsToSum[i].begin() + start2,
        interceptsToSum[i].begin() + ranks[start1 + j] + 1, Boundaries[child].intercepts[j]);

      std::fill(slopesToSum[i].begin() + start2,
        slopesToSum[i].begin() + ranks[start1 + j] + 1, Boundaries[child].slopes[j]);

      start2 = ranks[start1 + j] + 1;
    }
    std::fill(interceptsToSum[i].begin() + start2,
      interceptsToSum[i].end(), Boundaries[child].intercepts[bpLengths[i]]);

    std::fill(slopesToSum[i].begin() + start2,
      slopesToSum[i].end(), Boundaries[child].slopes[bpLengths[i]]);

    start1 += bpLengths[i];
  }

  std::deque<double> newIntercepts(totalLength + 1, 0.0);
  std::deque<double> newSlopes(totalLength + 1, 0.0);

  // Update the intercepts and slopes
  for (int j = 0; j < totalLength + 1; j++) {
    for (int i = 0; i < numChildren; i++) {
      newIntercepts[j] += interceptsToSum[i][j];
      newSlopes[j] += slopesToSum[i][j];
    }
  }

  Boundary b;
  b.boundaryPoints = mergedBoundaryPoints;
  b.intercepts = newIntercepts;
  b.slopes = newSlopes;

  return b;
}

// [[Rcpp::export]]
std::vector<double> computeTheta(
    const std::vector<double>& nodeValues,
    const double lambda,
    const std::vector<double>& vertices,
    const std::vector<double>& nodeTypes,
    const std::vector<double>& parents,
    const std::vector<double>& nodeWeights,
    const std::vector<double>& edgeWeights,
    const std::vector<std::vector<double>>& childrenList) {

  int numNodes = nodeValues.size();
  std::vector<Boundary> nodeBoundaries(numNodes);
  std::vector<double> lowerBounds(numNodes);
  std::vector<double> upperBounds(numNodes);

  std::vector<double> theta(numNodes, 0);
  Boundary updatedBoundary;
  double lb, ub, eta;

  // Forward pass
  for (int i = 0; i < numNodes; i++) {

    int vertex = vertices[i];
    double nodeType = nodeTypes[vertex];
    double nodeWeight = nodeWeights[vertex];
    double edgeWeight = edgeWeights[vertex];
    double nodeValue = nodeValues[vertex];

    if (nodeType == 0) { // leaves

      double wlambda = edgeWeight * lambda;
      std::tie(updatedBoundary, lb, ub) = updateBoundariesForLeaf(
        nodeValue, nodeWeight, wlambda, nodeBoundaries[vertex]
      );
      nodeBoundaries[vertex] = updatedBoundary;
      lowerBounds[i] = lb;
      upperBounds[i] = ub;

    }
    else if (nodeType == 1) { // single child

      int child = childrenList[vertex][0];
      nodeBoundaries[vertex] = nodeBoundaries[child];

      double wlambda = edgeWeight * lambda;
      std::tie(updatedBoundary, lb, ub) = updateBoundaries(
        nodeValue, nodeWeight, wlambda, nodeBoundaries[vertex]
      );
      nodeBoundaries[vertex] = updatedBoundary;
      lowerBounds[i] = lb;
      upperBounds[i] = ub;

    }
    else if (nodeType == 2) { // multiple children
      std::vector<double> children = childrenList[vertex];
      updatedBoundary = mergeQueues(children, nodeBoundaries);

      double wlambda = edgeWeight * lambda;
      std::tie(updatedBoundary, lb, ub) = updateBoundaries(
        nodeValue, nodeWeight, wlambda, updatedBoundary
      );
      nodeBoundaries[vertex] = updatedBoundary;
      lowerBounds[i] = lb;
      upperBounds[i] = ub;

    }
    else { // root
      if (nodeType == 3) { // root with single child
        int child = static_cast<int>(childrenList[vertex][0]);
        std::tie(updatedBoundary, eta) = updateRootTheta(
          nodeValue, nodeWeight, nodeBoundaries[child]
        );
        nodeBoundaries[child] = updatedBoundary;
        theta[vertex] = eta;
      }
      else { // root with multiple children
        std::vector<double> children = childrenList[vertex];

        updatedBoundary = mergeQueues(children, nodeBoundaries);

        std::tie(updatedBoundary, eta) = updateRootTheta(
          nodeValue, nodeWeight, updatedBoundary
        );
        nodeBoundaries[vertex] = updatedBoundary;
        theta[vertex] = eta;
      }
    }
  }

  // Backward pass
  for (int i = numNodes - 2; i >= 0; i--) {
    int vertex = static_cast<int>(vertices[i]);
    int parent = static_cast<int>(parents[vertex]);
    if (theta[parent] < lowerBounds[i]) {
      theta[vertex] = lowerBounds[i];
    }
    else if (theta[parent] > upperBounds[i]) {
      theta[vertex] = upperBounds[i];
    }
    else {
      theta[vertex] = theta[parent];
    }
  }

  return theta;
}
