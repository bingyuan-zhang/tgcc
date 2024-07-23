#include <RcppArmadillo.h>
#include "computeTheta.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

void updateBoundariesForLeaf(
    double nodeValue,
    double nodeWeight,
    double lambda,
    std::deque<double>& boundaryPoints,
    std::deque<double>& intercepts,
    std::deque<double>& slopes,
    double& lowerBound,
    double& upperBound) {

  // Calculate the lower and upper bounds
  lowerBound = nodeValue - lambda / nodeWeight;
  upperBound = nodeValue + lambda / nodeWeight;

  // Update boundary points (bp = c(l, u))
  boundaryPoints.push_front(lowerBound);
  boundaryPoints.push_back(upperBound);

  // Update intercepts (i = c(-lam, -nw*xnode, lam))
  intercepts.push_front(-lambda);
  intercepts.push_back(-nodeWeight * nodeValue);
  intercepts.push_back(lambda);

  // Update slopes (s = c(0, nw, 0))
  slopes.push_front(0);
  slopes.push_back(nodeWeight);
  slopes.push_back(0);
}


void updateRootTheta(
    double nodeValue,
    double nodeWeight,
    std::deque<double>& boundaryPoints,
    std::deque<double>& intercepts,
    std::deque<double>& slopes,
    double& theta) {

  // Calculate the initial value of theta
  theta = (nodeWeight * nodeValue - intercepts.front()) / (slopes.front() + nodeWeight);

  // Iterate through boundaryPoints to update theta
  while (boundaryPoints.front() < theta) {
    intercepts.pop_front();
    slopes.pop_front();
    boundaryPoints.pop_front();
    // Recalculate theta
    theta = (nodeWeight * nodeValue - intercepts.front()) / (slopes.front() + nodeWeight);
    // Break the loop if only one element remains in slopes
    if (slopes.size() == 1) break;
  }

}

void updateBoundaries(
    double nodeValue,
    double nodeWeight,
    double lambda,
    std::deque<double>& boundaryPoints,
    std::deque<double>& intercepts,
    std::deque<double>& slopes,
    double& lowerBound,
    double& upperBound) {

  // Calculate the initial lower bound
  lowerBound = (-lambda + nodeWeight * nodeValue -
    intercepts.front()) / (slopes.front() + nodeWeight);

  // Update lower bound by iterating through boundaryPoints
  while (boundaryPoints.front() < lowerBound) {
    intercepts.pop_front();
    slopes.pop_front();
    boundaryPoints.pop_front();
    lowerBound = (-lambda + nodeWeight * nodeValue -
      intercepts.front()) / (slopes.front() + nodeWeight);
    if (slopes.size() == 1) break;
  }

  // Calculate the initial upper bound
  if (slopes.size() == 1) {
    upperBound = (lambda + nodeWeight * nodeValue -
      intercepts.front()) / (slopes.front() + nodeWeight);
  } else {
    upperBound = (lambda + nodeWeight * nodeValue -
      intercepts.back()) / (slopes.back() + nodeWeight);
    while (boundaryPoints.back() > upperBound) {
      intercepts.pop_back();
      slopes.pop_back();
      boundaryPoints.pop_back();
      upperBound = (lambda + nodeWeight * nodeValue -
        intercepts.back()) / (slopes.back() + nodeWeight);
      if (slopes.size() == 1) break;
    }
  }

  // Update boundary points
  boundaryPoints.push_front(lowerBound);
  boundaryPoints.push_back(upperBound);

  // Update intercepts
  for (int i = 0; i < intercepts.size(); i++) {
    intercepts[i] -= nodeWeight * nodeValue;
  }
  intercepts.push_front(-lambda);
  intercepts.push_back(lambda);

  // Update slopes
  for (int i = 0; i < slopes.size(); i++) {
    slopes[i] += nodeWeight;
  }
  slopes.push_front(0);
  slopes.push_back(0);
}

void mergeQueues(
    const std::vector<double>& children,
    int parentIndex,
    std::vector<std::deque<double>>& boundaryPoints,
    std::vector<std::deque<double>>& intercepts,
    std::vector<std::deque<double>>& slopes) {

  int numChildren = children.size();
  int child;
  std::vector<double> bpLengths(numChildren);
  std::deque<double> mergedBoundaryPoints;

  // Concatenate all boundary points from children
  for (int i = 0; i < numChildren; i++) {
    child = children[i];
    bpLengths[i] = boundaryPoints[child].size();
    mergedBoundaryPoints.insert(mergedBoundaryPoints.end(),
      boundaryPoints[child].begin(), boundaryPoints[child].end());
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
  int start2 = 0;
  std::vector<std::vector<double>> interceptsToSum(numChildren, std::vector<double>(totalLength + 1));
  std::vector<std::vector<double>> slopesToSum(numChildren, std::vector<double>(totalLength + 1));
  for (int i = 0; i < numChildren; i++) {
    child = children[i];
    start2 = 0;

    for (int j = 0; j < bpLengths[i]; j++) {

      std::fill(interceptsToSum[i].begin() + start2,
        interceptsToSum[i].begin() + ranks[start1 + j] + 1, intercepts[child][j]);

      std::fill(slopesToSum[i].begin() + start2,
        slopesToSum[i].begin() + ranks[start1 + j] + 1, slopes[child][j]);

      start2 = ranks[start1 + j] + 1;
    }

    std::fill(interceptsToSum[i].begin() + start2,
      interceptsToSum[i].end(), intercepts[child][bpLengths[i]]);

    std::fill(slopesToSum[i].begin() + start2,
      slopesToSum[i].end(), slopes[child][bpLengths[i]]);

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

  // Update the boundary points, intercepts, and slopes
  boundaryPoints[parentIndex] = mergedBoundaryPoints;
  intercepts[parentIndex] = newIntercepts;
  slopes[parentIndex] = newSlopes;

}

// [[Rcpp::export]]
std::vector<double> computeTheta(
    const std::vector<double>& nodeValues,
    double lambda,
    const std::vector<double>& vertices,
    const std::vector<double>& nodeTypes,
    const std::vector<double>& parents,
    const std::vector<double>& nodeWeights,
    const std::vector<double>& edgeWeights,
    const std::vector<std::vector<double>>& childrenList) {

  int numNodes = nodeValues.size();
  std::vector<std::deque<double>> boundaryPoints(numNodes);
  std::vector<std::deque<double>> intercepts(numNodes);
  std::vector<std::deque<double>> slopes(numNodes);
  std::vector<double> lowerBounds(numNodes);
  std::vector<double> upperBounds(numNodes);

  int vertex;
  int nodeType;
  int parent;
  int child;

  double nodeWeight;
  double edgeWeight;
  double nodeValue;

  std::vector<double> children;
  std::vector<double> theta(numNodes, 0);

  // Forward pass
  for (int i = 0; i < numNodes; i++) {
    vertex = vertices[i];
    nodeType = nodeTypes[vertex];
    nodeWeight = nodeWeights[vertex];
    edgeWeight = edgeWeights[vertex];
    nodeValue = nodeValues[vertex];

    if (nodeType == 0) {

      updateBoundariesForLeaf(
        nodeValue, nodeWeight, edgeWeight * lambda,
        boundaryPoints[vertex], intercepts[vertex], slopes[vertex],
        lowerBounds[i], upperBounds[i]);

    } else if (nodeType == 1) {

      child = childrenList[vertex][0];
      boundaryPoints[vertex] = boundaryPoints[child];
      intercepts[vertex] = intercepts[child];
      slopes[vertex] = slopes[child];

      updateBoundaries(
        nodeValue, nodeWeight, edgeWeight * lambda,
        boundaryPoints[vertex], intercepts[vertex], slopes[vertex],
        lowerBounds[i], upperBounds[i]);

    } else if (nodeType == 2) {

      children = childrenList[vertex];
      mergeQueues(children, vertex, boundaryPoints, intercepts, slopes);

      updateBoundaries(
        nodeValue, nodeWeight, edgeWeight * lambda,
        boundaryPoints[vertex], intercepts[vertex], slopes[vertex],
        lowerBounds[i], upperBounds[i]);

    } else {
      if (nodeType == 3) {

        child = childrenList[vertex][0];

        updateRootTheta(nodeValue, nodeWeight,
          boundaryPoints[child], intercepts[child], slopes[child],
          theta[vertex]);

      } else {

        children = childrenList[vertex];
        mergeQueues(children, vertex, boundaryPoints, intercepts, slopes);

        updateRootTheta(nodeValue, nodeWeight,
          boundaryPoints[vertex], intercepts[vertex], slopes[vertex],
          theta[vertex]);

      }
    }
  }

  // Backward pass
  for (int i = numNodes - 2; i >= 0; i--) {
    vertex = vertices[i];
    parent = parents[vertex];
    if (theta[parent] < lowerBounds[i]) {
      theta[vertex] = lowerBounds[i];
    } else if (theta[parent] > upperBounds[i]) {
      theta[vertex] = upperBounds[i];
    } else {
      theta[vertex] = theta[parent];
    }
  }

  return theta;
}
