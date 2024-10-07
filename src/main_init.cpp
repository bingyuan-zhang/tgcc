#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

void findConnectedEdges(
    const std::vector<double>& sourceNodes,
    const std::vector<double>& targetNodes,
    const std::vector<double>& distances,
    std::vector<std::vector<double>>& adjList,
    std::vector<std::vector<double>>& distanceList) {

  int edgeCount = sourceNodes.size();

  for (int i = 0; i < edgeCount; i++) {
    int source = sourceNodes[i];
    int target = targetNodes[i];
    double distance = distances[i];

    adjList[source].push_back(target);
    adjList[target].push_back(source);
    distanceList[source].push_back(distance);
    distanceList[target].push_back(distance);
  }
}

void initializeParams(
    const std::vector<std::vector<double>>& adjList,
    const std::vector<std::vector<double>>& distanceList,
    std::vector<std::vector<double>>& childrenList,
    std::vector<double>& vertices,
    std::vector<double>& depths,
    std::vector<double>& types,
    std::vector<double>& edgeWeights,
    std::vector<double>& parents,
    std::vector<double>& partitions,
    double gamma) {

  int nodeCount = adjList.size();
  int seenCount = 1;
  int currentDepth = 0;
  std::vector<double> currentDepthVertices = adjList[0];
  std::vector<double> nextDepthVertices;

  // Initialize root parameters
  vertices[0] = 0;
  depths[0] = 0;
  parents[0] = 0;
  edgeWeights[0] = 0;
  childrenList[0] = adjList[0];
  types[0] = (adjList[0].size() == 1) ? 3 : 4;

  // Initialize root's children
  for (int i = 0; i < adjList[0].size(); ++i) {
    parents[adjList[0][i]] = 0;
    edgeWeights[adjList[0][i]] = std::exp(-std::pow(distanceList[0][i], 2) / gamma);
  }

  // BFS to initialize other nodes
  while (seenCount < nodeCount) {
    currentDepth += 1;
    for (int i = 0; i < currentDepthVertices.size(); ++i) {
      int vertex = currentDepthVertices[i];
      vertices[seenCount] = vertex;
      depths[seenCount] = currentDepth;
      for (int j = 0; j < adjList[vertex].size(); ++j) {
        int neighbor = adjList[vertex][j];
        if (neighbor != parents[vertex]) {
          // Find children of the vertex
          childrenList[vertex].push_back(neighbor);
          // Set the parent of the neighbor
          parents[neighbor] = vertex;
          // Set edge weights between vertex and its children
          edgeWeights[neighbor] = std::exp(-std::pow(distanceList[vertex][j], 2) / gamma);
        }
      }
      // Update types based on the number of children
      if (childrenList[vertex].empty()) {
        types[vertex] = 0; // Leaf node
      } else if (childrenList[vertex].size() == 1) {
        types[vertex] = 1; // Node with one child
      } else {
        types[vertex] = 2; // Node with more than one child
      }
      seenCount++;
      nextDepthVertices.insert(nextDepthVertices.end(), childrenList[vertex].begin(), childrenList[vertex].end());
    }
    currentDepthVertices = std::move(nextDepthVertices);
    nextDepthVertices.clear();
  }

  // Reverse the vertices and depths
  std::reverse(vertices.begin(), vertices.end());
  std::reverse(depths.begin(), depths.end());

  // Update partitions based on children
  for (int i = 0; i < nodeCount; ++i) {
    int vertex = vertices[i];
    if (types[vertex] != 0) {
      if (types[vertex] == 2) {
        for (const auto& child : childrenList[vertex]) {
          partitions[vertex] += partitions[child];
        }
      } else if (types[vertex] != 3) { // Type 1 or 4
        partitions[vertex] += partitions[childrenList[vertex][0]];
      }
    }
  }

  // Adjust partitions
  for (int i = 0; i < nodeCount; ++i) {
    int maxPartition = nodeCount - partitions[i];
    if (partitions[i] > maxPartition) {
      partitions[i] = maxPartition;
    }
  }
}

//[[Rcpp::export]]
Rcpp::List initPrepare (
    const std::vector<double>& from,
    const std::vector<double>& to,
    const std::vector<double>& dist,
    double gamma) {

  int n = from.size() + 1;

  std::vector<std::vector<double>> adjList(n);
  std::vector<std::vector<double>> distanceList(n);
  std::vector<std::vector<double>> childrenList(n);

  // Fill the adjacency list and distance list with edges
  findConnectedEdges(from, to, dist, adjList, distanceList);

  std::vector<double> vertices(n);
  std::vector<double> depths(n);
  std::vector<double> types(n);
  std::vector<double> edgeWeights(n);
  std::vector<double> nodeWeights(n, 1.0);
  std::vector<double> parents(n);
  std::vector<double> partitions(n, 1.0);

  // Initialize the node properties using BFS
  initializeParams(adjList, distanceList, childrenList, vertices,
    depths, types, edgeWeights, parents, partitions, gamma);

  Rcpp::List res = List::create(
    Named("Vertices")   = wrap(vertices),
    Named("Depths")     = wrap(depths),
    Named("Types")      = wrap(types),
    Named("EdgeWeights")= wrap(edgeWeights),
    Named("NodeWeights")= wrap(nodeWeights),
    Named("Parents")    = wrap(parents),
    Named("Partitions") = wrap(partitions),
    Named("Children")   = wrap(childrenList)
  );

  return (res);
}
