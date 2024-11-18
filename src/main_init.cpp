#include <RcppArmadillo.h>
#include <vector>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

struct edgeSet {
  std::vector<std::vector<double>> adjList;
  std::vector<std::vector<double>> distanceList;

  edgeSet(int numNodes) {
    adjList.resize(numNodes);
    distanceList.resize(numNodes);
  }
};

struct treeVar {
  std::vector<double> vertices;
  std::vector<double> depths;
  std::vector<double> types;
  std::vector<double> edgeWeights;
  std::vector<double> nodeWeights;
  std::vector<double> parents;
  std::vector<double> partitions;
  std::vector<std::vector<double>> childrenList;

  treeVar(int numNodes) {
    vertices.resize(numNodes, 0.0);
    depths.resize(numNodes, 0.0);
    types.resize(numNodes, 0.0);
    edgeWeights.resize(numNodes, 0.0);
    nodeWeights.resize(numNodes, 1.0);
    parents.resize(numNodes, 0.0);
    partitions.resize(numNodes, 1.0);
    childrenList.resize(numNodes, std::vector<double>());
  }
};

edgeSet findConnectedEdges(
    const std::vector<double>& sourceNodes,
    const std::vector<double>& targetNodes,
    const std::vector<double>& distances) {

  int edgeCount = sourceNodes.size();
  int nodeCount = edgeCount + 1; // For a tree，number of node = number of edge + 1

  edgeSet edges(nodeCount);

  for (int i = 0; i < edgeCount; i++) {
    int source = static_cast<int>(sourceNodes[i]);
    int target = static_cast<int>(targetNodes[i]);
    double distance = distances[i];

    edges.adjList[source].push_back(target);
    edges.adjList[target].push_back(source);
    edges.distanceList[source].push_back(distance);
    edges.distanceList[target].push_back(distance);
  }

  return edges;
}

treeVar initializeParams(
    const edgeSet& edges,
    double gamma) {

  int nodeCount = edges.adjList.size();
  treeVar tree(nodeCount);

  int seenCount = 1;
  int currentDepth = 0;
  std::vector<double> currentDepthVertices = edges.adjList[0];
  std::vector<double> nextDepthVertices;

  // Initialize root parameters
  tree.childrenList[0] = edges.adjList[0];
  tree.types[0] = (edges.adjList[0].size() == 1) ? 3 : 4;

  // Initialize root's children
  for (size_t i = 0; i < edges.adjList[0].size(); ++i) {
    int child = edges.adjList[0][i];
    tree.parents[child] = 0;
    tree.edgeWeights[child] = std::exp(-std::pow(edges.distanceList[0][i], 2) / gamma);
  }

  // BFS to initialize other nodes
  while (seenCount < nodeCount) {
    currentDepth += 1;
    for (size_t i = 0; i < currentDepthVertices.size(); ++i) {
      int vertex = currentDepthVertices[i];
      tree.vertices[seenCount] = vertex;
      tree.depths[seenCount] = currentDepth;

      for (size_t j = 0; j < edges.adjList[vertex].size(); ++j) {
        int neighbor = edges.adjList[vertex][j];
        if (neighbor != tree.parents[vertex]) {
          // Find children of the vertex
          tree.childrenList[vertex].push_back(neighbor);
          // Set the parent of the neighbor
          tree.parents[neighbor] = vertex;
          // Set edge weights between vertex and its children
          tree.edgeWeights[neighbor] = std::exp(-std::pow(edges.distanceList[vertex][j], 2) / gamma);
        }
      }

      // Update types based on the number of children
      if (tree.childrenList[vertex].empty()) {
        tree.types[vertex] = 0; // Leaf node
      } else if (tree.childrenList[vertex].size() == 1) {
        tree.types[vertex] = 1; // Node with one child
      } else {
        tree.types[vertex] = 2; // Node with multiple children
      }

      seenCount++;
      nextDepthVertices.insert(nextDepthVertices.end(), tree.childrenList[vertex].begin(), tree.childrenList[vertex].end());
    }
    currentDepthVertices = std::move(nextDepthVertices);
    nextDepthVertices.clear();
  }

  // Reverse the vertices and depths
  std::reverse(tree.vertices.begin(), tree.vertices.end());
  std::reverse(tree.depths.begin(), tree.depths.end());

  // Update partitions based on children
  for (int i = 0; i < nodeCount; ++i) {
    int vertex = tree.vertices[i];
    if (tree.types[vertex] != 0) {
      if (tree.types[vertex] == 2) {
        for (const auto& child : tree.childrenList[vertex]) {
          tree.partitions[vertex] += tree.partitions[child];
        }
      } else if (tree.types[vertex] != 3) {
          tree.partitions[vertex] += tree.partitions[tree.childrenList[vertex][0]];
      }
    }
  }

  // Adjust partitions
  for (int i = 0; i < nodeCount; ++i) {
    int maxPartition = nodeCount - tree.partitions[i];
    if (tree.partitions[i] > maxPartition) {
      tree.partitions[i] = maxPartition;
    }
  }

  return tree;
}

// [[Rcpp::export]]
List initPrepare (
    const std::vector<double>& from,
    const std::vector<double>& to,
    const std::vector<double>& dist,
    double gamma) {

  edgeSet edges = findConnectedEdges(from, to, dist);
  treeVar tree = initializeParams(edges, gamma);

  List res = List::create(
    Named("Vertices")    = wrap(tree.vertices),
    Named("Depths")      = wrap(tree.depths),
    Named("Types")       = wrap(tree.types),
    Named("EdgeWeights") = wrap(tree.edgeWeights),
    Named("NodeWeights") = wrap(tree.nodeWeights),
    Named("Parents")     = wrap(tree.parents),
    Named("Partitions")  = wrap(tree.partitions),
    Named("Children")    = wrap(tree.childrenList)
  );

  return res;
}
