#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

struct TreeParams {
  std::vector<double> vertices;
  std::vector<double> types;
  std::vector<double> parents;
  std::vector<double> edgeWeights;
  std::vector<double> nodeWeights;
  std::vector<std::vector<double>> childrenList;
  std::vector<double> pointer;
  std::vector<double> rowVertices;

  TreeParams(int numNodes) {
    vertices.resize(numNodes, 0.0);
    parents.resize(numNodes, 0.0);
    types.resize(numNodes, 0.0);
    edgeWeights.resize(numNodes, 0.0);
    nodeWeights.resize(numNodes, 0.0);
    childrenList.resize(numNodes, std::vector<double>());
    pointer.resize(numNodes, 0.0);
    rowVertices.resize(numNodes, 0.0);
  }
};

std::tuple<arma::mat, TreeParams> update(
    const arma::mat& Theta,
    const arma::mat& input,
    const TreeParams& params) {

  TreeParams newParams = params;

  int n = params.vertices.size();

  // new Params
  std::vector<double> newRowVertices;
  std::vector<double> newParents;
  std::vector<double> newTypes;
  std::vector<double> newNodeWeights;
  std::vector<double> newEdgeWeights;
  std::vector<double> newDepths;
  std::vector<double> siblings;
  std::vector<double> newChildren;
  std::vector<double> newPointer(n, 0);
  std::vector<double> parentDepths(n, 0.0);

  arma::mat newInput;

  int root = newParams.vertices[n - 1]; // Get the root node
  newRowVertices.push_back(root);
  newParents.push_back(newParams.parents[0]);
  newTypes.push_back(newParams.types[0]);
  newNodeWeights.push_back(newParams.nodeWeights[0]);
  newEdgeWeights.push_back(0);
  newDepths.push_back(0);
  parentDepths[0] = 0;
  newPointer[0] = 0;
  newInput.insert_rows(0, input.row(0));

  for (int i = n - 2; i >= 0; i--) {

    int node = newParams.vertices[i];
    int parent = newParams.parents[node];
    int parentIndex = newPointer[parent];

    // Calculate the distance between the node and its parent
    double distance = arma::norm(Theta.row(node) - Theta.row(parent), 2);

    if (distance != 0) {
      // If distance is not zero, update new structures
      newDepths.push_back(parentDepths[parent] + 1);
      newRowVertices.push_back(node);
      newParents.push_back(parent);
      newNodeWeights.push_back(newParams.nodeWeights[node]);
      newEdgeWeights.push_back(newParams.edgeWeights[node]);
      newTypes.push_back(newParams.types[node]);

      newPointer[node] = newRowVertices.size() - 1;
      parentDepths[node] = parentDepths[parent] + 1;
      newInput.insert_rows(newInput.n_rows, input.row(node));

    } else {
      // Merge nodes if distance is zero
      newInput.row(parentIndex) = (newNodeWeights[parentIndex] * newInput.row(parentIndex) +
        newParams.nodeWeights[node] * input.row(node)) / (newNodeWeights[parentIndex] + newParams.nodeWeights[node]);
      // Update node weight
      newNodeWeights[parentIndex] += newParams.nodeWeights[node];

      // Update the children of the merged node
      for (int j = 0; j < newParams.childrenList[node].size(); j++) {
        int child = newParams.childrenList[node][j];
        newParams.parents[child] = newParams.parents[node];
      }
      // Remove the merged node from its parent's children
      siblings = newParams.childrenList[parent];
      newChildren.clear();
      for (int j = 0; j < siblings.size(); j++) {
        if (siblings[j] != node) {
          newChildren.push_back(siblings[j]);
        }
      }
      newChildren.insert(newChildren.end(), newParams.childrenList[node].begin(), newParams.childrenList[node].end());
      newParams.childrenList[parent] = newChildren;

      // Update the type of the parent
      if (parent == root) {
        if (newChildren.size() == 1) {
          newTypes[parentIndex] = 3;
        } else {
          newTypes[parentIndex] = 4;
        }
      } else if (newChildren.size() == 0) {
        newTypes[parentIndex] = 0;
      } else if (newChildren.size() == 1) {
        newTypes[parentIndex] = 1;
      } else {
        newTypes[parentIndex] = 2;
      }

      newPointer[node] = parentIndex;
      parentDepths[node] = parentDepths[parent];
    }
  }

  // Map old node indices to new node indices
  std::vector<double> index(n, 0);
  for (int i = 0; i < newRowVertices.size(); i++) {
    index[newRowVertices[i]] = i;
  }

  // Update children with new indices
  std::vector<std::vector<double>> newChildrenList(newRowVertices.size());
  for (int i = 0; i < newRowVertices.size(); i++) {
    int node = newRowVertices[i];
    for (int j = 0; j < newParams.childrenList[node].size(); j++) {
      newChildrenList[i].push_back(index[newParams.childrenList[node][j]]);
    }
  }
  // Update parents with new indices
  for (int i = 0; i < newRowVertices.size(); i++) {
    newParents[i] = index[newParents[i]];
  }

  // Sort vertices by the order of newDepths
  std::vector<double> idx(newDepths.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::stable_sort(idx.begin(), idx.end(), [&newDepths](int i1, int i2) {
    return newDepths[i1] > newDepths[i2];
  });

  std::vector<double> newVertices(newRowVertices.size());
  for (int i = 0; i < newRowVertices.size(); i++) {
    newVertices[i] = idx[i];
  }

  // Return updated values
  newParams.vertices = newVertices;
  newParams.types = newTypes;
  newParams.parents = newParents;
  newParams.edgeWeights = newEdgeWeights;
  newParams.nodeWeights = newNodeWeights;
  newParams.childrenList = newChildrenList;
  newParams.pointer = newPointer;
  newParams.rowVertices = newRowVertices;

  return std::make_tuple(newInput, newParams);
}


// [[Rcpp::export]]
Rcpp::List updatenew(
    const arma::mat& Theta,
    const arma::mat& input,
    const std::vector<double>& vertices,
    const std::vector<double>& types,
    const std::vector<double>& parents,
    const std::vector<double>& nodeWeights,
    const std::vector<double>& edgeWeights,
    const std::vector<std::vector<double>>& childrenList) {

  int numNodes = vertices.size();
  TreeParams params(numNodes);
  TreeParams newParams(numNodes);

  params.vertices = vertices;
  params.types = types;
  params.parents = parents;
  params.edgeWeights = edgeWeights;
  params.nodeWeights = nodeWeights;
  params.childrenList = childrenList;

  arma::mat newInput;

  // Call the update function with the new variable names
  std::tie(newInput, newParams) = update(Theta, input, params);

  // Create and return the result list
  Rcpp::List res = Rcpp::List::create(
    Named("newInput") = wrap(newInput),
    Named("newVertices") = wrap(newParams.vertices),
    Named("newTypes") = wrap(newParams.types),
    Named("newParents") = wrap(newParams.parents),
    Named("newEdgeWeights") = wrap(newParams.edgeWeights),
    Named("newNodeWeights") = wrap(newParams.nodeWeights),
    Named("newChildrenList") = wrap(newParams.childrenList),
    Named("pointer") = wrap(newParams.pointer),
    Named("rowVertices") = wrap(newParams.rowVertices)
  );

  return res;
}
