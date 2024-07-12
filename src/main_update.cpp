#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// Function to update the tree structure and weights
void update (
    arma::mat Theta,
    arma::mat& input,
    std::vector<double>& vertices,
    std::vector<double>& types,
    std::vector<double>& parents,
    std::vector<double>& nodeWeights,
    std::vector<double>& edgeWeights,
    std::vector<std::vector<double>>& childrenList,
    std::vector<double>& pointer,
    std::vector<double>& rowVertices) {

  int n = vertices.size();
  std::vector<double> newRowVertices, newParents, newTypes, newNodeWeights, newDepths, newEdgeWeights;
  std::vector<double> newPointer(n, 0);
  std::vector<double> parentDepths(n, 0);
  arma::mat newInput;

  int root = vertices[n - 1]; // Get the root node
  int node, parent, child, parentIndex;
  double distance;
  std::vector<double> siblings, newChildren;

  newRowVertices.push_back(root);
  newParents.push_back(parents[0]);
  newTypes.push_back(types[0]);
  newNodeWeights.push_back(nodeWeights[0]);
  newEdgeWeights.push_back(0);
  newDepths.push_back(0);
  parentDepths[0] = 0;
  newPointer[0] = 0;
  newInput.insert_rows(0, input.row(0));

  for (int i = n - 2; i >= 0; i--) {

    node = vertices[i];
    parent = parents[node];
    parentIndex = newPointer[parent];

    // Calculate the distance between the node and its parent
    distance = arma::norm(Theta.row(node) - Theta.row(parent), 2);

    if (distance != 0) {
      // If distance is not zero, update new structures
      newDepths.push_back(parentDepths[parent] + 1);
      newRowVertices.push_back(node);
      newParents.push_back(parent);
      newNodeWeights.push_back(nodeWeights[node]);
      newEdgeWeights.push_back(edgeWeights[node]);
      newTypes.push_back(types[node]);

      newPointer[node] = newRowVertices.size() - 1;
      parentDepths[node] = parentDepths[parent] + 1;
      newInput.insert_rows(newInput.n_rows, input.row(node));

    } else {
      // Merge nodes if distance is zero
      newInput.row(parentIndex) = (newNodeWeights[parentIndex] * newInput.row(parentIndex) +
        nodeWeights[node] * input.row(node)) / (newNodeWeights[parentIndex] + nodeWeights[node]);
      // Update node weight
      newNodeWeights[parentIndex] += nodeWeights[node];

      // Update the children of the merged node
      for (int j = 0; j < childrenList[node].size(); j++) {
        child = childrenList[node][j];
        parents[child] = parents[node];
      }
      // Remove the merged node from its parent's children
      siblings = childrenList[parent];
      newChildren.clear();
      for (int j = 0; j < siblings.size(); j++) {
        if (siblings[j] != node) {
          newChildren.push_back(siblings[j]);
        }
      }
      newChildren.insert(newChildren.end(), childrenList[node].begin(), childrenList[node].end());
      childrenList[parent] = newChildren;

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
    node = newRowVertices[i];
    for (int j = 0; j < childrenList[node].size(); j++) {
      newChildrenList[i].push_back(index[childrenList[node][j]]);
    }
  }
  // Update parents with new indices
  for (int i = 0; i < newRowVertices.size(); i++) {
    newParents[i] = index[newParents[i]];
  }

  // Sort vertices by the order of newDepths
  std::vector<double> idx(newDepths.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::stable_sort(idx.begin(), idx.end(), [&newDepths](int i1, int i2) { return newDepths[i1] > newDepths[i2]; });

  std::vector<double> newVertices(newRowVertices.size());
  for (int i = 0; i < newRowVertices.size(); i++) {
    newVertices[i] = idx[i];
  }

  // Return updated values
  input = newInput;
  vertices = newVertices;
  types = newTypes;
  parents = newParents;
  edgeWeights = newEdgeWeights;
  nodeWeights = newNodeWeights;
  childrenList = newChildrenList;
  pointer = newPointer;
  rowVertices = newRowVertices;

}


// [[Rcpp::export]]
Rcpp::List updatenew(
    arma::mat Theta,
    arma::mat& input,
    std::vector<double>& vertices,
    std::vector<double>& types,
    std::vector<double>& parents,
    std::vector<double>& nodeWeights,
    std::vector<double>& edgeWeights,
    std::vector<std::vector<double>>& childrenList) {

  std::vector<double> pointer;
  std::vector<double> rowVertices;

  // Call the update function with the new variable names
  update(Theta, input, vertices, types, parents, nodeWeights, edgeWeights, childrenList, pointer, rowVertices);

  // Create and return the result list
  Rcpp::List res = Rcpp::List::create(
    Named("newInput") = wrap(input),
    Named("newVertices") = wrap(vertices),
    Named("newTypes") = wrap(types),
    Named("newParents") = wrap(parents),
    Named("newEdgeWeights") = wrap(edgeWeights),
    Named("newNodeWeights") = wrap(nodeWeights),
    Named("newChildrenList") = wrap(childrenList),
    Named("pointer") = wrap(pointer),
    Named("rowVertices") = wrap(rowVertices)
  );

  return res;
}



