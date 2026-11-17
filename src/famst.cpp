// famst.cpp -- C++ routines used by famst() (R/famst.R).
//
//   famstComponents : symmetric adjacency (CSR) of a kNN graph and its
//                     connected components.
//   famstKruskal    : Kruskal's minimum spanning tree with union-find.
//
// All node indices are 0-based.

#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <numeric>
using namespace Rcpp;

// ---- union-find helpers --------------------------------------------------
static int ufFind(std::vector<int>& parent, int x) {
  int root = x;
  while (parent[root] != root) root = parent[root];
  while (parent[x] != root) { int nx = parent[x]; parent[x] = root; x = nx; }
  return root;
}

static void ufUnion(std::vector<int>& parent, std::vector<int>& rank, int x, int y) {
  int rx = ufFind(parent, x), ry = ufFind(parent, y);
  if (rx == ry) return;
  if (rank[rx] < rank[ry]) std::swap(rx, ry);
  parent[ry] = rx;
  if (rank[rx] == rank[ry]) rank[rx]++;
}

// --------------------------------------------------------------------------
// Connected components of the undirected kNN graph.
//   neighbors : n x k matrix of 0-based neighbour indices (self excluded)
//   n         : number of nodes
// Returns the component id of each node, the CSR adjacency
// (adj_offsets, adj_indices) and the number of components.
// --------------------------------------------------------------------------
// [[Rcpp::export]]
List famstComponents(IntegerMatrix neighbors, int n) {
  int k = neighbors.ncol();

  // 1. Symmetric adjacency, then CSR with per-row sort + unique.
  std::vector< std::vector<int> > adj(n);
  for (int i = 0; i < n; ++i) {
    for (int c = 0; c < k; ++c) {
      int j = neighbors(i, c);
      if (j < 0 || j >= n) continue;
      adj[i].push_back(j);
      adj[j].push_back(i);
    }
  }
  IntegerVector offsets(n + 1);
  std::vector<int> indices;
  indices.reserve((size_t)n * k * 2);
  for (int i = 0; i < n; ++i) {
    std::vector<int>& a = adj[i];
    std::sort(a.begin(), a.end());
    a.erase(std::unique(a.begin(), a.end()), a.end());
    offsets[i + 1] = offsets[i] + (int)a.size();
    indices.insert(indices.end(), a.begin(), a.end());
  }

  // 2. Connected components via iterative DFS over the CSR adjacency.
  IntegerVector compId(n, -1);
  std::vector<int> stack;
  stack.reserve(n);
  int cid = 0;
  for (int s = 0; s < n; ++s) {
    if (compId[s] != -1) continue;
    stack.clear();
    stack.push_back(s);
    compId[s] = cid;
    while (!stack.empty()) {
      int v = stack.back(); stack.pop_back();
      for (int p = offsets[v]; p < offsets[v + 1]; ++p) {
        int w = indices[p];
        if (compId[w] == -1) { compId[w] = cid; stack.push_back(w); }
      }
    }
    ++cid;
  }

  return List::create(
    _["comp_id"]      = compId,
    _["adj_offsets"]  = offsets,
    _["adj_indices"]  = wrap(indices),
    _["n_components"] = cid
  );
}

// --------------------------------------------------------------------------
// Kruskal's minimum spanning tree over an edge list.
//   u, v : 0-based endpoints;  w : edge weights;  n : number of nodes (n >= 2)
// Returns the tree edges (u, v, w), the total weight and the number of
// connected components of the resulting forest.
// --------------------------------------------------------------------------
// [[Rcpp::export]]
List famstKruskal(IntegerVector u, IntegerVector v, NumericVector w, int n) {
  int m = w.size();
  std::vector<int> ord(m);
  std::iota(ord.begin(), ord.end(), 0);
  std::sort(ord.begin(), ord.end(),
            [&](int a, int b) { return w[a] < w[b]; });

  std::vector<int> parent(n), rank(n, 0);
  std::iota(parent.begin(), parent.end(), 0);

  IntegerVector mstU(n - 1), mstV(n - 1);
  NumericVector mstW(n - 1);
  int cnt = 0;
  double total = 0.0;
  for (int idx = 0; idx < m && cnt < n - 1; ++idx) {
    int e = ord[idx];
    int a = u[e], b = v[e];
    if (ufFind(parent, a) != ufFind(parent, b)) {
      ufUnion(parent, rank, a, b);
      mstU[cnt] = a; mstV[cnt] = b; mstW[cnt] = w[e];
      total += w[e];
      ++cnt;
    }
  }
  // Trim if the edge set did not connect all nodes.
  if (cnt < n - 1) {
    if (cnt == 0) {
      mstU = IntegerVector(0); mstV = IntegerVector(0); mstW = NumericVector(0);
    } else {
      mstU = mstU[Range(0, cnt - 1)];
      mstV = mstV[Range(0, cnt - 1)];
      mstW = mstW[Range(0, cnt - 1)];
    }
  }
  std::vector<int> roots;
  roots.reserve(n);
  for (int i = 0; i < n; ++i) roots.push_back(ufFind(parent, i));
  std::sort(roots.begin(), roots.end());
  int ncomp = (int)(std::unique(roots.begin(), roots.end()) - roots.begin());

  return List::create(
    _["u"] = mstU, _["v"] = mstV, _["w"] = mstW,
    _["weight"] = total, _["num_components"] = ncomp
  );
}
