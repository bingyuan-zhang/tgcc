#include <RcppArmadillo.h>
using namespace Rcpp;
// using namespace std;
// using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

void find_connected_edges (std::vector<double> from,
                           std::vector<double> to,
                           std::vector<double> dist,
                           std::vector<std::vector<double>>& L,
                           std::vector<std::vector<double>>& D ) {

  int edgenum = from.size();
  int i, j, k;
  double d;

  for (i = 0; i < edgenum; i++) {
    j = from[i];
    k = to[i];
    d = dist[i];

    L[j].push_back(k);
    L[k].push_back(j);
    D[j].push_back(d);
    D[k].push_back(d);
  }
}

// // Function to build adjacency and distance lists for the graph
// void find_connected_edges(
//     const std::vector<double>& sourceNodes,
//     const std::vector<double>& targetNodes,
//     const std::vector<double>& distances,
//     std::vector<std::vector<double>>& adjList,
//     std::vector<std::vector<double>>& distanceList) {
//
//   int edgeCount = sourceNodes.size();
//
//   for (int i = 0; i < edgeCount; i++) {
//     int source = sourceNodes[i];
//     int target = targetNodes[i];
//     double distance = distances[i];
//
//     adjList[source].push_back(target);
//     adjList[target].push_back(source);
//     distanceList[source].push_back(distance);
//     distanceList[target].push_back(distance);
//   }
// }
//

void initialize (std::vector<std::vector<double>> L,
                 std::vector<std::vector<double>> D,
                 std::vector<std::vector<double>>& Children,
                 std::vector<double>& Vertice,
                 std::vector<double>& Depth,
                 std::vector<double>& Type,
                 std::vector<double>& WeightE,
                 std::vector<double>& Parent,
                 std::vector<double>& Partition,
                 double gamma) {

  int n = L.size();
  int i, j, v;
  int seen = 1, depth = 0;
  std::vector<double> vertice_cur_depth = L[0];
  std::vector<double> vertice_next_depth;

  // find parameters for the root (always 0)
  Vertice[0] = 0;
  Depth[0] = 0;
  Parent[0] = 0;
  WeightE[0] = 0;
  Children[0] = L[0];
  if (L[0].size() == 1) {
    Type[0] = 3;
  } else {
    Type[0] = 4;
  }
  for (i=0; i<L[0].size(); i++) {
    // set root 0 as the parent root's children
    Parent[L[0][i]] = 0;
    // set edge weight between root and its children
    WeightE[L[0][i]] = std::exp(-std::pow(D[0][i], 2)/gamma);
  }

  while (seen < n) {
    depth = depth + 1;
    for (i=0; i < vertice_cur_depth.size(); i++) {
      v = vertice_cur_depth[i];
      Vertice[seen] = v;
      Depth[seen] = depth;
      for (j=0; j < L[v].size(); j++){
        if (L[v][j] != Parent[v]) {
          // find v's children
          Children[v].push_back(L[v][j]);
          // set v as the parent v's children
          Parent[L[v][j]] = v;
          // set edge weight between v and its children
          WeightE[L[v][j]] = std::exp(-std::pow(D[v][j], 2)/gamma);
        }
      }
      // std::cout <<"\n";
      if (Children[v].size() == 0) {
        Type[v] = 0;
      } else if (Children[v].size() == 1) {
        Type[v] = 1;
      } else {
        Type[v] = 2;
      }
      seen = seen + 1;
      vertice_next_depth.insert(vertice_next_depth.end(), Children[v].begin(), Children[v].end());
    }
    vertice_cur_depth.clear();
    vertice_cur_depth = vertice_next_depth;
    vertice_next_depth.clear();
  }

  std::reverse(Vertice.begin(), Vertice.end());
  std::reverse(Depth.begin(), Depth.end());

  for (i = 0; i < n; i++) {
    v = Vertice[i];
    if (Type[v] != 0) {
      if (Type[v] == 2) {
        for(j = 0; j < Children[v].size(); j++) {
          Partition[v] += Partition[Children[v][j]];
        }
      } else if (Type[v] != 3) { // Type = 1/4
        Partition[v] += Partition[Children[v][0]];
      }
    }
  }

  for (i = 0; i < n; i++) {
    int m = n - Partition[i];
    if (Partition[i] > m) {
      Partition[i] = m;
    }
  }
}

//[[Rcpp::export]]
Rcpp::List init_prepare (std::vector<double> from,
                         std::vector<double> to,
                         std::vector<double> dist,
                         double gamma) {

  int n = from.size() + 1;
  std::vector<std::vector<double>> L (n);
  std::vector<std::vector<double>> Children (n);
  std::vector<std::vector<double>> D (n);

  // The Connected Nodes.
  find_connected_edges(from, to, dist, L, D);

  std::vector<double> Vertice (n);
  std::vector<double> Depth (n);
  std::vector<double> Type (n);
  std::vector<double> WeightE (n);
  std::vector<double> WeightV (n, 1);
  std::vector<double> Parent (n);
  std::vector<double> Partition (n, 1);

  // Initialize the above parameters.
  initialize(L, D, Children, Vertice, Depth,  Type, WeightE, Parent, Partition, gamma);

  Rcpp::List res = List::create(
    Named("Vertice")   = wrap(Vertice),
    Named("Depth")     = wrap(Depth),
    Named("Type")      = wrap(Type),
    Named("WeightE")   = wrap(WeightE),
    Named("WeightV")   = wrap(WeightV),
    Named("Parent")    = wrap(Parent),
    Named("Partition") = wrap(Partition),
    Named("Children")  = wrap(Children)
  );

  return (res);
}

void update (arma::mat Theta,
             arma::mat& input,
             std::vector<double>& V,
             std::vector<double>& T,
             std::vector<double>& P,
             std::vector<double>& WV,
             std::vector<double>& WE,
             std::vector<std::vector<double>>& C,
             std::vector<double>& pointer,
             std::vector<double>& rowV) {

  int n = V.size();
  std::vector<double> newrowV, newP, newT, newWV, newD, newWE;
  std::vector<double> newpointer (n);
  std::vector<double> pdepth (n);
  arma::mat newinput;

  int root = V[n-1];
  int i, j, node, parent, child, pind;
  std::vector<double> brothers, newpc;
  double d;

  newrowV.push_back(V[n-1]);
  newP.push_back(P[0]);
  newT.push_back(T[0]);
  newWV.push_back(WV[0]);
  newWE.push_back(0);
  newD.push_back(0);
  pdepth[0] = 0;
  newpointer[0] = 0;
  newinput.insert_rows(0, input.row(0));

  for (i = n-2; i >= 0; i--) {

    node = V[i];
    parent = P[node];
    pind = newpointer[parent];

    d = arma::norm(Theta.row(node) - Theta.row(parent), 2);

    if (d != 0) {
      newD.push_back(pdepth[parent] + 1);
      newrowV.push_back(node);
      newP.push_back(parent);
      newWV.push_back(WV[node]);
      newWE.push_back(WE[node]);
      newT.push_back(T[node]);

      newpointer[node] = newrowV.size() - 1;
      pdepth[node] = pdepth[parent] + 1;
      newinput.insert_rows(newinput.n_rows, input.row(node));

    } else {

      // update input value
      newinput.row(pind) = (newWV[pind]*newinput.row(pind) + WV[node]*input.row(node))/(newWV[pind] + WV[node]);

      // update vertice weight
      newWV[pind] = newWV[pind] + WV[node];

      // update the children of the new node
      // update the parent of the nodes to be their grand parent.
      for (j = 0; j < C[node].size(); j++) {
        child = C[node][j];
        P[child] = P[node];
      }
      // remove the merged node from its parent's children
      brothers = C[parent];
      newpc.clear();
      for (j = 0; j < brothers.size(); j++) {
        if (brothers[j] != node) {
          newpc.push_back(brothers[j]);
        }
      }
      newpc.insert(newpc.end(), C[node].begin(), C[node].end());
      C[parent] = newpc;

      // update the class of the parent
      if (parent == root) {
        if (newpc.size() == 1) {
          newT[pind] = 3;
        } else {
          newT[pind] = 4;
        }
      } else if (newpc.size() == 0) {
        newT[pind] = 0;
      } else if (newpc.size() == 1) {
        newT[pind] = 1;
      } else {
        newT[pind] = 2;
      }

      newpointer[node] = pind;
      pdepth[node] = pdepth[parent];
    }
  }

  std::vector<double> index (n, 0);
  for (i = 0; i < newrowV.size(); i++) {
    index[newrowV[i]] = i;
  }

  std::vector<std::vector<double>> newC (newrowV.size());
  for (i = 0; i < newrowV.size(); i++) {
    node = newrowV[i];
    for (j = 0; j < C[node].size(); j++) {
      newC[i].push_back(index[C[node][j]]);
    }
  }

  for (i = 0; i < newrowV.size(); i++) {
    newP[i] = index[newP[i]];
  }

  // sort V by the order of D
  std::vector<double> idx (newD.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::stable_sort(idx.begin(), idx.end(),
    [&newD](int i1, int i2) {return newD[i1] > newD[i2];});

  std::vector<double> newV (newrowV.size());
  for (int i = 0; i < newrowV.size(); i++) {
    newV[i] = idx[i];
  }

  // return value
  input = newinput;
  V = newV;
  T = newT;
  P = newP;
  WE = newWE;
  WV = newWV;
  C = newC;
  pointer = newpointer;
  rowV = newrowV;

}


// [[Rcpp::export]]
Rcpp::List updatenew (arma::mat Theta,
                      arma::mat& input,
                      std::vector<double>& V,
                      std::vector<double>& T,
                      std::vector<double>& P,
                      std::vector<double>& WV,
                      std::vector<double>& WE,
                      std::vector<std::vector<double>>& C) {

  std::vector<double> pointer;
  std::vector<double> rowV;

  update (Theta, input, V, T, P, WV, WE, C, pointer, rowV);

  Rcpp::List res = List::create(
    Named("newinput") = wrap(input),
    Named("newV") = wrap(V),
    Named("newT") = wrap(T),
    Named("newP") = wrap(P),
    Named("newWE") = wrap(WE),
    Named("newWV") = wrap(WV),
    Named("newC") = wrap(C),
    Named("pointer") = wrap(pointer),
    Named("rowofV") = wrap(rowV)
  );

  return (res);
}



void forwardleaf (double xnode,
                  double wv,
                  double lam,
                  std::deque<double>& BP,
                  std::deque<double>& I,
                  std::deque<double>& S,
                  double& l,
                  double& u) {

  l = xnode - lam / wv;
  u = xnode + lam / wv;

  // bp = c(L,U)
  BP.push_front(l);
  BP.push_back(u);

  // I = c(-lam, -wv*xnode, lam)
  I.push_front(-lam);
  I.push_back(-wv*xnode);
  I.push_back(lam);

  // s = c(0, wv, 0)
  S.push_front(0);
  S.push_back(wv);
  S.push_back(0);

}

void forwardroot (double xnode,
                  double wv,
                  std::deque<double>& BP,
                  std::deque<double>& I,
                  std::deque<double>& S,
                  double& theta) {

  theta = (wv*xnode - I.front()) / (S.front() + wv);

  while (BP.front() < theta) {
    I.pop_front();
    S.pop_front();
    BP.pop_front();
    theta = (wv*xnode - I.front()) / (S.front() + wv);
    if(S.size() == 1) break;
  }

}

void forward (double xnode,
              double wv,
              double lam,
              std::deque<double>& BP,
              std::deque<double>& I,
              std::deque<double>& S,
              double& l,
              double& u) {

  l = (-lam + wv*xnode - I.front()) / (S.front() + wv);

  while (BP.front() < l) {
    I.pop_front();
    S.pop_front();
    BP.pop_front();
    l = (-lam + wv*xnode - I.front()) / (S.front() + wv);
    if(S.size() == 1) break;
  }

  if (S.size() == 1) {
    u = (lam + wv*xnode - I.front()) / (S.front() + wv);
  } else {
    u = (lam + wv*xnode - I.back()) / (S.back() + wv);
    while (BP.back() > u) {
      I.pop_back();
      S.pop_back();
      BP.pop_back();
      u = (lam + wv*xnode - I.back()) / (S.back() + wv);
      if (S.size() == 1) break;
    }
  }

  // BP = c(l, BP, u)
  BP.push_front(l);
  BP.push_back(u);

  // I = c(-lam, I -wv*xnode, lam)
  for (int i=0; i<I.size(); i++) {
    I[i] = I[i] - wv*xnode;
  }
  I.push_front(-lam);
  I.push_back(lam);

  // S = c(0, S+wv, 0)
  for (int i=0; i<S.size(); i++) {
    S[i] = S[i] + wv;
  }
  S.push_front(0);
  S.push_back(0);

}

void mergequeues (std::vector<double> children,
                  int index,
                  std::vector<std::deque<double>>& BPs,
                  std::vector<std::deque<double>>& Is,
                  std::vector<std::deque<double>>& Ss) {

  int ns = children.size();
  int child;
  std::vector<double> bp_len (ns);
  std::deque<double> new_bp;

  for (int i=0; i<ns; i++) {
    child = children[i];
    bp_len[i] = BPs[child].size();
    new_bp.insert(std::end(new_bp), std::begin(BPs[child]), std::end(BPs[child]));
  }
  int l = new_bp.size();

  // obtain the rank
  std::vector< std::pair<double, int> > bp_sort(l);
  std::vector<double> rank(l);
  for (int i=0; i<l; i++) {
    bp_sort[i] =  std::make_pair(new_bp[i], i);
  }
  std::sort(bp_sort.begin(), bp_sort.end());
  for (int i=0; i<l; i++) {
    rank[bp_sort[i].second] = i;
  }
  std::sort(new_bp.begin(), new_bp.end());

  // merge queues
  int s1 = 0;
  int s2 = 0;
  std::vector<std::vector<double>> i_tosum (ns, std::vector<double> (l + 1)); // included the extended sequence
  std::vector<std::vector<double>> s_tosum (ns, std::vector<double> (l + 1));
  for (int i = 0; i < ns; i++) {
    child = children[i];
    s2 = 0;
    for (int j = 0; j < bp_len[i]; j++) {
      std::fill(i_tosum[i].begin() + s2, i_tosum[i].begin() + rank[s1+j] + 1, Is[child][j]);
      std::fill(s_tosum[i].begin() + s2, s_tosum[i].begin() + rank[s1+j] + 1, Ss[child][j]);
      s2 = rank[s1+j] + 1;
    }
    std::fill(i_tosum[i].begin() + s2, i_tosum[i].end(), Is[child][bp_len[i]]);
    std::fill(s_tosum[i].begin() + s2, s_tosum[i].end(), Ss[child][bp_len[i]]);
    s1 = s1 + bp_len[i];
  }

  std::deque<double> new_i (l + 1, 0.0);
  std::deque<double> new_s (l + 1, 0.0);

  for (int j = 0; j < l+1; j++) {
    for (int i = 0; i < ns; i++) {
      new_i[j] = new_i[j] + i_tosum[i][j];
      new_s[j] = new_s[j] + s_tosum[i][j];
    }
  }

  BPs[index] = new_bp;
  Is[index] = new_i;
  Ss[index] = new_s;
}


// [[Rcpp::export]]
std::vector<double> dp (std::vector<double> x,
                        double lam,
                        std::vector<double> Vertice,
                        std::vector<double> Type,
                        std::vector<double> Parent,
                        std::vector<double> WeightV,
                        std::vector<double> WeightE,
                        std::vector<std::vector<double>> Children) {

  int n = x.size();
  std::vector<std::deque<double>> BPs (n); // branch points
  std::vector<std::deque<double>> Is (n); // intercepts
  std::vector<std::deque<double>> Ss (n); // slopes
  std::vector<double> L (n);
  std::vector<double> U (n);

  int v;
  int t;
  int p;
  double wv;
  double we;
  double xnode;

  int child;
  std::vector<double> children;
  std::vector<double> theta (n, 0);

  for (int i=0; i < n; i++) {

    v = Vertice[i];
    t = Type[v];
    wv = WeightV[v];
    we = WeightE[v];
    xnode = x[v];

    if (t == 0) {
      forwardleaf(xnode, wv, we*lam, BPs[v], Is[v], Ss[v], L[i], U[i]);
    } else if (t == 1) {
      child = Children[v][0];
      BPs[v] = BPs[child];
      Is[v] = Is[child];
      Ss[v] = Ss[child];
      forward(xnode, wv, we*lam, BPs[v], Is[v], Ss[v], L[i], U[i]);
    } else if (t == 2) {
      children = Children[v];
      mergequeues(children, v, BPs, Is, Ss);
      forward(xnode, wv, we*lam, BPs[v], Is[v], Ss[v], L[i], U[i]);
    } else {
      if (t == 3) {
        child = Children[v][0];
        forwardroot(xnode, wv, BPs[child], Is[child], Ss[child], theta[v]);
      } else {
        children = Children[v];
        mergequeues(children, v, BPs, Is, Ss);
        forwardroot(xnode, wv, BPs[v], Is[v], Ss[v], theta[v]);
      }
    }
  }

  //Backward pass
  for (int i=n-2; i>=0; i--){
    v = Vertice[i];
    p = Parent[v];
    if (theta[p] < L[i]){
      theta[v] = L[i];
    } else if (theta[p] > U[i]) {
      theta[v] = U[i];
    } else {
      theta[v] = theta[p];
    }
  }

  return(theta);
}




