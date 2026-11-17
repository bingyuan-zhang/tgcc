# tgcc 0.2.1

* The default FAMST settings of `tgCC()` are now shown in its signature:
  `famstControl = list(nNeighbors = 15, numRandomEdges = 3, metric = "euclidean")`.
  Results are unchanged.

# tgcc 0.2.0

* New function `famst()`: an R/C++ re-implementation of FAMST (fast
  approximate minimum spanning tree; Almansoori and Telek, 2025), with the
  approximate k-nearest-neighbour search provided by `rnndescent`.
* `tgCC()` gains the argument `mst`. `mst = "famst"` builds the tree by
  `famst()` instead of the exact Euclidean MST (`mlpack::emst`, the default),
  which removes the cost of exact MST construction for higher-dimensional
  data. The FAMST settings are passed through `famstControl`
  (defaults: `nNeighbors = 15`, `numRandomEdges = 3`, `metric = "euclidean"`).
* `tgCC()` gains the argument `keepTheta`. With `keepTheta = FALSE` the
  fitted matrices for every value of lambda are not returned, which reduces
  memory use for large data; `clusterLabel()` does not need them.
* The output of `tgCC()` now records the tree construction method (`mst`).
* `rnndescent` is a new dependency.

# tgcc 0.1

* Initial version.
