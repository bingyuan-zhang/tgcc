#ifndef tgcc_H
#define tgcc_H

#include <vector>
#include <tuple>

// Function computeTheta declaration
std::vector<double> computeTheta(
    const std::vector<double>& nodeValues,
    double lambda,
    const std::vector<double>& vertices,
    const std::vector<double>& nodeTypes,
    const std::vector<double>& parents,
    const std::vector<double>& nodeWeights,
    const std::vector<double>& edgeWeights,
    const std::vector<std::vector<double>>& childrenList);

#endif // tgcc_H
