#ifndef _tgcc_extension_H
#define _tgcc_extension_H

#define TGCC_SPCC_BICC


std::vector<double> dp (std::vector<double> x, double lam,
  std::vector<double> Vertice, std::vector<double> Type, std::vector<double> Parent,
  std::vector<double> WeightV, std::vector<double> WeightE,
  std::vector<std::vector<double>> Children);


#endif
