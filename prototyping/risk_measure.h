// This file contains the risk measure class, plus a sorting function

#include <iostream>
#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort

// Risk measure class
class RiskMeasure
{
protected:
    std::vector<double> lambdas;
    std::vector<double> betas;

public:
    // constructor
    RiskMeasure(std::vector<double> lambdas_new, std::vector<double> betas_new);

    // destructor
    ~RiskMeasure();

    // compute phi
    double phi(double p);

    // compute vector of weights
    std::vector<double> compute_weights(int num_scenarios);

    // compute rho
    double compute_rho(std::vector<double> scenario_values);

};


// Sorting function
template <typename T> // this indicates that "T" is a typename
std::vector<std::size_t> sort_indices(const std::vector<T> &v);