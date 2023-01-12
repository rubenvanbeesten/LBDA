#ifndef RISK_MEASURE_H
#define RISK_MEASURE_H

#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <sstream>

// This file contains the risk measure class, plus a sorting function and a string to vector function


// Risk measure class
class RiskMeasure
{
public:
    // main properties
    std::vector<double> lambdas;
    std::vector<double> betas;

    // constructor
    RiskMeasure(std::vector<double> lambdas_new, std::vector<double> betas_new);

    // destructor
    ~RiskMeasure();

    // sensibility check
    void sensibilityCheck();

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

#endif // RISK_MEASURE_H


// String to vector function
std::vector<double> str_to_double_vec(std::string s);
