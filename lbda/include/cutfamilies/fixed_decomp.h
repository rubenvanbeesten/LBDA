#ifndef FIXED_DECOMP_H
#define FIXED_DECOMP_H

#include "cutfamily.h"
#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort, std::stable_sort

#include "risk_measure.h"

class FixedDecomp : public CutFamily
{
    arma::vec const d_alpha;

    GRBConstr *d_constrs;
    GRBVar *d_vars;

    RiskMeasure d_risk_measure;

    // For each scenario, we store the basis matrices that we have visited
    // (encoded by vBasis, cBasis).
    // std::vector<std::vector<std::vector<int>>> d_visited;

    // For each visited basis matrix, we store the corresponding Gomory
    // objective value.
    // std::vector<std::vector<double>> d_objectives;

    /*double computeGomory(size_t scenario,
                         arma::vec &rhs,
                         arma::Col<int> const &vBasis,
                         arma::Col<int> const &cBasis);

    void update(arma::vec &rhs,
                arma::Col<int> const &vBasis,
                arma::Col<int> const &cBasis);
    */

public:
    /*
    LooseBenders(ProblemData const &problem,
                 RiskMeasure const &risk_measure,
                 arma::vec const alpha,
                 double timeLimit = 1e2);

    ~LooseBenders() override;

    CutFamily::Cut computeCut(arma::vec const &x) override;
    */
};

#endif  // LOOSEBENDERS_H
