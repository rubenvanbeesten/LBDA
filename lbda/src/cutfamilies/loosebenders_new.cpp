#include "cutfamilies/loosebenders_new.h"

#include "subproblem.h"

#include "risk_measure.h"

#include <algorithm>
#include <vector>


LooseBenders::LooseBenders(ProblemData const &problem,
                           RiskMeasure const &risk_measure,
                           arma::vec const alpha,
                           double timeLimit) :
    CutFamily(problem),
    d_alpha(alpha),
    d_risk_measure(risk_measure),
    d_visited(problem.nScenarios()),
    d_objectives(problem.nScenarios())
{
    auto const &Wmat = d_problem.Wmat();

    d_vars = d_model.addVars(d_problem.secondStageLowerBound().memptr(),
                             d_problem.secondStageUpperBound().memptr(),
                             problem.secondStageCoeffs().memptr(),
                             problem.secondStageVarTypes().memptr(),
                             nullptr,
                             problem.secondStageCoeffs().n_elem);

    GRBLinExpr lhs[Wmat.n_cols];

    for (auto iter = Wmat.begin(); iter != Wmat.end(); ++iter)
        lhs[iter.col()] += *iter * d_vars[iter.row()];  // Wy

    auto const &senses = d_problem.secondStageConstrSenses();
    arma::vec rhs = arma::zeros(senses.n_elem);

    d_constrs = d_model.addConstrs(lhs,
                                   senses.memptr(),
                                   rhs.memptr(),
                                   nullptr,
                                   rhs.n_elem);

    d_model.set(GRB_DoubleParam_TimeLimit, timeLimit);
    d_model.update();
}

LooseBenders::~LooseBenders()
{
    delete[] d_vars;
    delete[] d_constrs;
}

LooseBenders::Cut LooseBenders::computeCut(arma::vec const &x)
{
    auto const &Tmat = d_problem.Tmat(); 
    arma::vec Tx = Tmat.t() * x;

    // NOTE: a cut has the following form:
    // 

    // OLD: STORE IN A SINGLE "DUAL" AND "GAMMA"   
    /*
    arma::vec dual = arma::zeros(Tmat.n_cols); // this is lambda
    double gamma = 0;
    */

    // NEW: STORE CUTS IN VECTORS
    std::vector<double>    deltaVec(d_problem.nScenarios(), 0.0);                   // vector of delta constants; one for each scenario
    std::vector<arma::vec> betaVec(d_problem.nScenarios(), arma::zeros(x.n_elem));  // vector of beta coefficient vectors; one for each scenario
    std::vector<double>    uVec(d_problem.nScenarios(), 0.0);                       // vector of "u" values to be sorted (i.e., the value function objective values)


    for (size_t scenario = 0; scenario != d_problem.nScenarios(); ++scenario)
    {
        arma::vec omega = d_problem.scenarios().col(scenario);

        d_sub.updateRhs(omega - Tx);
        d_sub.solve();

        auto const basis = d_sub.basisInfo();
        auto const duals = d_sub.duals();
        //auto const curObj = d_sub.objective(); // objective value of the current scenario

        //double const prob = d_problem.scenarioProbability(scenario);

        // Gomory is lambda^T (omega - alpha) + psi(omega - alpha), so we add
        // lambda^T alpha.
        arma::vec rhs = omega - d_alpha;

        double res = computeGomory(scenario, rhs, basis.vBasis, basis.cBasis);
        res += arma::dot(duals.lambda, d_alpha);

        // extract cut coefficients of this scenario // NOTE: chedk if duals.lambdas is correct below
        double    cur_delta = res;                             // current value of delta (constant in the cut)
        arma::vec cur_beta = - Tmat * duals.lambda;            // current value of beta (x coefficient in the cut) //NOTE: minus sign taken care of
        double    cur_u = cur_delta + arma::dot(cur_beta, x);  // current value of "u", i.e., th value of the cut at current x

        // store values in vectors
        deltaVec[scenario] = cur_delta;  // store current delta value
        betaVec[scenario] = cur_beta;   // store current beta value
        uVec[scenario] = cur_u;         // store current cut value (to be sorted)
    }
  
    // TO DO:
    // 1. compute weights associated with risk measure (per scenario)                    --- function in risk measure class
    std::vector<double> weights = d_risk_measure.compute_weights(d_problem.nScenarios());

    // 2. sort uVec; store indices                                                     --- use function sort_indices

    // input vector: std::vector<double> uVec 

    std::vector<int> sorted_idx(uVec.size()); // vector with sorted indices corresponding to increasing values in uVec
    int idx_counter = 0;
    for (auto i: sort_indices(uVec)) { 
        sorted_idx.at(idx_counter) = i;
        idx_counter++;
    }

    // TEMPORARY: this is how one would output uVec in sorted order
    for (size_t i=0; i < sorted_idx.size(); i++){
        std::cout << sorted_idx.at(i) << ": " << uVec.at(sorted_idx.at(i)) << std::endl;
    }

    // 3. multiply weights by associated gamma and lambda values to create the cut
    
    // initialize
    arma::vec beta = arma::zeros(Tmat.n_rows);  // initialize as vector of zeros
    double    delta = 0.0;                      // initialize as zero

    // multiply scenario values by corresponding weights and aggregate
    for(size_t s = 0; s < d_problem.nScenarios(); s++) {
        double cur_weight = weights.at(sorted_idx.at(s));
        delta += cur_weight * deltaVec.at(s);
        for(size_t i=0; i < Tmat.n_rows; i++){ // loop over elements of beta
            beta[i] += cur_weight * betaVec.at(s).at(i);
        }
    }

    
    // 4. return the cut
    return Cut{beta, delta}; // we use "delta" instead of "gamma"

    // OLD:
    //return Cut{Tmat * dual, gamma};



}

double LooseBenders::computeGomory(size_t scenario,
                                   arma::vec &rhs,
                                   arma::Col<int> const &vBasis,
                                   arma::Col<int> const &cBasis)
{
    auto const &Wmat = d_problem.Wmat();

    // TODO make this more efficient?
    std::vector<int> basis(Wmat.n_rows + Wmat.n_cols);
    std::copy(vBasis.memptr(), vBasis.memptr() + Wmat.n_rows, basis.begin());
    std::copy(cBasis.memptr(),
              cBasis.memptr() + Wmat.n_cols,
              basis.begin() + Wmat.n_rows);

    auto &visited = d_visited[scenario];
    auto iterator = std::find(visited.begin(), visited.end(), basis);

    if (iterator != visited.end())
    {
        // Retrieve index and corresponding objective value
        size_t idx = std::distance(visited.begin(), iterator);
        return d_objectives[scenario][idx];
    }

    update(rhs, vBasis, cBasis);
    d_model.optimize();

    // Either optimal objective value, or best known lower bound (the latter
    // often in case of a time out).
    auto const objective = d_model.get(GRB_IntAttr_Status) == GRB_OPTIMAL
                               ? d_model.get(GRB_DoubleAttr_ObjVal)
                               : d_model.get(GRB_DoubleAttr_ObjBound);

    visited.emplace_back(basis);
    d_objectives[scenario].emplace_back(objective);

    return objective;
}

void LooseBenders::update(arma::vec &rhs,
                          arma::Col<int> const &vBasis,
                          arma::Col<int> const &cBasis)
{
    auto const &senses = d_problem.secondStageConstrSenses();

    for (size_t idx = 0; idx != senses.size(); ++idx)
    {
        if (cBasis[idx] != GRB_BASIC || senses[idx] == GRB_EQUAL)
            continue;

        if (senses[idx] == GRB_LESS_EQUAL)
            rhs[idx] = arma::datum::inf;

        if (senses[idx] == GRB_GREATER_EQUAL)
            rhs[idx] = -arma::datum::inf;
    }

    d_model.set(GRB_DoubleAttr_RHS, d_constrs, rhs.memptr(), rhs.n_elem);

    arma::vec lb = d_problem.secondStageLowerBound();
    arma::vec ub = d_problem.secondStageUpperBound();

    // Relax appropriate variable bounds if the bound is not tight.
    lb.elem(arma::find(vBasis != GRB_NONBASIC_LOWER)).fill(-arma::datum::inf);
    ub.elem(arma::find(vBasis != GRB_NONBASIC_UPPER)).fill(arma::datum::inf);

    d_model.set(GRB_DoubleAttr_LB, d_vars, lb.memptr(), lb.n_elem);
    d_model.set(GRB_DoubleAttr_UB, d_vars, ub.memptr(), ub.n_elem);
}
