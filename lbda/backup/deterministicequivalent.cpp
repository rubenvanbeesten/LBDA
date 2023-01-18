#include "deterministicequivalent.h"


DeterministicEquivalent::DeterministicEquivalent(ProblemData const &problem) :
    d_problem(problem), d_model(d_env)
{
    initFirstStage();
    initSecondStage();
}

void DeterministicEquivalent::initFirstStage()
{
    auto const &Amat = d_problem.Amat();

    // create names for x variabels
    std::string x_names[Amat.n_rows] = {""};
    for(size_t i=0; i < Amat.n_rows; i++){
        x_names[i] = "x_" + std::to_string(i);  // x_i
    }

    auto *xVars = d_model.addVars(d_problem.firstStageLowerBound().memptr(),
                                  d_problem.firstStageUpperBound().memptr(),
                                  d_problem.firstStageCoeffs().memptr(),
                                  d_problem.firstStageVarTypes().memptr(),
                                  x_names,
                                  d_problem.firstStageCoeffs().n_elem);

    GRBLinExpr lhs[Amat.n_cols];
    for (auto iter = Amat.begin(); iter != Amat.end(); ++iter)
        lhs[iter.col()] += *iter * xVars[iter.row()];

    auto *constrs = d_model
                        .addConstrs(lhs,
                                    d_problem.firstStageConstrSenses().memptr(),
                                    d_problem.firstStageRhs().memptr(),
                                    nullptr,
                                    Amat.n_cols);

    delete[] xVars;
    delete[] constrs;

    d_model.update();
}

void DeterministicEquivalent::initSecondStage()
{
    auto const &Tmat = d_problem.Tmat();
    auto const &Wmat = d_problem.Wmat();

    auto *xVars = d_model.getVars();

    for (size_t scenario = 0; scenario != d_problem.nScenarios(); ++scenario)
    {
        double const prob = d_problem.scenarioProbability(scenario);
        arma::vec const costs = prob * d_problem.secondStageCoeffs();

        GRBLinExpr lhs[Tmat.n_cols];

        for (auto iter = Tmat.begin(); iter != Tmat.end(); ++iter)
            lhs[iter.col()] += *iter * xVars[iter.row()];  // Tx

        // create names for y variabels
        std::string y_names[Wmat.n_rows] = {""};
        for(size_t i=0; i < Wmat.n_rows; i++){
            y_names[i] = "y_" + std::to_string(scenario) + "_" + std::to_string(i);  // y_s_i
        }

        // scenario-specific variables
        auto *yVars = d_model.addVars(d_problem.secondStageLowerBound().memptr(),
                                      d_problem.secondStageUpperBound().memptr(),
                                      costs.memptr(),
                                      d_problem.secondStageVarTypes().memptr(),
                                      y_names,
                                      Wmat.n_rows);

        for (auto iter = Wmat.begin(); iter != Wmat.end(); ++iter)
            lhs[iter.col()] += *iter * yVars[iter.row()];  // Wy

        auto const &senses = d_problem.secondStageConstrSenses();
        auto const rhs = d_problem.scenarios().colptr(scenario);

        auto *constrs = d_model.addConstrs(lhs,
                                           senses.memptr(),
                                           rhs,
                                           nullptr,
                                           Wmat.n_cols);

        delete[] yVars;
        delete[] constrs;
    }

    delete[] xVars;

    d_model.update();

    // write the model
    d_model.write("debug.lp");

}

std::unique_ptr<arma::vec> DeterministicEquivalent::solve(double timeLimit)
{
    d_model.set(GRB_DoubleParam_TimeLimit, timeLimit);
    d_model.optimize();

    auto status = d_model.get(GRB_IntAttr_Status);
    auto nSolutions = d_model.get(GRB_IntAttr_SolCount);

    // Time limit reached, and no feasible solution is available.
    if (status == GRB_TIME_LIMIT && nSolutions == 0)
        throw std::runtime_error("Time limit exceeded; no feasible solutions.");

    if (isOptimal() || nSolutions > 0)
    {
        auto const *vars = d_model.getVars();
        auto const numVars = d_problem.firstStageCoeffs().n_elem;

        auto const *xPtr = d_model.get(GRB_DoubleAttr_X, vars, numVars);
        auto result = std::make_unique<arma::vec>(xPtr, numVars);

        delete[] vars;
        delete[] xPtr;

        return result;
    }

    throw std::runtime_error("Deterministic equivalent is infeasible.");
}

double DeterministicEquivalent::firstStageObjective()
{
    auto const *vars = d_model.getVars();
    auto const numVars = d_problem.firstStageCoeffs().n_elem;

    auto const *xPtr = d_model.get(GRB_DoubleAttr_X, vars, numVars);
    arma::vec x(xPtr, numVars);

    delete[] vars;
    delete[] xPtr;

    return arma::dot(d_problem.firstStageCoeffs(), x);
}

double DeterministicEquivalent::secondStageObjective()
{
    return objective() - firstStageObjective();
}

bool DeterministicEquivalent::isOptimal() const
{
    return d_model.get(GRB_IntAttr_Status) == GRB_OPTIMAL;
}

double DeterministicEquivalent::mipGap() const
{
    return 100 * d_model.get(GRB_DoubleAttr_MIPGap);
}
