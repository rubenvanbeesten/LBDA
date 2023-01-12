#include "masterproblem.h"

MasterProblem::MasterProblem(ProblemData const &problem,
                             RiskMeasure const &risk_measure,
                             double lowerBound,
                             double upperBound) :
    d_problem(problem), d_model(d_env), d_risk_measure(risk_measure)
{
    d_model.addVar(lowerBound, upperBound, 1.0, GRB_CONTINUOUS, "theta");

    auto &Amat = d_problem.Amat();

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

    // write the model
    d_model.write("lbda_master_init_new.lp");
}

void MasterProblem::addCut(CutFamily::Cut const &cut)
{
    GRBVar const *vars = d_model.getVars();
    arma::vec const coeffs = -cut.beta;

    GRBLinExpr lhs;  // first variable is theta, and then all the x's.
    lhs.addTerms(coeffs.memptr(), vars + 1, coeffs.n_elem);
    lhs += vars[0];

    delete[] vars;

    d_model.addConstr(lhs, GRB_GREATER_EQUAL, cut.gamma);
}

std::unique_ptr<arma::vec> MasterProblem::solveWith(CutFamily &cutFamily,
                                                    double tol)
{
    while (true)
    {
        d_model.optimize();

        if (d_model.get(GRB_IntAttr_Status) != GRB_OPTIMAL)
            throw std::runtime_error("Master problem is infeasible.");

        size_t const numVars = d_problem.firstStageCoeffs().n_elem;

        auto *vars = d_model.getVars();

        auto *xValsPtr = d_model.get(GRB_DoubleAttr_X, vars + 1, numVars);
        double theta = vars[0].get(GRB_DoubleAttr_X);

        arma::vec xVals(xValsPtr, numVars);

        delete[] vars;
        delete[] xValsPtr;

        auto cut = cutFamily.computeCut(xVals);

        // Is the proposed cut violated by the current solution?
        if (cut.gamma + arma::dot(xVals, cut.beta) >= theta + tol)
            addCut(cut);
        else
            return std::make_unique<arma::vec>(xVals);
    }
}

double MasterProblem::firstStageObjective() const
{
    return objective() - secondStageObjective();
}

double MasterProblem::secondStageObjective() const
{
    // First variable is theta, and we need its value.
    return d_model.getVar(0).get(GRB_DoubleAttr_X);
};


std::vector<std::string> MasterProblem::getVarNames()
{
    auto variables = d_model.getVars();
    auto n_variables = d_model.get(GRB_IntAttr_NumVars);
    std::vector<std::string> varNames;
    for(int i = 0; i < n_variables; i++){
        std::string curName = variables[i].get(GRB_StringAttr_VarName);
        if(curName.at(0) == 'x' || curName.at(0) == 'u' || curName == "theta"){
            varNames.push_back(curName);
        }
    }

    return varNames;
}