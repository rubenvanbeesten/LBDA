#ifndef SCENARIOPROBLEM_H
#define SCENARIOPROBLEM_H

#include "problemdata.h"

#include <armadillo>
#include <gurobi_c++.h>
#include <iosfwd>
#include <memory>
#include <cmath>

// new:
#include "risk_measure.h"

/**
 * Second-stage scenario problem.
 */
class ScenarioProblem
{
    ProblemData const &d_problem;

    GRBEnv d_env = GRBEnv();
    GRBModel d_model;

    void initFirstStage();
    void initSecondStage();

public:
    ScenarioProblem(ProblemData const &problem);

    /**
     * Solves the scenario problem.
     *
     * @param timeLimit Time limit for the Gurobi solver, in seconds.
     * @return objective value
     */
    double solve(double timeLimit = arma::datum::inf);

    /**
     * @return True if the model was solved to optimality, false otherwise.
     */
    bool isOptimal() const;

    /**
     * @return Returns the optimality gap of the solution returned by
     *         <code>solve</code>.
     */
    //double mipGap() const;

    /**
     * @return First-stage objective value.
     */
    //double firstStageObjective();

    /*
     * @return expected cost-to-go of the second-stage.
     */
    //double secondStageObjective();

    /**
     * @return Objective value.
     */
    double objective() const
    {
        return d_model.get(GRB_DoubleAttr_ObjVal);
    }

    /**
     * @return variable names
     */
    std::vector<std::string> getVarNames();

    /**
     * Fix the first stage variables
     */
    //void fixFirstStageVariables(std::vector<std::string> fixedVarNames, std::vector<double> fixedVarValues);

    /**
     * Get the best lower bound
     */
    //double getBestBound(void);
    
};

#endif
