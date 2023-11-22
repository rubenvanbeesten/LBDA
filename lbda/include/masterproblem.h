#ifndef MASTERPROBLEM_H
#define MASTERPROBLEM_H

#include "cutfamilies/cutfamily.h"
#include "problemdata.h"

#include <armadillo>
#include <gurobi_c++.h>
#include <iosfwd>
#include <memory>

#include "risk_measure.h"

/**
 * First-stage (master) problem of the two-stage decomposition.
 */
class MasterProblem
{
    GRBEnv d_env = GRBEnv();

    ProblemData const &d_problem;
    GRBModel d_model;

    RiskMeasure const &d_risk_measure;

    /**
     * Adds cut <code>theta >= beta^T x + gamma</code>.
     */
    void addCut(CutFamily::Cut const &cut);

public:
    /**
     * Constructs the master problem from the data in the passed-in problem
     * instance.
     *
     * @param problem       ProblemData instance.
     * @param lowerBound    Lower bound for theta, the approximation of the
     *                      expected cost-to-go. Default 0.
     * @param upperBound    Upper bound for theta. Default +inf.
     */
    MasterProblem(ProblemData const &problem,
                  RiskMeasure const &risk_measure,
                  double lowerBound = 0.,
                  double upperBound = arma::datum::inf
                  );

    /**
     * Solves the master problem using the given cutting strategy. The
     * master problem is solved s.t. the optimality gap is smaller than tol.
     *
     * @param cutFamily The cutting strategy (family) to use.
     * @param tol       Maximum acceptable optimality gap.
     *
     * @return Vector of (near) optimal first-stage decisions.
     */
    std::unique_ptr<arma::vec> solveWith(CutFamily &cutFamily,
                                         double timeLimit = arma::datum::inf,
                                         double tol = 1e-3
                                         );

    /**
     * Returns the (exact) first stage objective value, that is, the value of
     * <code>c^T x*</code>.
     */
    double firstStageObjective() const;

    /**
     * Returns the approximation of the expected costs of the second stage
     * problems, that is, the value of <code>theta</code> in the master problem.
     */
    double secondStageObjective() const;

    /**
     * Returns the objective value of the master problem, that is, the sum
     * of the first and second stage objectives.
     */
    double objective() const
    {
        return d_model.get(GRB_DoubleAttr_ObjVal);
    }

    /**
     * Returns the variable names of the master problem
     */
    std::vector<std::string> getVarNames();
    
    
    // not used:
    double getBestBound(void);

};

#endif
