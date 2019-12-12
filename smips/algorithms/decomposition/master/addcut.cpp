#include "master.h"

bool Master::addCut(
    double *beta, double gamma, double *x, double theta, double tol)
{
    double betaxgamma = gamma;
    for (size_t var = 0; var != d_n1; ++var)
        betaxgamma += beta[var] * x[var];
    // cout << "gap: " << betaxgamma  - theta << '\n';

    if (betaxgamma > theta + tol)  // then add cut and return false
    {
        ++d_nSlacks;
        /*
            // theta >= beta^T x + gamma <==> -theta + beta^T x <=  - gamma
        GRBLinExpr lhs = GRBLinExpr(d_theta, -1.0); // expression with
        coefficient = -1.0 for theta lhs.addTerms(beta, d_xVars, d_n1); //
        expression is now of the form -theta + beta^T x d_model.addConstr(lhs,
        GRB_LESS_EQUAL, -gamma);
        */

        // adding the cut to d_cmodel
        GRBaddvar(d_cmodel,
                  0,
                  NULL,
                  NULL,
                  0,
                  0,
                  1e20,
                  GRB_CONTINUOUS,
                  NULL);  // slack

        size_t slackIdx = d_n1
                          + d_nSlacks;  // slack variable index (there are d_n1
                                        // + 1 + nSlacks variables in the model)

        int cind[d_n1 + 2];
        iota(cind, cind + d_n1 + 1, 0);


        cind[d_n1 + 1]
            = slackIdx;  // refers to the latest variable (i.e. the slack)
        double cval[d_n1 + 1];
        for (size_t var = 0; var != d_n1; ++var)
            cval[var + 1] = -beta[var];
        cval[0] = 1;
        cval[d_n1 + 1] = -1;  // >= constraint, so slack features with -1

        GRBaddconstr(d_cmodel, d_n1 + 2, cind, cval, GRB_EQUAL, gamma, NULL);

        // add cut to internal storage of master
        d_xcoefs.push_back(vector<double>(beta, beta + d_n1));
        d_cons.push_back(gamma);

        d_kappa.push_back(1);
        d_beta.push_back(vector<double>{beta, beta + d_n1});
        d_gamma.push_back(gamma);

        return false;
    }

    return true;  // betaxgamma >= theta, no cut added
}