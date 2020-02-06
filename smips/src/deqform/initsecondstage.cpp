#include "deqform.h"

void DeqForm::initSecondStage(size_t n1,
                              size_t n2,
                              size_t p2,
                              size_t m2,
                              size_t S,
                              size_t ss_leq,
                              size_t ss_geq,
                              double *lb,
                              double *ub,
                              double *probs,
                              double *q,
                              arma::mat &Tmat,
                              arma::mat &Wmat,
                              arma::mat &omega)
{
    GRBLinExpr Tx[m2];
    for (size_t conIdx = 0; conIdx != m2; ++conIdx)
    {
        double *row = Tmat.colptr(conIdx);
        Tx[conIdx].addTerms(row, d_xVars, n1);
    }

    // variable types
    char vTypes2[n2];
    std::fill_n(vTypes2, p2, GRB_INTEGER);
    std::fill_n(vTypes2 + p2, n2 - p2, GRB_CONTINUOUS);

    // constraint senses
    char senses2[m2];
    std::fill(senses2, senses2 + ss_leq, GRB_LESS_EQUAL);
    std::fill(senses2 + ss_leq, senses2 + ss_leq + ss_geq, GRB_GREATER_EQUAL);
    std::fill(senses2 + ss_leq + ss_geq, senses2 + m2, GRB_EQUAL);

    // for each scenario: add variables and constraints
    for (size_t s = 0; s != S; ++s)
    {
        double const prob = probs[s];
        double costs[n2];

        for (size_t var = 0; var != n2; ++var)
            costs[var] = prob * q[var];

        double *rhsOmega = omega.colptr(s);

        // add variables
        GRBVar *yVars = d_model.addVars(lb, ub, costs, vTypes2, nullptr, n2);

        // lhs expression of second-stage constraints, including Wy
        // TODO: I've removed a duplication for TxWy here, which did not compile
        for (size_t conIdx = 0; conIdx != m2; ++conIdx)
        {
            double *row = Wmat.colptr(conIdx);
            Tx[conIdx].addTerms(row, yVars, n2);
        }

        GRBConstr *constrs = d_model.addConstrs(Tx,
                                                senses2,
                                                rhsOmega,
                                                nullptr,
                                                m2);

        delete[] yVars;
        delete[] constrs;
    }
}