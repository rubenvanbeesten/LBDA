#include "cutfamilies/cutfamily.h"


CutFamily::CutFamily(ProblemData const &problem) :
    d_model(d_env), d_problem(problem), d_sub(problem)
{
    GRBModel *m = &d_model;
    m->set("OutputFlag", "0");
}
