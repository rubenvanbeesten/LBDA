#include "deterministicequivalent_new.h"


DeterministicEquivalent::DeterministicEquivalent(ProblemData const &problem, RiskMeasure const risk_measure) :
    d_problem(problem), d_model(d_env), d_risk_measure(risk_measure)
{
    initFirstStage();
    initSecondStage();
}

// deterministic equivalent formulation: 
//
// - old:   min_{x,                     c * x                                                           // first stage
//               y_s}                   + sum_s( q_s * y_s )                                            // second stage
//
//          s.t.                        x in X                                                          // first stage
//                                      W * y_s >= h_s                                                  // second stage
//                                      y_s >= 0                                                        // second stage
//
//              
// - new:   min_{x, u_i,                c * x + sum_i( lambda_i * u_i )                                 // first stage
//               y_s, tau_{i,s}}        + sum_s( sum_i( lambda_i * (1-beta_i)^(-1) * tau_{i,s} ) )      // second stage
//
//          s.t.                        x in X                                                          // first stage      (old)
//                                      W * y_s >= h_s                                                  // second stage     (old)
//                                      y_s >= 0                                                        // second stage     (old)
//                                      tau_{i,s} >= q * y_s - u_i                                      // second stage     (new)
//                                      tau_{i,s} >= 0                                                  // second stage     (new)


void DeterministicEquivalent::initFirstStage()
{   
    std::cout << "Initializing first stage..." << std::endl;

    // 1.A STANDARD VARIABLES AND CONSTRAINTS    
    //      variables:      x
    //      constraints:    x in X, i.e., A * x >= b

    // extract A matrix
    auto const &Amat = d_problem.Amat();

    // create names for x variabels
    std::string x_names[Amat.n_rows] = {""};
    for(size_t i=0; i < Amat.n_rows; i++){
        x_names[i] = "x_" + std::to_string(i);  // x_i
    }

    // create x variables based on d_problem
    auto *xVars = d_model.addVars(d_problem.firstStageLowerBound().memptr(),    // lower bound
                                  d_problem.firstStageUpperBound().memptr(),    // upper bound
                                  d_problem.firstStageCoeffs().memptr(),        // objective coefficients
                                  d_problem.firstStageVarTypes().memptr(),      // variable type
                                  x_names,                                      // variable names
                                  d_problem.firstStageCoeffs().n_elem);         // number of variables

    // initialize lhs expression
    GRBLinExpr lhs[Amat.n_cols];
    // fill lhs expression based on Amat
    for (auto iter = Amat.begin(); iter != Amat.end(); ++iter)
        lhs[iter.col()] += *iter * xVars[iter.row()];

    // add first stage constraints
    auto *constrs = d_model.addConstrs(lhs,
                                       d_problem.firstStageConstrSenses().memptr(),
                                       d_problem.firstStageRhs().memptr(),
                                       nullptr,
                                       Amat.n_cols);

    // memory management
    delete[] xVars;
    delete[] constrs;


    // 1.B RISK VARIABLES AND CONSTRAINTS      
    //      variables:      u_i
    //      constraints:    none

    // initialize information about u_i
    int n_cvars = d_risk_measure.lambdas.size();
    std::cout << "n_cvars: " << n_cvars << std::endl;
    double u_lb[n_cvars];
    double u_ub[n_cvars];
    double u_coef[n_cvars];
    char u_type[n_cvars];
    std::string u_names[n_cvars];
    for(int i = 0; i < n_cvars; i++){
        //u_lb[i] = -GRB_INFINITY;  // TEST
        u_lb[i] = - 10000;
        u_ub[i] = GRB_INFINITY;
        u_coef[i] = d_risk_measure.lambdas.at(i);
        u_type[i] = GRB_CONTINUOUS;
        u_names[i] = "u_" + std::to_string(i);
    }

    // add u_i to the model
    //auto *uVars = d_model.addVars(... // old: store in uVars 
    d_model.addVars(u_lb,     // lower bound
                    u_ub,     // upper bound
                    u_coef,   // objective coefficients
                    u_type,   // variable type
                    u_names,  // variable names
                    n_cvars); // number of variables
    
    // note: no constraints on u_i

    // memomry management
    //delete[] uVars; // not sure if we need this. probably a memory thing (old)

    // update model
    d_model.update();

    std::cout << "Finished initializing first stage." << std::endl;
}

void DeterministicEquivalent::initSecondStage()
{   
    std::cout << "Initializing second stage..." << std::endl;

    // extract T and W matrices
    auto const &Tmat = d_problem.Tmat();
    auto const &Wmat = d_problem.Wmat();

    // extract number of first-stage *x* variables (as opposed to u variables)
    int n_x_vars = d_problem.firstStageCoeffs().n_elem;

    // import x variables from model
    auto *xVars = d_model.getVars();    // NOTE: this also contains the u variables! Things still go right though, since we only use the first n elements of xVars below.

    // extract number of cvars
    int n_cvars = d_risk_measure.lambdas.size();

    // loop over scenarios
    for (size_t scenario = 0; scenario != d_problem.nScenarios(); ++scenario)
    {   
        // extract scenario probability
        double const prob = d_problem.scenarioProbability(scenario);

        // 2.A STANDARD VARIABLES AND CONSTRAINTS    
        //      variables:      y_s
        //      coefficients:   0                (NOTE: different from original)
        //      constraints:    W * y_s >= h_s   (and nonnegativity on y_s)

        // extract coefficients q
        //arma::vec const costs = prob * d_problem.secondStageCoeffs();   // prob[s] * q    (old)
        arma::vec const q_vec = d_problem.secondStageCoeffs();          // q

        // new coefficients for y in the subproblem (i.e., equal to zero)       // new
        arma::vec const y_coef = arma::zeros(Wmat.n_rows);                      // note: all matrices are transposed as compared to what I'm used to. so number of rows of W is length of vector y

        // initialize constraints lhs expression
        GRBLinExpr lhs_y[Wmat.n_cols];              // Wmat.n_cols: number of constraints in Tx + Wy >= h

        // add Tx to lhs expression
        for (auto iter = Tmat.begin(); iter != Tmat.end(); ++iter){ // loop over all elements of the matrix
            lhs_y[iter.col()] += *iter * xVars[iter.row()];         // Tx            
            // NOTE: we only use the first Tmat.end() elements of xVars (corresp. to x, not u), so things should still work.
        }

        // create names for y variabels
        std::string y_names[Wmat.n_rows] = {""};
        for(size_t i=0; i < Wmat.n_rows; i++){
            y_names[i] = "y_" + std::to_string(scenario) + "_" + std::to_string(i);  // y_s_i
        }

        // scenario-specific variables
        auto *yVars = d_model.addVars(d_problem.secondStageLowerBound().memptr(),
                                      d_problem.secondStageUpperBound().memptr(),
                                      y_coef.memptr(),                              // old: costs.memptr(),
                                      d_problem.secondStageVarTypes().memptr(),
                                      y_names,
                                      Wmat.n_rows);

        // add Wy to lhs expression
        for (auto iter = Wmat.begin(); iter != Wmat.end(); ++iter)
            lhs_y[iter.col()] += *iter * yVars[iter.row()];  // Wy

        // extract constraint senses (>=,<=,==) from problem
        auto const &senses_y = d_problem.secondStageConstrSenses();
        // extract rhs from problem
        auto const rhs_y = d_problem.scenarios().colptr(scenario);

        // add constraint: Tx + Wy >= h     (note: possible inequality in other direction)
        auto *constrs_y = d_model.addConstrs(lhs_y,               // lhs
                                             senses_y.memptr(),   // direction of inequality
                                             rhs_y,               // rhs
                                             nullptr,             // name
                                             Wmat.n_cols);        // number of constraints



        // 2.B RISK VARIABLES AND CONSTRAINTS      
        //      variables:     tau_{i,s} 
        //      constraints:   tau_{i,s} >= q * y_s - u_i,  tau_{i,s} >= 0

        // information for tauVars
        double tau_lb[n_cvars];
        double tau_ub[n_cvars];
        double tau_coef[n_cvars];
        char tau_type[n_cvars];
        std::string tau_names[n_cvars];

        for(int i = 0; i < n_cvars; i++){
            tau_lb[i] = 0.0;
            tau_ub[i] = GRB_INFINITY;
            tau_coef[i] = prob * d_risk_measure.lambdas.at(i) * (1.0/(1.0 - d_risk_measure.betas.at(i)));
            tau_type[i] = GRB_CONTINUOUS;
            tau_names[i] = "tau_" + std::to_string(scenario) + "_" + std::to_string(i);  // tau_s_i
        }

        // create scenario-specific tau_{i,s} variables     (note: s is a constant within the curret loop)
        auto *tauVars = d_model.addVars(tau_lb,     // lower bound
                                        tau_ub,     // upper bound
                                        tau_coef,   // objective coefficients
                                        tau_type,   // variable type
                                        tau_names,    // variable names
                                        n_cvars);   // number of variables
        
        // initialize constraints lhs expression for constraint tau_{i,s} >= q * y_s - u_i  (for current fixed s)
        // reformulation:   tau_{i,s} - q * y_s + u_i >= 0
        GRBLinExpr lhs_tau[n_cvars];

        // fill expression for lhs
        for(int i = 0; i < n_cvars; i++){
            // we will make expression for constraint i

            // add tau_{s,i}
            lhs_tau[i] += tauVars[i];                   // + tau_{s,i}
            
            // add (- q * y_s)
            for(size_t j = 0; j < q_vec.size(); j++){
                lhs_tau[i] -= q_vec.at(j) * yVars[j];   // - q_j * y_{s,j}
            }

            // add u_i
            lhs_tau[i] += xVars[n_x_vars + i];          // + u_i
        }
        
        // define rhs expression
        double rhs_tau[n_cvars] = {0.0}; // rhs equal to zero

        // define inequality senses
        char senses_tau[n_cvars];
        for(int i = 0; i < n_cvars; i++){
            senses_tau[i] = GRB_GREATER_EQUAL;  // lhs >= rhs
        }
    
        // add constraint: tau_{i,s} - q * y_s + u_i >= 0
        auto *constrs_tau = d_model.addConstrs(lhs_tau,     // lhs
                                               senses_tau,  // direction of inequality
                                               rhs_tau,     // rhs
                                               nullptr,     // name
                                               n_cvars);    // number of constraints

        


        // memory management 
        delete[] yVars;
        delete[] constrs_y;
        delete[] tauVars;
        delete[] constrs_tau;
    }

    // memory management
    delete[] xVars;

    
    //////////////

    // update model
    d_model.update();

    // write the model
    d_model.write("deq_new.lp");

    std::cout << "Finished initializing second stage." << std::endl;
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
        auto const numVars = d_problem.firstStageCoeffs().n_elem + d_risk_measure.lambdas.size(); // also add u_i

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

std::vector<std::string> DeterministicEquivalent::getVarNames()
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