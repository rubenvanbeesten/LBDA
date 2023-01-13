#include "main.h"

int main(int argc, char **argv)
try
{
    // parse arguments from command line
    auto arguments = parseArguments(argc, argv);

    // if option -h or --help, then just print explanation of program
    if (arguments.printUsage)
    {
        std::cout << USAGE << '\n';
        return EXIT_SUCCESS;
    }

    // read data from SMPS file
    auto problem = ProblemData::fromSmps(arguments.file);

    // test whether all the probabilities are equal
    for (size_t i = 0; i < problem.nScenarios(); i++){
        if(abs(problem.scenarioProbability(i) - (1.0/problem.nScenarios())) > 0.001){
            throw std::runtime_error("Probabilities are not equal.");
        }
    }
    
    // store risk measure
    auto risk_measure = RiskMeasure(str_to_double_vec(arguments.lambdaString), str_to_double_vec(arguments.betaString)); // convert lambda and beta strings to vector<double>

    
    // check the solution method
    switch (arguments.methodType)
    {
        // decomposition method
        case Arguments::DECOMPOSITION:
        {
            // define master problem and empty CutFamliy pointer
            MasterProblem master(problem, risk_measure, arguments.lb, arguments.ub);
            CutFamily *cutFamily;

            // check cut type and define cutFamily
            switch (arguments.cutType)
            {
                case Arguments::LP_DUAL:
                    cutFamily = new LpDual(problem);
                    break;
                case Arguments::STRONG_BENDERS:
                    cutFamily = new StrongBenders(problem);
                    break;
                case Arguments::LOOSE_BENDERS:
                default:
                    // TODO how to set alpha from the command line?
                    arma::vec alpha = arma::zeros(problem.Wmat().n_cols);
                    cutFamily = new LooseBenders(problem,
                                                 risk_measure,
                                                 alpha,
                                                 arguments.timeLimit);
                    break;
            }

            // solve problem
            auto start_time = std::chrono::high_resolution_clock::now(); // start time
            auto solution = master.solveWith(*cutFamily);   // solve model
            auto decisions = *solution;                     // extract solution
            auto stop_time = std::chrono::high_resolution_clock::now();  // stop time            
            auto sol_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time);

            delete cutFamily;

            // print solution
            printSolution(decisions, master, sol_time);

            // write solution to file
            writeSolution(decisions, master, sol_time, arguments);
            
            break;
        }

        // DEF
        case Arguments::DETERMINISTIC_EQUIVALENT:
        {
            //DeterministicEquivalent deq(problem);     // old
            DeterministicEquivalent deq(problem, risk_measure);

            // solve
            auto start_time = std::chrono::high_resolution_clock::now(); // start time
            auto solution = deq.solve(arguments.timeLimit);
            auto decisions = *solution;
            auto stop_time = std::chrono::high_resolution_clock::now();  // stop time
            auto sol_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time);


            // print solution
            printSolution(decisions, deq, sol_time);
            // report optimality gap
            std::cout << std::endl << "Gap (%) = " << deq.mipGap() << "%\n";

            // write solution to file
            writeSolution(decisions, deq, sol_time, arguments);

            break;
        }
    }
}
catch (GRBException const &e)
{
    std::cerr << e.getMessage() << '\n';
    return EXIT_FAILURE;
}
catch (std::exception const &e)
{
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
}
catch (...)
{
    std::cerr << "Something went wrong.\n";
    return EXIT_FAILURE;
}

argument_t parseArguments(int argc, char **argv)
{
    argument_t arguments;
    int option;

    while ((option = getopt(argc, argv, "hm:c:l:u:t:L:B:")) != -1)
    {
        switch (option)
        {
            case 'h': // help
                // If this flag is set, it is all we will process - immediately
                // after, the program ends. No need to continue processing.
                arguments.printUsage = true;
                return arguments;

            case 'm': // select solution method (deq, decomp). Default: decomp
                if (std::strcmp(optarg, "deq") == 0)
                    arguments.methodType = arguments.DETERMINISTIC_EQUIVALENT;

                break;

            case 'c': // select cut family (lbda, lp, sb). Default: lbda
                if (std::strcmp(optarg, "lp") == 0)
                    arguments.cutType = arguments.LP_DUAL;

                if (std::strcmp(optarg, "sb") == 0)
                    arguments.cutType = arguments.STRONG_BENDERS;

                break;

            case 'l': // set lower bound. Default: 0
                arguments.lb = std::stod(optarg);
                break;

            case 'u': // set upper bound. Default: +inf
                arguments.ub = std::stod(optarg);
                break;

            case 't': // set time limit. Default: none
                arguments.timeLimit = std::stod(optarg);
                break;
            
            case 'L': // set lambdas vector for risk measure. Default: {1.0}
                arguments.lambdaString = optarg;
                break;

            case 'B': // set betas vector for risk measure. Default: {0.0}
                arguments.betaString = optarg;
                break;

            default:
                break;  // unexpected option, but not the end of the world.
        }
    }

    if (argc == optind)
        throw std::runtime_error("Did not receive an SMPS file location.");

    arguments.file = argv[optind];
    return arguments;
}


template<class T> void printSolution(arma::vec const &decisions, T &method, std::chrono::milliseconds sol_time)
{
    // print solution
    std::cout << "First-stage decisions:" << std::endl;
    std::vector<std::string> varNames = method.getVarNames(); // first-stage variable names
    for(size_t i=0; i < varNames.size(); i++){
        std::cout << varNames.at(i) << " = " << decisions.at(i) << std::endl;
    }
    std::cout << std::endl;

    // print objective information
    std::cout << "Objective information:" << std::endl;
    std::cout << "Obj.    = " << method.objective() << '\n';            // objective value
    std::cout << "c^T x   = " << method.firstStageObjective() << '\n';  // first stage objective value
    std::cout << "Q(x)    = " << method.secondStageObjective() << '\n'; // second-stage objective value
    std::cout << std::endl;
    
    // print solution time
    std::cout << "Solution time: " << (sol_time.count()/1000.0) << " seconds" << std::endl;
}

template<class T> void writeSolution(arma::vec const &decisions, T &method, std::chrono::milliseconds sol_time, argument_t arguments)
{   
    // make output file name
    std::string problem_path = arguments.file;
    std::string problem_name = problem_path.substr(problem_path.find_last_of("/") + 1);
    std::string method_name;
    if(arguments.methodType == arguments.DETERMINISTIC_EQUIVALENT){
        method_name = "deq";
    }else{
        method_name = "decom";
    }

    mkdir("solution_reports", 0777); // make directory
    std::string outfilename = "solution_reports/" + problem_name + "_" + arguments.lambdaString + "_" + arguments.betaString + "_" + method_name + ".txt";

    // create/open a file with the right name
    std::ofstream outfile;
    outfile.open(outfilename);
    
    // write problem information
    outfile << "Problem: " << problem_name << std::endl << std::endl;

    // write solution method information
    outfile << "Solved with: " << method_name << std::endl << std::endl;

     // write solution time
    outfile << "Solution time: " << (sol_time.count()/1000.0) << " seconds" << std::endl << std::endl;

    // write objective information
    outfile << "Objective information:" << std::endl;
    outfile << "Obj.    = " << method.objective() << '\n';            // objective value
    outfile << "c^T x   = " << method.firstStageObjective() << '\n';  // first stage objective value
    outfile << "Q(x)    = " << method.secondStageObjective() << '\n'; // second-stage objective value
    outfile << std::endl;
   
    // write solution
    outfile << "First-stage decisions:" << std::endl;
    std::vector<std::string> varNames = method.getVarNames(); // first-stage variable names
    for(size_t i=0; i < varNames.size(); i++){
        outfile << varNames.at(i) << " = " << decisions.at(i) << std::endl;
    }
    outfile << std::endl;

    std::cin.get();

    // close the file
    outfile.close();    
}