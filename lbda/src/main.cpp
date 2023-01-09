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

    
    /*
    // check the solution method
    switch (arguments.methodType)
    {
        // decomposition method
        case Arguments::DECOMPOSITION:
        {
            // define master problem and empty CutFamliy pointer
            MasterProblem master(problem, arguments.lb, arguments.ub);
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
                                                 alpha,
                                                 arguments.timeLimit);
                    break;
            }

            // solve problem
            auto solution = master.solveWith(*cutFamily);
            auto decisions = *solution;

            delete cutFamily;

            // print solution
            printSolution(decisions, master);
            break;
        }

        // DEF
        case Arguments::DETERMINISTIC_EQUIVALENT:
        {
            DeterministicEquivalent deq(problem);

            auto solution = deq.solve(arguments.timeLimit);
            auto decisions = *solution;

            printSolution(decisions, deq);

            std::cout << "Gap (%) = " << deq.mipGap() << "%\n";
            break;
        }
    }
    */
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


template<class T> void printSolution(arma::vec const &decisions, T &method)
{
    std::cout << "x = \n" << decisions << '\n';
    std::cout << "Obj.    = " << method.objective() << '\n';
    std::cout << "c^T x   = " << method.firstStageObjective() << '\n';
    std::cout << "Q(x)    = " << method.secondStageObjective() << '\n';
}
