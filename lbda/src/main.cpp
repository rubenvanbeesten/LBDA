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
    auto risk_measure = RiskMeasure(arguments.lambdaString, arguments.betaString, problem.nScenarios()); 

    
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
                                                 arguments.subTimeLimit); // use subTimeLimit here
                    break;
            }

            // solve problem
            std::cout << "Starting LBDA..." << std::endl;
            auto start_time = std::chrono::high_resolution_clock::now(); // start time
            auto solution = master.solveWith(*cutFamily, arguments.timeLimit);   // solve model
            auto decisions = *solution;                     // extract solution
            auto stop_time = std::chrono::high_resolution_clock::now();  // stop time            
            auto sol_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time);
            std::cout << "LBDA finished." << std::endl;

            delete cutFamily;

            double mipGap = -1.0;           // note: mipGap = -1.0 used because it is irrelevant

            // print solution
            printSolution(decisions, master, sol_time, mipGap);

            // write solution report to file
            writeSolutionReport(decisions, master, sol_time, arguments);

            // write solution to results table
            writeResultsTable(master, sol_time, arguments, mipGap);   

            // write first stage solution to csv file
            writeFirstStageSolution(decisions, master, arguments);
            
            break;
        }

        // DEF
        case Arguments::DETERMINISTIC_EQUIVALENT:
        {
            DeterministicEquivalent deq(problem, risk_measure);

            // solve
            auto start_time = std::chrono::high_resolution_clock::now(); // start time
            auto solution = deq.solve(arguments.timeLimit);
            auto decisions = *solution;
            auto stop_time = std::chrono::high_resolution_clock::now();  // stop time
            auto sol_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time);


            // print solution
            printSolution(decisions, deq, sol_time, deq.mipGap());
            // report optimality gap
            std::cout << std::endl << "Gap (%) = " << deq.mipGap() << "%\n";

            // write solution report to file
            writeSolutionReport(decisions, deq, sol_time, arguments, deq.mipGap());

            // write solution to results table
            writeResultsTable(deq, sol_time, arguments, deq.mipGap());

            // write first stage solution to csv file
            writeFirstStageSolution(decisions, deq, arguments);

            break;
        }
        // FIXED x 
        case Arguments::FIXED_X:
        {   
            // Initialize DEf problem
            DeterministicEquivalent deq(problem, risk_measure);
            
            std::string problemPath = arguments.file;
            std::string problemName = problemPath.substr(problemPath.find_last_of("/") + 1);


            std::string solFileName = "first_stage_solution/" + problemName + "_" + arguments.lambdaString;
            if(arguments.lambdaString != "integral"){
                solFileName += "_" + arguments.betaString;
            }
            solFileName += "_decom.csv";

            std::cout << "solFileName = " << solFileName << std::endl;
            
            //std::string solFileName = "first_stage_solution/sizes3_1.0_0.0_decom.csv";  // TEMPORARY
            //std::string solFileName = "first_stage_solution/sizes3_1.0_0.0_deq.csv";  // TEMPORARY
            std::vector<std::string> fixedVarNames;     // initialize vector of variable names
            std::vector<double> fixedVarValues;         // initialize vector of variable values
            parseFirstStageSolution(solFileName, fixedVarNames, fixedVarValues);     // parse the data from the csv into the vectors

            //check (temp)
            std::cout << "Fixed solutions: " << std::endl;
            for(size_t i=0; i < fixedVarNames.size(); i++){
                std::cout << fixedVarNames[i] << " = " << fixedVarValues[i] << std::endl;
            }

            std::cout << "Fixing first stage variables...";
            deq.fixFirstStageVariables(fixedVarNames, fixedVarValues);
            std::cout << " done." << std::endl;

            // solve
            auto start_time = std::chrono::high_resolution_clock::now(); // start time
            auto solution = deq.solve(arguments.timeLimit);
            auto decisions = *solution;
            auto stop_time = std::chrono::high_resolution_clock::now();  // stop time
            auto sol_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time);


            // print solution
            printSolution(decisions, deq, sol_time, deq.mipGap());

            // report optimality gap
            std::cout << std::endl << "Gap (%) = " << deq.mipGap() << "%\n";

            // write solution report to file
            writeSolutionReport(decisions, deq, sol_time, arguments, deq.mipGap());

            // write solution to results table
            writeResultsTable(deq, sol_time, arguments, deq.mipGap());

            // write first stage solution to csv file
            writeFirstStageSolution(decisions, deq, arguments);

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

    while ((option = getopt(argc, argv, "hm:c:l:u:t:T:L:B:")) != -1)
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
                if (std::strcmp(optarg, "fixed") == 0)
                    arguments.methodType = arguments.FIXED_X;
                if (std::strcmp(optarg, "decom") == 0)
                    arguments.methodType = arguments.DECOMPOSITION;
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

            case 't': // set subproblem time limit. Default: none
                arguments.subTimeLimit = std::stod(optarg);
                break;
            
            case 'T': // set overall time limit. Default: none
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


template<class T> void printSolution(arma::vec const &decisions, T &method, std::chrono::milliseconds sol_time, double mipGap)
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

    if(mipGap != -1.0){
        std::cout << "MIP Gap: " << mipGap << "%" << std::endl << std::endl;
    }
    
    // print solution time
    std::cout << "Solution time: " << (sol_time.count()/1000.0) << " seconds" << std::endl;
}

template<class T> void writeSolutionReport(arma::vec const &decisions, T &method, std::chrono::milliseconds sol_time, argument_t arguments, double mipGap)
{   
    // make output file name
    std::string problem_path = arguments.file;
    std::string problem_name = problem_path.substr(problem_path.find_last_of("/") + 1);
    std::string method_name;
    if(arguments.methodType == arguments.DETERMINISTIC_EQUIVALENT){
        method_name = "deq";
    }else if(arguments.methodType == arguments.DECOMPOSITION){
        method_name = "decom";
    }else if(arguments.methodType == arguments.FIXED_X){
        method_name = "fixed";
    }

    mkdir("solution_reports", 0777); // make directory
    std::string outfilename = "solution_reports/" + problem_name + "_" + arguments.lambdaString;
    if(arguments.lambdaString != "integral"){
        outfilename = outfilename + "_" + arguments.betaString;
    }
    outfilename = outfilename + "_" + method_name + ".txt";

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
    if(mipGap != -1.0){
        outfile << "MIP Gap: " << mipGap << "%" << std::endl << std::endl;
    }
    if(arguments.methodType == arguments.DETERMINISTIC_EQUIVALENT){
        outfile << "Best lower bound: " << method.getBestBound() << std::endl << std::endl;
    }
   
    // write solution
    outfile << "First-stage decisions:" << std::endl;
    std::vector<std::string> varNames = method.getVarNames(); // first-stage variable names
    for(size_t i=0; i < varNames.size(); i++){
        outfile << varNames.at(i) << " = " << decisions.at(i) << std::endl;
    }
    outfile << std::endl;

    // close the file
    outfile.close();   
}


template<class T> void writeResultsTable(T &method, std::chrono::milliseconds sol_time, argument_t arguments, double mipGap)
{   
    // 1. Create file
    // file name
    std::string outfilename = "results_table.csv";

    // create file object
    std::ofstream outfile;
    // check if file already exists
    std::ifstream f(outfilename.c_str());
    if(f.good()){
        // if file exists, open file in append mode
        outfile.open(outfilename, std::ios_base::app);    
    }else{
        // if file does not exist, create it
        outfile.open(outfilename);
        // also add first line
        outfile << "instance;risk_measure;method;objective;mip_gap;solution_time" << std::endl;
    }

    // 2. Process information
    // problem (instance) name
    std::string problem_path = arguments.file;
    std::string problem_name = problem_path.substr(problem_path.find_last_of("/") + 1);
    // method name
    std::string method_name;
    if(arguments.methodType == arguments.DETERMINISTIC_EQUIVALENT){
        method_name = "DEQ";
    }else if(arguments.methodType == arguments.DECOMPOSITION){
        method_name = "decom";
    }else if(arguments.methodType == arguments.FIXED_X){
        method_name = "fixed";
    }
    // risk measure
    std::string risk_measure_name;
    if(arguments.lambdaString == "1.0" && arguments.betaString == "0.0"){
        risk_measure_name = "expectation";
    }else if(arguments.lambdaString == "1.0") {
        risk_measure_name = "CVaR_" + arguments.betaString;
    }else if(arguments.lambdaString == "0.5,0.5"){
        risk_measure_name = "mean-CVaR";
    }else if(arguments.lambdaString == "integral"){
        risk_measure_name = "integral";
    }else{
        risk_measure_name = "rho_" + arguments.lambdaString + "_" + arguments.betaString;
    }


    
    // 3. Write current run information to file
    outfile << problem_name << ";" << risk_measure_name << ";" << method_name << ";" << method.objective() << ";" << mipGap << ";" << (sol_time.count()/1000.0) << std::endl;

    // close the file
    outfile.close();   
}


template<class T> void writeFirstStageSolution(arma::vec const &decisions, T &method, argument_t arguments)
{   
    // make output file name
    std::string problem_path = arguments.file;
    std::string problem_name = problem_path.substr(problem_path.find_last_of("/") + 1);
    std::string method_name;
    if(arguments.methodType == arguments.DETERMINISTIC_EQUIVALENT){
        method_name = "deq";
    }else if(arguments.methodType == arguments.DECOMPOSITION){
        method_name = "decom";
    }else if(arguments.methodType == arguments.FIXED_X){
        method_name = "fixed";
    }

    mkdir("first_stage_solution", 0777); // make directory
    std::string solfilename = "first_stage_solution/" + problem_name + "_" + arguments.lambdaString;
    if(arguments.lambdaString != "integral"){
        solfilename = solfilename + "_" + arguments.betaString;
    }
    solfilename = solfilename + "_" + method_name + ".csv";

    // create/open a file with the right name
    std::ofstream solfile;
    solfile.open(solfilename, std::ofstream::trunc);
    
    // write solution
    std::vector<std::string> varNames = method.getVarNames(); // first-stage variable names
    solfile << "variable;value" << std::endl; // header
    for(size_t i=0; i < varNames.size(); i++){
        solfile << varNames.at(i) << ";" << decisions.at(i) << std::endl; // entry
    }

    // close the file
    solfile.close();   
}

void parseFirstStageSolution(std::string solFileName, std::vector<std::string> &fixedVarNames, std::vector<double> &fixedVarValues)
//void parseFirstStageSolution(std::string solFileName)
{
    // open the file
    std::fstream file(solFileName, std::ios::in);

    std::vector<std::vector<std::string>> content;  // not sure what this does
    std::vector<std::string> row;                   // not sure what this does

    std::string line;   // this will hold the contents of a line
    std::string entry;  // this will hold the contents of an entry
    if(file.is_open()){ // open the file
        while(getline(file, line)){  // read a line and write it into "line"
            row.clear();

            std::stringstream str(line);    // interpret line as a stringstream called "str"

            while(getline(str, entry, ';')){    // get all entries separated by ";" in "str"
                row.push_back(entry);           // push back the entry to "row"
            }
            content.push_back(row);
        }
    }
    else{
        std::cout << "Could not open the file." << std::endl;
    }

    // store everything in two vectors
    fixedVarNames.clear();
    fixedVarValues.clear();
    for(size_t i=1; i < content.size(); i++){ // loop over rows (start at row number 1)
        if(content[i].size() == 2){ // should be two entries in the row
            if(content[i][0].at(0) == 'x'){ // check if variable name starts with "x"
                fixedVarNames.push_back(content[i][0]);
                fixedVarValues.push_back(abs(stod(content[i][1])));    
            }           
        }else{
            throw std::runtime_error("Number of entries is not equal to two.");
        }
    }
}