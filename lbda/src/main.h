#ifndef MAIN_H
#define MAIN_H

#include "cutfamilies/loosebenders.h"
#include "cutfamilies/lpdual.h"
#include "cutfamilies/strongbenders.h"
#include "deterministicequivalent.h"
#include "masterproblem.h"

#include <cstring>
#include <gurobi_c++.h>
#include <iostream>
#include <stdexcept>
#include <unistd.h>

// new
#include <vector>
#include <chrono>
#include <fstream>
#include <string>
#include <sstream>
#include<sys/stat.h>    // mkdir
#include<sys/types.h>   // mkdir
#include "risk_measure.h"


auto const USAGE = R"(
LBDA+. A program for solving two-stage mixed-integer stochastic programs.

Usage:
    lbdaexec [-h] [-m method] [-c cut -l lb -u ub | -t time] <file>

Options:
    -h  Prints this help text.
    -m  Selects solution method. One of {deq, decomp}:
        * "deq". Solves the deterministic equivalent (extensive form).
        * "decomp". Solves a decomposition of the problem into a first- and
          second-stage problem. Default.
    -c  Cut family to use when solving the decomposition. One of {lbda, lp, sb}:
        * "lbda". Loose Benders approximation cuts. Default.
        * "lp". Uses cuts based on the LP dual of the second-stage problem.
        * "sb". Uses strengthened Benders' cuts, derived from the Lagrangian
          relaxation of the second-stage problem.
    -l  Lower bound on the expected cost-to-go in the second stage. Default 0.
    -u  Upper bound on the expected cost-to-go in the second stage. Default +inf.
    -t  Time limit (in seconds) to set when solving the deterministic equivalent,
        or time limit for the Gurobi relaxations when solving the decomposition
        using Loose Benders cuts.

Arguments:
    <file>  Location of the SMPS file triplet to solve. Should not contain any
            extensions.
)";

/**
 * Simple struct that gathers the command-line arguments.
 */
struct Arguments
{
    enum CutType
    {
        LP_DUAL,
        LOOSE_BENDERS,
        STRONG_BENDERS
    };

    enum MethodType
    {
        DETERMINISTIC_EQUIVALENT,
        DECOMPOSITION,
        FIXED_X
    };

    MethodType methodType = DECOMPOSITION;  // solution method to use
    CutType cutType = LOOSE_BENDERS;        // cut family to use
    double timeLimit = arma::datum::inf;    // max. solve time in seconds
    double subTimeLimit = arma::datum::inf; // max. subproblem solve time in seconds
    bool printUsage = false;                // print help text?
    double lb = 0;                          // lower bound on theta
    double ub = arma::datum::inf;           // upper bound on theta
    std::string file;                       // smps file location

    // new
    std::string lambdaString = "1.0";       // vector of lambdas for risk measure (initialized at expectation settings)
    std::string betaString = "0.0";         // vector of betas for risk measure (initialized at expectation settings)
};

using argument_t = struct Arguments;

/**
 * Parses command-line arguments. Heavily tied into the Posix/Unix way of doing
 * this, using <code>getopt</code>.
 *
 * @throws std::runtime_error when the arguments could not correctly be parsed.
 *
 * @param argc  Number of command-line arguments.
 * @param argv  Command-line arguments.
 * @return      Struct with all parsed arguments.
 */
argument_t parseArguments(int argc, char **argv);

/**
 * Prints the (near) optimal first-stage decisions passed in, and some objective
 * information derived from the method.
 *
 * @param decisions Near optimal first-stage decisions, as a vector.
 * @param method    Solution method used - the deterministic equivalent, or the
 *                  first-stage master problem.
 */
template<class T> void printSolution(arma::vec const &decisions, T &method, std::chrono::milliseconds sol_time);

/**
 * Writes the (near) optimal first-stage decisions passed in, and some objective
 * information derived from the method to a file.
 *
 * @param decisions Near optimal first-stage decisions, as a vector.
 * @param method    Solution method used - the deterministic equivalent, or the
 *                  first-stage master problem.
 */
template<class T> void writeSolutionReport(arma::vec const &decisions, T &method, std::chrono::milliseconds sol_time, argument_t arguments, double mipGap = -1.0);

/**
 * Writes the first-stage decisions to a csv file
 *
 * @param decisions Near optimal first-stage decisions, as a vector.
 * @param method    Solution method used - the deterministic equivalent, or the
 *                  first-stage master problem.
 */
template<class T> void writeFirstStageSolution(arma::vec const &decisions, T &method, argument_t arguments);


void parseFirstStageSolution(std::string solFileName, std::vector<std::string> &fixedVarNames, std::vector<double> &fixedVarValues);
//void parseFirstStageSolution(std::string solFileName);

#endif  // MAIN_H
