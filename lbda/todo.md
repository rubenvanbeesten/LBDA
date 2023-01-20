TODO:

- Check if we write all relevant information to file
    x Computation time
    x Objective value
    x MIP gap (for DEQ)
    - optimality gap (HOW?)

- plan experiments
    - make a list of risk measures
        - expectation, pure CVaR, mean-CVaR, integral
    - make a list of problems to use
        - all sslp instances

x automate running of problems
    x .sh file?





- do artificial parallel computing? (Niels didn't do so either)



DONE:

x check whether scenarioProbability are all equal (otherwise our sorting stuff doesn't work)
x read a risk measure from the command line or from a file (DONE)
x put the new files on the right place in the camkelists files (DONE)
x also make changes to the deterministic equivalent (or anything else we want to compare the LBDA with)
x find all places where objective function is implicitly computed using expected value
    x note: for loosebenders only used in subproblem, not in master
x write solution to file (and to screen)
    x measure solution time
    x change file name based on method
x figure out a way to set a lower bound on u_i
    x for now: just use 10e6
x test performance of lbda (does it work?)
    x compare new lbda with risk-neutral implementation (set risk measure to expectation)
    x remaining problem: mogelijk gaat er iets mis met het outputten van de oplossing (heeft te maken met dat xVars[0] = theta)
        x opgelost. probleem was dat in getVarNames ook theta mee werd genomen, terwijl in decisions alleen de x-variabelen zaten

