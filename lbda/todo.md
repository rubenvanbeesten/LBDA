TODO:

- save results in .csv file (like the excel file I'm maintaining)
- suppress Gurobi output (not yet working)

- NOTE: lbda is still very slow. why?
    - creating/destroying subproblems?
        - no
    - printing output?



TODO TO TEST ISSUE:

x run original lbda, compare results with 1.0_0.0 in my implementation
    x for sslp_10_50_50, result is the same
    x so probably, issue is in my "fixed" implementation
        - wait a second: LBDA in original implementation also does not find the optimal solution for 

- run "fixed" with deq solution. should give same result

- fix the -0.0 issue. We can do rounding, since first stage variables are integer anyway

- lower CVaR lower bound?
- extend subproblem time to 10 sec?
- run on a problem from the test set. these have been checked better I suppose
    - One that's in the test set: sslp_5_25_50

- also store 2nd stage solutions if I can't find the issue

- CHANGE LOWER BOUND -l IN RUNS!!!!!!!!
    - find a smart way to set the lower bound
    - probably use the deq solution for this purpose. this should speed up computation a bit


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

x check: do we make new subproblems every time or do we update the RHS?
    x answer: we update the RHS. See file loosebenders.cpp

x Check if we write all relevant information to file
    x Computation time
    x Objective value
    x MIP gap (for DEQ)
    x optimality gap (store best lower bound for DEQ, compute opt gap for lbda later)

x automate running of problems
    x .sh file?

x do artificial parallel computing? No.

x plan experiments
    x make a list of risk measures
        x expectation, pure CVaR, mean-CVaR, integral
    x make a list of problems to use
        x all sslp instances

