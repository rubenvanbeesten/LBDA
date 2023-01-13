TODO:

x check whether scenarioProbability are all equal (otherwise our sorting stuff doesn't work)
x read a risk measure from the command line or from a file (DONE)
x put the new files on the right place in the camkelists files (DONE)
x also make changes to the deterministic equivalent (or anything else we want to compare the LBDA with)
x find all places where objective function is implicitly computed using expected value
    x note: for loosebenders only used in subproblem, not in master
x write solution to file (and to screen)
    x measure solution time
    x change file name based on method
- figure out a way to set a lower bound on u_i
    - otherwise we get numerical issues
- test performance of lbda (does it work?)
    x compare new lbda with risk-neutral implementation (set risk measure to expectation)
    - remaining problem: mogelijk gaat er iets mis met het outputten van de oplossing (heeft te maken met dat xVars[0] = theta)
        - problem in sslp_5_25_100 -L 1.0 -B 0.0

