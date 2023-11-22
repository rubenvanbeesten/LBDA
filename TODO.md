# To do for revision

[Done] 1. Why are objectives of DEF and LBDA for CVaR different? [R2.4]
    x Answer: because they have different solutions. LBDA is an approximation.
    x We calculated the "out-of-sample" objective for the LBDA solution
        x By that I mean that we plugged the LBDA solution into the DEQ and computed the actual objective
    x So just explain this and no other actions necessary.


2. For instances not solved within time limit, report optimality gap [R2.3]
    x Rewrite code to store optimality gap
    - Run experiments again


3. Why do we have NA values for larger problems with spectral if they do solve within time limit? [E.4]
    - Answer: we use a DEQ with fixed first-stage variables (FIXED) to evaluate the solution quality of LBDA
        - DEQ is out of memory for the big instances
        - FIXED is also out of memory for these instances
        - So we cannot evaluate the quality of the LBDA solution
    - Question: how to deal with this?
        - We could use the objective function value found by LBDA itself
            - But this is misleading. It is often too optimistic
        - We could try to somehow run all the subproblems separately (the first-stage problem doesn't exist anyway)
            - This might be a lot of work, though
