# Test Suite
This is the first idea of trying to programm a test suite where we can check the robustness of a pFBA formulation. We use a standard dFBA process: for now this is my D-lactic acid producing fed-batch process.

We then change a couple of hyperparameters and run the IpOpt optimizer and see
* how long the optimization takes,
* if a local optimum is found,
* what the optimum value is,
* and ??? (I forgot).

With this information we can then try to rank the implementations depending on their accuracy, but also speed and usefulness while changing the bioprocess.


## Analysis

```check_simulations.ipynb``` is a notbook for analysis of the simulatons.

## List of Different Implementations

* version_01 ... pFBA every Collocation Point + Complementary Slackness in the Objective
