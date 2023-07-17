# Test Suite
This is the first idea of trying to programm a test suite where we can check the robustness of a dFBA formulation. We use a standard dFBA process: for now this is my D-lactic acid producing fed-batch process.

We then change a couple of hyperparameters and run the IpOpt optimizer and see
* how long the optimization takes,
* if a local optimum is found,
* what the optimum value is,
* and the size of the complementary slackness.

With this information we can then try to rank the implementations depending on their accuracy, but also speed and usefulness while changing the bioprocess.


## Analysis

```check_simulations.ipynb``` is a notbook for analysis of the simulatons.

## List of Different Implementations

| Version ID | Description |
|------------|-------------|
| version_01 | pFBA at every CP, CS in the OF
| version_02 | FBA at every CP, CS in the OF
| version_03 | pFBA at every FE, CS in the OF
| version_04 | FBA at every FE, CS in the OF
| version_05 | cp of v01, no lb on N_G


# List of Abbreviations

| Abbreviation | Meaning |
|---|---|
|CP | collocation point |
|CS | complementary slackness |
|OF | objective function |
|FE | finite element
|lb | lower bound
|N_G | mol of glucose in reactor
