[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# Regularized MIP Model for Integrating Energy Storage Systems and its Application for Solving a Trilevel Interdiction Problem

This project is distributed in association with the [INFORMS Journal on Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The software and data in this repository are a snapshot of the software and data that were used in the research reported on in the paper [Regularized MIP Model for Integrating Energy Storage Systems and its Application for Solving a Trilevel Interdiction Problem](https://doi.org/10.1287/ijoc.2024.0771) by Dahye Han, Nan Jiang, Santanu S. Dey and Weijun Xie.

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2024.0771

https://doi.org/10.1287/ijoc.2024.0771.cd

Below is the BibTex for citing this snapshot of the repository.

```
@misc{han2025regularized,
  title =     {Regularized MIP Model for Integrating Energy Storage Systems and its Application for Solving a Trilevel Interdiction Problem},
  author= {Han, Dahye and Jiang, Nan and Dey, Santanu S and Xie, Weijun},
  publisher = {INFORMS Journal on Computing},
  year =      {2025},
  doi =       {10.1287/ijoc.2024.0771.cd},
  url =       {https://github.com/INFORMSJoC/2024.0771},
  note =      {Available for download at https://github.com/INFORMSJoC/2024.0771},
}
```

## Description

This software provides the implementation of different models introduced for solving trilevel interdiction problem described in the paper, including the regularized MIP model, LP relaxation, and vertex-representation model.

## Content of the repository

This repository includes:

- the Julia source code to solve the trilevel application problem (directory src)
- the Julia source code to generate figures (directory src)
- results of the computational experiments (directory results)
- figures in the paper (directory results)
- the data to run trilevel application problem including demand, network information and battery parameters (directory data)

## Replicating

### Prerequisites

The experiments in the repository are based on:

- [Julia](https://julialang.org) version 1.7

- [Gurobi](https://www.gurobi.com) version 10.0

Additionally, it uses the following Julia packages: Dates, Printf, TOML, JSON, PowerModels, LinearAlgebra, DelimitedFiles, JuMP, Gurobi. To install packages in Julia, the following code can be run:
```bash
# In Julia, first press ']' to enter pkg REPL mode.
pkg> add package_name   
```

To replicate the results in Tables 3 - 7, run:
```bash
julia  src/3level.jl network b k scenario method
```
with `network` being one of: `ieee_14`, `ieee_73`, `pegase_89`, `ieee_118`, `ieee_162`, `ieee_300`, `pegase_1354`, `rte_1888`, number of batteries and number of interdictions pair `(b,k)` can be `(2,3)`, `(2,5)`, `(3,5)` or `(5,10)`, `scenarios` ranging from 1 to 10, method being one of `reg` (regularized MIP formulation -- corresponding to formulation (8) in the paper), `lpr` (lienar relaxation -- corresponding to formulation (9) in the paper), and `vch` (vertex-representation convex hull -- corresponding to formulation (10) in the paper).

For example, to solve trilevel interdiction problem for network pegase_1354 with 5 batteries and 10 interdictions for demand scenario 1 using regularized MIP model, we can run:
 ```bash
julia  src/3level.jl pegase_1354 5 10 1 reg
```

To produce Figure 2, run:
```bash
julia  src/figure.jl
```
