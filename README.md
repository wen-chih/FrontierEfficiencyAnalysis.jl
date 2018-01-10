# FrontierEfficiencyAnalysis.jl

#### Copyright © 2017 by Wen-Chih Chen.  Released under the MIT License.

[![Build Status](https://travis-ci.org/wen-chih/FrontierEfficiencyAnalysis.jl.svg?branch=master)](https://travis-ci.org/wen-chih/FrontierEfficiencyAnalysis.jl)

FrontierEfficiencyAnalysis.jl is a package for Frontier Efficiency Analysis (aka Data Envelopment Analysis, DEA) computation. It is embedded in the [Julia](https://julialang.org/) programming language, and is an extension to the [JuMP](https://github.com/JuliaOpt/JuMP.jl) modeling language. It is particularly designed to enhance large-scale DEA computation and to solve DEA problems by size-limited solvers.

**Disclaimer** : FrontierEfficiencyAnalysis is *not* developed or maintained by the JuMP developers.


## Installation
In Julia, call `Pkg.add("FrontierEfficiencyAnalysis")` to install FrontierEfficiencyAnalysis.


## Usage
DEA is a linear program (LP)-based method used to determine a firm’s relative efficiency. Users can use JuMP to model and solve the DEA problems (special LP problems). Rather than solve the LPs by calling `JuMP.solve()`, FrontierEfficiencyAnalysis.jl can solve the large-scale problems more efficiently and/or by a solver with size limitation (e.g. 300 variables).


Please refer to [Quick Start Guide of JuMP](https://jump.readthedocs.io/en/latest/quickstart.html) for modeling details. What needed is to call our FrontierEfficiencyAnalysis.jl function:

	solveDEA(model)

instead of calling

	JuMP.solve(model)


## Example


```julia
# The example determine the efficiencies for all DMUs based on the CRS input-oriented model (CCR model)
using JuMP
using Gurobi # Gurobi is used as the LP solver
using FrontierEfficiencyAnalysis

data = readcsv("example.csv") # input User's (.csv) data path
scale, dimension = size(data) # scale is the number of DMU, dimension is the total number of inputs and outputs

for t = 1 : scale
    ### Modeling section
    # Here is the CRS input-oriented model (CCR model) to evaluate DMU t
    model = Model(solver = GurobiSolver()) # Gurobi is used as the LP solver here. Users can choose their favorite solver.
    @variable(model, Lambda[1:scale] >= 0)
    @variable(model, Theta)
    @objective(model, Min, Theta)
    @constraint(model, inputCon[i=1:2], sum(Lambda[r]*data[r,i] for r = 1:scale) <= Theta*data[t,i])
    @constraint(model, outputCon[j=3:dimension], sum(Lambda[r]*data[r,j] for r = 1:scale) >= data[t,j])
    # uncomment to add the convexity constraint for the VRS model
    # @constraint(model, sum(Lambda[r] for r = 1:scale) == 1)

    ### Problem solving
    solveDEA(model)

    ### Display
    println("Results for DMU $t")
    println("The efficicncy: $(getobjectivevalue(model))")
    println("The efficicncy: $(getvalue(Theta))")
    println("The lambdas: $(getvalue(Lambda)))")


    ## Use the following if returning dual values and slacks is needed

    ### Problem solving
    # duals, slacks = solveDEA(model)
    ### Display
    # println("Results for DMU $t")
    # println("The efficicncy: $(getobjectivevalue(model))")
    # println("The efficicncy: $(getvalue(Theta))")
    # println("The lambdas: $(getvalue(Lambda)))")

    # #println("The dual values (weights): $duals") # a vector associated with the constraints you define from the top to the bottom
    # println("the slack values: $slacks") # a vector associated with the constraints you define from the top to the bottom
end
```

<br>

## Parameters

>
**incrementSize** : the incremental size to expand the sample ( default value: 100 ).

	solveDEA(model, incrementSize = 200) # set the incremental size to 200

>
**tol** : the solution tolerance for solving DEA problem (default value: 1e-6). It also resets the dual feasibility tolerance in the solver to the given value.
<br>

	solveDEA(model, tol = 10^-4) # set the solution tolerance to 1e-4

>
**lpUB** : the size limit of the LP, i.e. the limitation of number of variables in the LP (default value: Inf).
<br>

	solveDEA(model, lpUB = 300) # set the LP size limitation to 300 variables

>
**extremeValueSetFlag** : to enable (=1) or disable (=0) performing initial sampling by selecing extreme value in each input/output dimension (default value: 0).
<br>


	solveDEA(model, extremeValueSetFlag = 1) # enable




## Citation
If you find FrontierEfficiencyAnalysis useful in your work, we kindly request that you cite the following papers

	@article{ChenLai2017,
	author = {Wen-Chih Chen and Sheng-Yung Lai},
	title = {Determining radial efficiency with a large data set by solving small-size linear programs},
	journal = {Annals of Operations Research},
	volume = {250},
	number = {1},
	pages = {147-166},
	year = {2017},
	doi = {10.1007/s10479-015-1968-4},
	}
and

	@misc{chen2017b,
	Author = {Wen-Chih Chen and Yueh-Shan Chung},
	Title = {A generalized non-radial efficiency measure and its application in DEA computation},
	Year = {2017},
	Eprint = {http://dx.doi.org/10.2139/ssrn.2496847},
	}

## Acknowledgements
FrontierEfficiencyAnalysis has been developed under the financial support of the Ministry of Science and Technology, Taiwan (Grant No. 104-2410-H-009-026-MY2). The contributors include Yueh-Shan Chung and Hao-Yun Chen.
