# A Modeling Language Package for Solving Data Envelopment Analysis Problems
module FrontierEfficiencyAnalysis

# An algebraic modeling language for Julia
import JuMP
# Solver-independent functions and low-level interface for Mathematical Programming
import MathProgBase
# An example solver
import Gurobi

# function
export solveDEA

# solving dea problems by a smart algorithm
# it is computationally efficient and can solve large scale DEA problems via any solvers (even with LP size limitation)
include("solveDEA.jl")

end
