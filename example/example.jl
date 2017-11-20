# Copyright Wen-Chih Chen
# The example determine the efficiencies for all DMUs based on the CRS input-oriented, CCR, model (Charnes et al., 1978)

using JuMP
using Gurobi # Gurobi is used as the LP solver
using FrontierEfficiencyAnalysis

# read data from a csv file
data = readcsv("example.csv") # input User's (.csv) data path
# each row is for a DMU; each column is an input or output
# the first two columns are the inputs while the others are the outputs
# determine the size of the data
scale, dimension = size(data) # scale is the number of DMU, dimension is the total number of inputs and outputs

# compute efficiencies for all DMUs
for t = 1 : scale
    ### Modeling section
    # Here is the CRS input-oriented model (CCR model) to evaluate DMU t
    model = Model(solver = GurobiSolver()) # Gurobi is used as the LP solver here. Users can choose their favorite solver.
    @variable(model, Lambda[1:scale] >= 0)
    @variable(model, Theta)
    @objective(model, Min, Theta)
    @constraint(model, inputCon[i=1:2], sum(Lambda[r]*data[r,i] for r = 1:scale) <= Theta*data[t,i])
    @constraint(model, outputCon[j=3:dimension], sum(Lambda[r]*data[r,j] for r = 1:scale) >= data[t,j])
    ## uncomment the following line to add the convexity constraint for the VRS model
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

    # println("The dual values (weights): $duals") # a vector associated with the constraints you define from the top to the bottom
    # println("the slack values: $slacks") # a vector associated with the constraints you define from the top to the bottom

end
