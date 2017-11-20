# Copyright: Wen-Chih Chen
# The example detects outliers based on Chen and Johnson (2010)
# Chen, W. C., & Johnson, A. L. (2010). A unified model for detecting efficient and inefficient outliers
# in data envelopment analysis. Computers & Operations Research, 37(2), 417-425.
# it will write analysis result to output.csv

using JuMP
using Gurobi # Gurobi is used as the LP solver
# uncomment the following for a large-scale problem or with size-limited solvers
# using FrontierEfficiencyAnalysis

function outlierDetection()
    # ---------------------------------
    # data section
    # ---------------------------------
    # read data from a csv file
    println(">> read data")
    data = readcsv("outlierExample.csv") # input User's (.csv) data path
    # each row is for a DMU; each column is an input or output

    # you need to set the inputNo correspondingly
    inputNo = 1 # the first two columns are the inputs while the others are the outputs
    # determine the size of the data
    scale, dimension = size(data) # scale is the number of DMU, dimension is the total number of inputs and outputs
    inputs = data[:, 1:inputNo] # input data
    outputs = data[:, inputNo+1:dimension] # output data
    # end of data section
    # we have input data in inputs; output data in outputs
    # ---------------------------------

    println(">> solving ...")
    # output-oriented
    # outerImpact: vector for the total delta_o for removing the corresponding DMU
    # innerImpact: vector for the total delta_i for removing the corresponding DMU
    # outerImpactNo: vector for the no. of DMU impacted due to removing the corresponding DMU
    # innerImpactNo: vector for the no. of DMU impacted due to removing the corresponding DMU
    outerImpact, innerImpact, outerImpactNo, innerImpactNo = getDeltaOnY(inputs, outputs)
    totalImpact = outerImpact + innerImpact # eq. (8)

    # # for input-oriented models
    # outerImpact, innerImpact, outerImpactNo, innerImpactNo = getDeltaOnX(inputs, outputs)
    # totalImpact = outerImpact + innerImpact # eq. (15)

    # output section: export reslts to output.csv 
    writeResult("output.csv", [collect(1:scale)';outerImpact'; outerImpactNo'; innerImpact'; innerImpactNo']')

    println("done!")

end

# detect outlier measures in output-oriented models
# return (6) and (7) in Chen and Johnson (2010) in absolute value
function getDeltaOnY(inputs, outputs)
    scale = size(inputs)[1]
    outerDistWR = zeros(scale)
    innerDistWR = zeros(scale)
    for k = 1:scale
        outerDistWR[k] = outerDistOnY(inputs, outputs, inputs[k,:], outputs[k,:])
        innerDistWR[k] = innerDistOnY(inputs, outputs, inputs[k,:], outputs[k,:])
    end

    outerDistDelta = zeros(scale)
    innerDistDelta = zeros(scale)
    outerImpactNo = zeros(scale)
    innerImpactNo = zeros(scale)
    for delPt = 1:scale # DMU to be reomoved (R)
        set = collect(1:scale)
        deleteat!(set, delPt)

        # determine the total impact of removing delPt for each DMU except of delPt itself
        for k = 1:scale
            if k != delPt
                outerDiff = abs(outerDistWR[k]-outerDistOnY(inputs[set, :], outputs[set, :], inputs[k, :], outputs[k, :]))
                outerDistDelta[delPt] = outerDistDelta[delPt] + outerDiff
                if outerDiff > 0
                    outerImpactNo[delPt] = outerImpactNo[delPt] + 1
                end
                innerDiff = abs(innerDistWR[k]-innerDistOnY(inputs[set, :], outputs[set, :], inputs[k, :], outputs[k, :]))
                innerDistDelta[delPt] = innerDistDelta[delPt] + innerDiff
                if innerDiff > 0
                    innerImpactNo[delPt] = outerImpactNo[delPt] + 1
                end
            end
        end
    end
    return outerDistDelta, innerDistDelta, outerImpactNo, innerImpactNo
end

# detect outlier measures in input-oriented models
# return (13) and (14) in Chen and Johnson (2010) in absolute value
function getDeltaOnX(inputs, outputs)
    scale = size(inputs)[1]
    outerDistWR = zeros(scale)
    innerDistWR = zeros(scale)
    for k = 1:scale
        outerDistWR[k] = outerDistOnX(inputs, outputs, inputs[k,:], outputs[k,:])
        innerDistWR[k] = innerDistOnX(inputs, outputs, inputs[k,:], outputs[k,:])
    end

    outerDistDelta = zeros(scale)
    innerDistDelta = zeros(scale)
    outerImpactNo = zeros(scale)
    innerImpactNo = zeros(scale)
    for delPt = 1:scale # DMU to be reomoved (R)
        set = collect(1:scale)
        deleteat!(set, delPt)

        # determine the total impact of removing delPt for each DMU except of delPt itself
        for k = 1:scale
            if k != delPt
                outerDiff = abs(outerDistWR[k]-outerDistOnX(inputs[set, :], outputs[set, :], inputs[k, :], outputs[k, :]))
                outerDistDelta[delPt] = outerDistDelta[delPt] + outerDiff
                if outerDiff > 0
                    outerImpactNo[delPt] = outerImpactNo[delPt] + 1
                end
                innerDiff = abs(innerDistWR[k]-innerDistOnX(inputs[set, :], outputs[set, :], inputs[k, :], outputs[k, :]))
                innerDistDelta[delPt] = innerDistDelta[delPt] + innerDiff
                if innerDiff > 0
                    innerImpactNo[delPt] = outerImpactNo[delPt] + 1
                end
            end
        end
    end
    return outerDistDelta, innerDistDelta, outerImpactNo, innerImpactNo
end

# return the distance to the output boundary wrt to DMU t in an output-oriented analysis
# eq. (1) in Chen and Johnson (2010)
# inputs: input data; outputs: output data; t: DMU t under evaluation
function outerDistOnY(inputs, outputs, inputRHS, outputRHS)
    scale, inputNo = size(inputs)
    outputNo = size(outputs)[2]

    # modeling begins
    model = Model(solver = GurobiSolver(OutputFlag=0)) # Gurobi is used as the LP solver here. Users can choose their favorite solver.
    @variable(model, Lambda[1:scale] >= 0)
    @variable(model, Theta)
    @objective(model, Max, Theta)

    @constraint(model, inputCon[i=1:inputNo], sum(Lambda[r]*inputs[r,i] for r = 1:scale) == inputRHS[i])
    @constraint(model, outputCon[j=1:outputNo], sum(Lambda[r]*outputs[r,j] for r = 1:scale) == Theta*outputRHS[j])
    @constraint(model, sum(Lambda[r] for r = 1:scale) == 1)

    # use solveDEA(model) instead of solve(model) in the following for a large-scale problem or with size-limited solvers
    solve(model)
    # solveDEA(model)

    return getvalue(Theta)
end

# return the distance to the inner boundary wrt to DMU t in an output-oriented analysis
# eq. (2) in Chen and Johnson (2010)
# inputs: input data; outputs: output data; t: DMU t under evaluation
function innerDistOnY(inputs, outputs, inputRHS, outputRHS)
    scale, inputNo = size(inputs)
    outputNo = size(outputs)[2]
    # modeling begins
    model = Model(solver = GurobiSolver(OutputFlag=0)) # Gurobi is used as the LP solver here. Users can choose their favorite solver.
    @variable(model, Lambda[1:scale] >= 0)
    @variable(model, Theta)
    @objective(model, Min, Theta)
    @constraint(model, inputCon[i=1:inputNo], sum(Lambda[r]*inputs[r,i] for r = 1:scale) == inputRHS[i])
    @constraint(model, outputCon[j=1:outputNo], sum(Lambda[r]*outputs[r,j] for r = 1:scale) == Theta*outputRHS[j])
    @constraint(model, sum(Lambda[r] for r = 1:scale) == 1)

    # use solveDEA(model) instead of solve(model) in the following for a large-scale problem or with size-limited solvers
    solve(model)
    # solveDEA(model)

  return getvalue(Theta)
end

# return the distance to the output boundary wrt to DMU t in an output-oriented analysis
# eqs. (11) and (12) in Chen and Johnson (2010)
# inputs: input data; outputs: output data; t: DMU t under evaluation
function outerDistOnX(inputs, outputs, inputRHS, outputRHS)
    scale, inputNo = size(inputs)
    outputNo = size(outputs)[2]
    # modeling begins
    model = Model(solver = GurobiSolver(OutputFlag=0)) # Gurobi is used as the LP solver here. Users can choose their favorite solver.
    @variable(model, Lambda[1:scale] >= 0)
    @variable(model, Theta)
    @objective(model, Max, Theta)
    @constraint(model, inputCon[i=1:inputNo], sum(Lambda[r]*inputs[r,i] for r = 1:scale) == Theta*inputRHS[i])
    @constraint(model, outputCon[j=1:outputNo], sum(Lambda[r]*outputs[r,j] for r = 1:scale) == outputRHS[j])
    @constraint(model, sum(Lambda[r] for r = 1:scale) == 1)

    # use solveDEA(model) instead of solve(model) in the following for a large-scale problem or with size-limited solvers
    solve(model)
    # solveDEA(model)

    return getvalue(Theta)
end

# return the distance to the output boundary wrt to DMU t in an output-oriented analysis
# eqs. (9) and (10) in Chen and Johnson (2010)
# inputs: input data; outputs: output data; t: DMU t under evaluation
function innerDistOnX(inputs, outputs, inputRHS, outputRHS)
    scale, inputNo = size(inputs)
    outputNo = size(outputs)[2]
    # modeling begins
    model = Model(solver = GurobiSolver(OutputFlag=0)) # Gurobi is used as the LP solver here. Users can choose their favorite solver.
    @variable(model, Lambda[1:scale] >= 0)
    @variable(model, Theta)
    @objective(model, Min, Theta)
    @constraint(model, inputCon[i=1:inputNo], sum(Lambda[r]*inputs[r,i] for r = 1:scale) == Theta*inputRHS[i])
    @constraint(model, outputCon[j=1:outputNo], sum(Lambda[r]*outputs[r,j] for r = 1:scale) == outputRHS[j])
    @constraint(model, sum(Lambda[r] for r = 1:scale) == 1)

    # use solveDEA(model) instead of solve(model) in the following for a large-scale problem or with size-limited solvers
    solve(model)
    # solveDEA(model)

    return getvalue(Theta)
end

# write result to file
function writeResult(file, result)
    # write
    f = open(file, "w")
      write(f, "DMU, tol delta_O, no of DMUs impacted (O), tol delta_I, no of DMUs impacted (I)\r\n")
    close(f)

    for r = 1 : size(result)[1]
        f = open(file, "a")

        for c = 1 : size(result)[2]
            write(f, "$(result[r,c]),")
        end
        write(f, "\r\n")
        close(f)
    end
end

outlierDetection()
