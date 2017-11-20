# Copyright Wen-Chih Chen
# the code is for computing basic DEA models based on orientation (input or output)
# and returns-to-scale (crs or vrs)

using JuMP
using Gurobi # Gurobi is used as the LP solver
using FrontierEfficiencyAnalysis # for large scale problems and/or with size-limited solver

function basicDEAmodel()
  #----------------------------------------
  # data section
  data = readcsv("example.csv") # input User's (.csv) data path
  # each row is for a DMU; each column is an input or output
  # in the example the first two columns are the inputs while the others are the outputs
  inputs = data[:, 1:2]
  outputs = data[:, 3:end]

  # ---------------------------------------
  # to determine the efficiencies of all DMUs
  results = zeros(size(inputs)[1])
  getAllEfficiencies!(results, inputs, outputs) # get all DMU's efficiencies of input-oriented crs model
  ## or use the following to set the options of orientation, returns-to-scale and method
  # getAllEfficiencies!(results, inputs, outputs, orientation = "input", rts = "crs", method = 0)
  ## orientation = "input" (default) or "output"
  ## rts = "crs" (default) or "vrs"
  ## method = 0 (default) or 1 (by FrontierEfficiencyAnalysis for large scale problems and/or with size-limited solver)
  ## e.g. for output-oriented vsr model and solve by FrontierEfficiencyAnalysis do:
  # getAllEfficiencies!(results, inputs, outputs, orientation = "output", rts = "vrs", method = 1)
  println(results)

  # ---------------------------------------
  # to determine the efficiency of DMU k
  k = 2
  score = getEfficiency(inputs, outputs, k) # get DMU k's efficiency of input-oriented crs model
  ## or use the following to set the options of orientation, returns-to-scale and method
  # score = getEfficiency(inputs, outputs, k, orientation = "input", rts = "crs", method = 0)
  ## orientation = "input" (default) or "output"
  ## rts = "crs" (default) or "vrs"
  ## method = 0 (default) or 1 (by FrontierEfficiencyAnalysis for large scale problems and/or with size-limited solver)
  ## e.g. for output-oriented vsr model and solve by FrontierEfficiencyAnalysis do:
  # score = getEfficiency(inputs, outputs, k, orientation = "output", rts = "vrs", method = 1)
  println(score)

end

# compute DEA efficiencies for all DMUs in the data set and store results in scores
function getAllEfficiencies!(scores, inputData, outputData; orientation = "input", rts = "crs", method = 0)
  # parameters
  # scores: efficiency results
  # inputData: input data matrix (DMU size) x (input size)
  # outputData: output data matrix (DMU size) x (output size)
  # orientation: input- or output- oriented model (value= "input" or "output")
  # rts: returns to scale setting (value = "crs" or "vrs")
  # # method: how to solve LP (value = 0 (default) or 1 (by FrontierEfficiencyAnalysis.jl))

  if !(orientation == "input" || orientation == "output")
    error("orientation not supported")
  end

  if !(rts == "crs" || rts == "vrs")
    error("returns to scale not supported")
  end

  for t = 1:size(inputData)[1]
    scores[t] = getEfficiency(inputData, outputData, t, orientation = orientation, rts = rts, method = method)
  end
end

# compute DEA efficiency for DMU k
function getEfficiency(inputs, outputs, k; orientation = "input", rts = "crs", method = 0)
  # parameters
  # inputs: input data matrix (DMU size) x (input size)
  # outputs: output data matrix (DMU size) x (output size)
  # k: index of DMU to be evaluated
  # orientation: input- or output- oriented model (value= "input" or "output")
  # rts: returns to scale setting (value = "crs" or "vrs")
  # method: how to solve LP (value = 0 (default) or 1 (by FrontierEfficiencyAnalysis.jl))
  scale, inputNo = size(inputs)
  outputNo = size(outputs)[2]

  model = Model(solver = GurobiSolver(OutputFlag=0)) # Gurobi is used as the LP solver here. Users can choose their favorite solver.
  @variable(model, Lambda[1:scale] >= 0)
  @variable(model, Theta)

  if orientation == "input" # input-oriented model
    @objective(model, Min, Theta)
    @constraint(model, inputCon[i=1:inputNo], sum(Lambda[r]*inputs[r,i] for r = 1:scale) <= Theta*inputs[k,i])
    @constraint(model, outputCon[j=1:outputNo], sum(Lambda[r]*outputs[r,j] for r = 1:scale) >= outputs[k,j])
  elseif orientation == "output" # output-oriented model
    @objective(model, Max, Theta)
    @constraint(model, inputCon[i=1:inputNo], sum(Lambda[r]*inputs[r,i] for r = 1:scale) <= inputs[k,i])
    @constraint(model, outputCon[j=1:outputNo], sum(Lambda[r]*outputs[r,j] for r = 1:scale) >= Theta*outputs[k,j])
  else
    error("orientation not support")
  end

  if rts == "vrs"
    @constraint(model, sum(Lambda[r] for r = 1:scale) == 1)
  elseif rts == "crs"
  else
    error("returns to scale not supported")
  end

  if method == 0 # defalut
    solve(model)
  else  # using FrontierEfficiencyAnalysis for large scale problem or with size limited solver
    solveDEA(model)
  end

  return getvalue(Theta)
end

basicDEAmodel()
