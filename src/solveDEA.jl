using JuMP
using Gurobi
using MathProgBase

include("functions.jl")

# solveDEA: solving dea problems by a smart algorithm
function solveDEA(m::JuMP.Model ; incrementSize = 100 , tol = 10.0^-6 , lpUB = Inf , extremeValueSetFlag = 0)
  # parameter
  # m: model variable of JuMP which is entered by JuMP's scripts
  # incrementSize: smallLP increment size when its size under lpUP for resampling
  # tol: the accuracy for solving DEA problem
  # lpUB: the limit size of LP
  # extremeValueSetFlag: first sampling take extremeValueSet or not (default =  turn off)

  # data preprocessing
  # getting the data of full-size model which is entered by users
  fullsizeObjSense, fullsizeObjCoeff, fullsizeConstLB, fullsizeConstUB, fullsizeConstMatrix, fullsizeVarLB, fullsizeVarUB, fullsizeNumVar, fullsizeNumConst = getFullsizeModelData(m)
  denseFullsizeConstMatrix = full(fullsizeConstMatrix)

  # initial varibles related with data model
  iterations = 0
  dimension, scale = size(denseFullsizeConstMatrix)
  lpSize = 0
  additiveSize = 0
  precision = round(Int,round(log(10, 1/tol) ,1) )

  # user parameters simple exception catch
  lpUB = (lpUB > scale) ? scale : lpUB
  incrementSize = (incrementSize > lpUB) ? lpUB : incrementSize
  increment = [incrementSize , incrementSize]

  # create dictionary for fulldatamatrix mapping
  fulldatamatrix_dict = Dict([denseFullsizeConstMatrix[1:end,i] => i for i in 1:length(denseFullsizeConstMatrix[1,:])])

  # create critical component and data component from fulldatamatrix
  dataConstMatrix_dict = find(fullsizeObjCoeff.==0)
  criticalConstMatrix = denseFullsizeConstMatrix[:,find(fullsizeObjCoeff.!=0)]
  dataConstMatrix = denseFullsizeConstMatrix[:,dataConstMatrix_dict]

  # create map for fulldatamatrix mapping to dataConstMatrix
  dataConstMatrix_inv_dict = zeros(Int64, scale)
  for i in 1:length(dataConstMatrix_dict)
    dataConstMatrix_inv_dict[dataConstMatrix_dict[i] ] = i
  end

  # first sampling
  smallLP = criticalConstMatrix[1:end,:]
  critical_size = length(Array(criticalConstMatrix[1,:]))
  lpSize = lpSize + critical_size

  # create map for whether dataConstMatrix used or not
  dataConstMatrix_map = zeros(length(dataConstMatrix[1,:]))
  # record dataConstMatrix map sequence
  dataConstMatrix_map_seq = zeros(lpUB - critical_size)

  if (extremeValueSetFlag == 0)
    # create small LP (Critical data + first datamatrix random data in incrementSize)
    randNum = rand(1:length(dataConstMatrix[1,:]) , incrementSize)
    smallLP = [smallLP  dataConstMatrix[:,randNum[:]] ]
    dataConstMatrix_map[randNum] = 1
    dataConstMatrix_map_seq[1 : incrementSize] = randNum
    lpSize = lpSize + incrementSize
  else
    # create small LP (Critical data + extremeValueSet)
    extremeValueSet = zeros(Int64, dimension)
    for i in 1 : dimension
      extremeValueSet[i] = (fullsizeConstUB[i] == Inf) ? findmin(dataConstMatrix[i,:])[2] : findmax(dataConstMatrix[i,:])[2]
    end
    extremeValueSet =  sort(extremeValueSet)
    smallLP = [smallLP  dataConstMatrix[:,extremeValueSet[:]] ]
    dataConstMatrix_map[extremeValueSet] = 1
    dataConstMatrix_map_seq[1 : dimension] = extremeValueSet
    lpSize = lpSize + dimension
  end

  # generate other small condition matrix
  smallVarLB = zeros(lpSize)
  smallVarUB = zeros(lpSize)
  smallObjCoeff = zeros(lpSize)

  for i in 1 : lpSize
    smallVarLB[i] = fullsizeVarLB[fulldatamatrix_dict[smallLP[:,i]]]
    smallVarUB[i] = fullsizeVarUB[fulldatamatrix_dict[smallLP[:,i]]]
    smallObjCoeff[i] = fullsizeObjCoeff[fulldatamatrix_dict[smallLP[:,i]]]
  end

  smallConstLB = fullsizeConstLB
  smallConstUB = fullsizeConstUB
  smallSense = fullsizeObjSense

  # definition and initialization solving Status
  optStatus = false   # the optimality status of model
  iterations = 0      # the iterations that algorithm uses to solve model

  # JuMP solver initial define
  ### 20170924 need to add other solver that can be used in JuMP
  if solver_Outputflag_check(string(m.solver)) == "GurobiSolver"# check Output flag option for user catch
    setsolver(m , GurobiSolver(OutputFlag = 0))
  end

  lm = MathProgBase.LinearQuadraticModel(m.solver)  # model variable of MathProgBase, we use this to build small model
  stat = :Default   # model status
  m.objVal = NaN      # objective value of full-size model
  m.colVal = fill(NaN, fullsizeNumVar)    # values of all variables of full-size model
  m.linconstrDuals = Array(Float64, 0)    # dual values of full-size model

  # loop until model being solved to optimal or algorithm processing over 1000 iterations
  while optStatus == false
    iterations = iterations + 1
    if iterations > 1000
      warn("Iterations > 1000")
      stat = :Unsolved
      break
    end
    optStatus = true

    # load problems
    lm = MathProgBase.LinearQuadraticModel(m.solver)
    MathProgBase.loadproblem!(lm, smallLP, smallVarLB, smallVarUB, smallObjCoeff, smallConstLB, smallConstUB, smallSense)
    MathProgBase.optimize!(lm)
    stat = MathProgBase.status(lm)

    # get constduals
    constDuals = MathProgBase.getconstrduals(lm)
    # create check optimal array
    checkopt_arr = constDuals.'* fullsizeConstMatrix

    # addictive size calculate
    if lpSize < lpUB # there's room to increase sample size
      additiveSize = increment[min(iterations, length(increment))]
      if lpSize + additiveSize > lpUB
        additiveSize = lpUB - lpSize
      end
    else
        # if the lpSize is already equal to lpUB, then let additiveSize to be 1
        additiveSize = 1
    end

    # check optimal function
    infeasibleIdx = Array(Int, additiveSize)
    infeasibleValue = Array(Float64, additiveSize)

    # create notset for data index which didn't use in dataConstMatrix
    notset = find(dataConstMatrix_map.==0)
    for i in 1:length(notset)
      notset[i] = dataConstMatrix_dict[notset[i] ]
    end

    if smallSense == :Min
      if additiveSize == 1
        infeasibleValue, infeasibleIdx = findmax(checkopt_arr[notset])
        infeasibleIdx = dataConstMatrix_inv_dict[notset[infeasibleIdx] ]
        dataConstMatrix_map[infeasibleIdx] = 1
        optStatus = infeasibleValue[1] <= tol
      else
        infeasibleValue, infeasibleIdx = largestn(checkopt_arr[notset], additiveSize)
        for i in 1:length(infeasibleIdx)
          infeasibleIdx[i] = dataConstMatrix_inv_dict[notset[infeasibleIdx[i]]]
        end
        dataConstMatrix_map[infeasibleIdx] = 1
        optStatus = findmax(infeasibleValue)[1] <= tol
      end
    else
      if additiveSize == 1
        infeasibleValue, infeasibleIdx = findmin(checkopt_arr[notset])
        infeasibleIdx = dataConstMatrix_inv_dict[notset[infeasibleIdx] ]
        dataConstMatrix_map[infeasibleIdx] = 1
        optStatus = infeasibleValue[1] >= -tol
      else
        infeasibleValue, infeasibleIdx = smallestn(checkopt_arr[notset], additiveSize)
        infeasibleIdx = dataConstMatrix_inv_dict[notset[infeasibleIdx] ]
        for i in 1:length(infeasibleIdx)
          infeasibleIdx[i] = dataConstMatrix_inv_dict[notset[infeasibleIdx[i]]]
        end
        optStatus = findmin(infeasibleValue)[1] >= -tol
      end
    end

    # resampling and dynamic resampling
    if optStatus == false
      if additiveSize == 1
        # create delete index
        critical_length = length(Array(criticalConstMatrix[1,:]))

        deleteIdx = rand(1 : lpUB - critical_length)
        # marked used and not optimal data for "2"
        dataConstMatrix_map[dataConstMatrix_map_seq[deleteIdx]] = 2
        # swap now element at delete index
        dataConstMatrix_map_seq[deleteIdx] = infeasibleIdx
        # swap dataConstMatrix
        smallLP[:,critical_length + deleteIdx] = dataConstMatrix[:,infeasibleIdx]

        # swap other small condition matrix
        smallVarLB[critical_length + deleteIdx] = fullsizeVarLB[fulldatamatrix_dict[smallLP[:,critical_length + deleteIdx]]]
        smallVarUB[critical_length + deleteIdx] = fullsizeVarUB[fulldatamatrix_dict[smallLP[:,critical_length + deleteIdx]]]
        smallObjCoeff[critical_length + deleteIdx] = fullsizeObjCoeff[fulldatamatrix_dict[smallLP[:,critical_length + deleteIdx]]]
      else
        # increment data with incrementSize into smallLP
        smallLP = [smallLP  dataConstMatrix[:,infeasibleIdx] ]

        dataConstMatrix_map_seq[lpSize - critical_size + 1 : lpSize - critical_size + additiveSize] = infeasibleIdx

        # increment data with incrementSize into other small condition matrix
        tmp_smallVarLB = zeros(additiveSize)
        tmp_smallVarUB = zeros(additiveSize)
        tmp_smallObjCoeff = zeros(additiveSize)
        for i in lpSize + 1 : lpSize + additiveSize
          tmp_smallVarLB[i - lpSize] = fullsizeVarLB[fulldatamatrix_dict[smallLP[:,i]]]
          tmp_smallVarUB[i - lpSize] = fullsizeVarUB[fulldatamatrix_dict[smallLP[:,i]]]
          tmp_smallObjCoeff[i - lpSize] = fullsizeObjCoeff[fulldatamatrix_dict[smallLP[:,i]]]
        end

        smallVarLB = [smallVarLB ; tmp_smallVarLB]
        smallVarUB = [smallVarUB ; tmp_smallVarUB]
        smallObjCoeff = [smallObjCoeff ; tmp_smallObjCoeff]

        lpSize = lpSize + additiveSize

      end
    end
  end

  # Extracting solution
  if stat == :Optimal
    # recording objective value
    m.objVal = round(MathProgBase.getobjval(lm), precision)

    # recording clambda (DEA solution)
    m.colVal = fill(0.0, fullsizeNumVar)
    solution = MathProgBase.getsolution(lm)
    for i in 1 : length(solution)
      m.colVal[fulldatamatrix_dict[smallLP[:,i] ] ] = round(solution[i], precision)
    end

    # recording dual value
    m.linconstrDuals = try
        round(MathProgBase.getconstrduals(lm)[1:fullsizeNumConst], precision)
    catch
        fill(NaN, fullsizeNumConst)
    end
    duals = m.linconstrDuals

    # recording slack value
    slack = zeros(Float64,dimension)
    for k in 1:dimension
      for i in 1:length(find(m.colVal.!=0))
        slack[k] += m.colVal[find(m.colVal.!=0)[i]] * denseFullsizeConstMatrix[k,find(m.colVal.!=0)[i]]
      end
      if smallConstLB[k] == -Inf
        slack[k] = round(abs(smallConstUB[k] - slack[k]), precision)
      else
        slack[k] = round(abs(slack[k] - smallConstLB[k]), precision)
      end
    end
  end
  return duals,slack
end
