using JuMP
using Gurobi
using MathProgBase

# solver_Outputflag_check check if user didn't turn Original solver OutputFlag off
function solver_Outputflag_check(solver::ASCIIString)
  if contains(solver , "Any[<:OutputFlag,0>]") == false
    if contains(solver , "GurobiSolver")
      return "GurobiSolver"
    ###elseif other solver can be used in JuMP
    end
  end
end

# getFullsizeModelData: getting the data of full-size model which is entered by users
function getFullsizeModelData(m::JuMP.Model)
  # parameter
  # input:
  # m: the model entered by users using JuMP
  # output:
  # fullsizeObjSense: the objective sense of full-size model (Max or Min)
  # fullsizeObjCoeff: the coefficients of objective function of full-size model
  # fullsizeConstLB: the lower bound of each constraint of full-size model, which is a vector
  # fullsizeConstUB: the upper bound of each constraint of full-size model, which is a vector
  # fullsizeConstMatrix: the constraint matrix of full-size model
  # fullsizeVarLB: the lower bound of each varible of full-size model, which is a vector
  # fullsizeVarUB: the upper bound of each varible of full-size model, which is a vector
  # fullsizeNumVar: number of variables of full-size model
  # fullsizeNumConst: number of constraints of full-size model

  fullsizeObjSense = m.objSense
  fullsizeObjCoeff, fullsizeConstLB, fullsizeConstUB = JuMP.prepProblemBounds(m)
  fullsizeConstMatrix = JuMP.prepConstrMatrix(m)
  fullsizeVarLB = m.colLower
  fullsizeVarUB = m.colUpper

  fullsizeNumVar = m.numCols
  fullsizeNumConst = length(fullsizeConstLB)

  return fullsizeObjSense, fullsizeObjCoeff, fullsizeConstLB, fullsizeConstUB, fullsizeConstMatrix, fullsizeVarLB, fullsizeVarUB, fullsizeNumVar, fullsizeNumConst
end

# getBigger: getting the top n largest elements in vector v
function largestn(v::Vector{Float64}, n::Int64)
  # function getBigger(v::Vector{Float64}, n::Int)
  # parameters
  # input
  # v: the vector which is sorted to find top n largest elements
  # n: the number of elements which have to be finded from vector v
  # output
  # TopN_Idx: top n largest elements, which are the indices w.r.t the vector v
  # TopN_Value: values of top n largest elements

  # initialization
  TopN_Value = v[1:n]  # the top n largest value in vector v
  TopN_Idx = collect(1:n) # the top n largest elements in vector v

  minValue, minIdx = findmin(TopN_Value)

  for i = n+1:length(v)
    if v[i] > minValue
      TopN_Value[minIdx] = v[i]
      TopN_Idx[minIdx] = i

      minValue, minIdx = findmin(TopN_Value)
    end
  end

  return TopN_Value, TopN_Idx
end

# smallestn: getting the top n smallest elements in vector v
function smallestn(v::Vector{Float64}, n::Int64)
  # parameters
  # input
  # v: the vector which is sorted to find top n smallest elements
  # n: the number of elements which have to be finded from vector v
  # output
  # BottomN_Idx: top n smallest elements, which are the indices w.r.t the vector v
  # BottomN_Value: values of top n smallest elements

  # initialization
  BottomN_Value = v[1:n] # the top n smallest value in vector v
  BottomN_Idx = collect(1:n) # the top n smallest elements in vector v

  maxValue, maxIdx = findmax(BottomN_Value)

  for i = n+1:length(v)
    if v[i] < maxValue
      BottomN_Value[maxIdx] = v[i]
      BottomN_Idx[maxIdx] = i

      maxValue, maxIdx = findmax(BottomN_Value)
    end
  end

  return BottomN_Value, BottomN_Idx
end
