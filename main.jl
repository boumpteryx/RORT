using JuMP
using CPLEX
include("parser.jl")


function myModel(MyFileName::String)
  
  # Create the model
  m = Model(CPLEX.Optimizer)

  ## Variables

  ## Constraints

  ## Objective

  #resolution
  optimize!(m)
end
