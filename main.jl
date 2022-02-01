using JuMP
using CPLEX
include("parser.jl")


function myModel(MyFileName::String)
  cout_ouverture, Fct_commod, func_cost, func_capacity, nb_nodes, nb_arcs, nb_commodities, latency, nb_functions_per_node, commodity, nb_func, exclusion = read_instance(MyFileName)

  # Create the model
  m = Model(CPLEX.Optimizer)

  ## Variables


  ## Constraints

  ## Objective


  #resolution
  optimize!(m)
end
