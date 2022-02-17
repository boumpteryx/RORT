using JuMP
using CPLEX
include("parser.jl")

function PL_heuristic(cout_ouverture, Fct_commod, func_cost, func_capacity, nb_nodes, nb_arcs, nb_commodities, latency, nb_functions_per_node, commodity, nb_func, exclusion)
  # version normale pour iteration 1
  return JuMP.value.(x_i), JuMP.objective_value.(m),  JuMP.value.(x_fi)
end


function heuristic(MyFileName::String)
  cout_ouverture, Fct_commod, func_cost, func_capacity, nb_nodes, nb_arcs, nb_commodities, latency, nb_functions_per_node, commodity, nb_func, exclusion = read_instance(MyFileName)
  c = cout_ouverture
  cout_ouverture = Array{Int64,1}(zeros(nb_nodes))
  for i in 1:nb_nodes
    cout_ouverture[i] = c
  end
  sum_objectives = 0

  # TODO exclusion
  x_i, current_objective, x_fi = PL_heuristic(cout_ouverture, Fct_commod[1], func_cost, func_capacity, nb_nodes, nb_arcs, 1, latency, nb_functions_per_node, commodity[1], nb_func, exclusion)

  for i in 2:nb_commodities

    # update nb_functions_per_node
    for k in 1:nb_nodes
      for j in 1:nb_func
        nb_functions_per_node[k] = nb_functions_per_node[k] - x_fi[j,k]
      end
    end

    # update cout_ouverture
    for k in 1:nb_nodes
      cout_ouverture[k] = cout_ouverture[k] - c * x_i[k]
    end

    # update sum_objectives
    sum_objectives = sum_objectives + current_objective
    x_i, current_objective, x_fi = PL_heuristic(cout_ouverture, Fct_commod[i], func_cost, func_capacity, nb_nodes, nb_arcs, 1, latency, nb_functions_per_node, commodity[i], nb_func, exclusion)
  end

  return sum_objectives
end

heuristic("test")
