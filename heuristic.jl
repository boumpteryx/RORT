using JuMP
using CPLEX
include("parser.jl")

function PL_heuristic(cout_ouverture, Fct_commod, func_cost, func_capacity, nb_nodes, nb_arcs, nb_commodities, latency, nb_functions_per_node, commodity, nb_func, exclusion, y_fi, x_fi_prec)
  # version modifiee du PL normal
  # on ajoute y_fi qui decrit combien de place de fonction f au sommet i il reste de par les commodites precedentes. y_fi est une donnee du probleme TODO, optionnel mais important
  # on ajoute la contrainte de non-overlap entre les fonctions de x_fi_prec et celles posees dans cette instance DONE

  m = Model(CPLEX.Optimizer)

	# Var
	@variable(m, x_i[1:nb_nodes], Bin) # activation d'une fonction au sommet i
	@variable(m, x_fi[1:nb_func,1:nb_nodes], Int) # activation de la fonction f au sommet i
	@variable(m, x_ikf[1:nb_nodes,1:nb_commodities,1:nb_func], Bin) # activation de la fonction f au sommet i pour le traitement de la commodité k
	@variable(m, e[1:nb_nodes,1:nb_nodes,1:nb_commodities,1:nb_func+1], Bin) # passage par l'arc (i,j) par la commodité k lpour traitement par la fonction f, ou en direction du puit si on considère e[i,j,k,end]

	#Constraint on states
	for i in 1:nb_nodes
		@constraint(m,[f in 1:nb_func],x_fi[f,i]*func_capacity[f] >= sum(x_ikf[i,k,f]*commodity[k,3] for k in 1:nb_commodities), base_name = "n_"*string(f)*"_"*string(i)) # Nombre de fonctions f à placer en i
		@constraint(m,sum(x_fi[f,i] for f in 1:nb_func) <= node_capacity[i], base_name = "cap_"*string(i) ) #capacité de noeud
		for k in 1:nb_commodities
			@constraint(m,[f in 1:nb_func],x_i[i] >= x_ikf[i,k,f], base_name = "ouverture_"*string(i)*"_"*string(f)) #ouverture du noeud en i
			#@constraint(m,[f in 1:nb_func],x_ikf[i,k,f] <= length( findall( y -> y == f, Fct_commod[k] ))) #fixer à 0 x_ikf si fonction f n'est pas à appliquer à commodité k (donc pas sur le noeud i en particulier)
			if size(exclusion[k,:])[1]>0
				@constraint(m,sum(x_ikf[i,k,w] for w in exclusion[k,:]) <= 1, base_name = "exclusion_"*string(k)*"_"*string(i)) #exclusion
        @constraint(m,sum(x_ikf[i,k,w] for w in x_fi_prec) <= 1) #exclusion des precedents
			end
			for j in 1:nb_nodes
				if latency[i,j] == 0
					@constraint(m,[f in 1:nb_func],e[i,j,k,f]==0, base_name = "no_arc_"*string(i)*"_"*string(j)*"_"*string(f)) #absence d'arcs
				end
			end
		end
	end


	#Contrainte de flots successifs pour chaque commodité
	for k in 1:nb_commodities
		size_fk=length( findall( y -> y > 0, Fct_commod[k,:]))

		tab_fk=sortperm( Fct_commod[k,:])[1:size_fk]

		@constraint(m,sum(e[commodity[k,1],j,k,tab_fk[1]] - e[j,commodity[k,1],k,tab_fk[1]] for j in 1:nb_nodes)-x_ikf[commodity[k,1],k,tab_fk[1]]==1, base_name = "init_flot_"*string(k)) #initialisation du flot à la source
		@constraint(m,[i in 1:nb_nodes ; i!=commodity[k,1]],sum(e[i,j,k,tab_fk[1]] - e[j,i,k,tab_fk[1]]  for j in 1:nb_nodes)+x_ikf[i,k,tab_fk[1]]==0, base_name = "1er_flot_"*string(k)*"_"*string(i)) #flot pour traiter la première fonction de puis la source
		for t in 2:size_fk
			@constraint(m,[i in 1:nb_nodes],sum(e[i,j,k,tab_fk[t]] - e[j,i,k,tab_fk[t]] for j in 1:nb_nodes)-x_ikf[i,k,tab_fk[t-1]]+x_ikf[i,k,tab_fk[t]]==0, base_name = "mid_flow_"*string(t)*"_"*string(k)*"_"*string(i)) #flot d'un traitement au suivant
		end
		@constraint(m,[i in 1:nb_nodes ; i!=commodity[k,2]],sum(e[i,j,k,nb_func+1] - e[j,i,k,nb_func+1] for j in 1:nb_nodes)-x_ikf[i,k,tab_fk[end]]==0, base_name = "last_flow_"*string(k)*"_"*string(i)) #flot du dernier traitement vers le puit
		@constraint(m,sum(e[j,commodity[k,2],k,end] for j in 1:nb_nodes)+x_ikf[commodity[k,2],k,tab_fk[end]]==1, base_name = "end_flow_"*string(k)) #fin de flot sur le puit
	end

	for k in 1:nb_commodities
		for f in 1:nb_func
			@constraint(m,sum(x_ikf[i,k,f] for i in 1:nb_nodes)<=1, base_name = "concatenation_"*string(k)*"_"*string(f)) #concatenation des problemes de flot
		end
	end

	#Contraintes de latence
	@constraint(m, [k in 1:nb_commodities], sum( sum( sum(e[i,j,k,f]*latency[i,j] for f in 1:nb_func+1) for j in 1:nb_nodes) for i in 1:nb_nodes) <= commodity[k,4])

	#Objective
	@objective(m, Min, sum(x_i[i]*cout_ouverture + sum(x_fi[f,i]*func_cost[f,i] for f in 1:nb_func) for i in 1:nb_nodes) )

	optimize!(m)

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
  y_fi = Array{Int64,2}(zeros(nb_func,nb_nodes))
  x_fi_prec = Array{Int64,2}(zeros(nb_func,nb_nodes))

  x_i, current_objective, x_fi = PL_heuristic(cout_ouverture, Fct_commod[1], func_cost, func_capacity, nb_nodes, nb_arcs, 1, latency, nb_functions_per_node, commodity[1], nb_func, exclusion, y_fi, x_fi_prec)

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

    # update y_fi TODO


    # update x_fi_prec
    for k in 1:nb_nodes
      for j in 1:nb_func
        x_fi_prec[j,k] = x_fi_prec[j,k] + x_fi[j,k]
      end
    end

    # update sum_objectives
    sum_objectives = sum_objectives + current_objective
    x_i, current_objective, x_fi = PL_heuristic(cout_ouverture, Fct_commod[i], func_cost, func_capacity, nb_nodes, nb_arcs, 1, latency, nb_functions_per_node, commodity[i], nb_func, exclusion, y_fi, x_fi_prec)
  end

  return sum_objectives
end

heuristic("test")
