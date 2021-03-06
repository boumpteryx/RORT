using JuMP
using CPLEX
include("parser.jl")

export Static

function Static(fileName :: String,silent=true)
	cout_ouverture, Fct_commod, func_cost, func_capacity, nb_nodes, nb_arcs, nb_commodities, latency, node_capacity, commodity, nb_func, exclusion = read_instance(fileName)
	
	m = Model(CPLEX.Optimizer)
	if silent
		set_silent(m)
	end
	# Var
	@variable(m, x_i[1:nb_nodes], Bin) #activation d'une fonction au sommet i

	@variable(m, x_fi[1:nb_func,1:nb_nodes], Int) #activation de la fonction f au sommet i

	@variable(m, x_ikf[1:nb_nodes,1:nb_commodities,1:nb_func], Bin) #activation de la fonction f au sommet i pour le traitement de la commodité k

	@variable(m, e[1:nb_nodes,1:nb_nodes,1:nb_commodities,1:nb_func+1], Bin) #passage par l'arc (i,j) par la commodité k lpour traitement par la fonction f, ou en direction du puit si on considère e[i,j,k,end]

	#Constraint on states
	for i in 1:nb_nodes
		@constraint(m,[f in 1:nb_func],x_fi[f,i]*func_capacity[f] >= sum(x_ikf[i,k,f]*commodity[k,3] for k in 1:nb_commodities), base_name = "n_"*string(i)) #Nombre de fonctions f à placer en i
		@constraint(m,sum(x_fi[f,i] for f in 1:nb_func) <= node_capacity[i], base_name = "cap_"*string(i) ) #capacité de noeud
		for k in 1:nb_commodities
			@constraint(m,[f in 1:nb_func],x_i[i] >= x_ikf[i,k,f], base_name = "ouverture_"*string(i)) #ouverture du noeud en i
			if size(exclusion[k,:])[1]>0
				@constraint(m,sum(x_ikf[i,k,w] for w in exclusion[k,:] if w!=0) <= 1, base_name = "exclusion_"*string(k)*"_"*string(i)) #exclusion
			end
			for j in 1:nb_nodes
				if latency[i,j] == 0
					@constraint(m,[f in 1:nb_func+1],e[i,j,k,f]==0, base_name = "no_arc_"*string(i)*"_"*string(j)) #absence d'arcs
				end
			end
		end
	end


	#Contrainte de flots successifs pour chaque commodité
	for k in 1:nb_commodities
		size_fk=length( findall( y -> y > 0, Fct_commod[k,:]))

		tab_fk=sortperm( Fct_commod[k,:])[1:size_fk]

		@constraint(m,sum(e[commodity[k,1],j,k,tab_fk[1]] - e[j,commodity[k,1],k,tab_fk[1]] for j in 1:nb_nodes)-x_ikf[commodity[k,1],k,tab_fk[1]]==1, base_name = "init_flot_"*string(k)) #initialisation du flot à la source
		@constraint(m,[i in 1:nb_nodes ; i!=commodity[k,1]],sum(e[i,j,k,tab_fk[1]] - e[j,i,k,tab_fk[1]]  for j in 1:nb_nodes)+x_ikf[i,k,tab_fk[1]]==0, base_name = "1er_flot_"*string(k)) #flot pour traiter la première fonction de puis la source
		for t in 2:size_fk
			@constraint(m,[i in 1:nb_nodes],sum(e[i,j,k,tab_fk[t]] - e[j,i,k,tab_fk[t]] for j in 1:nb_nodes)-x_ikf[i,k,tab_fk[t-1]]+x_ikf[i,k,tab_fk[t]]==0, base_name = "mid_flow_"*string(t)*"_"*string(k)) #flot d'un traitement au suivant
		end
		@constraint(m,[i in 1:nb_nodes ; i!=commodity[k,2]],sum(e[i,j,k,nb_func+1] - e[j,i,k,nb_func+1] for j in 1:nb_nodes)-x_ikf[i,k,tab_fk[end]]==0, base_name = "last_flow_"*string(k)) #flot du dernier traitement vers le puit
		@constraint(m,sum(e[j,commodity[k,2],k,end] for j in 1:nb_nodes)+x_ikf[commodity[k,2],k,tab_fk[end]]==1, base_name = "end_flow_"*string(k)) #fin de flot sur le puit
	end

	for k in 1:nb_commodities
		for f in 1:nb_func
			@constraint(m,sum(x_ikf[i,k,f] for i in 1:nb_nodes)<=1, base_name = "concatenation_"*string(k)*"_"*string(f)) #concatenation des problemes de flot
		end
	end

	#Contraintes de latence
	@constraint(m, [k in 1:nb_commodities], sum( sum( sum(e[i,j,k,f]*latency[i,j] for f in 1:nb_func+1) for j in 1:nb_nodes) for i in 1:nb_nodes) <= commodity[k,4], base_name = "latence")

	#Objective
	@objective(m, Min, sum(x_i[i]*cout_ouverture + sum(x_fi[f,i]*func_cost[f,i] for f in 1:nb_func) for i in 1:nb_nodes) )

	optimize!(m)

	if !silent
		println(solution_summary(m))

		e=JuMP.value.(e)
		x_ikf=JuMP.value.(x_ikf)

		write_results(fileName,e,x_ikf)
	end
	status = termination_status(m)
	isOptimal = status == MOI.OPTIMAL
	return isOptimal, JuMP.value.(x_i), JuMP.value.(x_fi), JuMP.value.(x_ikf), JuMP.value.(e)
end

Static("di-yuan/di-yuan_1/")
