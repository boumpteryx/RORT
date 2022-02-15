using JuMP
using CPLEX
include("parser.jl")


export Relaxed

function Relaxed(fileName :: String)
	cout_ouverture, Fct_commod, func_cost, func_capacity, nb_nodes, nb_arcs, nb_commodities, latency, node_capacity, commodity, nb_func, exclusion = read_instance(fileName)

	tab=Vector{int64}(zeros(n_commodities))
	for k in 1:nb_commodities
		tab[k]=size(Fct_commod[k,:])
	end

	m = Model(CPLEX.Optimizer)

	# Var
	@variable(m, x_i[1:nb_nodes] >= 0) #activation d'une fonction au sommet i RELAXE

	@variable(m, x_fi[1:nb_func,1:nb_nodes], Int64) #activation de la fonction f au sommet i

	@variable(m, x_ikf[1:nb_nodes,1:nb_commodities,1:nb_func] >= 0) #activation de la fonction f au sommet i pour le traitement de la commodité k RELAXE

	@variable(m, e[1:nb_nodes,1:nb_nodes,1:nb_commodities,1:nb_func+1] >= 0) #passage par l'arc (i,j) par la commodité k lpour traitement par la fonction f, ou en direction du puit si on considère e[i,j,k,end] RELAXE

	# Constraint on states
  @constraint(m, [i in 1:nb_nodes], x_i[i] <= 1) # RELAXE
  @constraint(m, [i in 1:nb_nodes, k in 1:nb_commodities, f in 1:nb_func], x_ikf[i,k,f] >= 0) # RELAXE
  @constraint(m, [i in 1:nb_nodes, j in 1:nb_nodes, k in 1:nb_commodities, f in 1:nb_func+1], e[i,j,k,f] <= 1) # RELAXE


	for i in 1:nb_nodes
		@constraint(m,[f in 1:nb_func],x_fi[f,i]*func_capacity[f] >= sum(x_ikf[i,k,f]*commodity[k,3] for k in 1:nb_commodities)) #Nombre de fonctions f à placer en i
		@constraint(m,sum(x_fi[f,i] for f in 1:nb_func) <= node_capacity[i] ) #capacité de noeud
		for k in 1:nb_commodities
			@constraint(m,[f in 1:nb_func],x_i[i] >= x_ikf[i,k,f]) #ouverture du noeud en i
			@constraint(m,[f in 1:nb_func],x_ikf[i,k,f] <= length( findall( y -> y == f, Fct_commod[i] ))) #fixer à 0 x_ikf si fonction f n'est pas à appliquer à k
			@constraint(m,[(v,w) in exclusion[k]], x_ikf[i,k,v]+x_ikf[i,k,w] <= 1) #exclusion
			for j in 1:nb_nodes
				if latency[i,j] == 0
					@constraint(m,[f in 1:nb_func],e[i,j,k,f]==0) #absence d'arcs
				end
			end
		end
	end

	for k in 1:nb_commodities
		for f in 1:nb_func
			@constraint(m,sum(x_ikt[i,k,f] for i in 1:nb_nodes)==1) #concatenation des problemes de flot
		end
	end

	#Contrainte de flots successifs pour chaque commodité
	for k in 1:nb_commodities
		@constraint(m,sum(e[commodity[k,1],j,k,Fct_commod[k,1]] - e[j,commodity[k,1],k,Fct_commod[k,1]] for j in 1:nb_nodes)-x_ikf[commodity[k,1],k,Fct_commod[k,1]]==1) #initialisation du flot à la source
		@constraint(m,[i in 1:nb_nodes ; i!=commodity[k,1]],sum(e[i,j,k,Fct_commod[k,1]] - e[j,i,k,Fct_commod[k,1]]  for j in 1:nb_nodes)+x_ikf[i,k,Fct_commod[k,1]]==0) #flot pour traiter la première fonction de puis la source
		for t in 2:tab[k]
			@constraint(m,[i in 1:nb_nodes],sum(e[i,j,k,Fct_commod[k,t]] - e[j,i,k,Fct_commod[k,t]] for j in 1:nb_nodes)-x_ikf[i,k,Fct_commod[k,t-1]]+x_ikf[i,k,Fct_commod[k,t]]==0) #flot d'un traitement au suivant
		end
		@constraint(m,[i in 1:nb_nodes ; i!=commodity[k,2]],sum(e[i,j,k,Fct_commod[k,end]] - e[j,i,k,Fct_commod[k,end]] for j in 1:nb_nodes)-x_ikf[i,k,Fct_commod[k,end]]==0) #flot du dernier traitement vers le puit
		@constraint(m,sum(e[j,commodity[k,2],k,end] for j in 1:nb_nodes)+x_ikf[commodity[k,2],k,Fct_commod[k,end]]==1) #fin de flot sur le puit
	end

	#Contraintes de latence
	@constraint(m, [k in 1:nb_commodities], sum( sum( sum(e[i,j,k,f]*latency[i,j] for f in 1:nb_func+1) for j in 1:nb_nodes) for i in 1:nb_nodes) <= commodity[k,4])

	#Objective
	@objective(m, Min, sum(x_i[i]*cout_ouverture + sum(x_fi[f,i]*func_cost[f,i] for f in 1:nb_func) for i in 1:n) )

	optimize!(m)
	println(solution_summary(m))

	vX = JuMP.value.(x)

	status = termination_status(m)
	isOptimal = status == MOI.OPTIMAL

	return isOptimal
end
