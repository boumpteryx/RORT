using JuMP
using CPLEX

export Static

function Static(fileName :: String)
	cout_ouverture, Fct_commod, func_cost, func_capacity, nb_nodes, nb_arcs, nb_commodities, latency, node_capacity, commodity, nb_func, exclusion = read_instance(fileName)

	tab=Vector{int64}(zeros(n_commodities))	
	for k in 1:nb_commodities
		tab[k]=size(Fct_commod[k,:])
	end
	
	m = Model(CPLEX.Optimizer)

	# Var
	@variable(m, x_i[1:nb_nodes], Bin)
	
	@variable(m, x_fi[1:nb_func,1:nb_nodes], Int64)
	
	@variable(m, x_ikt[1:nb_nodes,1:nb_commodities,1:nb_func], Bin)
	
	@variable(m, e[1:nb_nodes,1:nb_nodes,1:nb_commodities,1:nb_func], Bin)
	
	#Constraint on states
	for i in 1:nb_nodes
		@constraint(m,[f in 1:nb_func],x_fi[f,i]*func_capacity[f] >= sum(  sum(x_ikt[i,k,t]*commodity[k,3] for k in 1:nb_commodities) for t in 1:tab[i])) #Nombre de fonctions f à placer en i
		@constraint(m,sum(x_fi[f,i] for f in 1:nb_func) <= node_capacity[i] ) #capacité de noeud
		for k in 1:nb_commodities
			@constraint(m,[t in 1:tab[k]],x_i[i] >= x_ikt[i,k,t]) #presence d'une fonction en i
			@constraint(m,[(v,w) in exclusion[k]], x_ikt[i,k,v]+x_ikt[i,k,w] <= 1) #exclusion
			for j in 1:nb_nodes
				if latency[i,j]==0
					@constraint(m,[t in 1:tab[i]],e[i,j,k,t]==0) #absence d'arcs
				end
			end
		end
	end
	
	for k in 1:nb_commodities
		for t in 1:tab[
			@constraint(m,sum(x_ikt[i,k,t] for i in 1:nb_nodes)==1) #concatenation des problemes de flot
		end
	end
	
	@constraint(m,[i in 1:n ; i!=s && i!=t],sum(x[k,i]*A[k,i] - x[i,k]*A[i,k] for k in 1:n)==0)
	for j in 1:n
		@constraint(m,[i in 1:n],x[i,j] <= A[i,j])
	end
	
	#Constraint on source and well
	@constraint(m,sum(x[s,k] - x[k,s] for k in 1:n)==1)
	@constraint(m,sum(x[k,t] - x[t,k] for k in 1:n)==1)

	#Constraints on vertex weights
	@constraint(m, sum(p[i]*sum(x[i,j] for j in 1:n) for i in 1:n) + p[t] <= S)
	
	#Objective
	@objective(m, Min, sum(x_i[i]*cout_ouverture + sum(x_fi[f,i]*func_cost[f,i] for f in 1:nb_func) for i in 1:n) )

	optimize!(m)
	println(solution_summary(m))

	vX = JuMP.value.(x)

	status = termination_status(m)
	isOptimal = status == MOI.OPTIMAL
		
	return isOptimal
end

