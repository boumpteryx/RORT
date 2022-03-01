using JuMP
using CPLEX
include("parser.jl")
include("relaxed.jl")
include("base.jl")

export DW

function DW(fileName :: String, X_i::Array{Float64,2},  X_fi::Array{Float64,3}, X_ikf::Array{Float64,4}, E::Array{Float64,5})
		
	#Data acquisition, model definition
	cout_ouverture, Fct_commod, func_cost, func_capacity, nb_nodes, nb_arcs, nb_commodities, latency, node_capacity, commodity, nb_func, exclusion = read_instance(fileName)

		
	#Column generation
	not_opti=true
	while not_opti
		m = Model(CPLEX.Optimizer)
		set_silent(m)
		
		#variable definition
		n=size(X_ikf,1)
		@variable(m, lbd[1:n] >= 0)
		
		#Convexity constraint
		@constraint(m,sum(lbd[w] for w in 1:n)==1, base_name="convex")
		
		#Constraint on states
		for i in 1:nb_nodes
			@constraint(m,[f in 1:nb_func],sum(lbd[w]*X_fi[w,f,i] for w in 1:n)*func_capacity[f] >=  sum( sum(lbd[w]*X_ikf[w,i,k,f] for w in 1:n) * commodity[k,3] for k in 1:nb_commodities), base_name = "n_"*string(i)) # (3) Nombre de fonctions f à placer en i
			for k in 1:nb_commodities
				@constraint(m,[f in 1:nb_func],sum(lbd[w]*X_i[w,i] for w in 1:n) >= sum(lbd[w]*X_ikf[w,i,k,f] for w in 1:n), base_name = "ouverture_"*string(i)*"_"*string(k)) # (2) ouverture du noeud en i
			end
		end
		
		#Objective
		@objective(m, Min, sum(sum(lbd[w]*X_i[w,i] for w in 1:n)*cout_ouverture + sum(sum(lbd[w]*X_fi[w,f,i] for w in 1:n)*func_cost[f,i] for f in 1:nb_func) for i in 1:nb_nodes) )
		optimize!(m)

		#Calcul dual values and solutions
		eta=abs(dual(constraint_by_name(m, "convex")))
		mu_2=Array{Float64,3}(zeros(nb_nodes,nb_commodities,nb_func))
		mu_3=Array{Float64,2}(zeros(nb_nodes,nb_func))
		for i in 1:nb_nodes
			for f in 1:nb_func
				mu_3[i,f]=dual(constraint_by_name(m,"n_"*string(i)*"["*string(f)*"]"))
				for k in 1:nb_commodities
					mu_2[i,k,f]=dual(constraint_by_name(m,"ouverture_"*string(i)*"_"*string(k)*"["*string(f)*"]"))
				end
			end
		end		
		lbd_sol=JuMP.value.(lbd)
		x_i_sol=sum(lbd_sol[w].*X_i[w,:] for w in 1:n)
		x_fi_sol=sum(lbd_sol[w].*X_fi[w,:,:] for w in 1:n)
		min_cr=-eta
		
		#setup for next column to add on first index
		X_i=vcat(X_i[1:1,:],X_i)
		X_fi=vcat(X_fi[1:1,:,:],X_fi)
		X_ikf=vcat(X_ikf[1:1,:,:,:],X_ikf)
		E=vcat(E[1:1,:,:,:,:],E)
		
		#state sub-problem resolution
		m0 = Model(CPLEX.Optimizer)
		set_silent(m0)
		@variable(m0, x_i[1:nb_nodes] >= 0)
		@variable(m0, x_fi[1:nb_func,1:nb_nodes]>=0)
		@constraint(m0,[i in 1:nb_nodes], x_i[i] <= 1) # Relaxed binary constraint
		#@constraint(m0,[i in 1:nb_nodes], sum(x_fi[f,i] for f in 1:nb_func) <= node_capacity[i], base_name = "cap" ) # (4) capacité de noeud
		@objective(m0, Min, sum(x_i[i]*(cout_ouverture-sum(sum(mu_2[i,k,f] for f in 1:nb_func) for k in nb_commodities))+sum(x_fi[f,i]*(func_cost[f,i]-mu_3[i,f]*func_capacity[f]) for f in 1:nb_func)  for i in 1:nb_nodes) )
		optimize!(m0)

		
		#Update on data
		min_cr=min_cr+objective_value(m0)
		X_i[1,:]=JuMP.value.(x_i)
		X_fi[1,:,:]=JuMP.value.(x_fi)

		#path sub-problems resolution
		for k in 1:nb_commodities
						
			#model definition
			mk = Model(CPLEX.Optimizer)
			set_silent(mk)
			@variable(mk, x_ikf[1:nb_nodes,1:nb_func] >= 0) 
			@variable(mk, e[1:nb_nodes,1:nb_nodes,1:nb_func+1] >= 0)
			
			#Constraints of exclusion and possible arc
			for i in 1:nb_nodes
				if size(exclusion[k,:])[1]>0
					@constraint(mk,sum(x_ikf[i,w] for w in exclusion[k,:]) <= 1, base_name = "exclusion_"*string(k)*"_"*string(i)) #exclusion
				end				
				for j in 1:nb_nodes
					if latency[i,j] == 0
						@constraint(mk,[f in 1:nb_func+1],e[i,j,f]==0, base_name = "no_arc_"*string(i)*"_"*string(j)) #absence d'arcs
					end
				end
				
			end
			
			#Contrainte de flots successifs pour chaque commodité
			size_fk=length( findall( y -> y > 0, Fct_commod[k,:]))
			tab_fk=sortperm( Fct_commod[k,:])[1:size_fk]
			@constraint(mk,sum(e[commodity[k,1],j,tab_fk[1]] - e[j,commodity[k,1],tab_fk[1]] for j in 1:nb_nodes)-x_ikf[commodity[k,1],tab_fk[1]]==1, base_name = "init_flot_"*string(k)) #initialisation du flot à la source
			@constraint(mk,[i in 1:nb_nodes ; i!=commodity[k,1]],sum(e[i,j,tab_fk[1]] - e[j,i,tab_fk[1]]  for j in 1:nb_nodes)+x_ikf[i,tab_fk[1]]==0, base_name = "1er_flot_"*string(k)) #flot pour traiter la première fonction de puis la source
			for t in 2:size_fk
				@constraint(mk,[i in 1:nb_nodes],sum(e[i,j,tab_fk[t]] - e[j,i,tab_fk[t]] for j in 1:nb_nodes)-x_ikf[i,tab_fk[t-1]]+x_ikf[i,tab_fk[t]]==0, base_name = "mid_flow_"*string(t)*"_"*string(k)) #flot d'un traitement au suivant
			end
			@constraint(mk,[i in 1:nb_nodes ; i!=commodity[k,2]],sum(e[i,j,nb_func+1] - e[j,i,nb_func+1] for j in 1:nb_nodes)-x_ikf[i,tab_fk[end]]==0, base_name = "last_flow_"*string(k)) #flot du dernier traitement vers le puit
			@constraint(mk,sum(e[j,commodity[k,2],end] for j in 1:nb_nodes)+x_ikf[commodity[k,2],tab_fk[end]]==1, base_name = "end_flow_"*string(k)) #fin de flot sur le puit
		
			for f in 1:nb_func
				@constraint(mk,sum(x_ikf[i,f] for i in 1:nb_nodes)<=1, base_name = "concatenation_"*string(k)*"_"*string(f)) #concatenation des problemes de flot
			end

			#Contraintes de latence
			@constraint(mk, sum( sum( sum(e[i,j,f]*latency[i,j] for f in 1:nb_func+1) for j in 1:nb_nodes) for i in 1:nb_nodes) <= commodity[k,4], base_name = "latence_"*string(k))
			
			#Objective
			v=Array{Float64,2}(zeros(nb_nodes,nb_func))
			for f in 1:nb_func
				if length( findall( y -> y == f, tab_fk))>0
					v[:,f]=-commodity[k,3].*mu_3[:,f]
				end
			end
			@objective(mk, Min, - sum(sum((v[i,f] - mu_2[i,k,f])*x_ikf[i,f] for i in 1:nb_nodes) for f in 1:nb_func) )
			optimize!(mk)
			
			#Update on data 
			min_cr = min_cr + objective_value(mk)
			X_ikf[1,:,k,:]=JuMP.value.(x_ikf)
			E[1,:,:,k,:]=JuMP.value.(e)
		end
		
		#Check if optimality
		if min_cr>=-1e-7
			not_opti=false
		end
		
		
		
	end
	
	#Last solve of model to get results
	m = Model(CPLEX.Optimizer)
	n=size(X_ikf,1)
	@variable(m, lbd[1:n] >= 0)
	@constraint(m,sum(lbd[w] for w in 1:n)==1, base_name="convex")  #Convexity constraint
	#Constraint on states
	for i in 1:nb_nodes
		@constraint(m,[f in 1:nb_func],sum(lbd[w]*X_fi[w,f,i] for w in 1:n)*func_capacity[f] >=  sum( sum(lbd[w]*X_ikf[w,i,k,f] for w in 1:n) * commodity[k,3] for k in 1:nb_commodities), base_name = "n_"*string(i)) # (3) Nombre de fonctions f à placer en i
		for k in 1:nb_commodities
			@constraint(m,[f in 1:nb_func],sum(lbd[w]*X_i[w,i] for w in 1:n) >= sum(lbd[w]*X_ikf[w,i,k,f] for w in 1:n), base_name = "ouverture_"*string(i)*"_"*string(k)) # (2) ouverture du noeud en i
		end
	end
	#Objective
	@objective(m, Min, sum(sum(lbd[w]*X_i[w,i] for w in 1:n)*cout_ouverture + sum(sum(lbd[w]*X_fi[w,f,i] for w in 1:n)*func_cost[f,i] for f in 1:nb_func) for i in 1:nb_nodes) )
	optimize!(m)

	#Results
	println(solution_summary(m))
	lbd_sol=JuMP.value.(lbd)
	x_i=sum(lbd_sol[w].*X_i[w,:] for w in 1:n)
	x_fi=sum(lbd_sol[w].*X_fi[w,:,:] for w in 1:n)
	x_ikf=sum(lbd_sol[w].*X_ikf[w,:,:,:] for w in 1:n)
	e=sum(lbd_sol[w].*E[w,:,:,:,:] for w in 1:n)
	
	return x_i, x_fi, x_ikf, e
end


io, x_i, x_fi, x_ikf, e = Relaxed("test1",true)
# bprintln(x_i, x_fi, x_ikf, e)
X_i=Array{Float64,2}(zeros(1,size(x_i,1)))
X_i[1,:]=x_i
X_fi=Array{Float64,3}(zeros(1,size(x_fi,1),size(x_fi,2)))
X_fi[1,:,:]=x_fi
X_ikf=Array{Float64,4}(zeros(1,size(x_ikf,1),size(x_ikf,2),size(x_ikf,3)))
X_ikf[1,:,:,:]=x_ikf
E=Array{Float64,5}(zeros(1,size(e,1),size(e,2),size(e,3),size(e,4)))
E[1,:,:,:,:]=e

a=DW("test1", X_i,X_fi,X_ikf,E)

