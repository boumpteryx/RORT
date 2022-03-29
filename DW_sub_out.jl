using JuMP
using CPLEX
include("parser.jl")
include("relaxed.jl")
include("base.jl")

export DW

function DW(fileName :: String, X_i::Array{Float64,2},  X_fi::Array{Float64,3}, X_ikf::Array{Float64,4}, E::Array{Float64,5})
		
	#Data acquisition, model definition
	cout_ouverture, Fct_commod, func_cost, func_capacity, nb_nodes, nb_arcs, nb_commodities, latency, node_capacity, commodity, nb_func, exclusion = read_instance(fileName)
		
	#state sub-problem definition
	m0 = Model(CPLEX.Optimizer)
	set_silent(m0)
	@variable(m0, x_i[1:nb_nodes] >= 0)
	@variable(m0, x_fi[1:nb_func,1:nb_nodes]>=0)
	@constraint(m0,[i in 1:nb_nodes], x_i[i] <= 1) # Relaxed binary constraint
	@constraint(m0,[i in 1:nb_nodes], sum(x_fi[f,i] for f in 1:nb_func) <= node_capacity[i], base_name = "cap" ) # (4) capacité de noeud

	#flow sub problem definition
	
	m_k=Dict{Int,JuMP.Model}()
	m_k[0]=m0
	
	for k in 1:nb_commodities

		#model definition
		mk = Model(CPLEX.Optimizer)
		set_silent(mk)
		@variable(mk, x_ikf[1:nb_nodes,1:nb_func]>=0, base_name="x_ikf_"*string(k)) 
		@variable(mk, e[1:nb_nodes,1:nb_nodes,1:nb_func+1] >= 0,base_name="e_"*string(k))

		for i in 1:nb_nodes
			@constraint(mk,[f in 1:nb_func], x_ikf[i,f] <= 1) # Relaxed binary constraint
			for j in 1:nb_nodes
				@constraint(mk,[f in 1:nb_func+1], e[i,j,f] <= 1) # Relaxed binary constraint
			end
		end

		#Constraints of exclusion and possible arc
		for i in 1:nb_nodes
			if exclusion[k,1]>0
				@constraint(mk,sum(x_ikf[i,w] for w in exclusion[k,:]) <= 1, base_name = "exclusion_"*string(k)*"_"*string(i)) #exclusion
			end				
			for j in 1:nb_nodes
				if latency[i,j] == 0
					@constraint(mk,[f in 1:nb_func+1],e[i,j,f]==0, base_name = "no_arc_"*string(i)*"_"*string(j)) #absence d'arcs
				end
			end
		end

		#Contrainte de flots successifs 
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
		
		m_k[k]=mk
	end

	
	#Column generation
	not_opti=true
	n=0
	while not_opti 
		
		#Probleme maitre
		m = Model(CPLEX.Optimizer)
		set_silent(m)
		n=size(X_ikf,1)
		@variable(m, lbd[1:n] >= 0)
		@constraint(m,sum(lbd[w] for w in 1:n)>=1, base_name="convex") #Convexity constraint
		for i in 1:nb_nodes
			@constraint(m,[f in 1:nb_func],sum(lbd[w]*X_fi[w,f,i] for w in 1:n)*func_capacity[f] >=  sum( sum(lbd[w]*X_ikf[w,i,k,f] for w in 1:n) * commodity[k,3] for k in 1:nb_commodities), base_name = "n_"*string(i)) # (3) Nombre de fonctions f à placer en i
			for k in 1:nb_commodities
				@constraint(m,[f in 1:nb_func],sum(lbd[w]*X_i[w,i] for w in 1:n) >= sum(lbd[w]*X_ikf[w,i,k,f] for w in 1:n), base_name = "ouverture_"*string(i)*"_"*string(k)) # (2) ouverture du noeud en i
			end
		end
		@objective(m, Min, sum(sum(lbd[w]*X_i[w,i] for w in 1:n)*cout_ouverture + sum(sum(lbd[w]*X_fi[w,f,i] for w in 1:n)*func_cost[f,i] for f in 1:nb_func) for i in 1:nb_nodes) )
		optimize!(m)
		
		#Get dual values and solutions
		#println(JuMP.value.(lbd))
		eta=dual(constraint_by_name(m, "convex"))
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

		
		#setup for next column to add on first index
		min_cr=-eta #Minimum ecart complementaire
		X_i=vcat(X_i,X_i[1:1,:])
		X_fi=vcat(X_fi,X_fi[1:1,:,:])
		X_ikf=vcat(X_ikf,X_ikf[1:1,:,:,:])
		E=vcat(E,E[1:1,:,:,:,:])
		
		#state sub-problem resolution
		@objective(m_k[0], Min, sum(x_i[i]*(cout_ouverture-sum(sum(mu_2[i,k,f] for f in 1:nb_func) for k in 1:nb_commodities))+sum(x_fi[f,i]*(func_cost[f,i]-mu_3[i,f]*func_capacity[f]) for f in 1:nb_func)  for i in 1:nb_nodes) )
		optimize!(m_k[0])
		min_cr=min_cr+objective_value(m_k[0])
		X_i[end,:]=JuMP.value.(x_i)
		X_fi[end,:,:]=JuMP.value.(x_fi)
		
		#flow sub-problem resolution
		for k in 1:nb_commodities
			x_ikf=object_dictionary(m_k[k])[:x_ikf]
			e=object_dictionary(m_k[k])[:e]
			@objective(m_k[k], Min, - sum(sum( (-mu_3[i,f]*commodity[k,3] - mu_2[i,k,f])*x_ikf[i,f] for i in 1:nb_nodes) for f in 1:nb_func) )
			optimize!(m_k[k])
			min_cr = min_cr + objective_value(m_k[k])
			X_ikf[end,:,k,:]=JuMP.value.(x_ikf)
			E[end,:,:,k,:]=JuMP.value.(e)		
		end
		
		#Check optimality
		if min_cr>=-1e-7
			not_opti=false
		end
	end
		
	#Last solve of model to get results
	m = Model(CPLEX.Optimizer)
	n=size(X_ikf,1)
	@variable(m, lbd[1:n] >= 0)
	@constraint(m, lbd[1] == 0)

	@constraint(m,sum(lbd[w] for w in 1:n)==1, base_name="convex")  #Convexity constraint
	for i in 1:nb_nodes
		@constraint(m,[f in 1:nb_func],sum(lbd[w]*X_fi[w,f,i] for w in 1:n)*func_capacity[f] >=  sum( sum(lbd[w]*X_ikf[w,i,k,f] for w in 1:n) * commodity[k,3] for k in 1:nb_commodities), base_name = "n_"*string(i)) # (3) Nombre de fonctions f à placer en i
		for k in 1:nb_commodities
			@constraint(m,[f in 1:nb_func],sum(lbd[w]*X_i[w,i] for w in 1:n) >= sum(lbd[w]*X_ikf[w,i,k,f] for w in 1:n), base_name = "ouverture_"*string(i)*"_"*string(k)) # (2) ouverture du noeud en i
		end
	end
	@objective(m, Min, sum(sum(lbd[w]*X_i[w,i] for w in 1:n)*cout_ouverture + sum(sum(lbd[w]*X_fi[w,f,i] for w in 1:n)*func_cost[f,i] for f in 1:nb_func) for i in 1:nb_nodes) )
	optimize!(m)
	println(JuMP.value.(lbd))
	mu_2_sol=Array{Float64,3}(zeros(nb_nodes,nb_commodities,nb_func))
	mu_3_sol=Array{Float64,2}(zeros(nb_nodes,nb_func))
	for i in 1:nb_nodes
		for f in 1:nb_func
			mu_3_sol[i,f]=dual(constraint_by_name(m,"n_"*string(i)*"["*string(f)*"]"))
			for k in 1:nb_commodities
				mu_2_sol[i,k,f]=dual(constraint_by_name(m,"ouverture_"*string(i)*"_"*string(k)*"["*string(f)*"]"))
			end
		end
	end		
	println(solution_summary(m))

	#with dual
	m = Model(CPLEX.Optimizer)
	set_silent(m)
	@variable(m,mu_2[1:nb_nodes,1:nb_commodities,1:nb_func]>=0)
	@variable(m,mu_3[1:nb_nodes,1:nb_func]>=0)
	#@constraint(m, [ i in 1:nb_nodes*nb_func], mu_3[i] == mu_3_sol[i])
	#@constraint(m, [ i in 1:nb_nodes*nb_func*nb_commodities], mu_2[i] == mu_2_sol[i])
	@constraint(m, [ i in 1:nb_nodes*nb_func], mu_3[i] <= 30000)
	@constraint(m, [ i in 1:nb_nodes*nb_func*nb_commodities], mu_2[i] <= 30000)

	@variable(m,eta)
	for idx in 2:size(X_ikf,1)
		@constraint(m, sum(cout_ouverture*X_i[idx,i]+sum(func_cost[f,i]*X_fi[idx,f,i] for f in 1:nb_func) for i in 1:nb_nodes) - sum( X_i[idx,i]*sum(sum(mu_2[i,k,f] for k in 1:nb_commodities) for f in 1:nb_func) + sum(func_capacity[f]*mu_3[i,f]*X_fi[idx,f,i] for f in 1:nb_func) for i in 1:nb_nodes) - sum(sum(sum( (-commodity[k,3]*mu_3[i,f]-mu_2[i,k,f])*X_ikf[idx,i,k,f] for f in 1:nb_func) for k in 1:nb_commodities) for i in 1:nb_nodes) - eta >=0, base_name="con_"*string(idx))
		#@constraint(m, sum(X_i[idx,i]*(cout_ouverture-sum(sum(mu_2[i,k,f] for f in 1:nb_func) for k in 1:nb_commodities))+sum(X_fi[idx,f,i]*(func_cost[f,i]-mu_3[i,f]*func_capacity[f]) for f in 1:nb_func)  for i in 1:nb_nodes) - sum( sum(sum( (-mu_3[i,f]*commodity[k,3] - mu_2[i,k,f])*X_ikf[idx,i,k,f] for i in 1:nb_nodes) for f in 1:nb_func) for k in 1:nb_commodities) - eta >=0, base_name="con_"*string(idx))
	end
	@objective(m, Max, eta )
	optimize!(m)
	for idx in 2:size(X_ikf,1)
		print(dual(constraint_by_name(m,"con_"*string(idx))))
		print(" ")
	end
	println(" ")

	
	#Results
	println(solution_summary(m))
	lbd_sol=JuMP.value.(lbd)
	x_i=sum(lbd_sol[w].*X_i[w,:] for w in 1:n)
	x_fi=sum(lbd_sol[w].*X_fi[w,:,:] for w in 1:n)
	x_ikf=sum(lbd_sol[w].*X_ikf[w,:,:,:] for w in 1:n)
	e=sum(lbd_sol[w].*E[w,:,:,:,:] for w in 1:n)
	
	return x_i, x_fi, x_ikf, e
end

io, x_i, x_fi, x_ikf, e = Static("di-yuan/di-yuan_1/",true)
X_i=Array{Float64,2}(zeros(1,size(x_i,1)))
X_i[1,:]=x_i
X_fi=Array{Float64,3}(zeros(1,size(x_fi,1),size(x_fi,2)))
X_fi[1,:,:]=x_fi
X_ikf=Array{Float64,4}(zeros(1,size(x_ikf,1),size(x_ikf,2),size(x_ikf,3)))
X_ikf[1,:,:,:]=x_ikf
E=Array{Float64,5}(zeros(1,size(e,1),size(e,2),size(e,3),size(e,4)))
E[1,:,:,:,:]=e


x_i, x_fi, x_ikf, e=DW("di-yuan/di-yuan_1/", X_i,X_fi,X_ikf,E)
println("fin")
#x_i, x_fi, x_ikf, e=Relaxed("di-yuan/di-yuan_1/")




