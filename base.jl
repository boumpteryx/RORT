using JuMP
using CPLEX

export Static


function Static(n :: Int64, s :: Int64, t :: Int64,  S :: Int64,  d1 :: Int64, d2 :: Int64, p :: Vector{Int64}, ph :: Vector{Int64}, A :: Array{Int64}, d :: Array{Int64}, D :: Array{Float64})

	m = Model(CPLEX.Optimizer)

	# Var
	@variable(m, x[1:n,1:n], Bin)
	
	#Constraint on flow
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
	@objective(m, Min, sum(x[i]*d[i] for i in 1:n*n) )

	optimize!(m)
	println(solution_summary(m))

	vX = JuMP.value.(x)

	status = termination_status(m)
	isOptimal = status == MOI.OPTIMAL
		
	return isOptimal
end




