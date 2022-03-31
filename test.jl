using JuMP
using CPLEX
include("relaxed.jl")
include("base.jl")
include("DW_mp_dual.jl")
include("heuristic.jl")

function test(str)
	Relaxed("di-yuan/di-yuan_1/")
	t=time()
	Relaxed(str)
	t = round(time()-t, digits=2)
	println(t)


	
	return 0
end

function func()
	dir_tab=["pdh","di-yuan","abilene","atlanta","dfn-bwin","dfn-gwin","newyork","nobel-germany","nobel-us","polska"]

	for dir in dir_tab
		
		files = readdir("Instances/"*dir)
		for file in files
			
			f = open("perf.txt","w")
			
			str=dir*"/"*file*"/"
			
			println(str)

			write(f,file*" & ")
			
			t = time()
			val1 = Relaxed(str)[1]
			t = round(time()-t, digits=2)
			write(file,string(t)*" & ")

			t = time()
			val2 = Static(str)[1]
			t = round(time()-t, digits=2)
			write(file,string(t)*" & ")

			if heuristic(str)
				
				t = time()
				heuristic(str)
				t = round(time()-t, digits=2)
				write(file,string(t)*" & ")

				ov, x_i, x_fi, x_ikf, e = heuristic(str,res=true)
				X_i=Array{Float64,2}(zeros(1,size(x_i,1)))
				X_i[1,:]=x_i
				X_fi=Array{Float64,3}(zeros(1,size(x_fi,1),size(x_fi,2)))
				X_fi[1,:,:]=x_fi
				X_ikf=Array{Float64,4}(zeros(1,size(x_ikf,1),size(x_ikf,2),size(x_ikf,3)))
				X_ikf[1,:,:,:]=x_ikf
				E=Array{Float64,5}(zeros(1,size(e,1),size(e,2),size(e,3),size(e,4)))
				E[1,:,:,:,:]=e

				t = time()
				DW_md(str,X_i,X_fi,X_ikf,E)
				t = round(time()-t, digits=2)
				write(file,string(t)*" & ")
			else
				write(file," \emptyset & \empty set & ")
			end
			write(file,string(val2)*" & "*string(val1)*'\n')
			
			close(f)

		end
	end

end

func()