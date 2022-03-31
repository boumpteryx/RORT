using JuMP
using CPLEX
include("parser.jl")
include("relaxed.jl")
include("base.jl")
include("DW_mp_dual.jl")
include("heuristic.jl")

function func()
	dir_tab=["abilene","atlanta","dfn-bwin","dfn-gwin","di-yuan","newyork","nobel-germany","nobel-us","pdh","polska"]

	f = open("perf.txt","w")


	for dir in dir_tab
		println(Static("di-yuan/di-yuan_1/")[2])

		
		files = readdir("Instances/"*dir)
		for file in files
			println(file)
			
			println(dir)
			println(Static("di-yuan/di-yuan_1/")[2])

			
			str=dir*"/"*file

			write(f,file*" & ")
			println(Static(str, silent=true))
			t = round(@elapsed Relaxed(str, silent=true), digits=2)
			write(file,string(t)*" & ")

			t = round(@elapsed Static(str, silent=true), digits=2)
			write(file,string(t)*" & ")

			if heuristic(dir*"/"*file)
				t = round(@elapsed heuristic(str), digit=2)
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

				t = round(@elapsed DW_md(str,X_i,X_fi,X_ikf,E), digit=2)
				write(file,string(t))
			else
				write(file," \emptyset & \empty set")
			end
			write(file,'\n')
		end
	end

	close(f)
end

func()