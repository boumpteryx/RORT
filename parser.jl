# parser

function read_instance(MyFileName::String)
  path = "./Instances/" * MyFileName
  path_Affinity = path * "_Affinity.txt"
  path_Commodity = path * "_Commodity.txt"
  path_Fct_commod = path * "_Fct_commod.txt"
  path_Functions = path * "_Functions.txt"
  path_Graph = path * "_Graph.txt"

  # Variables
  cout_ouverture = 1

  # Graph
  if isfile(path_Graph)
    myFile = open(path_Graph)
    readline(myFile)
    nb_nodes = parse(Int64, split(readline(myFile), " ")[2])
    nb_arcs = parse(Int64, split(readline(myFile), " ")[2])
    latency = Array{Int64,2}(zeros(nb_nodes,nb_nodes))
    nb_functions_per_node = Array{Int64,1}(zeros(nb_nodes))
    for i in 1:nb_arcs
      line = parse.(Int64, split(readline(myFile), " "))
      latency[line[1]+1,line[2]+1] = line[5]
      nb_functions_per_node[line[1]+1] = line[3]
      nb_functions_per_node[line[2]+1] = line[4]
    end
    close(myFile)
  end

  # Commodity
  if isfile(path_Commodity)
    myFile = open(path_Commodity)
    readline(myFile)
    nb_commodities = parse(Int64, split(readline(myFile), " ")[2])
    commodity = Array{Int64,2}(zeros(nb_commodities,4))
    for i in 1:nb_commodities
      line = parse.(Int64, split(readline(myFile), " "))
      commodity[i,1] = line[1] + 1
      commodity[i,2] = line[2] + 1
      commodity[i,3] = line[3]
      commodity[i,4] = line[4]
    end
    close(myFile)
  end

  
  # Functions
  if isfile(path_Functions)
    myFile = open(path_Functions)
    readline(myFile)
    nb_func = parse(Int64, split(readline(myFile), " ")[2])
    func_capacity = Array{Int64,1}(zeros(nb_func))
    func_cost = Array{Int64,2}(zeros(nb_func,nb_nodes))
    for i in 1:nb_func
      line = parse.(Int64, split(readline(myFile), " "))
      func_capacity[i] = line[1]
      popfirst!(line)
      for j in 1:nb_nodes
        func_cost[i,j] = line[j]
      end
    end
    close(myFile)
  end
  
  # Fct_commod
  if isfile(path_Fct_commod)
    myFile = open(path_Fct_commod)
    Fct_commod = Array{Int64,2}(zeros(nb_commodities,nb_func))
		for i in 1:nb_commodities
      line = parse.(Int64, split(readline(myFile), " "))
      cpt=1
      for f in line
        Fct_commod[i,f+1]=cpt
        cpt=cpt+1
      end
    end
    close(myFile)
  end


  # Affinity
  if isfile(path_Affinity)
    myFile = open(path_Affinity)
    exclusion = Array{Int64,2}(zeros(nb_commodities,2))
    for i in 1:nb_commodities
      line = parse.(Int64, split(readline(myFile), " "))
      if length(line) >=2
        exclusion[i,1] = line[1] + 1
        exclusion[i,2] = line[2] + 1
      end
    end
    close(myFile)
  end
  return cout_ouverture, Fct_commod, func_cost, func_capacity, nb_nodes, nb_arcs, nb_commodities, latency, nb_functions_per_node, commodity, nb_func, exclusion
end

function write_results(fileName,e,x_ikf)
  
  cout_ouverture, Fct_commod, func_cost, func_capacity, nb_nodes, nb_arcs, nb_commodities, latency, node_capacity, commodity, nb_func, exclusion = read_instance(fileName)
  
  if !isfile("./Resultats/"*fileName*".txt") 
		touch("./Resultats/"*fileName*".txt")
	end
	file = open("./Resultats/"*fileName*".txt","w")

	for k in 1:nb_commodities
    size_fk=length( findall( y -> y > 0, Fct_commod[k,:]))
		tab_fk=sortperm( Fct_commod[k,:])[1:size_fk]
		write(file,"commod "*string(k)*'\n')
		
    #write source to first function
		s=commodity[k,1]
		p= findall( y -> y == 1., x_ikf[:,k,tab_fk[1]])[1,1,1]
		sol=string(s)
		i=s
		while i!=p
			i=findall(y->y==1., e[i,:,k,tab_fk[1]])[1,1,1,1]
			sol=sol*" "*string(i)
		end
		sol=sol*'\n'
		write(file,sol)
		
    #write function to function
		for f in tab_fk[2:end]
			s=p
			p= findall( y -> y == 1., x_ikf[:,k,f])[1,1,1]
			sol=string(s)
			i=s
			while i!=p
				i=findall(y->y==1., e[i,:,k,f])[1,1,1,1]
				sol=sol*" "*string(i)
			end
			sol=sol*'\n'
			write(file,sol)
		end
    
    #write function to sink
    s=p
    p=commodity[k,2]
    sol=string(s)
    i=s
    while i!=p
      i=findall(y->y==1., e[i,:,k,end])[1,1,1,1]
      sol=sol*" "*string(i)
    end
    sol=sol*'\n'
    write(file,sol)
	end
	
	close(file)
end




