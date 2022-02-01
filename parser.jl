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
    latency = Array{Int64,2}(zeros(n,n))
    nb_functions_per_node = Array{Int64,1}(zeros(nb_nodes))
    for i in [1:nb_arcs]
      line = parse(Int64, split(readline(myFile), " "))
      latency[line[1],line[2]] = line[5]
      nb_functions_per_node[line[1]] = line[3]
      nb_functions_per_node[line[2]] = line[4]
    end

  # Commodity
  if isfile(path_Commodity)
    myFile = open(path_Commodity)
    readline(myFile)
    nb_commodities = parse(Int64, split(readline(myFile), " ")[2])
    commodity = Array{Int64,2}(zeros(nb_commodities,4))
    for i in [1:nb_commodities]
      line = parse(Int64, split(readline(myFile), " "))
      commodity[i] = line
    end

  # Fct_commod
  if isfile(path_Fct_commod)
    myFile = open(path_Fct_commod)
    Fct_commod = Array{Int64,2}(zeros(nb_commodities,1))
		for i in [1:nb_commodities]
      line = parse(Int64, split(readline(myFile), " "))
      Fct_commod[i,1] = line[1]
      for j in [2:length(line)]:
        Fct_commod[i] = vcat(Fct_commod[i], [line[j] + 1])
      end
    end

  # Functions
  if isfile(path_Functions)
    myFile = open(path_Functions)
    readline(myFile)
    nb_func = parse(Int64, split(readline(myFile), " ")[2])
    func_capacity = Array{Int64,1}(zeros(nb_func))
    func_cost = Array{Int64,2}(zeros(nb_func,nb_nodes))
    for i in [1:nb_func]
      line = parse(Int64, split(readline(myFile), " "))
      func_capacity[i] = line[1]
      popfirst!(line)
      func_cost[i] = line
    end

  # Affinity
  if isfile(path_Affinity)
    myFile = open(path_Affinity)
    exclusion = Array{Int64,2}(zeros(nb_commodities,2))
    for i in [1:nb_commodities]
      line = parse(Int64, split(readline(myFile), " "))
      if length(line) >=2
        exclusion[i,1] = line[1] + 1
        exclusion[i,2] = line[2] + 1
      end
    end


  return cout_ouverture, Fct_commod, func_cost, func_capacity, nb_nodes, nb_arcs, nb_commodities, latency, nb_functions_per_node, commodity, nb_func, exclusion
end
