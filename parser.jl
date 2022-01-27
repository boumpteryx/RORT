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
  commodity =

  # Graph
  if isfile(path_Graph)
    myFile = open(path_Graph)
    readline(myFile)
    nb_nodes = parse(Int64, split(readline(myFile), " ")[2])
    nb_arcs = parse(Int64, split(readline(myFile), " ")[2])
    for i in [1:nb_arcs]
      line = parse(Int64, split(readline(myFile), " "))

    end

  # Commodity
  if isfile(path_Commodity)
    myFile = open(path_Commodity)
    readline(myFile)
    nb_commodities = arse(Int64, split(readline(myFile), " ")[2])
    for i in [1:nb_commodities]

    end

  # Fct_commod
  if isfile(path_Fct_commod)
    myFile = open(path_Fct_commod)
    data = readlines(myFile)
		for datum in data

    end

  # Functions
  if isfile(path_Functions)
    myFile = open(path_Functions)

  # Affinity
  if isfile(path_Affinity)
    myFile = open(path_Affinity)

  return cout_ouverture,
end
