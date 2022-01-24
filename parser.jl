# parser

function read_instance(MyFileName::String)
  path = "./Instances/" * MyFileName
  # Si le fichier path existe
  if isfile(path) # exemple : "./Instances/grille2x3_Affinity.txt"
    # L’ouvrir
    myFile = open(path)

    # TODO

end
