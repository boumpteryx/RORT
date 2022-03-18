# RORT
Recherche Opérationnelle dans les Réseaux de Transports


Pour ce repo nous utilisons Julia 1.6.3 et CPLEX :

Récupérer et exécuter le fichier d'installation sur le site de julia.
Installer le package JuMP qui permet de modéliser des problèmes d'optimisation en julia :
- (Linux / Mac)
Ouvrer un terminal 
Taper la commande "julia" ;
- (Windows)
 Lancer l’application julia ;
 
Taper les commandes :
> import Pkg

> Pkg.add("JuMP")

Indiquer à JuMP où est installé CPLEX en fixant dans une fenêtre julia la variable ENV["CPLEX_STUDIO_BINARIES"]. 
Exemples à adapter en fonction de votre emplacement d’installation :
- Windows
> ENV["CPLEX_STUDIO_BINARIES"] = "C:\\Program Files\\CPLEX_Studio1210\\cplex\\bin\\x86-64_win\\"

- Mac
> ENV["CPLEX_STUDIO_BINARIES"] = "/Applications/CPLEX_Studio1210/cplex/bin/x86-64_osx/"

- Linux
> ENV["CPLEX_STUDIO_BINARIES"] = "/opt/CPLEX_Studio1210/cplex/bin/x86-64_linux/"

Installer le package CPLEX : 
> Pkg.add(“CPLEX”)

Tout d'abord, lancer la commande :

> cd("C:\\Users\\Antoine\\Documents\\MPRO\\ECMA\\projet_git\\ECMA") 

avec le chemin d'accès correspondant à votre ordinateur.
Puis éxécutez le fichier main :
> include("main.jl")

C'est dans le fichier main que vous pourrez choisir quelle méthode de résolution exécuter (n'oubliez pas de spécifier le nom de l'instance dans les arguments).
