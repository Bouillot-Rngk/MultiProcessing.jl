module Processdb

using LinearAlgebra, PosDefManifold, PosDefManifoldML, CovarianceEstimation,
      Dates, Distributions, PDMats, Revise, BenchmarkTools, Plots, GR, DataFrames,
      CSV
using EEGio, FileSystem, EEGpreprocessing, System, ERPs,
      DelphiApplications, Tyler, EEGtomography, Estimators,MPTools


export
      processDBP300,
      loadDBP300,
      multiProcessP300



#Charger une base
function loadDBP300(dbName)
#-----------------------------------------------------------------------------------#
#Permet de charger une base de donnees rapidement en utilisant uniquement le numero
#de la base (ordre alphabetique), ou en donnant directement son nom
#-----------------------------------------------------------------------------------#
#Input :
#     dbName::String || Int
#Output :
#     files::Vector{String} with N elements
      dbList = ["BI.EEG.2012-GIPSA","BI.EEG.2013-GIPSA","BI.EEG.2014a-GIPSA",
            "BI.EEG.2015a-GIPSA","BNCI2014008","BNCI2014009","BNCI2015003","EPFLP300"]
      Dir = homedir()*"/Documents/My Data/EEG data/npz"

      if dbName isa String && dbName in dbList
            dbSearch = Dir*"/P300/"*dbName;
      end
      if dbName isa Int
            dbSearch = Dir*"/P300/"*dbList[dbName];
      end
      try
            files = loadNYdb(dbSearch)
            return files
      catch e
            println("Base de donnees inexistante");
      end
end #loadDBP300

function processDBP300(dbName)
#-----------------------------------------------------------------------------------#
#Permet de charger une base de donnees rapidement en utilisant uniquement le numero
#de la base (ordre alphabetique), ou en donnant directement son nom
#-----------------------------------------------------------------------------------#
#Input :
#     dbLoaded => Vector{String} issu de loadDBP300
#Output :
#     meanA
#     sdA


dbList = ["BI.EEG.2012-GIPSA","BI.EEG.2013-GIPSA","BI.EEG.2014a-GIPSA",
      "BI.EEG.2015a-GIPSA","BNCI2014008","BNCI2014009","BNCI2015003","EPFLP300"]
Dir = homedir()*"/Documents/My Data/EEG data/npz"


files = loadDBP300(dbName)
meanA=Vector{Float64}(undef, length(files)); sdA = similar(meanA)
println("Methode TME")
for (i, file) ∈ enumerate(files)
      o=readNY(file; bandpass=(1, 16)) # read files and create the o structure

      print(rpad("$i.", 4), rpad("sj: $(o.subject), ss: $(o.session), run $(o.run): ", 26));
      ⌚ = now()

      #Case

      Clw = TMEP300(o)
      args=(shuffle=false, tol=1e-6, verbose=false)
      cvlw = cvAcc(ENLR(Fisher), Clw, o.y; args...)

      meanA[i] = cvlw.avgAcc
      sdA[i] = cvlw.stdAcc
      println(rpad(round(meanA[i]; digits=4), 6), " (", rpad(round(sdA[i]; digits=4), 6), ") done in ", now()-⌚)
end
println("")
return meanA, sdA
end #processData

function multiProcessP300(databases,estimators)::Bool
#-----------------------------------------------------------------------------------#
#Traitements de differentes bases de donnees avec differents estimateurs.
#Recuperation des donnees dans un csv.
#-----------------------------------------------------------------------------------#

#Variables utiles /!\ Comment eviter la demultiplication des inits de ces variables ?
#Functiom getdbList ? et getDir ? => a essayer
Dir, dbList, estimatorList = MPTools.init()
donnees = DataFrame(meanA = Float64[], sdA = Float64[], Time = DateTime[],
      Method = String[], Database = String[], Sujet_Session_Run = String[])

#Verifier que tout va bien

for (k,base) ∈ enumerate(databases)
      #affichage dans le REPL pour controle

      #chargement de la base
      files = loadDBP300(base)
      meanA=Vector{Float64}(undef, length(files)); sdA = similar(meanA)

      if typeof(base)==Int print("base ",dbList[base], "  w/ ", length(files)," elements \n"); base=dbList[base]
      else print("base", base, "  w/ ", length(files)," elements \n")
      end #end if

      #Choix de la methode  => Fonctionne
      for (j,method) ∈ enumerate(estimators)
            print("Methode ",method, "\n")
            count = 0
            for (i, file) ∈ enumerate(files)
                  o=readNY(file; bandpass=(1, 16)) # read files and create the o structure
                  print(method, " ")
                  print(rpad("$i.", 4), rpad("sj: $(o.subject), ss: $(o.session), run $(o.run): ", 26));
                  ⌚ = now()

                  if method == "SCM"
                        Clw = SCMP300(o)

                  elseif method == "TME"
                        Clw = TMEP300(o)

                  elseif method == "nrTME"
                        Clw = nrTMEP300(o)
                  else print("Estimator doesn't exist"); break

                  end #switch

                  args=(shuffle=false, tol=1e-6, verbose=false)
                  cvlw = cvAcc(ENLR(Fisher), Clw, o.y; args...)

                  meanA[i] = cvlw.avgAcc
                  sdA[i] = cvlw.stdAcc
                  time = now()-⌚
                  println(rpad(round(meanA[i]; digits=4), 6), " (", rpad(round(sdA[i]; digits=4), 6), ") done in ", time)

                  data = [cvlw.avgAcc, cvlw.stdAcc, time, method, base,
                        rpad("sj: $(o.subject), ss: $(o.session), run $(o.run): ", 26) ]
                  push!(donnees,data)
                  count += 1
                  if rem(count,3) == 0
                        CSV.write("data/output.csv",donnees)
                  end

            end #for file
            #Recup des meanA et sdA pour ecriture dans csv
            #Format de csv :


      end #for method
      #bar(sort(meanA, rev=true), ylim=(0.5, 1))
      #waste(files)
end #for database
print("Finished")
CSV.write("data/output_finished.csv", donnees)
end

return bool = true

end #module
