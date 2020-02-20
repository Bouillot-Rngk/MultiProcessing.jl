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



function loadDBP300(dbName)
#-----------------------------------------------------------------------------------#
#Load a npz database using the name of the database or the index in the dbList (see MPTools)
#corresponding to the alphabetical position in the folder
#-----------------------------------------------------------------------------------#
#Input :
#     dbName::String or Int
#Output :
#     files::Vector{String} with N elements

      Dir, dbList, t = MPTools.init()

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

function processDBP300(dbName) #obsolete
#-----------------------------------------------------------------------------------#
#Process 1 database with one specific method that need to be written in Estimators.jl
#and function must be modified according to the wanted covariance estimator
#-----------------------------------------------------------------------------------#
#Input :
#     dbName::String or Int
#Output :
#     meanA  mean computed by crossvalidatin
#     sdA    standart deviation computed by cross validation


      Dir, dbList, t = MPTools.init()


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
#Processing of length(databases) databases with length(estimators) covariance estimator
#The if-elseif must be update when a new estimator function is available
#Return a csv file with all meanA, sdA and time computed whilst crossvalidation that
#can be read and treated with MPTools.plotResults() function
#-----------------------------------------------------------------------------------#
#Input :
#     databases::Vector{String} containing names or index of to-be-processed DB
#     estimators::Vector{String} comtaining names of covariance estimator to be tested
#             The list of DB and estimators available is described in MPTools.jl
#Output :
#     bool::Bool = true if processing is successfully complete
#     output_finished.csv::CSV File stored in data/ folder
#     output.csv in data/ folder is a backup in case of any kind of crash during computing

Dir, dbList, estimatorList = MPTools.init()
donnees = DataFrame(meanA = Float64[], sdA = Float64[], Time = DateTime[],
      Method = String[], Database = String[], Sujet_Session_Run = String[])

#Maybe check if inputs are correct before launching code ?

for (k,base) ∈ enumerate(databases)
      #loading of 1 database
      files = loadDBP300(base)

      #Memory allocation
      meanA=Vector{Float64}(undef, length(files)); sdA = similar(meanA)

      #display in REPL for control => used base and number of elements
      if typeof(base)==Int print("base ",dbList[base], "  w/ ", length(files)," elements \n"); base=dbList[base]
      else print("base", base, "  w/ ", length(files)," elements \n")
      end #end if

      #Choice of covariance estimator
      for (j,method) ∈ enumerate(estimators)
            #REPL Display for control
            print("Methode ",method, "\n")
            count = 0 #var used for output.csv (back up csv)

            #Data processing
            for (i, file) ∈ enumerate(files)

                  o=readNY(file; bandpass=(1, 16)) # read files and create the o structure
                  print(method, " ")
                  print(rpad("$i.", 4), "/",length(files), " ", rpad("sj: $(o.subject), ss: $(o.session), run $(o.run): ", 26));
                  ⌚ = now()

                  #Choice of covariance estimator => is there any way to make it more flexible/optimized  ?

                  if method == "SCM"
                        Clw = SCMP300(o)

                  elseif method == "TME"
                        Clw = TMEP300(o)

                  elseif method == "nrTME"
                        Clw = nrTMEP300(o)

                  elseif method == "Wolf"
                        Clw = Wolf(o)

                  else print("Estimator doesn't exist"); break

                  end #switch-case


                  #beginning of crossvalidation
                  args=(shuffle=false, tol=1e-6, verbose=false)
                  cvlw = cvAcc(ENLR(Fisher), Clw, o.y; args...)

                  meanA[i] = cvlw.avgAcc
                  sdA[i] = cvlw.stdAcc


                  time = now()-⌚
                  #REPL Display for control
                  println(rpad(round(meanA[i]; digits=4), 6), " (", rpad(round(sdA[i]; digits=4), 6), ") done in ", time)

                  #Datastorage
                  data = [cvlw.avgAcc, cvlw.stdAcc, time, method, base,
                        rpad("sj: $(o.subject), ss: $(o.session), run $(o.run): ", 26) ]

                  push!(donnees,data)
                  count += 1
                  #Back up csv
                  if rem(count,3) == 0 || count==length(files)
                        CSV.write("Julia/MultiProcessing.jl/data/output-ssh.csv",donnees)
                  end

            end #for file
      end #for method
end #for database

print("Finished")
#Data recuperation (readable with MPTools.plotResults())
CSV.write("Julia/MultiProcessing.jl/data/output_finished-ssh.csv", donnees)
end



end #module
