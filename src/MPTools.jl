module MPTools
using DataFrames, CSV, Plots

export plotResults,
      loadResults,
      init,
      loadDBP300,
      loadDBTVA



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


function loadDBTVA(db)
      #-----------------------------------------------------------------------------------#
      #Load a npz database using the name of the database or the index in the dbList (see MPTools)
      #corresponding to the alphabetical position in the folder
      #-----------------------------------------------------------------------------------#
      #Input :
      #     dbName::String or Int
      #Output :
      #     files::Vector{String} with N elements

      Dir, dbList, t = MPTools.init()

      Dir = homedir()*"/Documents/My Data/EEG data/npz/BI.EEG.2013-Sorted"
      sub = ["/Sujet1","/Sujet2"]
      dbList = [readdir(Dir*s) for s ∈ sub]
      filter!(e->e∉[".DS_Store","._.DS_Store"],dbList)

      if db isa String && db in dbList  #To load a specific session dbname = "1-session1"
            dbSearch = Dir*"/Sujet"*dbName[first]*"/"*dbName*"/"
      end
      if db isa Int #to load a complete base : dbName = 1 will load Base 1
            dbSearch = Dir*"/Sujet$db/Base$db/"
      end
      try
            files = loadNYdb(dbSearch)
            return files
      catch e
            println("Base de donnees inexistante")
      end

end #loadDBPTVA


function loadResults(CSVFile)
#-----------------------------------------------------------------------------------#
#Formatting the data from output_finished.csv file (result of multiProcessP300)
#-----------------------------------------------------------------------------------#
#Input :
#     CSVFile::CSVFile result of multiProcessP300
#Output :
#     splitted::Vector{DataFrame} which contained one DataFrame per covariance estimator and database
#           used in multiProcessP300

~, dbList, estimatorList = init()


tot = DataFrame(CSV.File(CSVFile))
estimList = unique(tot[!, :Method])
v = Vector{DataFrame}()
for (i, base) ∈ enumerate(dbList)
      b = tot[in.(tot.Database, Ref([base])), :]
      if isempty(b) == false push!(v,b)
      end
end

splitted = Vector{DataFrame}()
for i = 1:length(v)
      for (j,method) ∈ enumerate(estimList)
            m = v[i][in.(v[i].Method, Ref([method])), :]
            if isempty(m) == false
                  push!(splitted,m)
            end
      end
end
return splitted

end
function plotResults(CSVFile)
#-----------------------------------------------------------------------------------#
#Formatting the data from output_finished.csv file (result of multiProcessP300)
#and creating plots and saving figure for data analysis
#-----------------------------------------------------------------------------------#
#Input :
#     CSVFile::CSVFile result of multiProcessP300
#Output :
#     splitted::Vector{DataFrame} which contained one DataFrame per covariance estimator and database
#           used in multiProcessP300
#     Computing of length(splitted) plots saved into data/ folder with appropriate name


      splitted = loadResults(CSVFile)

      savef="/nethome/bouillet/Julia/MultiProcessing.jl/plot/"
      for k = 1:length(splitted)
            sorted = sort!(splitted[k], rev= true)
            meanA = splitted[k][!, :meanA]
            sdA = splitted[k][!, :sdA]
            Base = splitted[k][!, :Database][1]
            Met = splitted[k][!, :Method][1]

            Plots.bar(meanA, ylim=(0.5, 1))
            Plots.savefig(savef*Base*"_"*Met*"_regu_mean.png")
      end

      return splitted
end

function getDBList()
#-----------------------------------------------------------------------------------#
#Get the list of database in the folder "Dir" (private function)
#-----------------------------------------------------------------------------------#
      Dir = getDir()
      dbList = readdir(Dir*"/P300")
      filter!(e->e∉[".DS_Store","._.DS_Store"],dbList)
      return dbList
end

function getDir()
#-----------------------------------------------------------------------------------#
#Get the Dir directory where database are stored (must be update) (private function)
#-----------------------------------------------------------------------------------#
      return Dir = homedir()*"/Documents/My Data/EEG data/npz"
end

function getEstimatorList()
#-----------------------------------------------------------------------------------#
#Get the list of available estimator (must be update if a new estimator is available)
#(private function)
#-----------------------------------------------------------------------------------#
      return estimatorList = ["SCM","TME","nrTME","Wolf","nrTMEFilt"]
end

function init()
#-----------------------------------------------------------------------------------#
#Get useful var to start the multiprocessing
#-----------------------------------------------------------------------------------#
      Dir = getDir()
      dbList = getDBList()
      estimatorList = getEstimatorList()
      return Dir, dbList, estimatorList
end

end #module
