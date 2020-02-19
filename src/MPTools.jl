module MPTools
using DataFrames, CSV, Plots

export plotResults,
      init

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
#     To implement => computing of length(splitted) plots saved into data/ folder with appropriate name


      ~, dbList, estimatorList = init()

      tot = DataFrame(CSV.File(CSVFile))
      v = Vector{DataFrame}()
      for (i, base) ∈ enumerate(dbList)
            b = tot[in.(tot.Database, Ref([base])), :]
            if isempty(b) == false push!(v,b)
            end
      end

      splitted = Vector{DataFrame}()
      for i = 1:length(v)
            for (j,method) ∈ enumerate(estimatorList)
                  m = v[i][in.(v[i].Method, Ref([method])), :]
                  if isempty(m) == false
                        push!(splitted,m)
                  end
            end
      end

      #WIP : Computing and saving plots
      for i=1:length(splitted)
            meanA = splitted[i][!, :meanA]
            bar(sort(meanA, rev=true), ylim=(0.5, 1))
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
      return estimatorList = ["SCM","TME","nrTME"]
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
