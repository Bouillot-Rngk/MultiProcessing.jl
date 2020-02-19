module MPTools
using DataFrames, CSV, Plots

export plotResults,
      init

function plotResults(fichierCsv)
#-----------------------------------------------------------------------------------#
#Traitement des donnees recuperees lors de multiProcessP300 pour analyse
#-----------------------------------------------------------------------------------#

      ~, dbList, estimatorList = init()

      tot = DataFrame(CSV.File(fichierCsv))
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

      for i=1:length(splitted)
            meanA = splitted[i][!, :meanA]
            bar(sort(meanA, rev=true), ylim=(0.5, 1))
      end

      return splitted
end

function getDBList()
      Dir = getDir()
      dbList = readdir(Dir*"/P300")
      filter!(e->e∉[".DS_Store","._.DS_Store"],dbList)
      return dbList
end

function getDir()
      return Dir = homedir()*"/Documents/My Data/EEG data/npz"
end

function getEstimatorList()
      return estimatorList = ["SCM","TME","nrTME"]
end

function init()
      Dir = getDir()
      dbList = getDBList()
      estimatorList = getEstimatorList()
      return Dir, dbList, estimatorList
end

end #module
