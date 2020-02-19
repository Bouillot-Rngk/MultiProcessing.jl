

push!(LOAD_PATH, joinpath(@__DIR__, "src"))

using LinearAlgebra, PosDefManifold, PosDefManifoldML, CovarianceEstimation,
      Dates, Distributions, PDMats, Revise, BenchmarkTools, Plots, GR
using EEGio, FileSystem, EEGpreprocessing, System, ERPs,
      DelphiApplications, Tyler, EEGtomography, Estimators, Processdb, MPTools


Dir, dbList, estimatorList = MPTools.init()
#                                                                                   #
#                                    Skeletton                                      #
#-----------------------------------------------------------------------------------#

multiProcessP300([1,4],["SCM","TME"])
bar(sort(meanA, rev=true), ylim=(0.5, 1))
                              Testing                                         #
#-----------------------------------------------------------------------------------#