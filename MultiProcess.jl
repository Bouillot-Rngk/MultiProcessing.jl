push!(LOAD_PATH, joinpath(@__DIR__, "src"))
#Resolving usual issues on launch
#Pkg.build(Rmath)
#Pkg.build(GLMnet)

using LinearAlgebra, PosDefManifold, PosDefManifoldML, CovarianceEstimation,
      Dates, Distributions, PDMats, Revise, BenchmarkTools, Plots, GR
using EEGio, FileSystem, EEGpreprocessing, System, ERPs,
      DelphiApplications, Tyler, EEGtomography, Estimators, Processdb, MPTools


Dir, dbList, estimatorList = MPTools.init()
#                                                                                   #
#                                    Skeletton                                      #
#-----------------------------------------------------------------------------------#

multiProcessP300([1],["Wolf"])
