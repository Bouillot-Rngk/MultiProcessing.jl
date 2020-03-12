push!(LOAD_PATH, joinpath(@__DIR__, "src"))
push!(LOAD_PATH, joinpath(@__DIR__, "/Julia/Multiprocessing.jl"))
#Resolving usual issues on launch
#Pkg.build(Rmath)
#Pkg.build(GLMnet)

using PosDefManifoldML, CovarianceEstimation, PDMats, Revise, BenchmarkTools, Plots, GR, EEGio, EEGpreprocessing, System, ERPs, DelphiApplications, Tyler,
    EEGtomography, Estimators, Processdb, MPTools, Diagonalizations


Dir, dbList, estimatorList = MPTools.init()
#                                                                                   #
#                                    Skeletton                                      #
#-----------------------------------------------------------------------------------#

#multiProcessCSV([6,7,8,9],["SCM","nrTME","TME","Wolf"],"testBNCI.csv",1)
#multiProcessCSV([3],["SCM"],"ssh1.csv",1)
#multiProcessCSV([3,4,5],["nrTME","nrTMEFilt"],"2013-2015Filter.csv",1)
