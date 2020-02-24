push!(LOAD_PATH, joinpath(@__DIR__, "src"))
push!(LOAD_PATH, joinpath(@__DIR__, "data"))
push!(LOAD_PATH,"/nethome/bouillet/.julia/packages/")


using Plots, GR
using MPTools, CSV, DataFrames


data = plotResults("data/output_finished-ssh.csv")


data1 = loadResults("data/output_B1_norm.csv")
data2 = loadResults("output_B1_nonorm.csv")
sorted1 = sort!(data1[k], rev= true)

sdA = data1[1][!, :sdA]

sorted2 = sort!(data2[1], rev= true)
meanA2 = sorted2[!, :meanA]
sdA = data2[k][!, :sdA]

xs =1:length(meanA)

      k=2
      sorted1 = sort!(data1[k], rev= true)
      meanA1 = sorted1[!, :meanA]
      sorted2 = sort!(data2[1], rev= true)
      meanA2 = sorted2[!, :meanA]

      diff = meanA1-meanA2

      Plots.bar(diff,ylim=(-0.1,0.1))
#Plots.savefig("/nethome/bouillet/Julia/MultiProcessing.jl/plot/GIPSA2012_SCM_rege4_mean.png")
