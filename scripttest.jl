push!(LOAD_PATH, joinpath(@__DIR__, "src"))
push!(LOAD_PATH, joinpath(@__DIR__, "data"))
push!(LOAD_PATH,"/nethome/bouillet/.julia/packages/")


using Plots, GR
using MPTools, CSV, DataFrames
using Processdb, EEGio

df =plotResults("data/2013-2015Filter.csv")



k=2

df[k][!, :Method] = ["TME_detnorm" for i=1:length(df[k][!, :Method])]
CSV.write("data/2014-2015Wolfnonorm-TME-norm-SCM-norm.csv",df)
tot = DataFrame(CSV.File("data/2014-2015Wolfnonorm-TME-norm-SCM-norm.csv"))



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

####

using NPZ, YAML, FileSystem, EEGpreprocessing, DSP
using Plots, LinearAlgebra, Processdb, MPTools, PosDefManifold
using EEGio, Diagonalizations,CovarianceEstimation


base = 2;
Dir, dbList, estimatorList = MPTools.init();
files = Processdb.loadDBP300(base);

#=  meanA=Vector{Float64}(undef, length(files)); sdA = similar(meanA);

  #display in REPL for control => used base and number of elements
  if typeof(base)==Int print("base ",dbList[base], "  w/ ", length(files)," elements \n"); base=dbList[base];
  else print("base", base, "  w/ ", length(files)," elements \n");
  end #end if
=#

file = files[1];
o=readNY(file; bandpass=(1, 16))
typeof(o.X)
size(o.X)
X = copy(o.X)
w=[[1/(norm(o.X[o.cstim[i][j]+o.offset:o.cstim[i][j]+o.offset+o.wl-1,:])^2) for j=1:length(o.cstim[i])] for i=1:o.nc]
Y=mean(o.X, o.wl, o.cstim; weights=w)[2]
Y=Y*eigvecs(cov(SimpleCovariance(), Y))[:, o.ne-3:o.ne];
wX = whitening([X Y];eVar=0.99)
Clw=ℍVector([Hermitian(cov(SimpleCovariance(), [X Y])) for X ∈ o.trials])
Cl=copy(Clw)
for i = 1:length(Clw)
      Cl[i]=Hermitian(wX.F'*Clw[i]*wX.F)
end
length(Clw)
C=Clw[1]
D=Cl[1]
Cmax=maximum(abs.(C));
h1 = heatmap(C, clim=(-Cmax, Cmax), yflip=true, c=:bluesreds, title="C");

h2 = heatmap(D, clim=(0, 1), yflip=true, c=:amp, title="F'*C*F");
📈=plot(h1,h2, size=(700, 300))


X = o.X
n, t = size(X)

iter, 😋, α, β, nt⁻¹ = 1, false, 0.0, 0.0, n / t
@inbounds x² = [x⋅x for x ∈ eachcol(X')]

Σx = sum(x²)

BPFilt = digitalfilter(Bandpass(1/(o.sr/2), 8/(o.sr/2)), Butterworth(2))
Xf = filtfilt(BPFilt, o.X)

x²f = [x⋅x for x ∈ eachcol(Xf')]
Σxf = sum(x²f)

x²fΣ = x²f*Σx/Σxf


Plots.plot(sort(x²fΣ))
Plots.plot!(sort(x²))

x
#=
 # read files and create the o structure

using NPZ, YAML, FileSystem, EEGpreprocessing, DSP;
filename = file;
bandpass = (1,16);
data = npzread(splitext(filename)[1]*".npz") # read data file;
info = YAML.load(open(splitext(filename)[1]*".yml")) # read info file;

      sr      = info["acquisition"]["samplingrate"]
      stim    = data["stim"]                  # stimulations
      (ns, ne)= size(data["data"])            # of sample, # of electrodes)
      os      = info["stim"]["offset"]        # offset for trial starting sample
      wl      = info["stim"]["windowlength"]  # trial duration
      nc      = info["stim"]["nclasses"]      # of classes
      BPfilter = digitalfilter(Bandpass(first(bandpass)/(sr/2), last(bandpass)/(sr/2)), Butterworth(2))
      X        = filtfilt(BPfilter, data["data"])


cstim=[[i+os for i in eachindex(stim) if stim[i]==j && i+os+wl<=ns] for j=1:nc]
trials=[X[cstim[i][j]:cstim[i][j]+wl-1,:] for i=1:nc for j=1:length(cstim[i])]
trials=nothing

EEG(
  info["id"],
  info["acquisition"],
  info["documentation"],
  info["formatversion"],

  info["id"]["database"],
  info["id"]["subject"],
  info["id"]["session"],
  info["id"]["run"],
  info["acquisition"]["sensors"],
  sr,
  ne,
  ns,
  wl,
  os, # trials offset
  nc,
  collect(keys(info["stim"]["labels"])), # clabels
  stim,
  cstim,
  [i for i=1:nc for j=1:length(cstim[i])], # y: all labels
  X, # whole EEG recording
  trials # all trials, by class
)

#mean = Statistics.mean(abs.(o.X))
#print(mean, " \n")
print(i+x-1, "/",length(files), " ", rpad("sj: $(o.subject), ss: $(o.session), run $(o.run): ", 26)," \n");


####
push!(LOAD_PATH, joinpath(@__DIR__, "src"))
using GLM, LinearAlgebra, DataFrames
using Processdb, MPTools, EEGio
using Plots, Statistics

base = 5;
Dir, dbList, estimatorList = MPTools.init();
files = Processdb.loadDBP300(base)

file = files[5];
o=readNY(file; bandpass=(1, 16));


gfp=[x⋅x for x ∈ eachrow(o.X)]

#Calcul du gfp + tri dans un DataFrame pour conserver l'ordre
lsgfp=log10.(gfp)
Plots.plot(lsgfp)
x = 1:length(lsgfp)
data=DataFrame(X = x, GFP = lsgfp, Deriv = missing)
datasorted = sort!(data, [:GFP, :X])
Plots.plot(datasorted[!, :GFP])


#Calcul de la derivée
derivative = Vector{Float64}(undef,length(lsgfp))
for x0 = 31:length(lsgfp)-30
  moy=Vector{Float64}(undef,29)
  for wl = 2:30
      moy[wl-1] = (datasorted[!, :GFP][x0+wl] - datasorted[!, :GFP][x0-wl])/(2*wl)
  end
  mean = Statistics.mean(moy)
  derivative[x0] = mean
end
for i=1:length(derivative)
  if derivative[i]<0.0000000000001 derivative[i] = 0 ; end
end
deriv = Statistics.mean(derivative)
Scalederivative = derivative/deriv
Plots.plot(Scalederivative)
for i=length(derivative)-50:length(derivative)
  Scalederivative[i] = 30
end
datasorted.Deriv = Scalederivative
dataArtifact = datasorted[1000:end,  :]

dataArtifact = dataArtifact[dataArtifact.Deriv .> 15, :]
Plots.plot(dataArtifact[!, :GFP])
indexArtif = copy(dataArtifact)
indexArtif = sort!(indexArtif, :X)
Plots.plot(indexArtif[!, :GFP])

stimArtifact = deleteat!(o.stim, indexArtif[!, :X])
cstim=[[i+o.os for i in eachindex(stimArtifact) if stimArtifact[i]==j && i+o.os+o.wl<=o.ns ] for j=1:o.nc]
