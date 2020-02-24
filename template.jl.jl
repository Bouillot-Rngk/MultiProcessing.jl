# This unit computes classification accuracy on a database
# and can be used to compare several cov matrix estimators,
# including shrinkage and tyler M-estimator.

DBDir = homedir()*"/Documents/My Data/EEG data/npz"

using LinearAlgebra, PosDefManifold, PosDefManifoldML, CovarianceEstimation,
      Dates, Distributions, PDMats, Revise, BenchmarkTools, Plots, GR

push!(LOAD_PATH, joinpath(@__DIR__, "Modules"))


using EEGio, FileSystem, EEGpreprocessing, System, ERPs,
      DelphiApplications, Tyler, EEGtomography

lse, nlse = LinearShrinkage, AnalyticalNonlinearShrinkage

function classify(o; shuffle=false)
    # multivariate regression ERP mean with data-driven weights
    w=[[1/(norm(o.X[o.cstim[i][j]+o.offset:o.cstim[i][j]+o.offset+o.wl-1,:])^2) for j=1:length(o.cstim[i])] for i=1:o.nc]

	Y=mean(o.X, o.wl, o.cstim; weights=w)[2]

    # PCA to keep only 4 components
    Y=Y*eigvecs(cov(SimpleCovariance(), Y))[:, o.ne-3:o.ne]

    # encoding

	# sample covariance matrix
    𝐂lw=ℍVector([ℍ(cov(SimpleCovariance(), [X Y])) for X ∈ o.trials])

	# linear shrinkage
    #𝐂lw=ℍVector([ℍ(cov(lse(ConstantCorrelation()), [X Y])) for X ∈ o.trials])

	# non-linear estimator (gives numerical problems)
	#𝐂lw=ℍVector([ℍ(cov(nlse(), [X Y])) for X ∈ o.trials])

	# TME estimator
	 𝐂lw=ℍVector([ℍ(Mest([X Y]')) for X ∈ o.trials])

	# nrTME estimator
	# 𝐂lw=ℍVector([ℍ(nrMest([X Y]'; reg=:LW)) for X ∈ o.trials])

	# regularization (necessary for non-linear shrinkage)
	#R=Hermitian(Matrix{eltype(𝐂lw[1])}(I, size(𝐂lw[1]))*0.001)
	#for C in 𝐂lw C+=R end

	# trace normalization
	#for i=1:length(𝐂lw) 𝐂lw[i]=𝐂lw[i]/tr(𝐂lw[i]) end

	# det normalization
	for i=1:length(𝐂lw) 𝐂lw[i]=det1(𝐂lw[i]) end

    # classification
    args=(shuffle=shuffle, tol=1e-6, verbose=false)
    cvlw = cvAcc(ENLR(Fisher), 𝐂lw, o.y; args...)

    return cvlw.avgAcc, cvlw.stdAcc
end


function processDB(dbName, Dir, paradigm)
    # find out the directory where the database is
    dbDir=DBDir*"/"*paradigm*"/"*dbName

    # get the complete list of file paths
    files=loadNYdb(dbDir)

    # get memory for storing results
    meanA=Vector{Float64}(undef, length(files)); sdA = similar(meanA)

    # read files, do computations and waste them, one by one
    println("\nProcessing DB "*dbName*"; # files: $(length(files))")
    for (i, file) ∈ enumerate(files)

        o=readNY(file; bandpass=(1, 16)) # read files and create the o structure

        print(rpad("$i.", 4), rpad("sj: $(o.subject), ss: $(o.session), run $(o.run): ", 26));
        ⌚ = now()
        meanA[i], sdA[i] = classify(o; shuffle=false)
        println(rpad(round(meanA[i]; digits=4), 6), " (", rpad(round(sdA[i]; digits=4), 6), ") done in ", now()-⌚)

        waste(o) # release memory
    end
    println("")
    return meanA, sdA
end

# Run Classification
using Plots
paradigm = "P300"
dbName = "BI.EEG.2013-Sorted"
dbDir=DBDir*"/"*paradigm*"/"*dbName
DirCS = dbDir*"/Sujet1/Base1"
files = loadNYdb(DirCS)
files = loadNYdbCS(dbDir)
meanA=Vector{Float64}(undef, length(files)); sdA = similar(meanA)
for (i, file) ∈ enumerate(files)

	o=readNY(file; bandpass=(1, 16)) # read files and create the o structure

	print(rpad("$i.", 4), rpad("sj: $(o.subject), ss: $(o.session), run $(o.run): ", 26));
	⌚ = now()
	meanA[i], sdA[i] = classify(o)
	println(rpad(round(meanA[i]; digits=4), 6), " (", rpad(round(sdA[i]; digits=4), 6), ") done in ", now()-⌚)

	waste(o) # release memory
end
println("")
return meanA, sdA


a, b = processDB("BI.EEG.2012-GIPSA", DBDir, "P300")
bar(sort(a, rev=true), ylim=(0.5, 1))

x

## Tools
#=

# get the data of one file
o=readNY(files[10]; bandpass=(1, 16))

# weighted multivariate regression ERP mean estimation
Z=mean(o.X, o.wl, o.cstim, true; weights=:a)
eegPlot(Z[2], o.sensors)

# send data to CSTP
erpPlot(o.X, o.sr, o.wl, o.offset, 16, o.sensors, o.stim, o.clabels)

# arithmetic mean
Y=fVec(mean, 𝕄Vector([o.trials[i] for i in eachindex(o.trials) if o.y[i]==2]))
eegPlot(Y, o.sensors)

# easy way to get the arithmetic mean
Z=mean(o.X, o.wl, o.cstim, false)


sr, stim, (ns, ne) = info["acquisition"]["samplingrate"], data["stim"], size(data["data"])

# offset for trial starting sample, trial duration, # of classes
offset, td, z = info["stim"]["offset"], info["stim"]["windowlength"], info["stim"]["nclasses"]

ICoNargs=(msg="\x1b[33m"*"Please wait while the Wave Editor runs", )

# PRE-PROCESSING
################
BPfilter = digitalfilter(Bandpass(1/(sr/2), 16/(sr/2)), Butterworth(2))
# Winsor standardize BP-filtered EEG and convert to Float64 if necessary
EEG	= standardizeEEG(filtfilt(BPfilter, data["data"]); prop=0.25)./200

# Write EEG Data file for visualization
showEEG && eegPlot(resample(EEG, 1//4), info["acquisition"]["sensors"]; ICoNargs...)

# EXTRACT TRIALS
################
# vectors of samples where the trials start for each class 1, 2,...
trig=[[i+offset for i in eachindex(stim) if stim[i]==j && i+offset+td<=ns] for j=1:z]

# Get Trials
𝐗=[EEG[trig[i][j]:trig[i][j]+td-1,:] for i=1:z for j=1:length(trig[i])]

# Get Labels
y=[i for i=1:z for j=1:length(trig[i])]

# ENCODING
##########
# arithmetic mean ERP of TargetTrial
Y=fVec(mean, 𝕄Vector([𝐗[i] for i in eachindex(𝐗) if y[i]==2]))
showEEG && eegPlot(resample(Y, 1//4), info["acquisition"]["sensors"]; ICoNargs...)

# Show topograpic maps of geometric means
if showTopoMaps
	Means=ℍVector([mean(Fisher, ℍVector([ℍ(𝐂[i][1:ne, 1:ne]) for i=1:length(𝐂) if y[i]==j])) for j=1:z])
	zlabels=Array{String, 1}(collect(keys(info["stim"]["labels"])))
	topoPlot(Means, info["acquisition"]["sensors"];
	    title="Class Means",
		mapLabels=zlabels,
	    monopolar=true,
	    scaleMode=:global)
end


# Minimum Distance to Mean
args=(shuffle=true, tol=1e-5)
cvMDM = cvAcc(MDM(Fisher), 𝐂, y; args...)
cvMDMscm = cvAcc(MDM(Fisher), 𝐂scm, y; args...)
cvMDMlw = cvAcc(MDM(Fisher), 𝐂lw, y; args...)
cvMDMnlw = cvAcc(MDM(Fisher), 𝐂nlw, y; args...)


# Lasso Logistic Regression (TS)
args=(shuffle=true, tol=1e-5, w=:b)
# Lasso Logistic Regression (TS) with reduced tangent vectors
cvLASSO = cvAcc(ENLR(Fisher), 𝐂, y; vecRange=1:ne, args...)
# Elastic-Net Logistic Regression (TS) with reduced tangent vectors
cvElNet05r = cvAcc(ENLR(Fisher), 𝐂, y; vecRange=1:ne, alpha=0.5, args...)

cvLASSOscm = cvAcc(ENLR(Fisher), 𝐂scm, y; vecRange=1:ne, args...)
cvLASSOlw = cvAcc(ENLR(Fisher), 𝐂lw, y; vecRange=1:ne, args...)
cvLASSOnl = cvAcc(ENLR(Fisher), 𝐂nl, y; vecRange=1:ne, args...)



# TArget geometric mean
# TAgm, iter, conv = gMean(ℍVector([𝐂[i] for i in eachindex(𝐂) if y[i]==2]); ⍰=true)

# Non Target geometric mean
# NTgm, iter, conv = gMean(ℍVector([𝐂[i] for i in eachindex(𝐂) if y[i]==1]); ⍰=true)

# using Plots
# heatmap(Matrix(TAgm); c=:bluesreds)


################ Example for working with a particular file, e.g., file 10

o=readNY(files[10]; bandpass=(1, 16), msg="read file: $(basename(file))")
# type `o` (letter o) and press the ENTER key in the REPL to see the structure

# multivariate regression ERP mean with data-driven weights
w=[[1/(norm(o.X[o.cstim[i][j]+o.offset:o.cstim[i][j]+o.offset+o.wl-1,:])^2) for j=1:length(o.cstim[i])] for i=1:o.nc]
Y=mean(o.X, o.wl, o.cstim; weights=w)[2]

# PCA to keep only 4 components
Y=Y*eigvecs(cov(SimpleCovariance(), Y))[:, o.ne-3:o.ne]

# encoding
covEstimator=LinearShrinkage(ConstantCorrelation())
𝐂lw=ℍVector([ℍ(cov(covEstimator, [X Y])) for X ∈ o.trials])


# MDM classification
args=(shuffle=true, tol=1e-5)
cvMDMlw = cvAcc(MDM(Fisher), 𝐂lw, o.y; args...)

# release memory
waste(o)

=#
