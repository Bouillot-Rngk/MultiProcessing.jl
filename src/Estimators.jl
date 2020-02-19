module Estimators
using LinearAlgebra, PosDefManifold, PosDefManifoldML, CovarianceEstimation,
      Dates, Distributions, PDMats, Revise, BenchmarkTools, Plots, GR
using EEGio, FileSystem, EEGpreprocessing, System, ERPs,
      DelphiApplications, Tyler, EEGtomography

# ? ¤ CONTENT ¤ ? #

# STRUCTURES
# EEG | holds data and metadata of an EEG recording

# FUNCTIONS:
# SCMP300 	 | read an EEG struct anc compute simple Sample Covariance Matrix
# TMEP300  	 | read an EEG struct anc compute Tyler's M-estimator Covariance Matrix
# nrTMEP300  | read an EEG struct anc compute normalize regularized Tyler's M-estimator Covariance Matrix

export
	SCMP300,
	TMEP300,
	nrTMEP300,
	nlseP300

#ajouter un regularisateur sur les estimateurs e^-6 et chercher un seuil

function SCMP300(o)
#-----------------------------------------------------------------------------------#
#calcul de la Sample Covariance Matrix simple de la base de donnees, avec moyenne
#ponderee
#-----------------------------------------------------------------------------------#
#Input :
#     o::EEG => Structure de EEGio.jl apres lecture de la bdd par readNY
#Output :
#     Clw : Matrice de covariance etendue (supertrials)

#Calcul des poids et moyenne + PCA pour ne conserver que 4 elements
      w=[[1/(norm(o.X[o.cstim[i][j]+o.offset:o.cstim[i][j]+o.offset+o.wl-1,:])^2) for j=1:length(o.cstim[i])] for i=1:o.nc]
	Y=mean(o.X, o.wl, o.cstim; weights=w)[2]
	Y=Y*eigvecs(cov(SimpleCovariance(), Y))[:, o.ne-3:o.ne]

#Calcul de la matrice de covariance ponderee
	Clw=ℍVector([ℍ(cov(SimpleCovariance(), [X Y])) for X ∈ o.trials])
	#regularization
	R=Hermitian(Matrix{eltype(Clw[1])}(I, size(Clw[1]))*0.0001)
	for C in Clw C+=R end
#Calculs annexes (normalization de det, etc)


	return Clw
end

function TMEP300(o)
#-----------------------------------------------------------------------------------#
#calcul de la Covariance Matrix par estimateur de Tyler de la base de donnees, avec
#moyenne ponderee
#-----------------------------------------------------------------------------------#
#Input :
#     o::EEG => Structure de EEGio.jl apres lecture de la bdd par readNY
#Output :
#     Clw : Matrice de covariance etendue, calcul TME (supertrials)
#Calcul des poids et moyenne + PCA pour ne conserver que 4 elements
	w=[[1/(norm(o.X[o.cstim[i][j]+o.offset:o.cstim[i][j]+o.offset+o.wl-1,:])^2) for j=1:length(o.cstim[i])] for i=1:o.nc]
	Y=mean(o.X, o.wl, o.cstim; weights=w)[2]
	Y=Y*eigvecs(cov(SimpleCovariance(), Y))[:, o.ne-3:o.ne]
#Calcul de la ;atrice de covariance par estimateur de Tyler
	Clw=ℍVector([ℍ(tme([X Y]')) for X ∈ o.trials])
	#regularization
	R=Hermitian(Matrix{eltype(Clw[1])}(I, size(Clw[1]))*0.0001)
	for C in Clw C+=R end

	return Clw
end

function nrTMEP300(o)
#-----------------------------------------------------------------------------------#
#calcul de la Covariance Matrix par estimateur de Tyler de la base de donnees, avec
#moyenne ponderee
#-----------------------------------------------------------------------------------#
#Input :
#     o::EEG => Structure de EEGio.jl apres lecture de la bdd par readNY
#Output :
#     Clw : Matrice de covariance etendue, calcul nrTME (supertrials)
#Calcul des poids et moyenne + PCA pour ne conserver que 4 elements
	w=[[1/(norm(o.X[o.cstim[i][j]+o.offset:o.cstim[i][j]+o.offset+o.wl-1,:])^2) for j=1:length(o.cstim[i])] for i=1:o.nc]
	Y=mean(o.X, o.wl, o.cstim; weights=w)[2]
	Y=Y*eigvecs(cov(SimpleCovariance(), Y))[:, o.ne-3:o.ne]
#Calcul de la matrice de covariance par estimateur de Tyler non regularized
	Clw=ℍVector([ℍ(nrtme([X Y]'; reg=:LW)) for X ∈ o.trials])
	#regularization
	R=Hermitian(Matrix{eltype(Clw[1])}(I, size(Clw[1]))*0.0001)
	for C in Clw C+=R end

	return Clw
end

function nlseP300(o)
#-----------------------------------------------------------------------------------#
#calcul de la Covariance Matrix par estimateur de Tyler de la base de donnees, avec
#moyenne ponderee
#-----------------------------------------------------------------------------------#
#Input :
#     o::EEG => Structure de EEGio.jl apres lecture de la bdd par readNY
#Output :
#     Clw : Matrice de covariance etendue, calcul nrTME (supertrials)
#Calcul des poids et moyenne + PCA pour ne conserver que 4 elements
	w=[[1/(norm(o.X[o.cstim[i][j]+o.offset:o.cstim[i][j]+o.offset+o.wl-1,:])^2) for j=1:length(o.cstim[i])] for i=1:o.nc]
	Y=mean(o.X, o.wl, o.cstim; weights=w)[2]
	Y=Y*eigvecs(cov(SimpleCovariance(), Y))[:, o.ne-3:o.ne]
#Calcul de la matrice de covariance par estimateur non linear shrinkage
	Clw=ℍVector([ℍ(cov(nlse(), [X Y])) for X ∈ o.trials])
#Regularization
	R=Hermitian(Matrix{eltype(Clw[1])}(I, size(Clw[1]))*0.001)
	for C in Clw C+=R end

	return Clw
end


end #module
