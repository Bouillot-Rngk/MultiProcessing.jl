#   Module "Tyler.jl" for Julia language
#
#   MIT License
#   Copyright (c) 2020
#   Marco Congedo, CNRS, UGA, Grenoble-INP, Grenoble, France
#   https://sites.google.com/site/marcocongedo/
#
# This module implements Tyler's M Estimator of covariance
# matrix 'shape' (Tyler, 1985) and the normalized regularized version
# of Zhang and Wiesel (2016)
#
#
# ? CONTENT
#
# FUNCTIONS:
# tme   | Tyler M-Estimator fixed point algorithm
# tme2  | Version with less memory consuption
# nrtme | normalized regularized Tyler's M-Estimator
#
# REFERENCES
# David E. Tyler (1987)
# A Distribution-Free M-Estimator of Multivariate Scatter
# The Annals of Statistics, 15(1), 234-251.
# https://projecteuclid.org/download/pdf_1/euclid.aos/1176350263

# Teng Zhang, Ami Wiesel (2016)
# Automatic diagonal loading for Tyler's robust covariance estimator
# IEEE Statistical Signal Processing Workshop (SSP), 1-5.
# https://sciences.ucf.edu/math/tengz/wp-content/uploads/sites/45/2016/08/automatic-diagonal-loading-3.pdf

module Tyler

using LinearAlgebra, Statistics, PosDefManifold

export
	tme,
	tme2,
	nrtme

function tme(X::AbstractArray{T};
			 tol 	:: Real = real(T)(1e-6),
			 maxiter :: Int = 200,
			 verbose :: Bool = false) where T<: Union{Real, Complex}
    n, t = size(X)
    R = Matrix{T}(I, n, n)
    💡 = Matrix{T}(undef, n, n)
    iter, 😋 = 1, false

	verbose && println("Iterating M-estimator fixed-point algorithm...")
    while true
        C=cholesky(R)
        fill!(💡, zero(T))
        @inbounds for i=1:t
			@views v = C.L\X[:, i]
			💡 += (X[:, i].*X[:, i]')./(v⋅v)
        end
        💡 *= inv(tr(💡))
        conv = norm(💡-R)/norm(R)
        R[:]=💡
        verbose && println("iteration: ", iter, "; convergence: ", conv)
        (overRun = iter == maxiter) && @warn("M-estimator reached the max number of iterations before convergence:", iter)
        (😋 = conv <= tol) || overRun==true ? break : iter += 1
    end # while
	verbose && @info("Convergence has "*(😋 ? "" : "not ")*"been attained.\n\n")
    return 💡
end


# version with less memory usage and slighter faster
function tme2(X::AbstractArray{T};
  			  tol 	:: Real = real(T)(1e-6),
			  maxiter :: Int = 100,
			  verbose :: Bool = false) where T<: Union{Real, Complex}
    n, t = size(X)
    R = Matrix{T}(I, n, n)
	iter, 😋, oldconv = 1, false, T(0.)

	verbose && println("Iterating M-estimator fixed-point algorithm...")
    while true
		C=cholesky(R)
        fill!(R, zero(T))
		temp = T(0.)
        @inbounds for i=1:t
			@views v = C.L\X[:, i]
			a = v⋅v
			temp += a
			R += (X[:, i].*X[:, i]')./a
        end
        R *= inv(tr(R))
		temp=sqrt(temp)/t; conv = abs(temp - oldconv); oldconv = temp
        verbose && println("iteration: ", iter, "; convergence: ", conv)
        (overRun = iter == maxiter) && @warn("M-estimator reached the max number of iterations before convergence:", iter)
        (😋 = conv <= tol) || overRun==true ? break : iter += 1
    end # while
	verbose && @info("Convergence has "*(😋 ? "" : "not ")*"been attained.\n\n")
    return R
end


function nrtme(X::AbstractArray{T};
  			   reg 	:: Symbol = :RMT,
			   tol 	:: Real = real(T)(1e-6),
			   maxiter :: Int = 200,
			   verbose :: Bool = false) where T<: Union{Real, Complex}
    n, t = size(X)
    R = Matrix{T}(I, n, n)
    💡 = zeros(T, n, n)
	x² = Vector{T}(undef, t)
	x = Matrix{T}(undef, n, 1)
	v = Vector{T}(undef, n)
    iter, 😋, α, β, nt⁻¹ = 1, false, 0., 0., n/t

	@inbounds @views for i=1:t x²[i] = X[:, i]⋅X[:, i] end

	if reg == :RMT
		@inbounds for i=1:t
			x[:] = X[:, i]									# |
			BLAS.gemm!('N', 'T', inv(x²[i]), x, x, 1., 💡) 	# | instead of 💡 += (X[:, i].*X[:, i]')./x²[i]
		end
		ζ = n*tr((💡./t)^2)-nt⁻¹-1.
	else
		scm = (X'*X).*inv(n)
		ζ = (n*tr(scm^2)/(tr(scm))^2)-1.
	end
	α = clamp(inv(t) * ((ζ + 1 + n) / (ζ + nt⁻¹)), 0., 1.)
	β = 1. - α
	αn⁻¹ = α/n

	verbose && println("Iterating nrM-estimator fixed-point algorithm...")
    while true
		if n ≤ 40
        	L, U⁻¹ = choInv(R)
			trR⁻¹ = sumOfSqr(LowerTriangular(U⁻¹')) # tr(R⁻¹) = tr(Diagonal(U⁻¹*U⁻¹'))
		else
			L = cholesky(R)
			trR⁻¹ = sumOfSqr(inv(L)) # tr(R⁻¹) = tr(Diagonal(L⁻¹'*L⁻¹))
		end
        fill!(💡, zero(T))
		g(x, β) = BLAS.gemm('N', 'T', β, x, x)
		for i = 1:t
			x[:] = X[:, i]
			v[:] = L\x
			#💡 += (β*(x.*x')+(αn⁻¹*x²[i])*I) ./ (β*(v⋅v)+αn⁻¹*trR⁻¹*x²[i])
			💡 += (g(x, β)+(αn⁻¹*x²[i])*I) ./ (β*(v⋅v)+αn⁻¹*trR⁻¹*x²[i])
        end
        💡 *= (inv(tr(💡)))
        conv = norm(💡-R)/norm(R)
        R[:] = 💡

        verbose && println("iteration: ", iter, "; convergence: ", conv)
        (overRun = iter == maxiter) && @warn("nrM-estimator reached the max number of iterations before convergence:", iter)
        (😋 = conv <= tol) || overRun==true ? break : iter += 1
    end # while
	verbose && @info("Convergence has "*(😋 ? "" : "not ")*"been attained.\n\n")
    return 💡
end

end # module


# test the M-estimators
#
# create data drawn randomly from a multivariate t-student
# distribution with `df` degrees of freedom
# and check how far the estimated shape is different from
# `trueC`

#=
using BenchmarkTools, Distributions

n, 512, df = 30, 3.0
trueC=randP(n)
trueC=trueC/tr(trueC)
tdist=MvTDist(3., zeros(n), PDMat(Matrix(trueC)))
X=rand(tdist, t)

tme(X, verbose=true)
tme2(X, verbose=true)
nrtme(X, verbose=true)
nrtme(X, reg=:LW, verbose=true)

# run 100 simulations
# and check the Fisher distance between true and estimated shape
t=1_000
for i=1:100
	tdist=MvTDist(df, zeros(n), PDMat(Matrix(trueC)))
	X=rand(tdist, t);
	#println(cov(td))
	#scm=cov(X')
	#scm=scm/tr(scm)
	M = (tme(X))
	d1[i]=distance(Fisher, Hermitian(trueC), Hermitian(scm))
	d2[i]=distance(Fisher, Hermitian(trueC), Hermitian(M))
end
10*log10(mean(d1))
10*log10(mean(d2))

# benchmark algorithms
@benchmark(tme(X))
@benchmark(tme2(X))
@benchmark(nrtme(X))

=#
