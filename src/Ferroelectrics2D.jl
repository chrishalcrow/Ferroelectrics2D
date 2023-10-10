module Ferroelectrics2D

using Roots
using LinearAlgebra
using Makie, CairoMakie, LaTeXStrings
using FileIO
using StaticArrays
using FFTW

mutable struct Parameter
	PV::Vector{Float64}
	V0::Float64
	A2::Matrix{Float64}
	A4::Array{Float64, 4}
	G4::Array{Float64,4}
	A2t::Matrix{Float64}
	A4t::Array{Float64, 4}
	A2s::Matrix{Float64}
	A4s::Array{Float64, 4}
	rescaling::Float64
	length_scale::Float64
    is_electrostatic::Bool
end


mutable struct PolarisationField2D
    field::Array{Float64, 3}
    lp::Vector{Int64}
    ls::Vector{Float64}
    x::Vector{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}
    parameters::Parameter
end

export PolarisationField2D

include("Parameters.jl")

include("Properties.jl")
export energy

include("Derivatives.jl")


include("Diff.jl")
export gradient_flow!

PolarisationField2D(lp::Vector{Int64},ls::Vector{Float64},R2,material; is_electrostatic=false) = PolarisationField2D( zeros(lp[1],lp[2],3), lp, ls, [-ls[a]*(lp[a] - 1)/2.0 : ls[a] : ls[a]*(lp[a] - 1)/2.0 for a in 1:2 ], set_parameters_lithium!(R2,material, is_electrostatic=is_electrostatic) )

#=

include("Parameters.jl")
export set_parameters_lithium!


PolarisationField2D(lp,ls,R2,material,b) = PolarisationField2D( zeros(lp[1],lp[2],3), lp, ls, -ls*(lp - 1)/2.0 : ls : ls*(lp - 1)/2.0, set_parameters_lithium!(R2,material) , false, b)


function setgrid(N,dx)
	
	println("grid has ", N, " points and goes from ", -0.5*dx*(N-1), " to ", 0.5*dx*(N-1))	
	return -0.5*dx*(N-1):dx:0.5*dx*(N-1)
	
end

function xyzfsrt(P,R2)
	return R2'*P
end

function srtfxyz(P,R2)
	return R2*P
end

function inbfsrt(P,R2)
	return R2*P
end

function srtfinb(P,R2)
	return R2*P
end

=#
	
end

