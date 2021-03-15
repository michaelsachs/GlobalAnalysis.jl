

#define composite type to hold data
struct DataStruct
    name::String
    x::Vector{Float64}
    y::Vector{Float64}
    z::Array{Union{Float64,Missing},2}
    var::Array{String}
    val::Array{Float64}
    unit::Array{String}
    file::Array{String}
end

#=
#julia 1.6.0rc1 crashes when unsing interpolations if x or y are allowed missing
struct DataStruct{F<:Float64, M<:Missing, S<:String}
    name::S
    x::Vector{Union{F,M}}
    y::Vector{Union{F,M}}
    z::Array{Union{F,M},2}
    var::Array{S}
    val::Array{F}
    unit::Array{S}
    file::Array{S}
end
=#