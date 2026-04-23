
# define composite type to hold data
struct DataStruct{T}
    name::String
    x::Vector{Float64}
    y::Vector{Float64}
    z::Array{T,2}
    var::Array{String}
    val::Array{Float64}
    unit::Array{String}
    file::Array{String}
end

"""
Precomputed dataset metadata reused by `paramToSSR` during kinetic optimisation.
"""
struct SSRMetaData{T}
    dataNorm::T
    hasNaN::Bool
end