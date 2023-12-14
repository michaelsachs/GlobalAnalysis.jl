
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
