# set BLAS threads to 1 to avoid conflicts with Julia threading
using LinearAlgebra
LinearAlgebra.BLAS.set_num_threads(1)

using Catalyst
using GlobalAnalysis
using Metaheuristics
using Plots