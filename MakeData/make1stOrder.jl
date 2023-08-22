using Catalyst, Plots, DifferentialEquations
using Statistics 
include("../src/IRFConvolution.jl")

# define kinetic model
rs = @reaction_network begin
    k1, A --> B
    k2, B --> C
    k3, C --> 0
end
  
# starting conditions
u0 = [:A => 1, :B => 0, :C => 0]
# parameters
p = (:k1 => 1.5, :k2 => 1, :k3 => 0.5)
# time
time = 0:0.1:10
tpos = time[time.>0]
tspan = (0, 10)

# solve problem
prob = ODEProblem(rs, u0, tspan, p; saveat=tpos)
sol  = solve(prob, Tsit5())
sol_array = Array(sol)

# plot kinetic traces 
plot(sol, legend=:none, yaxis="Population", xaxis="Time", linewidth = 2, left_margin=5Plots.mm)

# IRF parameters
μ=0.3
σ=0.08

# define IRF as a gaussian distribution 
IRF = pdf.(Normal(μ, σ), time)
# convolve kinetic traces
ConvKin1 = GlobalAnalysis.convolveIRFReg(time, sol_array[1,:], IRF)
ConvKin2 = GlobalAnalysis.convolveIRFReg(time, sol_array[2,:], IRF)
ConvKin3 = GlobalAnalysis.convolveIRFReg(time, sol_array[3,:], IRF)

plot(time, ConvKin1, title="Kinetics", xlabel="Time", ylabel="Population", linewidth = 2, xguidefontsize=10, yguidefontsize=10, label="A(t)")
plot!(time, ConvKin2, label="B(t)", linewidth = 2)
plot!(time, ConvKin3, label="C(t)", linewidth = 2)


# Define a range of wavelengths for the spectral signatures 
len = length(ConvKin1)
step = (800 - 400) / (len - 1) 
wave_range = 400:step:800

# define three functions for the three signatures
gaussian_func(wave) = 100 * exp(-((wave - 600)^2) / (2 * 20^2))
gaussian_vals = gaussian_func.(wave_range)

gaussian_func1(wave) = 70 * exp(-((wave - 500)^2) / (2 * 20^2))
gaussian_vals1 = gaussian_func1.(wave_range) 

gaussian_func2(wave) = 50 * exp(-((wave - 700)^2) / (2 * 20^2))
gaussian_vals2 = gaussian_func2.(wave_range)

# Plot spectral signatures
plot(wave_range, gaussian_vals, title="Spectral Signatures",
    xlabel="Wavelength", legend=:none, linewidth=2, left_margin=5Plots.mm, label="A(t)", ylabel="Δ Absorbance")
plot!(wave_range, gaussian_vals1, linewidth=2,  label="B(t)")
plot!(wave_range, gaussian_vals2, linewidth=2,  label="C(t)")


# Assign each spectral signature to a kinetic trace by multiplying them together
trace_a = gaussian_vals .* transpose(ConvKin1)
trace_b = gaussian_vals1 .* transpose(ConvKin2)
trace_c = gaussian_vals2 .* transpose(ConvKin3)

# Create the final fake data by summing the individual traces
fake_data = trace_a .+ trace_b .+ trace_c
surface(fake_data, c=:thermal, gridalpha=0.3)

# Define the noise level (adjust as needed)
noise = randn(size(fake_data))
noise_level = 0.2
# Add the noise to the fake_data
fake_data_noise = fake_data + noise * noise_level
# Plot the fake_data with noise
heatmap(time, wave_range, fake_data_noise, title="Fake Data with Noise", xlabel="Time", ylabel="Wavelength (nm)",
colorbar_title="\n\nΔ Absorbance", right_margin=15Plots.mm, left_margin=10Plots.mm,
xguidefontsize=10, yguidefontsize=10, c=:thermal)


# convert to CSV, first row should be the time vector and first column should
# be the spectral range
t_ax = collect(t)
freq_ax = collect(freq_range)
dataInclAxes = [vcat(NaN,freq_ax) vcat(transpose(t_ax),fake_data_noise)]

# Convert the data matrix to a DataFrame
df = DataFrame(dataInclAxes[2:end, 1:end], :auto)
column_names = string.(dataInclAxes[1, 1:end])
rename!(df, Symbol.(column_names))
# CSV.write("/path/to/desired/location/1stOrderData.csv", df)
