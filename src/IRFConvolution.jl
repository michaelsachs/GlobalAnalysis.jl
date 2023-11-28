using Distributions
using Interpolations 
using DSP
using Plots

#=
interpolation with gaussian IRF uses window around zero with width adapting depending on gaussian width,
then reinterpolates onto original time vector and inserts data outside interpolation window to stitch
trace back together

To address for variable time spacing:
- undefined if less negative time exists than assumed: might not be a problem for conv, just need to cut
    all negative data with no timepoints after conv
- signal currently needs to onset around 0
- what if both signal and IRF are shifted (TCSPC)? or if just signal is shifted (slow TAS)?
    need to shift gaussian to match signal onset? option of per WL IRF?
- take into account μ shift during ODE such that final conv kin is sampled like expKin (time-μ)

- shift interpolation window gives same result as long as centred around gaussian as is indep of time vector
=#

#= some test parameters
μ = 0.3
σ = 0.1
timePos = time[time .> 0]
simKin = exp.(-0.1 .* timePos)
=#

function defineIRFtime(time, μ, σ)
    timeDiff = diff(time)
    #time = time .+ μ
    # check if time points are evenly spaced; all(x,itr) can short circuit, 
    #   thus performance hit should be small
    if all(x->x==timeDiff[1], timeDiff) == false
        # time step in IRF window is smallest time step in dataset
        IRFStep = minimum(timeDiff)
        # IRF window width in pos/neg direction
        IRFWindow = 30*σ
        # evenly spaced time vector for IRF convolution
        timeIRFPos = collect(IRFStep/2:IRFStep:IRFWindow)
        # convolved trace will be cut here as it will have a sharp drop where kinetic data ends
        IRFWindowCut = IRFWindow*0.5
        
        # combine IRF and post-IRF times in array for single call to ODE solver
        # shift post-IRF time by -μ to evaluate correct points for original time vector in ODE
        timeODE = sort(unique([timeIRFPos; time[time .> IRFWindowCut] .- μ]))

        # keep track of which time in ODE time vector is IRF/post-IRF (some may be in both)
        timeODEBool = [in.(timeODE, Ref(Set(timeIRFPos))), in.(timeODE, Ref(Set(time .- μ)))]
        
        return timeODE, timeIRFPos, IRFWindowCut, timeODEBool
    end
end


# also works for more than one kinetic in vector
# time includeds negative data, kin is positive data out of simulation only
# TODO: small offset btw original and convkin in convKin. check 
function convolveIRF(time, kin, μ, σ, timeIRFPos, IRFWindowCut, timeODEBool)
    # check if time points are evenly spaced; all(x,itr) can short circuit, 
    #   thus performance hit should be small
    timeDiff = diff(time)
    if all(x->x==timeDiff[1], timeDiff) == false
        # zero pre-pad simulated data
        kinIRF = [zeros(length(timeIRFPos),size(kin,2)); kin[timeODEBool[1],:]]

        # mirror array around zero, pos/neg symmetric
        timeIRF = [reverse(-timeIRFPos); timeIRFPos]

        # gaussian IRF
        IRF = pdf.(Normal(μ, σ), timeIRF)
        # discrete convolution, hence multiply by time step
        kinIRFConv = conv(kinIRF, IRF) .* (timeIRFPos[2]-timeIRFPos[1])

        # make conv data same length as IRF data to restore time correspondence
        idxFirst = Int64(1 + floor(length(timeIRF) / 2))
        idxLast = idxFirst + size(kinIRF,1) - 1       
        #kinConvSame = @view kinIRFConv[idxFirst:idxLast]
        kinConvSame = kinIRFConv[idxFirst:idxLast,:]

        # for interpolation back onto original time vector
        # no interpolation along second dimension 
        if size(kin,2) == 1
            kinConvSame = [kinConvSame...]
            itpKinConv = interpolate((timeIRF,), kinConvSame, Gridded(Linear()))
        else 
            itpKinConv = interpolate((timeIRF,collect(1:size(kin,2)),), kinConvSame, 
                (Gridded(Linear()),NoInterp()))
        end

    

        

        # final data contains: zeros pre-IRF, interpolated cut IRF data, unaltered post-IRF data
        #timePos = @view time[time .> 0]
        #kinConv = [zeros(length(@view time[time .< timeIRF[1]])); 
        ##    itpKinConv(@view time[timeIRF[1] .≤ time .≤ IRFWindowCut]); 
        #    @view kin[timeODEBool[2]]]
        #kinConv = deepcopy(kinConv)

        timePos = time[time .> 0]
        if size(kin,2) == 1
            kinConv = [zeros(length(time[time .< timeIRF[1]])); 
                itpKinConv(time[timeIRF[1] .≤ time .≤ IRFWindowCut]); 
                kin[timeODEBool[2]]]
        else
            kinConv = [zeros(length(time[time .< timeIRF[1]]), size(kin,2)); 
                itpKinConv(time[timeIRF[1] .≤ time .≤ IRFWindowCut], collect(1:size(kin,2))); 
                kin[timeODEBool[2],:]]
        end
        kinConv = deepcopy(kinConv)


        #plot(timeIRF,kinIRF, linewidth=2)
        #plot!(time,kinConv, xlim=(-5,20))
        
        return kinConv

    end

end


# time includeds negative data, kin is positive data out of simulation only
function convolveIRF(time, kin, μ, σ)
    # check if time points are evenly spaced; all(x,itr) can short circuit, 
    #   thus performance hit should be small
    timeDiff = diff(time)
    # equal spacing: no IRF window needed
    if all(x->x==timeDiff[1], timeDiff) == true
        # zero pre-pad simulated data
        timePos = @view time[time .> 0]
        kinIRF = [zeros(length(time)-length(timePos)); kin]

        # gaussian IRF
        IRF = pdf.(Normal(μ, σ), time)
        # discrete convolution, hence multiply by time step
        kinIRFConv = conv(kinIRF, IRF) * (time[2]-time[1])

        # make conv data same length as IRF data to restore time correspondence
        idxFirst = Int64(1 + floor(length(timePos) / 2))
        idxLast = idxFirst + length(kinIRF) - 1       
        #kinConvSame = @view kinIRFConv[idxFirst:idxLast]
        kinConvSame = kinIRFConv[idxFirst:idxLast]

        #plot(time,kinConvSame, xlim=(-5,20), linewidth=2, label= "μ = "*string(μ))
       
        return kinConv

    end

end

"""
Performs numerical integration of `y` over `x` via trapezoidal method.
"""
function trapezIntegration(x, y) 
    # Check same length of X and Y
    @assert length(x) == length(y)
    out = 0.0
    for n in 2:length(x)
      out += 0.5*(x[n] - x[n-1])*(y[n] + y[n-1])
    end
    return out
end

# even spacing, measured IRF
function defineODEtime(time::Array{<:Real,1})
    timeDiff = diff(time)
    # check if time points are evenly spaced; all(x,itr) can short circuit, 
    # thus performance hit should be small
    if all(x->x==timeDiff[1], timeDiff) == true

        # time vector with same number of time steps as IRF but no negative zeros
        timeODE = collect(0:time[2]-time[1]:time[end]-time[1])
        
        return timeODE
    end
end



# even spacing, measured IRF
# time includes negative data, kin is positive data out of simulation only
# kin must have same number of time steps as IRF but no negative zeros to ensure 
# that convolved trace has enough time points before dropoff (this implementation 
# shoudl work as long as IRF has at least one zero point before rise)
#function convolveIRFReg(time::Array{<:Real,1}, kin::Array{<:Real,1}, IRF::Array{<:Real,1}) 
function convolveIRFReg(time, kin, IRF) 
    # check if time points are evenly spaced; all(x,itr) can short circuit, 
    #   thus performance hit should be small
    timeDiff = diff(time)
    # equal spacing required
    #if all(x->x==timeDiff[1], timeDiff) == true

        #normalise IRF by its area
        IRF = IRF ./ trapezIntegration(time,IRF)

        # discrete convolution, hence multiply by time step
        kinIRFConv = DSP.conv(kin, IRF) .* (time[2]-time[1])

        # make conv data same length as IRF data to restore time correspondence
        if size(kin,2) == 1
            kinConvSame = kinIRFConv[1:length(IRF)]
        else
            kinConvSame = kinIRFConv[1:length(IRF),:]
        end

        return kinConvSame

    #end

end




