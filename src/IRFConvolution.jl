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
- take into account μ shift during ODE such that final conv kin is sampled like expKin (Time-μ)

- shift interpolation window gives same result as long as centred around gaussian as is indep of time vector
=#

#= some test parameters
μ = 0.3
σ = 0.1
TimePos = Time[Time .> 0]
simKin = exp.(-0.1 .* TimePos)
=#

function defineIRFTime(Time, μ, σ)
    TimeDiff = diff(Time)
    #Time = Time .+ μ
    # check if time points are evenly spaced; all(x,itr) can short circuit, 
    #   thus performance hit should be small
    if all(x->x==TimeDiff[1], TimeDiff) == false
        # time step in IRF window is smallest time step in dataset
        IRFStep = minimum(TimeDiff)
        # IRF window width in pos/neg direction
        IRFWindow = 30*σ
        # evenly spaced time vector for IRF convolution
        TimeIRFPos = collect(IRFStep/2:IRFStep:IRFWindow)
        # convolved trace will be cut here as it will have a sharp drop where kinetic data ends
        IRFWindowCut = IRFWindow*0.5
        
        # combine IRF and post-IRF times in array for single call to ODE solver
        # shift post-IRF time by -μ to evaluate correct points for original time vector in ODE
        TimeODE = sort(unique([TimeIRFPos; Time[Time .> IRFWindowCut] .- μ]))

        # keep track of which time in ODE time vector is IRF/post-IRF (some may be in both)
        TimeODEBool = [in.(TimeODE, Ref(Set(TimeIRFPos))), in.(TimeODE, Ref(Set(Time .- μ)))]
        
        return TimeODE, TimeIRFPos, IRFWindowCut, TimeODEBool
    end
end


# also works for more than one kinetic in vector
# Time includeds negative data, kin is positive data out of simulation only
# TODO: small offset btw original and convkin in convKin. check 
function convolveIRF(Time, Kin, μ, σ, TimeIRFPos, IRFWindowCut, TimeODEBool)
    # check if time points are evenly spaced; all(x,itr) can short circuit, 
    #   thus performance hit should be small
    TimeDiff = diff(Time)
    if all(x->x==TimeDiff[1], TimeDiff) == false
        # zero pre-pad simulated data
        KinIRF = [zeros(length(TimeIRFPos),size(Kin,2)); Kin[TimeODEBool[1],:]]

        # mirror array around zero, pos/neg symmetric
        TimeIRF = [reverse(-TimeIRFPos); TimeIRFPos]

        # gaussian IRF
        IRF = pdf.(Normal(μ, σ), TimeIRF)
        # discrete convolution, hence multiply by time step
        KinIRFConv = conv(KinIRF, IRF) .* (TimeIRFPos[2]-TimeIRFPos[1])

        # make conv data same length as IRF data to restore time correspondence
        idxFirst = Int64(1 + floor(length(TimeIRF) / 2))
        idxLast = idxFirst + size(KinIRF,1) - 1       
        #KinConvSame = @view KinIRFConv[idxFirst:idxLast]
        KinConvSame = KinIRFConv[idxFirst:idxLast,:]

        # for interpolation back onto original time vector
        # no interpolation along second dimension 
        if size(Kin,2) == 1
            KinConvSame = [KinConvSame...]
            itpKinConv = interpolate((TimeIRF,), KinConvSame, Gridded(Linear()))
        else 
            itpKinConv = interpolate((TimeIRF,collect(1:size(Kin,2)),), KinConvSame, 
                (Gridded(Linear()),NoInterp()))
        end

    

        

        # final data contains: zeros pre-IRF, interpolated cut IRF data, unaltered post-IRF data
        #TimePos = @view Time[Time .> 0]
        #KinConv = [zeros(length(@view Time[Time .< TimeIRF[1]])); 
        ##    itpKinConv(@view Time[TimeIRF[1] .≤ Time .≤ IRFWindowCut]); 
        #    @view Kin[TimeODEBool[2]]]
        #KinConv = deepcopy(KinConv)

        TimePos = Time[Time .> 0]
        if size(Kin,2) == 1
            KinConv = [zeros(length(Time[Time .< TimeIRF[1]])); 
                itpKinConv(Time[TimeIRF[1] .≤ Time .≤ IRFWindowCut]); 
                Kin[TimeODEBool[2]]]
        else
            KinConv = [zeros(length(Time[Time .< TimeIRF[1]]), size(Kin,2)); 
                itpKinConv(Time[TimeIRF[1] .≤ Time .≤ IRFWindowCut], collect(1:size(Kin,2))); 
                Kin[TimeODEBool[2],:]]
        end
        KinConv = deepcopy(KinConv)


        #plot(TimeIRF,KinIRF, linewidth=2)
        #plot!(Time,KinConv, xlim=(-5,20))
        
        return KinConv

    end

end


# Time includeds negative data, kin is positive data out of simulation only
function convolveIRF(Time, Kin, μ, σ)
    # check if time points are evenly spaced; all(x,itr) can short circuit, 
    #   thus performance hit should be small
    TimeDiff = diff(Time)
    # equal spacing: no IRF window needed
    if all(x->x==TimeDiff[1], TimeDiff) == true
        # zero pre-pad simulated data
        TimePos = @view Time[Time .> 0]
        KinIRF = [zeros(length(Time)-length(TimePos)); Kin]

        # gaussian IRF
        IRF = pdf.(Normal(μ, σ), Time)
        # discrete convolution, hence multiply by time step
        KinIRFConv = conv(KinIRF, IRF) * (Time[2]-Time[1])

        # make conv data same length as IRF data to restore time correspondence
        idxFirst = Int64(1 + floor(length(TimePos) / 2))
        idxLast = idxFirst + length(KinIRF) - 1       
        #KinConvSame = @view KinIRFConv[idxFirst:idxLast]
        KinConvSame = KinIRFConv[idxFirst:idxLast]

        #plot(Time,KinConvSame, xlim=(-5,20), linewidth=2, label= "μ = "*string(μ))
       
        return KinConv

    end

end

#numerical integration via trapezoidal method
function trapezIntegration(X, Y) 
    # Check same length of X and Y
    @assert length(X) == length(Y)
    out = 0.0
    for n in 2:length(X)
      out = out + 0.5*(X[n] - X[n-1])*(Y[n] + Y[n-1])
    end
    return out
end

# even spacing, measured IRF
function defineODETime(Time::Array{<:Real,1})
    TimeDiff = diff(Time)
    # check if time points are evenly spaced; all(x,itr) can short circuit, 
    # thus performance hit should be small
    if all(x->x==TimeDiff[1], TimeDiff) == true

        # time vector with same number of time steps as IRF but no negative zeros
        TimeODE = collect(0:Time[2]-Time[1]:Time[end]-Time[1])
        
        return TimeODE
    end
end



# even spacing, measured IRF
# Time includes negative data, kin is positive data out of simulation only
# Kin must have same number of time steps as IRF but no negative zeros to ensure 
# that convolved trace has enough time points before dropoff (this implementation 
# shoudl work as long as IRF has at least one zero point before rise)
#function convolveIRFReg(Time::Array{<:Real,1}, Kin::Array{<:Real,1}, IRF::Array{<:Real,1}) 
function convolveIRFReg(Time, Kin, IRF) 
    # check if time points are evenly spaced; all(x,itr) can short circuit, 
    #   thus performance hit should be small
    TimeDiff = diff(Time)
    # equal spacing required
    #if all(x->x==TimeDiff[1], TimeDiff) == true

        #normalise IRF by its area
        IRF = IRF ./ trapezIntegration(Time,IRF)

        # discrete convolution, hence multiply by time step
        KinIRFConv = DSP.conv(Kin, IRF) .* (Time[2]-Time[1])

        # make conv data same length as IRF data to restore time correspondence
        if size(Kin,2) == 1
            KinConvSame = KinIRFConv[1:length(IRF)]
        else
            KinConvSame = KinIRFConv[1:length(IRF),:]
        end

        return KinConvSame

    #end

end




