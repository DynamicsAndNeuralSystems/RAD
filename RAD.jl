using StatsBase
"""
    RAD(x, τ=1, doAbs=true)
Compute the rescaled auto-density, a metric for inferring the
distance to criticality that is insensitive to uncertainty in the noise strength.
Calibrated to experiments on the radial compenent of the Hopf bifurcation with 
variable and unknown measurement noise.

Inputs:
    x:      The input time series (vector).
    doAbs:  Whether to centre the time series at 0 then take absolute values (logical flag)
    τ:      The embedding and differencing delay in units of the timestep (integer)

Outputs:
    f:      The RAD feature value
"""
function RAD(z, τ::Integer = 1, doAbs::Bool = true)
    if doAbs
        z = z .- median(z)
        z = abs.(z)
    end

    y = @view z[(τ + 1):end]
    x = @view z[1:(end - τ)]

    # Median split
    subMedians = x .< median(x)
    superMedianSD = std(x[.!subMedians])
    subMedianSD = std(x[subMedians])

    # Properties of the auto-density
    sigma_dx = std(y - x)
    densityDifference = 1 / superMedianSD - 1 / subMedianSD

    f = sigma_dx * densityDifference
end
