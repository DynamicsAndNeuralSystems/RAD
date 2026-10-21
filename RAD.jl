using StatsBase
"""
    RAD(x, τ=1, doCentre=true)
Compute the rescaled auto-density, a metric for inferring the
distance to criticality that is insensitive to uncertainty in the noise strength.
Calibrated to experiments on the radial component of the Hopf bifurcation with
variable and unknown dynamical noise.

Inputs:
    x:        The input time series (vector).
    τ:        The embedding and differencing delay in units of the timestep (positive integer)
    doCentre: Whether to centre the time series at its median then take absolute values;
              use this unless the time series represents radial values (logical flag)

Outputs:
    f:      The RAD feature value, or NaN (with a warning) when either half of the
            median split has fewer than two values or no spread, as happens for short
            or heavily tied (e.g., integer-valued) time series
"""
function RAD(z::AbstractVector, τ::Integer = 1, doCentre::Bool = true)
    τ ≥ 1 || throw(ArgumentError("τ must be a positive integer"))
    if doCentre
        z = z .- median(z)
        z = abs.(z)
    end

    y = @view z[(τ + 1):end]
    x = @view z[1:(end - τ)]

    # Median split
    subMedians = x .< median(x)
    if count(subMedians) < 2 || count(.!subMedians) < 2
        @warn "RAD: a half of the median split has fewer than two values; returning NaN"
        return NaN
    end
    superMedianSD = std(x[.!subMedians])
    subMedianSD = std(x[subMedians])
    if superMedianSD == 0 || subMedianSD == 0
        @warn "RAD: a half of the median split has no spread (tied values); returning NaN"
        return NaN
    end

    # Properties of the auto-density
    sigma_dx = std(y - x)
    densityDifference = 1 / superMedianSD - 1 / subMedianSD

    f = sigma_dx * densityDifference
end
