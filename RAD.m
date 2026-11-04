function f = RAD(x, tau, doCentre)
    % RAD Compute the rescaled auto-density, a metric for inferring the
    % distance to criticality that is insensitive to uncertainty in the noise strength
    % INPUTS
    %   x:        A vector, the input time series
    %   tau:      A positive integer; the embedding and differencing delay in
    %             units of the timestep (default: 1)
    %   doCentre: Whether to center the time series at its median then take
    %             absolute values; use this unless the time series represents
    %             radial values (default: true)
    % OUTPUTS
    %   f:        The RAD feature value, or NaN (with a warning) when either half
    %             of the median split has fewer than two values or no spread, as
    %             happens for short or heavily tied (e.g., integer-valued) series
    %
    % Note: the argument order was RAD(x, center, tau) before it was unified
    % with the Python and Julia versions.
    if nargin < 2 || isempty(tau)
        tau = 1;
    end
    if nargin < 3 || isempty(doCentre)
        doCentre = true;
    end

    if ~isscalar(tau) || tau < 1 || tau ~= round(tau)
        error('RAD:tau', 'tau must be a positive integer (argument order: RAD(x, tau, doCentre))');
    end
    if ~isscalar(doCentre) || ~(islogical(doCentre) || doCentre == 0 || doCentre == 1)
        error('RAD:doCentre', 'doCentre must be true or false (argument order: RAD(x, tau, doCentre))');
    end
    if ~isvector(x)
        error('RAD:x', 'x must be a vector');
    end

    x = x(:);

    if doCentre
        x = x - median(x);
        x = abs(x);
    end

    % Delay embed at interval tau, m = 2
    y = x(tau + 1:end);
    x = x(1:end - tau);

    % Median split
    subMedians = x < median(x);
    if sum(subMedians) < 2 || sum(~subMedians) < 2
        warning('RAD:degenerate', 'A half of the median split has fewer than two values; returning NaN');
        f = NaN;
        return
    end
    superMedianSD = std(x(~subMedians));
    subMedianSD = std(x(subMedians));
    if superMedianSD == 0 || subMedianSD == 0
        warning('RAD:degenerate', 'A half of the median split has no spread (tied values); returning NaN');
        f = NaN;
        return
    end

    % Properties of the auto-density
    sigma_dx = std(y - x);
    densityDifference = 1 ./ superMedianSD - 1 ./ subMedianSD;

    f = sigma_dx .* densityDifference;
end
