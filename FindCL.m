function [CL,ARL,SDRL,info] = FindCL(X,targetARL0,profile_start)
% Robust empirical UCL calibration for stepwise ARL curves.
% Allows negative thresholds and avoids fminbnd on a discontinuous objective.
%
% profile_start:
%   first index at which an alarm is allowed / evaluated.
%   default = 1.

if nargin < 3
    profile_start = 1;
end

X(isnan(X)) = inf;
finiteX = X(isfinite(X));
if isempty(finiteX)
    error('FindCL_fixed:NoFiniteValues', 'Test matrix contains no finite values.');
end

lo = min(finiteX) - max(1, abs(min(finiteX))) * eps;
hi = max(finiteX);

[ARL_lo, ~] = ComputeARL_fixed(X, lo, profile_start);
[ARL_hi, ~] = ComputeARL_fixed(X, hi, profile_start);

info = struct();
info.lower_bound = lo;
info.upper_bound = hi;
info.ARL_at_lower = ARL_lo;
info.ARL_at_upper = ARL_hi;
info.profile_start = profile_start;

if targetARL0 <= ARL_lo
    warning('FindCL_fixed:TargetBelowRange', ...
        ['Target ARL0=%.4f is below the minimum empirical ARL=%.4f. ', ...
         'Returning the smallest feasible threshold.'], ...
        targetARL0, ARL_lo);
    CL = lo;
    [ARL,SDRL] = ComputeARL_fixed(X, CL, profile_start);
    return;
end

if targetARL0 >= ARL_hi
    warning('FindCL_fixed:TargetAboveRange', ...
        ['Target ARL0=%.4f is above the maximum empirical ARL=%.4f. ', ...
         'Returning the largest feasible threshold.'], ...
        targetARL0, ARL_hi);
    CL = hi;
    [ARL,SDRL] = ComputeARL_fixed(X, CL, profile_start);
    return;
end

for iter = 1:80
    mid = (lo + hi) / 2;
    ARL_mid = ComputeARL_fixed(X, mid, profile_start);
    if ARL_mid < targetARL0
        lo = mid;
    else
        hi = mid;
    end
end

[ARL_lo, SDRL_lo] = ComputeARL_fixed(X, lo, profile_start);
[ARL_hi, SDRL_hi] = ComputeARL_fixed(X, hi, profile_start);

if abs(ARL_lo - targetARL0) <= abs(ARL_hi - targetARL0)
    CL = lo;
    ARL = ARL_lo;
    SDRL = SDRL_lo;
else
    CL = hi;
    ARL = ARL_hi;
    SDRL = SDRL_hi;
end

info.final_low = lo;
info.final_high = hi;
info.final_ARL_low = ARL_lo;
info.final_ARL_high = ARL_hi;
end

function [ARL,SDRL] = ComputeARL_fixed(Test,UCL,profile_start)
[rep,T_sum] = size(Test);
profile_start = max(1, min(T_sum, profile_start));

detectP = T_sum * ones(rep,1);
for s = 1:rep
    temp = find(Test(s, profile_start:end) > UCL, 1);
    if ~isempty(temp)
        detectP(s) = temp + profile_start - 1;
    end
end

runlen = detectP - profile_start + 1;
ARL = mean(runlen);
SDRL = std(runlen);
end
