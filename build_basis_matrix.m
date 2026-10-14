function Phi = build_basis_matrix(M, d)
% BUILD_BASIS_MATRIX  Construct a B-spline basis matrix.
%
%   Phi = build_basis_matrix(M, d)
%
%   Inputs
%   ------
%   M : number of design (evaluation) points on [0, 1]
%   d : number of B-spline basis functions
%       For order-p splines with n_knot interior knots: d = n_knot + p
%       Default order is 4 (cubic).  To get d bases you need d-p+2
%       equally-spaced breakpoints (including the two endpoints).
%
%   Output
%   ------
%   Phi : M-by-d matrix; Phi(i,k) = k-th basis function evaluated at
%         the i-th design point.
%
%   Example (replicates the paper setting: M=200, d=43)
%   -------
%       Phi = build_basis_matrix(200, 43);
%       figure; plot(Phi); title('B-spline bases');
%
%   Notes
%   -----
%   * Uses 4th-order (cubic) B-splines, matching the paper description
%     "four-order B-spline bases with 41 equally spaced knots".
%   * The knot vector has clamped (repeated) boundary knots so that the
%     basis spans the full domain and sums to 1 at every interior point.
%   * Design points are equally spaced on [0, 1].

p      = 4;                          % spline order (degree = p-1 = 3)
t      = linspace(0, 1, M)';        % M evaluation points in [0,1]

% ---- build the knot vector ----------------------------------------
% d bases of order p  =>  n_interior = d - p  interior knots
% total unique breakpoints = n_interior + 2  (including 0 and 1)
n_int  = d - p;                      % number of interior knots
breaks = linspace(0, 1, n_int + 2); % equally spaced breakpoints

% Clamped knot vector: repeat first and last knot p times
knots  = [repmat(breaks(1), 1, p), ...
          breaks(2:end-1), ...
          repmat(breaks(end), 1, p)];
% length of knots = p + n_int + p = 2p + (d-p) = d + p  ✓

% ---- evaluate using the Cox-de Boor recursion ---------------------
Phi = bspline_eval(t, knots, d, p);

end % ================================================================


function B = bspline_eval(t, knots, n, p)
% Evaluate n B-spline basis functions of order p at points t.
% Uses the stable Cox-de Boor recursion.
%
%   B(i, k) = N_{k,p}(t_i)

M = numel(t);
B = zeros(M, n);

% ----- order-1 (piecewise constant) basis -------------------------
for k = 1:n
    B(:, k) = double(t >= knots(k) & t < knots(k+1));
end
% Handle the right boundary: the last point belongs to the last span
B(t == knots(end), end) = 1;

% ----- raise order from 1 up to p via recursion -------------------
for r = 2:p
    B_new = zeros(M, n);
    for k = 1:n
        % left term
        denom1 = knots(k + r - 1) - knots(k);
        if denom1 > 0
            B_new(:, k) = B_new(:, k) + ...
                ((t - knots(k)) / denom1) .* B(:, k);
        end
        % right term  (needs basis k+1 from previous order)
        if k + 1 <= n + 1          % guard: n+1-th column needed
            denom2 = knots(k + r) - knots(k + 1);
            if denom2 > 0 && k < n
                B_new(:, k) = B_new(:, k) + ...
                    ((knots(k + r) - t) / denom2) .* B(:, k + 1);
            elseif denom2 > 0 && k == n
                % last basis: right neighbor is zero by convention
                B_new(:, k) = B_new(:, k);
            end
        end
    end
    B = B_new;
end

end
