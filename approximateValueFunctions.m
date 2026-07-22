function valueFunc = approximateValueFunctions(solutionSP, gridSf, gridSm, gridP, param)
    % APPROXIMATEVALUEFUNCTIONS - Approximate value functions with Chebyshev polynomials
    % 
    % This function takes the solved value functions on grids and approximates them
    % using Chebyshev polynomials, returning function handles that can be evaluated
    % at any point within the state space.
    %
    % Inputs:
    %   solutionSP - Solution structure containing value functions
    %   gridSf     - Single female grid
    %   gridSm     - Single male grid  
    %   gridP      - Pair/couple grid
    %   param      - Parameter structure
    %
    % Outputs:
    %   valueFunc  - Structure containing approximated value function handles
    
    % fprintf('Starting approximation of value functions...\n');
    
    % Initialize output structure
    valueFunc = struct();
    
    
    %% SINGLE FEMALE VALUE FUNCTIONS
    % fprintf('Approximating single female value functions...\n');

    valueFunc.single_female = struct();
    
    % Employed (ability_f, wage_f)
    if isfield(solutionSP, 'JflValue')
        valueFunc.single_female.employed_nokid = ...
            approximate2D_qspline(solutionSP.JflValue, gridSf.a.a, gridSf.w.w);
    end
    if isfield(solutionSP, 'JfcValue')
        valueFunc.single_female.employed_kid = ...
            approximate2D_qspline(solutionSP.JfcValue, gridSf.a.a, gridSf.w.w);
    end
    
    % Unemployed / OLF (ability_f)
    if isfield(solutionSP, 'UflValue')
        valueFunc.single_female.unemployed_nokid = ...
            approximate1D_qspline(solutionSP.UflValue, gridSf.a.a);
    end
    if isfield(solutionSP, 'UfcValue')
        valueFunc.single_female.unemployed_kid = ...
            approximate1D_qspline(solutionSP.UfcValue, gridSf.a.a);
    end
    if isfield(solutionSP, 'OflValue')
        valueFunc.single_female.olf_nokid = ...
            approximate1D_qspline(solutionSP.OflValue, gridSf.a.a);
    end
    if isfield(solutionSP, 'OfcValue')
        valueFunc.single_female.olf_kid = ...
            approximate1D_qspline(solutionSP.OfcValue, gridSf.a.a);
    end
    
    %% SINGLE MALE VALUE FUNCTIONS
    % fprintf('Approximating single male value functions...\n');
    
    valueFunc.single_male = struct();
    
    % Employed (ability_m, wage_m)
    if isfield(solutionSP, 'JmlValue')
        valueFunc.single_male.employed_nokid = ...
            approximate2D_qspline(solutionSP.JmlValue, gridSm.a.a, gridSm.w.w);
    end
    if isfield(solutionSP, 'JmcValue')
        valueFunc.single_male.employed_kid = ...
            approximate2D_qspline(solutionSP.JmcValue, gridSm.a.a, gridSm.w.w);
    end
    
    % Unemployed / OLF (ability_m)
    if isfield(solutionSP, 'UmlValue')
        valueFunc.single_male.unemployed_nokid = ...
            approximate1D_qspline(solutionSP.UmlValue, gridSm.a.a);
    end
    if isfield(solutionSP, 'UmcValue')
        valueFunc.single_male.unemployed_kid = ...
            approximate1D_qspline(solutionSP.UmcValue, gridSm.a.a);
    end
    if isfield(solutionSP, 'OmlValue')
        valueFunc.single_male.olf_nokid = ...
            approximate1D_qspline(solutionSP.OmlValue, gridSm.a.a);
    end
    if isfield(solutionSP, 'OmcValue')
        valueFunc.single_male.olf_kid = ...
            approximate1D_qspline(solutionSP.OmcValue, gridSm.a.a);
    end
    
    %% COUPLE VALUE FUNCTIONS
    % fprintf('Approximating couple value functions...\n');
    
    valueFunc.couple = struct();
    
    % Both employed: JJ (ability_f, wage_f, ability_m, wage_m)
    if isfield(solutionSP, 'JJlValue')
        valueFunc.couple.both_employed_nokid = ...
            approximate4D_qspline(solutionSP.JJlValue, ...
                gridP.af.a, gridP.wf.w, gridP.am.a, gridP.wm.w);
    end
    if isfield(solutionSP, 'JJcValue')
        valueFunc.couple.both_employed_kid = ...
            approximate4D_qspline(solutionSP.JJcValue, ...
                gridP.af.a, gridP.wf.w, gridP.am.a, gridP.wm.w);
    end
    
    % Female employed, male unemployed: JU (ability_f, wage_f, ability_m)
    if isfield(solutionSP, 'JUlValue')
        valueFunc.couple.f_employed_m_unemployed_nokid = ...
            approximate3D_qspline(solutionSP.JUlValue, ...
                gridP.af.a, gridP.wf.w, gridP.am.a);
    end
    if isfield(solutionSP, 'JUcValue')
        valueFunc.couple.f_employed_m_unemployed_kid = ...
            approximate3D_qspline(solutionSP.JUcValue, ...
                gridP.af.a, gridP.wf.w, gridP.am.a);
    end
    
    % Female unemployed, male employed: UJ (ability_f, ability_m, wage_m)
    if isfield(solutionSP, 'UJlValue')
        valueFunc.couple.f_unemployed_m_employed_nokid = ...
            approximate3D_qspline(solutionSP.UJlValue, ...
                gridP.af.a, gridP.am.a, gridP.wm.w);
    end
    if isfield(solutionSP, 'UJcValue')
        valueFunc.couple.f_unemployed_m_employed_kid = ...
            approximate3D_qspline(solutionSP.UJcValue, ...
                gridP.af.a, gridP.am.a, gridP.wm.w);
    end
    
    % Both unemployed: UU (ability_f, ability_m)
    if isfield(solutionSP, 'UUlValue')
        valueFunc.couple.both_unemployed_nokid = ...
            approximate2D_qspline(solutionSP.UUlValue, gridP.af.a, gridP.am.a);
    end
    if isfield(solutionSP, 'UUcValue')
        valueFunc.couple.both_unemployed_kid = ...
            approximate2D_qspline(solutionSP.UUcValue, gridP.af.a, gridP.am.a);
    end
    
    % Female OLF, male employed: OJ (ability_f, ability_m, wage_m)
    if isfield(solutionSP, 'OJlValue')
        valueFunc.couple.f_olf_m_employed_nokid = ...
            approximate3D_qspline(solutionSP.OJlValue, ...
                gridP.af.a, gridP.am.a, gridP.wm.w);
    end
    if isfield(solutionSP, 'OJcValue')
        valueFunc.couple.f_olf_m_employed_kid = ...
            approximate3D_qspline(solutionSP.OJcValue, ...
                gridP.af.a, gridP.am.a, gridP.wm.w);
    end
    
    % Female employed, male OLF: JO (ability_f, wage_f, ability_m)
    if isfield(solutionSP, 'JOlValue')
        valueFunc.couple.f_employed_m_olf_nokid = ...
            approximate3D_qspline(solutionSP.JOlValue, ...
                gridP.af.a, gridP.wf.w, gridP.am.a);
    end
    if isfield(solutionSP, 'JOcValue')
        valueFunc.couple.f_employed_m_olf_kid = ...
            approximate3D_qspline(solutionSP.JOcValue, ...
                gridP.af.a, gridP.wf.w, gridP.am.a);
    end
    
    % Female OLF, male unemployed: OU (ability_f, ability_m)
    if isfield(solutionSP, 'OUlValue')
        valueFunc.couple.f_olf_m_unemployed_nokid = ...
            approximate2D_qspline(solutionSP.OUlValue, gridP.af.a, gridP.am.a);
    end
    if isfield(solutionSP, 'OUcValue')
        valueFunc.couple.f_olf_m_unemployed_kid = ...
            approximate2D_qspline(solutionSP.OUcValue, gridP.af.a, gridP.am.a);
    end
    
    % Female unemployed, male OLF: UO (ability_f, ability_m)
    if isfield(solutionSP, 'UOlValue')
        valueFunc.couple.f_unemployed_m_olf_nokid = ...
            approximate2D_qspline(solutionSP.UOlValue, gridP.af.a, gridP.am.a);
    end
    if isfield(solutionSP, 'UOcValue')
        valueFunc.couple.f_unemployed_m_olf_kid = ...
            approximate2D_qspline(solutionSP.UOcValue, gridP.af.a, gridP.am.a);
    end
    
    % Both OLF: OO (ability_f, ability_m)
    if isfield(solutionSP, 'OOlValue')
        valueFunc.couple.both_olf_nokid = ...
            approximate2D_qspline(solutionSP.OOlValue, gridP.af.a, gridP.am.a);
    end
    if isfield(solutionSP, 'OOcValue')
        valueFunc.couple.both_olf_kid = ...
            approximate2D_qspline(solutionSP.OOcValue, gridP.af.a, gridP.am.a);
    end

    
    % fprintf('Chebyshev approximation completed successfully.\n');
    % fprintf('Total value functions approximated: %d singles + %d couples\n', ...
    %     length(fieldnames(valueFunc.single_female)) + length(fieldnames(valueFunc.single_male)), ...
    %     length(fieldnames(valueFunc.couple)));
end

%% APPROXIMATION FUNCTIONS BY DIMENSION
function f = approximate1D_qspline(V_data, x1)
% Quadratic (degree 2) spline on 1D uniform grid
% V_data: [n1,1] values on x1
% x1    : vector of grid points (uniform, increasing)

    x1 = x1(:).';                 % row vector for spapi convenience
    V  = V_data(:).';             % row vector

    try
        % order = degree+1 = 3  -> quadratic spline
        sp = spapi(3, x1, V);     % Spline Toolbox
        f  = @(xq) fnval(sp, xq);
    catch
        % Fallback (no Spline Toolbox): use high-quality cubic (makima)
        F = griddedInterpolant(x1, V, 'makima', 'none');
        f  = @(xq) F(xq);
    end
end
function f = approximate2D_qspline(V_data, x1, x2)
% Quadratic (degree 2) tensor-product spline on 2D uniform grid
% V_data: [n1,n2] with first dim  x1, second  x2
% x1,x2 : vectors of grid points (uniform, increasing)

    x1 = x1(:).'; x2 = x2(:).';   % row vectors for spapi
    V  = V_data;                  % keep [n1,n2] layout (ndgrid convention)

    try
        sp = spapi({3,3}, {x1,x2}, V);  % quadratic in each dim
        % Evaluate at paired points or grids:
        f  = @(x1q,x2q) local_eval_tp(sp, x1q, x2q);
    catch
        F = griddedInterpolant({x1,x2}, V, 'spline', 'none'); % cubic fallback
        f  = @(x1q,x2q) F(x1q, x2q);
    end
end
function f = approximate3D_qspline(V_data, x1, x2, x3)
% Quadratic (degree 2) tensor-product spline on 3D uniform grid
% V_data: [n1,n2,n3] with dims  x1,x2,x3
% x1,x2,x3: vectors (uniform, increasing)

    x1 = x1(:).'; x2 = x2(:).'; x3 = x3(:).';
    V  = V_data;

    try
        sp = spapi({3,3,3}, {x1,x2,x3}, V);
        f  = @(x1q,x2q,x3q) local_eval_tp(sp, x1q, x2q, x3q);
    catch
        F = griddedInterpolant({x1,x2,x3}, V, 'spline', 'none');
        f  = @(x1q,x2q,x3q) F(x1q, x2q, x3q);
    end
end
function f = approximate4D_qspline(V_data, x1, x2, x3, x4)
% Quadratic (degree 2) tensor-product spline on 4D uniform grid
% V_data: [n1,n2,n3,n4] with dims  x1,x2,x3,x4
% x1..x4: vectors (uniform, increasing)

    x1 = x1(:).'; x2 = x2(:).'; x3 = x3(:).'; x4 = x4(:).';
    V  = V_data;

    try
        sp = spapi({3,3,3,3}, {x1,x2,x3,x4}, V);
        f  = @(x1q,x2q,x3q,x4q) local_eval_tp(sp, x1q, x2q, x3q, x4q);
    catch
        F = griddedInterpolant({x1,x2,x3,x4}, V, 'spline', 'none');
        f  = @(x1q,x2q,x3q,x4q) F(x1q, x2q, x3q, x4q);
    end
end




%% UTILITY FUNCTIONS
function Vq = local_eval_tp(sp, varargin)
% Evaluate tensor-product spline 'sp' at query points.
% If all inputs have the same size, treat them as paired points.
% Else, treat them as grid vectors and return a tensor grid.

    sz = cellfun(@size, varargin, 'uni', 0);
    sameShape = all(cellfun(@(s) isequal(s, sz{1}), sz));

    if sameShape
        % Paired points: stack as columns [d x N]
        d = numel(varargin);
        N = numel(varargin{1});
        X = zeros(d, N);
        for k = 1:d
            X(k,:) = varargin{k}(:).';
        end
        Vq = fnval(sp, X);
        Vq = reshape(Vq, size(varargin{1}));
    else
        % Tensor grid: pass as cell, fnval returns grid array
        Vq = fnval(sp, varargin);
    end
end


function x_mapped = mapFromUnitInterval(x_unit, bounds)
    % Map from [-1,1] to [bounds(1), bounds(2)]
    x_mapped = bounds(1) + (bounds(2) - bounds(1)) * (x_unit + 1) / 2;
end

function x_unit = mapToUnitInterval(x, bounds)
    % Map from [bounds(1), bounds(2)] to [-1,1]
    x_unit = 2 * (x - bounds(1)) / (bounds(2) - bounds(1)) - 1;
end

%% CHEBYSHEV FITTING FUNCTIONS

function coeffs = fitChebyshev1D(V_data, order)
    % Fit 1D Chebyshev polynomial coefficients using discrete cosine transform
    coeffs = dct(V_data(1:order)) / (order/2);
    coeffs(1) = coeffs(1) / 2;
end

function coeffs = fitChebyshev2D(V_data, order)
    % Fit 2D Chebyshev polynomial coefficients
    coeffs = zeros(order, order);

    % precompute Chebyshev nodes (length "order")
    x_nodes = cos(pi * (0:order-1) / (order-1));
    y_nodes = x_nodes;

    for i = 1:order
        for j = 1:order
            sum_val = 0;
            norm_factor = 0;

            for k = 1:order
                for l = 1:order
                    T_i = chebyshevT(i-1, x_nodes(k));
                    T_j = chebyshevT(j-1, y_nodes(l));

                    % --- MATLAB-friendly endpoint weights (Clenshaw�Curtis) ---
                    if k==1 || k==order
                        weight_k = 0.5;
                    else
                        weight_k = 1.0;
                    end
                    if l==1 || l==order
                        weight_l = 0.5;
                    else
                        weight_l = 1.0;
                    end
                    % -----------------------------------------------------------

                    w = weight_k * weight_l;
                    sum_val   = sum_val   + V_data(k,l) * T_i * T_j * w;
                    norm_factor = norm_factor + (T_i*T_i) * (T_j*T_j) * w;
                end
            end

            coeffs(i,j) = sum_val / norm_factor;
        end
    end
end

function coeffs = fitChebyshev3D(V_data, order)
    % Simplified 3D fitting using tensor product
    coeffs = V_data(1:order, 1:order, 1:order);
end

function coeffs = fitChebyshev4D(V_data, order)
    % Simplified 4D fitting using tensor product
    coeffs = V_data(1:order, 1:order, 1:order, 1:order);
end

%% CHEBYSHEV EVALUATION FUNCTIONS

function V = evaluateChebyshev1D(x, coeffs, bounds)
    % Evaluate 1D Chebyshev approximation
    u = mapToUnitInterval(x, bounds);
    
    V = 0;
    for i = 1:length(coeffs)
        V = V + coeffs(i) * chebyshevT(i-1, u);
    end
end

function V = evaluateChebyshev2D(x1, x2, coeffs, bounds1, bounds2)
    % Evaluate 2D Chebyshev approximation
    u1 = mapToUnitInterval(x1, bounds1);
    u2 = mapToUnitInterval(x2, bounds2);
    
    [order1, order2] = size(coeffs);
    V = 0;
    
    for i = 1:order1
        for j = 1:order2
            T_i = chebyshevT(i-1, u1);
            T_j = chebyshevT(j-1, u2);
            V = V + coeffs(i,j) * T_i * T_j;
        end
    end
end

function V = evaluateChebyshev3D(x1, x2, x3, coeffs, bounds1, bounds2, bounds3)
    % Evaluate 3D Chebyshev approximation (simplified)
    u1 = mapToUnitInterval(x1, bounds1);
    u2 = mapToUnitInterval(x2, bounds2);
    u3 = mapToUnitInterval(x3, bounds3);
    
    % Use linear interpolation of coefficients as approximation
    V = interp3(coeffs, u1, u2, u3, 'linear', 0);
end

function V = evaluateChebyshev4D(x1, x2, x3, x4, coeffs, bounds)
    % Evaluate 4D Chebyshev approximation (simplified)
    u1 = mapToUnitInterval(x1, bounds{1});
    u2 = mapToUnitInterval(x2, bounds{2});
    u3 = mapToUnitInterval(x3, bounds{3});
    u4 = mapToUnitInterval(x4, bounds{4});
    
    % Use linear interpolation of coefficients as approximation
    V = interpn(coeffs, u1, u2, u3, u4, 'linear', 0);
end

function T = chebyshevT(n, x)
    % Evaluate Chebyshev polynomial of the first kind using recurrence relation
    if n == 0
        T = ones(size(x));
    elseif n == 1
        T = x;
    else
        T_prev2 = ones(size(x));
        T_prev1 = x;
        for k = 2:n
            T_current = 2 * x .* T_prev1 - T_prev2;
            T_prev2 = T_prev1;
            T_prev1 = T_current;
        end
        T = T_prev1;
    end
end