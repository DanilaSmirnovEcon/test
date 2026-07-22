function [gridS, pg] = BuildSingleGrid(param, pg)
%==========================================================================
% BUILDSINGLE GRID - Construct state space grids for single household problems
%==========================================================================
% This function builds the state space grids and transition matrices for
% single agent problems. The state space consists of idiosyncratic productivity
% and current wage offers.
%
% Inputs:
%   param - General model parameters
%   pg    - Gender-specific parameters (pf for females, pm for males)
%
% Outputs:
%   gridS - State space grids and derived objects
%   pg    - Updated parameter structure with grid information
%
% The function discretizes the continuous-time productivity process:
%   Employed:   da = ( - a)dt + dW  
%   Unemployed: da = _u(a_min - a)dt + dW
%
% using an implicit finite difference scheme with upwind differencing.
%
% Reference: Morazzoni, M. and Smirnov, D. (2025)
%==========================================================================

%==========================================================================
%% GRID SETUP
%==========================================================================

% Grid dimensions (refined by Sdimmult factor)
gridS.a.n = (param.a.n - 1) * param.Sdimmult + 1;  % Productivity grid points
gridS.w.n = (param.w.n - 1) * param.Sdimmult + 1;  % Wage grid points

%--------------------------------------------------------------------------
% Productivity grid
%--------------------------------------------------------------------------
gridS.a.min = param.a.min;
gridS.a.max = param.a.max;
gridS.a.a = linspace(gridS.a.min, gridS.a.max, gridS.a.n)';  % Grid vector
gridS.a.aa = gridS.a.a * ones(1, gridS.w.n);                % Grid matrix (a � w)
gridS.a.da = (gridS.a.max - gridS.a.min) / (gridS.a.n - 1); % Grid spacing
gridS.a.da2 = gridS.a.da^2;                                  % Grid spacing squared

%--------------------------------------------------------------------------
% Wage grid  
%--------------------------------------------------------------------------
gridS.w.min = param.w.min;
gridS.w.max = param.w.max;
gridS.w.w = linspace(gridS.w.min, gridS.w.max, gridS.w.n);   % Grid vector
gridS.w.ww = ones(gridS.a.n, 1) * gridS.w.w;                % Grid matrix (a � w)
gridS.w.dw = (gridS.w.max - gridS.w.min) / (gridS.w.n - 1); % Grid spacing

%==========================================================================
%% PRODUCTIVITY PROCESS DISCRETIZATION
%==========================================================================
% The productivity follows different processes based on employment status:
% - When employed: mean reversion toward pg.amean with drift pg.mu
% - When unemployed: drift toward minimum with faster rate pg.mu_u

% Drift terms (from Ito's lemma)
gridS.a.mu = pg.mu * (pg.amean - gridS.a.aa);         % Employed drift
gridS.a.mu_u = pg.mu_u * (param.a.min - gridS.a.aa); % Unemployed drift

% Variance term (constant across states)
gridS.a.var = (pg.a_sigma^2) / 2;  % Variance coefficient for finite difference

%==========================================================================
%% TRANSITION MATRIX DIMENSIONS
%==========================================================================

% State space sizes
pg.NS = gridS.a.n * gridS.w.n;     % Employed states (productivity � wage)
pg.NSall = pg.NS + gridS.a.n;      % Total states (employed + unemployed)

%==========================================================================
%% BUILD TRANSITION MATRICES
%==========================================================================
% These matrices capture the evolution of the productivity process and
% exogenous transitions between employment states

% Initialize sparse transition matrices
pg.CJ = sparse(zeros(pg.NS, pg.NS));           % Employed-to-employed transitions
pg.CU = sparse(zeros(gridS.a.n, gridS.a.n));  % Unemployed-to-unemployed transitions  
pg.CJ_U = sparse(zeros(pg.NS, gridS.a.n));    % Employed-to-unemployed (job loss)

%--------------------------------------------------------------------------
% Build unemployed-to-unemployed transition matrix
%--------------------------------------------------------------------------

for ai = 1:gridS.a.n
    % Current productivity level
    current_a = gridS.a.a(ai);
    
    % Finite difference coefficients for productivity evolution
    % Uses upwind scheme: backward differences for negative drift, forward for positive
    variance_term = gridS.a.var;
    backward_coeff = min(0, gridS.a.mu_u(ai, 1));  % Negative drift coefficient
    forward_coeff = max(0, gridS.a.mu_u(ai, 1));   % Positive drift coefficient
    
    % Fill transition matrix using finite difference approximation
    if current_a > gridS.a.min  % Not at lower boundary
        % Backward difference term
        pg.CU(ai, ai-1) = pg.CU(ai, ai-1) + ...
                          (-backward_coeff/gridS.a.da + variance_term/gridS.a.da2);
        pg.CU(ai, ai) = pg.CU(ai, ai) + ...
                        (backward_coeff/gridS.a.da - variance_term/gridS.a.da2);
    end
    
    if current_a < gridS.a.max  % Not at upper boundary
        % Forward difference term
        pg.CU(ai, ai+1) = pg.CU(ai, ai+1) + ...
                          (forward_coeff/gridS.a.da + variance_term/gridS.a.da2);
        pg.CU(ai, ai) = pg.CU(ai, ai) + ...
                        (-forward_coeff/gridS.a.da - variance_term/gridS.a.da2);
    end
    
    % Retirement flow (exit from the model)
    pg.CU(ai, ai) = pg.CU(ai, ai) - param.zeta;
    
    %----------------------------------------------------------------------
    % Build employed-to-employed transition matrix
    %----------------------------------------------------------------------
    
    for ri = 1:gridS.w.n
        % Linear index for employed state (productivity, wage)
        employed_index = (ri - 1) * gridS.a.n + ai;
        
        % Finite difference coefficients for employed productivity evolution
        backward_coeff = min(0, gridS.a.mu(ai, ri));
        forward_coeff = max(0, gridS.a.mu(ai, ri));
        variance_term = gridS.a.var;
        
        if current_a > gridS.a.min  % Not at lower boundary
            % Backward difference term
            pg.CJ(employed_index, employed_index-1) = ...
                pg.CJ(employed_index, employed_index-1) + ...
                (-backward_coeff/gridS.a.da + variance_term/gridS.a.da2);
            pg.CJ(employed_index, employed_index) = ...
                pg.CJ(employed_index, employed_index) + ...
                (backward_coeff/gridS.a.da - variance_term/gridS.a.da2);
        end
        
        if current_a < gridS.a.max  % Not at upper boundary
            % Forward difference term  
            pg.CJ(employed_index, employed_index+1) = ...
                pg.CJ(employed_index, employed_index+1) + ...
                (forward_coeff/gridS.a.da + variance_term/gridS.a.da2);
            pg.CJ(employed_index, employed_index) = ...
                pg.CJ(employed_index, employed_index) + ...
                (-forward_coeff/gridS.a.da - variance_term/gridS.a.da2);
        end
        
        % Retirement flow (exit from the model)
        pg.CJ(employed_index, employed_index) = ...
            pg.CJ(employed_index, employed_index) - param.zeta;
            
        % Exogenous job separation (employed to unemployed)
        pg.CJ(employed_index, employed_index) = ...
            pg.CJ(employed_index, employed_index) - pg.delta;
        pg.CJ_U(employed_index, ai) = ...
            pg.CJ_U(employed_index, ai) + pg.delta;
    end
end

pg.CO = pg.CU;
%==========================================================================
%% GRID SUMMARY INFORMATION
%==========================================================================

% Store key grid information for later use
gridS.info.n_productivity = gridS.a.n;
gridS.info.n_wages = gridS.w.n; 
gridS.info.n_employed_states = pg.NS;
gridS.info.n_unemployed_states = gridS.a.n;
gridS.info.n_total_states = pg.NSall;

% Display grid information
% fprintf('   Built single agent grid:\n');
% fprintf('     - Productivity points: %d\n', gridS.a.n);
% fprintf('     - Wage points: %d\n', gridS.w.n);
% fprintf('     - Total employed states: %d\n', pg.NS);
% fprintf('     - Total unemployed states: %d\n', gridS.a.n);
% fprintf('     - Grid spacing (productivity): %.4f\n', gridS.a.da);
% fprintf('     - Grid spacing (wages): %.4f\n', gridS.w.dw);

end