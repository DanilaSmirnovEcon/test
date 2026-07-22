function param = Param_Gen()
%==========================================================================
% PARAM_GEN - Generate baseline model parameters
%==========================================================================
% This function initializes all baseline parameters for the joint search
% model with endogenous marriage decisions.
%
% Output:
%   param - Structure containing all model parameters
%
% Parameter Categories:
%   - Preferences and demographics
%   - Labor market parameters  
%   - Marriage market parameters
%   - Fiscal policy parameters
%   - Numerical solution parameters
%   - State space grid parameters
%
% Reference: Morazzoni, M. and Smirnov, D. (2025)
%==========================================================================

%==========================================================================
%% PREFERENCES AND DEMOGRAPHICS
%==========================================================================

% Time preference and risk aversion
param.rho = 0.0034;%-log(1-0.05)/12;    % Monthly discount rate (paper Table 2)
param.nu = 2.00;                   % Coefficient of relative risk aversion
param.zeta = (1/40)/12;         % Monthly retirement probability (40-year careers)

% Life-cycle and demographics
param.T = 0;                    % Terminal condition parameter

%==========================================================================
%% LABOR MARKET PARAMETERS
%==========================================================================

% Firm parameters
param.F = 1;                    % Mass of vacant firms (normalized)

% Wage offer distribution (Pareto)
param.w.min = 10;%1.5;              % Minimum wage offer
param.w.max = 30;%4.0;              % Maximum wage offer
param.w.n = 9;                  % Number of wage grid points
ParetoFlat = 15;%5;  % Shape parameter (smaller = steeper dropoff)

% Productivity parameters
param.a.min = 0.00;%0.75;             % Minimum productivity level
param.a.max = 3.00;%3.5;              % Maximum productivity level  
param.a.n = 5;                  % Number of productivity grid points

% Unemployment benefits (home production)
param.bmin = 15.85;%4.0;              % Minimum unemployment benefit
param.bmax = 30.0;%9.00;              % Maximum unemployment benefit
param.bscale = 1.05;%0.25;            % Benefit scaling parameter
param.b = @(x) min(param.bmin + param.bscale.*exp(x), param.bmax);

% Additional benefits by household type
param.bc = 0;                   % Extra unemployment benefit if has children
param.bm = 0;                   % Extra unemployment benefit if married
param.bs = 0;                   % Extra unemployment benefit if single
param.t = 0;                    % General transfer parameter

% Out-of-labor-force benefits by household type
param.oc = 0;%3.00;                % Extra OLF benefit if has children (= cmin)
param.om = -1.80;%2.80;                % Extra OLF benefit if married
param.os = +1.20;%1.70;                % Extra OLF benefit if single

%==========================================================================
%% MARRIAGE MARKET PARAMETERS
%==========================================================================

% Marriage process
param.m = 0.375;%0.00375;              % Monthly marriage meeting intensity (paper Table 2)
param.always_marry = false;     % Flag for forced marriages (debugging)
param.tremble = 1e-1;         % Trembling hand parameter for marriage decisions

% Fertility parameters
param.pcp = 0.005;%0.005;              % Monthly probability of children for couples
param.partner_kid_prob_male=0.0005;
param.partner_kid_prob_female=0.002;

% Spousal complementarity
param.l = @(af, am) -0.0.*abs(exp(af) - exp(am)).^1.0; % Utility from spousal match quality

%==========================================================================
%% CHILD-RELATED PARAMETERS
%==========================================================================

% Child costs
param.cscale = 0.35;%0.50;            % Child cost scaling parameter
param.cmin = 7.0;%3.00;              % Minimum child cost
param.c = @(x) min(param.cmin + param.cscale.*exp(x), +Inf); % Child cost function
param.fkidwpenalty = 0.86;      % Female wage penalty with children

%==========================================================================
%% FISCAL POLICY PARAMETERS
%==========================================================================

% Tax system parameters
param.lambda_tax = 0.906;%0.906;       % Tax progressivity parameter
param.tau_tax = 0.18;          % Tax curvature parameter
param.averagewage = 61;         % Average wage for tax calculations

% Tax function: progressive tax schedule
param.tax = @(x, param) max((param.lambda_tax .* ...
    ((x./param.averagewage).^(1-param.tau_tax)) .* param.averagewage), 0.5.*x);

% Tax credits parameters
param.pi1_n0 = 7.65/100;        % Tax credit rate (no children, low income)
param.pi1_n1 = 34/100;          % Tax credit rate (with children, low income)
param.pi2_n0 = 7.65/100;        % Tax credit rate (no children, high income)
param.pi2_n1 = 15.98/100;       % Tax credit rate (with children, high income)

% Tax credit thresholds
param.adjustW      = param.averagewage/25;% wage rescaling
% param.adjustH      = 40/param.adjustW;          % hours rescaling
% param.E1_n0 = (6750/(48*40*param.adjustH));  % Lower income threshold (no children)
% param.E1_n1 = (10150/(48*40*param.adjustH)); % Lower income threshold (with children)
% param.E2_n0 = (14250/(48*40*param.adjustH)); % Upper income threshold (no children)
% param.E2_n1 = (24350/(48*40*param.adjustH)); % Upper income threshold (with children)
% param.E3_n0 = (20950/(48*40*param.adjustH)); % Max income threshold (no children)
% param.E3_n1 = (46010/(48*40*param.adjustH)); % Max income threshold (with children)
param.adjustWScale = param.averagewage/45000;
param.E1_n0 = (6750*param.adjustWScale);  % Lower income threshold (no children)
param.E1_n1 = (10150*param.adjustWScale); % Lower income threshold (with children)
param.E2_n0 = (14200*param.adjustWScale); % Upper income threshold (no children)
param.E2_n1 = (24350*param.adjustWScale); % Upper income threshold (with children)
param.E3_n0 = (20950*param.adjustWScale); % Max income threshold (no children)
param.E3_n1 = (46010*param.adjustWScale); % Max income threshold (with children)


% Tax rebate functions
param.tax_rebait_kid = @(x, param) tax_rebait_kid(x, param);
param.tax_rebait_nokid = @(x, param) tax_rebait_nokid(x, param);

% Disposable income functions
param.nokid_tax = @(x, a, param) param.tax(x, param) + param.tax_rebait_nokid(x, param);
param.onekid_tax = @(x, a, param) param.tax(x, param) + param.tax_rebait_kid(x, param);  
param.no_tax = @(x, a, param) x;  % No tax function (for comparison)

%==========================================================================
%% EMPIRICAL TARGETS FOR CALIBRATION
%==========================================================================

% Marriage and demographics
param.moment.mfrac = 0.61047;       % Target fraction of married individuals
param.moment.unemp_f = 0.055695;    % Target female unemployment rate
param.moment.unemp_m = 0.055695;    % Target male unemployment rate

% Labor market dynamics
param.moment.tenure = 1.5;          % Target average job tenure
param.moment.w5010 = 5;             % Target 50th/10th wage percentile ratio
param.moment.w5010entr = 3;         % Target entry wage dispersion
param.moment.dwvar = 0.5;           % Target wage change variance

%==========================================================================
%% NUMERICAL SOLUTION PARAMETERS
%==========================================================================

% Value function iteration
param.updMult = 0.15;               % Update speed for value functions (conservative)
param.updMult_s = 0.15;             % Update speed for single agent distributions
param.updMult_d = 0.5;              % Update speed for couple distributions
param.updMultf = 0.5;               % Update speed for flows between states
param.updMultLCP = 0.85;            % Update speed for linear complementarity problems

% Convergence criteria
param.maxiter = 100;                % Maximum iterations for outer loops
param.crit = 1e-4;                  % General convergence criterion
param.critC = 1e-6;                 % Strict convergence criterion
param.Delta = 100;                  % Time step parameter

% Grid refinement
param.Sdimmult = 1;                 % Grid refinement multiplier for single agents
param.minValue = 1e-5;              % Minimum value function threshold

%==========================================================================
%% WAGE OFFER DISTRIBUTION SETUP
%==========================================================================

% Construct wage grid for baseline problem
wvec = linspace(param.w.min, param.w.max, param.w.n);
dw = (param.w.max - param.w.min)/(param.w.n - 1);

% Pareto distribution for wage offers
param.f = gppdf(wvec, 0, ParetoFlat, 0);  % Generalized Pareto PDF
param.f = param.f .* dw;                  % Adjust for grid spacing
param.f(end) = param.f(end);             % Ensure proper tail behavior
param.f = param.f./sum(param.f);

% Extended wage grid for refined single agent problems
wvec_long = linspace(param.w.min, param.w.max, (param.w.n-1)*param.Sdimmult+1);
dw_long = (param.w.max - param.w.min)/((param.w.n-1)*param.Sdimmult);
param.f_long = gppdf(wvec_long, 0, ParetoFlat, 0);
param.f_long = param.f_long .* dw_long;
param.f_long = param.f_long ./ sum(param.f_long) .* sum(param.f);
param.f_long(end) = param.f_long(end);
param.f_long = param.f_long./sum(param.f_long);

%==========================================================================
%% NESTED FUNCTIONS FOR TAX REBATES
%==========================================================================

function rebait = tax_rebait_nokid(x, param)
    % Tax rebate schedule for households without children
    rebait = param.pi1_n0 .* x .* (x < param.E1_n0) + ...
             param.pi1_n0 .* param.E1_n0 .* (x >= param.E1_n0) .* (x < param.E2_n0) + ...
             max(0, param.pi1_n0 .* param.E1_n0 - param.pi2_n0 .* (x - param.E2_n0)) .* (x >= param.E2_n0)+...
             0.*(x>=param.E3_n0);
end

function rebait = tax_rebait_kid(x, param)
    % Tax rebate schedule for households with children  
    rebait = param.pi1_n1 .* x .* (x < param.E1_n1) + ...
             param.pi1_n1 .* param.E1_n1 .* (x >= param.E1_n1) .* (x < param.E2_n1) + ...
             max(0, param.pi1_n1 .* param.E1_n1 - param.pi2_n1 .* (x - param.E2_n1)) .* (x >= param.E2_n1)+...
             0.*(x>=param.E3_n1);
end

end