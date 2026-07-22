function [param, gridSf] = Param_F(param)
%==========================================================================
% PARAM_F - Generate parameters specific to female single agents
%==========================================================================
% This function sets parameters specific to female agents and builds
% the corresponding state space grids.
%
% Inputs:
%   param - Base parameter structure from Param_Gen()
%
% Outputs:
%   param  - Parameter structure with added female-specific parameters (param.pf)
%   gridSf - State space grids for female single agents
%
% Female-specific features:
%   - Different productivity evolution parameters
%   - Different job separation rates
%   - Different search intensities by employment status
%   - Different out-of-labor-force benefits (reflecting childcare roles)
%
% Reference: Morazzoni, M. and Smirnov, D. (2025)
%==========================================================================

%==========================================================================
%% PRODUCTIVITY PARAMETERS
%==========================================================================

% Productivity process parameters
pf.a_sigma = 0.055;%0.085;             % Standard deviation of productivity shocks
pf.mu = 0.0095;%0.007;                  % Productivity drift intensity when employed
pf.mu_u = pf.mu * 2.0;%pf.mu * 2.5;          % Productivity drift intensity when unemployed (faster depreciation)
pf.amean = 1.350;%3.0 * 0.75;         % Mean productivity level (lower than males)

% Entry-level productivity distribution
pf.aergmean = 0.0;%1.00;             % Mean of entry productivity (log scale)
pf.aergvar = 0.805;%0.205;              % Standard deviation of entry productivity

%==========================================================================
%% LABOR MARKET PARAMETERS
%==========================================================================

% Job destruction and search
pf.delta = 0.015;               % Monthly job separation rate (paper Table 2: higher than males)

% Search intensity by employment status
pf.lambda = 0;                  % Search intensity when employed (on-the-job search)
pf.lambda_u = 0.45;               % Search intensity when unemployed (paper Table 2)
pf.lambda_o = 50;              % Search intensity when out of labor force
pf.lambda_oE = 0.105;            % OLF job arrival rate (paper Table 2: higher than males)

%==========================================================================
%% OUT-OF-LABOR-FORCE BENEFITS
%==========================================================================
% These parameters reflect the gendered nature of childcare responsibilities
% and home production, with females receiving higher OLF benefits when 
% they have children (reflecting childcare value)

pf.boutmin = param.bmin*1.80;
pf.boutmax = param.bmin*2;
pf.outscale=0.10;
pf.bout = @(x) min(pf.boutmin + pf.outscale.*exp(x), pf.boutmax);%2.30;                
pf.boutminc = pf.boutmin - 0.50;%1.0;      
pf.boutc = @(x) min(pf.boutminc + pf.outscale.*exp(x), pf.boutmax);%2.30;                

% pf.bout = param.bmin+9.50;%3.6;               
% pf.boutc = pf.bout + 3.0;%2.6;    

%==========================================================================
%% FERTILITY PARAMETERS
%==========================================================================

% Fertility rates for single females
pf.pcs = 0.0004;%0.0004;                % Monthly probability of having children (single females)

%==========================================================================
%% BUILD FEMALE SINGLE AGENT GRIDS
%==========================================================================

% Build state space grids using the single grid builder
[gridSf, pf] = BuildSingleGrid(param, pf);

%==========================================================================
%% ENTRY DISTRIBUTION
%==========================================================================

% Distribution of productivity for newly entering female agents
% Uses normal distribution in levels (not logs) centered at aergmean
pf.a.erg = normpdf(gridSf.a.a, pf.aergmean, sqrt(pf.aergvar)); 
pf.a.erg = pf.a.erg ./ sum(pf.a.erg);  % Normalize to sum to 1

%==========================================================================
%% STORE IN PARAMETER STRUCTURE
%==========================================================================

% Add all female parameters to the main parameter structure
param.pf = pf;

end