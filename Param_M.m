function [param, gridSm] = Param_M(param)
%==========================================================================
% PARAM_M - Generate parameters specific to male single agents
%==========================================================================
% This function sets parameters specific to male agents and builds
% the corresponding state space grids.
%
% Inputs:
%   param - Base parameter structure from Param_Gen()
%
% Outputs:
%   param  - Parameter structure with added male-specific parameters (param.pm)
%   gridSm - State space grids for male single agents
%
% Male-specific features:
%   - Higher average productivity levels
%   - Lower job separation rate than females (paper Table 2)
%   - Different search intensities by employment status
%   - Lower out-of-labor-force benefits (reflecting lower childcare responsibilities)
%   - Gender-specific fertility rates
%
% Reference: Morazzoni, M. and Smirnov, D. (2025)
%==========================================================================

%==========================================================================
%% PRODUCTIVITY PARAMETERS
%==========================================================================

% Productivity process parameters
pm.a_sigma = 0.055;%0.085;             % Standard deviation of productivity shocks
pm.mu = 0.0095;%0.007;                  % Productivity drift intensity when employed  
pm.mu_u = pm.mu * 2.0;%pm.mu * 2.5;          % Productivity drift intensity when unemployed (faster depreciation)
pm.amean = 1.350;%3.0 * 0.90;          % Mean productivity level (higher than females)

% Entry-level productivity distribution
pm.aergmean = 0.0;%1.00;             % Mean of entry productivity (log scale)
pm.aergvar = 0.805;%0.205;              % Variance of entry productivity (log scale)

%==========================================================================
%% LABOR MARKET PARAMETERS
%==========================================================================

% Job destruction and search
pm.delta = 0.012;               % Monthly job separation rate (paper Table 2: lower than females)

% Search intensity by employment status
pm.lambda = 0;                  % Search intensity when employed (on-the-job search)
pm.lambda_u = 0.55;               % Search intensity when unemployed (paper Table 2: higher than females)
pm.lambda_o = 50;              % Search intensity when out of labor force
pm.lambda_oE = 0.095;            % OLF job arrival rate (paper Table 2)

%==========================================================================
%% OUT-OF-LABOR-FORCE BENEFITS
%==========================================================================
% These parameters reflect that males typically have lower childcare
% responsibilities, resulting in lower OLF benefits overall

pm.boutmin = param.bmin*1.385;
pm.boutmax = param.bmin*2;
pm.outscale=0.10;
pm.bout = @(x) min(pm.boutmin + pm.outscale.*exp(x), pm.boutmax);%2.30;                
pm.boutminc = pm.boutmin - 0.60;%1.0;      
pm.boutc = @(x) min(pm.boutminc + pm.outscale.*exp(x), pm.boutmax);%2.30;                

% pm.bout = param.bmin+6.0;%2.30;                
% pm.boutc = pm.bout + 1.0;%1.0;      

%==========================================================================
%% FERTILITY PARAMETERS
%==========================================================================

% Fertility rates for single males
pm.pcs = 0.0004;%0.0004;                % Monthly probability of having children (single males)

%==========================================================================
%% BUILD MALE SINGLE AGENT GRIDS
%==========================================================================

% Build state space grids using the single grid builder
[gridSm, pm] = BuildSingleGrid(param, pm);

%==========================================================================
%% ENTRY DISTRIBUTION
%==========================================================================

% Distribution of productivity for newly entering male agents
% Log-normal, a0 ~ logN(aergmean, aergvar), consistent with BuildDoubleGrid.m
pm.a.erg = lognpdf(gridSm.a.a + gridSm.a.da, pm.aergmean, sqrt(pm.aergvar));
pm.a.erg = pm.a.erg ./ sum(pm.a.erg);  % Normalize to sum to 1

%==========================================================================
%% STORE IN PARAMETER STRUCTURE
%==========================================================================

% Add all male parameters to the main parameter structure
param.pm = pm;

end