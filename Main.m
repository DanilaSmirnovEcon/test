%==========================================================================
% MAIN.M - Labor and Family Dynamics in a Joint-Search Framework
%==========================================================================
% This is the main file that solves the joint search problem on the labor 
% market with endogenous marriage decisions.
%
% The model features:
% - Single agents and couples who search for jobs when unemployed/employed
% - Endogenous marriage decisions for single agents
% - Productivity accumulation on the job and depletion when unemployed
% - Fertility decisions affecting household costs
%
% Reference: Morazzoni, M. and Smirnov, D. (2025)
% "Labor and Family Dynamics in a Joint-Search Framework"
%
% Last updated: August 2025
%==========================================================================

clear all; close all; %clc;

% Display startup message
fprintf('\n');
fprintf('==========================================================================\n');
% fprintf('Labor and Family Dynamics in a Joint-Search Framework\n');
% fprintf('Morazzoni & Smirnov (2025)\n');
fprintf('==========================================================================\n');
% fprintf('Starting model solution...\n\n');

%==========================================================================
%% 1. PARAMETER INITIALIZATION
%==========================================================================
% fprintf('Step 1: Initializing parameters...\n');

% Generate baseline parameters
param = Param_Gen();

% Build grids for female and male single agents
[param, gridSf] = Param_F(param);
[param, gridSm] = Param_M(param);

% Build joint grids for couples
[param, gridP] = BuildDoubleGrid(param, true);

% Set total single agents count
param.NSall = param.pf.NSall + param.pm.NSall;

% fprintf('   - Parameter structures initialized\n');
% fprintf('   - State space grids constructed\n');

%==========================================================================
%% 2. SOLVE SINGLE AGENT PROBLEMS
%==========================================================================
% fprintf('\nStep 2: Solving single agent problems...\n');

% Set convergence parameters for single agent problems
param.updMult = 0.90;          % Update multiplier for value function iteration
param.critC = 1e-6;            % Convergence criterion
param.maxiter = 5000;          % Maximum iterations

% Solve for female single agents
% fprintf('   - Solving female single agent problem...\n');
solutionSf = SolveSingleProblem(param,param.pf, gridSf);

% Solve for male single agents  
% fprintf('   - Solving male single agent problem...\n');
solutionSm = SolveSingleProblem(param,param.pm, gridSm);

% fprintf('   - Single agent problems solved successfully\n');

%==========================================================================
%% 3. SOLVE ENDOGENOUS MARRIAGE PROBLEM
%==========================================================================
fprintf('\nStep 3: Solving endogenous marriage problem...\n');

% Set parameters for marriage market equilibrium
param.critC = 1e-6;            % Convergence criterion for marriage market
param.updMult = 0.75;           % Update multiplier for value functions
param.updMult_s = 0.5;%0.95;        % Update multiplier for single agent distributions
param.updMult_d = 0.75;%1.0;         % Update multiplier for couple distributions  
param.updMultf = 0.75;          % Update multiplier for flows
param.maxiter = 200;           % Maximum iterations for marriage equilibrium

% Option to load previous solution (set to false for clean run)
load_previous_solution = false;
dispiter = true;

% Solve the endogenous marriage problem
% fprintf('   - Computing marriage market equilibrium...\n');
tic;
solutionSP = SolveEndogProblem(solutionSf, solutionSm, load_previous_solution, ...
                               dispiter, param, gridSf, gridSm, gridP);
toc;

% Save results
% fprintf('   - Saving solution...\n');
save('solutionSP.mat', 'solutionSP', 'param', 'gridSf', 'gridSm', 'gridP','-v7.3');
% fprintf('   - Solution saved to solutionSP.mat\n');

%%
fprintf('\nStep 4: Simulate the distribution...\n');
valueFunc = approximateValueFunctions(solutionSP, gridSf, gridSm, gridP, param);
% Set simulation parameters
simParam = struct( ...
  'N',                 1000, ...
  'seed',              12345, ...
  'nWorkers',          [], ...           % [] => default pool
  'coarse_periods',    69, ...           % years
  'coarse_step_months',12, ...           % 12 months per step
  'burn_periods',      200, ...          % steps (check that the simul plots not trended over the simulation)
  'burn_step_months',  0.25, ...         % 1 week per step
  'sim_periods',       480, ...           % steps
  'sim_step_months',   0.25 ...          % 1 week per step
);
tic;
agentPanel= SimulatePanel(solutionSP, valueFunc, gridSf, gridSm, gridP, param, simParam);
toc;

%==========================================================================
%% 5. DISPLAY RESULTS
%==========================================================================
fprintf('\nStep 5: Model Results Summary\n');
fprintf('==========================================================================\n');
moments_stats = compute_moments(agentPanel, simParam, param);


% Display targeted moments (used in calibration)
display_targeted_moments(moments_stats, param);

% Display untargeted moments (model predictions)
display_untargeted_moments(moments_stats, param);

% Generate plots and figures
fprintf('   - Generating plots and figures...\n');
% PlotSimul(agentPanel, param, valueFunc);
% PlotSimul_Extras(agentPanel, param);  

fprintf('\n==========================================================================\n');
fprintf('Model solution completed successfully!\n');
fprintf('==========================================================================\n\n');
close all
%==========================================================================
%% HELPER FUNCTIONS FOR DISPLAYING RESULTS
%==========================================================================

function display_targeted_moments(moments_stats,param)
    if nargin < 2 || ~isstruct(param), param = struct(); end
    fprintf('\nTARGETED MOMENTS (Calibration Targets)\n');
    fprintf('--------------------------------------\n');
    
    % Marriage and demographics
    fprintf('Married share of population:           %6.2f%%  (target: 50%%)\n', ...
            100*moments_stats.married.all);
    
    fprintf('\n');
    
    % Labor market transitions
    fprintf('UE transition rate - Female:           %6.2f%%  (target: 23%%)\n', ...
            100*moments_stats.transitions.female.UE);
    fprintf('UE transition rate - Male:             %6.2f%%  (target: 26%%)\n', ...
            100*moments_stats.transitions.male.UE);
    fprintf('EU transition rate - Female:           %6.2f%%  (target: 1.2%%)\n', ...
            100*moments_stats.transitions.female.EU);
    fprintf('EU transition rate - Male:             %6.2f%%  (target: 1.5%%)\n', ...
            100*moments_stats.transitions.male.EU);
    fprintf('OE transition rate - Female:           %6.2f%%  (target: 4.0%%)\n', ...
            100*moments_stats.transitions.female.OE);
    fprintf('OE transition rate - Male:             %6.2f%%  (target: 5.3%%)\n', ...
            100*moments_stats.transitions.male.OE);
    
    fprintf('\n');
    
    % Out of labor force shares
    % fprintf('Out of LF - Single W K:                   %6.2f%%  (target: 12.7%%)\n', ...
    %         100*moments_stats.never_married.out_of_lf/moments_stats.never_married.all)
    % fprintf('Out of LF - Single W NK:                   %6.2f%%  (target: 13.1%%)\n', ...
    %         100*moments_stats.married.out_of_lf/moments_stats.married.all)
    % fprintf('Out of LF - Single M K:                    %6.2f%%  (target: 20.9%%)\n', ...
    %         100*moments_stats.women.out_of_lf/moments_stats.women.all)
    % fprintf('Out of LF - Single M NK:                      %6.2f%%  (target: 5.0%%)\n', ...
    %         100*moments_stats.men.out_of_lf/moments_stats.men.all)
    % fprintf('Out of LF - Married W K:                   %6.2f%%  (target: 12.7%%)\n', ...
    %         100*moments_stats.never_married.out_of_lf/moments_stats.never_married.all)
    % fprintf('Out of LF - Married W NK:                   %6.2f%%  (target: 13.1%%)\n', ...
    %         100*moments_stats.married.out_of_lf/moments_stats.married.all)
    % fprintf('Out of LF - Married M K:                    %6.2f%%  (target: 20.9%%)\n', ...
    %         100*moments_stats.women.out_of_lf/moments_stats.women.all)
    % fprintf('Out of LF - Married M NK:                      %6.2f%%  (target: 5.0%%)\n', ...
    %         100*moments_stats.men.out_of_lf/moments_stats.men.all)
    % 
    % fprintf('\n');
    % 
    fprintf('Out of LF - Singles:                   %6.2f%%  (target: 12.7%%)\n', ...
            100*moments_stats.never_married.out_of_lf)
    fprintf('Out of LF - Married:                   %6.2f%%  (target: 13.1%%)\n', ...
            100*moments_stats.married.out_of_lf)
    fprintf('Out of LF - Female:                    %6.2f%%  (target: 20.9%%)\n', ...
            100*moments_stats.women.out_of_lf)
    fprintf('Out of LF - Male:                      %6.2f%%  (target: 5.0%%)\n', ...
            100*moments_stats.men.out_of_lf)
    fprintf('Out of LF - W/ Kids:                   %6.2f%%  (target: 14.9%%)\n', ...
            100*moments_stats.haschild.out_of_lf)
    fprintf('Out of LF - W/O Kids:                  %6.2f%%  (target: 10.1%%)\n', ...
            100*moments_stats.nochild.out_of_lf)

    fprintf('\n');
    

    % Wage dynamics
    fprintf('log Wage var (entry):      %6.3f   (target: 0.35)\n', ...
            moments_stats.wageentryvar);
    fprintf('Wage progression (entry/average):      %6.3f   (target: 0.5-0.55)\n', ...
            moments_stats.wageprogression);
    fprintf('Small wage changes (<20%% annually):    %6.2f%%  (target: 60-65%%)\n', ...
            100*moments_stats.wagechangeyear);
    fprintf('Log wage dispersion:                   %6.3f   (target: 0.55-0.65)\n', ...
            moments_stats.wagedisp);
    fprintf('Wage loss per month unemployed:        %6.2f%%  (target: -1.0%%)\n', ...
            100*moments_stats.wagechangeunemployment_per_period_logpct);
    fprintf('CV Job offers:                         %6.2f%%  (target: 30-40%%)\n', ...
            100*moments_stats.stats_wages.cvOffers);
   
    fprintf('\n');

    % Cost of a child (low/high by median income) ---
    cc = moments_stats.child_cost;
    fprintf('Child cost on average = %.3f\n', cc.mean_share);
    fprintf('  Low-income:  mean share cost = %6.2f%%  (target: 20-30%%) \n', ...
            100*cc.low.mean_share);
    fprintf('  High-income: mean share cost = %6.2f%%  (target: 10-15%%)  \n', ...
            100*cc.high.mean_share);

    % Unemp Benefit (low/high by median income) ---
    fprintf('Unemp benefit/income:      %6.3f   (target: 0.35-0.4)\n', ...
            moments_stats.unemp_benefits.all_share_of_median_income);
    fprintf('Unemp benefit/income (low):      %6.3f   (target: 0.45)\n', ...
            moments_stats.unemp_benefits.lowinc_mean_share);
    fprintf('Unemp benefit/income (high):      %6.3f   (target: 0.35)\n', ...
            moments_stats.unemp_benefits.highinc_mean_share);
end

function display_untargeted_moments(moments_stats,param)
    if nargin < 2 || ~isstruct(param), param = struct(); end
    if ~isfield(param,'adjustW') || isempty(param.adjustW), param.adjustW = 1; end

    fprintf('\nUNTARGETED MOMENTS (Model Predictions)\n');
    fprintf('--------------------------------------\n');
    
    % Additional labor market transitions
    fprintf('EO transition rate - Female:           %6.2f%%  (target: 2.0%%)\n', ...
            100*moments_stats.transitions.female.EO);
    fprintf('EO transition rate - Male:             %6.2f%%  (target: 1.5%%)\n', ...
            100*moments_stats.transitions.male.EO);
    fprintf('UO transition rate - Female:           %6.2f%%  (target: 22-25%%)\n', ...
            100*moments_stats.transitions.female.UO);
    fprintf('UO transition rate - Male:             %6.2f%%  (target: 18-20%%)\n', ...
            100*moments_stats.transitions.male.UO);
    
    fprintf('\n');
    
    % Couple employment patterns
    fprintf('Couples - Both Employed (EE):          %6.2f%%  (target: 55%%)\n', ...
            100*moments_stats.couplesemploymentshares.EE);
    fprintf('Couples - Mixed Employment (UE/EU):    %6.2f%%  (target: 3-4%%)\n', ...
            100*moments_stats.couplesemploymentshares.EU);
    fprintf('Couples - Both Unemployed (UU):        %6.2f%%  (target: 0.39%%)\n', ...
            100*moments_stats.couplesemploymentshares.UU);
    fprintf('Couples - Mixed w/ Out of LF (EO/OE):  %6.2f%%  (target: 30-35%%)\n', ...
            100*moments_stats.couplesemploymentshares.EO);
    fprintf('Couples - Both Out of LF (OO):         %6.2f%%  (target: 7-9%%)\n', ...
            100*moments_stats.couplesemploymentshares.OO);
    fprintf('Couples - Mixed Unemp/OutLF (UO/OU):   %6.2f%%  (target: 0.8%%)\n', ...
            100*moments_stats.couplesemploymentshares.UO);
    
    fprintf('\n');

    % % Single employment patterns
    % fprintf('Single Female - Employed (E):                %6.2f%%  (target: 73%%)\n', ...
    %         100*moments_stats.singlesemploymentshares.female.E);
    % fprintf('Single Female - Unemployed (U):              %6.2f%%  (target: 5%%)\n', ...
    %         100*moments_stats.singlesemploymentshares.female.U);
    % fprintf('Single Female - Out of LF (O):               %6.2f%%  (target: 22%%)\n', ...
    %         100*moments_stats.singlesemploymentshares.female.O);    
    % fprintf('Single Male - Employed (E):                  %6.2f%%  (target: 85%%)\n', ...
    %         100*moments_stats.singlesemploymentshares.male.E);
    % fprintf('Single Male - Unemployed (U):                %6.2f%%  (target: 5%%)\n', ...
    %         100*moments_stats.singlesemploymentshares.male.U);
    % fprintf('Single Male - Out of LF (O):                 %6.2f%%  (target: 10%%)\n', ...
    %         100*moments_stats.singlesemploymentshares.male.O);
    % fprintf('\n');
    
    
    % Unemployment rates by group
    fprintf('Unemployment rate - Singles:           %6.2f%%  (target: 5-6%%)\n', ...
            100*moments_stats.never_married.unemployed);%
    fprintf('Unemployment rate - Married:           %6.2f%%  (target: 3.5%%)\n', ...
            100*moments_stats.married.unemployed);%
    fprintf('Unemployment rate - Female:            %6.2f%%  (target: 3.7%%)\n', ...
            100*moments_stats.women.unemployed);%
    fprintf('Unemployment rate - Male:              %6.2f%%  (target: 5.1%%)\n', ...
            100*moments_stats.men.unemployed);%
    fprintf('Unemployment rate - W/ Kids:           %6.2f%%  (target: 3.9%%)\n', ...
            100*moments_stats.haschild.unemployed);%
    fprintf('Unemployment rate - W/O Kids:          %6.2f%%  (target: 5.2%%)\n', ...
            100*moments_stats.nochild.unemployed);%

    fprintf('\n');
    
    % Average wages by group  
    fprintf('Average wage - All:                    %6.2f   (target: 25.0)\n', ...
            moments_stats.stats_wages.all/param.adjustW);
    fprintf('Average wage - Singles:                %6.2f   (target: 23.6)\n', ...
            moments_stats.stats_wages.never_married/param.adjustW);
    fprintf('Average wage - Married:                %6.2f   (target: 26.7)\n', ...
            moments_stats.stats_wages.married/param.adjustW);
    fprintf('Average wage - Female:                 %6.2f   (target: 22.5)\n', ...
            moments_stats.stats_wages.women/param.adjustW);
    fprintf('Average wage - Male:                   %6.2f   (target: 27.5)\n', ...
            moments_stats.stats_wages.men/param.adjustW);
    fprintf('Average wage - W/ Kids:                %6.2f   (target: 24)\n', ...
            moments_stats.stats_wages.haschild/param.adjustW);
    fprintf('Average wage - W/O Kids:               %6.2f   (target: 26.5)\n', ...
            moments_stats.stats_wages.nochild/param.adjustW);

    fprintf('\n');
    
    % Assortative matching
    fprintf('Wage correlation in couples:           %6.2f%%  (target: 35-45%%)\n', ...
            100*moments_stats.corr_wages);
    fprintf('Ability correlation in couples:        %6.2f%%  (model prediction)\n', ...
            100*moments_stats.corr_ability);
 
    
    % --- Average age at (first) marriage ---
    am = moments_stats.age_at_marriage;
    fprintf('Average age at (first) marriage (years):  Women = %.2f,  Men = %.2f,  All = %.2f\n', ...
            am.women_years, am.men_years, am.all_years);
    fprintf('=================================================================\n');


end