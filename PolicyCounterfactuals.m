%==========================================================================
% POLICYCOUNTERFACTUALS.M - Section 6 policy experiments
%==========================================================================
% Sets up the two fiscal counterfactuals from Section 6 of Morazzoni &
% Smirnov, "Labor and Family Dynamics in a Joint-Search Framework":
%   1. ub_plus15      - 15% higher unemployment-benefit baseline rate b1
%   2. credit_plus20  - 20% higher income caps for the child tax credit,
%                       households with a dependent child only (k=1)
%
% For each policy (plus the baseline, solved fresh here for a
% consistent comparison point), this script solves the model once and
% simulates it twice:
%   - agentPanel_endog: endogenous, bilateral marriage decisions (SimulatePanel)
%   - agentPanel_exog:  the exogenous-marriage counterfactual used for the
%                       "without endogenous marriage" comparison in
%                       Table 5/tab:decom (SimulatePanelCounter)
%
% This script is NOT called from Main.m; run it manually. It only wires
% up the parameter changes and the solve/simulate pipeline and reports
% compute_moments for each of the 3 policies x 2 marriage regimes = 6
% runs. It does NOT compute the welfare (consumption-equivalent) or
% fiscal-expenditure comparisons in Tables 6-9 - those were computed
% separately and are not reproduced here.
%
% Expect this to take roughly 3x as long as a single run of Main.m.
%==========================================================================

clear all; close all;

policies = struct( ...
    'name',  {'baseline',  'ub_plus15',        'credit_plus20'}, ...
    'apply', {@(p) p,      @apply_ub_plus15,   @apply_credit_plus20} ...
);

simParam = struct( ...
  'N',                 1000, ...
  'seed',              12345, ...
  'nWorkers',          [], ...
  'coarse_periods',    69, ...
  'coarse_step_months',12, ...
  'burn_periods',      200, ...
  'burn_step_months',  0.25, ...
  'sim_periods',       480, ...
  'sim_step_months',   0.25 ...
);

results = struct();

for k = 1:numel(policies)
    name = policies(k).name;
    fprintf('\n==========================================================================\n');
    fprintf('Policy: %s\n', name);
    fprintf('==========================================================================\n');

    % ---- 1. Parameters and grids ----
    param = Param_Gen();
    param = policies(k).apply(param);
    [param, gridSf] = Param_F(param);
    [param, gridSm] = Param_M(param);
    [param, gridP]  = BuildDoubleGrid(param, true);
    param.NSall = param.pf.NSall + param.pm.NSall;

    % ---- 2. Single-agent problems ----
    param.updMult = 0.90;
    param.critC   = 1e-6;
    param.maxiter = 5000;
    solutionSf = SolveSingleProblem(param, param.pf, gridSf);
    solutionSm = SolveSingleProblem(param, param.pm, gridSm);

    % ---- 3. Endogenous marriage equilibrium ----
    param.critC     = 1e-6;
    param.updMult   = 0.75;
    param.updMult_s = 0.5;
    param.updMult_d = 0.75;
    param.updMultf  = 0.75;
    param.maxiter   = 200;
    solutionSP = SolveEndogProblem(solutionSf, solutionSm, false, true, ...
                                   param, gridSf, gridSm, gridP);
    valueFunc  = approximateValueFunctions(solutionSP, gridSf, gridSm, gridP, param);

    % ---- 4a. Simulate with endogenous marriage ----
    agentPanel_endog = SimulatePanel(solutionSP, valueFunc, gridSf, gridSm, gridP, param, simParam);
    moments_endog = compute_moments(agentPanel_endog, simParam, param);

    % ---- 4b. Simulate with the exogenous-marriage counterfactual ----
    agentPanel_exog = SimulatePanelCounter(solutionSP, valueFunc, gridSf, gridSm, gridP, param, simParam);
    moments_exog = compute_moments(agentPanel_exog, simParam, param);

    % ---- 5. Store and save ----
    results.(name).param         = param;
    results.(name).moments_endog = moments_endog;
    results.(name).moments_exog  = moments_exog;

    save(sprintf('solutionSP_%s.mat', name), 'solutionSP', 'param', 'gridSf', 'gridSm', 'gridP', '-v7.3');
end

save('PolicyCounterfactuals_results.mat', 'results', '-v7.3');

fprintf('\nDone. Results saved to PolicyCounterfactuals_results.mat.\n');
fprintf('See results.<policy>.moments_endog / .moments_exog for wage, unemployment,\n');
fprintf('and marriage-rate moments under each policy and marriage regime.\n');
fprintf('Welfare and fiscal-expenditure comparisons (Tables 6-9) are computed\n');
fprintf('separately and are not part of this script.\n');

%==========================================================================
%% POLICY PARAMETER OVERRIDES
%==========================================================================

function param = apply_ub_plus15(param)
    % Section 6, Table policy1: 15% higher unemployment-benefit baseline
    % rate b1. In this model b(x) = min(bmin + bscale.*exp(x), bmax); b1
    % is the intercept, param.bmin. Other benefit parameters unchanged.
    param.bmin = param.bmin * 1.15;
end

function param = apply_credit_plus20(param)
    % Section 6, Table policy2: 20% higher income caps for the child tax
    % credit, households with a dependent child only (k=1): the
    % phase-in end (E1_n1) and phase-out start (E2_n1) thresholds.
    param.E1_n1 = param.E1_n1 * 1.20;
    param.E2_n1 = param.E2_n1 * 1.20;
end
