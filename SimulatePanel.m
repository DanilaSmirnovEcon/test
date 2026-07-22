function [agentPanel, agentPanel_full] = SimulatePanel(solutionSP, valueFunc, gridSf, gridSm, gridP, param, simParam)
% SimulatePanel  � agent-outer, time-inner parallel simulation
%   � Three phases with custom step lengths:
%       1) coarse burn-in, 2) fine burn-in, 3) actual simulation
%   � Hazards/drifts scale with step length h (in MONTHS)
%   � Output is trimmed to phase (3) and tagged with months_per_step
%   � Transition rates in simResults are reported in MONTHLY units
%
% Inputs:
%   solutionSP, valueFunc, gridSf, gridSm, gridP, param  � model objects
%   simParam: struct with fields (defaults below)
%       N, seed, nWorkers
%       coarse_periods, coarse_step_months
%       burn_periods,   burn_step_months
%       sim_periods,    sim_step_months
%
% Output:
%   agentPanel: panel for the FINAL phase only (actual simulation)
%     - includes agentPanel.months_per_step = simParam.sim_step_months
%   simResults: time series + steady-state summaries (monthly rates)

% ---------------- Defaults ----------------
if ~exist('simParam','var'), simParam = struct(); end
if ~isfield(simParam,'N'),                simParam.N = 10000; end
if ~isfield(simParam,'seed'),             simParam.seed = 12345; end
if ~isfield(simParam,'nWorkers'),         simParam.nWorkers = []; end

if ~isfield(simParam,'coarse_periods'),   simParam.coarse_periods = 49; end
if ~isfield(simParam,'coarse_step_months'), simParam.coarse_step_months = 12; end

if ~isfield(simParam,'burn_periods'),     simParam.burn_periods = 24; end
if ~isfield(simParam,'burn_step_months'), simParam.burn_step_months = 0.25; end  % 1 week

if ~isfield(simParam,'sim_periods'),      simParam.sim_periods = 120; end        % 10 years monthly
if ~isfield(simParam,'sim_step_months'),  simParam.sim_step_months = 0.25; end   % 1 week

if ~isfield(param,'fkidwpenalty'), param.fkidwpenalty = 0; end
param.fpem = 0;%max(0, min(1, param.fkidwpenalty));     % clamp to [0,1]

if ~isfield(param,'entry_age_years'), param.entry_age_years = 25; end
age0_m = param.entry_age_years * 12;   % months

% ---- Acceptance softening (defaults) ----
if ~isfield(param,'accept_model'), param.accept_model = 'probit'; end  % 'probit' or 'logit'
if ~isfield(param,'accept_sigma'), param.accept_sigma = 0.00;     end  % noise scale
if ~isfield(param,'accept_tau'),   param.accept_tau   = 0.00;     end  % surplus tolerance (>=0 stricter, <0 laxer)


% ---- Build step schedule (months per step for each phase) ----
stepMonths = [ ...
    repmat(simParam.coarse_step_months, simParam.coarse_periods, 1); ...
    repmat(simParam.burn_step_months,   simParam.burn_periods,   1); ...
    repmat(simParam.sim_step_months,    simParam.sim_periods,    1) ...
];
dt_vec  = stepMonths(:)';          % 1 x nSteps, in months
nSteps  = numel(stepMonths);
T_total = nSteps + 1;              % including initial column t=1
N       = simParam.N;

% ---- CDFs for wage offers (inverse-CDF draws) ----
if isfield(param, 'f_female'), f_female = param.f_female(:); else, f_female = param.f_long(:); end
if isfield(param, 'f_male'),   f_male   = param.f_male(:);   else, f_male   = param.f_long(:); end

wgrid_f = gridSf.w.w(:);   wgrid_m = gridSm.w.w(:);
cdf_f = cumsum(max(f_female,0)); cdf_f = cdf_f / max(sum(f_female),eps);
cdf_m = cumsum(max(f_male,  0)); cdf_m = cdf_m / max(sum(f_male),  eps);
edges_f = [0; cdf_f(:)];
edges_m = [0; cdf_m(:)];

% ---- Partner samplers (outside marriage market) ----
shadowF = solutionSP.solutionSf_s;
shadowM = solutionSP.solutionSm_s;
samplerPartnerF = buildShadowPartnerSampler(gridSf.a.a(:), gridSf.w.w(:), shadowF);
samplerPartnerM = buildShadowPartnerSampler(gridSm.a.a(:), gridSm.w.w(:), shadowM);

% ---- Pre-allocate (plain arrays; parfor-safe) ----
Z  = NaN(N, T_total);  O = zeros(N, T_total);  NZ = NaN(N, T_total);

ability          = Z;
wage             = Z;
employed         = O;
olf              = O;
married          = O;
has_kids         = O;

partner_ability  = NZ;
partner_wage     = NZ;
partner_employed = NZ;
partner_olf      = NZ;

value_now        = Z;
job_offer        = O;
marriage_offer   = O;
divorce          = O;             %#ok<NASGU>
birth            = O;

age              = O;             % months since (re)birth
age_chron_full   = O;
newborn          = O;             % =1 at (re)birth

% ---- Fixed agent attributes ----
rng(simParam.seed,'twister');
female = rand(N,1) < 0.5;

% Initial ability draw by gender
a0 = zeros(N,1);
iF = female; iM = ~female;
a0(iF) = clip(randn(sum(iF),1).*sqrt(param.pf.aergvar) + param.pf.aergmean, param.a.min, param.a.max);
a0(iM) = clip(randn(sum(iM),1).*sqrt(param.pm.aergvar) + param.pm.aergmean, param.a.min, param.a.max);

% t=1 initialization
ability(:,1) = a0;
employed(:,1) = 0;
olf(:,1)      = 0;   % <<< CHANGED: newborns start OLF, not unemployed
married(:,1)  = 0;
has_kids(:,1) = 0;
wage(:,1)     = 0;
age(:,1)        = 0;        % tenure since (re)birth
age_chron_full(:,1) = age0_m; % chronological age at entry (e.g., 25y)
newborn(:,1)    = 1;

% Initial value function (OLF, no kids)
for i = 1:N
    if female(i)
        value_now(i,1) = valueFunc.single_female.unemployed_nokid(ability(i,1));  % <<< CHANGED
    else
        value_now(i,1) = valueFunc.single_male.unemployed_nokid(ability(i,1));    % <<< CHANGED
    end
end

% ---- Pre-generate RNG for reproducibility ----
U_sep     = rand(N,T_total);
U_offer   = rand(N,T_total);
U_olf     = rand(N,T_total);
U_birth_s = rand(N,T_total);
U_birth_m = rand(N,T_total);
U_meet    = rand(N,T_total);
U_death   = rand(N,T_total);

Up_sep    = rand(N,T_total);
Up_offer  = rand(N,T_total);
Up_olf    = rand(N,T_total);

Z_mu_self    = randn(N,T_total);
Z_mu_u_self  = randn(N,T_total);
Z_mu_partner = randn(N,T_total);
Z_mu_u_part  = randn(N,T_total);

U_w_self     = rand(N,T_total);
U_w_partner  = rand(N,T_total);

% ---- Acceptance shocks (for soft acceptance) ----
U_acc_self  = rand(N,T_total);    % for logit
U_acc_part  = rand(N,T_total);
Z_acc_self  = randn(N,T_total);   % for probit
Z_acc_part  = randn(N,T_total);

% ---- Parallel pool ----
if isempty(gcp('nocreate'))
    if isempty(simParam.nWorkers)
        parpool;
    else
        parpool(simParam.nWorkers);
    end
end

% Kept-window bounds (compute BEFORE parfor so inside loop can see them)
B  = simParam.coarse_periods + simParam.burn_periods;
K  = simParam.sim_periods;
T0 = B + 1;
T1 = B + K;

% Per-agent counters (parfor-safe: each i writes its own row)
meet_cnt         = zeros(N,1);
accept_cnt       = zeros(N,1);
self_reject_cnt  = zeros(N,1);
part_reject_cnt  = zeros(N,1);

% (optional) partner characteristics at meets
meet_partnerKids_cnt = zeros(N,1);
meet_partnerW_sum    = zeros(N,1);

% ======================== Core simulation ================================
parfor i = 1:N
    % local paths (row vectors)
    a  = ability(i,:); e = employed(i,:); o = olf(i,:); m = married(i,:);
    k  = has_kids(i,:); w = wage(i,:);
    pa = partner_ability(i,:);  pw = partner_wage(i,:);
    pe = partner_employed(i,:); po = partner_olf(i,:);
    V  = value_now(i,:);  jo = job_offer(i,:); mo = marriage_offer(i,:);
    ag = age(i,:);        nb = newborn(i,:); agc = age_chron_full(i,:); 
    br = birth(i,:);   % <<< ADD THIS (local buffer for birth events)

    isF = female(i);

    for t = 2:T_total
        h = dt_vec(t-1);  % months in this step

        % shorthand uniforms/shocks
        u_sep   = U_sep(i,t);      u_off  = U_offer(i,t);     u_olf = U_olf(i,t);
        u_b_s   = U_birth_s(i,t);  u_b_m  = U_birth_m(i,t);   u_me  = U_meet(i,t);
        u_die   = U_death(i,t);
        up_sep  = Up_sep(i,t);     up_off = Up_offer(i,t);    up_olf = Up_olf(i,t);
        z_s     = Z_mu_self(i,t);  z_su   = Z_mu_u_self(i,t);
        z_p     = Z_mu_partner(i,t); z_pu = Z_mu_u_part(i,t);
        u_w_s   = U_w_self(i,t);   u_w_p  = U_w_partner(i,t);

        u_acc_s = U_acc_self(i,t);  u_acc_p = U_acc_part(i,t);
        z_acc_s = Z_acc_self(i,t);  z_acc_p = Z_acc_part(i,t);

        % advance age
        ag(t)  = ag(t-1)  + h;    % tenure
        agc(t) = agc(t-1) + h;    % chronological
        nb(t)  = 0;
        br(t)  = 0;

        if m(t-1)==0
            % -------- SINGLE evolution --------
        [a(t), e(t), o(t), w(t), k(t), V(t), br(t), jo(t)] = ...
            step_single_h(a(t-1), e(t-1), o(t-1), k(t-1), isF, w(t-1), ...
                          u_sep, u_off, u_olf, u_b_s, z_s, z_su, ...
                          valueFunc, param, h, wgrid_f, edges_f, wgrid_m, edges_m, u_w_s);


            % Marriage market meeting
            meet = u_me < 1 - exp(-param.m*h); % (u_me > exp(-param.m * h));
            in_keep = (t >= T0) && (t <= T1);
            if meet
                if in_keep
                meet_cnt(i) = meet_cnt(i) + 1;
                end
                if isF
                    pm = samplerPartnerM(1);
                    % track partner features at the meet
                    p_kid_partner = param.partner_kid_prob_male;
                    pm.has_kids = pm.has_kids || (rand < p_kid_partner);
                    hk = k(t) || pm.has_kids;
                    if in_keep
                        meet_partnerKids_cnt(i) = meet_partnerKids_cnt(i) + double(pm.has_kids);
                        meet_partnerW_sum(i)    = meet_partnerW_sum(i)    + max(0, pm.wage);
                    end

                    if param.always_marry
                    % Force marriage: accept partner without value comparisons
                    % --- compute these before overwriting k(t) ---
                    k_pre = k(t);                 % after step_single_h (includes births this step)                    
                    m(t)=1; mo(t)=1; k(t)=hk;
                    pa(t)=pm.ability; pe(t)=pm.employed; po(t)=pm.olf; pw(t)=pm.wage;

                    % Adoption penalty: only if she had no kids before, partner had kids,
                    % and she is currently employed.
                    if isF && e(t)==1 && (k_pre==0) && pm.has_kids
                        w(t) = w(t) * (1 - param.fpem);
                    end
                    else
                    V_self_s = getSingleValue(a(t), e(t), w(t), o(t), k(t), true,  valueFunc);
                    V_part_s = getSingleValue(pm.ability, pm.employed, pm.wage, pm.olf, pm.has_kids, false, valueFunc);
                    V_self_m = getMarriedValue(a(t), e(t), w(t), pm.ability, pm.employed, pm.wage, hk, true,  valueFunc, o(t), pm.olf);
                    V_part_m = getMarriedValue(pm.ability, pm.employed, pm.wage, a(t), e(t), w(t), hk, false, valueFunc, pm.olf, o(t));
                    %self_ok = (V_self_m > V_self_s);
                    %part_ok = (V_part_m > V_part_s);

                    % ----- compute margins -----
                    d_self = (V_self_m - V_self_s) - param.accept_tau;
                    d_part = (V_part_m - V_part_s) - param.accept_tau;
                    
                    % ----- softened acceptance (choose ONE style) -----
                    % binding-side probit (robust)
                    d_bind = min(d_self, d_part);
                    if strcmpi(param.accept_model,'probit')
                        accept = (d_bind + param.accept_sigma * z_acc_s) > 0;
                    else % 'logit'
                        p = 1 ./ (1 + exp(- d_bind / max(1e-8, param.accept_sigma)));
                        accept = (u_acc_s < p);
                    end
                    
                    % ----- update counters USING FINAL DECISION -----
                    if in_keep
                        meet_cnt(i) = meet_cnt(i) + 1;
                        if accept
                            accept_cnt(i) = accept_cnt(i) + 1;
                        else
                            if d_self <= d_part, self_reject_cnt(i) = self_reject_cnt(i) + 1;
                            else,                 part_reject_cnt(i) = part_reject_cnt(i) + 1;
                            end
                        end
                    end

                    if accept %V_self_m > V_self_s && V_part_m > V_part_s
                        % capture pre-marriage kids flag (after step_single_h; births this step already counted)
                        k_pre = k(t);
                        m(t)=1; mo(t)=1; k(t)=hk;
                        pa(t)=pm.ability; pe(t)=pm.employed; po(t)=pm.olf; pw(t)=pm.wage;
                        % Adoption penalty ONLY: she had no kids, partner had kids, and she’s employed
                        if e(t)==1 && (k_pre==0) && pm.has_kids
                            w(t) = w(t) * (1 - param.fpem);
                        end
                    end

                    end
                else
                    pf = samplerPartnerF(1);
                    p_kid_partner = param.partner_kid_prob_female;
                    pf.has_kids = pf.has_kids || (rand < p_kid_partner);
                    hk = k(t) || pf.has_kids;
                    if in_keep
                        meet_partnerKids_cnt(i) = meet_partnerKids_cnt(i) + double(pf.has_kids);
                        meet_partnerW_sum(i)    = meet_partnerW_sum(i)    + max(0, pf.wage);
                    end
                    % if in_keep
                    %     meet_partnerKids_cnt(i) = meet_partnerKids_cnt(i) + double(pf.has_kids);
                    %     meet_partnerW_sum(i)    = meet_partnerW_sum(i)    + max(0, pf.wage);
                    % end
                    % p_kid_partner = param.partner_kid_prob_female; % calibrate this
                    % pf.has_kids = pf.has_kids || (rand < p_kid_partner);
                    % hk = k(t) || pf.has_kids;
                    if param.always_marry
                    % Force marriage: accept partner without value comparisons
                    % --- compute these before overwriting k(t) ---
                    k_pre = k(t);                 % after step_single_h (includes births this step)
                    m(t)=1; mo(t)=1; k(t)=hk;
                    pa(t)=pf.ability; pe(t)=pf.employed; po(t)=pf.olf; pw(t)=pf.wage;
                    % Adoption penalty: she had no kids, man had kids (k_pre==1), and she’s employed
                    if pe(t)==1 && (pf.has_kids==0) && (k_pre==1)
                        pw(t) = pw(t) * (1 - param.fpem);
                    end
                    else
                    V_self_s = getSingleValue(a(t), e(t), w(t), o(t), k(t), false, valueFunc);
                    V_part_s = getSingleValue(pf.ability, pf.employed, pf.wage, pf.olf, pf.has_kids, true,  valueFunc);
                    V_self_m = getMarriedValue(a(t), e(t), w(t), pf.ability, pf.employed, pf.wage, hk, false, valueFunc, o(t), pf.olf);
                    V_part_m = getMarriedValue(pf.ability, pf.employed, pf.wage, a(t), e(t), w(t), hk, true,  valueFunc, pf.olf, o(t));
                    %self_ok = (V_self_m > V_self_s);
                    %part_ok = (V_part_m > V_part_s);

                    % ----- compute margins -----
                    d_self = (V_self_m - V_self_s) - param.accept_tau;
                    d_part = (V_part_m - V_part_s) - param.accept_tau;
                    
                    % ----- softened acceptance (choose ONE style) -----
                    % binding-side probit (robust)
                    d_bind = min(d_self, d_part);
                    if strcmpi(param.accept_model,'probit')
                        accept = (d_bind + param.accept_sigma * z_acc_s) > 0;
                    else % 'logit'
                        p = 1 ./ (1 + exp(- d_bind / max(1e-8, param.accept_sigma)));
                        accept = (u_acc_s < p);
                    end
                    
                    % ----- update counters USING FINAL DECISION -----
                    if in_keep
                        meet_cnt(i) = meet_cnt(i) + 1;
                        if accept
                            accept_cnt(i) = accept_cnt(i) + 1;
                        else
                            if d_self <= d_part, self_reject_cnt(i) = self_reject_cnt(i) + 1;
                            else,                 part_reject_cnt(i) = part_reject_cnt(i) + 1;
                            end
                        end
                    end

                    if accept %V_self_m > V_self_s && V_part_m > V_part_s
                        k_pre = k(t);
                        m(t)=1; mo(t)=1; k(t)=hk;
                        pa(t)=pf.ability; pe(t)=pf.employed; po(t)=pf.olf; pw(t)=pf.wage;
                        % Adoption penalty for female partner ONLY:
                        % she had no kids, man had kids (k_pre==1), and she’s employed
                        if pe(t)==1 && (pf.has_kids==0) && (k_pre==1)
                            pw(t) = pw(t) * (1 - param.fpem);
                        end
                    end

                    end
                end
            end

        else
            % -------- MARRIED evolution (self + partner) --------
        [a(t), e(t), o(t), w(t), k(t), V(t), br(t), ...
         pa(t), pe(t), po(t), pw(t), jo(t)] = ...
            step_married_h(a(t-1), e(t-1), o(t-1), k(t-1), isF, w(t-1), ...
                           pa(t-1), pe(t-1), po(t-1), pw(t-1), ...
                           u_sep, u_off, u_olf, u_b_m, z_s, z_su, z_p, z_pu, ...
                           up_sep, up_off, up_olf, ...
                           valueFunc, param, h, wgrid_f, edges_f, wgrid_m, edges_m, u_w_s, u_w_p);
            m(t) = 1;  % stays married unless attrition resets
        end

        % ---- Attrition (death) ----
        if (u_die > exp(-param.zeta * h))
            % reset to newborn
            ag(t) = 0; nb(t) = 1; agc(t) = age0_m;
            m(t) = 0; e(t) = 0; o(t) = 0;  % <<< CHANGED: newborns start unemployed
            w(t) = 0; k(t) = 0;
            pa(t)=NaN; pe(t)=NaN; po(t)=NaN; pw(t)=NaN; mo(t)=0; jo(t)=0;
            if isF
                a(t) = clip(randn*sqrt(param.pf.aergvar) + param.pf.aergmean, param.a.min, param.a.max);
                value_now_i = valueFunc.single_female.unemployed_nokid(a(t));  % <<< CHANGED
            else
                a(t) = clip(randn*sqrt(param.pm.aergvar) + param.pm.aergmean, param.a.min, param.a.max);
                value_now_i = valueFunc.single_male.unemployed_nokid(a(t));    % <<< CHANGED
            end
            V(t) = value_now_i;
        end

    end

    % write back (parfor-friendly)
    ability(i,:)          = a;     employed(i,:)         = e;
    olf(i,:)              = o;     married(i,:)          = m;
    has_kids(i,:)         = k;     wage(i,:)             = w;
    partner_ability(i,:)  = pa;    partner_employed(i,:) = pe;
    partner_olf(i,:)      = po;    partner_wage(i,:)     = pw;
    value_now(i,:)        = V;     job_offer(i,:)        = jo;
    marriage_offer(i,:)   = mo;    age(i,:)              = ag;
    newborn(i,:)          = nb;    birth(i,:)            = br; 
    age_chron_full(i,:)   = agc;
end

% ======================== DIAGNOSTICS ON MEETINGS ================================

isW = female; isM = ~female;

meetW = sum(meet_cnt(isW));  accW = sum(accept_cnt(isW));
rSelfW = sum(self_reject_cnt(isW)); rPartW = sum(part_reject_cnt(isW));

meetM = sum(meet_cnt(isM));  accM = sum(accept_cnt(isM));
rSelfM = sum(self_reject_cnt(isM)); rPartM = sum(part_reject_cnt(isM));

fprintf('MEETS:   W=%d,  M=%d\n', meetW, meetM);
fprintf('ACCEPTS: W=%d (rate given meet=%.3f),  M=%d (rate=%.3f)\n', ...
        accW, accW/max(1,meetW), accM, accM/max(1,meetM));
fprintf('REJECT (given meet): Women self %.3f, partner %.3f | Men self %.3f, partner %.3f\n', ...
        rSelfW/max(1,meetW), rPartW/max(1,meetW), rSelfM/max(1,meetM), rPartM/max(1,meetM));

% Optional partner mix at meets
avgKidsAtMeet_W = sum(meet_partnerKids_cnt(isW)) / max(1,meetW);
avgKidsAtMeet_M = sum(meet_partnerKids_cnt(isM)) / max(1,meetM);
fprintf('Partner kids at meet (avg prob): W sees %.3f, M sees %.3f\n', avgKidsAtMeet_W, avgKidsAtMeet_M);

% ======================== PRE_TRIM PANEL ================================
agentPanel_full = struct();
agentPanel_full.id      = (1:N)';
agentPanel_full.female  = female;

agentPanel_full.married = married;          % N x T_total
agentPanel_full.newborn = newborn;          % N x T_total
agentPanel_full.age_chron = age_chron_full; % N x T_total (chronological age, months)

% Useful to know varying step sizes over the full path:
agentPanel_full.step_months = dt_vec;       % 1 x (T_total-1) months per step

% ========= TRIM TO KEPT WINDOW (after parfor, before packaging) =========
B  = simParam.coarse_periods + simParam.burn_periods;  % steps to burn
K  = simParam.sim_periods;                             % steps to keep
T0 = B + 1;                         % first kept column (1-based)
T1 = B + K;                         % last kept column (inclusive)

if T1 > T_total
    error('Trim window exceeds simulated horizon: last=%d, T_total=%d', T1, T_total);
end
cols = T0:T1;  % exactly K columns

% Slice raw matrices to the kept window
ability          = ability(:, cols);
employed         = employed(:, cols);
olf              = olf(:, cols);
married          = married(:, cols);
has_kids         = has_kids(:, cols);
wage             = wage(:, cols);

partner_ability  = partner_ability(:, cols);
partner_employed = partner_employed(:, cols);
partner_olf      = partner_olf(:, cols);
partner_wage     = partner_wage(:, cols);

value_now        = value_now(:, cols);
job_offer        = job_offer(:, cols);
marriage_offer   = marriage_offer(:, cols);
birth            = birth(:, cols);
age              = age(:, cols);
age_chron_full   = age_chron_full(:, cols);
newborn          = newborn(:, cols);

% If you don't track divorces explicitly, keep zeros with correct size:
divorce          = zeros(N, numel(cols));

% ====================== ASSEMBLE agentPanel ======================
agentPanel = struct();
agentPanel.id               = (1:N)';
agentPanel.female           = female;

agentPanel.ability          = ability;
agentPanel.employed         = employed;
agentPanel.olf              = olf;
agentPanel.married          = married;
agentPanel.has_kids         = has_kids;
agentPanel.wage             = wage;

agentPanel.partner_ability  = partner_ability;
agentPanel.partner_employed = partner_employed;
agentPanel.partner_olf      = partner_olf;
agentPanel.partner_wage     = partner_wage;

agentPanel.value            = value_now;
agentPanel.job_offer        = job_offer;
agentPanel.marriage_offer   = marriage_offer;
agentPanel.birth            = birth;
agentPanel.divorce          = divorce;
agentPanel.age              = age;
agentPanel.age_chron        = age_chron_full;
agentPanel.newborn          = newborn;

% Expose calendar unit for the kept window (used by compute_moments)
agentPanel.months_per_step  = simParam.sim_step_months;

% ====================== STATS ======================
%simResults = calculatePanelStatistics(agentPanel, param);

end

% ========================================================================
% ============================== HELPERS =================================
% ========================================================================

function x = clip(x, xmin, xmax)
x = min(max(x, xmin), xmax);
end

function [a_t, emp_t, olf_t, w_t, kids_t, V_now, birth_t, job_offer_t] = ...
    step_single_h(ability, employed, olf, has_kids, is_female, wage_tm1, ...
                  u_sep, u_offer, u_olf, u_birth, z_mu, z_mu_u, ...
                  valueFunc, param, h, wgrid_f, edges_f, wgrid_m, edges_m, u_w)

% OU-like drift scaled by h; shocks  sqrt(h)
alpha  = @(mu,hh) 1 - exp(-mu*hh);    % always in [0,1), no complex values
sigmah = @(sig,hh) sqrt(hh).*sig;

emp_t = employed; olf_t = olf; kids_t = has_kids; w_t = wage_tm1; birth_t = 0; job_offer_t = 0;

if is_female
    vF = valueFunc.single_female;
    lam_u = param.pf.lambda_u; lam_o = param.pf.lambda_o; delta = param.pf.delta; lam_oE = param.pf.lambda_oE; %%
    a_sig = param.pf.a_sigma;  mu    = param.pf.mu;        mu_u  = param.pf.mu_u;  amean = param.pf.amean;
    draw_idx = @(u) max(1, discretize(min(max(u, realmin), 1 - eps), edges_f));
    draw_w   = @(u) wgrid_f(draw_idx(min(max(u, realmin), 1 - eps)));
    % draw_idx_f  = @(u) max(1, discretize(min(max(u, realmin), 1 - eps), edges_f));
    % grid_now_f  = wgrid_f .* (1 - param.fpem*double(has_kids~=0));
    % draw_w      = @(u) grid_now_f(draw_idx_f(min(max(u, realmin), 1 - eps)));
else
    vF = valueFunc.single_male;
    lam_u = param.pm.lambda_u; lam_o = param.pm.lambda_o; delta = param.pm.delta; lam_oE = param.pm.lambda_oE;%%
    a_sig = param.pm.a_sigma;  mu    = param.pm.mu;        mu_u  = param.pm.mu_u;  amean = param.pm.amean;
    draw_idx = @(u) max(1, discretize(min(max(u, realmin), 1 - eps), edges_m));
    draw_w   = @(u) wgrid_m(draw_idx(min(max(u, realmin), 1 - eps)));
end
p_birth = is_female * param.pf.pcs + (~is_female) * param.pm.pcs;

sep_hit   = (u_sep   < 1 - exp(-delta * h));
offer_hit = (u_offer < 1 - exp(-lam_u * h));
olf_hit   = (u_olf   < 1 - exp(-lam_o * h));
birth_hit = (u_birth < 1 - exp(-p_birth * h));

if employed
    if sep_hit
        emp_t=0; w_t=0;
        a_t = ability + alpha(mu_u,h)*(param.a.min-ability) + sigmah(a_sig,h)*z_mu_u;
    else
        a_emp = ability + alpha(mu,h)*(amean-ability) + sigmah(a_sig,h)*z_mu; % if stay E
        if olf_hit
            V_e = has_kids*vF.employed_kid(a_emp,w_t) + (~has_kids)*vF.employed_nokid(a_emp,w_t);
            V_o = has_kids*vF.olf_kid(a_emp)          + (~has_kids)*vF.olf_nokid(a_emp);
            if V_o > V_e
                emp_t=0; olf_t=1; w_t=0;
                a_t = ability + alpha(mu_u,h)*(param.a.min-ability) + sigmah(a_sig,h)*z_mu_u;
            else
                a_t = a_emp;
            end
        else
            a_t = a_emp;
        end
    end

elseif ~olf
    a_rej = ability + alpha(mu_u,h)*(param.a.min-ability) + sigmah(a_sig,h)*z_mu_u; % stay U
    if offer_hit
        job_offer_t = 1;
        w_offer = draw_w(u_w);
        a_acc   = ability + alpha(mu,h)*(amean-ability) + sigmah(a_sig,h)*z_mu;      % take job

        V_acc = has_kids*vF.employed_kid(a_acc,w_offer) + (~has_kids)*vF.employed_nokid(a_acc,w_offer);
        V_rej = has_kids*vF.unemployed_kid(a_rej)       + (~has_kids)*vF.unemployed_nokid(a_rej);
        if V_acc > V_rej, emp_t=1; w_t=w_offer; a_t=a_acc; else, a_t=a_rej; end
    else
        a_t=a_rej;
    end
    if (~emp_t) && olf_hit
        V_u = has_kids*vF.unemployed_kid(a_t) + (~has_kids)*vF.unemployed_nokid(a_t);
        V_o = has_kids*vF.olf_kid(a_t)        + (~has_kids)*vF.olf_nokid(a_t);
        if V_o > V_u, olf_t=1; w_t=0; end
    end
else
    % ---------- OLF state (no OLF -> U; OLF -> E offers at rate lambda_oE) ----------
    a_olf = ability + alpha(mu_u,h)*(param.a.min-ability) + sigmah(a_sig,h)*z_mu_u;
    w_t   = 0;

    offer_o_hit = (u_offer < 1- exp(-lam_oE * h));  % OLF→E offer arrival
    if offer_o_hit
        job_offer_t = 1;
        w_offer = draw_w(u_w);
        a_acc   = ability + alpha(mu,h)*(amean-ability) + sigmah(a_sig,h)*z_mu;

        V_acc  = has_kids*vF.employed_kid(a_acc,w_offer) + (~has_kids)*vF.employed_nokid(a_acc,w_offer);
        V_stay = has_kids*vF.olf_kid(a_olf)              + (~has_kids)*vF.olf_nokid(a_olf);

        if V_acc > V_stay
            emp_t = 1; olf_t = 0; w_t = w_offer; a_t = a_acc;
        else
            a_t = a_olf;
        end
    else
        a_t = a_olf;
    end
end

a_t = clip(a_t, param.a.min, param.a.max);

% value + fertility
if ~has_kids && birth_hit
    kids_t  = 1;
    birth_t = 1;
else
    kids_t = has_kids;
end

if is_female && emp_t && (birth_t==1)
    w_t = w_t * (1 - param.fpem);
end

% --- Clamp E/OLF invariants & sanitize wage ------------------------------
emp_t = double(emp_t~=0);
olf_t = double(olf_t~=0);

olf_t = (~emp_t) & (olf_t ~= 0);           % OLF only if not employed
w_t   = w_t .* double(emp_t);              % zero-out wage when not employed

% (optional) keep offer flag tidy
job_offer_t = double(job_offer_t~=0);

% Ensure V_now matches the resulting state (after potential birth)
if emp_t
    V_now = kids_t*vF.employed_kid(a_t,w_t) + (~kids_t)*vF.employed_nokid(a_t,w_t);
elseif ~olf_t
    V_now = kids_t*vF.unemployed_kid(a_t)   + (~kids_t)*vF.unemployed_nokid(a_t);
else
    V_now = kids_t*vF.olf_kid(a_t)          + (~kids_t)*vF.olf_nokid(a_t);
end
% if emp_t
%     V_now = has_kids*vF.employed_kid(a_t,w_t) + (~has_kids)*vF.employed_nokid(a_t,w_t);
% elseif ~olf_t
%     V_now = has_kids*vF.unemployed_kid(a_t)   + (~has_kids)*vF.unemployed_nokid(a_t);
% else
%     V_now = has_kids*vF.olf_kid(a_t)          + (~has_kids)*vF.olf_nokid(a_t);
% end
% 
% if ~has_kids && birth_hit
%     if emp_t
%         Vk = vF.employed_kid(a_t,w_t); Vn = vF.employed_nokid(a_t,w_t);
%     elseif ~olf_t
%         Vk = vF.unemployed_kid(a_t);   Vn = vF.unemployed_nokid(a_t);
%     else
%         Vk = vF.olf_kid(a_t);          Vn = vF.olf_nokid(a_t);
%     end
%     if Vk > Vn, kids_t=1; birth_t=1; V_now=Vk; else, kids_t=has_kids; end
% else
%     kids_t = has_kids;
% end
end

function [a_t, emp_t, olf_t, w_t, kids_t, V_now, birth_t, ...
          p_a_t, p_e_t, p_o_t, p_w_t, job_offer_t] = ...
    step_married_h( ability, employed, olf, has_kids, is_female, wage_tm1, ...
                    p_ability, p_employed, p_olf, p_wage, ...
                    u_sep, u_offer, u_olf, u_birth, z_mu, z_mu_u, z_pm, z_pm_u, ...
                    up_sep, up_offer, up_olf, ...
                    valueFunc, param, h, wgrid_f, edges_f, wgrid_m, edges_m, u_w_self, u_w_part)

alpha  = @(mu,hh) 1 - exp(-mu*hh);    % always in [0,1), no complex values
sigmah = @(sig,hh) sqrt(hh).*sig;
clipa  = @(x) clip(x, param.a.min, param.a.max);

emp_t = employed; olf_t = olf; w_t = wage_tm1; kids_t = has_kids; birth_t = 0; job_offer_t = 0;
p_e_t = p_employed; p_o_t = p_olf; p_w_t = p_wage; p_a_t = p_ability;
if isnan(p_e_t), p_e_t=0; end; if isnan(p_o_t), p_o_t=0; end
if isnan(p_w_t), p_w_t=0; end; if isnan(p_a_t), p_a_t = (is_female)*param.pm.amean + (~is_female)*param.pf.amean; end

if is_female
    lam_u = param.pf.lambda_u; lam_o = param.pf.lambda_o; delta = param.pf.delta; lam_oE = param.pf.lambda_oE;%%
    a_sig = param.pf.a_sigma;  mu    = param.pf.mu;        mu_u = param.pf.mu_u;  amean = param.pf.amean;
    draw_idx_self = @(u) max(1, discretize(min(max(u, realmin), 1 - eps), edges_f));
    draw_w_self   = @(u) wgrid_f(draw_idx_self(min(max(u, realmin), 1 - eps)));
    % draw_idx_self = @(u) max(1, discretize(min(max(u, realmin), 1 - eps), edges_f));
    % grid_self_now = wgrid_f .* (1 - param.fpem*double(kids_t~=0));  % <- kids_t is the couple’s current kids flag
    % draw_w_self   = @(u) grid_self_now(draw_idx_self(min(max(u, realmin), 1 - eps)));

    p_lam_u = param.pm.lambda_u; p_lam_o = param.pm.lambda_o; p_delta = param.pm.delta; p_lam_oE = param.pm.lambda_oE;%%
    p_a_sig = param.pm.a_sigma;  p_mu    = param.pm.mu;       p_mu_u  = param.pm.mu_u; p_amean = param.pm.amean;
    draw_idx_p   = @(u) max(1, discretize(min(max(u, realmin), 1 - eps), edges_m));
    draw_w_part  = @(u) wgrid_m(draw_idx_p(min(max(u, realmin), 1 - eps)));
    % draw_idx_p   = @(u) max(1, discretize(min(max(u, realmin), 1 - eps), edges_m));
    % grid_part_now = wgrid_m;
    % draw_w_part  = @(u) grid_part_now(draw_idx_p(min(max(u, realmin), 1 - eps)));
else
    lam_u = param.pm.lambda_u; lam_o = param.pm.lambda_o; delta = param.pm.delta; lam_oE = param.pm.lambda_oE;%%
    a_sig = param.pm.a_sigma;  mu    = param.pm.mu;        mu_u = param.pm.mu_u;  amean = param.pm.amean;
    draw_idx_self = @(u) max(1, discretize(min(max(u, realmin), 1 - eps), edges_m));
    draw_w_self   = @(u) wgrid_m(draw_idx_self(min(max(u, realmin), 1 - eps)));
    % draw_idx_self = @(u) max(1, discretize(min(max(u, realmin), 1 - eps), edges_m));
    % grid_self_now = wgrid_m;
    % draw_w_self   = @(u) grid_self_now(draw_idx_self(min(max(u, realmin), 1 - eps)));

    p_lam_u = param.pf.lambda_u; p_lam_o = param.pf.lambda_o; p_delta = param.pf.delta; p_lam_oE = param.pf.lambda_oE;%%
    p_a_sig = param.pf.a_sigma;  p_mu    = param.pf.mu;       p_mu_u  = param.pf.mu_u; p_amean = param.pf.amean;
    draw_idx_p   = @(u) max(1, discretize(min(max(u, realmin), 1 - eps), edges_f));
    draw_w_part  = @(u) wgrid_f(draw_idx_p(min(max(u, realmin), 1 - eps)));
    % draw_idx_p    = @(u) max(1, discretize(min(max(u, realmin), 1 - eps), edges_f));
    % grid_part_now = wgrid_f .* (1 - param.fpem*double(kids_t~=0));
    % draw_w_part   = @(u) grid_part_now(draw_idx_p(min(max(u, realmin), 1 - eps)));
end
p_birth = param.pcp;

sep_hit   = (u_sep   < 1 - exp(-delta * h));
offer_hit = (u_offer < 1 - exp(-lam_u * h));
olf_hit   = (u_olf   < 1 - exp(-lam_o * h));
birth_hit = (u_birth < 1 - exp(-p_birth * h));

p_sep_hit   = (up_sep   < 1 - exp(-p_delta * h));
p_offer_hit = (up_offer < 1 - exp(-p_lam_u * h));
p_olf_hit   = (up_olf   < 1 - exp(-p_lam_o * h));

% --- Partner dynamics
if p_e_t
    if p_sep_hit
        p_e_t=0; p_w_t=0;
        p_a_t = p_ability + alpha(p_mu_u,h)*(param.a.min - p_ability) + sigmah(p_a_sig,h)*z_pm_u;
    else
        p_a_emp = p_ability + alpha(p_mu,h)*(p_amean - p_ability) + sigmah(p_a_sig,h)*z_pm;
        if p_olf_hit
            V_e = getMarriedValue(ability, employed, w_t, p_a_emp, 1, p_w_t, kids_t, is_female, valueFunc, olf_t, 0 );
            V_o = getMarriedValue(ability, employed, w_t, p_a_emp, 0, 0,     kids_t, is_female, valueFunc, olf_t, 1 );
            if V_o > V_e
                p_e_t=0; p_o_t=1; p_w_t=0;
                p_a_t = p_ability + alpha(p_mu_u,h)*(param.a.min - p_ability) + sigmah(p_a_sig,h)*z_pm_u;
            else
                p_a_t = p_a_emp;
            end
        else
            p_a_t = p_a_emp;
        end
    end
elseif ~p_o_t
    p_a_rej = p_ability + alpha(p_mu_u,h)*(param.a.min - p_ability) + sigmah(p_a_sig,h)*z_pm_u;
    if p_offer_hit
        w_offer_p = draw_w_part(u_w_part);
        p_a_acc   = p_ability + alpha(p_mu,h)*(p_amean - p_ability) + sigmah(p_a_sig,h)*z_pm;

        V_acc = getMarriedValue(ability, employed, w_t, p_a_acc, 1, w_offer_p, kids_t, is_female, valueFunc, olf_t, 0 );
        V_rej = getMarriedValue(ability, employed, w_t, p_a_rej, 0, 0,         kids_t, is_female, valueFunc, olf_t, p_o_t);
        if V_acc > V_rej
            p_e_t=1; p_w_t=w_offer_p; p_a_t=p_a_acc;
        else
            p_a_t=p_a_rej;
        end
    else
        p_a_t=p_a_rej;
    end
    if (~p_e_t) && p_olf_hit
        V_u = getMarriedValue(ability, employed, w_t, p_a_t, 0, 0, kids_t, is_female, valueFunc, olf_t, 0 );
        V_o = getMarriedValue(ability, employed, w_t, p_a_t, 0, 0, kids_t, is_female, valueFunc, olf_t, 1 );
        if V_o > V_u, p_o_t=1; end
    end
else
    % ---------- Partner is OLF (no OLF -> U; OLF -> E offers at rate p_lambda_oE) ----------
    p_w_t  = 0;
    p_a_olf = p_ability + alpha(p_mu_u,h)*(param.a.min - p_ability) + sigmah(p_a_sig,h)*z_pm_u;

    offer_o_hit_p = (up_offer > exp(-p_lam_oE * h));
    if offer_o_hit_p
        w_offer_p = draw_w_part(u_w_part);
        p_a_acc   = p_ability + alpha(p_mu,h)*(p_amean - p_ability) + sigmah(p_a_sig,h)*z_pm;

        V_acc  = getMarriedValue(ability, employed, w_t, p_a_acc, 1, w_offer_p, kids_t, is_female, valueFunc, olf_t, 0);
        V_stay = getMarriedValue(ability, employed, w_t, p_a_olf, 0, 0,         kids_t, is_female, valueFunc, olf_t, 1);

        if V_acc > V_stay
            p_e_t = 1; p_o_t = 0; p_w_t = w_offer_p; p_a_t = p_a_acc;
        else
            p_a_t = p_a_olf;
        end
    else
        p_a_t = p_a_olf;
    end
end
p_a_t = clipa(p_a_t);

% --- Self dynamics
if employed
    if sep_hit
        emp_t=0; w_t=0;
        a_t = ability + alpha(mu_u,h)*(param.a.min-ability) + sigmah(a_sig,h)*z_mu_u;
    else
        a_emp = ability + alpha(mu,h)*(amean-ability) + sigmah(a_sig,h)*z_mu;
        if olf_hit
            V_e = getMarriedValue(a_emp, 1, w_t, p_a_t, p_e_t, p_w_t, kids_t, is_female, valueFunc, 0, p_o_t);
            V_o = getMarriedValue(a_emp, 0, 0,   p_a_t, p_e_t, p_w_t, kids_t, is_female, valueFunc, 1, p_o_t);
            if V_o > V_e
                emp_t=0; olf_t=1; w_t=0;
                a_t = ability + alpha(mu_u,h)*(param.a.min-ability) + sigmah(a_sig,h)*z_mu_u;
            else
                a_t = a_emp;
            end
        else
            a_t = a_emp;
        end
    end
else
    if ~olf
        a_rej = ability + alpha(mu_u,h)*(param.a.min-ability) + sigmah(a_sig,h)*z_mu_u;
        if offer_hit
            job_offer_t=1;
            w_offer = draw_w_self(u_w_self);
            a_acc   = ability + alpha(mu,h)*(amean-ability) + sigmah(a_sig,h)*z_mu;

            V_acc = getMarriedValue(a_acc, 1, w_offer, p_a_t, p_e_t, p_w_t, kids_t, is_female, valueFunc, 0, p_o_t);
            V_rej = getMarriedValue(a_rej, 0, 0,       p_a_t, p_e_t, p_w_t, kids_t, is_female, valueFunc, 0, p_o_t);
            if V_acc > V_rej, emp_t=1; w_t=w_offer; a_t=a_acc; else, a_t=a_rej; end
        else
            a_t=a_rej;
        end
        if (~emp_t) && olf_hit
            V_u = getMarriedValue(a_t, 0, 0, p_a_t, p_e_t, p_w_t, kids_t, is_female, valueFunc, 0, p_o_t);
            V_o = getMarriedValue(a_t, 0, 0, p_a_t, p_e_t, p_w_t, kids_t, is_female, valueFunc, 1, p_o_t);
            if V_o > V_u, olf_t=1; w_t=0; end
        end
   else
    % ---------- Self is OLF (no OLF -> U; allow OLF -> E at rate lambda_o) ----------
    w_t = 0;
    a_olf = ability + alpha(mu_u,h)*(param.a.min-ability) + sigmah(a_sig,h)*z_mu_u;

    offer_o_hit = (u_offer > exp(-lam_oE * h));
    if offer_o_hit
        job_offer_t = 1;
        w_offer = draw_w_self(u_w_self);
    
        a_acc   = ability + alpha(mu,h)*(amean-ability) + sigmah(a_sig,h)*z_mu;

        V_acc  = getMarriedValue(a_acc, 1, w_offer, p_a_t, p_e_t, p_w_t, kids_t, is_female, valueFunc, 0, p_o_t);
        V_stay = getMarriedValue(a_olf, 0, 0,       p_a_t, p_e_t, p_w_t, kids_t, is_female, valueFunc, 1, p_o_t);

        if V_acc > V_stay
            emp_t = 1; olf_t = 0; w_t = w_offer; a_t = a_acc;
        else
            a_t = a_olf;  % remain OLF
        end
    else
        a_t = a_olf;      % remain OLF
    end
end
end
a_t = clipa(a_t);

if ~kids_t && birth_hit
    kids_t  = 1;
    birth_t = 1;
end

if birth_t==1
    if is_female && emp_t
        w_t = w_t * (1 - param.fpem);
    elseif (~is_female) && p_e_t
        p_w_t = p_w_t * (1 - param.fpem);
    end
end

% --- Clamp partner -------------------------------------------------------
p_e_t = double(p_e_t~=0);
p_o_t = double(p_o_t~=0);

p_o_t = (~p_e_t) & (p_o_t ~= 0);           % OLF only if not employed
p_w_t   = p_w_t .* double(p_e_t);          % zero-out wage when not employed

% (optional) keep offer flag tidy
job_offer_t = double(job_offer_t~=0);


% --- Now compute couple value with clamped flags -------------------------
V_now = getMarriedValue(a_t, emp_t, w_t, ...
                        p_a_t, p_e_t, p_w_t, ...
                        kids_t, is_female, valueFunc, ...
                        olf_t, p_o_t);

% V_now = getMarriedValue(a_t, emp_t, w_t, p_a_t, p_e_t, p_w_t, kids_t, is_female, valueFunc, olf_t, p_o_t);
end

function V = getSingleValue(ability, employed, wage, olf, has_kids, is_female, valueFunc)
if is_female, vFunc = valueFunc.single_female; else, vFunc = valueFunc.single_male; end
if employed
    V = has_kids * vFunc.employed_kid(ability, wage) + (~has_kids) * vFunc.employed_nokid(ability, wage);
elseif ~olf
    V = has_kids * vFunc.unemployed_kid(ability)     + (~has_kids) * vFunc.unemployed_nokid(ability);
else
    V = has_kids * vFunc.olf_kid(ability)            + (~has_kids) * vFunc.olf_nokid(ability);
end
end

function V = getMarriedValue(ability_self, employed_self, wage_self, ...
                             ability_partner, employed_partner, wage_partner, ...
                             has_kids, is_female, valueFunc, olf_self, olf_partner)
if nargin < 10, olf_self = 0; end
if nargin < 11, olf_partner = 0; end

ef = ~isnan(employed_self)    & (employed_self~=0);
em = ~isnan(employed_partner) & (employed_partner~=0);
of = ~isnan(olf_self)         & (olf_self~=0);
om = ~isnan(olf_partner)      & (olf_partner~=0);
hk = ~isnan(has_kids)         & (has_kids~=0);

if isnan(wage_self),    wage_self = 0; end
if isnan(wage_partner), wage_partner = 0; end

vC = valueFunc.couple;

% map to female/male roles expected by couple value functions
if is_female
    af = ability_self; wf = wage_self; ef0 = ef; of0 = of;
    am = ability_partner; wm = wage_partner; em0 = em; om0 = om;
else
    af = ability_partner; wf = wage_partner; ef0 = em; of0 = om;
    am = ability_self;    wm = wage_self;    em0 = ef; om0 = of;
end

% choose function name base by state
if ef0 && ~of0 && em0 && ~om0
    nm = 'both_employed';                 args = {af,wf,am,wm};
elseif ef0 && ~of0 && ~em0 && ~om0
    nm = 'f_employed_m_unemployed';       args = {af,wf,am};
elseif ~ef0 && ~of0 && em0 && ~om0
    nm = 'f_unemployed_m_employed';       args = {af,am,wm};
elseif ~ef0 && ~of0 && ~em0 && ~om0
    nm = 'both_unemployed';               args = {af,am};
elseif ef0 && ~of0 && ~em0 &&  om0
    nm = 'f_employed_m_olf';              args = {af,wf,am};
elseif ~ef0 &&  of0 && em0 && ~om0
    nm = 'f_olf_m_employed';              args = {af,am,wm};
elseif ~ef0 &&  of0 && ~em0 &&  om0
    nm = 'both_olf';                      args = {af,am};
elseif ~ef0 && ~of0 && ~em0 &&  om0
    nm = 'f_unemployed_m_olf';            args = {af,am};
elseif ~ef0 &&  of0 && ~em0 && ~om0
    nm = 'f_olf_m_unemployed';            args = {af,am};
else
    nm = 'both_unemployed';               args = {af,am};
end

suf = ternary(hk,'_kid','_nokid');
fn  = [nm suf];

% if missing, fall back to nearest olf form
if ~isfield(vC, fn) || ~isa(vC.(fn), 'function_handle')
    switch nm
        case {'f_employed_m_olf','f_employed_m_unemployed'}
            fn = ['f_employed_m_unemployed' suf]; args = {af,wf,am};
        case {'f_olf_m_employed','f_unemployed_m_employed'}
            fn = ['f_unemployed_m_employed' suf]; args = {af,am,wm};
        case {'both_olf','both_unemployed','f_unemployed_m_olf','f_olf_m_unemployed'}
            fn = ['both_unemployed' suf];         args = {af,am};
        otherwise
            fn = ['both_unemployed' suf];         args = {af,am};
    end
end
V = vC.(fn)(args{:});
end

function y = ternary(cond,a,b), if cond, y=a; else, y=b; end, end

% ====================== STATS (monthly) ================================
function simResults = calculatePanelStatistics(agentPanel)
    simResults = struct();
    [N, T] = size(agentPanel.employed);
    h = agentPanel.months_per_step;                      % months per step in the KEPT window
    steps_per_month = max(1e-12, 1/h);

    females = agentPanel.female;
    males   = ~females;

    for t = 1:T
        emp_t = agentPanel.employed(:,t)==1;
        olf_t = agentPanel.olf(:,t)==1;

        simResults.emp_rate_female(t,1) = mean(emp_t(females));
        simResults.emp_rate_male(t,1)   = mean(emp_t(~females));
        simResults.emp_rate_all(t,1)    = mean(emp_t);

        in_lf = (emp_t | ~olf_t);
        simResults.lfp_rate_female(t,1) = mean(in_lf(females));
        simResults.lfp_rate_male(t,1)   = mean(in_lf(~females));
        simResults.lfp_rate_all(t,1)    = mean(in_lf);

        simResults.marriage_rate(t,1)   = mean(agentPanel.married(:,t));
        simResults.fertility_rate(t,1)  = mean(agentPanel.has_kids(:,t));

        if any(emp_t)
            rw = exp(agentPanel.ability(emp_t,t)) .* agentPanel.wage(emp_t,t);
            simResults.avg_wage(t,1) = mean(rw);
            simResults.wage_std(t,1) = std(rw);

            emp_f = emp_t & females;
            emp_m = emp_t & males;
            simResults.avg_wage_female(t,1) = iff(any(emp_f), mean(exp(agentPanel.ability(emp_f,t)).*agentPanel.wage(emp_f,t)), NaN);
            simResults.avg_wage_male(t,1)   = iff(any(emp_m), mean(exp(agentPanel.ability(emp_m,t)).*agentPanel.wage(emp_m,t)), NaN);
        else
            simResults.avg_wage(t,1) = NaN; simResults.wage_std(t,1) = NaN;
            simResults.avg_wage_female(t,1) = NaN; simResults.avg_wage_male(t,1) = NaN;
        end
    end

    % Transitions over last half (stability) � converted to monthly probability
    T_start = max(2, floor(T/2));
    X = agentPanel.employed(:,T_start:end)~=0;  % logical employed sequence
    simResults.job_finding_rate    = transition_rate_monthly(X, 0, 1, h);
    simResults.job_separation_rate = transition_rate_monthly(X, 1, 0, h);

    % Marriage formation (monthly)
    M = agentPanel.married(:,T_start:end)~=0;
    simResults.marriage_formation_rate = transition_rate_monthly(M, 0, 1, h);

    % Convenience: steady-state summaries over last 20%
    idx = max(1, floor(0.8*T)) : T;
    fields = {'emp_rate_female','emp_rate_male','emp_rate_all', ...
              'lfp_rate_female','lfp_rate_male','lfp_rate_all', ...
              'marriage_rate','fertility_rate','avg_wage'};
    for f = 1:numel(fields)
        if isfield(simResults, fields{f})
            simResults.([fields{f} '_ss']) = nanmean(simResults.(fields{f})(idx));
        end
    end

    % Final cross-section by gender & marital status
    females = agentPanel.female;
    married = agentPanel.married(:,end)==1;
    simResults.final_emp_by_gender_marital = zeros(2,2);
    simResults.final_emp_by_gender_marital(1,1) = mean(agentPanel.employed(females & ~married, end));
    simResults.final_emp_by_gender_marital(1,2) = mean(agentPanel.employed(females &  married, end));
    simResults.final_emp_by_gender_marital(2,1) = mean(agentPanel.employed(~females & ~married, end));
    simResults.final_emp_by_gender_marital(2,2) = mean(agentPanel.employed(~females &  married, end));

    simResults.months_per_step = h;
end

function y = iff(c,a,b), if c, y=a; else, y=b; end, end

function rate_month = transition_rate_monthly(X, from_state, to_state, h_months)
    % X: logical matrix (#agents x #time), state==1 if "on"
    prev = (X(:,1:end-1) == (from_state~=0));
    next = (X(:,2:end)   == (to_state~=0));
    opp  = sum(prev(:));
    if opp==0
        p_step = NaN;
    else
        p_step = sum(prev(:) & next(:)) / opp;   % probability per step
    end
    % Convert per-step probability to monthly probability:
    % p_month = 1 - (1 - p_step)^(1/h)
    if ~isfinite(p_step), rate_month = NaN; return; end
    if h_months <= 0,      rate_month = NaN; return; end
    rate_month = 1 - (1 - p_step)^(1/h_months);
end

% =================== Samplers (unchanged) ====================
function sampler = buildShadowPartnerSampler(a_grid, w_grid, S)
    a = a_grid(:); w = w_grid(:);
    if ~issorted(a) || ~issorted(w), error('Grids must be strictly increasing.'); end

    c_hl = max(S.hl(:),0);  c_hc = max(S.hc(:),0);
    c_ol = max(S.ol(:),0);  c_oc = max(S.oc(:),0);
    samp_hl = buildInvSampler(a, c_hl);  samp_hc = buildInvSampler(a, c_hc);
    samp_ol = buildInvSampler(a, c_ol);  samp_oc = buildInvSampler(a, c_oc);
    mass_hl = sum(c_hl); mass_hc = sum(c_hc); mass_ol = sum(c_ol); mass_oc = sum(c_oc);

    EL = max(S.el,0); EC = max(S.ec,0);
    [samp_el, mass_el] = build2DMixtureCellSampler(a, w, EL);
    [samp_ec, mass_ec] = build2DMixtureCellSampler(a, w, EC);

    comp_mass = [mass_el, mass_ec, mass_hl, mass_hc, mass_ol, mass_oc];
    if all(comp_mass==0), error('Shadow distribution has zero mass.'); end
    comp_p   = comp_mass / sum(comp_mass);
    comp_cdf = cumsum(comp_p);

    sampler = @draw_one;
    function p = draw_one(nDraw)
        if nargin<1, nDraw = 1; end
        p(nDraw,1) = struct('ability',[],'employed',[],'olf',[],'wage',[],'has_kids',[]);
        for dd=1:nDraw
            u = rand;
            if u <= comp_cdf(1)
                [a_draw, w_draw] = samp_el(1);
                p(dd) = pack(a_draw, 1, 0, w_draw, 0);
            elseif u <= comp_cdf(2)
                [a_draw, w_draw] = samp_ec(1);
                p(dd) = pack(a_draw, 1, 0, w_draw, 1);
            elseif u <= comp_cdf(3)
                a_draw = samp_hl(1);
                p(dd) = pack(a_draw, 0, 0, 0, 0);
            elseif u <= comp_cdf(4)
                a_draw = samp_hc(1);
                p(dd) = pack(a_draw, 0, 0, 0, 1);
            elseif u <= comp_cdf(5)
                try 
                    a_draw = samp_ol(1);
                catch
                    a_draw = samp_hl(1);
                end
                p(dd) = pack(a_draw, 0, 1, 0, 0);
            else
                try
                    a_draw = samp_oc(1);
                catch
                    a_draw = samp_hc(1);
                end
                p(dd) = pack(a_draw, 0, 1, 0, 1);
            end
        end
    end

    function s = pack(a_, emp_, olf_, w_, kid_)
        s.ability = a_; s.employed = emp_; s.olf = olf_; s.wage = w_; s.has_kids = kid_;
    end
end

function [sampler, total_mass] = build2DMixtureCellSampler(x, y, Z)
    x = x(:); y = y(:);
    [nX,nY] = size(Z);
    if nX ~= numel(x) || nY ~= numel(y), error('build2DMixtureCellSampler: size mismatch.'); end
    Z = max(Z,0);
    dx = diff(x); dy = diff(y);
    Zavg = 0.25*(Z(1:end-1,1:end-1)+Z(2:end,1:end-1)+Z(1:end-1,2:end)+Z(2:end,2:end));
    Area = dx * dy.';
    M = Zavg .* Area;
    mass_vec = M(:);
    total_mass = sum(mass_vec);
    if total_mass <= 0
        sampler = @(n) deal(nan(n,1), nan(n,1));
        return;
    end
    pmf = mass_vec / total_mass;
    cdf = cumsum(pmf);
    nCx = numel(dx);
    sampler = @draw2d;
    function [xdraw, ydraw] = draw2d(nDraw)
        if nargin<1, nDraw = 1; end
        u = rand(nDraw,1);
        k = discretize(u, [0; cdf]); k = max(k,1);
        j = ceil(k / nCx);
        i = k - (j-1)*nCx;
        rx = rand(nDraw,1); ry = rand(nDraw,1);
        xdraw = x(i) + rx.*dx(i);
        ydraw = y(j) + ry.*dy(j);
    end
end

function sampler = buildInvSampler(x_grid, mass, fallback_value)
    if nargin < 3, fallback_value = NaN; end
    x_grid = x_grid(:); mass = max(mass(:), 0);
    if ~issorted(x_grid), error('buildInvSampler: x_grid must be sorted ascending.'); end
    total = sum(mass);
    if total <= 0
        sampler = @(nDraw) repmat(fallback_value, max(1,nDraw), 1);
        return;
    end
    pmf = mass / total;
    cdf = cumsum(pmf); cdf(end) = 1;           % ensure exact 1
    edges = [0; cdf(:)];
    sampler = @sample_once;
    function x = sample_once(nDraw)
        if nargin < 1 || isempty(nDraw), nDraw = 1; end
        u = rand(nDraw,1);
        k = discretize(u, edges);
        k(~isfinite(k) | k < 1) = 1;
        x = x_grid(k);
    end
end
