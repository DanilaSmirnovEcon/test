function moments = compute_moments(agentPanel, simParam, param)
% COMPUTE_MOMENTS
% Build a moments struct from the (post-burn-in) panel.
% Includes:
%   • Person-time shares (married / never-married, women/men, has child)
%   • Out-of-LF and Unemployment rates using person-time denominators
%   • Transition rates (monthly), couple/single employment shares
%   • Wage stats, wage dispersion/progression, assortative matching
%   • Child cost by income bin (if param.c provided)
%   • Avg. age at (first) marriage (years)

if nargin < 3, param = struct(); end

% ---------- Basics ----------
[N, T] = size(agentPanel.employed);

% Focal gender (row's person). Accept Nx1 or NxT, use first column if needed.
% ---------- Gender masks (focal) ----------
if isfield(agentPanel,'female') && ~isempty(agentPanel.female)
    fem = logical(agentPanel.female(:));         % N×1
elseif isfield(agentPanel,'is_female') && ~isempty(agentPanel.is_female)
    fem = logical(agentPanel.is_female(:));      % fallback field name
else
    error('compute_moments: missing agentPanel.female (or is_female).');
end
male = ~fem;

% Expand to N×T masks for focal person
gF = repmat(fem ,1,T);   % focal is female
gM = repmat(male,1,T);   % focal is male

if isfield(agentPanel,'female') && ~isempty(agentPanel.female)
    pmale = logical(agentPanel.female(:));         % N×1
elseif isfield(agentPanel,'is_female') && ~isempty(agentPanel.is_female)
    pmale= logical(agentPanel.is_female(:));      % fallback field name
else
    error('compute_moments: missing gender partner .');
end
pfem = ~pmale;
% Expand to N×T masks for partner person
pgF = repmat(pfem ,1,T);   % partner is female
pgM = repmat(pmale,1,T);   % partner is male

% ---------- Calendar units ----------
if isfield(agentPanel,'months_per_step') && ~isempty(agentPanel.months_per_step)
    months_per_step = agentPanel.months_per_step;
elseif isfield(simParam,'months_per_step') && ~isempty(simParam.months_per_step)
    months_per_step = simParam.months_per_step;
elseif isfield(simParam,'fine_step_months') && ~isempty(simParam.fine_step_months)
    months_per_step = simParam.fine_step_months;
else
    months_per_step = 1; % default: step = 1 month
end
steps_per_year = max(1, round(12 / months_per_step));

% ---------- States ----------
E  = agentPanel.employed == 1;
O  = agentPanel.olf      == 1;
U  = ~E & ~O;

M   = agentPanel.married   == 1;
NM  = ~M;
KK  = agentPanel.has_kids  == 1;
KN  = ~KK;

% Partner states (robust to missing partner info)
pE     = agentPanel.partner_employed == 1;
pO     = agentPanel.partner_olf      == 1;
pKnown = ~isnan(agentPanel.partner_employed) | ~isnan(agentPanel.partner_olf);
pU     = (~pE & ~pO) & pKnown;

% Realized wages / abilities (self & partner)
A   = agentPanel.ability;
W   = exp(agentPanel.ability) .* agentPanel.wage;                 % self
W   = W.*param.fkidwpenalty.*(agentPanel.female) + W.*(1-agentPanel.female);
pW  = exp(agentPanel.partner_ability) .* agentPanel.partner_wage; % partner
pW  = pW.*param.fkidwpenalty.*(1-agentPanel.female) + pW.*(agentPanel.female);
WO  = agentPanel.wage;     % offer wage (self)
pWO = agentPanel.partner_wage;

% Unemployment benefits (as fn of ability)
if isfield(param,'b') && isa(param.b,'function_handle')
    B  = param.b(agentPanel.ability);
    pB = param.b(agentPanel.partner_ability);
else
    B  = NaN(size(U));
    pB = NaN(size(U));
end

% Helpers
avg   = @(x) mean(x(:),'omitnan');
moments = struct();

%% ====================== 1) Shares + Rates ===============================
SMALL = 1e-12;
single  = NM;    % never married person-time
married = M;     % married person-time

% ---- Person-time denominators (each "person" counted once; married counts twice)
den_all   = sum(single,'all') + 2*sum(married,'all');                         den_all   = max(den_all,SMALL);
den_nm    = sum(single,'all');                                                den_nm    = max(den_nm ,SMALL);
den_m     = 2*sum(married,'all');                                             den_m     = max(den_m  ,SMALL);
den_women = sum(gF & single,'all') + sum(gF & married,'all') + sum(pgF & married,'all');  den_women = max(den_women,SMALL);
den_men   = sum(gM & single,'all') + sum(gM & married,'all') + sum(pgM & married,'all');  den_men   = max(den_men  ,SMALL);
den_k     = sum(KK & single ,'all') + 2*sum(KK & married ,'all');             den_k     = max(den_k    ,SMALL);
den_nk    = sum(KN & single ,'all') + 2*sum(KN & married ,'all');             den_nk    = max(den_nk   ,SMALL);

% ---- Overall shares (person-time, person-level accounting)
allpeople     = den_all;
women_count   = sum(gF,'all') + sum(pgF & married,'all');
men_count     = sum(gM,'all') + sum(pgM & married,'all');
haschild_cnt  = sum(KK & single,'all') + 2*sum(KK & married,'all');

moments.married.all        = (1*sum(married,'all')) / (sum(single,'all') + 1*sum(married,'all'));%allpeople
moments.never_married.all  = 1-moments.married.all ;%(  sum(single ,'all')) / allpeople;

moments.women.all          = women_count / allpeople;
moments.men.all            = men_count   / allpeople;

moments.haschild.all       = haschild_cnt / allpeople;
moments.nochild.all        = 1 - moments.haschild.all;

% ---- Numerators: person-time in each state
% Overall
numO_all = sum(O & single,'all') + sum(O & married,'all') + sum(pO & married,'all');
numU_all = sum(U & single,'all') + sum(U & married,'all') + sum(pU & married,'all');
% Never-married (focal only)
numO_nm  = sum(O & single,'all');
numU_nm  = sum(U & single,'all');
% Married (count both spouses)
numO_m   = sum(O & married,'all') + sum(pO & married,'all');
numU_m   = sum(U & married,'all') + sum(pU & married,'all');
% Women
numO_w   = sum(O & single  & gF,'all') + sum(O  & married & gF ,'all') + sum(pO  & married & pgF,'all');
numU_w   = sum(U & single  & gF,'all') + sum(U  & married & gF ,'all') + sum(pU  & married & pgF,'all');
% Men
numO_men = sum(O & single  & gM,'all') + sum(O  & married & gM ,'all') + sum(pO  & married & pgM,'all');
numU_men = sum(U & single  & gM,'all') + sum(U  & married & gM ,'all') + sum(pU  & married & pgM,'all');
% Has child
numO_k   = sum(O & single  & KK,'all') + sum(O  & married & KK,'all') + sum(pO  & married & KK,'all');
numU_k   = sum(U & single  & KK,'all') + sum(U  & married & KK,'all') + sum(pU  & married & KK,'all');
% No child
numO_nk  = sum(O & single  & KN,'all') + sum(O  & married & KN,'all') + sum(pO  & married & KN,'all');
numU_nk  = sum(U & single  & KN,'all') + sum(U  & married & KN,'all') + sum(pU  & married & KN,'all');

% ---- Rates: state | group
moments.all.out_of_lf            = numO_all / den_all;
moments.never_married.out_of_lf  = numO_nm  / den_nm;
moments.married.out_of_lf        = numO_m   / den_m;
moments.women.out_of_lf          = numO_w   / den_women;
moments.men.out_of_lf            = numO_men / den_men;
moments.haschild.out_of_lf       = numO_k   / den_k;
moments.nochild.out_of_lf        = numO_nk  / den_nk;

moments.all.unemployed           = numU_all / den_all;
moments.never_married.unemployed = numU_nm  / den_nm;
moments.married.unemployed       = numU_m   / den_m;
moments.women.unemployed         = numU_w   / den_women;
moments.men.unemployed           = numU_men / den_men;
moments.haschild.unemployed      = numU_k   / den_k;
moments.nochild.unemployed       = numU_nk  / den_nk;

%% ========== 2) Transitions by gender (MONTHLY) ==========================
tr_step_f = transitions_by_gender(fem,  U, E, O, T);
tr_step_m = transitions_by_gender(male, U, E, O, T);

clip01 = @(x) min(max(x,0), 1-1e-12);
scale_to_month = @(p) 1 - (1 - clip01(p)).^(1/max(months_per_step,eps));

moments.transitions.female.UE = scale_to_month(tr_step_f.UE);
moments.transitions.female.EU = scale_to_month(tr_step_f.EU);
moments.transitions.female.EO = scale_to_month(tr_step_f.EO);
moments.transitions.female.UO = scale_to_month(tr_step_f.UO);
moments.transitions.female.OE = scale_to_month(tr_step_f.OE);

moments.transitions.male.UE = scale_to_month(tr_step_m.UE);
moments.transitions.male.EU = scale_to_month(tr_step_m.EU);
moments.transitions.male.EO = scale_to_month(tr_step_m.EO);
moments.transitions.male.UO = scale_to_month(tr_step_m.UO);
moments.transitions.male.OE = scale_to_month(tr_step_m.OE);

% tr_step_pf = transitions_by_gender(pfem,  pU, pE, pO, T);
% tr_step_pm = transitions_by_gender(pmale, pU, pE, pO, T);
% 
% clip01 = @(x) min(max(x,0), 1-1e-12);
% scale_to_month = @(p) 1 - (1 - clip01(p)).^(1/max(months_per_step,eps));
% 
% moments.transitions.pfemale.UE = scale_to_month(tr_step_pf.UE);
% moments.transitions.pfemale.EU = scale_to_month(tr_step_pf.EU);
% moments.transitions.pfemale.EO = scale_to_month(tr_step_pf.EO);
% moments.transitions.pfemale.UO = scale_to_month(tr_step_pf.UO);
% 
% moments.transitions.pmale.UE = scale_to_month(tr_step_pm.UE);
% moments.transitions.pmale.EU = scale_to_month(tr_step_pm.EU);
% moments.transitions.pmale.EO = scale_to_month(tr_step_pm.EO);
% moments.transitions.pmale.UO = scale_to_month(tr_step_pm.UO);

%% ========== 3) Couple employment pattern shares (married & known) =======
mar_known = (M == 1) & pKnown;
denom = sum(mar_known(:));
if denom == 0
    cee = NaN; ceu = NaN; cuu = NaN; ceo = NaN; coo = NaN; cuo = NaN;
else
    cee = sum( mar_known &  (E &  pE) , 'all') / denom;
    ceu = sum( mar_known & ((E &  pU) | (U &  pE)) , 'all') / denom;
    cuu = sum( mar_known &  (U &  pU) , 'all') / denom;
    ceo = sum( mar_known & ((E &  pO) | (O &  pE)) , 'all') / denom;
    coo = sum( mar_known &  (O &  pO) , 'all') / denom;
    cuo = sum( mar_known & ((U &  pO) | (O &  pU)) , 'all') / denom;
end
moments.couplesemploymentshares = struct('EE',cee,'EU',ceu,'UU',cuu,'EO',ceo,'OO',coo,'UO',cuo);

%% ========== 3b) Singles’ employment shares ==============================
singles = NM;

% females
singF = singles & gF;
denF  = sum(singF(:));
if denF==0
    fE=NaN; fU=NaN; fO=NaN;
else
    fE = sum((singF & E), 'all')/denF;
    fU = sum((singF & U), 'all')/denF;
    fO = sum((singF & O), 'all')/denF;
end
moments.singlesemploymentshares.female = struct('E', fE, 'U', fU, 'O', fO);

% males
singM = singles & gM;
denM  = sum(singM(:));
if denM==0
    mE=NaN; mU=NaN; mO=NaN;
else
    mE = sum((singM & E), 'all')/denM;
    mU = sum((singM & U), 'all')/denM;
    mO = sum((singM & O), 'all')/denM;
end
moments.singlesemploymentshares.male = struct('E', mE, 'U', mU, 'O', mO);

%% ========== 4) Wages: means, dispersion, yearly change ==================
W_pos  = W;  W_pos(~E)  = NaN;
pW_pos = pW; pW_pos(~pE)= NaN;

% Offer wage dispersion (combine self & partner offers symmetrically)
moments.stats_wages.cvOffers  = std(WO(:),'omitnan')./mean(WO(:),'omitnan');
moments.stats_wages.pcvOffers = std(pWO(:),'omitnan')./mean(pWO(:),'omitnan');
moments.stats_wages.cvOffers  = 0.5*(moments.stats_wages.cvOffers + moments.stats_wages.pcvOffers);

% Simple, person-based wage averages
moments.stats_wages.women          = mean(W_pos(gF), 'omitnan');
moments.stats_wages.men            = mean(W_pos(gM), 'omitnan');
moments.stats_wages.never_married  = mean(W_pos(NM), 'omitnan');
moments.stats_wages.married        = mean(W_pos(M ), 'omitnan');
moments.stats_wages.haschild       = mean(W_pos(KK), 'omitnan');
moments.stats_wages.nochild        = mean(W_pos(KN), 'omitnan');
moments.stats_wages.all            = mean(W_pos(:),  'omitnan');

% Yearly dispersion/change (self)
ppy = steps_per_year;
if T > ppy
    E12  = E(:,1:T-ppy) & E(:,1+ppy:T);
    yearlyW = zeros(size(W,1),size(W,2)-ppy);
    yearlyE = zeros(size(W,1),size(W,2)-ppy);
    for j = 1:ppy
        yearlyW = yearlyW + W(:,j:T-ppy+j-1)./ppy;
        yearlyE = yearlyE + E(:,j:T-ppy+j-1)./ppy;
    end
    yearlyW = yearlyW*48*40;                 % to annualized hours if desired
    base = W(:,1:T-ppy);  nxt = W(:,1+ppy:T);
    pct   = abs(nxt - base) ./ base;
    valid = E12 & isfinite(pct) & (base > 0);
    den   = sum(valid(:));
    if den>0
        moments.wagechangeyear = sum(pct(valid) <= 0.20) / den;
    else
        moments.wagechangeyear = NaN;
    end
else
    moments.wagechangeyear = NaN;
end
yearlyW(yearlyE==0)=NaN;
moments.wagedisp = var(log(yearlyW(:)),'omitnan');
moments.wagechangeyear = moments.wagechangeyear;

% Partner yearly dispersion/change (for reference)
if T > ppy
    pE12  = pE(:,1:T-ppy) & pE(:,1+ppy:T);
    yearlypW = zeros(size(pW,1),size(pW,2)-ppy);
    yearlypE = zeros(size(pW,1),size(pW,2)-ppy);
    for j = 1:ppy
        yearlypW = yearlypW + pW(:,j:T-ppy+j-1)./ppy;
        yearlypE = yearlypE + pE(:,j:T-ppy+j-1)./ppy;
    end
    yearlypW = yearlypW*48*40;
    pbase = pW(:,1:T-ppy);  pnxt = pW(:,1+ppy:T);
    ppct   = abs(pnxt - pbase) ./ pbase;
    pvalid = pE12 & isfinite(ppct) & (pbase > 0);
    pden   = sum(pvalid(:));
    if pden>0
        moments.pwagechangeyear = sum(pct(valid) <= 0.20) / den;
    else
        moments.pwagechangeyear = NaN;
    end
else
    moments.pwagechangeyear = NaN;
end
yearlypW(yearlypE==0)=NaN;
moments.pwagedisp = var(log(yearlypW(:)),'omitnan');

% Blend self/partner dispersion (optional)
resp = (sum(M==0,'all') + sum(M==1,'all')) / allpeople; % person share of "self" records
moments.wagedisp  = resp*moments.wagedisp + (1-resp)*moments.pwagedisp;

%% ========== 5) Wage progression (entry wage / overall avg wage) =========
if isfield(agentPanel,'newborn')
    NEWB = agentPanel.newborn == 1;
else
    NEWB = false(size(E)); NEWB(:,1) = true;
end
entry_wages_newborn = [];
for i = 1:N
    t0s = find(NEWB(i,:)); if isempty(t0s), continue; end
    t0s = [t0s, T+1];
    for k = 1:numel(t0s)-1
        t_start = t0s(k); t_end = t0s(k+1)-1;
        if t_start > t_end, continue; end
        rel = find(E(i,t_start:t_end),1,'first');
        if ~isempty(rel)
            w_entry = W(i, t_start+rel-1);
            if isfinite(w_entry) && w_entry>0
                entry_wages_newborn(end+1,1) = w_entry; %#ok<AGROW>
            end
        end
    end
end
moments.wageentryvar        = var(log(entry_wages_newborn),'omitnan');
moments.wageprogression     = mean(entry_wages_newborn,'omitnan') / moments.stats_wages.all;

enter_any   = E(:,2:end) & ~E(:,1:end-1);
entryW_any  = W(:,2:end); entryW_any(~enter_any) = NaN;
moments.wageprogression_any = mean(entryW_any(:),'omitnan') / moments.stats_wages.all;

% ====== 5b) Wage change after unemployment (robust) ======================
if isfield(simParam,'reemployment_horizon_months') && ~isempty(simParam.reemployment_horizon_months)
    hzn_months = simParam.reemployment_horizon_months;
else
    hzn_months = 12;
end
% Need >= 2 steps because after a separation at t, E(t+1)=0 by construction
max_hzn_steps = max(2, ceil(hzn_months / months_per_step));

spell_pct  = [];
per_step   = [];
durs_steps = [];

n_sep = 0; n_censored = 0; n_found = 0; n_kept = 0;

for i = 1:N
    % separations: E=1 then 0 next period
    sep_t = find(E(i,1:end-1) & ~E(i,2:end));
    for k = 1:numel(sep_t)
        n_sep = n_sep + 1;
        t  = sep_t(k);

        % wage at separation
        w0 = W(i,t);
        if ~(isfinite(w0) && w0>0), continue; end

        % search window: start at t+2 (earliest possible re-employment)
        search_start = t + 2;
        t_end = min(T, t + max_hzn_steps);
        if search_start > t_end
            n_censored = n_censored + 1;
            continue;
        end

        % find first re-employment inside [t+2, t_end]
        relpos = find(E(i,search_start:t_end), 1, 'first');  % position within the subarray
        if isempty(relpos)
            n_censored = n_censored + 1;
            continue;
        end

        t2 = search_start + relpos - 1;       % absolute time of re-employment
        w1 = W(i,t2);
        if ~(isfinite(w1) && w1>0), continue; end

        % unemployment duration in steps: from (t+1) to (t2-1) inclusive
        d_u = t2 - (t+1);                      % strictly >= 1 by construction
        if d_u <= 0, continue; end             % safety, should not trigger

        % store outcomes
        durs_steps(end+1,1) = d_u;                  %#ok<AGROW>
        per_step(end+1,1)   = (log(w1) - log(w0)) / d_u;  %#ok<AGROW>
        spell_pct(end+1,1)  = (w1 - w0) / w0;            %#ok<AGROW>
        n_found = n_found + 1;
    end
end

% aggregate
if isempty(per_step)
    moments.wagechangeunemployment                       = NaN;
    moments.wagechangeunemployment_spell_pct             = NaN;
    moments.wagechangeunemployment_per_period_logpct     = NaN;
    moments.wagechangeunemployment_per_month_logpct      = NaN;
    moments.wagechangeunemployment_per_year_logpct       = NaN;
else
    log_spell = per_step .* durs_steps;  % = log(w1) - log(w0)
    keep = isfinite(per_step) & isfinite(spell_pct) & (durs_steps > 0) & (abs(log_spell) <= log(5));
    pp   = per_step(keep);
    spct = spell_pct(keep);

    n_kept = sum(keep);
    moments.wagechangeunemployment           = mean(spct,'omitnan');
    moments.wagechangeunemployment_spell_pct = moments.wagechangeunemployment;
    moments.wagechangeunemployment_per_period_logpct = mean(pp,'omitnan');
    moments.wagechangeunemployment_per_month_logpct  = (1/months_per_step) * mean(pp,'omitnan');
    moments.wagechangeunemployment_per_year_logpct   = 12 * moments.wagechangeunemployment_per_month_logpct;
end

% (optional) quick diagnostics:
fprintf('Spells: sep=%d, reemp_found=%d, kept=%d, censored=%d\n', n_sep, n_found, n_kept, n_censored);
if ~isempty(durs_steps)
    fprintf('Unemp duration (steps): min=%g, p50=%g, mean=%g, max=%g\n', ...
        min(durs_steps), median(durs_steps), mean(durs_steps), max(durs_steps));
end

%% ========== 6) Assortative matching (within couples) ====================
mask_couple_we = M & E & pE;
if any(mask_couple_we(:))
    x = log(W);  y = log(pW);
    x = x(mask_couple_we);
    x = nonzeros(x);
    y = y(mask_couple_we);
    y = nonzeros(y);
    C = corrcoef(x, y, 'Rows','complete');%C = corrcoef(x(mask_couple_we), y(mask_couple_we), 'Rows','complete');
    moments.corr_wages = C(1,2);
else
    moments.corr_wages = NaN;
end


mask_couple = M;
if any(mask_couple(:))
    xa = agentPanel.ability(mask_couple);
    ya = agentPanel.partner_ability(mask_couple);
    C  = corrcoef(xa, ya, 'Rows','complete');
    moments.corr_ability = C(1,2);
else
    moments.corr_ability = NaN;
end

%% ========== 7) Child cost by income bin (optional) ======================
% If you have a function param.c(ability) returning the (annual) cost,
% we scale by realized income y for households with kids.
if isfield(param,'c') && isa(param.c,'function_handle')
    R_lin  = W(:);
    K_vec  = KK(:);
    A_all  = A(:);
    valid  = isfinite(R_lin) & (R_lin>0) & K_vec;
    y      = R_lin(valid);
    a      = A_all(valid);
    cost_abs = param.c(a)./y;
    % Split by median income of valid obs
    med_inc = median(y, 'omitnan');
    low  = (y <= med_inc);
    high = (y >  med_inc);
    moments.child_cost.median_income  = med_inc;
    moments.child_cost.low.mean_share = mean(cost_abs(low),  'omitnan');
    moments.child_cost.high.mean_share= mean(cost_abs(high), 'omitnan');
    moments.child_cost.mean_share     = mean(cost_abs, 'omitnan');
else
    moments.child_cost = struct('median_income',NaN,'low',struct('mean_share',NaN),...
                                'high',struct('mean_share',NaN),'mean_share',NaN);
end

%% ========== 8) Average age at (first) marriage (years) ==================
moments.age_at_marriage.women_years = NaN;
moments.age_at_marriage.men_years   = NaN;
moments.age_at_marriage.all_years   = NaN;

if isfield(agentPanel,'age')
    AgeM = agentPanel.age;  % months
    if isfield(agentPanel,'newborn')
        NEWB = agentPanel.newborn == 1;
    else
        NEWB = false(size(AgeM));
        NEWB(:,1) = true;
        NEWB(:,2:end) = AgeM(:,2:end) < AgeM(:,1:end-1);
    end

    ages_w = []; ages_m = [];
    for i = 1:N
        t0s = find(NEWB(i,:)); if isempty(t0s), t0s = 1; end
        t0s = [t0s, T+1];
        for s = 1:numel(t0s)-1
            t0 = t0s(s); t1 = t0s(s+1)-1;
            if t0>=t1, continue; end
            trans = find( (~M(i,t0:t1-1)) & M(i,t0+1:t1), 1, 'first');
            if ~isempty(trans)
                t_mar = t0 + trans;
                age_mar_years = AgeM(i,t_mar)/12;
                if isfinite(age_mar_years) && age_mar_years>=0
                    if fem(i), ages_w(end+1,1) = age_mar_years; %#ok<AGROW>
                    else       
                        ages_m(end+1,1) = age_mar_years; %#ok<AGROW>
                    end
                end
            end
        end
    end
    if ~isempty(ages_w), moments.age_at_marriage.women_years = median(ages_w,'omitnan'); end
    if ~isempty(ages_m), moments.age_at_marriage.men_years   = median(ages_m,'omitnan'); end
    both = [ages_w; ages_m];
    if ~isempty(both),   moments.age_at_marriage.all_years   = median(both,'omitnan');   end
end
% %% ========== 8) Age at (first) marriage: naive + left-truncation aware ==========
% % Output fields (years)
% moments.age_at_marriage.women_years = NaN;
% moments.age_at_marriage.men_years   = NaN;
% moments.age_at_marriage.all_years   = NaN;
% 
% moments.age_at_marriage_km.women_median_years = NaN;
% moments.age_at_marriage_km.men_median_years   = NaN;
% moments.age_at_marriage_km.all_median_years   = NaN;
% 
% moments.age_at_marriage_expMLE.women_mean_years = NaN;
% moments.age_at_marriage_expMLE.men_mean_years   = NaN;
% moments.age_at_marriage_expMLE.all_mean_years   = NaN;
% 
% % ---- Inputs (prefer chronological age; fallback if needed) ----------------------
% if isfield(agentPanel,'age_chron') && ~isempty(agentPanel.age_chron)
%     AgeM = agentPanel.age_chron;          % months, chronological
% elseif isfield(agentPanel,'age') && ~isempty(agentPanel.age)
%     % Fallback: treat 'age' as chronological if no explicit 'age_chron'
%     if exist('param','var') && isfield(param,'entry_age_years')
%         AgeM = agentPanel.age + param.entry_age_years*12;  % best-effort
%     else
%         AgeM = agentPanel.age;                              % last resort
%     end
% else
%     warning('compute_moments:age','No age field found; skipping age-at-marriage moments.');
%     goto_end_age_section = true; %#ok<NASGU>
% end
% 
% if exist('goto_end_age_section','var')
%     % nothing to do
% else
%     if isfield(agentPanel,'newborn') && ~isempty(agentPanel.newborn)
%         NEWB = agentPanel.newborn ~= 0;
%     else
%         % Detect (re)births by age resets
%         NEWB = false(size(AgeM));
%         NEWB(:,1) = true;
%         NEWB(:,2:end) = AgeM(:,2:end) < AgeM(:,1:end-1);
%     end
%     M   = agentPanel.married ~= 0;
%     fem = agentPanel.female  ~= 0;
%     [N,T] = size(AgeM);
% 
%     % --------------------- (i) Naive within-window medians -----------------------
%     ages_w = []; ages_m = [];   % realized ages (years) for spells with an event
%     for i = 1:N
%         t0s = find(NEWB(i,:)); if isempty(t0s), t0s = 1; end
%         t0s = [t0s, T+1];
%         for s = 1:numel(t0s)-1
%             t0 = t0s(s); t1 = t0s(s+1)-1;  if t0>=t1, continue; end
%             % first 0->1 transition (unmarried -> married) within [t0+1, t1]
%             trans = find(~M(i,t0:t1-1) & M(i,t0+1:t1), 1, 'first');
%             if ~isempty(trans)
%                 age_mar_years = AgeM(i, t0+trans)/12;
%                 if isfinite(age_mar_years) && age_mar_years>=0
%                     if fem(i), ages_w(end+1,1) = age_mar_years; %#ok<AGROW>
%                     else       ages_m(end+1,1) = age_mar_years; %#ok<AGROW>
%                     end
%                 end
%             end
%         end
%     end
%     if ~isempty(ages_w), moments.age_at_marriage.women_years = median(ages_w,'omitnan'); end
%     if ~isempty(ages_m), moments.age_at_marriage.men_years   = median(ages_m,'omitnan'); end
%     both_naive = [ages_w; ages_m];
%     if ~isempty(both_naive), moments.age_at_marriage.all_years = median(both_naive,'omitnan'); end
% 
%     % ------- (ii-a) KM medians on AGE (years), with LEFT TRUNCATION -------------
%     % Build left-truncated spells on AGE scale
%     [entryY, exitY, status, sexID] = buildMarriageSpells_AgeLT(M, NEWB, AgeM, fem);
% 
%     % Women
%     isW = (sexID==1);
%     [tw, Sw] = KM_leftTrunc_Age(entryY(isW), exitY(isW), status(isW));
%     moments.age_at_marriage_km.women_median_years = kmMedianFromCurve(tw, Sw);
% 
%     % Men
%     isM = (sexID==0);
%     [tm, Sm] = KM_leftTrunc_Age(entryY(isM), exitY(isM), status(isM));
%     moments.age_at_marriage_km.men_median_years = kmMedianFromCurve(tm, Sm);
% 
%     % Pooled
%     [ta, Sa] = KM_leftTrunc_Age(entryY, exitY, status);
%     moments.age_at_marriage_km.all_median_years = kmMedianFromCurve(ta, Sa);
% 
%     % --- (ii-b) Exponential MLE means on AGE (years), with LEFT TRUNCATION ------
%     % Hazard per YEAR with left truncation: lambda = (#events) / sum(exposure years)
%     expMLE_age = @(a,b,s) sum(s==1) / max(eps, sum(b - a));   % a=entryY, b=exitY, s∈{0,1}
%     lam_w = expMLE_age(entryY(isW), exitY(isW), status(isW));
%     lam_m = expMLE_age(entryY(isM), exitY(isM), status(isM));
%     lam_a = expMLE_age(entryY,       exitY,      status);
% 
%     moments.age_at_marriage_expMLE.women_mean_years = mean(entryY(isW),'omitnan') + 1/max(eps,lam_w);
%     moments.age_at_marriage_expMLE.men_mean_years   = mean(entryY(isM),'omitnan') + 1/max(eps,lam_m);
%     moments.age_at_marriage_expMLE.all_mean_years   = mean(entryY,'omitnan')      + 1/max(eps,lam_a);
% end

% %% ========== 8) Age at (first) marriage: full-path (unbiased) vs kept-window ==========
% % Existing quick medians (within-window, events only) can stay as you had them earlier
% % -- we keep them for comparison in `moments.age_at_marriage.*` (years).
% 
% % -------------------------- FULL-PATH (preferred) --------------------------
% use_full = exist('agentPanel_full','var') && ~isempty(agentPanel_full) ...
%            && isfield(agentPanel_full,'married') && isfield(agentPanel_full,'age_chron');
% 
% if use_full
%     APf   = agentPanel_full;
%     M_f   = APf.married ~= 0;          % N x T_total
%     NEWBf = APf.newborn ~= 0;          % N x T_total (1 at (re)birth/spell start)
%     AgeMf = APf.age_chron;             % N x T_total (months, chronological)
%     FEMf  = APf.female  ~= 0;          % N x 1 (logical)
% 
%     % Build "from-(re)birth" at-risk spells on the AGE scale (years)
%     [entryY, exitY, status, sexID] = buildMarriageSpells_Full(M_f, NEWBf, AgeMf, FEMf);
% 
%     % --- KM medians on AGE (years), no left truncation needed on the full path
%     % Women
%     isW = (sexID==1);
%     [tw, Sw] = KM_noToolbox_onAge(exitY(isW), status(isW));
%     moments_full.age_at_marriage_km.women_median_years = kmMedianFromCurve(tw, Sw);
% 
%     % Men
%     isM = (sexID==0);
%     [tm, Sm] = KM_noToolbox_onAge(exitY(isM), status(isM));
%     moments_full.age_at_marriage_km.men_median_years = kmMedianFromCurve(tm, Sm);
% 
%     % Pooled
%     [ta, Sa] = KM_noToolbox_onAge(exitY, status);
%     moments_full.age_at_marriage_km.all_median_years = kmMedianFromCurve(ta, Sa);
% 
%     % --- Exponential-MLE mean AGE (years) on full path
%     % Hazard per YEAR: lambda = (#events) / (total exposure years)
%     hazardY = @(a,b,s) sum(s==1) / max(eps, sum(b - a));   % a=entryY, b=exitY, s∈{0,1}
% 
%     lam_w = hazardY(entryY(isW), exitY(isW), status(isW));
%     lam_m = hazardY(entryY(isM), exitY(isM), status(isM));
%     lam_a = hazardY(entryY,       exitY,      status);
% 
%     moments_full.age_at_marriage_expMLE.women_mean_years = mean(entryY(isW),'omitnan') + 1/max(eps,lam_w);
%     moments_full.age_at_marriage_expMLE.men_mean_years   = mean(entryY(isM),'omitnan') + 1/max(eps,lam_m);
%     moments_full.age_at_marriage_expMLE.all_mean_years   = mean(entryY,'omitnan')      + 1/max(eps,lam_a);
% else
%     % ------------------ Fallback: kept-window, left-truncated KM ------------------
%     % Uses your already-built kept-window variables: agentPanel.*, NEWB (or detect),
%     % AgeM (prefer age_chron), and the left-trunc helpers you already have.
%    moments.age_at_marriage.women_years = NaN;
% moments.age_at_marriage.men_years   = NaN;
% moments.age_at_marriage.all_years   = NaN;
% 
% % ---- Inputs (prefer chronological age; fallback if needed) ----------------------
% if isfield(agentPanel,'age_chron') && ~isempty(agentPanel.age_chron)
%     AgeM = agentPanel.age_chron;          % months, chronological
% elseif isfield(agentPanel,'age') && ~isempty(agentPanel.age)
%     % Fallback: treat 'age' as chronological if no explicit 'age_chron'
%     if exist('param','var') && isfield(param,'entry_age_years')
%         AgeM = agentPanel.age + param.entry_age_years*12;  % best-effort
%     else
%         AgeM = agentPanel.age;                              % last resort
%     end
% else
%     warning('compute_moments:age','No age field found; skipping age-at-marriage moments.');
%     goto_end_age_section = true; %#ok<NASGU>
% end
% 
% if exist('goto_end_age_section','var')
%     % nothing to do
% else
%     if isfield(agentPanel,'newborn') && ~isempty(agentPanel.newborn)
%         NEWB = agentPanel.newborn ~= 0;
%     else
%         % Detect (re)births by age resets
%         NEWB = false(size(AgeM));
%         NEWB(:,1) = true;
%         NEWB(:,2:end) = AgeM(:,2:end) < AgeM(:,1:end-1);
%     end
%     M   = agentPanel.married ~= 0;
%     fem = agentPanel.female  ~= 0;
%     [N,T] = size(AgeM);
% 
%     % --------------------- (i) Naive within-window medians -----------------------
%     ages_w = []; ages_m = [];   % realized ages (years) for spells with an event
%     for i = 1:N
%         t0s = find(NEWB(i,:)); if isempty(t0s), t0s = 1; end
%         t0s = [t0s, T+1];
%         for s = 1:numel(t0s)-1
%             t0 = t0s(s); t1 = t0s(s+1)-1;  if t0>=t1, continue; end
%             % first 0->1 transition (unmarried -> married) within [t0+1, t1]
%             trans = find(~M(i,t0:t1-1) & M(i,t0+1:t1), 1, 'first');
%             if ~isempty(trans)
%                 age_mar_years = AgeM(i, t0+trans)/12;
%                 if isfinite(age_mar_years) && age_mar_years>=0
%                     if fem(i), ages_w(end+1,1) = age_mar_years; %#ok<AGROW>
%                     else       ages_m(end+1,1) = age_mar_years; %#ok<AGROW>
%                     end
%                 end
%             end
%         end
%     end
%     if ~isempty(ages_w), moments.age_at_marriage.women_years = median(ages_w,'omitnan'); end
%     if ~isempty(ages_m), moments.age_at_marriage.men_years   = median(ages_m,'omitnan'); end
%     both_naive = [ages_w; ages_m];
%     if ~isempty(both_naive), moments.age_at_marriage.all_years = median(both_naive,'omitnan'); end
% 
%         % Build left-truncated spells on AGE (years)
%         [entryY, exitY, status, sexID] = buildMarriageSpells_AgeLT(M, NEWB, AgeM, FEM);
% 
%         % KM (left-trunc.) medians on AGE
%         isW = (sexID==1); isM = ~isW;
%         [tw, Sw] = KM_leftTrunc_Age(entryY(isW), exitY(isW), status(isW));
%         [tm, Sm] = KM_leftTrunc_Age(entryY(isM), exitY(isM), status(isM));
%         [ta, Sa] = KM_leftTrunc_Age(entryY,       exitY,      status);
%         moments.age_at_marriage_km.women_median_years = kmMedianFromCurve(tw, Sw);
%         moments.age_at_marriage_km.men_median_years   = kmMedianFromCurve(tm, Sm);
%         moments.age_at_marriage_km.all_median_years   = kmMedianFromCurve(ta, Sa);
% 
%         % Exp-MLE mean AGE with left truncation
%         hazardY = @(a,b,s) sum(s==1) / max(eps, sum(b - a));
%         lam_w = hazardY(entryY(isW), exitY(isW), status(isW));
%         lam_m = hazardY(entryY(isM), exitY(isM), status(isM));
%         lam_a = hazardY(entryY,       exitY,      status);
%         moments.age_at_marriage_expMLE.women_mean_years = mean(entryY(isW),'omitnan') + 1/max(eps,lam_w);
%         moments.age_at_marriage_expMLE.men_mean_years   = mean(entryY(isM),'omitnan') + 1/max(eps,lam_m);
%         moments.age_at_marriage_expMLE.all_mean_years   = mean(entryY,'omitnan')      + 1/max(eps,lam_a);
%     end
% end


%% ====== 9) Unemployment Benefits (overall and by pre-separation income) ======

% Keep your original “mask-at-source” variables
B_pos  = B;   B_pos(~U)  = NaN;
pB_pos = pB;  pB_pos(~pU) = NaN;

% Also define wages with NaNs when not employed
W_pos  = W;   W_pos(~E)  = NaN;
pW_pos = pW;  pW_pos(~pE) = NaN;

% --- Overall means (report both own, partner, and pooled) ---
moments.unemp_benefits.all_individuals = mean(B_pos(:),  'omitnan');
moments.unemp_benefits.all_partners    = mean(pB_pos(:), 'omitnan');
moments.unemp_benefits.all             = mean([B_pos(:); pB_pos(:)], 'omitnan');  % pooled

R_lin  = W(:);
y      = R_lin(valid);
med_inc = median(y, 'omitnan');
low  = (y <= med_inc);
high = (y >  med_inc);

moments.unemp_benefits.median_income              = med_inc;
moments.unemp_benefits.lowinc_mean_share          = moments.unemp_benefits.all/mean(y(low),  'omitnan');
moments.unemp_benefits.highinc_mean_share         = moments.unemp_benefits.all/mean(y(high), 'omitnan');
moments.unemp_benefits.all_share_of_median_income = moments.unemp_benefits.all/med_inc;

end % ===== end compute_moments =====


% ========================== Helpers (local functions) ============================
% ============================ Helpers (add once) ===========================
function [entryY, exitY, status, sexID] = buildMarriageSpells_Full(M, NEWB, AgeM, FEM)
% From-(re)birth spells over the FULL path; AGE scale (years).
% Each NEWB==1 starts a new at-risk spell. If M==1 at entry, skip that spell.
[N,T] = size(M);
entryY=[]; exitY=[]; status=[]; sexID=[];
for i=1:N
    t0s = find(NEWB(i,:)); if isempty(t0s), t0s = 1; end
    t0s = sort(unique([t0s, T+1]));      % sentinel
    for s=1:numel(t0s)-1
        t0 = t0s(s); t1 = t0s(s+1)-1; if t0>=t1, continue; end
        if M(i,t0), continue; end        % already married at entry -> not at risk of first marriage
        trans = find(~M(i,t0:t1-1) & M(i,t0+1:t1), 1, 'first');   % first 0->1
        eY = AgeM(i,t0)/12;
        if ~isempty(trans), xY = AgeM(i,t0+trans)/12; st=1;
        else,               xY = AgeM(i,t1)/12;       st=0;
        end
        entryY(end+1,1)=eY; exitY(end+1,1)=xY; status(end+1,1)=st; sexID(end+1,1)=FEM(i);
    end
end
end

function [t,S] = KM_noToolbox_onAge(eventAgeY, status)
% Product-limit survivor on AGE; eventAgeY (years) at exit; status=1 if event, 0 if censored.
eventAgeY = eventAgeY(:); status = status(:)~=0;
evt = unique(eventAgeY(status));
t=[]; S=[]; if isempty(evt), return; end
S_now = 1;
for k=1:numel(evt)
    tau = evt(k);
    n_risk = sum(eventAgeY >= tau);             % at risk just before tau
    d_tau  = sum(eventAgeY==tau & status);
    if n_risk>0 && d_tau>0
        S_now = S_now * (1 - d_tau/n_risk);
        t(end+1,1)=tau; S(end+1,1)=S_now; %#ok<AGROW>
    end
end
end

function [entryY, exitY, status, sexID] = buildMarriageSpells_AgeLT(M, NEWB, AgeM, FEM)
% Left-truncated spells within the KEPT window; AGE scale (years).
[N,T] = size(M);
entryY=[]; exitY=[]; status=[]; sexID=[];
for i=1:N
    t0s = find(NEWB(i,:));
    if isempty(t0s)
        if ~M(i,1), t0s = 1; else, t0s = []; end
    else
        if ~M(i,1) && t0s(1) > 1, t0s = [1, t0s]; end
    end
    if isempty(t0s), continue; end
    t0s = sort(unique([t0s, T+1]));
    for s=1:numel(t0s)-1
        t0 = t0s(s); t1 = t0s(s+1)-1; if t0>=t1, continue; end
        if M(i,t0), continue; end
        trans = find(~M(i,t0:t1-1) & M(i,t0+1:t1), 1, 'first');
        eY = AgeM(i,t0)/12;
        if ~isempty(trans), xY = AgeM(i,t0+trans)/12; st=1;
        else,               xY = AgeM(i,t1)/12;       st=0;
        end
        entryY(end+1,1)=eY; exitY(end+1,1)=xY; status(end+1,1)=st; sexID(end+1,1)=FEM(i);
    end
end
end

function [t,S] = KM_leftTrunc_Age(entryY, exitY, status)
% Left-truncated Kaplan–Meier on AGE (years).
entryY = entryY(:); exitY = exitY(:); status = status(:)~=0;
ok = isfinite(entryY) & isfinite(exitY) & (exitY >= entryY);
entryY = entryY(ok); exitY = exitY(ok); status = status(ok);
evt = unique(exitY(status));
t=[]; S=[]; if isempty(evt), return; end
S_now = 1;
for k=1:numel(evt)
    tau = evt(k);
    n_risk = sum( (entryY < tau) & (exitY >= tau) );
    d_tau  = sum( (exitY == tau) & status );
    if n_risk>0 && d_tau>0
        S_now = S_now * (1 - d_tau/n_risk);
        t(end+1,1)=tau; S(end+1,1)=S_now; %#ok<AGROW>
    end
end
end

function med = kmMedianFromCurve(t,S)
if isempty(t) || isempty(S), med = NaN; return; end
k = find(S <= 0.5, 1, 'first');
if isempty(k), med = NaN; else, med = t(k); end
end

function tr = transitions_by_gender(gmask, U, E, O, T)
    G = repmat(gmask,1,T-1);
    prevU = U(:,1:end-1) & G;    nextE = E(:,2:end) & G;
    prevE = E(:,1:end-1) & G;    nextU = U(:,2:end) & G;
    prevO = O(:,1:end-1) & G;    nextO = O(:,2:end) & G;

    tr.UE = deal_rate(prevU, nextE);
    tr.EU = deal_rate(prevE, nextU);
    tr.EO = deal_rate(prevE, nextO);
    tr.UO = deal_rate(prevU, nextO);
    tr.OE = deal_rate(prevO, nextE);
end

function r = deal_rate(prev, next)
    SMALL = 1e-12;
    flows   = prev & next;
    at_risk = prev;
    r = sum(flows,'all') / max(sum(at_risk,'all'), SMALL);
end
