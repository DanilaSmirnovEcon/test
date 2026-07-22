function PlotSimul(agentPanel, param,valueFunc, opts)
% PlotSimul  Recreate model graphs using *simulated panel* data.
% Usage: PlotSimul(agentPanel, param)
%        PlotSimul(agentPanel, param, struct('outdir','Plots','nBinsAbility',25))
%
% Inputs
%   agentPanel : struct from SimulatePanel (burn-in already trimmed)
%   param      : struct with tax handles (nokid_tax, onekid_tax), optional periods_per_year
%   opts       : struct (optional)
%       .outdir         default 'Plots'
%       .nBinsAbility   default 25
%       .pp_change_hzn  default param.periods_per_year (or 12)

% ------------------- options & setup -------------------
if nargin<4, opts = struct(); end
if ~isfield(opts,'outdir'),        opts.outdir = 'Plots'; end
if ~isfield(opts,'nBinsAbility'),  opts.nBinsAbility = 25; end
if ~isfield(opts,'nBinsWages'),  opts.nBinsWages = 25; end

ppy = 12;
if isfield(param,'periods_per_year') && ~isempty(param.periods_per_year)
    ppy = param.periods_per_year;
end
if ~isfield(opts,'pp_change_hzn') || isempty(opts.pp_change_hzn)
    opts.pp_change_hzn = ppy;
end

if ~exist(opts.outdir, 'dir'), mkdir(opts.outdir); end
SMALL = 1e-10;

% ------------------- pull panel data -------------------
[N,T] = size(agentPanel.employed);

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
pA   = agentPanel.partner_ability;
pW  = exp(agentPanel.partner_ability) .* agentPanel.partner_wage; % partner
pW  = pW.*param.fkidwpenalty.*(1-agentPanel.female) + pW.*(agentPanel.female);
WO  = agentPanel.wage;     % offer wage (self)
pWO = agentPanel.partner_wage;

% Real wage (only meaningful when employed)
R = W;           % NT
R(~E) = NaN;
pR = pW;           % NT
pR(~pE) = NaN;

% Ability bins (same for all figs)
Aall = A(:); Aall = Aall(isfinite(Aall));
if isempty(Aall), error('No valid ability data.'); end
edges = linspace(min(Aall), max(Aall), opts.nBinsAbility+1);
centers = 0.5*(edges(1:end-1)+edges(2:end));

Wall = W(:); Wall = Wall(isfinite(Wall));
if isempty(Wall), error('No valid wage data.'); end
edgesW = linspace(min(Wall), max(Wall), opts.nBinsWages+1);
centersW = 0.5*(edgesW(1:end-1)+edgesW(2:end));

% bin for own ability (already had this)
bin  = @(X,mask,agg,varargin) bin_by_ability(A,  X, mask, edges, agg, varargin{:});   % accepts denom mask
% bin for partner ability
pbin = @(X,mask,agg,varargin) bin_by_ability(pA, X, mask, edges, agg, varargin{:});   % <-- partner/spouse

% =====================================================================
% 1) PRODUCTIVITY DISTRIBUTIONS (singles vs married incl. spouses)
% =====================================================================
dens_single        = bin([], ~M, 'share');         % singles over own ability
dens_married_own   = bin([],  M, 'share');         % married over own ability
dens_married_spous = pbin([], M, 'share');         % married over partner ability

% Combine to represent "married incl. spouses"
% (average keeps total area = 1 because each side integrates to 1)
dens_married_both = 0.5 * (dens_married_own + dens_married_spous);

fig = figure('Color','w');
plot(centers, dens_married_both, 'LineWidth',2); hold on;
plot(centers, dens_single,       'LineWidth',2);
formatPlot('Productivity', 'Density', 'Distribution over Productivities');
legend({'Married', 'Singles'}, 'Location','best','Interpreter','none');
saveas(fig, fullfile(opts.outdir,'EndogProductivities.png'));

% gender-specific (married, incl. spouses)
% own ability distributions
dens_f_mar_own = bin([], M & gF, 'share');
dens_m_mar_own = bin([], M & gM, 'share');
% spouse ability distributions
dens_f_mar_spouse = pbin([], M & pgF, 'share');
dens_m_mar_spouse = pbin([], M & pgM, 'share');
% combine own + spouse
dens_f_mar_both = 0.5 * (dens_f_mar_own + dens_f_mar_spouse);
dens_m_mar_both = 0.5 * (dens_m_mar_own + dens_m_mar_spouse);
% plot
fig = figure('Color','w');
plot(centers, dens_f_mar_both, 'LineWidth',2); hold on;
plot(centers, dens_m_mar_both, 'LineWidth',2);
formatPlot('Productivity','Density','Married: Distribution over Productivities');
legend({'Female','Male'}, 'Location','best','Interpreter','none');
saveas(fig, fullfile(opts.outdir,'EndogProductivitiesGenderMarried.png'));

% gender-specific (single)
dens_f_sing = bin([], ~M & gF, 'share');
dens_m_sing = bin([], ~M & gM, 'share');
fig = figure('Color','w');
plot(centers, dens_f_sing,'LineWidth',2); hold on;
plot(centers, dens_m_sing,'LineWidth',2);
formatPlot('Productivity','Density','Single: Distribution over Productivities');
legend({'Female','Male'}, 'Location','best','Interpreter','none');
saveas(fig, fullfile(opts.outdir,'EndogProductivitiesGenderSingle.png'));

% =====================================================================
% 2) WAGE MARITAL PREMIUMS BY ABILITY (female, male, combined)
%    Including spouses (own + spouse, averaged)
% =====================================================================

% --- FEMALE ---
% Married females: own wages by own ability
mean_wage_f_M_own   = bin(R,  E & M  & gF, 'mean');
% Wages of spouses who are female, placed into the partner's ability bins
mean_wage_f_M_spous = pbin(pR, pE & M  & pgF, 'mean');
% Married (incl. spouses): average keeps scale comparable to singles
mean_wage_f_M_both  = 0.5 * (mean_wage_f_M_own + mean_wage_f_M_spous);
% Singles (no spouse contribution)
mean_wage_f_S = bin(R,  E & ~M & gF, 'mean');

% --- MALE ---
mean_wage_m_M_own   = bin(R,  E & M  & gM, 'mean');
mean_wage_m_M_spous = pbin(pR, pE & M  & pgM, 'mean');
mean_wage_m_M_both  = 0.5 * (mean_wage_m_M_own + mean_wage_m_M_spous);
mean_wage_m_S = bin(R,  E & ~M & gM, 'mean');

% --- COMBINED (ability-bin weighted by employed mass), incl. spouses ---
% Masses for weighting (own + spouse, averaged)
mass_f_M_own   = bin([], E & M  & gF, 'mass');
mass_f_M_spous = pbin([], pE & M  & pgF, 'mass');
mass_f_M_both  = 0.5 * (mass_f_M_own + mass_f_M_spous);

mass_m_M_own   = bin([], E & M  & gM, 'mass');
mass_m_M_spous = pbin([], pE & M  & pgM, 'mass');
mass_m_M_both  = 0.5 * (mass_m_M_own + mass_m_M_spous);

% Singles masses unchanged
mass_f_S = bin([], E & ~M & gF, 'mass');
mass_m_S = bin([], E & ~M & gM, 'mass');

% Married means incl. spouses (already computed above): mean_wage_*_M_both
% Singles means: mean_wage_*_S

M_all = (mean_wage_m_M_both.*mass_m_M_both + mean_wage_f_M_both.*mass_f_M_both) ./ ...
        max(mass_m_M_both + mass_f_M_both, SMALL);

S_all = (mean_wage_m_S.*mass_m_S + mean_wage_f_S.*mass_f_S) ./ ...
        max(mass_m_S + mass_f_S, SMALL);

% Female premiums
prem_f = 100 * (mean_wage_f_M_both - mean_wage_f_S) ./ max(mean_wage_f_S, SMALL);
prem_f_smooth = smoothdata(prem_f, 'movmean', 3);

% Male premiums
prem_m = 100 * (mean_wage_m_M_both - mean_wage_m_S) ./ max(mean_wage_m_S, SMALL);
prem_m_smooth = smoothdata(prem_m, 'movmean', 3);

% Plot female + male together
fig = figure('Color','w');
plot(centers, prem_f_smooth, 'LineWidth',2); hold on;
plot(centers, prem_m_smooth, 'LineWidth',2);
formatPlot('Productivity','% Wage difference','Wage Marital Premiums (M−S)');
legend({'Female','Male'}, 'Location','best','Interpreter','none');
saveas(fig, fullfile(opts.outdir,'EndogWageDiff_FvsM.png'));

% Combined
prem_all = 100 * (M_all - S_all) ./ max(S_all, SMALL);
prem_all_smooth = smoothdata(prem_all, 'movmean', 3);

fig = figure('Color','w');
plot(centers, prem_all_smooth, 'LineWidth',2);
formatPlot('Productivity','% Wage difference','Average Wage Marital Premium (M−S)');
saveas(fig, fullfile(opts.outdir,'EndogWageDiff_All.png'));

% =====================================================================
% 2b) AVERAGE RAW (GRID) WAGE BY ABILITY  �  E[w | a-bin], employed only
% =====================================================================
W_pos = W;                 % use the stored grid wage, not exp(a).*w
W_pos(~E) = NaN;            % only meaningful when employed

avg_w_all = bin(W_pos, E, 'mean');  % mean(w) in each ability bin among employed

fig = figure('Color','w');
plot(centers, avg_w_all, 'LineWidth', 2);
formatPlot('Productivity','Average grid wage w','Average wage by ability (raw w)');
saveas(fig, fullfile(opts.outdir,'AvgGridWage_by_Ability.png'));

% (optional) by gender:
avg_w_f = bin(W_pos, E & gF, 'mean');
avg_w_m = bin(W_pos, E & gM, 'mean');
fig = figure('Color','w');
plot(centers, avg_w_f, 'LineWidth', 2); hold on;
plot(centers, avg_w_m, 'LineWidth', 2);
formatPlot('Productivity','Average grid wage w','Average wage by ability (raw w) � by gender');
legend({'Female','Male'}, 'Location','best','Interpreter','none');
saveas(fig, fullfile(opts.outdir,'AvgGridWage_by_Ability_Gender.png'));


% =====================================================================
% 2.5) WAGES & UNEMPLOYMENT BY ABILITY (levels)
%      - Wages: mean real wage among employed in each ability bin
%      - Unemployment: share unemployed within the labor force (E|U)
% =====================================================================

% --- Average wage by ability (employed only) ---
avg_w_all = bin(R,  E,                      'mean');
avg_w_f   = bin(R,  E & gF,      'mean');
avg_w_m   = bin(R,  E & gM,      'mean');

avg_pw_all = bin(pR,  pE,                      'mean');
avg_pw_f   = bin(pR,  pE & pgF,      'mean');
avg_pw_m   = bin(pR,  pE & pgM,      'mean');

avg_w_all  = 0.5.*(avg_w_all+avg_pw_all);
avg_w_f    = 0.5.*(avg_w_f+avg_pw_f);
avg_w_m    = 0.5.*(avg_w_m+avg_pw_m);

fig = figure('Color','w'); plot(centers, avg_w_all,'LineWidth',2);
formatPlot('Productivity','Average wage (real)','Average Wage by Ability (employed only)');
saveas(fig, fullfile(opts.outdir,'AvgWage_byAbility_All.png'));

fig = figure('Color','w'); plot(centers, avg_w_f,'LineWidth',2);
formatPlot('Productivity','Average wage (real)','Average Wage by Ability � Women (employed only)');
saveas(fig, fullfile(opts.outdir,'AvgWage_byAbility_Female.png'));

fig = figure('Color','w'); plot(centers, avg_w_m,'LineWidth',2);
formatPlot('Productivity','Average wage (real)','Average Wage by Ability � Men (employed only)');
saveas(fig, fullfile(opts.outdir,'AvgWage_byAbility_Male.png'));

% --- Unemployment rate by ability (within the labor force) ---
inLF    = E | U;  % labor force indicator
ur_all  = 100 * bin([], U & inLF,                        'share_in_group', inLF);
ur_f    = 100 * bin([], U & inLF & gF,        'share_in_group', inLF & gF);
ur_m    = 100 * bin([], U & inLF & gM,        'share_in_group', inLF & gM);

pinLF    = pE | pU;  % labor force indicator
pur_all  = 100 * pbin([], pU & pinLF,                        'share_in_group', pinLF);
pur_f    = 100 * pbin([], pU & pinLF & pgF,        'share_in_group', pinLF & pgF);
pur_m    = 100 * pbin([], pU & pinLF & pgM,        'share_in_group', pinLF & pgM);

ur_all  = 0.5.*(ur_all+pur_all);
ur_f    = 0.5.*(ur_f+pur_f);
ur_m    = 0.5.*(ur_m+pur_m);

fig = figure('Color','w'); plot(centers, ur_all,'LineWidth',2);
formatPlot('Productivity','Unemployment rate (%, LF)','Unemployment by Ability (all)');
saveas(fig, fullfile(opts.outdir,'UnempRate_byAbility_All.png'));

fig = figure('Color','w'); plot(centers, ur_f,'LineWidth',2);
formatPlot('Productivity','Unemployment rate (%, LF)','Unemployment by Ability � Women');
saveas(fig, fullfile(opts.outdir,'UnempRate_byAbility_Female.png'));

fig = figure('Color','w'); plot(centers, ur_m,'LineWidth',2);
formatPlot('Productivity','Unemployment rate (%, LF)','Unemployment by Ability � Men');
saveas(fig, fullfile(opts.outdir,'UnempRate_byAbility_Male.png'));

% =====================================================================
% 3) UNEMPLOYMENT MARITAL GAPS BY ABILITY (female, male, combined)
% =====================================================================
% ur_f_M = bin([], U & M  & gF, 'share_in_group', M  & gF);
% ur_f_S = bin([], U & ~M & gF, 'share_in_group', ~M & gF);
% ur_m_M = bin([], U & M  & gM, 'share_in_group', M  & gM);
% ur_m_S = bin([], U & ~M & gM, 'share_in_group', ~M & gM);
% 
% gap_f = -100*(ur_f_M - ur_f_S);   % (S - M)
% gap_m = -100*(ur_m_M - ur_m_S);
% 
% fig = figure('Color','w'); plot(centers, gap_f,'LineWidth',2);
% formatPlot('Productivity','pp','Female Unemployment Marital Gap (S-M)');
% saveas(fig, fullfile(opts.outdir,'EndogFemaleUnemploymentDiff.png'));
% 
% fig = figure('Color','w'); plot(centers, gap_m,'LineWidth',2);
% formatPlot('Productivity','pp','Male Unemployment Marital Gap (S-M)');
% saveas(fig, fullfile(opts.outdir,'EndogMaleUnemploymentDiff.png'));
% 
% % combined (weighted across genders by group mass in bin)
% mass_f_all = bin([], (M|~M) & gF, 'mass');   % all female obs by bin
% mass_m_all = bin([], (M|~M) & gM, 'mass');
% ur_M_all = (ur_m_M.*mass_m_all + ur_f_M.*mass_f_all) ./ max(mass_m_all+mass_f_all,SMALL);
% ur_S_all = (ur_m_S.*mass_m_all + ur_f_S.*mass_f_all) ./ max(mass_m_all+mass_f_all,SMALL);
% gap_all = -100*(ur_M_all - ur_S_all);
% 
% fig = figure('Color','w'); plot(centers, gap_all,'LineWidth',2);
% formatPlot('Productivity','pp','Unemployment Marital Gap (S-M)');
% saveas(fig, fullfile(opts.outdir,'EndogUnemploymentDiff.png'));

% =====================================================================
% 3) UNEMPLOYMENT MARITAL GAPS BY ABILITY (female, male, combined)
%     Partner-inclusive (own + spouse), by ability bin
% =====================================================================

% Labor-force / unemployment masks (focal)
LF = ~O;             % in labor force
U  = LF & ~E;        % unemployed

% Partner-side masks (use pO/pE if available; else best-effort fallback)
if exist('pO','var') && exist('pE','var')
    pLF = ~pO;
    pU  = pLF & ~pE;
elseif exist('pE','var')
    % If no partner OLF flag, we can only mark employed partners; treat others as non-unemployed
    pLF = pE;        % crude fallback
    pU  = false(size(pE));
else
    error('Need partner-side employment/OLF masks (pE and ideally pO) to include partner unemployment.');
end

% -----------------------------
% FEMALE
% -----------------------------
% Married females (own + partner is female)
num_f_M_A  = bin([], U & M & gF,         'mass');
den_f_M_A  = bin([],     M & gF,         'mass');
num_f_M_pA = pbin([], pU & M & pgF,      'mass');
den_f_M_pA = pbin([],     M & pgF,       'mass');

ur_f_M_both = (num_f_M_A + num_f_M_pA) ./ max(den_f_M_A + den_f_M_pA, SMALL);

% Single females (no partner contribution)
num_f_S = bin([], U & ~M & gF, 'mass');
den_f_S = bin([],     ~M & gF, 'mass');
ur_f_S  = num_f_S ./ max(den_f_S, SMALL);

gap_f = -100*(ur_f_M_both - ur_f_S);   % (S - M), in percentage points

fig = figure('Color','w'); plot(centers, gap_f,'LineWidth',2);
formatPlot('Productivity','pp','Female Unemployment Marital Gap (incl. spouse; S−M)');
saveas(fig, fullfile(opts.outdir,'EndogFemaleUnemploymentDiff.png'));

% -----------------------------
% MALE
% -----------------------------
num_m_M_A  = bin([], U & M & gM,         'mass');
den_m_M_A  = bin([],     M & gM,         'mass');
num_m_M_pA = pbin([], pU & M & pgM,      'mass');
den_m_M_pA = pbin([],     M & pgM,       'mass');

ur_m_M_both = (num_m_M_A + num_m_M_pA) ./ max(den_m_M_A + den_m_M_pA, SMALL);

num_m_S = bin([], U & ~M & gM, 'mass');
den_m_S = bin([],     ~M & gM, 'mass');
ur_m_S  = num_m_S ./ max(den_m_S, SMALL);

gap_m = -100*(ur_m_M_both - ur_m_S);   % (S - M)

fig = figure('Color','w'); plot(centers, gap_m,'LineWidth',2);
formatPlot('Productivity','pp','Male Unemployment Marital Gap (incl. spouse; S−M)');
saveas(fig, fullfile(opts.outdir,'EndogMaleUnemploymentDiff.png'));

% -----------------------------
% COMBINED (composition-weighted across genders, by bin)
% Match your previous weighting approach but with partner-inclusive masses
% -----------------------------
mass_f_all_both = bin([], gF, 'mass') + pbin([], pgF, 'mass');
mass_m_all_both = bin([], gM, 'mass') + pbin([], pgM, 'mass');

ur_M_all = (ur_m_M_both .* mass_m_all_both + ur_f_M_both .* mass_f_all_both) ./ ...
           max(mass_m_all_both + mass_f_all_both, SMALL);

ur_S_all = (ur_m_S .* mass_m_all_both + ur_f_S .* mass_f_all_both) ./ ...
           max(mass_m_all_both + mass_f_all_both, SMALL);

gap_all = -100*(ur_M_all - ur_S_all);  % (S - M)

fig = figure('Color','w'); plot(centers, gap_all,'LineWidth',2);
formatPlot('Productivity','pp','Unemployment Marital Gap (incl. spouse; S−M)');
saveas(fig, fullfile(opts.outdir,'EndogUnemploymentDiff.png'));

% =====================================================================
% 4) MARRIAGE RATES BY ABILITY (count own ability + spouse ability)
% =====================================================================

% Denominators: total mass in each ability bin (own and spouse ability)
mass_all_A  = bin([], true(size(A)),  'mass');   % everyone by own ability
mass_all_pA = pbin([], true(size(pA)), 'mass');  % everyone by partner ability

% Numerators: married mass in each ability bin (own and spouse ability)
mass_M_A    = bin([], M,  'mass');               % married, binned by own ability
mass_M_pA   = pbin([], M, 'mass');               % married, binned by spouse ability

% Combine by summing numerators and denominators, then take the share
mar_rate_both = 100 * (mass_M_A + mass_M_pA) ./ max(mass_all_A + mass_all_pA, SMALL);

% Plot
fig = figure('Color','w');
plot(centers, mar_rate_both, 'LineWidth', 2);
formatPlot('Ability','% Married','Share Married by Ability (own + spouse ability)');
saveas(fig, fullfile(opts.outdir,'EndogMarriedShare.png'));

% % Group masks (focal)
% LF = ~O;          % in labor force
% U  = LF & ~E;     % unemployed (in LF but not employed)
% 
% % Partner-side masks (use what you have)
% if exist('pO','var')
%     pLF = ~pO;
%     pU  = pLF & ~pE;
% else
%     % Fallbacks if pO not available (best-effort)
%     % If you only have pE, treat partners in LF ≈ pE (crude)
%     pLF = pE;
%     pU  = false(size(pE));  % cannot distinguish unemployment without pO/pU
% end
% 
% % =========================
% % Married share among the LF (by ability bin), counting partners
% % Numerator: married in group (own + partner bins)
% % Denominator: group mass (own + partner bins)
% % =========================
% num_A_LF  = bin([], M & LF,  'mass');
% den_A_LF  = bin([], LF,      'mass');
% 
% num_pA_LF = pbin([], M & pLF, 'mass');
% den_pA_LF = pbin([], pLF,     'mass');
% 
% mar_rate_LF_both = 100 * (num_A_LF + num_pA_LF) ./ max(den_A_LF + den_pA_LF, SMALL);
% 
% fig = figure('Color','w');
% plot(centers, mar_rate_LF_both, 'LineWidth',2);
% formatPlot('Ability','% Married','Married Share by Ability (Labor Force; own + spouse)');
% saveas(fig, fullfile(opts.outdir,'EndogMarriedShare_LF_partners.png'));

% =====================================================================
% 5) WAGE CHANGE AFTER UNEMPLOYMENT (by initial ability) — with partners
%     Plot 1: Singles vs Married   (smoothed only)
%     Plot 2: With kids vs Without kids   (smoothed only)
% =====================================================================

chg   = [];        % event-level % change in wage after unemployment
abil0 = [];        % initial OWN ability at separation
pabil0= [];        % initial PARTNER ability at separation (NaN if none)
isMar = [];        % marital status at separation (0/1)
hasKid= [];        % child indicator at separation (0/1)

for i = 1:N
    sep_t = find(E(i,1:end-1) & ~E(i,2:end));  % separation between t and t+1
    for k = 1:numel(sep_t)
        t  = sep_t(k);
        w0 = R(i,t);
        if ~isfinite(w0) || w0 <= 0, continue; end

        % first re-employment within horizon
        tt = find(E(i, t+1:min(T, t+opts.pp_change_hzn)), 1, 'first');
        if isempty(tt), continue; end
        t2 = t + tt;

        w1 = R(i,t2);
        if ~isfinite(w1) || w1 <= 0, continue; end

        % record event
        chg(end+1,1)    = (w1 - w0) / w0;              %#ok<AGROW>
        abil0(end+1,1)  = A(i,t);                       %#ok<AGROW>
        if exist('pA','var') && numel(pA)>=i && isfinite(pA(i,t))
            pabil0(end+1,1) = pA(i,t);                 %#ok<AGROW>
        else
            pabil0(end+1,1) = NaN;                     %#ok<AGROW>
        end
        isMar(end+1,1)  = M(i,t);                      %#ok<AGROW>
        hasKid(end+1,1) = (KK(i,t)==1);                %#ok<AGROW>
    end
end

if ~isempty(chg)
    nB = opts.nBinsAbility;

    % Helper: binned mean (returns nB×1 with NaN where no data)
    mean_by_bin = @(x, idx) accumarray(idx(idx>0), x(idx>0), [nB,1], @mean, NaN);

    % Precompute bin indices
    [~,~,idx_own] = histcounts(abil0,  edges);
    [~,~,idx_pnr] = histcounts(pabil0, edges);  % partner bins (NaN for singles)

    % Smoothing window (in bins)
    win = 5;  % increase for stronger smoothing; use odd numbers when possible

    % ============ Plot 1: Singles vs Married ============
    mS = (isMar==0);
    g_S = mean_by_bin(100*chg(mS), idx_own(mS));  % % change

    mM = (isMar==1);
    g_M_own = mean_by_bin(100*chg(mM), idx_own(mM));
    mMp = mM & (idx_pnr>0);
    g_M_pnr = mean_by_bin(100*chg(mMp), idx_pnr(mMp));
    g_M_both = 0.5*(g_M_own + g_M_pnr);
    miss = isnan(g_M_both) & ~isnan(g_M_own);
    g_M_both(miss) = g_M_own(miss);

    % Smooth
    g_S_s      = smoothdata(g_S(:),      'movmean', win);
    g_M_both_s = smoothdata(g_M_both(:), 'movmean', win);

    % Plot (smoothed only)
    fig = figure('Color','w');
    plot(centers(:), g_S_s,      'LineWidth',2); hold on;
    plot(centers(:), g_M_both_s, 'LineWidth',2);
    formatPlot('Initial productivity (ability)','% wage change after unemployment', ...
               'Return-to-work wage change by ability');
    legend({'Singles','Married (own + spouse bins)'}, 'Location','best','Interpreter','none');
    saveas(fig, fullfile(opts.outdir,'WageChangeAfterUnemployment_S_vs_M.png'));

    % ============ Plot 2: With kids vs Without kids ============
    k1 = (hasKid==1);
    k0 = (hasKid==0);

    g_K1_own = mean_by_bin(100*chg(k1), idx_own(k1));
    k1p = k1 & (idx_pnr>0);
    g_K1_pnr = mean_by_bin(100*chg(k1p), idx_pnr(k1p));
    g_K1_both = 0.5*(g_K1_own + g_K1_pnr);
    miss = isnan(g_K1_both) & ~isnan(g_K1_own);
    g_K1_both(miss) = g_K1_own(miss);

    g_K0_own = mean_by_bin(100*chg(k0), idx_own(k0));
    k0p = k0 & (idx_pnr>0);
    g_K0_pnr = mean_by_bin(100*chg(k0p), idx_pnr(k0p));
    g_K0_both = 0.5*(g_K0_own + g_K0_pnr);
    miss = isnan(g_K0_both) & ~isnan(g_K0_own);
    g_K0_both(miss) = g_K0_own(miss);

    % Smooth
    g_K0_both_s = smoothdata(g_K0_both(:), 'movmean', win);
    g_K1_both_s = smoothdata(g_K1_both(:), 'movmean', win);

    % Plot (smoothed only)
    fig = figure('Color','w');
    plot(centers(:), g_K0_both_s, 'LineWidth',2); hold on;
    plot(centers(:), g_K1_both_s, 'LineWidth',2);
    formatPlot('Initial productivity (ability)','% wage change after unemployment', ...
               'Return-to-work wage change by ability');
    legend({'No kids','Has kids'}, 'Location','best','Interpreter','none');
    saveas(fig, fullfile(opts.outdir,'WageChangeAfterUnemployment_Kids.png'));
end


% chg = []; abil0 = [];
% for i=1:N
%     sep_t = find(E(i,1:end-1) & ~E(i,2:end));   % t where separation between t and t+1
%     for k=1:numel(sep_t)
%         t = sep_t(k);
%         w0 = R(i,t); if ~isfinite(w0) || w0<=0, continue; end
%         % find first re-employment within horizon
%         tt = find(E(i,t+1:min(T,t+opts.pp_change_hzn)), 1,'first');
%         if isempty(tt), continue; end
%         t2 = t + tt;
%         w1 = R(i,t2); if ~isfinite(w1) || w1<=0, continue; end
%         chg(end+1,1)  = (w1 - w0) / w0; %#ok<AGROW>
%         abil0(end+1,1)= A(i,t);         %#ok<AGROW>
%     end
% end
% if ~isempty(chg)
%     [~,~,binIdx] = histcounts(abil0, edges);
%     g = accumarray(binIdx(binIdx>0), chg(binIdx>0), [opts.nBinsAbility,1], @mean, NaN);
%     fig = figure('Color','w'); plot(centers, 100*g,'-o','LineWidth',2,'MarkerSize',5);
%     formatPlot('Initial productivity (ability)','% wage change after unemployment', ...
%                'Return-to-work wage change by initial ability');
%     saveas(fig, fullfile(opts.outdir,'WageChangeAfterUnemployment_Aggregate.png'));
% end

% =====================================================================
% 6) ABILITY CORRELATION HEATMAPS (all married; employed couples)
% =====================================================================
mask_pair_all = M & isfinite(A) & isfinite(pA);
[A_f, A_m] = role_map(A, pA, gF);
A_fv = A_f(mask_pair_all); A_mv = A_m(mask_pair_all);
if ~isempty(A_fv)
    [H,~,~] = histcounts2(A_fv, A_mv, edges, edges);
    fig = figure('Color','w');
    imagesc(centers, centers, H'); axis xy; colorbar;
    formatPlot('Female Ability','Male Ability','Distribution of Married Couples');
    saveas(fig, fullfile(opts.outdir,'HeatmapAbilitiesMarried.png'));
end

mask_pair_emp = M & E & pE & isfinite(A) & isfinite(pA);
A_fv = A_f(mask_pair_emp); A_mv = A_m(mask_pair_emp);
if ~isempty(A_fv)
    [H,~,~] = histcounts2(A_fv, A_mv, edges, edges);
    fig = figure('Color','w');
    imagesc(centers, centers, H'); axis xy; colorbar;
    formatPlot('Female Ability','Male Ability','Distribution of Married Couples (Employed Only)');
    saveas(fig, fullfile(opts.outdir,'HeatmapAbilitiesMarriedNoOLF.png'));
end

% =====================================================================
% 7) TAX SCHEDULES & INCOME DISTRIBUTIONS
% =====================================================================
Wall = R(:); Wall = Wall(isfinite(Wall) & Wall>0);
if isempty(Wall), Wall = linspace(1,50,100)'; end
wgrid = linspace(min(Wall), max(Wall), 400)';

if isfield(param,'nokid_tax') && isfield(param,'onekid_tax')
    inc_pre = wgrid;
    inc_post_nokid = param.nokid_tax(inc_pre, 0, param);
    inc_post_kid   = param.onekid_tax(inc_pre, 0, param);
    atr_nokid = 1 - inc_post_nokid./max(inc_pre,SMALL);
    atr_kid   = 1 - inc_post_kid  ./max(inc_pre,SMALL);

    fig = figure('Color','w');
    yyaxis left
    plot(inc_pre, inc_pre,'k:','LineWidth',1.2); hold on;
    plot(inc_pre, inc_post_nokid,'-','LineWidth',2);
    plot(inc_pre, inc_post_kid,'--','LineWidth',2);
    ylabel('Disposable income (USD)','Interpreter','none');
    ylim([0 1.05*max(inc_pre)]); grid on;
    yyaxis right
    plot(inc_pre, 100*atr_nokid,'-','LineWidth',1.5);
    plot(inc_pre, 100*atr_kid,'--','LineWidth',1.5);
    ylabel('Average tax rate (percent)','Interpreter','none');
    xlabel('Pre-tax hourly income (USD)','Interpreter','none');
    title('Tax schedule & average tax rates','Interpreter','none');
    legend({'45 line','Disp. (no kids)','Disp. (with kids)','ATR (no kids)','ATR (with kids)'}, ...
           'Location','SouthEast','Interpreter','none');
    saveas(fig, fullfile(opts.outdir,'TaxSchedule_DisposableIncome_ATR.png'));
end

% Income distribution (employed periods)
inc = R(:); inc = inc(isfinite(inc) & inc>0);
if ~isempty(inc)
    inc_sorted = sort(inc);
    cdf_emp = (1:numel(inc_sorted))' / numel(inc_sorted);

    epsl = 1e-6;
    y = min(max(cdf_emp,epsl),1-epsl);
    logit_y = log(y./(1-y));
    deg = 5; pL = polyfit(inc_sorted, logit_y, deg);
    xgrid = linspace(min(inc_sorted), max(inc_sorted), 400);
    ell = polyval(pL, xgrid);
    cdf_fit = 1./(1+exp(-ell));
    ellp = polyval(polyder(pL), xgrid);
    pdf_fit = ellp .* exp(-ell) ./ (1+exp(-ell)).^2;
    pdf_fit(pdf_fit<0) = 0;

    % PDF
    fig = figure('Color','w');
    plot(xgrid, pdf_fit,'LineWidth',2); grid on;
    xlabel('Pre-tax hourly income (USD)','Interpreter','none');
    ylabel('Density','Interpreter','none');
    title('Fitted PDF of pre-tax hourly income','Interpreter','none');
    saveas(fig, fullfile(opts.outdir,'PDF_PreTaxIncome.png'));

    % CDF compare
    fig = figure('Color','w');
    stairs(inc_sorted, cdf_emp,'k','LineWidth',1); hold on;
    plot(xgrid, cdf_fit,'r--','LineWidth',2); grid on;
    xlabel('Pre-tax hourly income (USD)','Interpreter','none');
    ylabel('CDF','Interpreter','none');
    legend({'Empirical CDF','Fitted CDF'},'Location','SouthEast','Interpreter','none');
    title('Goodness-of-fit: empirical vs. fitted CDF','Interpreter','none');
    saveas(fig, fullfile(opts.outdir,'CDF_PreTaxIncome_Fit.png'));
end
% =====================================================================
% 8) AGE PROFILES: ability vs age; income vs age; wage vs age
% =====================================================================
if ~isfield(agentPanel,'age')
    warning('agentPanel.age not found. Skipping age-profile plots.');
else
    AgeM = agentPanel.age;              % N x T, months since (re)birth
    ppy  = 12;
    if isfield(param,'periods_per_year') && ~isempty(param.periods_per_year)
        ppy = param.periods_per_year;
    end

    % Helper: compute mean(Y) by integer age (months). Returns age grid in years.
    mean_by_age = @(Y) local_mean_by_age(AgeM, Y, ppy);

    % --- (a) Ability vs age (use everyone; ability always defined) ---
    [ageY, muAbility] = mean_by_age(A);
    fig = figure('Color','w');
    plot(ageY, muAbility, 'LineWidth',2);
    formatPlot('Age (years)','Average ability','Ability vs age');
    saveas(fig, fullfile(opts.outdir,'AgeProfile_Ability.png'));

    % --- (b) Income vs age: exp(a)*w, employed only ---
    R_age = W;     % income measure consistent with the rest
    R_age(~E) = NaN;          % only meaningful when employed
    [ageY, muIncome] = mean_by_age(R_age);
    fig = figure('Color','w');
    plot(ageY, muIncome, 'LineWidth',2);
    formatPlot('Age (years)','Average income (exp(a) * w)','Income vs age (employed only)');
    saveas(fig, fullfile(opts.outdir,'AgeProfile_Income.png'));

    % --- (c) Wage vs age: grid wage w, employed only ---
    W_age = W;
    W_age(~E) = NaN;
    [ageY, muWage] = mean_by_age(W_age);
    fig = figure('Color','w');
    plot(ageY, muWage, 'LineWidth',2);
    formatPlot('Age (years)','Average wage (w)','Wage vs age (employed only)');
    saveas(fig, fullfile(opts.outdir,'AgeProfile_Wage.png'));
end
%% ================= Figure: Value of Marrying an Average Productivity Person =================
% Productivity grid on the x-axis (P = exp(a))
Pgrid  = 0.5:0.5:4;                        % as in the paper
agrid  = log(Pgrid);                       % ability points behind those productivities

% Convenience
[N,T] = size(agentPanel.employed);
fem   = logical(agentPanel.female(:));
male  = ~fem;
maskF = repmat(fem,1,T);
maskM = repmat(male,1,T);
E     = agentPanel.employed==1;

% Typical wages (use observed wages among employed; robust fallbacks)
wf_obs = agentPanel.wage(E & maskF);
wm_obs = agentPanel.wage(E & maskM);

w_ref_f = median(wf_obs,'omitnan');
if isempty(w_ref_f) || ~isfinite(w_ref_f)
    w_ref_f = median(agentPanel.wage(:),'omitnan');
end

w_ref_m = median(wm_obs,'omitnan');
if isempty(w_ref_m) || ~isfinite(w_ref_m)
    % fall back to female ref or overall if needed
    tmp = median(agentPanel.wage(:),'omitnan');
    if isfinite(tmp), w_ref_m = tmp; else, w_ref_m = w_ref_f; end
end

% Average abilities by sex (for the "average partner")
a_bar_f = mean(agentPanel.ability(maskF),'omitnan');
a_bar_m = mean(agentPanel.ability(maskM),'omitnan');

% ---------- (A) Female perspective: own ability varies, partner fixed at avg male ----------
V_pairWW_F = zeros(numel(agrid),1);
V_pairWU_F = zeros(numel(agrid),1);
V_singleW_F= zeros(numel(agrid),1);
for j = 1:numel(agrid)
    a_self = agrid(j);
    % pair: (Working, Working)
    V_pairWW_F(j) = getMarriedValue(a_self, 1, w_ref_f, a_bar_m, 1, w_ref_m, 0, true,  valueFunc, 0, 0);
    % pair: (Working, Unemployed partner)
    V_pairWU_F(j) = getMarriedValue(a_self, 1, w_ref_f, a_bar_m, 0, 0,       0, true,  valueFunc, 0, 0);
    % single: (Working)
    V_singleW_F(j)= getSingleValue (a_self, 1, w_ref_f, 0, 0, true, valueFunc);
end

figure; hold on; box on; grid on;
plot(Pgrid, V_pairWW_F, '-s', 'LineWidth',1.5, 'MarkerSize',5, 'DisplayName','Pair (Working,Working)');
plot(Pgrid, V_singleW_F,'--', 'LineWidth',1.5,               'DisplayName','Single (Working)');
plot(Pgrid, V_pairWU_F, '-o', 'LineWidth',1.5, 'MarkerSize',5, 'DisplayName','Pair (Working,Unemployed)');
xlabel('Productivity'); ylabel('Value');
title('Value of Marrying an Average Productivity Person — Female');
legend('Location','best');

% ---------- (B) Male perspective: own ability varies, partner fixed at avg female ----------
V_pairWW_M = zeros(numel(agrid),1);
V_pairWU_M = zeros(numel(agrid),1);
V_singleW_M= zeros(numel(agrid),1);
for j = 1:numel(agrid)
    a_self = agrid(j);
    V_pairWW_M(j) = getMarriedValue(a_self, 1, w_ref_m, a_bar_f, 1, w_ref_f, 0, false, valueFunc, 0, 0);
    V_pairWU_M(j) = getMarriedValue(a_self, 1, w_ref_m, a_bar_f, 0, 0,       0, false, valueFunc, 0, 0);
    V_singleW_M(j)= getSingleValue (a_self, 1, w_ref_m, 0, 0, false, valueFunc);
end

figure; hold on; box on; grid on;
plot(Pgrid, V_pairWW_M, '-s', 'LineWidth',1.5, 'MarkerSize',5, 'DisplayName','Pair (Working,Working)');
plot(Pgrid, V_singleW_M,'--', 'LineWidth',1.5,               'DisplayName','Single (Working)');
plot(Pgrid, V_pairWU_M, '-o', 'LineWidth',1.5, 'MarkerSize',5, 'DisplayName','Pair (Working,Unemployed)');
xlabel('Productivity'); ylabel('Value');
title('Value of Marrying an Average Productivity Person — Male');
legend('Location','best');

%% =========== Additional Figure: Average respondent vs partner ability  ===========
% Fix the respondent at their group average ability; vary the partner’s ability (productivity)
Pgrid_p  = Pgrid;                   % reuse the same axis
agrid_p  = agrid;

% --- Female respondent fixed at a_bar_f; vary male partner ability
V_pairWW_F_fix = zeros(numel(agrid_p),1);
V_pairWU_F_fix = zeros(numel(agrid_p),1);
V_singleW_F_fix= getSingleValue(a_bar_f, 1, w_ref_f, 0, 0, true, valueFunc) * ones(numel(agrid_p),1);
for j = 1:numel(agrid_p)
    a_part = agrid_p(j);
    V_pairWW_F_fix(j) = getMarriedValue(a_bar_f, 1, w_ref_f, a_part, 1, w_ref_m, 0, true,  valueFunc, 0, 0);
    V_pairWU_F_fix(j) = getMarriedValue(a_bar_f, 1, w_ref_f, a_part, 0, 0,       0, true,  valueFunc, 0, 0);
end

figure; hold on; box on; grid on;
plot(Pgrid_p, V_pairWW_F_fix, '-s', 'LineWidth',1.5, 'MarkerSize',5, 'DisplayName','Pair (Working,Working)');
plot(Pgrid_p, V_singleW_F_fix,'--', 'LineWidth',1.5,               'DisplayName','Single (Working)');
plot(Pgrid_p, V_pairWU_F_fix, '-o', 'LineWidth',1.5, 'MarkerSize',5, 'DisplayName','Pair (Working,Unemployed)');
xlabel('Partner Productivity'); ylabel('Value');
title('Average Female Respondent: Value vs Partner Ability');
legend('Location','best');

% --- Male respondent fixed at a_bar_m; vary female partner ability
V_pairWW_M_fix = zeros(numel(agrid_p),1);
V_pairWU_M_fix = zeros(numel(agrid_p),1);
V_singleW_M_fix= getSingleValue(a_bar_m, 1, w_ref_m, 0, 0, false, valueFunc) * ones(numel(agrid_p),1);
for j = 1:numel(agrid_p)
    a_part = agrid_p(j);
    V_pairWW_M_fix(j) = getMarriedValue(a_bar_m, 1, w_ref_m, a_part, 1, w_ref_f, 0, false, valueFunc, 0, 0);
    V_pairWU_M_fix(j) = getMarriedValue(a_bar_m, 1, w_ref_m, a_part, 0, 0,       0, false, valueFunc, 0, 0);
end

figure; hold on; box on; grid on;
plot(Pgrid_p, V_pairWW_M_fix, '-s', 'LineWidth',1.5, 'MarkerSize',5, 'DisplayName','Pair (Working,Working)');
plot(Pgrid_p, V_singleW_M_fix,'--', 'LineWidth',1.5,               'DisplayName','Single (Working)');
plot(Pgrid_p, V_pairWU_M_fix, '-o', 'LineWidth',1.5, 'MarkerSize',5, 'DisplayName','Pair (Working,Unemployed)');
xlabel('Partner Productivity'); ylabel('Value');
title('Average Male Respondent: Value vs Partner Ability');
legend('Location','best');

%% === ONE FIGURE: Unemployment marital gap (women vs men) ===
% gap_* are in percentage points (S − M)
fig = figure('Color','w'); hold on; box on; grid on;
plot(centers, gap_f, '-','LineWidth',2);    % women
plot(centers, gap_m, '--','LineWidth',2);   % men
formatPlot('Productivity','pp','Unemployment Marital Gap (Single − Married)');
legend({'Women','Men'}, 'Location','best','Interpreter','none');
saveas(fig, fullfile(opts.outdir,'EndogUnemploymentDiff_FvM.png'));

%% === ONE FIGURE: Wage marital premium (women vs men) ===
% prem_* are in percent ((M − S)/S * 100)
fig = figure('Color','w'); hold on; box on; grid on;
plot(centers, prem_f, '-','LineWidth',2);   % women
plot(centers, prem_m, '--','LineWidth',2);  % men
formatPlot('Productivity','% Wage difference','Wage Marital Premium (Married − Single)');
legend({'Women','Men'}, 'Location','best','Interpreter','none');
saveas(fig, fullfile(opts.outdir,'EndogWageDiff_FvM.png'));


disp('PlotSimul: all figures saved.');
end

% ==================== helpers ====================
function out = bin_by_ability(A, X, mask, edges, agg, denom_mask)
% A: NT ability
% X: NT target (NaN ok) or [] when agg='share'/'mass'
% mask: logical NT selecting records
% agg: 'mean' | 'share' | 'mass' | 'share_in_group'
% denom_mask: (optional) group mask for 'share_in_group'
    if nargin<6, denom_mask = []; end
    if isempty(mask), mask = true(size(A)); end
    a = A(:); m = mask(:);
    [~,~,b] = histcounts(a, edges);   % 0 outside bins
    nb = numel(edges)-1;
    out = nan(nb,1);

    switch agg
        case 'mean'
            x = X(:);
            valid = m & isfinite(x);
            idx = b(valid); idx = idx(idx>0);
            val = x(valid); val = val(idx>0); %#ok<NASGU>
            out = accumarray(idx, x(valid & b>0), [nb,1], @mean, NaN);
        case 'share'
            valid = m & b>0;
            out = accumarray(b(valid), 1, [nb,1], @sum, 0);
            out = out / max(sum(valid),1);
        case 'mass'
            valid = m & b>0;
            out = accumarray(b(valid), 1, [nb,1], @sum, 0);
        case 'share_in_group'
            if isempty(denom_mask), error('denom_mask needed for share_in_group'); end
            dm = denom_mask(:);
            num = accumarray(b(m & b>0), 1, [nb,1], @sum, 0);
            den = accumarray(b(dm & b>0), 1, [nb,1], @sum, 0);
            out = num ./ max(den,1);  % fraction within group per bin
        otherwise
            error('Unknown agg');
    end
end

function [A_f, A_m] = role_map(A, pA, femaleMask)
% Map abilities to female/male roles for married pairs
    A_f = NaN(size(A)); A_m = NaN(size(A));
    f = femaleMask; m = ~f;
    A_f(f) = A(f);    A_f(m) = pA(m);
    A_m(f) = pA(f);   A_m(m) = A(m);
end

function formatPlot(xlabel_text, ylabel_text, title_text)
% Standardized plot formatting without LaTeX (avoids %/^ issues)
    set(gca,'FontSize',12,'FontName','Times');
    xlabel(xlabel_text,'Interpreter','none','FontSize',16);
    ylabel(ylabel_text,'Interpreter','none','FontSize',16);
    title(title_text,'Interpreter','none','FontSize',16);
    grid on;
end
function [age_years, m] = local_mean_by_age(AgeM, Y, ppy)
% local_mean_by_age  Mean of Y by integer age in months; returns age in years.
% Inputs:
%   AgeM : N x T integer ages (months)
%   Y    : N x T values (NaN OK)
%   ppy  : periods per year (e.g., 12)
% Outputs:
%   age_years : vector of unique ages (years)
%   m         : mean(Y) at each age (same length)

    a = AgeM(:); y = Y(:);
    valid = isfinite(a) & isfinite(y);
    if ~any(valid)
        age_years = []; m = [];
        return;
    end

    % Ensure integer month bins (ages produced by the sim should already be ints)
    a = round(a(valid));
    y = y(valid);

    % Shift to 1-based indexing for accumarray
    a_idx = a + 1;
    maxA  = max(a_idx);
    sumY  = accumarray(a_idx, y, [maxA 1], @nansum, 0);
    cntY  = accumarray(a_idx, 1, [maxA 1], @sum,    0);

    keep  = cntY > 0;
    m     = sumY(keep) ./ cntY(keep);
    age_months = (find(keep)-1);          % back to 0-based months
    age_years  = age_months / max(ppy,1); % convert to years
end
% ===================== Local value helpers =====================
function V = getSingleValue(ability, employed, wage, olf, has_kids, is_female, valueFunc)
    % Reads from valueFunc.single_female / single_male
    if is_female
        vF = valueFunc.single_female;
    else
        vF = valueFunc.single_male;
    end

    if employed
        if has_kids
            V = vF.employed_kid(ability, wage);
        else
            V = vF.employed_nokid(ability, wage);
        end
    elseif ~olf
        if has_kids
            V = vF.unemployed_kid(ability);
        else
            V = vF.unemployed_nokid(ability);
        end
    else
        if has_kids
            V = vF.olf_kid(ability);
        else
            V = vF.olf_nokid(ability);
        end
    end
end

function V = getMarriedValue(ability_self, employed_self, wage_self, ...
                             ability_partner, employed_partner, wage_partner, ...
                             has_kids, is_female, valueFunc, olf_self, olf_partner)
    % Default OLF flags if omitted
    if nargin < 10 || isempty(olf_self),    olf_self    = 0; end
    if nargin < 11 || isempty(olf_partner), olf_partner = 0; end

    % Robust logicals (treat NaN/0 as false)
    ef = isfinite(employed_self)    && (employed_self~=0);
    em = isfinite(employed_partner) && (employed_partner~=0);
    of = isfinite(olf_self)         && (olf_self~=0);
    om = isfinite(olf_partner)      && (olf_partner~=0);
    hk = isfinite(has_kids)         && (has_kids~=0);

    % Replace NaN wages with 0
    if ~isfinite(wage_self),    wage_self = 0; end
    if ~isfinite(wage_partner), wage_partner = 0; end

    vC = valueFunc.couple;

    % Map arguments to female/male roles expected by couple functions
    if is_female
        af = ability_self;  wf = wage_self;  ef0 = ef;  of0 = of;
        am = ability_partner; wm = wage_partner; em0 = em;  om0 = om;
    else
        af = ability_partner; wf = wage_partner; ef0 = em;  of0 = om;
        am = ability_self;    wm = wage_self;    em0 = ef;  om0 = of;
    end

    % Determine the couple state name base
    if     ef0 && ~of0 && em0 && ~om0
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

    % Kids suffix
    if hk, suf = '_kid'; else, suf = '_nokid'; end
    fn = [nm suf];

    % Fallbacks if a specific function is missing in valueFunc.couple
    if ~isfield(vC, fn) || ~isa(vC.(fn),'function_handle')
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

