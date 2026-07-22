function PlotSimul_Extras(agentPanel, param, opts)
% PlotSimul_Extras  Additional, paper-style figures from the simulated panel.
% Usage: PlotSimul_Extras(agentPanel, param)
%        PlotSimul_Extras(agentPanel, param, struct('outdir','Plots','nBinsAbility',25))
%
% Figures produced (saved under opts.outdir):
%   1) AvgWage_MarriedVsSingle.png          — Figure 5 analog (wage levels)
%   2) UnempRate_MarriedVsSingle.png        — Figure 5 analog (unemployment levels)
%   3) MarriedShare_byWageQuintile.png      — Appendix A.1 analog (by wage quintiles)
%   4) RankRank_Wages_EmployedCouples.png   — Appendix A.1 analog (rank–rank wages in couples)
%   5) SpouseEmploymentByAbility.png        — diagnostic: P(spouse employed | own ability)
%   6) EventStudy_AroundMarriage.png        — dynamics of wage & unemployment around 1st marriage
%   7) UnempDuration_byMaritalAtStart.png   — distribution of unemployment spell lengths by status
%
% This file is additive: it does not change PlotSimul.m. To fix the
% "share of married by ability" panel inside PlotSimul.m as requested,
% see the one-line patch provided in the main chat.

if nargin<3, opts=struct(); end
if ~isfield(opts,'outdir'), opts.outdir = 'PlotsExtras'; end
if ~isfield(opts,'nBinsAbility'), opts.nBinsAbility = 25; end

if ~exist(opts.outdir,'dir'), mkdir(opts.outdir); end
SMALL = 1e-10;

% ------------------- unpack panel -------------------
[N,T] = size(agentPanel.employed);
fem   = logical(agentPanel.female(:));
male  = ~fem;

A  = agentPanel.ability;                       % N×T
Wg = agentPanel.wage;                          % N×T (grid wage)
E  = agentPanel.employed==1;
O  = agentPanel.olf==1;
U  = ~E & ~O;
M  = agentPanel.married==1;

% partner states (robust)
pE = agentPanel.partner_employed==1; pE(isnan(agentPanel.partner_employed)) = false;

% real wages (only for employed)
R = exp(A) .* Wg;    R(~E) = NaN;

% ability bins (like PlotSimul)
Aall = A(:); Aall = Aall(isfinite(Aall));
edgesA   = linspace(min(Aall), max(Aall), opts.nBinsAbility+1);
centersA = 0.5*(edgesA(1:end-1) + edgesA(2:end));

% helper handles
maskF = repmat(fem,1,T);
maskM = repmat(male,1,T);
binA  = @(X,mask,agg,varargin) bin_by_ability(A, X, mask, edgesA, agg, varargin{:});

%% =====================================================================
%% 1) Figure 5 analog: Average WAGES and UNEMPLOYMENT by marital status
%% =====================================================================
% Average real wage among employed, by current marital status
W_emp = R; % already NaN where not employed
wM = mean(W_emp(M), 'omitnan');
wS = mean(W_emp(~M), 'omitnan');
fig = figure('Color','w');
bar([wS, wM]); set(gca,'XTickLabel',{'Single','Married'});
formatPlot('Marital status','Average real wage','Average wages by marital status');
saveas(fig, fullfile(opts.outdir,'AvgWage_MarriedVsSingle.png'));
close(fig);

% Unemployment rate by marital status
uM = mean(U(M));
uS = mean(U(~M));
fig = figure('Color','w');
bar(100*[uS, uM]); set(gca,'XTickLabel',{'Single','Married'});
formatPlot('Marital status','Unemployment rate (pp)','Unemployment by marital status');
saveas(fig, fullfile(opts.outdir,'UnempRate_MarriedVsSingle.png'));
close(fig);

%% =====================================================================
%% 2) Appendix analog: Married share by wage quintile & wage sorting
%% =====================================================================
% Married share by *own* wage quintile (among employed periods)
Wvec = R(:); Wvec = Wvec(isfinite(Wvec) & Wvec>0);
if ~isempty(Wvec)
    q = quantile(Wvec, [0 .2 .4 .6 .8 1]);
    WQ = discretize(R, q);   % NaN where not employed
    married_share_q = nan(5,1); mass_q = nan(5,1);
    for qk=1:5
        inq = (WQ==qk);
        mass_q(qk) = nnz(inq);
        if mass_q(qk)>0
            married_share_q(qk) = nnz(inq & M) / mass_q(qk);
        end
    end
    fig = figure('Color','w');
    bar(100*married_share_q); set(gca,'XTickLabel',{'Q1','Q2','Q3','Q4','Q5'});
    formatPlot('Own wage quintile','% married','Married share by wage quintile (employed)');
    saveas(fig, fullfile(opts.outdir,'MarriedShare_byWageQuintile.png'));
    close(fig);
end

% Rank–rank wages in employed couples (both employed same period)
maskCoupleWE = M & E & pE & isfinite(R) & isfinite(agentPanel.partner_wage) & isfinite(agentPanel.partner_ability);
if any(maskCoupleWE(:))
    % construct spouses' real wages
    pR = exp(agentPanel.partner_ability).*agentPanel.partner_wage;
    x = R(maskCoupleWE); y = pR(maskCoupleWE);
    % percentile ranks
    xr = tiedrank(x)/numel(x);
    yr = tiedrank(y)/numel(y);
    % bin own rank, plot mean partner rank
    nb = 20; edgesR = linspace(0,1,nb+1); centersR = 0.5*(edgesR(1:end-1)+edgesR(2:end));
    [~,~,bx] = histcounts(xr, edgesR);
    rr = accumarray(bx(bx>0), yr(bx>0), [nb,1], @mean, NaN);
    fig = figure('Color','w');
    plot(centersR, rr, '-o','LineWidth',2,'MarkerSize',5); hold on;
    plot([0 1],[0 1],'k:');
    formatPlot('Own wage rank','Partner wage rank','Rank–rank wages in employed couples');
    ylim([0 1]); xlim([0 1]);
    saveas(fig, fullfile(opts.outdir,'RankRank_Wages_EmployedCouples.png'));
    close(fig);
end

%% =====================================================================
%% 3) Diagnostic: P(spouse employed | own ability) by gender
%% =====================================================================
pe_f = binA([], M & maskF & pE, 'share_in_group', M & maskF);
pe_m = binA([], M & maskM & pE, 'share_in_group', M & maskM);
fig = figure('Color','w');
plot(centersA, 100*pe_f, '-','LineWidth',2); hold on;
plot(centersA, 100*pe_m, '--','LineWidth',2);
formatPlot('Own ability','% spouse employed','Spouse employment by own ability (married only)');
legend({'Female','Male'},'Location','best','Interpreter','none');
saveas(fig, fullfile(opts.outdir,'SpouseEmploymentByAbility.png'));
close(fig);

%% =====================================================================
%% 4) Event study around first marriage (wage & unemployment)
%% =====================================================================
% window K periods before/after first transition into marriage
K = min(6, max(3, round(T/10)));  % up to ±6 or 10% of T
first_mar_t = nan(N,1);
for i=1:N
    m = find(~M(i,1:end-1) & M(i,2:end), 1, 'first');
    if ~isempty(m), first_mar_t(i) = m+1; end
end
sel = find(isfinite(first_mar_t));
if ~isempty(sel)
    Wrel = nan(numel(sel), 2*K+1);
    Urel = nan(numel(sel), 2*K+1);
    for kk=1:numel(sel)
        i = sel(kk); t0 = first_mar_t(i);
        tL = max(1, t0-K); tR = min(T, t0+K);
        idx = (tL:tR);
        Wrel(kk, (K+1)-(t0-tL) : (K+1)+(tR-t0)) = R(i, idx);
        Urel(kk, (K+1)-(t0-tL) : (K+1)+(tR-t0)) = U(i, idx);
    end
    muW = mean(Wrel, 1, 'omitnan');
    muU = mean(Urel, 1, 'omitnan');
    tau  = (-K:K);
    % wages (conditional on being employed)
    fig = figure('Color','w');
    plot(tau, muW, '-o','LineWidth',2,'MarkerSize',4);
    formatPlot('Event time (0 = first marriage)','Average real wage','Event study: wages around marriage');
    xline(0,'k:');
    saveas(fig, fullfile(opts.outdir,'EventStudy_AroundMarriage_Wage.png'));
    close(fig);
    % unemployment
    fig = figure('Color','w');
    plot(tau, 100*muU, '-o','LineWidth',2,'MarkerSize',4);
    formatPlot('Event time (0 = first marriage)','Unemployment rate (pp)','Event study: unemployment around marriage');
    xline(0,'k:');
    saveas(fig, fullfile(opts.outdir,'EventStudy_AroundMarriage_Unemp.png'));
    close(fig);
end

%% =====================================================================
%% 5) Unemployment duration distribution by marital status at spell start
%% =====================================================================
[durS, durM] = unemployment_spell_lengths(U, M);
if ~isempty(durS) || ~isempty(durM)
    maxD = max([durS; durM; 0]);
    e = 1:max(1,maxD);
    hS = histcounts(durS, [e, inf], 'Normalization','probability');
    hM = histcounts(durM, [e, inf], 'Normalization','probability');
    fig = figure('Color','w');
    stairs(e, 100*hS, 'LineWidth',1.8); hold on;
    stairs(e, 100*hM, '--','LineWidth',1.8);
    formatPlot('Unemployment spell length (periods)','Share of spells (pp)','Unemployment durations by marital status at start');
    legend({'Single at start','Married at start'},'Location','northeast','Interpreter','none');
    saveas(fig, fullfile(opts.outdir,'UnempDuration_byMaritalAtStart.png'));
    close(fig);
end

end % function

%% ==================== helpers (local) ====================
function out = bin_by_ability(A, X, mask, edges, agg, denom_mask)
% bin_by_ability  (copied from PlotSimul for consistency)
    if nargin<6, denom_mask = []; end
    if isempty(mask), mask = true(size(A)); end
    a = A(:); m = mask(:);
    [~,~,b] = histcounts(a, edges);
    nb = numel(edges)-1; out = nan(nb,1);
    switch agg
        case 'mean'
            x = X(:); valid = m & isfinite(x) & b>0;
            out = accumarray(b(valid), x(valid), [nb,1], @mean, NaN);
        case 'share'
            valid = m & b>0; out = accumarray(b(valid), 1, [nb,1], @sum, 0);
            out = out / max(sum(valid),1);
        case 'mass'
            valid = m & b>0; out = accumarray(b(valid), 1, [nb,1], @sum, 0);
        case 'share_in_group'
            if isempty(denom_mask), error('denom_mask needed for share_in_group'); end
            dm = denom_mask(:);
            num = accumarray(b(m & b>0), 1, [nb,1], @sum, 0);
            den = accumarray(b(dm & b>0), 1, [nb,1], @sum, 0);
            out = num ./ max(den,1);
        otherwise, error('Unknown agg');
    end
end

function formatPlot(xlabel_text, ylabel_text, title_text)
    set(gca,'FontSize',12,'FontName','Times');
    xlabel(xlabel_text,'Interpreter','none','FontSize',16);
    ylabel(ylabel_text,'Interpreter','none','FontSize',16);
    title(title_text,'Interpreter','none','FontSize',16);
    grid on;
end

function [durS, durM] = unemployment_spell_lengths(U, M)
% Extract spell lengths in unemployment; attribute to marital status at spell start.
    [N,T] = size(U); durS = []; durM = [];
    for i=1:N
        u = U(i,:);
        start_idx = find([false, ~u(1:end-1) & u(2:end)]);
        end_idx   = find([u(1:end-1) & ~u(2:end), false]);
        % ensure matching lengths
        k = min(numel(start_idx), numel(end_idx));
        start_idx = start_idx(1:k); end_idx = end_idx(1:k);
        for j=1:k
            d = end_idx(j) - start_idx(j) + 1;  % inclusive
            if d>0
                if M(i, start_idx(j))
                    durM(end+1,1) = d; %#ok<AGROW>
                else
                    durS(end+1,1) = d; %#ok<AGROW>
                end
            end
        end
    end
end
