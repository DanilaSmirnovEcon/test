function createSummaryPlots(solutionSP, param, gridSf, gridSm, gridP)
    % CREATESUMMARYPLOTS - Generate key results plots illustrating economic mechanisms
    % This function creates diverse plots highlighting the paper's main mechanisms
    
    close all;
    
    % Create plots directory if it doesn't exist
    if ~exist('Plots', 'dir')
        mkdir('Plots');
    end
    
    % Calculate aggregate statistics
    stats = calculateAggregateStats(solutionSP, param, gridSf, gridSm, gridP);
    
    %% 1. CORE LABOR MARKET OUTCOMES - Bar charts
    createLaborMarketSummary(stats);
    
    %% 2. MARRIAGE MARKET SORTING - Scatter plot
    createAssortativeMatingPlot(solutionSP, gridP);
    
    %% 3. WAGE CORRELATION IN COUPLES - Scatter plot
    createWageCorrelationPlot(solutionSP, gridP, stats);
    
    %% 4. JOB SEARCH INTENSITY BY MARITAL STATUS - Line plot
    createJobSearchIntensityPlot(stats);
    
    %% 5. EMPLOYMENT TRANSITIONS - Grouped bar chart
    createTransitionComparisonPlot(stats);
    
    %% 6. WAGE PROGRESSION BY MARRIAGE STATUS - Line plot
    createWageProgressionPlot(solutionSP, gridSf, gridSm, gridP);
    
    close all;
end

function stats = calculateAggregateStats(solutionSP, param, gridSf, gridSm, gridP)
    % Calculate comprehensive aggregate statistics
    
    SMALL_NUMBER = 1e-10;
    
    % Basic population shares
    total_pop = sum(solutionSP.fNew);
    
    % Singles
    single_f = sum(solutionSP.efl, 'all') + sum(solutionSP.hfl, 'all') + sum(solutionSP.ofl, 'all') + ...
               sum(solutionSP.efc, 'all') + sum(solutionSP.hfc, 'all') + sum(solutionSP.ofc, 'all');
    single_m = sum(solutionSP.eml, 'all') + sum(solutionSP.hml, 'all') + sum(solutionSP.oml, 'all') + ...
               sum(solutionSP.emc, 'all') + sum(solutionSP.hmc, 'all') + sum(solutionSP.omc, 'all');
    
    % Married (count each person once)
    married_total = total_pop - single_f - single_m;
    
    % Employment status calculations
    single_f_emp = sum(solutionSP.efl, 'all') + sum(solutionSP.efc, 'all');
    single_f_unemp = sum(solutionSP.hfl, 'all') + sum(solutionSP.hfc, 'all');
    single_m_emp = sum(solutionSP.eml, 'all') + sum(solutionSP.emc, 'all');
    single_m_unemp = sum(solutionSP.hml, 'all') + sum(solutionSP.hmc, 'all');
    
    married_f_emp = sum(solutionSP.eel, 'all') + sum(solutionSP.ehl, 'all') + sum(solutionSP.eol, 'all') + ...
                    sum(solutionSP.eec, 'all') + sum(solutionSP.ehc, 'all') + sum(solutionSP.eoc, 'all');
    married_f_unemp = sum(solutionSP.hel, 'all') + sum(solutionSP.hhl, 'all') + sum(solutionSP.hol, 'all') + ...
                      sum(solutionSP.hec, 'all') + sum(solutionSP.hhc, 'all') + sum(solutionSP.hoc, 'all');
    married_m_emp = sum(solutionSP.eel, 'all') + sum(solutionSP.hel, 'all') + sum(solutionSP.oel, 'all') + ...
                    sum(solutionSP.eec, 'all') + sum(solutionSP.hec, 'all') + sum(solutionSP.oec, 'all');
    married_m_unemp = sum(solutionSP.ehl, 'all') + sum(solutionSP.hhl, 'all') + sum(solutionSP.ohl, 'all') + ...
                      sum(solutionSP.ehc, 'all') + sum(solutionSP.hhc, 'all') + sum(solutionSP.ohc, 'all');
    
    % Calculate key rates
    stats.employment_rates = struct();
    stats.employment_rates.single_f = 100 * single_f_emp / (single_f_emp + single_f_unemp + SMALL_NUMBER);
    stats.employment_rates.single_m = 100 * single_m_emp / (single_m_emp + single_m_unemp + SMALL_NUMBER);
    stats.employment_rates.married_f = 100 * married_f_emp / (married_f_emp + married_f_unemp + SMALL_NUMBER);
    stats.employment_rates.married_m = 100 * married_m_emp / (married_m_emp + married_m_unemp + SMALL_NUMBER);
    
    stats.unemployment_rates = struct();
    stats.unemployment_rates.single_f = 100 * single_f_unemp / (single_f_emp + single_f_unemp + SMALL_NUMBER);
    stats.unemployment_rates.single_m = 100 * single_m_unemp / (single_m_emp + single_m_unemp + SMALL_NUMBER);
    stats.unemployment_rates.married_f = 100 * married_f_unemp / (married_f_emp + married_f_unemp + SMALL_NUMBER);
    stats.unemployment_rates.married_m = 100 * married_m_unemp / (married_m_emp + married_m_unemp + SMALL_NUMBER);
    
    % Marriage rate
    stats.marriage_rate = 100 * married_total / total_pop;
    
    % Store raw numbers for other calculations
    stats.raw = struct();
    stats.raw.single_f_emp = single_f_emp;
    stats.raw.single_m_emp = single_m_emp;
    stats.raw.married_f_emp = married_f_emp;
    stats.raw.married_m_emp = married_m_emp;
    stats.raw.total_pop = total_pop;
    stats.raw.married_total = married_total;
end

function createLaborMarketSummary(stats)
    % Create summary bar charts for key labor market outcomes
    
    %% Employment Rates Comparison
    fig = figure('Color', 'w', 'Position', [100, 100, 800, 600]);
    
    categories = {'Single Women', 'Married Women', 'Single Men', 'Married Men'};
    employment_data = [stats.employment_rates.single_f, stats.employment_rates.married_f, ...
                      stats.employment_rates.single_m, stats.employment_rates.married_m];
    
    bar_colors = [0.8, 0.3, 0.3; 0.3, 0.3, 0.8; 0.8, 0.3, 0.3; 0.3, 0.3, 0.8];
    
    b = bar(employment_data, 'FaceColor', 'flat');
    b.CData = bar_colors;
    
    set(gca, 'XTickLabel', categories, 'FontSize', 14, 'FontName', 'Times');
    ylabel('Employment Rate (\%)', 'Interpreter', 'latex', 'FontSize', 18, 'FontName', 'Times');
    title('Employment Rates by Gender and Marital Status', 'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Times');
    grid on;
    ylim([0, max(employment_data) * 1.1]);
    
    % Add value labels on bars
    for i = 1:length(employment_data)
        text(i, employment_data(i) + max(employment_data) * 0.02, ...
             sprintf('%.1f', employment_data(i)), ...
             'HorizontalAlignment', 'center', 'FontSize', 12);
    end
    
    saveas(fig, 'Plots/Summary_EmploymentRates.png');
    
    %% Marriage Rate
    fig = figure('Color', 'w', 'Position', [100, 100, 400, 600]);
    
    b = bar(stats.marriage_rate, 'FaceColor', [0.2, 0.6, 0.2]);
    
    set(gca, 'XTickLabel', {'Population'}, 'FontSize', 14, 'FontName', 'Times');
    ylabel('Marriage Rate (\%)', 'Interpreter', 'latex', 'FontSize', 18, 'FontName', 'Times');
    title('Overall Marriage Rate', 'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Times');
    grid on;
    ylim([0, stats.marriage_rate * 1.2]);
    
    text(1, stats.marriage_rate + stats.marriage_rate * 0.05, ...
         sprintf('%.1f', stats.marriage_rate), ...
         'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
    
    saveas(fig, 'Plots/Summary_MarriageRate.png');
end

function createAssortativeMatingPlot(solutionSP, gridP)
    % Create scatter plot showing assortative mating patterns
    
    % Calculate joint distribution of married couples by ability
    joint_dist = squeeze(sum(solutionSP.eel, [2,4])) + ...
                 squeeze(sum(solutionSP.ehl, 2)) + ...
                 squeeze(sum(solutionSP.eol, 2)) + ...
                 squeeze(sum(solutionSP.hel, 3)) + ...
                 squeeze(sum(solutionSP.oel, 3)) + ...
                 solutionSP.hhl + solutionSP.hol + solutionSP.ohl + solutionSP.ool + ...
                 squeeze(sum(solutionSP.eec, [2,4])) + ...
                 squeeze(sum(solutionSP.ehc, 2)) + ...
                 squeeze(sum(solutionSP.eoc, 2)) + ...
                 squeeze(sum(solutionSP.hec, 3)) + ...
                 squeeze(sum(solutionSP.oec, 3)) + ...
                 solutionSP.hhc + solutionSP.hoc + solutionSP.ohc + solutionSP.ooc;
    
    % Create mesh grids
    [Female_Ability, Male_Ability] = meshgrid(gridP.af.a, gridP.am.a);
    
    % Flatten for scatter plot
    female_abilities = Female_Ability(:);
    male_abilities = Male_Ability(:);
    weights = joint_dist(:);
    
    % Remove zero weights
    non_zero = weights > 1e-10;
    female_abilities = female_abilities(non_zero);
    male_abilities = male_abilities(non_zero);
    weights = weights(non_zero);
    
    fig = figure('Color', 'w', 'Position', [100, 100, 700, 700]);
    
    % Create scatter plot with size proportional to weight
    scatter(female_abilities, male_abilities, weights * 50000, weights, 'filled');
    
    % Add 45-degree line
    hold on;
    min_ability = min([gridP.af.a; gridP.am.a]);
    max_ability = max([gridP.af.a; gridP.am.a]);
    plot([min_ability, max_ability], [min_ability, max_ability], 'r--', 'LineWidth', 2);
    
    xlabel('Female Productivity', 'Interpreter', 'latex', 'FontSize', 18, 'FontName', 'Times');
    ylabel('Male Productivity', 'Interpreter', 'latex', 'FontSize', 18, 'FontName', 'Times');
    title('Assortative Mating by Productivity', 'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Times');
    
    colorbar;
    colormap(parula);
    grid on;
    axis equal;
    xlim([min_ability, max_ability]);
    ylim([min_ability, max_ability]);
    
    legend('Married Couples', '45$^{\circ}$ Line', 'Location', 'southeast', 'Interpreter', 'latex', 'FontSize', 12);
    
    saveas(fig, 'Plots/Summary_AssortativeMating.png');
end

function createWageCorrelationPlot(solutionSP, gridP, stats)
    % Create scatter plot showing wage correlation in couples
    
    % Calculate wages for employed couples only
    wm_f = exp(gridP.af.aa) .* (gridP.wf.ww);
    wm_m = exp(gridP.am.aa) .* (gridP.wm.ww);
    
    % Get employed couple distributions
    employed_couples = solutionSP.eel + solutionSP.eec;  % Both employed
    
    [af_idx, am_idx, wf_idx, wm_idx] = ind2sub(size(employed_couples), find(employed_couples > 1e-10));
    
    female_wages = [];
    male_wages = [];
    couple_weights = [];
    
    for i = 1:length(af_idx)
        female_wages(end+1) = wm_f(af_idx(i), wf_idx(i));
        male_wages(end+1) = wm_m(am_idx(i), wm_idx(i));
        couple_weights(end+1) = employed_couples(af_idx(i), am_idx(i), wf_idx(i), wm_idx(i));
    end
    
    if ~isempty(female_wages)
        fig = figure('Color', 'w', 'Position', [100, 100, 700, 700]);
        
        scatter(female_wages, male_wages, couple_weights * 1000, 'filled');
        
        % Add trend line
        hold on;
        p = polyfit(female_wages, male_wages, 1);
        x_trend = linspace(min(female_wages), max(female_wages), 100);
        y_trend = polyval(p, x_trend);
        plot(x_trend, y_trend, 'r-', 'LineWidth', 2);
        
        % Calculate correlation
        correlation = corr(female_wages', male_wages');
        
        xlabel('Female Wage', 'Interpreter', 'latex', 'FontSize', 18, 'FontName', 'Times');
        ylabel('Male Wage', 'Interpreter', 'latex', 'FontSize', 18, 'FontName', 'Times');
        title(sprintf('Wage Correlation in Married Couples ($\\rho = %.2f$)', correlation), ...
              'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Times');
        
        grid on;
        
        saveas(fig, 'Plots/Summary_WageCorrelation.png');
    end
end

function createJobSearchIntensityPlot(stats)
    % Create line plot showing job search differences by marital status
    
    fig = figure('Color', 'w', 'Position', [100, 100, 800, 600]);
    
    % Simulate job finding rates (would need actual transition matrices for precise calculation)
    marital_status = {'Single', 'Married'};
    
    % Use unemployment rates as proxy for job search intensity (lower = more intense search)
    female_rates = [stats.unemployment_rates.single_f, stats.unemployment_rates.married_f];
    male_rates = [stats.unemployment_rates.single_m, stats.unemployment_rates.married_m];
    
    plot(1:2, female_rates, 'o-', 'LineWidth', 3, 'MarkerSize', 8, 'Color', [0.8, 0.3, 0.3]);
    hold on;
    plot(1:2, male_rates, 's-', 'LineWidth', 3, 'MarkerSize', 8, 'Color', [0.3, 0.3, 0.8]);
    
    set(gca, 'XTickLabel', marital_status, 'FontSize', 14, 'FontName', 'Times');
    ylabel('Unemployment Rate (\%)', 'Interpreter', 'latex', 'FontSize', 18, 'FontName', 'Times');
    title('Job Search Outcomes by Marital Status', 'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Times');
    
    legend('Women', 'Men', 'Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 14);
    grid on;
    
    % Add value labels
    for i = 1:2
        text(i, female_rates(i) + 2, sprintf('%.1f', female_rates(i)), ...
             'HorizontalAlignment', 'center', 'Color', [0.8, 0.3, 0.3], 'FontSize', 11);
        text(i, male_rates(i) + 2, sprintf('%.1f', male_rates(i)), ...
             'HorizontalAlignment', 'center', 'Color', [0.3, 0.3, 0.8], 'FontSize', 11);
    end
    
    saveas(fig, 'Plots/Summary_JobSearchIntensity.png');
end

function createTransitionComparisonPlot(stats)
    % Create grouped bar chart for employment transitions
    
    fig = figure('Color', 'w', 'Position', [100, 100, 800, 600]);
    
    % Use model output values for transition rates
    UE_rates = [19.01, 3.37];  % Female, Male UE rates from model output
    EU_rates = [1.64, 2.34];   % Female, Male EU rates from model output
    
    x = 1:2;
    width = 0.35;
    
    b1 = bar(x - width/2, UE_rates, width, 'FaceColor', [0.3, 0.7, 0.3]);
    hold on;
    b2 = bar(x + width/2, EU_rates, width, 'FaceColor', [0.8, 0.3, 0.3]);
    
    set(gca, 'XTickLabel', {'Women', 'Men'}, 'FontSize', 14, 'FontName', 'Times');
    ylabel('Monthly Transition Rate (\%)', 'Interpreter', 'latex', 'FontSize', 18, 'FontName', 'Times');
    title('Labor Market Transition Rates by Gender', 'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Times');
    legend('Unemployment $\rightarrow$ Employment', 'Employment $\rightarrow$ Unemployment', ...
           'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 14);
    grid on;
    
    % Add value labels on bars
    for i = 1:length(UE_rates)
        text(i - width/2, UE_rates(i) + 0.5, sprintf('%.1f', UE_rates(i)), ...
             'HorizontalAlignment', 'center', 'FontSize', 12);
        text(i + width/2, EU_rates(i) + 0.5, sprintf('%.1f', EU_rates(i)), ...
             'HorizontalAlignment', 'center', 'FontSize', 12);
    end
    
    saveas(fig, 'Plots/Summary_TransitionRates.png');
end

function createWageProgressionPlot(solutionSP, gridSf, gridSm, gridP)
    % Create line plot showing wage progression by marital status
    
    try
        % Use single grids for productivity levels
        productivity_levels_f = gridSf.a.a;
        productivity_levels_m = gridSm.a.a;
        
        % Single wages
        ws_f = exp(gridSf.a.aa) .* (gridSf.w.ww);
        ws_m = exp(gridSm.a.aa) .* (gridSm.w.ww);
        
        single_f_avg_wages = zeros(size(productivity_levels_f));
        single_m_avg_wages = zeros(size(productivity_levels_m));
        
        for i = 1:length(productivity_levels_f)
            single_f_avg_wages(i) = mean(ws_f(i, :));
        end
        
        for i = 1:length(productivity_levels_m)
            single_m_avg_wages(i) = mean(ws_m(i, :));
        end
        
        % Married wages - use couple grid dimensions
        wm_f = exp(gridP.af.aa) .* (gridP.wf.ww);
        wm_m = exp(gridP.am.aa) .* (gridP.wm.ww);
        
        % Get dimensions of couple grids
        nAf = size(gridP.af.a, 1);
        nAm = size(gridP.am.a, 1);
        
        married_f_avg_wages = zeros(nAf, 1);
        married_m_avg_wages = zeros(nAm, 1);
        
        for i = 1:nAf
            % Average across all male partner types and wage states for female i
            married_f_avg_wages(i) = mean(wm_f(i, :, :), 'all');
        end
        
        for i = 1:nAm
            % Average across all female partner types and wage states for male i
            married_m_avg_wages(i) = mean(wm_m(:, i, :), 'all');
        end
        
        fig = figure('Color', 'w', 'Position', [100, 100, 800, 600]);
        
        plot(productivity_levels_f, single_f_avg_wages, 'o-', 'LineWidth', 3, 'MarkerSize', 6, 'Color', [0.8, 0.3, 0.3]);
        hold on;
        plot(gridP.af.a, married_f_avg_wages, 's--', 'LineWidth', 3, 'MarkerSize', 6, 'Color', [0.8, 0.3, 0.3]);
        plot(productivity_levels_m, single_m_avg_wages, 'o-', 'LineWidth', 3, 'MarkerSize', 6, 'Color', [0.3, 0.3, 0.8]);
        plot(gridP.am.a, married_m_avg_wages, 's--', 'LineWidth', 3, 'MarkerSize', 6, 'Color', [0.3, 0.3, 0.8]);
        
        xlabel('Productivity Level', 'Interpreter', 'latex', 'FontSize', 18, 'FontName', 'Times');
        ylabel('Average Wage', 'Interpreter', 'latex', 'FontSize', 18, 'FontName', 'Times');
        title('Wage Progression by Productivity and Marital Status', 'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Times');
        
        legend('Single Women', 'Married Women', 'Single Men', 'Married Men', ...
               'Location', 'northwest', 'Interpreter', 'latex', 'FontSize', 12);
        grid on;
        
        saveas(fig, 'Plots/Summary_WageProgression.png');
        
    catch ME
        warning('Wage progression plot failed: %s', ME.message);
        fprintf('Skipping wage progression plot due to error.\n');
    end
end