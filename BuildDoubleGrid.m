function [param, gridP] = BuildDoubleGrid(param, buildC)
%==========================================================================
% BUILDDOUBLEGRID - Construct state space grids for couples
%==========================================================================
% This function builds the discretized state space grids for married couples
% and constructs transition matrices for the joint stochastic productivity
% process. The state space is 4-dimensional: (female productivity, 
% female wage, male productivity, male wage).
%
% Inputs:
%   param  - Main parameter structure
%   buildC - Boolean flag to build transition matrices (computationally intensive)
%
% Outputs:
%   param  - Updated parameter structure with couple grid information
%   gridP  - State space grids for couples
%
% The function creates:
% - 4D grids for couple productivity and wage states
% - Transition matrices for joint productivity evolution
% - Index mappings for different employment combinations
%
% Reference: Morazzoni & Smirnov (2025)
%==========================================================================

%--------------------------------------------------------------------------
% 1. GRID DIMENSIONS AND SETUP
%--------------------------------------------------------------------------

% Set grid dimensions (same as base parameter grids for computational efficiency)
gridP.af.n = param.a.n;     % Female productivity grid points
gridP.wf.n = param.w.n;     % Female wage offer grid points  
gridP.am.n = param.a.n;     % Male productivity grid points
gridP.wm.n = param.w.n;     % Male wage offer grid points

%--------------------------------------------------------------------------
% 2. FEMALE PRODUCTIVITY AND WAGE GRIDS
%--------------------------------------------------------------------------

% Female productivity grid
gridP.af.min = param.a.min;
gridP.af.max = param.a.max;
gridP.af.a = linspace(gridP.af.min, gridP.af.max, gridP.af.n)';

% Create 4D female productivity array: (af, wf, am, wm)
aa_female = repmat(gridP.af.a, 1, gridP.wf.n, gridP.am.n, gridP.wm.n);
gridP.af.aa = permute(aa_female, [1,2,3,4]);
gridP.af.da = (gridP.af.max - gridP.af.min) / (gridP.af.n - 1);
gridP.af.da2 = gridP.af.da^2;

% Female wage offer grid
gridP.wf.min = param.w.min;
gridP.wf.max = param.w.max;
gridP.wf.w = linspace(gridP.wf.min, gridP.wf.max, gridP.wf.n)';

% Create 4D female wage array: (af, wf, am, wm)
ww_female = repmat(gridP.wf.w, 1, gridP.af.n, gridP.am.n, gridP.wm.n);
gridP.wf.ww = permute(ww_female, [2,1,3,4]);
gridP.wf.dw = (gridP.wf.max - gridP.wf.min) / (gridP.wf.n - 1);

% Female productivity process parameters
gridP.af.mu = param.pf.mu * (param.pf.amean - gridP.af.aa);     % Employed drift
gridP.af.mu_u = param.pf.mu_u * (param.a.min - gridP.af.aa);   % Unemployed drift
gridP.af.var = (param.pf.a_sigma^2) / 2;                       % Diffusion coefficient

%--------------------------------------------------------------------------
% 3. MALE PRODUCTIVITY AND WAGE GRIDS
%--------------------------------------------------------------------------

% Male productivity grid
gridP.am.min = param.a.min;
gridP.am.max = param.a.max;
gridP.am.a = linspace(gridP.am.min, gridP.am.max, gridP.am.n)';

% Create 4D male productivity array: (af, wf, am, wm)
aa_male = repmat(gridP.am.a, 1, gridP.af.n, gridP.wf.n, gridP.wm.n);
gridP.am.aa = permute(aa_male, [2,3,1,4]);
gridP.am.da = (gridP.am.max - gridP.am.min) / (gridP.am.n - 1);
gridP.am.da2 = gridP.am.da^2;

% Male wage offer grid
gridP.wm.min = param.w.min;
gridP.wm.max = param.w.max;
gridP.wm.w = linspace(gridP.wm.min, gridP.wm.max, gridP.wm.n)';

% Create 4D male wage array: (af, wf, am, wm)
ww_male = repmat(gridP.wm.w, 1, gridP.af.n, gridP.wf.n, gridP.am.n);
gridP.wm.ww = permute(ww_male, [2,3,4,1]);
gridP.wm.dw = (gridP.wm.max - gridP.wm.min) / (gridP.wm.n - 1);

% Male productivity process parameters
gridP.am.mu = param.pm.mu * (param.pm.amean - gridP.am.aa);     % Employed drift
gridP.am.mu_u = param.pm.mu_u * (param.a.min - gridP.am.aa);   % Unemployed drift
gridP.am.var = (param.pm.a_sigma^2) / 2;                       % Diffusion coefficient

%--------------------------------------------------------------------------
% 4. ERGODIC DISTRIBUTIONS FOR NEW COUPLES
%--------------------------------------------------------------------------

% Female entry distribution (used when forming new couples)
param.pf.a.Perg = lognpdf(gridP.af.a + gridP.af.da, param.pf.aergmean, ...
                          sqrt(param.pf.aergvar)); 
param.pf.a.Perg = param.pf.a.Perg ./ sum(param.pf.a.Perg);

% Male entry distribution (used when forming new couples)
param.pm.a.Perg = lognpdf(gridP.am.a + gridP.am.da, param.pm.aergmean, ...
                          sqrt(param.pm.aergvar)); 
param.pm.a.Perg = param.pm.a.Perg ./ sum(param.pm.a.Perg);

%--------------------------------------------------------------------------
% 5. STATE SPACE DIMENSIONS
%--------------------------------------------------------------------------

% Total number of states for different employment combinations
param.NP = gridP.af.n * gridP.wf.n * gridP.am.n * gridP.wm.n;  % Both employed (JJ)
param.NPall = param.NP + param.NP/gridP.wm.n + param.NP/gridP.wf.n + ...
              gridP.af.n * gridP.am.n;  % All employment combinations

%--------------------------------------------------------------------------
% 6. TRANSITION MATRIX CONSTRUCTION (OPTIONAL)
%--------------------------------------------------------------------------

% Build transition matrices for couple productivity evolution (computationally intensive)
if buildC
    % fprintf('   - Building couple transition matrices (this may take time)...\n');
    
    % Initialize sparse transition matrices for all employment combinations
    param.CJJ = sparse(zeros(param.NP, param.NP));                    % Both employed
    param.CUU = sparse(zeros(gridP.af.n*gridP.am.n, gridP.af.n*gridP.am.n)); % Both unemployed
    param.CJU = sparse(zeros(param.NP/gridP.wm.n, param.NP/gridP.wm.n));     % Female employed, male unemployed
    param.CUJ = sparse(zeros(param.NP/gridP.wf.n, param.NP/gridP.wf.n));     % Female unemployed, male employed
    
    % Job separation transition matrices  
    param.CJU_UU = sparse(zeros(param.NP/gridP.wm.n, gridP.af.n*gridP.am.n));
    param.CUJ_UU = sparse(zeros(param.NP/gridP.wf.n, gridP.af.n*gridP.am.n));
    param.CJJ_UJ = sparse(zeros(param.NP, param.NP/gridP.wf.n));
    param.CJJ_JU = sparse(zeros(param.NP, param.NP/gridP.wm.n));
    
    % Build finite difference transition matrices using original working logic
    for afi = 1:gridP.af.n
        for ami = 1:gridP.am.n    
            
            %--------------------------------------------------------------
            % Case 1: Both unemployed (UU state)
            %--------------------------------------------------------------
            i = (ami-1)*gridP.af.n + afi;
            
            % Female productivity evolution (unemployed drift toward minimum)
            sb = min(0, gridP.af.mu_u(afi,1,ami,1));
            sf = max(0, gridP.af.mu_u(afi,1,ami,1));
            runninglength = 1;
            ss = gridP.af.var;        
            if (gridP.af.a(afi) > gridP.af.min)   
                param.CUU(i, i-runninglength) = param.CUU(i, i-runninglength) + ...
                                               (-sb/gridP.af.da + ss/gridP.af.da2);
                param.CUU(i, i) = param.CUU(i, i) + ...
                                 (sb/gridP.af.da - ss/gridP.af.da2);
            end
            if (gridP.af.a(afi) < gridP.af.max)
                param.CUU(i, i+runninglength) = param.CUU(i, i+runninglength) + ...
                                               (sf/gridP.af.da + ss/gridP.af.da2);
                param.CUU(i, i) = param.CUU(i, i) + ...
                                 (-sf/gridP.af.da - ss/gridP.af.da2);
            end
            
            % Male productivity evolution (unemployed drift toward minimum)
            sb = min(0, gridP.am.mu_u(afi,1,ami,1));
            sf = max(0, gridP.am.mu_u(afi,1,ami,1));
            runninglength = gridP.af.n;
            ss = gridP.am.var;
            if (gridP.am.a(ami) > gridP.am.min)   
                param.CUU(i, i-runninglength) = param.CUU(i, i-runninglength) + ...
                                               (-sb/gridP.am.da + ss/gridP.am.da2);
                param.CUU(i, i) = param.CUU(i, i) + ...
                                 (sb/gridP.am.da - ss/gridP.am.da2);
            end
            if (gridP.am.a(ami) < gridP.am.max)
                param.CUU(i, i+runninglength) = param.CUU(i, i+runninglength) + ...
                                               (sf/gridP.am.da + ss/gridP.am.da2);
                param.CUU(i, i) = param.CUU(i, i) + ...
                                 (-sf/gridP.am.da - ss/gridP.am.da2);
            end
            param.CUU(i, i) = param.CUU(i, i) - param.zeta;  % Retirement

            %--------------------------------------------------------------
            % Case 2: Female employed, male unemployed (JU state)
            %--------------------------------------------------------------
            for wfi = 1:gridP.wf.n
                i = (ami-1)*gridP.af.n*gridP.wf.n + (wfi-1)*gridP.af.n + afi;

                % Female productivity evolution (employed drift toward mean)
                runninglength = 1;
                sb = min(0, gridP.af.mu(afi,wfi,ami,1));
                sf = max(0, gridP.af.mu(afi,wfi,ami,1));
                ss = gridP.af.var;            
                if (gridP.af.a(afi) > gridP.af.min)   
                    param.CJU(i, i-runninglength) = param.CJU(i, i-runninglength) + ...
                                                   (-sb/gridP.af.da + ss/gridP.af.da2);
                    param.CJU(i, i) = param.CJU(i, i) + ...
                                     (sb/gridP.af.da - ss/gridP.af.da2);
                end
                if (gridP.af.a(afi) < gridP.af.max)
                    param.CJU(i, i+runninglength) = param.CJU(i, i+runninglength) + ...
                                                   (sf/gridP.af.da + ss/gridP.af.da2);
                    param.CJU(i, i) = param.CJU(i, i) + ...
                                     (-sf/gridP.af.da - ss/gridP.af.da2);
                end

                % Male productivity evolution (unemployed drift)
                sb = min(0, gridP.am.mu_u(afi,wfi,ami,1));
                sf = max(0, gridP.am.mu_u(afi,wfi,ami,1));
                runninglength = gridP.af.n*gridP.wf.n;
                ss = gridP.am.var;
                if (gridP.am.a(ami) > gridP.am.min)   
                    param.CJU(i, i-runninglength) = param.CJU(i, i-runninglength) + ...
                                                   (-sb/gridP.am.da + ss/gridP.am.da2);
                    param.CJU(i, i) = param.CJU(i, i) + ...
                                     (sb/gridP.am.da - ss/gridP.am.da2);
                end
                if (gridP.am.a(ami) < gridP.am.max)
                    param.CJU(i, i+runninglength) = param.CJU(i, i+runninglength) + ...
                                                   (sf/gridP.am.da + ss/gridP.am.da2);
                    param.CJU(i, i) = param.CJU(i, i) + ...
                                     (-sf/gridP.am.da - ss/gridP.am.da2);
                end

                param.CJU(i, i) = param.CJU(i, i) - param.zeta;  % Retirement
                param.CJU(i, i) = param.CJU(i, i) - param.pf.delta;  % Female job loss
                param.CJU_UU(i, (ami-1)*gridP.af.n+afi) = param.CJU_UU(i, (ami-1)*gridP.af.n+afi) + param.pf.delta;   
            end

            %--------------------------------------------------------------
            % Case 3: Female unemployed, male employed (UJ state)  
            %--------------------------------------------------------------
            for wmi = 1:gridP.wm.n
                i = (wmi-1)*gridP.af.n*gridP.am.n + (ami-1)*gridP.af.n + afi;            

                % Female productivity evolution (unemployed drift)
                runninglength = 1;
                sb = min(0, gridP.af.mu_u(afi,1,ami,wmi));
                sf = max(0, gridP.af.mu_u(afi,1,ami,wmi));
                ss = gridP.af.var;
                if (gridP.af.a(afi) > gridP.af.min)   
                    param.CUJ(i, i-runninglength) = param.CUJ(i, i-runninglength) + ...
                                                   (-sb/gridP.af.da + ss/gridP.af.da2);
                    param.CUJ(i, i) = param.CUJ(i, i) + ...
                                     (sb/gridP.af.da - ss/gridP.af.da2);
                end
                if (gridP.af.a(afi) < gridP.af.max)
                    param.CUJ(i, i+runninglength) = param.CUJ(i, i+runninglength) + ...
                                                   (sf/gridP.af.da + ss/gridP.af.da2);
                    param.CUJ(i, i) = param.CUJ(i, i) + ...
                                     (-sf/gridP.af.da - ss/gridP.af.da2);
                end

                % Male productivity evolution (employed drift)
                runninglength = gridP.af.n;
                sb = min(0, gridP.am.mu(afi,1,ami,wmi));
                sf = max(0, gridP.am.mu(afi,1,ami,wmi));
                ss = gridP.am.var;            
                if (gridP.am.a(ami) > gridP.am.min)   
                    param.CUJ(i, i-runninglength) = param.CUJ(i, i-runninglength) + ...
                                                   (-sb/gridP.am.da + ss/gridP.am.da2);
                    param.CUJ(i, i) = param.CUJ(i, i) + ...
                                     (sb/gridP.am.da - ss/gridP.am.da2);
                end
                if (gridP.am.a(ami) < gridP.am.max)
                    param.CUJ(i, i+runninglength) = param.CUJ(i, i+runninglength) + ...
                                                   (sf/gridP.am.da + ss/gridP.am.da2);
                    param.CUJ(i, i) = param.CUJ(i, i) + ...
                                     (-sf/gridP.am.da - ss/gridP.am.da2);
                end

                param.CUJ(i, i) = param.CUJ(i, i) - param.zeta;  % Retirement
                param.CUJ(i, i) = param.CUJ(i, i) - param.pm.delta;  % Male job loss
                param.CUJ_UU(i, (ami-1)*gridP.af.n+afi) = param.CUJ_UU(i, (ami-1)*gridP.af.n+afi) + param.pm.delta;   
            end

            %--------------------------------------------------------------
            % Case 4: Both employed (JJ state)
            %--------------------------------------------------------------
            for wfi = 1:gridP.wf.n
                for wmi = 1:gridP.wm.n
                    i = (wmi-1)*gridP.am.n*gridP.af.n*gridP.wf.n + (ami-1)*gridP.af.n*gridP.wf.n + (wfi-1)*gridP.af.n + afi;        
                    
                    % Female productivity evolution (employed drift)
                    runninglength = 1;
                    sb = min(0, gridP.af.mu(afi,wfi,ami,wmi));
                    sf = max(0, gridP.af.mu(afi,wfi,ami,wmi));
                    ss = gridP.af.var;            
                    if (gridP.af.a(afi) > gridP.af.min)   
                        param.CJJ(i, i-runninglength) = param.CJJ(i, i-runninglength) + ...
                                                       (-sb/gridP.af.da + ss/gridP.af.da2);
                        param.CJJ(i, i) = param.CJJ(i, i) + ...
                                         (sb/gridP.af.da - ss/gridP.af.da2);
                    end
                    if (gridP.af.a(afi) < gridP.af.max)
                        param.CJJ(i, i+runninglength) = param.CJJ(i, i+runninglength) + ...
                                                       (sf/gridP.af.da + ss/gridP.af.da2);
                        param.CJJ(i, i) = param.CJJ(i, i) + ...
                                         (-sf/gridP.af.da - ss/gridP.af.da2);
                    end
                    
                    % Male productivity evolution (employed drift)
                    runninglength = gridP.af.n*gridP.wf.n;
                    sb = min(0, gridP.am.mu(afi,wfi,ami,wmi));
                    sf = max(0, gridP.am.mu(afi,wfi,ami,wmi));
                    ss = gridP.am.var;            
                    if (gridP.am.a(ami) > gridP.am.min)   
                        param.CJJ(i, i-runninglength) = param.CJJ(i, i-runninglength) + ...
                                                       (-sb/gridP.am.da + ss/gridP.am.da2);
                        param.CJJ(i, i) = param.CJJ(i, i) + ...
                                         (sb/gridP.am.da - ss/gridP.am.da2);
                    end
                    if (gridP.am.a(ami) < gridP.am.max)
                        param.CJJ(i, i+runninglength) = param.CJJ(i, i+runninglength) + ...
                                                       (sf/gridP.am.da + ss/gridP.am.da2);
                        param.CJJ(i, i) = param.CJJ(i, i) + ...
                                         (-sf/gridP.am.da - ss/gridP.am.da2);
                    end

                    param.CJJ(i, i) = param.CJJ(i, i) - param.zeta;  % Retirement
                    param.CJJ(i, i) = param.CJJ(i, i) - param.pf.delta - param.pm.delta;  % Job losses
                    param.CJJ_JU(i, (ami-1)*gridP.af.n*gridP.wf.n+(wfi-1)*gridP.af.n+afi) = ...
                        param.CJJ_JU(i, (ami-1)*gridP.af.n*gridP.wf.n+(wfi-1)*gridP.af.n+afi) + param.pm.delta;   
                    param.CJJ_UJ(i, (wmi-1)*gridP.af.n*gridP.am.n+(ami-1)*gridP.af.n+afi) = ...
                        param.CJJ_UJ(i, (wmi-1)*gridP.af.n*gridP.am.n+(ami-1)*gridP.af.n+afi) + param.pf.delta;   
                end
            end                
        end
    end
    
    % fprintf('   - Couple transition matrices completed\n');
end

end