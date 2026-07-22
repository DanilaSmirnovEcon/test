function welfaregrid=searchPolicyGrid(bcvec, bmvec, solutionSPdist, param, gridS, gridP)

% this function solves the problem on the given policy grid and returns the
% set of corresponding welfare
disp("Search policy on the grid");

%% calculate ability-equalizing transfers

% Aggregate income
solutionSP=solutionSPdist;

welfareSS = sum(solutionSPdist.e.*solutionSPdist.JValue, 'all')+...
                      sum(solutionSPdist.h.*solutionSPdist.UValue, 'all')+...
                      sum(solutionSPdist.ec.*solutionSPdist.JcValue, 'all')+...
                      sum(solutionSPdist.hc.*solutionSPdist.UcValue, 'all')+...
                      sum(solutionSPdist.ee.*solutionSPdist.JJValue, 'all')+...
                      sum(solutionSPdist.eh.*solutionSPdist.JUValue, 'all')+...
                      sum(solutionSPdist.he.*solutionSPdist.UJValue, 'all')+...
                      sum(solutionSPdist.hh.*solutionSPdist.UUValue, 'all')+...
                      sum(solutionSPdist.eec.*solutionSPdist.JJcValue, 'all')+...
                      sum(solutionSPdist.ehc.*solutionSPdist.JUcValue, 'all')+...
                      sum(solutionSPdist.hec.*solutionSPdist.UJcValue, 'all')+...
                      sum(solutionSPdist.hhc.*solutionSPdist.UUcValue, 'all')+...
                      sum(solutionSPdist.fOut.*solutionSPdist.constr, 'all');

welfaregrid = ones(length(bcvec), length(bmvec)).*welfareSS.*NaN;
                       
param_s=param;
param_s.Sdimmult=1;
[param_s, gridS_s, ~] = param_bprog(param_s);
solutionS_s = SolveSingleProblem(param_s, gridS_s);
solutionS_s.h = solutionS_s.h+solutionS_s.o;
solutionS_s.hc = solutionS_s.hc+solutionS_s.oc;

%% calculate policy grid
bprogr=1/1;
param.Ball=0;
% for i = 1:length(bcvec)
%     param.bc = bcvec(i);
for j = 1:length(bmvec)
      param.bs = bmvec(j);
     
%       param.Ball = bmvec(j);
%       param.b = @(x) param.bnew.*exp(x.*bprogr) + param.Ball;
      
%     param.updMult  = 0.8;
%     solutionS = SolveSingleGE(0.02, param, gridS);
%     solutionSO = SolveSingleOutProblem(param, gridP);
%     solutionP = SolveDoubleProblem(solutionSO, false, true, param, gridP);%load, display    
%     param.updMult  = 0.99;
%     solutionSP = SolveEndogProblem(solutionS, solutionSO, solutionP, false, true, param, gridS, gridP);              
    
%     for j = 1:length(bmvec)
%         param.bs = bmvec(j);  
    for i = 1:length(bcvec)
        param.bc = bcvec(i);

        param.nokid = @(x, a) x;
        param.onekid = @(x, a) x - param.c;%.*(x+1)./2;%max(x-param.c, 0.01);
        param.onekidU = @(x, a) x + param.bc - param.c;%.*(x+1)./2;
        
        %param.b = @(x) param.bnew.*exp(x.*bprogr) + param.bc;
        %param.onekidU = @(x, a) x+param.bc.*exp(a)-param.c;%.*(x+1)./2;
        try
            for k=1:6
                [solutionSP, diff, param]=fingEqB(solutionSPdist, solutionSP, solutionS_s, param, gridS, gridP);
                disp(num2str(diff));
                if diff<0.001
                    break
                end
            end
            if diff>0.01
                error("err")
            end
            disp("delta benefit kid = " + num2str(param.bc) +...
                 "; delta benefit single = " + num2str(param.bs) +...
                 "; transfer unemployed = " + num2str(param.Ball));

            welfare = sum(solutionSP.e.*solutionSP.JValue, 'all')+...
                      sum(solutionSP.h.*solutionSP.UValue, 'all')+...
                      sum(solutionSP.ec.*solutionSP.JcValue, 'all')+...
                      sum(solutionSP.hc.*solutionSP.UcValue, 'all')+...
                      sum(solutionSP.ee.*solutionSP.JJValue, 'all')+...
                      sum(solutionSP.eh.*solutionSP.JUValue, 'all')+...
                      sum(solutionSP.he.*solutionSP.UJValue, 'all')+...
                      sum(solutionSP.hh.*solutionSP.UUValue, 'all')+...
                      sum(solutionSP.eec.*solutionSP.JJcValue, 'all')+...
                      sum(solutionSP.ehc.*solutionSP.JUcValue, 'all')+...
                      sum(solutionSP.hec.*solutionSP.UJcValue, 'all')+...
                      sum(solutionSP.hhc.*solutionSP.UUcValue, 'all')+...
                      sum(solutionSP.fOut.*solutionSP.constr, 'all');
                  
            if welfare>0
                error("err")
            end
            welfaregrid(i,j) = welfare;
            
        catch
            welfare = NaN;
            welfaregrid(i,j) = welfare;

        end
        disp("Welfare = " + num2str(welfare))
        disp(num2str(round((length(bmvec)*(i-1)+j)/length(bcvec)/length(bmvec)*100)) + "% done");
        disp("")     
        
        figure(1)                          
        heatmap(bmvec, bcvec, (welfaregrid-welfareSS)./abs(welfareSS).*100);
        title('Welfare');
        ylabel('kid benefit');
        xlabel('single benefit');
    end
end

disp(welfaregrid);

%save('WelfareGridMarrKid','welfaregrid','bcvec','bmvec','param','gridS','gridP');


                  
fig = heatmap(bmvec, bcvec, (welfaregrid-welfareSS)./abs(welfareSS));
title('Welfare');
ylabel('kid benefit');
xlabel('pair benefit');
saveas(fig, 'Plots/WelfarePolicyMarrKid.png');
end

function [solutionSP, diff, param]=fingEqB(solutionSPdist, solutionSP, solutionS_s, param, gridS, gridP)
    tau=0.02;
    diff=1;
    bprogr=1;        
            
    try      
        solutionS = SolveSingleProblem(param, gridS);
        %[solutionS, param] = SolveSingleGE(0.02, param, gridS, false, Ball);
        solutionSO = SolveSingleOutProblem(param, gridP);
        solutionP = SolveDoubleProblem(solutionSO, true, false, param, gridP);%load, display
        save('NoAsolutionP','solutionP');
        solutionSP = SolveEndogProblemDistFixed(solutionSPdist, solutionS, solutionS_s, solutionSO, solutionP, false, false, param, gridS, gridP);       
        %save('solutionSP','solutionSP');
        %solutionSP = SolveEndogProblem(solutionS, solutionSO, solutionP, false, true, param, gridS, gridP); %load, display
        
%         unempbasem = sum((exp(gridP.a1.a.*bprogr).*(sum(((solutionSP.hh+solutionSP.hhc)), [2])+...
%                                                   sum(((solutionSP.ho+solutionSP.hoc)), [2])+...
%                                                   sum(((solutionSP.he+solutionSP.hec)), [2,3]))),'all');

%         unempbasec = sum((exp(gridP.a1.a.*bprogr).*(sum(((solutionSP.hhc)), [2])+...
%                                                   sum(((solutionSP.hoc)), [2])+...
%                                                   sum(((solutionSP.hec)), [2,3]))),'all')+...
%                      sum((exp(gridS.a.a.*bprogr).*(solutionSP.hc)),'all');
        unempbasec = sum(((sum(((solutionSP.hhc)), [2])+...
                                              sum(((solutionSP.hoc)), [2])+...
                                              sum(((solutionSP.hec)), [2,3]))),'all')+...
                 sum(((solutionSP.hc)),'all');
        
        basec = sum(((sum(((solutionSP.hhc)), 'all')+...
                      sum(((solutionSP.hoc)), 'all')+...
                      sum(((solutionSP.hec)), 'all'))),'all')+...
                      sum(((solutionSP.hc)),'all')+...
                sum(((sum(((solutionSP.ehc)), 'all')+...
                      sum(((solutionSP.eoc)), 'all')+...
                      sum(((solutionSP.eec)), 'all'))),'all')+...
                      sum(((solutionSP.ec)),'all')+...
                sum(solutionSP.omc,'all');
        %basec = sum(((solutionSP.hc)),'all')+sum(((solutionSP.ec)),'all');
                 
        unempbase = sum((exp(gridP.a1.a.*bprogr).*(sum(((solutionSP.hh+solutionSP.hhc)), [2])+...
                                                  sum(((solutionSP.ho+solutionSP.hoc)), [2])+...
                                                  sum(((solutionSP.he+solutionSP.hec)), [2,3]))),'all')+...
                     sum((exp(gridS.a.a.*bprogr).*(solutionSP.h+solutionSP.hc)),'all');
                 
        unemp = sum(((sum(((solutionSP.hh+solutionSP.hhc)), [2])+...
                                                  sum(((solutionSP.ho+solutionSP.hoc)), [2])+...
                                                  sum(((solutionSP.he+solutionSP.hec)), [2,3]))),'all')+...
                     sum(((solutionSP.h+solutionSP.hc)),'all');
                 
        unemps = sum(((solutionSP.h+solutionSP.hc)),'all');
        
        w=(exp(gridP.a1.aa).*(gridP.r1.rr));
        wh = squeeze(w(:,:,:,1));
        wage_s= sum((exp(gridS.a.aa).*(gridS.r.rr)).*(solutionSP.e+solutionSP.ec),'all')+...
                sum(w.*(solutionSP.ee+solutionSP.eec),'all')+...
                sum(wh.*((solutionSP.eh+solutionSP.ehc)+(solutionSP.eo+solutionSP.eoc)),'all');
        emp= sum((solutionSP.e+solutionSP.ec),'all')+...
                sum((solutionSP.ee+solutionSP.eec),'all')+...
                sum(((solutionSP.eh+solutionSP.ehc)+(solutionSP.eo+solutionSP.eoc)),'all');

        B = wage_s/(1-tau)*tau;

        %bnew = (B-param.bm*unempbasem-param.bc*basec)/unempbase;  
        %tnew = (B - param.bs*unemps - param.bc*unempbasec - param.bnew.*unempbase)/emp;
%         tnew = (B - param.Ball*unemp - param.bnew.*unempbase)/emp;       
%         diff=abs(param.t-tnew);
%         param.t = 0.85*(tnew)+(1-0.85)*param.t;
        Bnew = (B - param.bs*unemps - param.bc*unempbasec - param.bnew.*unempbase)/unemp; 
        Bnew = ( - param.bs*unemps - param.bc*unempbasec)/unemp; 
        diff=abs(param.Ball-Bnew);
        param.Ball = 0.85*(Bnew)+(1-0.85)*param.Ball;
        param.b = @(x) param.bnew.*exp(x.*bprogr) + param.Ball;
        
%         diff=abs(param.bnew-bnew);
%         
%         bnew = 0.15*(param.bnew)+(1-0.15)*bnew;
%         param.bnew=bnew;
%         param.b = @(x) param.bnew.*exp(x.*bprogr);
    catch 
        solutionSP=solutionSPdist;
        disp("")
        error("err")
    end
end

