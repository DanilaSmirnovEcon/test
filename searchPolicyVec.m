function welfarevec=searchPolicyVec(bvec, solutionSPdist)
global param gridS gridP
% this function solves the problem on the given policy vector and returns the
% set of corresponding welfare
disp("Search policy on the grid");
welfarevec = zeros(1,length(bvec));

bbar = param.bnew;
B = bbar*(sum(solutionSPdist.h+solutionSPdist.hc, 'all')+...
                2*sum(solutionSPdist.hh+solutionSPdist.hhc, 'all')+...
                    sum(solutionSPdist.he+solutionSPdist.hec, 'all')+...
                    sum(solutionSPdist.eh+solutionSPdist.ehc, 'all'));

for i = 1:length(bvec)
    param.bm = bvec(i);       
    Tbm = (param.bm)*(2*sum(solutionSPdist.hh+solutionSPdist.hhc, 'all')+...
                    sum(solutionSPdist.he+solutionSPdist.hec, 'all')+...
                    sum(solutionSPdist.eh+solutionSPdist.ehc, 'all'));
    Bres = B - Tbm;
    param.bnew = bbar*Bres/B;
    param.b = @(x) param.bnew;%
    
    disp("Benefit single = " + num2str(param.b(0))+"; Benefit married = " +num2str(param.b(0)+param.bm));
    try
        solutionS = SolveSingleProblem();
        solutionSO = SolveSingleOutProblem();
        solutionP = SolveDoubleProblem(solutionSO, false, false); %load, display
        save('NoAsolutionP','solutionP','param','gridP');
        solutionS.h = solutionS.h+solutionS.o;
        solutionS.hc = solutionS.hc+solutionS.oc;
        %solutionSP = SolveEndogProblem(solutionS, solutionSO, solutionP, false, false); %load, display
        solutionSP = SolveEndogProblemDistFixed(solutionSPdist, solutionS, solutionSO, solutionP, false, false);
        
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
                  sum(solutionSP.hhc.*solutionSP.UUcValue, 'all');
        welfarevec(i) = welfare; 
        %solutionSPdist = solutionSP;
        save('NoAsolutionSP','solutionSP','param','gridS','gridP');           
    catch 
        welfarevec(i) = NaN;
    end
    disp("Welfare = " + num2str(welfare))
    disp("")
end

fig = plot(bvec, welfarevec);
title('Welfare');
ylabel('welfare');
xlabel('benefit');
saveas(fig, 'Plots/WelfarePolicyVec.png');
end