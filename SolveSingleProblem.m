function solutionS = SolveSingleProblem(param, pg, gridS)
% SolveSingleProblem — singles' problem (female or male, based on pg)
%
% Returns: solutionS with value functions, A-blocks, and stationary masses.

    % ---------- housekeeping ----------
    tremble = 0;                 % small >0 only if you need tie smoothing
    tol     = 1e-10;             % tie-breaking tolerance for argmax switches

    % Offer distribution: must be a prob. vector of length gridS.w.n
    if ~isfield(param,'f_long') || isempty(param.f_long)
        error('param.f_long must be provided and sum to 1');
    end
    param.f_long = param.f_long(:);
    s = sum(param.f_long);
    if ~(isfinite(s) && s>0), error('param.f_long must sum to a positive finite value'); end
    param.f_long = param.f_long / s;
    if numel(param.f_long) ~= gridS.w.n
        error('length(param.f_long) (%d) must equal gridS.w.n (%d)', numel(param.f_long), gridS.w.n);
    end

    % Precompute ability × wage matrices
    AA    = repmat(gridS.a.a, 1, gridS.w.n);            % (a × r)
    if ~isequal(size(gridS.w.ww), [gridS.a.n, gridS.w.n])
        error('gridS.w.ww must be size [gridS.a.n × gridS.w.n]');
    end
    W_eff = exp(AA) .* gridS.w.ww;                       % realized earnings on the job (a × r)
    c_a   = param.c(gridS.a.a);                          % child cost vector (a × 1)

    % ---------- initial guesses ----------
    solutionS.JlValue = util(param, param.nokid_tax(W_eff, AA, param));
    solutionS.UlValue = util(param, param.nokid_tax(param.b(gridS.a.a), gridS.a.a, param));
    solutionS.OlValue = util(param, param.no_tax(  pg.bout(gridS.a.a),  gridS.a.a, param));

    solutionS.JcValue = util(param, param.onekid_tax(W_eff, AA, param) - repmat(c_a,1,gridS.w.n));
    solutionS.UcValue = util(param, param.onekid_tax(param.b(gridS.a.a), gridS.a.a, param) - c_a);
    solutionS.OcValue = util(param, param.no_tax(  pg.boutc(gridS.a.a),  gridS.a.a, param) - c_a);

    % ---------- iterate on values ----------
    % diff = Inf;
    % nIter = max(1, param.maxiter*10);
    % for it = 1:nIter
    %     [solutionS, diff] = updateValue(solutionS);
    %     if diff < param.critC, break; end
    % end
    % ---------- Value iteration (actually linear system updates) ----------
    diff = Inf;
    for j = 1:(param.maxiter*10)
        [solutionS, diff] = updateValue(solutionS);
        if diff < param.critC
            % fprintf('Single converged in %d iterations; diff %.3e\n', j, diff);
            break
        end
    end
    if diff > param.critC
        % If you want a hard stop, error here. For replication, just notify.
        % fprintf('SolveSingleProblem: not fully converged; final diff=%.3e\n', diff);
    end
    
    % ---------- assemble full A and stationary distribution ----------
    AL_O = solutionS.AL_O;
    AO_L = solutionS.AO_L;
    Aall = [solutionS.AL, AL_O; AO_L, solutionS.AO];
    solutionS.Aall = Aall;

    % Stationary mass: solve Aall' f = rhs with normalization ∑f = 1
    nVar = size(Aall,1);
    rhs  = zeros(nVar+1,1);

    % Optional demographic injection into Ul rows (if fields exist)
    if isfield(param,'zeta') && isfield(pg,'a') && isfield(pg.a,'erg') ...
       && numel(pg.a.erg)==gridS.a.n
        ul_lo =  (pg.NS) + 1;                 % Jl block size = pg.NS
        ul_hi =  (pg.NS) + gridS.a.n;
        rhs(ul_lo:ul_hi) = rhs(ul_lo:ul_hi) - param.zeta * pg.a.erg(:);
    end

    % normalization row
    Msys     = [Aall.'; ones(1,nVar)];
    rhs(end) = 1;

    fNew = Msys \ rhs;

    % unpack (ordering inside AL: [Jl ; Ul ; Jc ; Uc], then AO: [Ol ; Oc])
    cur = 0;
    solutionS.el = reshape(fNew(cur+(1:pg.NS)), size(solutionS.JlValue)); cur = cur + pg.NS;
    solutionS.hl =             fNew(cur+(1:gridS.a.n));                   cur = cur + gridS.a.n;
    solutionS.ec = reshape(fNew(cur+(1:pg.NS)), size(solutionS.JcValue)); cur = cur + pg.NS;
    solutionS.hc =             fNew(cur+(1:gridS.a.n));                   cur = cur + gridS.a.n;
    solutionS.ol =             fNew(cur+(1:gridS.a.n));                   cur = cur + gridS.a.n;
    solutionS.oc =             fNew(cur+(1:gridS.a.n));

    % ================== nested generators ==================
    function [AJ_J, AJ_U, AJ_O] = generateAJl(S) % employed, no kid
        AJ_J = pg.CJ;                                 % (NS×NS) exogenous component
        AJ_U = pg.CJ_U;                               % (NS×a)
        AJ_O = spalloc(pg.NS, gridS.a.n, 0);         % (NS×a) endogenous quits

        for ai = 1:gridS.a.n
            for ri = 1:gridS.w.n
                i = (ri-1)*gridS.a.n + ai;

                % quit J -> O if O strictly better
                if S.OlValue(ai) > S.JlValue(ai,ri) + tol
                    AJ_J(i,i)  = AJ_J(i,i)  - pg.lambda_o;
                    AJ_O(i,ai) = AJ_O(i,ai) + pg.lambda_o;
                elseif tremble > 0
                    t = min(pg.lambda_o, tremble * exp(S.OlValue(ai) - S.JlValue(ai,ri)));
                    AJ_J(i,i)  = AJ_J(i,i)  - t;
                    AJ_O(i,ai) = AJ_O(i,ai) + t;
                end

                % accept better offers J -> J'
                for r2 = 1:gridS.w.n
                    if S.JlValue(ai,ri) + tol < S.JlValue(ai,r2)
                        w  = pg.lambda * param.f_long(r2);
                        i2 = (r2-1)*gridS.a.n + ai;
                        AJ_J(i,i)  = AJ_J(i,i)  - w;
                        AJ_J(i,i2) = AJ_J(i,i2) + w;
                    end
                end
            end
        end
        AJ_J = sparse(AJ_J);
    end

    function [AU_J, AU_U, AU_O] = generateAUl(S) % unemployed, no kid
        AU_J = spalloc(gridS.a.n, pg.NS,     0);
        AU_U = sparse(gridS.a.n,  gridS.a.n) + pg.CU;
        AU_O = spalloc(gridS.a.n, gridS.a.n, 0);

        for ai = 1:gridS.a.n
            % U -> O if O strictly better
            if S.OlValue(ai) > S.UlValue(ai) + tol
                AU_U(ai,ai) = AU_U(ai,ai) - pg.lambda_o;
                AU_O(ai,ai) = AU_O(ai,ai) + pg.lambda_o;
            elseif tremble > 0
                t = min(pg.lambda_o, tremble * exp(S.OlValue(ai) - S.UlValue(ai)));
                AU_U(ai,ai) = AU_U(ai,ai) - t;
                AU_O(ai,ai) = AU_O(ai,ai) + t;
            end

            % U -> J(ai,r) via offers
            for r2 = 1:gridS.w.n
                i2 = (r2-1)*gridS.a.n + ai;
                if S.UlValue(ai) + tol < S.JlValue(ai,r2)
                    w = pg.lambda_u * param.f_long(r2);
                    AU_U(ai,ai) = AU_U(ai,ai) - w;
                    AU_J(ai,i2) = AU_J(ai,i2) + w;
                elseif tremble > 0
                    t = min(pg.lambda_u*param.f_long(r2), ...
                            tremble*exp(S.UlValue(ai)-S.JlValue(ai,r2))*param.f_long(r2));
                    AU_U(ai,ai) = AU_U(ai,ai) - t;
                    AU_J(ai,i2) = AU_J(ai,i2) + t;
                end
            end
        end
    end

    function [AO_J, AO_O] = generateAOl(S) % OLF, no kid
        % Direct OLF -> Employed job draws at rate lambda_oE (paper: nu_g),
        % accepted if the job beats staying OLF. No OLF -> Unemployed path:
        % the paper does not model OLF agents "starting to search" as a
        % separate step, only a direct job-offer draw from f(w).
        AO_J = spalloc(gridS.a.n, pg.NS, 0);
        AO_O = sparse(gridS.a.n,  gridS.a.n) + pg.CU; % stay-in-O components
        for ai = 1:gridS.a.n
            for wi_new = 1:gridS.w.n
                i_new = (wi_new-1)*gridS.a.n + ai;
                if S.OlValue(ai) + tol < S.JlValue(ai, wi_new)
                    w = pg.lambda_oE * param.f_long(wi_new);
                    AO_O(ai,ai)    = AO_O(ai,ai)    - w;
                    AO_J(ai,i_new) = AO_J(ai,i_new) + w;
                elseif tremble > 0
                    t = min(pg.lambda_oE, tremble * exp(S.OlValue(ai) - S.JlValue(ai,wi_new))) * param.f_long(wi_new);
                    AO_O(ai,ai)    = AO_O(ai,ai)    - t;
                    AO_J(ai,i_new) = AO_J(ai,i_new) + t;
                end
            end
        end
    end

    function [AJ_J, AJ_U, AJ_O] = generateAJc(S) % employed, kid
        AJ_J = pg.CJ; AJ_U = pg.CJ_U; AJ_O = spalloc(pg.NS,gridS.a.n,0);
        for ai = 1:gridS.a.n
            for ri = 1:gridS.w.n
                i = (ri-1)*gridS.a.n + ai;
                if S.OcValue(ai) > S.JcValue(ai,ri) + tol
                    AJ_J(i,i)  = AJ_J(i,i)  - pg.lambda_o;
                    AJ_O(i,ai) = AJ_O(i,ai) + pg.lambda_o;
                elseif tremble > 0
                    t = min(pg.lambda_o, tremble * exp(S.OcValue(ai) - S.JcValue(ai,ri)));
                    AJ_J(i,i)  = AJ_J(i,i)  - t;
                    AJ_O(i,ai) = AJ_O(i,ai) + t;
                end
                for r2 = 1:gridS.w.n
                    if S.JcValue(ai,ri) + tol < S.JcValue(ai,r2)
                        w  = pg.lambda * param.f_long(r2);
                        i2 = (r2-1)*gridS.a.n + ai;
                        AJ_J(i,i)  = AJ_J(i,i)  - w;
                        AJ_J(i,i2) = AJ_J(i,i2) + w;
                    end
                end
            end
        end
        AJ_J = sparse(AJ_J);
    end

    function [AU_J, AU_U, AU_O] = generateAUc(S) % unemployed, kid
        AU_J = spalloc(gridS.a.n, pg.NS,     0);
        AU_U = sparse(gridS.a.n,  gridS.a.n) + pg.CU;
        AU_O = spalloc(gridS.a.n, gridS.a.n, 0);
        for ai = 1:gridS.a.n
            if S.OcValue(ai) > S.UcValue(ai) + tol
                AU_U(ai,ai) = AU_U(ai,ai) - pg.lambda_o;
                AU_O(ai,ai) = AU_O(ai,ai) + pg.lambda_o;
            elseif tremble > 0
                t = min(pg.lambda_o, tremble * exp(S.OcValue(ai) - S.UcValue(ai)));
                AU_U(ai,ai) = AU_U(ai,ai) - t;
                AU_O(ai,ai) = AU_O(ai,ai) + t;
            end
            for r2 = 1:gridS.w.n
                i2 = (r2-1)*gridS.a.n + ai;
                if S.UcValue(ai) + tol < S.JcValue(ai,r2)
                    w = pg.lambda_u * param.f_long(r2);
                    AU_U(ai,ai) = AU_U(ai,ai) - w;
                    AU_J(ai,i2) = AU_J(ai,i2) + w;
                elseif tremble > 0
                    t = min(pg.lambda_u*param.f_long(r2), ...
                            tremble*exp(S.UcValue(ai)-S.JcValue(ai,r2))*param.f_long(r2));
                    AU_U(ai,ai) = AU_U(ai,ai) - t;
                    AU_J(ai,i2) = AU_J(ai,i2) + t;
                end
            end
        end
    end

    function [AO_J, AO_O] = generateAOc(S) % OLF, kid
        % Same direct OLF -> Employed mechanism as generateAOl, for parents.
        AO_J = spalloc(gridS.a.n, pg.NS, 0);
        AO_O = sparse(gridS.a.n,  gridS.a.n) + pg.CU;
        for ai = 1:gridS.a.n
            for wi_new = 1:gridS.w.n
                i_new = (wi_new-1)*gridS.a.n + ai;
                if S.OcValue(ai) + tol < S.JcValue(ai, wi_new)
                    w = pg.lambda_oE * param.f_long(wi_new);
                    AO_O(ai,ai)    = AO_O(ai,ai)    - w;
                    AO_J(ai,i_new) = AO_J(ai,i_new) + w;
                elseif tremble > 0
                    t = min(pg.lambda_oE, tremble * exp(S.OcValue(ai) - S.JcValue(ai,wi_new))) * param.f_long(wi_new);
                    AO_O(ai,ai)    = AO_O(ai,ai)    - t;
                    AO_J(ai,i_new) = AO_J(ai,i_new) + t;
                end
            end
        end
    end

    % ---------- instantaneous utilities (explicit broadcasting) ----------
    function u = instUtilityOl
        cons = param.no_tax(pg.bout(gridS.a.a), gridS.a.a, param)+param.os;
        u    = util(param, cons) .* ones(gridS.a.n,1);
    end
    function u = instUtilityUl
        cons = param.nokid_tax(param.b(gridS.a.a), gridS.a.a, param)+param.bs;
        u    = util(param, cons) .* ones(gridS.a.n,1);
    end
    function u = instUtilityJl
        cons = param.nokid_tax(W_eff, AA, param)+param.t;
        u    = util(param, cons);
    end
    function u = instUtilityOc
        cons = param.no_tax(pg.boutc(gridS.a.a), gridS.a.a, param) - c_a +param.os+param.oc;
        u    = util(param, cons) .* ones(gridS.a.n,1);
    end
    function u = instUtilityUc
        cons = param.onekid_tax(param.b(gridS.a.a), gridS.a.a, param) - c_a +param.bs+param.bc;
        u    = util(param, cons) .* ones(gridS.a.n,1);
    end
    function u = instUtilityJc
        cons = param.onekid_tax(W_eff, AA, param).*param.fkidwpenalty - repmat(c_a,1,gridS.w.n) +param.t;
        u    = util(param, cons);
    end

    % ---------- single value update step ----------
    function [S, dmax] = updateValue(S)
        % LFP (no kid)
        [S.AJl_Jl, S.AJl_Ul, S.AJl_Ol] = generateAJl(S);
        [S.AUl_Jl, S.AUl_Ul, S.AUl_Ol] = generateAUl(S);
        % LFP (kid)
        [S.AJc_Jc, S.AJc_Uc, S.AJc_Oc] = generateAJc(S);
        [S.AUc_Jc, S.AUc_Uc, S.AUc_Oc] = generateAUc(S);
        % OLF
        [S.AOl_Jl, S.AOl_Ol]           = generateAOl(S);
        [S.AOc_Jc, S.AOc_Oc]           = generateAOc(S);

        % subtract pcs on “no kid” diagonals only; add coupling to kid
        ALl   = [S.AJl_Jl, S.AJl_Ul;
                 S.AUl_Jl, S.AUl_Ul] - pg.pcs .* speye(pg.NSall);

        ALc   = [S.AJc_Jc, S.AJc_Uc;
                 S.AUc_Jc, S.AUc_Uc];

        ALl_c = pg.pcs .* speye(pg.NSall);
        ALc_l = spalloc(pg.NSall, pg.NSall, 0);
        S.AL  = [ALl,  ALl_c;
                 ALc_l, ALc];

        % OLF: subtract pcs on Ol diag only (Ol -> Oc coupling below)
        AO_Ol = S.AOl_Ol - pg.pcs .* speye(gridS.a.n);
        AO_Oc = S.AOc_Oc;
        AOl_c = pg.pcs .* speye(gridS.a.n);
        AOc_l = spalloc(gridS.a.n, gridS.a.n, 0);
        S.AO  = [AO_Ol, AOl_c;
                 AOc_l, AO_Oc];

        % Cross L <-> O
        % rows follow [Jl; Ul; Jc; Uc], columns [Ol, Oc]
        leftOl  = [ S.AJl_Ol; ...
                    S.AUl_Ol; ...
                    spalloc(pg.NS,      gridS.a.n, 0); ...
                    spalloc(gridS.a.n,  gridS.a.n, 0) ];
        rightOc = [ spalloc(pg.NS+gridS.a.n, gridS.a.n, 0); ...
                    S.AJc_Oc; ...
                    S.AUc_Oc ];
        S.AL_O  = [sparse(leftOl), sparse(rightOc)];   % (2*NSall) × (2*a)

        % rows [Ol; Oc], columns [Jl; Ul; Jc; Uc]
        % OLF re-entry lands directly in the employed (Jl/Jc) block, not Ul/Uc.
        topOl = [ S.AOl_Jl, ...
                  spalloc(gridS.a.n, gridS.a.n, 0), ...
                  spalloc(gridS.a.n, pg.NS, 0), ...
                  spalloc(gridS.a.n, gridS.a.n, 0) ];
        botOc = [ spalloc(gridS.a.n, pg.NS, 0), ...
                  spalloc(gridS.a.n, gridS.a.n, 0), ...
                  S.AOc_Jc, ...
                  spalloc(gridS.a.n, gridS.a.n, 0) ];
        S.AO_L = [sparse(topOl); sparse(botOc)];       % (2*a) × (2*NSall)

        % full operator and Bellman linear system
        ALO = [S.AL, S.AL_O;
               S.AO_L, S.AO];

        uJl = instUtilityJl;  uUl = instUtilityUl;
        uJc = instUtilityJc;  uUc = instUtilityUc;
        uOl = instUtilityOl;  uOc = instUtilityOc;

        uL  = [reshape(uJl, pg.NS,1); uUl; reshape(uJc, pg.NS,1); uUc];
        uO  = [uOl; uOc];
        u   = [uL; uO];

        B = (param.rho) .* speye(size(ALO)) - ALO;
        V = B \ u;

        % unpack values
        cur = 0;
        JlNew = reshape(V(cur+(1:pg.NS)), size(S.JlValue));  cur = cur + pg.NS;
        UlNew = V(cur+(1:gridS.a.n));                       cur = cur + gridS.a.n;
        JcNew = reshape(V(cur+(1:pg.NS)), size(S.JcValue));  cur = cur + pg.NS;
        UcNew = V(cur+(1:gridS.a.n));                        cur = cur + gridS.a.n;
        OlNew = V(cur+(1:gridS.a.n));                        cur = cur + gridS.a.n;
        OcNew = V(cur+(1:gridS.a.n));                        cur = cur + gridS.a.n;

        dmax = max([ max(abs(JlNew-S.JlValue),[],'all') , ...
                     max(abs(UlNew-S.UlValue))           , ...
                     max(abs(JcNew-S.JcValue),[],'all') , ...
                     max(abs(UcNew-S.UcValue))           , ...
                     max(abs(OlNew-S.OlValue))           , ...
                     max(abs(OcNew-S.OcValue)) ]);

        upd = param.updMult;
        S.JlValue = (1-upd).*S.JlValue + upd.*JlNew;
        S.UlValue = (1-upd).*S.UlValue + upd.*UlNew;
        S.JcValue = (1-upd).*S.JcValue + upd.*JcNew;
        S.UcValue = (1-upd).*S.UcValue + upd.*UcNew;
        S.OlValue = (1-upd).*S.OlValue + upd.*OlNew;
        S.OcValue = (1-upd).*S.OcValue + upd.*OcNew;

        % expose assembled blocks for the outer scope
        S.AL   = sparse(S.AL);
        S.AO   = sparse(S.AO);
        S.AL_O = sparse(S.AL_O);
        S.AO_L = sparse(S.AO_L);
    end
end
