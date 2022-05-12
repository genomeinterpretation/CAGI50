function [lrPlus_sup, lrPlus_mod, lrPlus_st, lrPlus_vst, C, fails] = prior2lrThresholds(prior, original, strict)
%% Computes the lr+ thresholds for the supporting (sup), moderate (mod), strong (st)
%% and very strong (vst) levels of evidence from the prior under Tavtigan's framework [1].
% Note that [1] derives thresholds for "odds of pathogenicity (OP)", which from its 
% relationship to the posterior probability of pathogenicity (given in the
% paper), can be inferred to be the positive likelihood ratio. Note that the 
% positive likelihood ratio is equal to the ratio of the posterior odds of 
% pathogencity to the prior odds of pathogenciity. 
% In our work, we interpret the lr+ thresholds as the thresholds for local
% lr+, instead of the global LR+. See the CAGI flagship paper for details.
% [1] Tavtigian, Sean V., et al. "Modeling the ACMG/AMP variant classification 
% guidelines as a Bayesian classification framework." Genetics in Medicine 
% 20.9 (2018): 1054-1060.
%
% Inputs:
% prior: (numeric scalar in (0,1)) the proprtion of positives in the
%   reference population.
% original: (boolean, defualt false) set to true for using the original 
%   ACMG/AMP combing rules. Set to false if modified rules are to be used.
% strict: (boolean, defualt false) set to true if a likely pathogenic or
%   a likely benign rule should fail when the posterior probability goes
%   above 0.99 or below 0.01, respectively. Set to false if a rule should
%   not fail in those cases.
%
% Output:
% lrPlus_xxx: (numeric scalar) the lr+ threshold at xxx evidence level.
% C: (numeric scalar) The constant O_PVSt in [1] 
% fails: (numeric scalar) The number of rules failed when using C.
if nargin < 3
    strict = false;
end
if nargin < 2
    original = false;
end
% Search over 1:30000 for C with the lowest number of rules failed.
Cs = 1:30000;
% Number of pathogenic rules failed at each candidate C because the posterior probability was 
% not above 0.99. 
paths = arrayfun(@(x) sum(failsIfLessThan0(pathogenicRulesPosterior(x, prior, original)-0.99)), Cs);
% Number of likely pathogenic rules failed at each candidate C because the posterior 
% probability was not above 0.9. 
lpaths = arrayfun(@(x) sum(failsIfLessThan0(likelypathogenicRulesPosterior(x, prior, original)-0.90)), Cs);
% Number of likely pathogenic rules failed at each candidate C because the posterior 
% probability was not below 0.99.
lpathsStrict = arrayfun(@(x) sum(failsIfLessThan0(0.99-likelypathogenicRulesPosterior(x, prior, original))), Cs);
% Number of benign rules failed at each candidate C because the posterior 
% probability was not below 0.01.
benigns = arrayfun(@(x) sum(failsIfLessThan0(0.01 - benignRulesPosterior(x, prior))), Cs);
% Number of likely benign rules failed at each candidate C because the posterior 
% probability was not below 0.1.
lbenigns = arrayfun(@(x) sum(failsIfLessThan0(0.1 - likelybenignRulesPosterior(x, prior, original))), Cs);
% Number of likely benign rules failed at each candidate C because the posterior 
% probability was not above 0.01.
lbenignsStrict = arrayfun(@(x) sum(failsIfLessThan0(likelybenignRulesPosterior(x, prior, original)-0.01)), Cs);

% Add the counts of failiures from diferent rule categories for each
% candidate C.
if strict 
    fails = lpaths + paths + benigns + lbenigns+ lpathsStrict +lbenignsStrict; 
else
    % the likely pathogenic and likely benign rules achieving posterior probability 
    % above 0.99 and below 0.01, respectively, are not counted as
    % failiures.
    fails = lpaths + paths + benigns + lbenigns;
end
% Find the C among the candidates with the least number of failiures. When the
% minimum is achieved at multiple candidates, pick the smallest one.
[~,ix] = min(fails);
fails = fails(ix);
C = Cs(ix);
% Compute the lrPlus thresholds from C.
lrPlus_sup = C^0.125;
lrPlus_mod = C^0.25;
lrPlus_st= C^0.5;
lrPlus_vst= C; 

    function failiure = failsIfLessThan0(delta)
        % Assigns a failiure count of 1 to all indices of delta where delta<0,
        % otherwise assign 0.
        ixx = delta < 0;
        failiure = zeros(length(delta),1);
        failiure(ixx) = 1;
    end
end

