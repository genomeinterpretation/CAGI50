function posterior = pathogenicRulesPosterior(C, prior, original)
%% Computes the posterior values for the ACMG/AMP pathogenic rules under 
%% Tavtigan's framework [1] for a candidate C (O_PVSt) and a prior.
% [1] Tavtigian, Sean V., et al. "Modeling the ACMG/AMP variant classification 
% guidelines as a Bayesian classification framework." Genetics in Medicine 
% 20.9 (2018): 1054-1060.
% Inputs:
% C: (numeric scalar) The constant O_PVSt in [1]. 
% prior: (numeric scalar in (0,1)) the proprtion of positives in the
%   reference population.
% original: (boolean, defualt false) set to true for using the original 
%   ACMG/AMP combing rules. Set to false if modified rules are to be used.
%
% Output:
% posterior: (numeric vector of length 8 (equal to number of pathogenic rules))
%   the posterior computed for each rule.

% fracs contains the exponents of C to be used to compute the positive likelihood 
% ratio (lr+) for a single line of supporting, moderate, strong and very strong 
% evidence towards pathogenicity. For example, the the lr+ for a single line of 
% supporting evidence is C^(2^-3). 
fracs = [2^-3, 2^-2, 2^-1, 1]';

% The combined lr+ of a rule is obtained by multiplying the lr+ for evidences 
% supported in the rule. For example, for the fifth rule, requiring 1 strong 
% and 3 moderate lines of evidence, combined lr+ is C^(3*2^-2 + 2^-1) =
% (C^(3*2^-2)) * (C^(1*2^-1)).
% The posterior for the rule is computed from the prior and the combined lr+.
% The vector of counts below specify the number of evidence of each kind
% supported in the rule. For example, the fifth rule supports 3 moderate
% and 1 strong evidence.
posterior(1) = locallrPlus2Posterior(C^(dotprod([0,0,1,1],fracs)), prior);
posterior(2) = locallrPlus2Posterior(C^(dotprod([0,2,0,1],fracs)), prior);
posterior(3) = locallrPlus2Posterior(C^(dotprod([1,1,0,1],fracs)), prior);
posterior(4) = locallrPlus2Posterior(C^(dotprod([2,0,0,1],fracs)), prior);
posterior(5) = locallrPlus2Posterior(C^(dotprod([0,3,1,0],fracs)), prior);
posterior(6) = locallrPlus2Posterior(C^(dotprod([2,2,1,0],fracs)), prior);
posterior(7) = locallrPlus2Posterior(C^(dotprod([4,1,1,0],fracs)), prior);
%posterior(8) = llrp2Posterior(C^(dotprod([0,1,0,1],fracs)), prior);
if original
    % The original rules consider 2 strong lines of evidence as pathogenic. 
    posterior(8) = locallrPlus2Posterior(C^(dotprod([0,0,2,0],fracs)), prior);
else
    % The modified rules consider 1 moderate and 1 very strong lines as
    % pathogenic. This replaces the 2 strong lines of evidence rule, which
    % is moved to the likely pathogenic rules.
    posterior(8) = locallrPlus2Posterior(C^(dotprod([0,1,0,1],fracs)), prior);
end
end

