function v = likelybenignRulesPosterior(e, alpha, original)
%% Computes the posterior values for the ACMG/AMP likely benign rules under 
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
% posterior: (numeric vector of length 1 (number of benign rules))
%   The posterior computed for each rule.

% fracs contains the exponents of C to be used to compute the positive likelihood 
% ratio (lr+) for a single line of supporting, moderate, strong and very strong 
% evidence towards benignity. For example, the the lr+ for a single line of 
% supporting evidence is C^(-2^-3). 
fracs = -[2^-3, 2^-2, 2^-1, 1]';

% The combined lr+ of a rule is obtained by multiplying the lr+ for evidences 
% supported in the rule. For example, for the rule requiring 1 supporting and 
% 1 strong evidence, combined lr+ is C^(-2^-3 -2^-1).
% The posterior for the rule is computed from the prior and the combined lr+.
% The vector of counts below specify the number of evidence of each kind
% supported in the rule. For example, the rule supports
% and 1 supporting and 1 strong evidence.
v(1) = locallrPlus2Posterior(e^(dotprod([1,0,1,0],fracs)), alpha);
if original
     % The original rules consider 2 supporting lines of evidence as likely benign.
     % This rule is not included in the modified rules.
    v(2) = locallrPlus2Posterior(e^(dotprod([2,0,0,0],fracs)), alpha);
end
end