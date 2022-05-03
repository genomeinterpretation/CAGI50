function [LLRPlus, RR, scores] = PRSFromClinical(classes,scores, prior)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin < 3
    prior = mean(classes);
end

[scores,ix] = randomSort(scores);
classes = classes(ix);
opts.alpha =  prior;
opts.EvidenceNotRequired = true;
opts.presorted = true;
out = clinicalMetrics(classes,scores, opts);
LLRPlus = out.llrp_raw';
RR = LLRPlus./(prior*LLRPlus + (1-prior));

end

