function [llrPlus, RR, scores] = PRSFromClinical(classes, scores, prior)
% Computes the llrPlus and RR used in the PRS figure for Crohns dataset
if nargin < 3
    prior = mean(classes);
end
[scores,ix] = randomSort(scores);
classes = classes(ix);
opts.EvidenceNotRequired = true;
opts.presorted = true;
out = clinicalMetrics(classes,scores, prior, opts);
llrPlus = out.llrPlus';
RR = llrPlus./(prior*llrPlus + (1-prior));
end

