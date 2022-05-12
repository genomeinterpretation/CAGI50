function posterior = locallrPlus2Posterior(llrp, prior)
%% Converts local lr+ to posterior adjusted to prior.
ratio = llrp*(prior)/(1-prior);
posterior = ratio./(1 + ratio);
posterior(isnan(posterior)) = 1;
end

