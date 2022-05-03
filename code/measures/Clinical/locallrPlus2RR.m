function RR = locallrPlus2RR(llrPlus,prior)
%% Converts local lr+ to Relative Risk w.r.t. the entire reference population.
RR = llrPlus./((llrPlus-1)*prior + 1);
end

