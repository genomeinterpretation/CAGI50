function ix = isClinvar(sig_str)
% Returns the index
ix = ~(strcmpi(sig_str, 'DM?') | strcmpi(sig_str, 'DM'));
end

