function ix = isClinvarPathogenic(sig_str)
% Returns the index
ix = contains(sig_str, 'pathogenic', 'IgnoreCase', true);
end
