function ix = isBenign(sig_str)
%returns a boolean indicating id the variant belongs to benign class
ix = contains(sig_str, 'benign', 'IgnoreCase', true);
end

