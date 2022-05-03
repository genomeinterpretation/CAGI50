function ix = isPathogenic(sig_str)
%returns a boolean indicating id the variant belongs to pathogenic class
ix_path = contains(sig_str, 'pathogenic', 'IgnoreCase', true);
ix_dm  = contains(sig_str, 'DM' , 'IgnoreCase', true);
ix = ix_path | ix_dm;
end

