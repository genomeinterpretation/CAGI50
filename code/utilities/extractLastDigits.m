function digit_strs = extractLastDigits(strs)
digit_strs = cellfun(@(celll) celll{end}, regexp(strs, '\d*', 'match'), 'UniformOutput', false);
end
