function digit_strs = extractFrontDigits(strs)
digit_strs = cellfun(@(celll) celll{1}, regexp(strs, '\d*', 'match'), 'UniformOutput', false);
end

