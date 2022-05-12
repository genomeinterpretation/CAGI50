function v = averageStrings(strstr)
v = nan(length(strstr),1);
for i = 1: length(strstr)
    str = strstr{i};
    if ~isempty(str)
        strstr2 = regexpi(str,',|;', 'split');
        strstr2 = cellfun(@(str) extractFirstDecimal(str), strstr2, 'UniformOutput', false);
        if iscell(strstr2)
            ix =  cellfun(@(str) ~isempty(str2num(str)), strstr2);
            v(i) =  mean(cellfun(@(str) str2num(str), strstr2(ix)));
        else
            %if the str is made up of a single number.
            %applying mean would convert a non-number to a NAN.
            v(i) = mean(str2num(strstr2));
        end
    end
end
    function strd = extractFirstDecimal(str)
        str = regexpi(str, '-*\d*\.?\d*', 'match');
        strd = 'nan';
        if ~isempty(str)
            strd = str{1};
        end
    end
end

