function grp2Color = group2Color(perfs)
    [BLUE, GREEN, RED, PURPLE, ORANGE, GREY, DBLUE] = mainColors();
    groups = {perfs.group};
    nGroups = length(groups);
    ix_baseline = find([perfs.isBaseline], 1);
    ix_expMax = find([perfs.isExpMax], 1);
    if ~isempty(ix_baseline)
        grp2Color.(groups{ix_baseline}) = ORANGE;
    end
    if ~isempty(ix_expMax)
        grp2Color.(groups{ix_expMax}) = GREY;
    end
    fourColors = [BLUE;GREEN;RED;PURPLE];
    ix_main = setdiff(1:nGroups, [ix_baseline, ix_expMax]);
    if length(ix_main) > 4
        error(['No more than four main methods, in addition to the baseline and ' ...
            'Experimental-Max, can be assigned colors with this code.' ...
            'Please modify the code to have more method colors'])
    end
    for m = 1:length(ix_main)
        grp2Color.(groups{m}) = fourColors(m,:);
    end
end