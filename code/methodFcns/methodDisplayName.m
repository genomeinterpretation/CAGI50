function dispName = methodDisplayName(group, method)
%% returns the display name of a method
Grps = {'ExpMax', 'Mora', 'Licht', 'Moult', 'Dunbr', 'Jones', 'Bromb', 'Bolog', 'PolyPhen'};
DispNames = {'Experimental-Max', 'MutPred2', 'Evolutionary Action', 'Maryland', 'FoxChase', 'London', 'Bromberg', 'Bologna', 'PolyPhen-2'};

ix  = find(strcmpi(group, Grps));
if ~isempty(ix)
    dispName = DispNames{ix};
else
    dispName = group;
end
end