function [col_sup, col_mod, col_st] = colorsEvidence()
addpath('utilities/cbrewer/');
colors_unique = cbrewer('qual', 'Set1', 3);
col_sup = colors_unique(3,:);
col_mod = colors_unique(2,:);
col_st = colors_unique(1,:);
end