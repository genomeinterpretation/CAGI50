function  perfs = selectMethodsAndBaselines(perfs_all, ds, k)
%% Sorts methods in perfs_all based on the ranking measures applicable to the dataset. 
%% Select the best performer (representative) from each group and return it in 'perfs' 
%% in the sorted order. Any method with coverage below 0.9 is listed at the end. 
%% The best performing baseline and Experimental-Max, if present in 'perfs_all' are 
%% listed after all the methods.
% Inputs
% perfs_all: a structure array, each entry of which contains the
%   measures for a method computed using 'perfMetrics()'.
% ds: The Dataset object for the dataset.
% k: (numeric) Gives the number of main methods (excluding baselines and
%   Experimental-Max) to be included in the output variable 'perfs'. If the 
%   number of main methods is less than k, then all main methods are included.
%   If perfs_all also contains meta methods, a total 2k (k meta and k non meta)
%   main methods can be included, provided they are present in perfs_all. 
% Output
% perfs: a structure array containg subset of methods from perfs_all. 


rMeasures = rankingMeasures(ds);


if nargin < 3
    grps_uniq = unique({perfs_all.group});
    k = length(grps_uniq);
end


% Sepearate the methods with coverage above and below 0.9
% and sort them separately
ixCov = arrayfun(@(perf) perf.coverage > 0.9, perfs_all);
perfs1 = perfs_all(ixCov);
perfs_srt1 = sortPerfs(perfs1, rMeasures);
perfs2 = perfs_all(~ixCov);
perfs_srt2 = sortPerfs(perfs2, rMeasures);
% Put the mehods with less than 0.9 coverage at the end.
perfs_srt = [perfs_srt1, perfs_srt2];

% Get the index of all the meta methods, baselines and Experimental-Max.
ix_eMax = find([perfs_srt.isExpMax]);
ix_bl = find([perfs_srt.isBaseline]);
ix_meta = find([perfs_srt.isMeta]);

% Get the index of non meta main methods.
ix_noMeta = setdiff(1:length(perfs_srt), [ix_meta, ix_bl, ix_eMax]);

% Extract the top k meta and non meta methods. Only one method is included
% from each group. And arrange the extracted main methods in the sorted
% order.
[~, ix_noMeta2] = firstByGrp(perfs_srt(ix_noMeta), k);
[~, ix_meta2] = firstByGrp(perfs_srt(ix_meta), k);
srt_ix_main = sort([ix_meta(ix_meta2), ix_noMeta(ix_noMeta2)]);
perfs_main = perfs_srt(srt_ix_main);

% Extract the top performing baseline and Experimental-Max
perfs_eMax = perfs_srt(ix_eMax);
perfs_bl = perfs_srt(ix_bl);
if ~isempty(perfs_bl)
    perfs_bl = perfs_bl(1);
end
if ~isempty(perfs_eMax)
    perfs_eMax = perfs_eMax(1);
end

% Put the selected baseline and Experiemntal-Max after the selected main
% methods.
perfs = [perfs_main, perfs_bl, perfs_eMax];



end