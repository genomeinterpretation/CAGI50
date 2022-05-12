function [perfs_srt, index_srt, rank_srt, val_srt] = sortPerfs(perfs, rankingMeasures)
%% Sorts the methods in 'perfs' based on the 'rankingMeasures' 
    nMetrics = length(rankingMeasures);
    nPerfs = length(perfs);
    values = nan(nMetrics, nPerfs);
    ranks = nan(nMetrics, nPerfs);
    avgVal = zeros(1, nPerfs);
    for  j = 1:length(rankingMeasures)
        values(j,:) = arrayfun(@(perf) perf.specs2property(rankingMeasures(j), 'val'), perfs);
        if rankingMeasures(j).ismax
            ranks(j,:) = tiedrank(-values(j,:));
            avgVal = avgVal + values(j,:);
        else
            ranks(j,:) = tiedrank(values(j,:));
            avgVal = avgVal - values(j,:);
        end
    end
    % avg_val is used to break ties if avgRank are tied
    avgVal = avgVal/nMetrics;
    avgRank = mean(ranks,1);
    groups = arrayfun(@(perf) upper(perf.group), perfs, 'UniformOutput', false);
    index = 1:nPerfs;
    
    index_srt = sortIndexByRank(index);
    perfs_srt = perfs(index_srt);
    rank_srt = avgRank(index_srt);
    val_srt = avgVal(index_srt);
    
    function [index_srt, rank] = sortIndexByRank(index)
        % if there are ties in rank use the avg_val to break rank ties
        unique_rank = unique(avgRank(index));
        index2 = [];
        for rank = unique_rank
            ixx = find(avgRank(index) == rank);
            [index_srt, ~] = sortIndexByVal(index(ixx));
            index2 = [index2, index_srt];
        end
        if length(index2) ~= length(index)
            ixx = find(isnan(avgRank(index)));
            index2 = [index2, index(ixx)]; 
            if length(index2) ~= length(index)
                error('index2 and index should have same lengths')
            end
        end
        [~, ixx] = sort(avgRank(index2));
        index_srt = index2(ixx);
        rank = avgRank(index_srt);
     end
    
    function [index_srt, val] = sortIndexByVal(index)
        [~, ixx] = sort(avgVal(index), 'descend');
        index_srt = index(ixx);
        val = avgVal(index_srt);
    end

end

