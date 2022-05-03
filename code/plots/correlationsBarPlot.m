function correlationsBarPlot(perfs, ds, showPearson, showKendall, showSpearman)
%% Creates the correlation bar plots for all methods in 'perfs'.
% Inputs:
% perfs: a structure array, each entry of which contains the
%   measures for a method computed using 'perfMetrics()'.
% ds: The Dataset object for the dataset.
% showPearson: (boolean) Tells if pearson correlation should be included in
%   the figure.
% showKendall: (boolean) Tells if kendall's tau should be included in
%   the figure.
% showSpearman: (boolean) Tells if psearman's correlation should be included in
%   the figure.


    corrStrs = {};
    corrLabels = {};
    if showPearson
        corrStrs = [corrStrs, {'pearson'}];
        corrLabels = [corrLabels, {"Pearson's corr."}];
    end
    if showKendall
        corrStrs = [corrStrs, {'kendall'}];
        corrLabels = [corrLabels, {"Kendall's tau"}];
    end
    if showSpearman
        corrStrs = [corrStrs, {'spearman'}];
        corrLabels = [corrLabels, {"Spearman's corr."}];
    end
    [BLUE, GREEN, RED, PURPLE, ORANGE, GREY, DBLUE] = mainColors();
    
    nPerfs = length(perfs);
    grp2Color = group2Color(perfs);
    methodNames = {};
    for c = 1:length(corrStrs)
        corrStr = corrStrs{c};
        if c > 1
            xline((nPerfs+1)*(c-1), 'LineStyle','-', 'LineWidth',2, 'Color', [0.5,0.5,0.5])
        end
        for p = 1:length(perfs)
            perf = perfs(p);
            corr = perf.(corrStr).val;
            cor95 = perf.(corrStr).bt.p95;
            cor5 = perf.(corrStr).bt.p5;
            pos = (c-1)*(nPerfs+1) + p;
            b(pos) = bar(pos, corr);
            hold on;
            if c == 1
                legendEntries(pos) = b(pos);
                name = methodDisplayName(perf.group, perf.method);
                methodNames = [methodNames, {name}];
            end
            b(pos).FaceColor = 'flat';
            col = grp2Color.(perf.group);
            b(pos).CData(1,:) = col;
            text(pos, 0.07, numToStr(corr),'FontSize', 12, 'HorizontalAlignment','center');
            plot([pos,pos], [cor5, cor95], '-k', 'LineWidth', 2);
            plot([pos-0.2,pos+0.2], [cor5, cor5], '-k', 'LineWidth', 2);
            plot([pos-0.2,pos+0.2], [cor95, cor95], '-k', 'LineWidth', 2);
        end
    end
    
    ax = gca;
    
    xticklabels(repmat({''}, length(ax.XTick)))
    ylim([0,1])
    ax.Box = 'on';
    ax.BoxStyle = 'full';
    ax.LineWidth = 1.5;
    set(gca,'XTick',[])
    hold off;
    set(gca,'XTick',(1:2:(2*length(corrStrs)))*(nPerfs+1)/2)
    xticklabels(corrLabels);
    ylabel('Correlation', 'Rotation', 90)
    set(gca, 'FontSize', 14)
    title(ds.displayName)
    lgd = legend(legendEntries, methodNames, 'Location', 'southoutside', 'Orientation','horizontal', NumColumns=2);
    lgd.FontSize = 14;
    legend boxoff            
end