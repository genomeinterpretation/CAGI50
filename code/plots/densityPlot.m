function densityPlot(perf, ds)
    [BLUE, GREEN, RED, PURPLE, ORANGE, GREY, DBLUE] = mainColors();
    scores  = perf.scores;
    classes = perf.classes;
    scores0 = scores(classes==0);
    scores1 = scores(classes==1);
    [scores, ix] = sort(scores);
    classes = classes(ix);
    [f0, x0] = ksdensity(scores0 , scores);
    p(1) = plot(x0, f0, 'Color', GREEN, 'DisplayName', 'negatives', 'LineWidth' ,3);
    hold on
    [f1, x1] = ksdensity(scores1, scores);
    p(2) = plot(x1, f1, 'Color', RED, 'DisplayName', 'positives', 'LineWidth', 3);
    legend(p,{'negatives', 'positives'});
    hold off;
    ylabel('Density', 'FontSize', 20);
    xlab = [ds.classScoreLabel, ' by ', methodDisplayName(perf.group, perf.method)];
    xlabel(xlab)
    ax=gca;
    ax.FontSize = 16;
    lgd = legend('Location', 'northeast');
    lgd.FontSize = 20;
    legend boxoff;  
end