function RRPlot(perf, ds, prior, smoothing)
    [BLUE, GREEN, RED, PURPLE, ORANGE, GREY, DBLUE] = mainColors();
    prior_ix = find(abs(prior - 0.013) < 10^-3, 1);
    rr = perf.clinical.RR.val{prior_ix};
    scores = perf.clinical.scores.val;
    plot(scores, log2(rr), 'LineStyle', ':', 'LineWidth', 2, 'Color', BLUE)
    hold on;
    if smoothing
        smoothrr = smoothCurve(scores,rr,2);
        plot(scores, log2(smoothrr), 'LineStyle', '-', 'LineWidth', 2, 'Color', BLUE)
    end
    ax = gca;
    ax.YTick = [-2, -1, 0, 1,2,3];
    ax.YTickLabel =  arrayfun(@(n) num2str(n), 2.^ax.YTick, 'UniformOutput', false);
    ylabel('Relative Risk', 'FontSize', 16)
    xlab = [ds.classScoreLabel, ' by ', methodDisplayName(perf.group, perf.method)];
    xlabel(xlab, 'FontSize', 16)
end