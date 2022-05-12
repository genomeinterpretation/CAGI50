function ROCPlot(perfs, ds, type, prior, changeCurveWidth)
%% Displays the ROC plots for all methods in 'perfs' and the random classifier.
%% The function can be used to display the 'standard', 'truncated' and the 'log-log' 
%% ROC plots. The area under the curve is also given. The clinical thresholds corresponding 
%% to the class prior in 'prior' is also displayed for the truncated ROC curve.
% Inputs
% perfs: a structure array, each entry of which contains the
%   measures for a method computed using 'perfMetrics()'.
% ds: The Dataset object for the dataset.
% type: (string: 'standard', 'truncated', 'log-log') Specified the type of
%   ROC curve to plot.
% prior: (numeric scalar) contains the class prior for which the clinical thresholds are
%   displayed for the truncated ROC curve. The prior should be from the ones used while 'perf' was computed. 
%   If prior is empty the clinical results are not displayed.
% changeCurveWidth: (boolean, default: false) When the ROC curves of the methods have
%   significant overlap, set 'changeCurveWidth' to true so that ROC curves
%   are plotted in decreasing width order and can be seen in spite of the
%   overlap.
    showClinical = true;
    if isempty(prior) || nargin <4
        showClinical = false;
    end

    if nargin <5
        changeCurveWidth = false;
    end
    fprLabel = 'false positive rate';
    tprLabel = 'true positive rate';
    isStandard = false;
    isTrunc = false;
    isLog = false;
    if strcmpi(type, 'standard')
        isStandard = true;
        aucStr = 'roc auc';
        fprStr = 'roc fpr';
        tprStr = 'roc tpr';
    
    end
    if strcmpi(type, 'truncated')
        isTrunc = true;
        aucStr = 'roc auc_trunc';
        fprStr = 'roc fpr_trunc';
        tprStr = 'roc tpr_trunc';
    end
    if strcmpi(type, 'log-log')
        isLog = true;
        aucStr = 'roc auc_log';
        fprStr = 'roc fpr_log';
        tprStr = 'roc tpr_log';
        tprLabel = ['log_1_0 ', tprLabel];
        fprLabel = ['log_1_0 ', fprLabel];
    end
    
    aucSpec = struct('specs', aucStr);
    fprSpec = struct('specs', fprStr);
    tprSpec = struct('specs', tprStr);
    nPerfs = length(perfs);
    grp2Color = group2Color(perfs);
    [col_sup, col_mod, col_st] = colorsEvidence();
    lgdLabels = {};
    yRange = 1;
    if isStandard
        auc_random = 0.5;
        plts(nPerfs + 1) = plot([0,1], [0,1], 'LineStyle','--','Color',[0.5,0.5,0.5]);
    elseif isTrunc
        auc_random = 0.1;
        plts(nPerfs + 1) = plot([0, 0.2], [0, 0.2], 'LineStyle','--','Color',[0.5,0.5,0.5]);
    elseif isLog
        n1 = sum(perfs(1).allClasses);
        n0 = sum(1-perfs(1).allClasses);
        tpr_log_min = log10(1/n1);
        fpr_log_min = log10(1/n0) - log10(2);
        auc_random = randomAUC_log(n1,n0);
        plts(nPerfs + 1) = plot([fpr_log_min, 0] ,[tpr_log_min,0], 'LineStyle','--','Color',[0.5,0.5,0.5]);
        xlim([fpr_log_min, 0]);
        ylim([tpr_log_min, 0]);
        yRange = -tpr_log_min;
    end
    hold on
    label = ['random', ' (', numToStr(auc_random),')'];
    lgdLabels = [{label}, lgdLabels];

    for i = length(perfs):-1:1
        perf = perfs(i);
        methodName = methodDisplayName(perf.group, perf.method);
        auc = perf.specs2property(aucSpec, 'val');
        bt = perf.specs2property(aucSpec, 'bt');
        std = bt.std;
        fpr = perf.specs2property(fprSpec, 'val');
        tpr = perf.specs2property(tprSpec, 'val');
        label = [methodName, ' (', numToStr(auc),'±',numToStr(1.96*std),')'];
        lgdLabels = [{label}, lgdLabels];
        col = grp2Color.(perf.group);
        if changeCurveWidth
            plts(i) = plot(fpr, tpr, 'Color', col , 'LineWidth', i+2);
        else
            plts(i) = plot(fpr, tpr, 'Color', col , 'LineWidth', 3);
        end
    end
    
    if showClinical
        prior_ix = find(abs(perfs(1).priors-prior)<10^-3, 1);
        y = perfs(1).clinical.tpr_sup.val(prior_ix); x=perfs(1).clinical.fpr_sup.val(prior_ix);
        if y > 0 
            if isLog
               x = max(log10(x), fpr_log_min); y = max(log10(y), tpr_log_min);
            end
            scatter(x, y,70, col_sup, 'filled', 'o', 'MarkerEdgeColor', 'k')
            text(x, y+yRange*0.04, 'Sup.', 'HorizontalAlignment','center','FontSize',12)
        end
        y = perfs(1).clinical.tpr_mod.val(prior_ix); x=perfs(1).clinical.fpr_mod.val(prior_ix);
        if y > 0 
            if isLog
                x = max(log10(x), fpr_log_min); y = max(log10(y), tpr_log_min);
            end
            scatter(x, y,70, col_mod, 'filled', 'o', 'MarkerEdgeColor', 'k')
            text(x, y+ yRange*0.04, 'Mod.', 'HorizontalAlignment','center','FontSize',12)
        end
        y = perfs(1).clinical.tpr_st.val(prior_ix); x =perfs(1).clinical.fpr_st.val(prior_ix);
        if y > 0
            if isLog
                x = max(log10(x), fpr_log_min); y = max(log10(y), tpr_log_min);
            end
            scatter(x, y,70, col_st, 'filled', 'o', 'MarkerEdgeColor', 'k');
            text(x, y+ yRange*0.04, 'Str.', 'HorizontalAlignment','center','FontSize',12);
        end
    end
    
    xlabel(fprLabel)
    ylabel(tprLabel, 'Rotation', 90)
    title(ds.displayName)
    
    lgd = legend(plts, lgdLabels, 'Location', 'best','Orientation','horizontal','NumColumns',2);
    lgd.FontSize = 14;
    legend boxoff
    ax = gca;
    ax.Box = 'on';
    ax.BoxStyle = 'full';
    ax.LineWidth = 1.5;
    set(gca, 'FontSize', 14)
    
end