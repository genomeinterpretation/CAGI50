function scatterPlot(perf, ds, prior) 
%% Applicable when both regression and classification analysis are performed and the 
%% classification scores are created by using the experimental value predictions as is or 
%% after negating. The figure contains the scatter plot of experimental values versus 
%% its prediction for the method evaluated in 'perf'. The points belonging to the positive
%% and negative classes are colored purple and yellow, respectively. Any remaining points are
%% colored grey. The clinical threshold along with their 95% confidence intervald are also 
%% displayed if a clinical analysis is applicable. If the ground truth class labels were 
%%   created by applying a class boundary to the Experiemental values, the boundary 
%%   is be displayed on the figure as a horizontal  and vertical line separating the 
%%   experimental values and the predictions respectively.

% Inputs:
% perf: contains the measures evaluated for a method computed using 'perfMetrics()'.
% ds: The Dataset object for the dataset.
% prior: (numeric scalar) contains the class prior for which the clinical thresholds are
%   displayed. The prior should be from the ones used while 'perf' was computed. 
%   If prior is empty the clinical results are not displayed.


    
   expValLabel = ds.expValLabel;
    
   classBoundary = ds.classBoundary;
    
    if length(classBoundary) > 2
        error('the class boundary cannot contain more than 2 entries')
    end
    if nargin < 3 || isempty(prior)
        showClinThr = false;
    else 
        showClinThr = true;
    end
    
    yLabel = ['observed ', expValLabel];
    xLabel = ['predicted ', expValLabel, ' by ' methodDisplayName(perf.group, perf.method)];
    figTitle = ds.displayName;
    
    
    % Extract the experimental values and the predictions
    preds = perf.predictions;
    expVals = perf.expVals;
    n = length(expVals);
    
    % Extract index of postives, negatives and neutral experimental values.
    class_index = perf.class_index;
    pos_index = class_index(perf.classes==1);
    neg_index = class_index(perf.classes==0);
    neutral_index = setdiff(1:n, class_index);
    % If variant positions are well defined, then the figure should only show
    % a single variant for each position. To this end, a random index is picked
    % from all indices having the same position
    if ~isempty(perf.positions) && ~any(isnan(perf.positions))
        positions = perf.positions;
        pos_unique = unique(positions);
        index = arrayfun(@(pos) randsample(find(positions == pos), 1), pos_unique);
        pos_index = pos_index(ismember(pos_index, index));
        neg_index = neg_index(ismember(neg_index, index));
        neutral_index = neutral_index(ismember(neutral_index, index));
    end
    
    % Pick colors for positive, negative and neutral points
    [BLUE, GREEN, RED, PURPLE, ORANGE, GREY, DBLUE] = mainColors();
    col_pos = PURPLE;
    col_neg = ORANGE;
    col_neut = [0.3,0.3,0.3];
    
    % Plot the experimental values agains the predicitons, colored to reflect
    % the classes.
    scatter(preds(pos_index), expVals(pos_index), 50, col_pos, 'LineWidth', 2);
    hold on
    scatter(preds(neg_index), expVals(neg_index), 50, col_neg, 'LineWidth', 2);
    scatter(preds(neutral_index), expVals(neutral_index), 50, col_neut, 'LineWidth', 2);
    

%     mn_x = min(preds);
%     mx_x = max(preds);
%     xx = [mn_x, mx_x];
%     plot(xx, xx, '-', 'LineWidth', 1, 'Color', [0.75,0.75,0.75]);
    
    % The limits of the y-axis are rounded to an appropriate decimal place 
    mn_y = min(expVals);
    mx_y = max(expVals);
    y95 = prctile(expVals, 95);
    y5 = prctile(expVals, 5);
    yRange = y95 -y5;
    y_lim = [myFloor(mn_y, yRange), myCeil(mx_y, yRange)];
    ylim(y_lim);
    
    
    % Plot the clinical thresholds and its 95% percent condifdence interval.
    if showClinThr
        % Pick colors for supporting, moderate and strong evidence
        [col_sup, col_mod, col_st] = colorsEvidence();
    
        clinical = perf.clinical;  
    
        % If the classification scores were obtained by negating the
        % predicitons, the scores and the clinical thresholds should be negated
        % back to have the original sign. 
        scores = clinical.scores.val;
        invertScores = perf.negated;
        sign = (1-2*invertScores);
        scores = sign*scores;
        
        % Find if the clinical results were computed w.r.t. the prior. 
        % If yes, extract the prior index.
        priors = clinical.priors.val;
        prior_ix = find(abs(prior - priors)< 10^-3, 1);
       
        
        % Determine the width of the panel beneath te scatter plot that will
        % contain the clinical thresholds.
        fudge = (y_lim(2) - y_lim(1))*0.05;
        width1 = (y_lim(2) - y_lim(1))*0.3;
        width = width1/3;
        ylim([y_lim(1)-(width1 + fudge), y_lim(2)])
        ylimFinal = ylim;
        
    
        % Display the prior.    
        x_lim = xlim;
        letterWidth = (x_lim(2)- x_lim(1))*0.02;
        if invertScores
            text(x_lim(2) - letterWidth, y_lim(1)- 1.75*width, ['prior =', '\rm ', num2str(prior,2)], 'FontSize', 14, 'HorizontalAlignment', 'right')
            ev_pos_left = false;
        else
            text(x_lim(1) + letterWidth, y_lim(1)- 1.75*width, ['prior =', '\rm ', num2str(prior,2)], 'FontSize', 14)
            ev_pos_left = true;
        end
        

        % Plot the clinical thresholds and their 95% confidence intervals
        if clinical.tpr_sup.val(prior_ix) > 0.0
            thr_sup = sign*clinical.thr_sup.val(prior_ix);
            thr_5 = sign*clinical.thr_sup.bt.p5(prior_ix);
            thr_95 = sign*clinical.thr_sup.bt.p95(prior_ix);
            plotThresholds(thr_sup, thr_5, thr_95, col_sup, 1, width, fudge, 'Sup.');
        end
        if clinical.tpr_mod.val(prior_ix) > 0.0
            thr_mod = sign*clinical.thr_mod.val(prior_ix);
            thr_5 = sign*clinical.thr_mod.bt.p5(prior_ix);
            thr_95 = sign*clinical.thr_mod.bt.p95(prior_ix);
            plotThresholds(thr_mod, thr_5, thr_95, col_mod, 2, width, fudge, 'Mod.');
        end
        if clinical.tpr_st.val(prior_ix) > 0.0
            thr_st = sign*clinical.thr_st.val(prior_ix);
            thr_5 = sign*clinical.thr_st.bt.p5(prior_ix);
            thr_95 = sign*clinical.thr_st.bt.p95(prior_ix);
            plotThresholds(thr_st, thr_5, thr_95, col_st, 3, width, fudge, 'Str.');
        end
    end

    % Draw the line separating the scatter plot from clinical threshold
    % panel
    yline(y_lim(1), 'LineWidth',2)
    
    
    % Show the class boundaries defined on the experimental values for
    % both the experimental values and the predictions.
    if ~isempty(classBoundary)
        % classBoundary can contain one or two numbers. Two class boundaries are required 
        % when neutral variants are also defined.
        for cb = classBoundary
            x_lim = xlim;
            xRange = x_lim(2) - x_lim(1);
            % Plot the class boundary for experimental values
            yline(cb, 'LineStyle','--', 'Color', [0.3,0.3,0.3], 'LineWidth',2);
            % Plot the class boundary for predictions
            plot([cb, cb], [y_lim(1), ylimFinal(2)], 'LineStyle','--', 'Color', [0.5,0.5,0.5], 'LineWidth',2);
            % Extend the x-axis limits if the class boundaries are not in
            % the range of predictions.
            if x_lim(1) > cb- xRange*0.05
                xlim([cb - xRange*0.05, x_lim(2)])
            elseif cb + xRange*0.05 > x_lim(2)
                xlim([x_lim(1), cb + xRange*0.05])
            end
        end
    end

    ax = gca;
    noTicklabels_ix = find(ax.YTick < y_lim(1));
    ax.YTickLabel(noTicklabels_ix) = repmat({''}, 1, length(noTicklabels_ix));

    % Plot a line at 45 degree anlge.
    xx = [y_lim(1), min(ylimFinal(2), x_lim(2))];
    plot(xx, xx, '-', 'LineWidth', 1, 'Color', [0.75,0.75,0.75]);

    
    % Display the title, x-axis and y-axis labels
    title(figTitle, 'FontSize', 14);
    xlabel(xLabel,'FontSize', 14);
    ylabel(yLabel, 'FontSize', 14);
    
    ax.Box = 'on';
    ax.BoxStyle = 'full';
    ax.LineWidth = 1.5;
    set(gca, 'FontSize', 14)


    function plotThresholds(thr, thr_5, thr_95, col, pos, width, fudge, evidence)
        % Plot light vertical lines at the clinical thresholds spanning the
        % the scatter plot and the clinical threshold pannel
        plot([thr,thr], [ylimFinal(1), ylimFinal(2)], '-', 'LineWidth', 1, 'Color', [0.8,0.8,0.8]);
        
        % Plot the clinical threshold and its 95% confidence intervals
        y_pos = y_lim(1) - ((pos-1) + 0.5)*width - fudge;
        scatter(thr, y_pos, 50, col, 'filled', 'MarkerEdgeColor', 'k');
        plot([thr_5,thr_95], [y_pos, y_pos], 'Color', col, 'LineWidth', 2);
        plot([thr_5,thr_5], [y_pos-width/2, y_pos+width/2], 'Color', col, 'LineWidth', 3);
        plot([thr_95,thr_95], [y_pos-width/2, y_pos+width/2], 'Color', col,'LineWidth', 3);
       
        % Display the evidence descriptors (Sup., Mod., Str.).
        x_lim = xlim;
        if invertScores
            % If there is more empty space to the left of the
            % confidence interval than the right, then plot the 
            % evidnce descriptors to the left
            if pos == 1 && x_lim(2) - thr_5 < thr_95 - x_lim(1)
                ev_pos_left = true;
            end
            if ~ev_pos_left
                text(thr_5+ letterWidth, y_pos, evidence, 'FontSize', 12)
            else
                text(thr_95 - letterWidth, y_pos, evidence, 'FontSize', 12, 'HorizontalAlignment','right')
            end
        else
            % If there is more empty space to the right of the
            % confidence interval than the left, then plot the 
            % evidnce descriptors to the right
            if pos == 1 && x_lim(2) - thr_95 > thr_5 - x_lim(1)
                ev_pos_left = false;
            end
            if ~ev_pos_left
                text(thr_5 - letterWidth, y_pos, evidence, 'FontSize', 12, 'HorizontalAlignment','right')
            else
                text(thr_95 + letterWidth, y_pos, evidence, 'FontSize', 12)
            end
        end
    end

end

