function posteriorAndllrPlusPlot(perf, ds, priors, show_llrp, smoothing)
%% Creates the posterior and/or local lr+ plots for the method in 'perf' as function of 
%% the classification score. The clinical thresholds along with their 95% confidence intervals 
%% are also diplayed. The percent of positive predictions (ppp) at the thresholds is also displayed.
% Inputs:
% perf: contains the measures evaluated for a method computed using 'perfMetrics()'.
% ds: The Dataset object for the dataset.
% priors: (numeric array) contains all the class priors for which the
%   posteriors should be displayed. If empty, the posteriors, clinical thresholds and ppp are
%   are not displayed. The prior should be from the ones used while 'perf' was computed. 
%   Note that the displayed ppp values also depend on the prior. 
% show_llrp: (boolean, default: true) Tells if the local lr+ should be included in
%   the figure. In figures with more than one prior (and posterior), it is
%   best to not include posterior and lr+ isn the same figure for clarity.
% smoothing: (boolean, defualt: true) Tells if the smoothed versions of the posterior and 
%   local lr+ curve should be displayed. 
% clsScoreLabel: The x-axis label to be used in the figure. If not provided
%   it is computed by 'ds2clsScoreLabel()' function


yLabel_llrp = 'local lr+';
yLabelPosterior = 'posterior probability of pathogenicity';

if perf.negated
    priorLegendLoc = 'east';
    evidenceLegendLoc = 'northeast';
else
    priorLegendLoc = 'west';
    evidenceLegendLoc = 'northwest';
end
    



xLabel = [ds.classScoreLabel, ' by ' methodDisplayName(perf.group, perf.method)];

if nargin < 5
    smoothing = true;
end

if nargin < 4
    show_llrp =true;
end
showPosterior = true;
if nargin < 3 || isempty(priors)
    showPosterior = false;
else
    nPriors = length(priors);
    if nPriors > 2
        warning('the figure will be too busy with more than two posterior curves')
    end
    if nPriors > 3
        warning('Cannot plot more than 3 posterior curves')
    end
end


[BLUE, GREEN, RED, PURPLE, ORANGE, GREY, DBLUE] = mainColors();
[col_sup, col_mod, col_st] = colorsEvidence();

postColor = RED;

llrpColor = BLUE;

grey = [0.8,0.8,0.8];

invertScores = perf.negated;
sign = (1-2*invertScores);
clinical = perf.clinical;
scores = clinical.scores.val;
scores = sign*scores;
thr_sup_all = sign*clinical.thr_sup.val;
thr_mod_all = sign*clinical.thr_mod.val;
thr_st_all = sign*clinical.thr_st.val;
post_sup_all = clinical.post_sup.val;
post_mod_all = clinical.post_mod.val;
post_st_all = clinical.post_st.val;
llrp_sup_all = clinical.llrp_sup.val;
llrp_mod_all = clinical.llrp_mod.val;
llrp_st_all = clinical.llrp_st.val;

lx = min(scores);
ux = max(scores);

fig = gcf;
ax = gca;
ax.Color = 'none';


if smoothing
    rawCurveLS = ':';
else
    rawCurveLS = '-';
end

if showPosterior
    prior_str={};
    post_curve = [];
    priors = sort(priors);
    posteriors =  clinical.posterior_adj.val;
    max_posterior = max(arrayfun(@(pr) max(posteriors{find(abs(pr-perf.priors)<10^-3, 1)}), priors));
    for i = 1:nPriors
        postColor_i = i*postColor/nPriors;
        prior = priors(i);
        %prior_str = [prior_str, ['prior = ', num2str(prior, 2)]];
        prior_str = [prior_str, [num2str(prior, 2)]];
        prior_ix = find(abs(prior-perf.priors)<10^-3, 1);
        posterior = posteriors{prior_ix};
        
        thr_sup = thr_sup_all(prior_ix);
        thr_mod = thr_mod_all(prior_ix);
        thr_st = thr_st_all(prior_ix);
        post_sup = post_sup_all(prior_ix);
        post_mod = post_mod_all(prior_ix);
        post_st = post_st_all(prior_ix);
        llrp_sup = llrp_sup_all(prior_ix);
        llrp_mod = llrp_mod_all(prior_ix);
        llrp_st = llrp_st_all(prior_ix);
      
        %ax1=axes('position',get(ax,'position'),'visible','off');
        ax1 = ax;
        hold on
        %uistack(ax1,"bottom");
        xlim(ax1, [lx,ux])
        if ~isnan(thr_sup)
            if show_llrp
                yyaxis right
                plot([ux,thr_sup], repmat(llrp_sup,1,2), 'LineStyle', '-', 'LineWidth', 0.25, 'Color',grey);
                hold on
                yyaxis left
            end
             plot([lx,thr_sup], repmat(post_sup,1,2), 'LineWidth', 0.25, 'LineStyle', '-', 'Color',grey);
             plot([thr_sup, thr_sup], [0,1], 'LineStyle', '-', 'LineWidth', 0.25, 'Color',grey)
        end
        if ~isnan(thr_mod)
            if show_llrp
                yyaxis right
                plot([ux,thr_mod], repmat(llrp_mod,1,2), 'LineStyle', '-', 'LineWidth', 0.25, 'Color',grey);
                yyaxis left
            end
            plot([lx,thr_mod],repmat(post_mod,1,2), 'LineWidth', 0.25, 'LineStyle', '-', 'Color', grey);
            plot([thr_mod, thr_mod], [0,1], 'LineStyle', '-', 'LineWidth', 0.25, 'Color',grey)
        end
        if ~isnan(thr_st)
            if show_llrp
                yyaxis right
                plot([ux,thr_st], repmat(llrp_st,1,2), 'LineStyle', '-', 'LineWidth', 0.25, 'Color',grey);
                yyaxis left
            end
            plot([lx,thr_st],repmat(post_st,1,2), 'LineWidth', 0.25, 'LineStyle', '-', 'Color', grey);
            plot([thr_st, thr_st], [0,1], 'LineStyle', '-', 'LineWidth', 0.25, 'Color',grey)
        end
        set(ax1, 'FontSize', 14)
        fig.CurrentAxes = ax;
    
        
        %---------------------Plot posterior----------------------
        
        post_curve(i) = plot(scores, posterior, rawCurveLS, 'Color', postColor_i, 'LineWidth', 3);
        hold on
        if smoothing
            post_smooth = smoothCurve(scores, posterior, 2);
            plot(scores, post_smooth, '-', 'Color', postColor_i, 'LineWidth', 3, 'Marker', 'none');
        end
        ylabel(yLabelPosterior, 'FontSize', 14)
        ylim([0, myCeil(max_posterior, max_posterior)])
     
        
        %--------------------Put markers at thresholds---------------
        s = [];
        if ~isnan(thr_sup)
            s(1) = scatter(thr_sup, post_sup, 50, 'MarkerFaceColor', col_sup, 'LineWidth', 1, 'MarkerEdgeColor','k', 'DisplayName','Sup.' );
        end
        if ~isnan(thr_mod)
            s(2) = scatter(thr_mod, post_mod, 50, 'MarkerFaceColor', col_mod, 'LineWidth', 1, 'MarkerEdgeColor','k', 'DisplayName','Mod.');
        end
        if ~isnan(thr_st)
            s(3) = scatter(thr_st, post_st, 50, 'MarkerFaceColor', col_st, 'LineWidth', 1, 'MarkerEdgeColor','k', 'DisplayName','Str.');
        end
        
        
        % -----------------PPP--------------------------------------
        offset_y = max_posterior*0.05;
        %offset_x = -sign*(ux-lx)*0.04;
        offset_x = 0;
        if sign == -1
            align = 'left';
        else
            align = 'right';
            offset_x = 0;
        end
       
        if ~isnan(thr_sup)
            ppp_sup = [num2str(clinical.ppp_sup.val(prior_ix)*100, 2),'%'];
            text(thr_sup + offset_x, post_sup + offset_y, ppp_sup, 'FontSize', 14, 'HorizontalAlignment', align)
        end
        if ~isnan(thr_mod)
            ppp_mod = [num2str(clinical.ppp_mod.val(prior_ix)*100, 2),'%'];
            text(thr_mod + offset_x, post_mod+offset_y,ppp_mod, 'FontSize', 14, 'HorizontalAlignment',align)
        end
        if ~isnan(thr_st)
           ppp_st = [num2str(clinical.ppp_st.val(prior_ix)*100, 2),'%'];
           text(thr_st + offset_x, post_st+offset_y,ppp_st, 'FontSize', 14, 'HorizontalAlignment', align)
        end
    
    
    
    
   
    if ~isnan(thr_sup)
        hold on
        thr_sup5 = sign*clinical.thr_sup.bt.p5(prior_ix);
        thr_sup95 = sign*clinical.thr_sup.bt.p95(prior_ix);
        plotCI(thr_sup5, thr_sup95,  post_sup, col_sup);
    end
    if ~isnan(thr_mod)
        thr_mod5 = sign*clinical.thr_mod.bt.p5(prior_ix);
        thr_mod95 = sign*clinical.thr_mod.bt.p95(prior_ix);
        plotCI(thr_mod5, thr_mod95, post_mod, col_mod);
    end
    if ~isnan(thr_st)
        thr_st5 = sign*clinical.thr_st.bt.p5(prior_ix);
        thr_st95 = sign*clinical.thr_st.bt.p95(prior_ix);
        plotCI(thr_st5, thr_st95, post_st, col_st);
    end
   
    end
end

    %---------------------Plot lr+----------------------
   
    if show_llrp
        if showPosterior
            fig.CurrentAxes = ax1;
            yyaxis right
        end
        xlim([lx,ux])
        hold on
        llrPlus = clinical.llrPlus.val;
        plot(scores, llrPlus, rawCurveLS, 'Color', llrpColor, 'LineWidth', 3);
        if smoothing
            if ~showPosterior
                smooth_llrp = smoothCurve(scores, llrPlus, 2);
            else
                smooth_llrp = (post_smooth./(1-post_smooth))*(1-prior)/prior;
            end
            plot(scores,  smooth_llrp, '-','Color', llrpColor, 'LineWidth', 3);
        end
        llrp_lim =  max(llrPlus);
        if showPosterior
            llrp_lim = llrp_lim*1.5;
        end

        ylim2 = [0, llrp_lim + 1];
        ylim(ylim2)
        ylabel(yLabel_llrp, 'FontSize', 14)
        if showPosterior
            yyaxis left
            fig.CurrentAxes = ax;
            yyaxis right
            ylim(ylim2);
            yyaxis left
        end
    end

xlim([lx,ux])
ax.YColor = [0 0 0];
ax.Color = 'none';


x_limits = xlim;
if show_llrp && showPosterior
    ax.YAxis(1).Color = [0 0 0];
    ax.YAxis(2).Color = [0 0 0];
    xline(x_limits(1), 'Color', postColor, 'LineWidth', 3)
    xline(x_limits(2), 'Color', llrpColor, 'LineWidth', 3)
end

if showPosterior
    lgd = legend(post_curve, prior_str);
    legend('Location', priorLegendLoc);
    title(lgd, 'Prior');
    title(ds.displayName, 'FontSize', 14)
    lgd.FontSize = 14;
    legend('boxoff');
    ax = gca;
    if ~isempty(s)
        ax2=axes('position',get(ax,'position'),'visible','off');
        lgd2 = legend(ax2, s,'Location',evidenceLegendLoc, 'Orientation','vertical');
        lgd2.FontSize = 14;
        lgd2.ItemTokenSize(1) = 10;
        legend('Box','off');
        set(ax2, 'FontSize', 14)
    end
end





fig.CurrentAxes = ax;
xlabel(xLabel, 'FontSize', 14)
ax.Box = 'on';
ax.BoxStyle = 'full';
ax.LineWidth = 1.5;
set(ax, 'FontSize', 14)



    function plotCI(thr_5, thr_95, y_pos, col)
        y_lim = ylim;
        width = (y_lim(2)-y_lim(1))*0.1;
        plot([thr_5,thr_95], [y_pos, y_pos], 'LineWidth', 2, 'Color',col, 'LineStyle', '-',  'Marker', 'none');
        plot([thr_5,thr_5], [y_pos-width/2, y_pos+width/2], 'Color', col, 'LineWidth', 2, 'LineStyle', '-',  'Marker', 'none');
        plot([thr_95,thr_95], [y_pos-width/2, y_pos+width/2], 'Color', col,'LineWidth', 2, 'LineStyle', '-', 'Marker', 'none');
    end

end


