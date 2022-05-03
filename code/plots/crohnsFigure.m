function crohnsFigure(Perfs, saveFig)
    addpath('subtightplot/subtightplot/');
    addpath('metrics/');
    BLUE = [0, 0.45, 0.74];
    GREEN = [0.47, 0.67, 0.19];
    RED = [0.85, 0.33, 0.10];
    PURPLE = [0.49, 0.18, 0.56];
    ORANGE = [0.93, 0.69, 0.13];
    GREY = [0.5, 0.5, 0.5];
    %ref_prc = [0.1,0.25,0.5,0.75,0.9,1];
    %ColorsRR = [BLUE; GREEN; RED; PURPLE; ORANGE; GREY];
    ref_prc = 1;
    ColorsRR = BLUE;
    pix =1;
    %Bins = {'0%-25%','25%-50%', '50%-75%', '75%-100%'};
    %Bins = {'0%-10%','10%-20%','20%-30%','30%-40%', '40%-50%', '50%-60%', '60%-70%','70%-80%', '80%-90%', '90%-100%'};
    nboot = 1000;
    if ~saveFig
        subplot = @(i) subtightplot (2,2, i, [0.15 0.1], [0.1 0.05], [0.1 0.01]);
        subplot(pix);
        pix = pix+1;
    else
        figure;
    end
    ds = fieldnames(Perfs);
    perfs = Perfs.(ds{1}).perfs;
    metric = struct('specs1', 'roc auc', 'ismax', true);
    
    [perfs_srt,index_srt, rank_srt, val_srt] = sortPerfs(perfs, metric);
    topTwo = firstByGrp(perfs_srt, 2);
    top = topTwo(1);
    second = topTwo(2);
    classes = top.classes;
    scores = top.scores;
    prior = top.clinical.alpha.val;
    %[fpr,tpr,t]=perfcurve(classes, scores, 1);
    fpr = top.roc.fpr.val;
    tpr = top.roc.tpr.val;
    auc = top.roc.auc.val;
    auc = num2str(auc,2);
    std = 1.96*top.roc.auc.bt.std;
    std = num2str(std,2);
    label = [top.group, '(' ,num2str(auc,2),'±',std(2:end),')'];
    plot(fpr,tpr, 'color' , BLUE, 'LineWidth',2,  'DisplayName', label);
    hold on
    classes2 =second.classes;
    scores2 = second.scores;
    fpr = second.roc.fpr.val;
    tpr = second.roc.tpr.val;
    auc = second.roc.auc.val;
    auc = num2str(auc,2);
    std =  1.96*second.roc.auc.bt.std;
    std = num2str(std,2);
    label = [second.group, '(' ,num2str(auc,2),'±',std(2:end),')'];
    plot(fpr,tpr, 'color' , GREEN, 'LineWidth',2,  'DisplayName', label);
    %title('ROC','FontSize', 20);
    xlabel('fpr');
    ylabel('tpr');
%     a = get(gca,'XTickLabel');  
%     set(gca,'XTickLabel',a,'fontsize',14)
    ax=gca;
    ax.FontSize = 16;
    lgd = legend('Location', 'southeast');
    lgd.FontSize = 20;
    legend boxoff  
    hold off
    
    
    if ~saveFig
        subplot(pix);
        pix = pix +1;
    else
        saveas(gcf, fullfile('..', 'plotsTables','Crohns_ROC' ))
        figure;
    end
    
    
    scores0 = scores(classes==0);
    scores1 = scores(classes==1);
    [scores, ix] = sort(scores);
    classes = classes(ix);
    [f0, x0] = ksdensity(scores0 , scores);
%     par0 = betafit(scores0);
%     f0 = betapdf(scores, par0(1), par0(2));
    p(1) = plot(x0, f0, 'Color', GREEN, 'DisplayName', 'negatives', 'LineWidth' ,3);
    hold on
    %histogram(scores0, 15, 'Normalization', 'pdf', 'FaceAlpha', 0.5)
    [f1, x1] = ksdensity(scores1, scores);
%     par1 = betafit(scores1);
%     f1 = betapdf(scores, par1(1), par1(2));
    p(2) = plot(x1, f1, 'Color', RED, 'DisplayName', 'positives', 'LineWidth', 3);
    %histogram(scores1, 15, 'Normalization', 'pdf','FaceAlpha', 0.5)
    legend(p,{'negatives', 'positives'});
    hold off;
    %title('Density', 'FontSize', 20);
    ylabel('Density', 'FontSize', 20);
    xlabel('method score')
    %xlim([0,1.01]);
    %a = get(gca,'XTickLabel');  
    %set(gca,'XTickLabel',a,'fontsize',14);
    ax=gca;
    ax.FontSize = 16;
    lgd = legend('Location', 'northeast');
    lgd.FontSize = 20;
   legend boxoff;  
   
    
    if ~saveFig
        subplot(pix);
        pix = pix +1;
    else
         saveas(gcf, fullfile('..', 'plotsTables','Crohns_densities' ))
         figure;
    end
    [LLRPlus_bt, RR_bt, scores_bt] = PRSFromClinical(classes, scores, prior);
%     for i = 1:nboot
%         pos_bt = bootstrp(1, @(x) x, scores1);
%         neg_bt = bootstrp(1, @(x) x, scores0);
%         cls = [ones(length(pos_bt),1); zeros(length(neg_bt), 1)];
%         scs = [pos_bt'; neg_bt'];
%         [LRPlus(i,:), LRMinus(i,:), DOR(i,:), thr(i,:), tpr(i,:), fpr(i,:), tnr(i,:), fnr(i,:)]= PRS(cls, scs, 0.013, p);
%     end

    if nboot > 1
        for i = 1:nboot
            pos_bt = bootstrp(1, @(x) x, scores1);
            neg_bt = bootstrp(1, @(x) x, scores0);
            cls = [ones(length(pos_bt),1); zeros(length(neg_bt), 1)];
            scs = [pos_bt'; neg_bt'];
            [LLRPlus_bt(:,i), RR_bt(:,i), scores_bt(:,i)]= PRSFromClinical(cls, scs, prior);
        end
    end
     LLRPlus = nan(length(scores),1);
     RR = nan(length(scores),1);
    if nboot > 1
        for i = 1:length(scores)
            ixx = scores_bt == scores(i);
            LLRPlus(i) = mean(LLRPlus_bt(ixx));
            %LRMinus(i) = mean(LRMinus_bt(ixx)); 
            %RR(i) = mean(RR_bt(ixx));
        end
    else
         LLRPlus = LLRPlus_bt;
         %RR = RR_bt;
    end
    
    s = reshape(scores, 1, []);

%     plot(scores, PPV, 'color', BLUE, 'LineStyle', ':', 'LineWidth', 2);
%     hold on
%     ppv = reshape(PPV, 1, []);
%     smooth_ppv = posteriorSmoothNN(s,ppv,2);
%     plot(s,  smooth_ppv, 'color', BLUE, 'LineWidth', 2);
%     title('PPV','FontSize', 20)
%     ax=gca;
%     ax.FontSize = 16;
%     xlabel('method score');
%     if ~saveFig
%         subplot(pix);
%         pix = pix +1;
%     else
%          saveas(gcf, fullfile('..', 'plotsTables','Crohns_PPV' ))
%          figure;
%     end
   
    plot(scores, LLRPlus, 'color', BLUE, 'LineStyle', ':', 'LineWidth', 2);
    hold on
    LR = reshape(LLRPlus, 1, []);
    smoothLRP = posteriorSmoothNN(s,LR,2);
    plot(s,  smoothLRP, 'color', BLUE, 'LineWidth', 2);
    %title('LR^+','FontSize', 20)
    ylabel('local lr^+','FontSize', 20)
    ax=gca;
    ax.FontSize = 16;
    xlabel('method score');
    if ~saveFig
        subplot(pix);
        pix = pix +1;
    else
         saveas(gcf, fullfile('..', 'plotsTables','Crohns_LRPlus' ))
         figure;
    end
    
   
%     
%     boxplot(LRMinus)
%     title('LRMinus')
%     xticklabels(Bins)
%     if ~saveFig
%          subplot(4);
%     else
%          saveas(gcf, fullfile('..', 'plotsTables','Crohns_LRMinus' ))
%          figure;
%     end
    
%    plot(scores, DOR, 'color', BLUE, 'LineStyle', ':', 'LineWidth', 2);
%    hold on
%    DR = reshape(DOR, 1, []);
%    smoothDOR = posteriorSmoothNN(s,DR,2);
%    plot(s,  smoothDOR, 'color', BLUE, 'LineWidth', 2);
%    title('DOR','FontSize', 20);
%    xlabel('method score')
%    ax=gca;
%    ax.FontSize = 16 ;
%     if ~saveFig
%          subplot(5);
%     else
%          saveas(gcf, fullfile('..', 'plotsTables','Crohns_DOR' ))
%          figure;
%     end
   prc = {};
   RR = localLR2RR(LLRPlus,prior);
   prc = {[num2str(100),'%']};
   plt = plot(scores, log2(RR), 'color', ColorsRR(1,:), 'LineStyle', ':', 'LineWidth', 2);
   hold on
         rr = reshape(RR(:,1), 1, []);
         smoothrr = posteriorSmoothNN(s,rr,2);
%    smoothrr = localLR2RR(smoothLRP,prior);
    plot(s,  log2(smoothrr), 'color', ColorsRR(1,:), 'LineWidth', 2);
  
  
   xlabel('method score')
   ylabel('Relative Risk')
   ax=gca;
   ax.FontSize = 16 ;
   %lgd=legend('Lowest 10%','','100%','');
   %lgd=legend(prc{:});
   if length(ref_prc) > 1
        lgd=legend(plt, prc);
        lgd.Title.String= 'Reference population';
        legend('Location', 'northwest');
        legend('boxoff');
        %title('RR','FontSize', 20);
   else
       %title(['RR @ ',prc{1}],'FontSize', 20);
       %title('RR','FontSize', 20);
   end
   %When plotting log2 RR use
   ax.YTick = [-2, -1, 0, 1,2,3];
   ax.YTickLabel =  arrayfun(@(n) num2str(n), 2.^ax.YTick, 'UniformOutput', false);
    if ~saveFig
        subplot(pix);
        pix = pix +1;
    else
         saveas(gcf, fullfile('..', 'plotsTables','Crohns_RR' ))
         figure;
    end
    
  
%    p = find(ref_prc == 1, 1);
%    prc = [num2str(100),'%'];
%    plt = plot(scores, RR(:,p), 'color', ColorsRR(p,:), 'LineStyle', ':', 'LineWidth', 2);
%    hold on
%    rr = reshape(RR(:,p), 1, []);
%    smoothrr = posteriorSmoothNN(s,rr,2);
%    plot(s,  smoothrr, 'color', ColorsRR(p,:), 'LineWidth', 2);
%    title('RR','FontSize', 20);
%    xlabel('method score')
%    ax=gca;
%    ax.FontSize = 16 ;
%    %lgd=legend('Lowest 10%','','100%','');
%    %lgd=legend(prc{:});
%    lgd=legend(plt, prc);
%    
%    lgd.Title.String= 'Reference population';
%    legend('Location', 'northwest');
%     legend('boxoff');
%     if ~saveFig
%          subplot(6);
%     else
%          saveas(gcf, fullfile('..', 'plotsTables','Crohns_RR' ))
%          figure;
%     end

%  prc = {};
%    for p = 1:length(ref_prc)
%        prc = [prc, {[num2str(100*ref_prc(p)),'%']}];
%        plt(p) = plot(scores, log2(RR(:,p)), 'color', ColorsRR(p,:), 'LineStyle', ':', 'LineWidth', 2);
%         hold on
%         rr = reshape(RR(:,p), 1, []);
%         smoothrr = posteriorSmoothNN(s,rr,2);
%         plot(s,  log2(smoothrr), 'color', ColorsRR(p,:), 'LineWidth', 2);
%    end
%    title('log2(RR)','FontSize', 20);
%    xlabel('method score')
%    ax=gca;
%    ax.FontSize = 16 ;
%    %lgd=legend('Lowest 10%','','100%','');
%    %lgd=legend(prc{:});
%    lgd=legend(plt, prc);
%    
%    lgd.Title.String= 'Reference population';
%    legend('Location', 'northwest');
%     legend('boxoff');
%     if ~saveFig
%          subplot(7);
%     else
%          saveas(gcf, fullfile('..', 'plotsTables','Crohns_RR' ))
%          figure;
%     end


%    p2=plot(scores, RR2, 'color', PURPLE, 'LineStyle', ':', 'LineWidth', 2);
%    hold on
%    RR2 = reshape(RR2, 1, []);
%    smoothRR2 = posteriorSmoothNN(s,RR2,2);
%    plot(s,  smoothRR2, 'color', PURPLE, 'LineWidth', 2);
%    title('RR','FontSize', 20);
%    xlabel('method score')
%    ax=gca;
%    ax.FontSize = 16 ;
   
   
    
   
%     boxplot(tpr)
%     title('TPR')
%     xticklabels(Bins )
%     if ~saveFig
%          subplot(6);
%     else
%          saveas(gcf, fullfile('..', 'plotsTables','Crohns_TPR' ))
%          figure;
%     end
%     boxplot(fpr)
%     title('FPR')
%     xticklabels(Bins )
%     if ~saveFig
%          subplot(7);
%     else
%          saveas(gcf, fullfile('..', 'plotsTables','Crohns_FPR' ))
%          figure;
%     end
%     boxplot(tnr)
%     title('TNR')
%     xticklabels(Bins )
%     if ~saveFig
%          subplot(8);
%     else
%          saveas(gcf, fullfile('..', 'plotsTables','Crohns_TNR' ))
%          figure;
%     end
%     boxplot(fnr)
%     title('FNR')
%     xticklabels(Bins )
%     if ~saveFig
%          subplot(9);
%     else
%          saveas(gcf, fullfile('..', 'plotsTables','Crohns_FNR' ))
%          figure;
%     end
end

