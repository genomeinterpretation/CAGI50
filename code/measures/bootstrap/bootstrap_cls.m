function out = bootstrap_cls(fnc, classes, scores, arg, nboot, cheaper)
%% Bootstrap function for the measures computed on the classification dataset.
%% Use this function to bootstrap all classification and clinical measures.
if nargin < 6
    cheaper = true;
end
[~,ncol] = size(scores);

% Result using the entire dataset without bootstrap. For Experimental-Max
% it contains the result for the first column in 'scores'.
n = length(classes);
nScore_u = length(unique(scores(:,1)));
is_discrete = nScore_u <= 2;
[scores_srt, ix_srt] = randomSort(scores(:,1));
classes_srt = classes(ix_srt);
arg1 = arg;
arg1.presorted = true;
res1 = fnc(classes_srt, scores_srt(:,1), arg1);

% The scalar Metrics is used to collect the bootstrap results of the scalar outputs
% of the called function.
metrics = struct();
metrics_v =struct();

if isstruct(res1)
    allFields = fieldnames(res1);
end
index = 1:n;
if ncol == 1
    ix_pos = index(classes_srt==1);
    ix_neg = index(classes_srt==0);
    i = 1;
    while i <= nboot
        ix_pos_bt = bootstrp(1, @(x) x, ix_pos);
        ix_neg_bt = bootstrp(1, @(x) x, ix_neg);
        ix_bt = [ix_pos_bt';ix_neg_bt'];
        
        scores_bt = scores_srt(ix_bt);
        %classes_bt = [ones(length(ix_pos_bt),1); zeros(length(ix_neg_bt), 1)];
        classes_bt = classes_srt(ix_bt);
        [scores_bt, ix_sort] = randomSort(scores_bt);
        nScore_bt_u = length(unique(scores_bt));
        if (~is_discrete  || nScore_bt_u == nScore_u) && (is_discrete  || nScore_bt_u > 2) 
            classes_bt = classes_bt(ix_sort);
            ix_bt = ix_bt(ix_sort);
            arg.presorted = true;
            res = fnc(classes_bt, scores_bt, arg);
            updateMetrics(res,i, nboot);
            %updateVectorMetrics(res,i, nboot, ix_bt, res1);
            i = i + 1;
        end
    end
else
    for i = 1:ncol
       arg.presorted = false;
       res = fnc(classes, scores(:,i), arg);
       updateMetrics(res,i, ncol);
    end
end

out = struct();

%Create the output datastructure containing summary statistics for the
%scalar results (obtained from the bootstraps). For the nonscalar values
%the results in res1 are stored. When the score column contains multiple
%columns the mean result is stored in the val field of scalar metrics.
if isstruct(res1)
scalarFields = fieldnames(metrics);
%vectorFields = fieldnames(metrics_v);
for j = 1:length(allFields)
    fld = allFields{j};
    out.(fld).val = res1.(fld);
    if any(strcmp(fld, scalarFields))
        out.(fld).bt.mean = mean(metrics.(fld), 'omitnan');
        out.(fld).bt.median = median(metrics.(fld));
        out.(fld).bt.p5 = prctile(metrics.(fld),5);
        out.(fld).bt.p95 = prctile(metrics.(fld),95);
        out.(fld).bt.std = std(metrics.(fld), 'omitnan');
        if ~cheaper
            out.(fld).bt.bt = metrics.(fld);
        end
        if ncol > 1
           out.(fld).val =  out.(fld).bt.mean;
        end
%     elseif any(strcmp(fld, vectorFields))
%         out.(fld).mean = metrics_v.(fld).mean;
%         out.(fld).cnt = metrics_v.(fld).cnt;
    end
end
else
    out.bt.mean = mean(metrics);
    out.bt.median = median(metrics);
    out.bt.p5 = prctile(metrics,5);
    out.bt.p95 = prctile(metrics,95);
    out.bt.std = std(metrics);
    if ~cheaper
        out.bt.bt = metrics;
    end
    if ncol == 1
        out.val = res1;
    else
        out.val = out.bt.mean;
    end
end




    function updateMetrics(res, b, B)
        %records the bootstrap samples from the scalar metric values
        %contained in res1
        if isstruct(res)
            for ii = 1:length(allFields)
                fld1 = allFields{ii};
                r = res.(fld1);
                if  ~iscell(r) && (isscalar(r) || (length(r) < n)) 
                    if b == 1
                        metrics.(fld1) = nan(B,length(r));
                    end
                    %disp(fld1)
                    if ~isfield(metrics, fld1) || length(metrics.(fld1)(b,:)) - length(r) ~=0
                        error('dimensions expected to match')
                    end
                    metrics.(fld1)(b,:) = r;
                end
            end
        else
            if b == 1
                metrics = nan(B,length(res));
            end
            metrics(b,:) = res;
        end
    end

    
end