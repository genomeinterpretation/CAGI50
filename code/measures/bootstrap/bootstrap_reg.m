function out = bootstrap_reg(fnc, y, predictions, nboot)
% Number of bootstrap samples
[~,ncol] = size(predictions);
cheaper =true;
%Result using the entire dataset without bootstrap contained in the first
res1 = fnc(y, predictions(:,1));

% The scalar Metrics is used to collect the bootstrap results of the scalar outputs
% of the called function.
metrics = struct();

if isstruct(res1)
    allFields = fieldnames(res1);
end

if ncol == 1
    for i = 1:nboot
        ix = bootstrp(1, @(x) x, 1:length(y));
        yy = y(ix);
        preds = predictions(ix);
        res = fnc(yy, preds);
        updateMetrics(res,i, nboot);
    end
else
    for i = 1:ncol
       res = fnc(y, predictions(:,i));
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
    end
end
else
    out.bt.mean = mean(metrics, 'omitnan');
    out.bt.median = median(metrics);
    out.bt.p5 = prctile(metrics,5);
    out.bt.p95 = prctile(metrics,95);
    out.bt.std = std(metrics, "omitnan");
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
                if isscalar(r) || ((length(r) < length(y)) && ~iscell(r))
                    if b == 1
                        metrics.(fld1) = nan(B,length(r));
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