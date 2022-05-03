s = load('../results/Crohns.mat');
perfs = s.Results.Crohns.perfs;
ds = s.Results.Crohns.ds;
rMeasures = rankingMeasures(ds);
perfs = sortPerfs(perfs,rMeasures);
cls = perfs(1).classes;
scs = perfs(1).scores;
prior = mean(cls);
for ii = 1:1000    
    [scs, ix] = randomSort(scs);
    cls = cls(ix);
    p = posteriorWindow(cls, scs);
    llrPlus = (p./(1-p))*(1-prior)/prior;
    plot(scs, llrPlus)
    hold on
    disp(llrPlus(end))
    if llrPlus(end) < 6
        disp('hi')
    end
end
hold off