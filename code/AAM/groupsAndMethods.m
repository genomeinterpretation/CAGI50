function [M1, G1, M2, G2, invert] = groupsAndMethods()

GG2 = struct();
invert={'SIFT', 'PROVEAN', 'FATHMM'};
groups = {'SIFT', 'PROVEAN', 'PolyPhen2', 'LRT', 'MutationTaster', 'MutationAssessor',...
    'FATHMM', 'fathmm_MKL', 'CADD', 'VEST', 'fitCons', 'DANN', 'MetaSVM', 'MetaLR', ...
	'GenoCanyon', 'Eigen' 'M_CAP', 'REVEL', 'MutPred','GERP', ...
    'phyloP', 'phastCons', 'SiPhy'};

GG2.PolyPhen2 = {'Polyphen2_HDIV_score', 'Polyphen2_HVAR_score'};
GG2.fathmm_MKL = {'fathmm_MKL_coding_score'};
GG2.CADD = {'CADD_raw'};
GG2.fitCons = {'integrated_fitCons_score', 'GM12878_fitCons_score', 'H1_hESC_fitCons_score', 'HUVEC_fitCons_score'};
GG2.Eigen= {'Eigen_raw', 'Eigen_PC_raw'};
GG2.VEST= {'VEST3_score'};
GG2.SiPhy= {'SiPhy_29way_logOdds'};
GG2.phyloP = {'phyloP20way_mammalian','phyloP100way_vertebrate'};
GG2.phastCons = {'phastCons100way_vertebrate', 'phastCons20way_mammalian'};
GG2.GERP = {'GERP___RS'};


M2 = {};
G2 = {};
for g = groups
    grp = g{1};
    if ~isfield(GG2, grp)
       M2 = [M2, {[grp,'_score']}];
       G2 = [G2, {grp}];
    else
       M2 = [M2, GG2.(grp)];    
       G2 = [G2, repmat({grp}, 1, length(GG2.(grp)))];
    end
end

M1 =  {'VEST4', 'MutPred', 'Bologna', 'Turkey', 'SIFT', 'PolyPhen', 'Condel', 'CADD_v1_4_RawScore', ...
    'DANN_score', 'Eigen_PC_raw_coding', 'Eigen_raw_coding', 'fathmm_MKL_coding_score', 'GenoCanyon_score', ...
    'LRT_score', 'M_CAP_score', 'MetaLR_score', 'MetaSVM_score', 'REVEL_score','pph2_prob'};
GG1 = struct();
GG1.VEST = {'VEST4'};
GG1.PolyPhen2 = {'PolyPhen','pph2_prob'};
G1 = {};

gg1 = fieldnames(GG1);
for m  = M1
    ix = find(cellfun(@(g) contains(m{1}, g), groups));
    if ~isempty(ix)
        G1 = [G1, groups(ix)];
    else
        groupFound = false;
        for ff = 1:length(gg1)
            f =  gg1{ff};
            if any(strcmpi(m{1}, GG1.(f)))
                groupFound = true;
                G1 = [G1, {f}];
            end
        end
        if ~groupFound
            G1 = [G1, m];
        end
    end
end

end