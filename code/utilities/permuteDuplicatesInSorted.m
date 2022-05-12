function [x,ix] = permuteDuplicatesInSorted(x)
%Assumes x is sorted. Permutes the duplicates in sorted x
i = 1;
ix = 1:length(x);
while i <= length(x)
    j=1;
    while i+j <= length(x) && x(i+j) == x(i)
        j = j+1;
    end
    if j > 1
        ixx = i:(i+j-1);
        ix(ixx) = ixx(randperm(length(ixx)));
    end
    i=i+j;
end
end

