function y = myCeil(x, range)
p = -(round(log10(range))-1);
y = ceil(x*(10^p))/(10^p);  
end