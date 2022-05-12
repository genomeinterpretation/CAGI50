function y = myFloor(x, range)
p = -(round(log10(range))-1);
y = floor(x*(10^p))/(10^p);  
end