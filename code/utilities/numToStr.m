function str = numToStr(num)
str = num2str(round(num,2),2);
num = str2num(str);
if num > 0 && num < 1
    str = str(2:end);
elseif num<0 && num > -1
    str = ['-',str(3:end)];
end
end