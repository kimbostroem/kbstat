function str = capitalize(str)

str = char(str);
idx = regexp([' ' str],'(?<=\s+)\S','start')-1;
str(idx) = upper(str(idx));

end