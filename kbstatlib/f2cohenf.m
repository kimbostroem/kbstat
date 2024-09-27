function cohenf = f2cohenf(F,df1,df2)
% calculate partial eta square from F-value and degrees-of-freedom

cohenf = sqrt((F .* df1) ./ df2);

end