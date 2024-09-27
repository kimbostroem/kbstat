function smd = f2smd(F,df1,df2)
%% Calculate Standardized Mean Difference (SMD) from F-value and degrees-of-freedom

smd = 2 * sqrt((F .* df1) ./ df2);

end