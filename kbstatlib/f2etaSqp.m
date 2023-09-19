function etaSqP = f2etaSqp(F,df1,df2)
% calculate partial eta square from F-value and degrees-of-freedom

etaSqP = (F .* df1) ./ (F .* df1 + df2);

end