function [index] = getindex2(D)
%Given a monomial x^\alpha, compute the index of the corresponding moment
%Input: 
%   l: the exponent alpha in the monomial x^\alpha
%Ouput:
%   index: the index (position) of the corresponding moment in the moment
%   sequence

s2=factorial(size(D,2)+sum(D,2))./(factorial(sum(D,2))*factorial(size(D,2)));


index2=zeros(size(s2));
for i=size(D,2):-1:2
    di=sum(D,2)-sum(D(:,i+1:size(D,2)),2);
    index2=index2+factorial(i-1+di)./(factorial(i-1)*factorial(di))-factorial(i-1+di-D(:,i))./(factorial(di-D(:,i))*factorial(i-1));
end
index=round(s2-index2);


end

