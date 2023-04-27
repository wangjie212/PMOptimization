function degs = deglist(nvar, mindeg, maxdeg)
%Given the number nvar of variables, 
%compute the list of exponets of monomials
%of degree from min to max



if mindeg < 1
    first_ind = 0;
else
    first_ind = nchoosek(mindeg-1+nvar, nvar);
end

if nvar == 1
    degs = (mindeg:maxdeg)';
    return;
end

degs = zeros(nchoosek(maxdeg+nvar,nvar)-first_ind, nvar);
degs(end-maxdeg:end,1:2) = [(0:maxdeg);(maxdeg:-1:0)]';

s = maxdeg + 1;
for i=2:nvar
    t = s; 
    for j=maxdeg-1:-1:((i==nvar)*mindeg)
        t = t*(j+1)/(i+j); 
        degs(end-s-t+1:end-s,1:i-1) = degs(end-s+1:end-s+t,1:i-1);
        degs(end-s-t+1:end-s,i) = degs(end-s+1:end-s+t,i) - 1;
        if i < nvar
            degs(end-s-t+1:end-s,i+1) = maxdeg - j;
        end
        s = s + t;
    end
end



   
    
