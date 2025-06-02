function C = Kn(a,p)
% C = KN(a,p) gives click probabilities for an n-fold partition
% The partition is defined by the quantity p.part{p.nobserve}
% This is a cell array of index vectors
% The count range needed is defined by p.xk{p.nobserve}
%*************************************************************

if p.phase ~= 1                                  %check phase-space is +P
    error('xqsim supports click probabilities only for +P method = 1');
end
obs = p.nobserve;                                %store observe index
partition = p.part{obs};                         %store partition
nd = length(partition);                          %store partition dimension
ensemb = size(a,2);                              %first ensemble size
xkcell = p.xk{obs};                              %count numbers for graph
if isempty(xkcell)                               %if no count numbers input
    fprint('\nWarning: no axis numbers input for graph %d\n', obs);
    xkcell{nd+1} = [];
end
C = 1;                                          %initial count probability
for axis = 1:nd                                  %loop over the axis number
    part  = partition{axis};                     %initialise axis partition
    xk = xkcell{axis};
    if isempty(xk)                               %if no count numbers input
       xk = 0:length(part);                      %set counts to partition
    end                                          %end if no counts
    ncounts = xk(end)-xk(1)+1;
    sz = [ones(1,axis-1),ncounts,ones(1,nd-axis),ensemb];
    C1 = K1(a,part,xk,p);
    C1 = reshape(C1,sz);
    C = C.*C1;
end
end                                              %end kn