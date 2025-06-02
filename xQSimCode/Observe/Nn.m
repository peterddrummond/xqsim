function P = Nn(a,p)
% P = NN(a,p) generates count probabilities for an nd-fold partition
% The partition is defined by the quantity p.part{p.nobserve}
% This is  a cell array of index vectors
% The count range needed is defined by p.xk{p.nobserve}
%*************************************************************

if p.phase ~= 1                                  %check phase-space is +P
    error('Nn requires +P method = 1');
end
obs = p.nobserve;                                %store observe index
partition = p.part{obs};                         %store partition
nd = length(partition);                          %store partition dimension
ensemb = size(a,2);                              %first ensemble size
xkcell = p.xk{obs};                              %count numbers for graph
if isempty(xkcell) || isempty(xkcell{1})         %if no count numbers input
    xkcell{nd+1} = [];
end
P = 1;                                           %initial count probability
for axis = 1:nd                                  %loop over the axis number
    part  = partition{axis};                     %initialise axis partition
    count = xkcell{axis};
    if isempty(count)                            %if no count numbers input
       count = 0:length(part);                   %set counts to partition
    end                                          %end if no counts
    ncounts = count(end)-count(1)+1;
    sz = [ones(1,axis-1),ncounts,ones(1,nd-axis),ensemb];
    P = P.*reshape(N1(a,part,count,p),sz);
end
end                                              %end kn