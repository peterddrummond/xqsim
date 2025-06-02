function P = Mnn(a,p)
% P = MNN(a,p) generates count probabilities for an nd-fold partition
% Uses coherent projection matrix-P methods
% Requires a pure squeezed-state input for validity
% Needs:   modevec = p.part{p.nobserve};
% and      counts  = p.xk{p.nobserve};
% Returns probability of counts on list for all modes counted 
% Output is: P(m)    = g(n,m)Prod(1/m_i!)(np_i^m_i)exp(-np_i)
% where      g(n,m)  = (1+par(m)exp(2(n-np))/(1+exp(-2np)).
% Here par(m) = parity(m) = Prod Par(m_i) = Prod(1-2*rem(m_i,2))
% where np is the output photon number, n the input number
%*************************************************************

if ~isequal(p.phase,1)                                  %check phase-space is +P
    error('Mnn requires p.phase = 1, not %d',p.phase);
end
if ~(isequal(p.thermal,0)||isequal(p.thermal, zeros(1,p.modes)))                               
    error('Mnn requires p.thermal = 0, not %d',p.thermal(1)');
end
if ~isequal(p.alpha,0)                                  %check no coheemt input
    error('Mnn requires p.alpha = 0, not %d',p.alpha');
end
obs = p.nobserve;                                %store observe index
partition = p.part{obs};                         %store partition
nd = length(partition);                          %store partition dimension
ensemb = size(a,2);                              %vector ensemble size
xkcell = p.xk{obs};                              %count numbers for graph
if isempty(xkcell) || isempty(xkcell{1})         %if no count numbers input
    xkcell{nd+1} = [];
end
s = [ones(1,nd),ensemb];
n = reshape(sum(a(1+2*p.modes:3*p.modes,:),1),s);%input boson numbers
P = 1;
par = 1;                                         %initial count probability
np  = 0;
for axis = 1:nd                                  %loop over the axis number
    part  = partition{axis};                     %initialise axis partition
    count = xkcell{axis};
    if isempty(count)                            %if no count numbers input
       count = 0:length(part);                   %set counts to partition
    end                                          %end if no counts   
    ncounts = count(end)-count(1)+1;             %range for counts 
    np = np+sum(a(part,:).*a(p.modes+part,:),1); %photon numbers
    sz = [ones(1,axis-1),ncounts,ones(1,nd-axis),ensemb];
    par1 = 1-2*mod(count,2);                     %single partition parity
    if nd > 1
        par1 = reshape(par1,sz(1:nd));           %multi partition parity
    else
        par1 = reshape(par1,[sz(1),1]);           %multi partition parity
    end
    P = P.*reshape(N1(a,part,count,p),sz);
    par = par.*par1;
end
n  = reshape(n,s);
np = reshape(np,s);
P = P.*(1+par.*exp(2*(np-n)))./(1+exp(-2*n));
end                                              %end Mnn