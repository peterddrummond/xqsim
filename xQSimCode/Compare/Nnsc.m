function C = Nnsc(p)
% C = NNSC(p); generates comparison count probabilities for an exact
% lossless pure squeezed state n-fold partition. Only valid for pure
% squeezed states with uniform input squeezing + Identity transmission
% matrix.
if p.thermal ~= 0
    error('nnsc uses uniform pure squeezed state inputs');
elseif p.tr(1) ~= 1
    error('nnsc is for lossless distributions only');
end
partition = p.part{p.noutput};                   %current graph partition

obs = p.noutput;
xkcell = p.xk{obs};                              %count numbers for graph
if isempty(xkcell)                               %if no count numbers input
    xkcell{1} = 0:p.modes;                       %set counts to modes
end
maxm = 1;
nd = length(partition);                          %partition dimension
for i =1:nd
    maxm = max(maxm,1+max(xkcell{i}));
end
sz = ones(1,nd)*maxm;                            %size of probability data                                     
C = zeros([sz,1]);                               %initial total probability
epC = zeros(maxm,1);                             %initial even count prod.
opC = zeros(maxm,1) + log(p.cutoff);             %initial odd count prod.
p1 = zeros([1,nd]);                              %initialize numbers
M  = zeros([1,nd]);                              %initialize sizes

if iscell(partition)                             %check if cell partition
    for j =1:nd                                  %loop on partitions
        i = partition{j}(1);                     %first partition index
        p1(j)  = 1/(1+(sinh(p.sqz(i))).^2);      %thermal numbers
        M(j)  = length(partition{j});            %partition length
    end                                          %end loop on partitions
else                                             %if partition vector
    i = 1;                                       %first detector index
    M = partition;                               %partition lengths
    for j=1:nd                                   %loop on vector index         
        p1(j) = 1/(1+(sinh(p.sqz(i))).^2);       %thermal numbers
        i=i+M(j);                                %update detector index
    end                                          %end loop on vector
end                                              %end check if cell

M = M*0.5;                                       %Change mode number to M/2
for j =1:nd                                      %loop on partitions
    n1=maxm^(j-1);                               %array size below j
    n2=maxm^(nd-j);                              %array size above j
    epC(1,:) = epC(1,:) + M(j)*log(p1(j));       %m = 0 bin. coeff. product
    C = reshape(C,n1,maxm,n2);                   %reshape log prob.
    C(:,1,:) = C(:,1,:) + epC(1,:);              %Add m = 0 count 
    lp2 = log(1-p1(j));                          %get log(1-n) for j
    for m = 1:(maxm)*0.5                         %start loop on counts
        epC(m+1,:) = epC(m,:)+log(1+(M(j)-1)/m)+lp2;%ven bin. prod.
        C(:,2*m,:) = C(:,2*m,:) + opC(m,:);       %add odd count prob.
        if m < (maxm)*0.5
            C(:,2*m+1,:) = C(:,2*m+1,:)+ epC(m+1,:);%add even count prob.
        end
    end
epC = zeros(maxm,1);
end                                              %start loop on counts

C = reshape(C,[sz,1]);                           %reshape log prob.
C = exp(C);                                      %probability
testC = size(C)
in = [p.xk{p.noutput}{:},{0},{0},{0},{0},{0}];
testin = in{1}

C = C(1+in{1}(1):1+in{1}(end),1+in{2}(1):1+in{2}(end),...
    1+in{3}(1):1+in{3}(end),1+in{4}(1):1+in{4}(end),...
    1+in{5}(1):1+in{5}(end),1+in{6}(1):1+in{6}(end));%up to 6D binning
end                                              %end nnc
