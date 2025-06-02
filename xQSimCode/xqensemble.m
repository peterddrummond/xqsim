function d = xqensemble(r,input)
% [d] = xqensemble(r,input) simulates one xqsim serial ensemble
% Data output is compatible with the xGRAPH graphics program
% Includes multiple cycles if p.cyc>1

% INITIALIZE CELLS 

sequence  =   length(input);                     %sequence length
p   =   input{1};                                %first sequence parameters
ensp =  p.ensembles(3);                          %parallel ensemble
ens =   p.ensembles(1);                          %vector ensemble
d  =    cell(1,sequence);                        %cell for average data
for seq = 1:sequence                             %loop over input sequence
p = input{seq};                                  %set the parameter struct
    for k = p.averages                           %loop over data types
      d{seq}{k} = zeros(p.sz{seq}{k});           %initialize averages
    end                                          %end loop over data types
end

% LOOP OVER  SERIAL ENSEMBLE  FOR ENSEMBLE AVERAGES

for r1 = 1:p.ensembles(2)                        %loop on serial ensembles
  if p.verbose > 0
    fprintf('Ensemble %d\n',(r1-1)*ensp+r);      %ensemble numbers
  end  
            
% LOOP OVER SEQUENCE AND CYCLES

  for seq = 1:sequence                           %loop over input sequence
    p = input{seq};                              %set the parameter struct
    p.seq = seq;                                 %store sequence number
    if p.phase > 1
      a = zeros(p.modes,p.ensembles(1));         %set amplitudes to zero
      if p.re > 0                                %if not normal & recycled
        q =   sqrt((p.phase-1)/4);               %quantum noise coefficient
        x  =  q*randn(p.modes,ens);              %x-quadrature noise
        y  =  q*randn(p.modes,ens);              %y-quadrature noise
        a(1:p.modes,:) = (x+1i*y);               %input vacuum noise
      end
    else
      a = zeros(3*p.modes,p.ensembles(1));       %set amplitudes to zero
    end                                          %end if not normal
    for cy = 1:p.cyc                             %loop over cycles
      p.cy = cy;                                 %store cycle number
      a = p.gen(a,p);                            %get samples
           
% LOOP OVER THE OBSERVABLES FOR ENSEMBLE AVERAGES

      for k = p.averages                         %loop over data types
        p.nobserve = k;                          %store graph number
        D1 = xqdata(a,k,p);
        D = reshape(real(D1),p.els1{seq}{k});      %reshape data 
        d{seq}{k} = reshape(d{seq}{k},p.els{seq}{k});       %reshape output 
        d{seq}{k}(cy,:,1) = d{seq}{k}(cy,:,1) + D/p.rep;    %mean data
        d{seq}{k}(cy,:,2) = 0;                              %reserved
        d{seq}{k}(cy,:,3) = d{seq}{k}(cy,:,3) + D.^2/p.rep; %mean square
        d{seq}{k} = reshape(d{seq}{k},p.sz{seq}{k});        %reshape data
      end                                        %end loop over data type
    end                                          %end cycle loop
  end                                            %end sequence loop
end                                              %end serial loop 
end                                              %end ensemble function