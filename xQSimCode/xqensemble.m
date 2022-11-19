function d = xqensemble(r,input)
% [d] = xqensemble(r,input) simulates one qcpsim serial ensemble
% Data output is compatible with the xGRAPH graphics program
% Includes multiple cycles if p.cyc>1

% INITIALIZE CELLS 

ls  =   length(input);                           %sequence length
p   =   input{1};                                %first sequence parameters
ensp =  p.ensembles(3);                          %parallel ensemble
ens =   p.ensembles(1);                          %vector ensemble
d  =    cell(ls);                                %cell for average data
for s = 1:ls                                     %loop over input sequence
p = input{s};                                    %set the parameter struct
    for k = 1:p.graphs                           %loop over data types
      d{s}{k} = zeros(p.sz{s}{k});               %initialize averages
    end                                          %end loop over data types
end

% LOOP OVER  SERIAL ENSEMBLE  FOR ENSEMBLE AVERAGES

for r1 = 1:p.ensembles(2)                        %loop on serial ensembles
  if p.print > 0
    fprintf('Ensemble %d\n',(r1-1)*ensp+r);      %print ensemble numbers
  end
  a = zeros(2*p.M,p.ensembles(1));               %set amplitudes to zero
  if p.method > 1 && p.re > 0                    %if not normal & recycled
      q=sqrt((p.method-1)/4);                    %quantum noise coefficient
      x  =  q*randn(p.M,ens);                    %x-quadrature noise
      y  =  q*randn(p.M,ens);                    %y-quadrature noise
      a(1:p.M,:) = (x+1i*y);                     %input vacuum noise
      a(1+p.m:2*p.M,:) = (x-1i*y);               %conjugate vacuum noise
  end                                            %end if not normal 
            
% LOOP OVER SEQUENCE AND CYCLES

  for s = 1:ls                                   %loop over input sequence
    p = input{s};                                %set the parameter struct
    p.seq = s;                                   %store sequence number
    for cy = 1:p.cyc                             %loop over cycles
      p.cy = cy;                                 %store cycle number
      a = p.gen(a,p);                            %get samples
           
% LOOP OVER THE OBSERVABLES FOR ENSEMBLE AVERAGES

      for k = 1:p.graphs                         %loop over data types
        p.k = k;                                 %store graph number
        D1 = p.observe{k}(a,p);                  %get vector average
        D = reshape(real(D1),p.els1{s}{k});      %reshape data 
        d{s}{k} = reshape(d{s}{k},p.els{s}{k});  %reshape output 
        d{s}{k}(cy,:,1) = d{s}{k}(cy,:,1) + D/p.rep;     %mean data
        d{s}{k}(cy,:,2) = 0;                             %reserved
        d{s}{k}(cy,:,3) = d{s}{k}(cy,:,3) + D.^2/p.rep;  %mean square
        d{s}{k} = reshape(d{s}{k},p.sz{s}{k});   %reshape data
      end                                        %end loop over data type
    end                                          %end cycle loop
  end                                            %end sequence loop
end                                              %end serial loop 
end                                              %end ensemble function
