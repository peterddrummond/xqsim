function [] =  xmultiplot(n,datan,param,g) 
%   [] = xmultiplot(n,data,axis,in) plots multidimensional data files
%   The first or time axis is plotted from data in 'axis' 
%   Input:  number 'n', data & axis array 'datan', 'axis',parameters 'g'.
%   Output: graphs, with multiple plots giving different dimensional views.
%   Called by: xgraph
%   Needs:     xreduce, ximageplot, xtransverseplot, xplot3, xplot, xheader
%   xGRAPH functions are licensed by Peter D. Drummond, (2021) - see License
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  REDUCE THE INPUT DATA TO 3D
%                          
set (0, 'defaultaxesfontsize', g.font{n})
set (0, 'defaulttextfontsize', g.font{n})     
head = g.headers{n};                         %%set graph name    
[datan,param,ax,nx,x,xlab] =  xreduce(n,datan,param,g); %%reduce data
grd = length(nx)-2;                          %%graph dimension = grd
tlabel = xlab{1};
tcoord = x{1};
olab = g.olabels{n};
if grd == 0
        fprintf('xGRAPH warning: nothing to plot in function{%d}\n',n);
        return
end
mx = 1+floor(nx/2);                          %%adjust midpoint: space 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PLOT DATA IN 3D
%
if grd > 2                                   %Check if dimension > 2
    ximageplot(n,datan,nx,x,xlab,olab,ax,g);%%Image plot
    sz = size(datan);
    fprintf('Lattice reduced: %s = %d\n',xlab{3},x{3}(mx(4))); 
    datan=reshape(datan(:,:,:,mx(4),:),[sz(1:3),sz(5)]);%%Central point
    sz = size(param);
    param=reshape(param(:,:,:,mx(4),:),[sz(1:3),sz(5)]);%%Central point 
end                                          %%End check if dimension > 2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PLOT DATA IN 2D
%        
if grd  > 1                                  %%Check plot dimension >1
    xtransverseplot(n,datan,param,nx,x,xlab,olab,ax(2:end),g) %%Transverse
    im = reshape(datan(1,:,:,1),nx(2),nx(3));%%Image vs t,x
    fprintf('3D plot in graph %d\n',n);
    xplot3(x{1},x{2},im',xlab{1},xlab{2},olab,head,n,ax,g);
    sz = size(datan);
    fprintf('Lattice reduced: %s = %d\n',xlab{2},x{2}(mx(3)));
    datan=reshape(datan(:,:,mx(3),:),[sz(1:2),sz(4)]);%%Reduce
    sz = size(param);
    param=reshape(param(:,:,mx(3),:),[sz(1:2),sz(4)]);
end                                          %%End check if dimension > 1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  PLOT DATA IN 1D
%            
if g.parametric{n}(2) == 1
    tlabel = g.olabels{g.parametric{n}(1)};
    tcoord = param ;                         %%Axis data 
end
fprintf('2D plot in graph %d\n',n);
xplot(tcoord,datan,n,ax,g);                  %%xspde 2D plot  
xheader(head,tlabel,olab,' ');               %%xspde 2D plot title          
end                                          %%end function
