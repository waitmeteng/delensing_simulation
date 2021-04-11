function [map_lens]=lensing1(resol,angu_size,dx_angle,dy_angle,map)
% lensing (ray-tracing)
% -------------------------------------------------------------------------
% INPUT:
% .resol                 : resolution of lensed maps.
%                          This value should be less than the resolution
%                          of input unlensed map. 
% .angu_size             : angular size of maps (degree)
% .dx_angle              : deflection angle in x direction
% .dy_angle              : deflection angle in y direction
% .tmap                  : unlensed map
%--------------------------------------------------------------------------
% OUTPUT:
% .map_lens             : lensed map
%--------------------------------------------------------------------------
% Wei-Hsiang Teng, NTU, 2011
%--------------------------------------------------------------------------
reso=size(dx_angle,1);
delta_x=angu_size/reso*pi/180;
if (size(dy_angle,1)~=reso || size(map,1)~=reso)
   disp('The size of input maps are not the same!'); return 
end

if (resol>reso)
   disp('warning:The pixel number of lensed map is supposed to be less than unlensed map!\n');
   disp('the pixel number of lensed map will be the same with that of unlensed map.');
   resoll=reso;
else
   resoll=resol; 
end

%--------------------------------------------------------------------------    
    % If reso is too large suffering memory problem use below part...
%--------------------------------------------------------------------------  

    if (floor(reso/resoll)~=reso/resoll)
       disp('error: the ratio of lensed map pixel number and unlensed map should\n')
       disp('be integer.');return 
    end
    
    add=16;
    
    % create the larger maps corresponding to input unlensed maps
    % using periodic boundary condition.
    map_unlens_=zeros(reso+add,reso+add);
    
    map_unlens_(1+add/2:reso+add/2,1+add/2:reso+add/2)=map;
    
    map_unlens_(1+add/2:reso+add/2,1:add/2)=map(:,reso-add/2+1:reso);
    
    map_unlens_(1+add/2:reso+add/2,reso+1+add/2:reso+add)=map(:,1:add/2);
    
    map_unlens_(1:add/2,:)=map_unlens_(reso+1:reso+add/2,:);
    
    map_unlens_(reso+add/2+1:reso+add,:)=map_unlens_(add/2+1:add,:);
    
    clear map 
        
    map_lens=zeros(resoll,resoll);
    [xinterp,yinterp]=meshgrid(1:reso/resoll:reso/2);   
    [xinterp1,yinterp1]=meshgrid(-(add/2-1):1:reso/2+add/2);
    
    for ii=1:4
    
    switch ii
        case 1 % leftup part
            temp_x=xinterp-dx_angle(1:reso/resoll:reso/2,1:reso/resoll:reso/2)/delta_x; 
            temp_y=yinterp-dy_angle(1:reso/resoll:reso/2,1:reso/resoll:reso/2)/delta_x;

            map_lens(1:resoll/2,1:resoll/2)= ...
                interp2(xinterp1,yinterp1,map_unlens_(1:reso/2+add,1:reso/2+add),temp_x,temp_y,'spline');
            
        case 2 % rightup part
            temp_x=xinterp-dx_angle(1:reso/resoll:reso/2,reso/2+1:reso/resoll:reso)/delta_x; 
            temp_y=yinterp-dy_angle(1:reso/resoll:reso/2,reso/2+1:reso/resoll:reso)/delta_x;
            
            map_lens(1:resoll/2,resoll/2+1:resoll)= ...
                interp2(xinterp1,yinterp1,map_unlens_(1:reso/2+add,reso/2+1:reso+add),temp_x,temp_y,'spline');
            
        case 3 % leftdown part
            temp_x=xinterp-dx_angle(reso/2+1:reso/resoll:reso,1:reso/resoll:reso/2)/delta_x; 
            temp_y=yinterp-dy_angle(reso/2+1:reso/resoll:reso,1:reso/resoll:reso/2)/delta_x;
            
            map_lens(resoll/2+1:resoll,1:resoll/2)= ...
                interp2(xinterp1,yinterp1,map_unlens_(reso/2+1:reso+add,1:reso/2+add),temp_x,temp_y,'spline');
            
        case 4 % rightdown part
            temp_x=xinterp-dx_angle(reso/2+1:reso/resoll:reso,reso/2+1:reso/resoll:reso)/delta_x; 
            temp_y=yinterp-dy_angle(reso/2+1:reso/resoll:reso,reso/2+1:reso/resoll:reso)/delta_x;
            
            map_lens(resoll/2+1:resoll,resoll/2+1:resoll)= ...
                interp2(xinterp1,yinterp1,map_unlens_(reso/2+1:reso+add,reso/2+1:reso+add),temp_x,temp_y,'spline');
           
    end
    end
