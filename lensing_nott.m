function [qmap_lens,umap_lens]=lensing_nott(resol,number, angu_size,dx_angle,dy_angle,qmap,umap)
% lensing (ray-tracing)
% -------------------------------------------------------------------------
% INPUT:
% .resol                 : resolution of lensed maps.
%                          This value should be less than the resolution
%                          of input unlensed map. 
% .angu_size             : angular size of maps (degree)
% .dx_angle              : deflection angle in x direction
% .dy_angle              : deflection angle in y direction
% .tmap                  : unlensed T map
% .qmap                  : unlensed Q map
% .umap                  : unlensed U map
%--------------------------------------------------------------------------
% OUTPUT:
% .tmap_lens             : lensed T map
% .qmap_lens             : lensed Q map
% .umap_lens             : lensed U map
%--------------------------------------------------------------------------
% Wei-Hsiang Teng, NTU, 2010
%--------------------------------------------------------------------------
reso=size(dx_angle,1);
delta_x=angu_size/reso*pi/180;
if (size(dy_angle,1)~=reso || size(qmap,1)~=reso || size(umap,1)~=reso)
   disp('The size of input maps are not the same!'); return 
end

if (resol>reso)
   disp('warning:The pixel number of lensed map is supposed to be less than unlensed map!\n');
   disp('the pixel number of lensed map will be the same with that of unlensed map.');
   resoll=reso;
else
   resoll=resol; 
end
% if (floor(reso/resoll)==reso/resoll)
%     [xinterp,yinterp]=meshgrid(1:reso/resoll:reso);
%     temp_x=xinterp-dx_angle(1:reso/resoll:reso,1:reso/resoll:reso)/delta_x; 
%     temp_y=yinterp-dy_angle(1:reso/resoll:reso,1:reso/resoll:reso)/delta_x;
%     clear xinterp yinterp
%     add=16;
%     
%     tmap_unlens_=zeros(reso+add,reso+add);
%     qmap_unlens_=zeros(reso+add,reso+add);
%     umap_unlens_=zeros(reso+add,reso+add);
%     
%     tmap_unlens_(1+add/2:reso+add/2,1+add/2:reso+add/2)=tmap;
%     qmap_unlens_(1+add/2:reso+add/2,1+add/2:reso+add/2)=qmap;
%     umap_unlens_(1+add/2:reso+add/2,1+add/2:reso+add/2)=umap;
%     
%     tmap_unlens_(1+add/2:reso+add/2,1:add/2)=tmap(:,reso-add/2+1:reso);
%     qmap_unlens_(1+add/2:reso+add/2,1:add/2)=qmap(:,reso-add/2+1:reso);
%     umap_unlens_(1+add/2:reso+add/2,1:add/2)=umap(:,reso-add/2+1:reso);
%     
%     tmap_unlens_(1+add/2:reso+add/2,reso+1+add/2:reso+add)=tmap(:,1:add/2);
%     qmap_unlens_(1+add/2:reso+add/2,reso+1+add/2:reso+add)=qmap(:,1:add/2);
%     umap_unlens_(1+add/2:reso+add/2,reso+1+add/2:reso+add)=umap(:,1:add/2);
%     
%     tmap_unlens_(1:add/2,:)=tmap_unlens_(reso+1:reso+add/2,:);
%     qmap_unlens_(1:add/2,:)=qmap_unlens_(reso+1:reso+add/2,:);
%     umap_unlens_(1:add/2,:)=umap_unlens_(reso+1:reso+add/2,:);
%     
%     tmap_unlens_(reso+add/2+1:reso+add,:)=tmap_unlens_(add/2+1:add,:);
%     qmap_unlens_(reso+add/2+1:reso+add,:)=qmap_unlens_(add/2+1:add,:);
%     umap_unlens_(reso+add/2+1:reso+add,:)=umap_unlens_(add/2+1:add,:);
%     [xinterp1,yinterp1]=meshgrid(-(add/2-1):1:reso+add/2);
%     clear tmap qmap umap
% 
%     %==========================================================================
%     tmap_lens=interp2(xinterp1,yinterp1,tmap_unlens_,temp_x,temp_y,'spline');
%     qmap_lens=interp2(xinterp1,yinterp1,qmap_unlens_,temp_x,temp_y,'spline');
%     umap_lens=interp2(xinterp1,yinterp1,umap_unlens_,temp_x,temp_y,'spline');
% else
%     disp('error: the ratio of lensed map pixel number and unlensed map should\n')
%     disp('be integer.');return 
% end

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
%     tmap_unlens_=zeros(reso+add,reso+add);
    qmap_unlens_=zeros(reso+add,reso+add);
    umap_unlens_=zeros(reso+add,reso+add);
    
%     tmap_unlens_(1+add/2:reso+add/2,1+add/2:reso+add/2)=tmap;
    qmap_unlens_(1+add/2:reso+add/2,1+add/2:reso+add/2)=qmap;
    umap_unlens_(1+add/2:reso+add/2,1+add/2:reso+add/2)=umap;
    
%     tmap_unlens_(1+add/2:reso+add/2,1:add/2)=tmap(:,reso-add/2+1:reso);
    qmap_unlens_(1+add/2:reso+add/2,1:add/2)=qmap(:,reso-add/2+1:reso);
    umap_unlens_(1+add/2:reso+add/2,1:add/2)=umap(:,reso-add/2+1:reso);
    
%     tmap_unlens_(1+add/2:reso+add/2,reso+1+add/2:reso+add)=tmap(:,1:add/2);
    qmap_unlens_(1+add/2:reso+add/2,reso+1+add/2:reso+add)=qmap(:,1:add/2);
    umap_unlens_(1+add/2:reso+add/2,reso+1+add/2:reso+add)=umap(:,1:add/2);
    
    clear qmap umap
    
%     tmap_unlens_(1:add/2,:)=tmap_unlens_(reso+1:reso+add/2,:);
    qmap_unlens_(1:add/2,:)=qmap_unlens_(reso+1:reso+add/2,:);
    umap_unlens_(1:add/2,:)=umap_unlens_(reso+1:reso+add/2,:);
    
%     tmap_unlens_(reso+add/2+1:reso+add,:)=tmap_unlens_(add/2+1:add,:);
    qmap_unlens_(reso+add/2+1:reso+add,:)=qmap_unlens_(add/2+1:add,:);
    umap_unlens_(reso+add/2+1:reso+add,:)=umap_unlens_(add/2+1:add,:);
            
%     tmap_lens=zeros(resoll,resoll);
    qmap_lens=zeros(resoll,resoll);
    umap_lens=zeros(resoll,resoll);
    
   
    [xinterp,yinterp]=meshgrid(1:reso/resoll:reso/number);   
    [xinterp1,yinterp1]=meshgrid(-(add/2-1):1:reso/number+add/2);
    
    for ii=1:number
        for jj=1:number
    
    
            temp_x=xinterp-dx_angle(reso/number*(ii-1)+1:reso/resoll:reso/number*ii,reso/number*(jj-1)+1:reso/resoll:reso/number*jj)/delta_x; 
            temp_y=yinterp-dy_angle(reso/number*(ii-1)+1:reso/resoll:reso/number*ii,reso/number*(jj-1)+1:reso/resoll:reso/number*jj)/delta_x;

            qmap_lens(reso/number*(ii-1)+1:reso/resoll:reso/number*ii,reso/number*(jj-1)+1:reso/resoll:reso/number*jj)= ...
                interp2(xinterp1,yinterp1,qmap_unlens_(reso/number*(ii-1)+1:reso/number*ii+add,reso/number*(jj-1)+1:reso/number*jj+add),temp_x,temp_y,'spline');
            umap_lens(reso/number*(ii-1)+1:reso/resoll:reso/number*ii,reso/number*(jj-1)+1:reso/resoll:reso/number*jj)= ...
                interp2(xinterp1,yinterp1,umap_unlens_(reso/number*(ii-1)+1:reso/number*ii+add,reso/number*(jj-1)+1:reso/number*jj+add),temp_x,temp_y,'spline'); 
            
        
        end
    end
