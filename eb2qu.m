function [ qmap, umap ] = eb2qu(angu_size, emap, bmap)
% this code is to tranfer E, B maps into Q, U maps. You can choose if you 
% want to use Wiener filter according to the input power spectra and noise 
% level. 
% -------------------------------------------------------------------------
% INPUT:
% .angu_size                : angular size of E, B maps
% .emap                     : input E map
% .bmap                     : input B map
% OUTPUT:
% .qmap                     : output Q map
% .umap                     : output U map
%--------------------------------------------------------------------------
reso=size(emap,1);
if (size(bmap,1)~=reso)
   disp('The size of input maps are not the same!'); return 
end
delta_x=angu_size/reso*pi/180;
%--------------------------------------------------------------------------
cos2phi_m=zeros(reso/2+1,reso);
sin2phi_m=zeros(reso/2+1,reso);

for i=1:reso
   for j=1:reso/2+1
       kx=i-(reso/2+1);
       ky=j-(reso/2+1); 
       
       phi=atan2(ky,kx);
       cos2phi_m(j,i)=cos(2*phi); 
       sin2phi_m(j,i)=sin(2*phi);       
       
   end
end 
qmap_l=zeros(reso,reso);
umap_l=zeros(reso,reso);
emap_l=fftshift(fft2(emap))*delta_x^2; 
bmap_l=fftshift(fft2(bmap))*delta_x^2;

qmap_l(1:reso/2+1,:)= ...
    (-bmap_l(1:reso/2+1,:).*sin2phi_m+emap_l(1:reso/2+1,:).*cos2phi_m);
umap_l(1:reso/2+1,:)= ...
    (bmap_l(1:reso/2+1,:).*cos2phi_m+emap_l(1:reso/2+1,:).*sin2phi_m);

qmap_l(reso/2+1,reso/2+1)=0;
qmap_l(1,1)=real(qmap_l(1,1));
qmap_l(1,reso/2+1)=real(qmap_l(1,reso/2+1));
qmap_l(reso/2+1,1)=real(qmap_l(reso/2+1,1));
qmap_l(reso/2+2:reso,1)=conj(qmap_l(reso/2:-1:2,1));
qmap_l(1,reso/2+2:reso)=conj(qmap_l(1,reso/2:-1:2));
qmap_l(reso/2+1:reso,reso/2+1:reso)=conj(qmap_l(reso/2+1:-1:2,reso/2+1:-1:2));
qmap_l(reso/2+2:reso,2:reso/2)=conj(qmap_l(reso/2:-1:2,reso:-1:reso/2+2));
qmap=ifft2(ifftshift(qmap_l))./(delta_x^2);

umap_l(reso/2+1,reso/2+1)=0;
umap_l(1,1)=real(umap_l(1,1));
umap_l(1,reso/2+1)=real(umap_l(1,reso/2+1));
umap_l(reso/2+1,1)=real(umap_l(reso/2+1,1));
umap_l(reso/2+2:reso,1)=conj(umap_l(reso/2:-1:2,1));
umap_l(1,reso/2+2:reso)=conj(umap_l(1,reso/2:-1:2));
umap_l(reso/2+1:reso,reso/2+1:reso)=conj(umap_l(reso/2+1:-1:2,reso/2+1:-1:2));
umap_l(reso/2+2:reso,2:reso/2)=conj(umap_l(reso/2:-1:2,reso:-1:reso/2+2));
umap=ifft2(ifftshift(umap_l))./(delta_x^2);
