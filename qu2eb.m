function [ emap, bmap ] = qu2eb(angu_size, lmax, qmap, umap, fwhm)
% this code is to tranfer Q, U maps into E, B maps
% -------------------------------------------------------------------------
% INPUT:
% .angu_size                : angular size of Q, U maps
% .qmap                     : input Q map
% .umap                     : input U map
% .fwhm                     : if fwhm=0, no deconvolution is done.
% OUTPUT:
% .emap                     : output E map
% .bmap                     : output B map
%--------------------------------------------------------------------------
reso=size(qmap,1);
if (size(umap,1)~=reso)
   disp('The size of input maps are not the same!'); return 
end
delta_x=angu_size/reso*pi/180;
fwhm_rad=fwhm/60*pi/180;
dl=360/angu_size;
%--------------------------------------------------------------------------
cos2phi_m=zeros(reso/2+1,reso);
sin2phi_m=zeros(reso/2+1,reso);
beam_m=zeros(reso/2+1,reso);
for i=1:reso
   for j=1:reso/2+1
       kx=i-(reso/2+1);
       ky=j-(reso/2+1);      
       k=sqrt(kx^2+ky^2);
       ll=k*dl;
       phi=atan2(ky,kx);
       cos2phi_m(j,i)=cos(2*phi); 
       sin2phi_m(j,i)=sin(2*phi);
       if (ll<=lmax)
       beam_m(j,i)=exp(ll.*(ll+1)*(fwhm_rad^2)/16/log(2));
       end
   end
end 
emap_l=zeros(reso,reso);
bmap_l=zeros(reso,reso);
qmap_l=fftshift(fft2(qmap))*delta_x^2; 
umap_l=fftshift(fft2(umap))*delta_x^2;
bmap_l(1:reso/2+1,:)= ...
    (-qmap_l(1:reso/2+1,:).*sin2phi_m+umap_l(1:reso/2+1,:).*cos2phi_m).*beam_m;
emap_l(1:reso/2+1,:)= ...
    (qmap_l(1:reso/2+1,:).*cos2phi_m+umap_l(1:reso/2+1,:).*sin2phi_m).*beam_m;

emap_l(reso/2+1,reso/2+1)=0;
emap_l(1,1)=real(emap_l(1,1));
emap_l(1,reso/2+1)=real(emap_l(1,reso/2+1));
emap_l(reso/2+1,1)=real(emap_l(reso/2+1,1));
emap_l(reso/2+2:reso,1)=conj(emap_l(reso/2:-1:2,1));
emap_l(1,reso/2+2:reso)=conj(emap_l(1,reso/2:-1:2));
emap_l(reso/2+1:reso,reso/2+1:reso)=conj(emap_l(reso/2+1:-1:2,reso/2+1:-1:2));
emap_l(reso/2+2:reso,2:reso/2)=conj(emap_l(reso/2:-1:2,reso:-1:reso/2+2));
emap=ifft2(ifftshift(emap_l))./(delta_x^2);

bmap_l(reso/2+1,reso/2+1)=0;
bmap_l(1,1)=real(bmap_l(1,1));
bmap_l(1,reso/2+1)=real(bmap_l(1,reso/2+1));
bmap_l(reso/2+1,1)=real(bmap_l(reso/2+1,1));
bmap_l(reso/2+2:reso,1)=conj(bmap_l(reso/2:-1:2,1));
bmap_l(1,reso/2+2:reso)=conj(bmap_l(1,reso/2:-1:2));
bmap_l(reso/2+1:reso,reso/2+1:reso)=conj(bmap_l(reso/2+1:-1:2,reso/2+1:-1:2));
bmap_l(reso/2+2:reso,2:reso/2)=conj(bmap_l(reso/2:-1:2,reso:-1:reso/2+2));
bmap=ifft2(ifftshift(bmap_l))./(delta_x^2);
