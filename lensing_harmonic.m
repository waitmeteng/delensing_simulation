function [bmap_lens] = lensing_harmonic(angu_size, lmax, bmap, emap, defmap)
% this code is to generate delensed B map according to the delensing
% estimator suggested by Smith et al. 
% arXiv:1010.0048 Eq. (9)
% Here we use flat approximation.
%--------------------------------------------------------------------------
% INPUT:
% .angu_size                     : angular size of map
% .lmax                          : LMAX
% .bmap                          : input lensed B map
% .emap                          : input Wiener filtered lensed E map 
% .defmap                        : input reconstructed deflection map
%                                  (should Wiener filtered)
% OUTPUT:
% .bmap_delens                   : delensed B map
%--------------------------------------------------------------------------
reso=size(bmap,1);
if (size(emap,1)~=reso || size(defmap,1)~=reso)
   disp('The size of input maps are not the same!'); return 
end
delta_x=angu_size/reso*pi/180;
lreso=360/angu_size;
%--------------------------------------------------------------------------
emap_l=fftshift(fft2(emap))*delta_x^2;   
bmap_l=fftshift(fft2(bmap))*delta_x^2; 
defmap_l=fftshift(fft2(defmap))*delta_x^2; 

lx_m=zeros(reso/2+1,reso);
ly_m=zeros(reso/2+1,reso);
phi_l=zeros(reso/2+1,reso);
cos2phi_m=zeros(reso/2+1,reso);
sin2phi_m=zeros(reso/2+1,reso);
for i=1:reso
    for j=1:reso/2+1
        kx=i-(reso/2+1);
        ky=j-(reso/2+1);
        k=sqrt(kx^2+ky^2);
        ll=k*lreso;
        lx_m(j,i)=kx*lreso;
        ly_m(j,i)=ky*lreso;
        phi=atan2(ky,kx);
        cos2phi_m(j,i)=cos(2*phi);
        sin2phi_m(j,i)=sin(2*phi);
        
         if (ll<=lmax)
            phi_l(j,i)=defmap_l(j,i)/ll;           
         else 
            emap_l(j,i)=0;
            bmap_l(j,i)=0;
         end
    end
end
emap_l=emap_l(1:reso/2+1,:);
bmap_l=bmap_l(1:reso/2+1,:);
%--------------------------------------------------------------------------
% e111_l
%       =1i*emap_l*sin(2phi)*lx 
e111_l=1i*emap_l.*sin2phi_m.*lx_m;
e111_l(reso/2+1,reso/2+1)=0;
e111_l(1,1)=real(e111_l(1,1));
e111_l(1,reso/2+1)=real(e111_l(1,reso/2+1));
e111_l(reso/2+1,1)=real(e111_l(reso/2+1,1));
e111_l(reso/2+2:reso,1)=conj(e111_l(reso/2:-1:2,1));
e111_l(1,reso/2+2:reso)=conj(e111_l(1,reso/2:-1:2));
e111_l(reso/2+1:reso,reso/2+1:reso)=conj(e111_l(reso/2+1:-1:2,reso/2+1:-1:2));
e111_l(reso/2+2:reso,2:reso/2)=conj(e111_l(reso/2:-1:2,reso:-1:reso/2+2));
e111=ifft2(ifftshift(e111_l))./(delta_x^2);
%--------------------------------------------------------------------------
% e222_l
%       =1i*emap_l*sin(2phi)*ly
e222_l=1i*emap_l.*sin2phi_m.*ly_m;
e222_l(reso/2+1,reso/2+1)=0;
e222_l(1,1)=real(e222_l(1,1));
e222_l(1,reso/2+1)=real(e222_l(1,reso/2+1));
e222_l(reso/2+1,1)=real(e222_l(reso/2+1,1));
e222_l(reso/2+2:reso,1)=conj(e222_l(reso/2:-1:2,1));
e222_l(1,reso/2+2:reso)=conj(e222_l(1,reso/2:-1:2));
e222_l(reso/2+1:reso,reso/2+1:reso)=conj(e222_l(reso/2+1:-1:2,reso/2+1:-1:2));
e222_l(reso/2+2:reso,2:reso/2)=conj(e222_l(reso/2:-1:2,reso:-1:reso/2+2));
e222=ifft2(ifftshift(e222_l))./(delta_x^2);
%--------------------------------------------------------------------------
% e121_l
%       =1i*emap_l*cos(2phi)*lx
e121_l=1i*emap_l.*cos2phi_m.*lx_m;
e121_l(reso/2+1,reso/2+1)=0;
e121_l(1,1)=real(e121_l(1,1));
e121_l(1,reso/2+1)=real(e121_l(1,reso/2+1));
e121_l(reso/2+1,1)=real(e121_l(reso/2+1,1));
e121_l(reso/2+2:reso,1)=conj(e121_l(reso/2:-1:2,1));
e121_l(1,reso/2+2:reso)=conj(e121_l(1,reso/2:-1:2));
e121_l(reso/2+1:reso,reso/2+1:reso)=conj(e121_l(reso/2+1:-1:2,reso/2+1:-1:2));
e121_l(reso/2+2:reso,2:reso/2)=conj(e121_l(reso/2:-1:2,reso:-1:reso/2+2));
e121=ifft2(ifftshift(e121_l))./(delta_x^2);
%--------------------------------------------------------------------------
% e122_l
%       =1i*emap_l*cos(2phi)*ly
e122_l=1i*emap_l.*cos2phi_m.*ly_m;
e122_l(reso/2+1,reso/2+1)=0;
e122_l(1,1)=real(e122_l(1,1));
e122_l(1,reso/2+1)=real(e122_l(1,reso/2+1));
e122_l(reso/2+1,1)=real(e122_l(reso/2+1,1));
e122_l(reso/2+2:reso,1)=conj(e122_l(reso/2:-1:2,1));
e122_l(1,reso/2+2:reso)=conj(e122_l(1,reso/2:-1:2));
e122_l(reso/2+1:reso,reso/2+1:reso)=conj(e122_l(reso/2+1:-1:2,reso/2+1:-1:2));
e122_l(reso/2+2:reso,2:reso/2)=conj(e122_l(reso/2:-1:2,reso:-1:reso/2+2));
e122=ifft2(ifftshift(e122_l))./(delta_x^2);
%--------------------------------------------------------------------------
% b12_l
%      =1i*phi_l*lx
b12_l=1i*phi_l.*lx_m;
b12_l(reso/2+1,reso/2+1)=0;
b12_l(1,1)=real(b12_l(1,1));
b12_l(1,reso/2+1)=real(b12_l(1,reso/2+1));
b12_l(reso/2+1,1)=real(b12_l(reso/2+1,1));
b12_l(reso/2+2:reso,1)=conj(b12_l(reso/2:-1:2,1));
b12_l(1,reso/2+2:reso)=conj(b12_l(1,reso/2:-1:2));
b12_l(reso/2+1:reso,reso/2+1:reso)=conj(b12_l(reso/2+1:-1:2,reso/2+1:-1:2));
b12_l(reso/2+2:reso,2:reso/2)=conj(b12_l(reso/2:-1:2,reso:-1:reso/2+2));
b12=ifft2(ifftshift(b12_l))./(delta_x^2);
%--------------------------------------------------------------------------
% b22_l
%      =1i*phi_l*ly
b22_l=1i*phi_l.*ly_m;
b22_l(reso/2+1,reso/2+1)=0;
b22_l(1,1)=real(b22_l(1,1));
b22_l(1,reso/2+1)=real(b22_l(1,reso/2+1));
b22_l(reso/2+1,1)=real(b22_l(reso/2+1,1));
b22_l(reso/2+2:reso,1)=conj(b22_l(reso/2:-1:2,1));
b22_l(1,reso/2+2:reso)=conj(b22_l(1,reso/2:-1:2));
b22_l(reso/2+1:reso,reso/2+1:reso)=conj(b22_l(reso/2+1:-1:2,reso/2+1:-1:2));
b22_l(reso/2+2:reso,2:reso/2)=conj(b22_l(reso/2:-1:2,reso:-1:reso/2+2));
b22=ifft2(ifftshift(b22_l))./(delta_x^2);
%--------------------------------------------------------------------------
bmap_lens_l=zeros(reso,reso);
g1=e111.*b12+e222.*b22;
g2=e121.*b12+e122.*b22;
g1_l=fftshift(fft2(g1))*delta_x^2; 
g2_l=fftshift(fft2(g2))*delta_x^2;
bmap_diff_l= ...
    g1_l(1:reso/2+1,:).*cos2phi_m-g2_l(1:reso/2+1,:).*sin2phi_m;
bmap_lens_l(1:reso/2+1,:)=bmap_l+bmap_diff_l;
bmap_lens_l(reso/2+1,reso/2+1)=0;
bmap_lens_l(1,1)=real(bmap_lens_l(1,1));
bmap_lens_l(1,reso/2+1)=real(bmap_lens_l(1,reso/2+1));
bmap_lens_l(reso/2+1,1)=real(bmap_lens_l(reso/2+1,1));
bmap_lens_l(reso/2+2:reso,1)=conj(bmap_lens_l(reso/2:-1:2,1));
bmap_lens_l(1,reso/2+2:reso)=conj(bmap_lens_l(1,reso/2:-1:2));
bmap_lens_l(reso/2+1:reso,reso/2+1:reso)=conj(bmap_lens_l(reso/2+1:-1:2,reso/2+1:-1:2));
bmap_lens_l(reso/2+2:reso,2:reso/2)=conj(bmap_lens_l(reso/2:-1:2,reso:-1:reso/2+2));
bmap_lens=ifft2(ifftshift(bmap_lens_l))./(delta_x^2);
