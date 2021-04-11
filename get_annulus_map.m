function [ map_w ] = get_annulus_map( angu_size, lmin, lmax, map)
% this code is to extract some part of annulus in fourier space .
%--------------------------------------------------------------------------
% INPUT:
% .angu_size                  : angular size of map
% .lmax                       : outer of annulus
% .lmin                       : inner of annulus
% .map                        : input map
% OUTPUT                      : output map
%--------------------------------------------------------------------------
reso=size(map,1);
dl=360/angu_size;
delta_x=angu_size/reso*pi/180;
%--------------------------------------------------------------------------
filter_m=zeros(reso/2+1,reso);

for i=1:reso
   for j=1:reso/2+1
       kx=i-(reso/2+1);
       ky=j-(reso/2+1); 
       k=sqrt(kx^2+ky^2);
       ll=round(k*dl);       
       if (ll<lmax && ll>=lmin)
          filter_m(j,i)=1; 
          
       end       
   end
end 
map_l=fftshift(fft2(map))*delta_x^2; 
map_l(1:reso/2+1,:)=map_l(1:reso/2+1,:).*filter_m;

map_l(reso/2+1,reso/2+1)=0;
map_l(1,1)=real(map_l(1,1));
map_l(1,reso/2+1)=real(map_l(1,reso/2+1));
map_l(reso/2+1,1)=real(map_l(reso/2+1,1));
map_l(reso/2+2:reso,1)=conj(map_l(reso/2:-1:2,1));
map_l(1,reso/2+2:reso)=conj(map_l(1,reso/2:-1:2));
map_l(reso/2+1:reso,reso/2+1:reso)=conj(map_l(reso/2+1:-1:2,reso/2+1:-1:2));
map_l(reso/2+2:reso,2:reso/2)=conj(map_l(reso/2:-1:2,reso:-1:reso/2+2));
map_w=ifft2(ifftshift(map_l))./(delta_x^2);
end