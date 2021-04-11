function [ map_w ] = wiener_filter( angu_size, lmax, map, cl_s, cl_n)
% this code is to produce Wiener filtered map for visual reason.
%--------------------------------------------------------------------------
% INPUT:
% .angu_size                  : angular size of map
% .lmax                       : LMAX
% .map                        : input map
% .cl_s                       : power of signal
% .cl_n                       : power of noise
% OUTPUT                      : Wiener filtered map
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
       if (ll<=lmax)
          filter_m(j,i)=cl_s(ll+1)/(cl_s(ll+1)+cl_n(ll+1)); 
          
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

