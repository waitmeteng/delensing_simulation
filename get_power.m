function [l_array,power,count]=get_power(angu_size,fwhm,map,bin_n,bin1)
% This code is to compute certain power corresponding to certain map
% and certain bin. 
%--------------------------------------------------------------------------
% INPUT:
% .angu_size          : angular size of map (degree)
% .fwhm               : FWHM (arcmin)
% .map                : input map
% .bin_n              : number of bin
% .bin1               : bin
% OUTPUT:
% .l_array            : l for each bin (average of l in each bin)
% .power              : output power
%--------------------------------------------------------------------------
% Wei-Hsiang Teng, NTU, 2010
%--------------------------------------------------------------------------

reso=size(map,1);

area=(angu_size*pi/180)^2;
delta_x=angu_size/reso*pi/180;
fwhm_rad=fwhm/60*pi/180;
dl=360/angu_size; 

l_m=zeros(reso/2+1,reso);
beam_m=zeros(reso/2+1,reso);
for i=1:reso
   for j=1:reso/2+1
       kx=i-(reso/2+1);
       ky=j-(reso/2+1);             
       k=sqrt(kx^2+ky^2);
       ll=k*dl;
       l_m(j,i)=round(ll);
       beam_m(j,i)=exp(ll.*(ll+1)*(fwhm_rad^2)/8/log(2));           
   end
end 

map_l=fftshift(fft2(map))*delta_x^2; 
power_m=l_m.*(l_m+1)/2/pi.*beam_m.*(abs(map_l(1:reso/2+1,:)).^2);

power=zeros(1,bin_n);
count=zeros(1,bin_n);
l_array=zeros(1,bin_n);
for ix=1:reso
    for iy=1:reso/2+1
        kx=ix-(reso/2+1);
        ky=iy-(reso/2+1); 
        k=sqrt(kx^2+ky^2);
        
        ll=k*dl;
        temp=find(bin1>ll);
        if (numel(temp)~=0 && temp(1)~=1)
        bin=temp(1)-1;
        l_array(bin)=l_array(bin)+ll;
        count(bin)=count(bin)+1;
        power(bin)=power(bin)+power_m(iy,ix);
        end
    end
end

power=power./count/area;
l_array=l_array./count;

