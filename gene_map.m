function [ map ] = gene_map( reso, angu_size, seed, lmax, l, cl )
% this code is to generate one map for certain power 

dl=360/angu_size;
angu_size_rad=angu_size*pi/180;
area = angu_size_rad*angu_size_rad;
delta_x=angu_size/reso*pi/180;
resocut=2*ceil(lmax/dl);
if (floor(resocut/2)~=resocut/2)
   resocut=resocut+1; 
end

if (2*lmax >= reso/2*dl)
   disp('warning: the resolution is too small that may induce aliasing problem!');  
end

number_l=lmax-l+1;
for ln=1:number_l
    if (cl(ln)<0)
        disp('The input power spectra are wrong!!!'); return 
    end
    
end
s = RandStream('mcg16807', 'Seed',seed);
RandStream.setGlobalStream(s);

if (resocut>reso)
    resocut=reso;
end
l_m=zeros(resocut/2+1,resocut);
for i=1:resocut
    for j=1:resocut/2+1
        kx=i-(resocut/2+1);
        ky=j-(resocut/2+1);
        k=sqrt(kx^2+ky^2);
       ll=round(k*dl);
       if (ll<=lmax)
           l_m(j,i)=ll;
       end
       
    end
end
cl_m=cl(l_m-l+1);

number1=(randn(resocut/2+1,resocut)+1i*randn(resocut/2+1,resocut))/sqrt(2);

map_l=zeros(reso,reso);
lbegin=reso/2+1-(resocut/2+1)+1;
lendy=lbegin+resocut/2;
lendx=lbegin+resocut-1;
map_l(lbegin:lendy,lbegin:lendx)=number1.*sqrt(area).*(cl_m.^0.5);

map_l(reso/2+1,reso/2+1)=0;
map_l(1,1)=real(map_l(1,1));
map_l(1,reso/2+1)=real(map_l(1,reso/2+1));
map_l(reso/2+1,1)=real(map_l(reso/2+1,1));
map_l(reso/2+2:reso,1)=conj(map_l(reso/2:-1:2,1));
map_l(1,reso/2+2:reso)=conj(map_l(1,reso/2:-1:2));
map_l(reso/2+1:reso,reso/2+1:reso)=conj(map_l(reso/2+1:-1:2,reso/2+1:-1:2));
map_l(reso/2+2:reso,2:reso/2)=conj(map_l(reso/2:-1:2,reso:-1:reso/2+2));
map=ifft2(ifftshift(map_l))./(delta_x^2);


end

