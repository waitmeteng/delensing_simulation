function [tmap,qmap,umap] = gene_tqu(reso,angu_size,seed_n,lmax_cut,l,cl_tt,cl_ee,cl_bb,cl_te)
% this code is to produce T,Q,U maps using flat approximation 
% according to the input power spectra.
% -------------------------------------------------------------------------
% INPUT:
% .reso        : resolution of map
% .angu_size   : angular size of map (degree)
% .seed_n      : random seed of map
% .lmax_cut    : maximum l that not equal to zero
% .l           : minimum l corresponds to input power spectra
% .cl_tt       : TT power
% .cl_ee       : EE power
% .cl_bb       : BB power
% .cl_te       : TE power
%
% OUTPUT:
% .tmap        : T map
% .qmap        : Q map
% .umap        : U map
%--------------------------------------------------------------------------
% Wei-Hsiang Teng, NTU, 2010
%--------------------------------------------------------------------------
dl=360/angu_size;
angu_size_rad=angu_size*pi/180;
area = angu_size_rad*angu_size_rad;
delta_x=angu_size/reso*pi/180;
resocut=2*ceil(lmax_cut/dl);
if (floor(resocut/2)~=resocut/2)
   resocut=resocut+1; 
end
% l_max=ceil(sqrt(2)*dl*reso/2);
% if (l(size(l,2))<l_max)
%     disp('The number of power spectra are not enough!!!'); return
% end
% if (l(1)>dl)
%     disp('The multipole number of first power spectrum should be less than 360/(map of size)!!!'); return
% end

if (2*lmax_cut >= reso/2*dl)
   disp('warning: the resolution is too small that may induce aliasing problem!');  
end

number_l=lmax_cut-l+1;
for ln=1:number_l
    if (cl_tt(ln)*cl_ee(ln)<cl_te(ln)^2 || cl_tt(ln)<0 || cl_ee(ln)<0 || cl_bb(ln)<0)
        disp('The input power spectra are wrong!!!'); return 
    end
    
end
s = RandStream('mcg16807', 'Seed',seed_n);
RandStream.setGlobalStream(s);

if (resocut>reso)
    resocut=reso;
end
l_m=zeros(resocut/2+1,resocut);
cos2phi_m=zeros(resocut/2+1,resocut);
sin2phi_m=zeros(resocut/2+1,resocut);
for i=1:resocut
    for j=1:resocut/2+1
        kx=i-(resocut/2+1);
        ky=j-(resocut/2+1);
        k=sqrt(kx^2+ky^2);
        ll=round(k*dl);
       if (ll<=lmax_cut)
           l_m(j,i)=ll;
       end
        phi=atan2(ky,kx);
        cos2phi_m(j,i)=cos(2*phi); 
        sin2phi_m(j,i)=sin(2*phi); 
    end
end
cl_ttm=cl_tt(l_m-l+1);
cl_eem=cl_ee(l_m-l+1);
cl_bbm=cl_bb(l_m-l+1);
cl_tem=cl_te(l_m-l+1);

number1=(randn(resocut/2+1,resocut)+1i*randn(resocut/2+1,resocut))/sqrt(2);
number2=(randn(resocut/2+1,resocut)+1i*randn(resocut/2+1,resocut))/sqrt(2);
number3=(randn(resocut/2+1,resocut)+1i*randn(resocut/2+1,resocut))/sqrt(2);
tmap_l=zeros(reso,reso);
qmap_l=zeros(reso,reso);
umap_l=zeros(reso,reso);
lbegin=reso/2+1-(resocut/2+1)+1;
lendy=lbegin+resocut/2;
lendx=lbegin+resocut-1;
tmap_l(lbegin:lendy,lbegin:lendx)=number1.*sqrt(area).*(cl_ttm.^0.5);
emap_l=number1.*sqrt(area).*cl_tem./(cl_ttm.^0.5)+...
      number2.*sqrt(area).*((cl_eem-cl_tem.^2./cl_ttm).^0.5);
emap_l(resocut/2+1,resocut/2+1)=0;
for i=1:resocut
    for j=1:resocut/2+1
        kx=i-(resocut/2+1);
        ky=j-(resocut/2+1);
        k=sqrt(kx^2+ky^2);
        ll=round(k*dl);        
        if (ll>lmax_cut)
           emap_l(j,i)=0;
        end
    end
end
bmap_l=number3.*sqrt(area).*(cl_bbm.^0.5);

qmap_l(lbegin:lendy,lbegin:lendx)=emap_l.*cos2phi_m-bmap_l.*sin2phi_m;
umap_l(lbegin:lendy,lbegin:lendx)=emap_l.*sin2phi_m+bmap_l.*cos2phi_m;

tmap_l(reso/2+1,reso/2+1)=0;
tmap_l(1,1)=real(tmap_l(1,1));
tmap_l(1,reso/2+1)=real(tmap_l(1,reso/2+1));
tmap_l(reso/2+1,1)=real(tmap_l(reso/2+1,1));
tmap_l(reso/2+2:reso,1)=conj(tmap_l(reso/2:-1:2,1));
tmap_l(1,reso/2+2:reso)=conj(tmap_l(1,reso/2:-1:2));
tmap_l(reso/2+1:reso,reso/2+1:reso)=conj(tmap_l(reso/2+1:-1:2,reso/2+1:-1:2));
tmap_l(reso/2+2:reso,2:reso/2)=conj(tmap_l(reso/2:-1:2,reso:-1:reso/2+2));
tmap=ifft2(ifftshift(tmap_l))./(delta_x^2);

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



