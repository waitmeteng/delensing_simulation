function [defmap]=ebee_estimatorclnn(angu_size,fwhm,l_max,l_min,l_min0,al,qmap,umap,l,cl_eeg,cl_ee,cl_bbg,cl_bb,cl_dd,cl_enp,cl_bnp,wiener,cl_nn)
% Hi! Nice to meet you.
% This function only supports the periodic boundary condition.

reso=size(qmap,1);
fwhm_rad=fwhm/60*pi/180;
dl=360/angu_size;
if (size(umap,1)~=reso)
   disp('The size of input maps are not the same!'); return 
end
delta_x=angu_size/reso*pi/180;
qmap_l=fftshift(fft2(qmap))*delta_x^2;   
umap_l=fftshift(fft2(umap))*delta_x^2; 

l_m=zeros(reso/2+1,reso);
l_m1=zeros(reso/2+1,reso);
lx_m=zeros(reso/2+1,reso);
ly_m=zeros(reso/2+1,reso);
cl_eem=zeros(reso/2+1,reso);
cl_bbm=zeros(reso/2+1,reso);
cl_bbgm=zeros(reso/2+1,reso);
cl_eegm=zeros(reso/2+1,reso);
cl_ddm=zeros(reso/2+1,reso);
cl_ddm1=zeros(reso/2+1,reso);
cos2phi_m=zeros(reso/2+1,reso);
sin2phi_m=zeros(reso/2+1,reso);
beam_m=zeros(reso/2+1,reso);
for i=1:reso
    for j=1:reso/2+1
        kx=i-(reso/2+1);
        ky=j-(reso/2+1);
        k=sqrt(kx^2+ky^2);
        ll=k*dl;
        l_m(j,i)=round(ll);
        lx_m(j,i)=kx*dl;
        ly_m(j,i)=ky*dl;
        phi=atan2(ky,kx);
        cos2phi_m(j,i)=cos(2*phi); 
        sin2phi_m(j,i)=sin(2*phi);
        if (ll<=l_max)
        l_m1(j,i)=round(ll);    
        
        cl_eem(j,i)=1/(cl_ee(round(ll)-l+1)+cl_enp(round(ll)-l+1));
        cl_bbm(j,i)=1/(cl_bb(round(ll)-l+1)+cl_bnp(round(ll)-l+1));
        cl_eegm(j,i)=cl_eeg(round(ll)-l+1);
        cl_bbgm(j,i)=cl_bbg(round(ll)-l+1);
        cl_ddm(j,i)=cl_dd(round(ll)-l+1)/(cl_dd(round(ll)-l+1)+cl_nn(round(ll)-l+1));
        cl_ddm1(j,i)=sqrt(cl_dd(round(ll)-l+1)/(cl_dd(round(ll)-l+1)+cl_nn(round(ll)-l+1)));
        end
        if (ll<=l_max && ll>=l_min )
           beam_m(j,i)=exp(ll.*(ll+1)*(fwhm_rad^2)/16/log(2)); 
        end
        if (ll<l_min0)
           beam_m(j,i)=exp(ll.*(ll+1)*(fwhm_rad^2)/16/log(2)); 
        end
    end
end

cl_eem(reso/2+1,reso/2+1)=0;
cl_bbm(reso/2+1,reso/2+1)=0;
cl_eegm(reso/2+1,reso/2+1)=0;
cl_ddm(reso/2+1,reso/2+1)=0;
al_m=al(l_m1-l+1);

% cl_nnm=cl_nn(l_m-l(1)+1);

emap_l=(qmap_l(1:reso/2+1,:).*cos2phi_m+umap_l(1:reso/2+1,:).*sin2phi_m).*beam_m.*cl_eegm.*cl_eem;
bmap_l=(-qmap_l(1:reso/2+1,:).*sin2phi_m+umap_l(1:reso/2+1,:).*cos2phi_m).*beam_m.*cl_bbgm.*cl_bbm;
emap_l1=(qmap_l(1:reso/2+1,:).*cos2phi_m+umap_l(1:reso/2+1,:).*sin2phi_m).*beam_m.*cl_eem;
bmap_l1=(-qmap_l(1:reso/2+1,:).*sin2phi_m+umap_l(1:reso/2+1,:).*cos2phi_m).*beam_m.*cl_bbm;
%-----------------------start to generate the deflection field-------------
q1_l=zeros(reso,reso);
q1_l(1:reso/2+1,:)=1i*(lx_m).*(emap_l.*cos2phi_m-bmap_l.*sin2phi_m);
q1_l(reso/2+1,reso/2+1)=0;
q1_l(1,1)=real(q1_l(1,1));
q1_l(1,reso/2+1)=real(q1_l(1,reso/2+1));
q1_l(reso/2+1,1)=real(q1_l(reso/2+1,1));
q1_l(reso/2+2:reso,1)=conj(q1_l(reso/2:-1:2,1));
q1_l(1,reso/2+2:reso)=conj(q1_l(1,reso/2:-1:2));
q1_l(reso/2+1:reso,reso/2+1:reso)=conj(q1_l(reso/2+1:-1:2,reso/2+1:-1:2));
q1_l(reso/2+2:reso,2:reso/2)=conj(q1_l(reso/2:-1:2,reso:-1:reso/2+2));
q1=ifft2(ifftshift(q1_l))./(delta_x^2);

q2_l=zeros(reso,reso);
q2_l(1:reso/2+1,:)=1i*(ly_m).*(emap_l.*cos2phi_m-bmap_l.*sin2phi_m);
q2_l(reso/2+1,reso/2+1)=0;
q2_l(1,1)=real(q2_l(1,1));
q2_l(1,reso/2+1)=real(q2_l(1,reso/2+1));
q2_l(reso/2+1,1)=real(q2_l(reso/2+1,1));
q2_l(reso/2+2:reso,1)=conj(q2_l(reso/2:-1:2,1));
q2_l(1,reso/2+2:reso)=conj(q2_l(1,reso/2:-1:2));
q2_l(reso/2+1:reso,reso/2+1:reso)=conj(q2_l(reso/2+1:-1:2,reso/2+1:-1:2));
q2_l(reso/2+2:reso,2:reso/2)=conj(q2_l(reso/2:-1:2,reso:-1:reso/2+2));
q2=ifft2(ifftshift(q2_l))./(delta_x^2);

u1_l=zeros(reso,reso);
u1_l(1:reso/2+1,:)=1i*(lx_m).*(emap_l.*sin2phi_m+bmap_l.*cos2phi_m);
u1_l(reso/2+1,reso/2+1)=0;
u1_l(1,1)=real(u1_l(1,1));
u1_l(1,reso/2+1)=real(u1_l(1,reso/2+1));
u1_l(reso/2+1,1)=real(u1_l(reso/2+1,1));
u1_l(reso/2+2:reso,1)=conj(u1_l(reso/2:-1:2,1));
u1_l(1,reso/2+2:reso)=conj(u1_l(1,reso/2:-1:2));
u1_l(reso/2+1:reso,reso/2+1:reso)=conj(u1_l(reso/2+1:-1:2,reso/2+1:-1:2));
u1_l(reso/2+2:reso,2:reso/2)=conj(u1_l(reso/2:-1:2,reso:-1:reso/2+2));
u1=ifft2(ifftshift(u1_l))./(delta_x^2);

u2_l=zeros(reso,reso);
u2_l(1:reso/2+1,:)=1i*(ly_m).*(emap_l.*sin2phi_m+bmap_l.*cos2phi_m);
u2_l(reso/2+1,reso/2+1)=0;
u2_l(1,1)=real(u2_l(1,1));
u2_l(1,reso/2+1)=real(u2_l(1,reso/2+1));
u2_l(reso/2+1,1)=real(u2_l(reso/2+1,1));
u2_l(reso/2+2:reso,1)=conj(u2_l(reso/2:-1:2,1));
u2_l(1,reso/2+2:reso)=conj(u2_l(1,reso/2:-1:2));
u2_l(reso/2+1:reso,reso/2+1:reso)=conj(u2_l(reso/2+1:-1:2,reso/2+1:-1:2));
u2_l(reso/2+2:reso,2:reso/2)=conj(u2_l(reso/2:-1:2,reso:-1:reso/2+2));
u2=ifft2(ifftshift(u2_l))./(delta_x^2);

q_l=zeros(reso,reso);
q_l(1:reso/2+1,:)=(emap_l1.*cos2phi_m-bmap_l1.*sin2phi_m);
q_l(reso/2+1,reso/2+1)=0;
q_l(1,1)=real(q_l(1,1));
q_l(1,reso/2+1)=real(q_l(1,reso/2+1));
q_l(reso/2+1,1)=real(q_l(reso/2+1,1));
q_l(reso/2+2:reso,1)=conj(q_l(reso/2:-1:2,1));
q_l(1,reso/2+2:reso)=conj(q_l(1,reso/2:-1:2));
q_l(reso/2+1:reso,reso/2+1:reso)=conj(q_l(reso/2+1:-1:2,reso/2+1:-1:2));
q_l(reso/2+2:reso,2:reso/2)=conj(q_l(reso/2:-1:2,reso:-1:reso/2+2));
q=ifft2(ifftshift(q_l))./(delta_x^2);

u_l=zeros(reso,reso);
u_l(1:reso/2+1,:)=(emap_l1.*sin2phi_m+bmap_l1.*cos2phi_m);
u_l(reso/2+1,reso/2+1)=0;
u_l(1,1)=real(u_l(1,1));
u_l(1,reso/2+1)=real(u_l(1,reso/2+1));
u_l(reso/2+1,1)=real(u_l(reso/2+1,1));
u_l(reso/2+2:reso,1)=conj(u_l(reso/2:-1:2,1));
u_l(1,reso/2+2:reso)=conj(u_l(1,reso/2:-1:2));
u_l(reso/2+1:reso,reso/2+1:reso)=conj(u_l(reso/2+1:-1:2,reso/2+1:-1:2));
u_l(reso/2+2:reso,2:reso/2)=conj(u_l(reso/2:-1:2,reso:-1:reso/2+2));
u=ifft2(ifftshift(u_l))./(delta_x^2);
%--------------------------------------------------------------------------
g1=(q1.*q+u1.*u);
g2=(q2.*q+u2.*u);
g1_l=fftshift(fft2(g1))*delta_x^2; 
g2_l=fftshift(fft2(g2))*delta_x^2;
defmap_l=zeros(reso,reso);
if (strcmp(wiener,'wiener')==1)
     defmap_l(1:reso/2+1,:)=1i*al_m./l_m.*(lx_m.*g1_l(1:reso/2+1,:) ...
         +ly_m.*g2_l(1:reso/2+1,:)).*cl_ddm;
elseif (strcmp(wiener,'wiener_modified')==1)
     defmap_l(1:reso/2+1,:)=1i*al_m./l_m.*(lx_m.*g1_l(1:reso/2+1,:) ...
         +ly_m.*g2_l(1:reso/2+1,:)).*cl_ddm1;
else
     defmap_l(1:reso/2+1,:)=1i*al_m./l_m.*(lx_m.*g1_l(1:reso/2+1,:) ...
         +ly_m.*g2_l(1:reso/2+1,:));
end
defmap_l(reso/2+1,reso/2+1)=0;
defmap_l(1,1)=real(defmap_l(1,1));
defmap_l(1,reso/2+1)=real(defmap_l(1,reso/2+1));
defmap_l(reso/2+1,1)=real(defmap_l(reso/2+1,1));
defmap_l(reso/2+2:reso,1)=conj(defmap_l(reso/2:-1:2,1));
defmap_l(1,reso/2+2:reso)=conj(defmap_l(1,reso/2:-1:2));
defmap_l(reso/2+1:reso,reso/2+1:reso)=conj(defmap_l(reso/2+1:-1:2,reso/2+1:-1:2));
defmap_l(reso/2+2:reso,2:reso/2)=conj(defmap_l(reso/2:-1:2,reso:-1:reso/2+2));
defmap=ifft2(ifftshift(defmap_l))./(delta_x^2);
