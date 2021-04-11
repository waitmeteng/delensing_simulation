function [al_ee] = gene_al_ee(reso,angu_size,lmax,cl_eeg,cl_ee,cl_np)

% generate renormalized factor defined in
% arXiv:astro-ph/0111606 Eq.(11)
% ATTENSION: this code is only for EE estimator!!!
%--------------------------------------------------------------------------
%
%  Wei-Hsiang Teng, NTU Jan 2011
%
%--------------------------------------------------------------------------
delta_x=angu_size/reso*pi/180;
lreso=360/angu_size;
%--------------------------------------------------------------------------
% Define the needed matrix
%--------------------------------------------------------------------------
l_m=zeros(reso/2+1,reso);
lx_m=zeros(reso/2+1,reso);
ly_m=zeros(reso/2+1,reso);
cos2phi_m=zeros(reso/2+1,reso);
cosphi_m=zeros(reso/2+1,reso);
sin2phi_m=zeros(reso/2+1,reso);
cl_eem=zeros(reso/2+1,reso);
cl_eem1=zeros(reso/2+1,reso);
cl_eem2=zeros(reso/2+1,reso);
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
        cosphi_m(j,i)=cos(phi);
        
         if (ll<=lmax)
             l_m(j,i)=round(ll);
            cl_eem(j,i)=cl_eeg(round(ll)+1)^2 ...
                /(cl_ee(round(ll)+1)+cl_np(round(ll)+1));
            cl_eem1(j,i)=1/(cl_ee(round(ll)+1)+cl_np(round(ll)+1)); 
            cl_eem2(j,i)=cl_eeg(round(ll)+1)/(cl_ee(round(ll)+1)+cl_np(round(ll)+1)); 
         end
    end
end
%--------------------------------------------------------------------------
% calculate the needed matrix in Fourier space.
% there are totally 9 matries.
%--------------------------------------------------------------------------
% e11s_l=
% cl_ee^2/(cl_ee+cl_np)*lx^2*sin(2phi)^2
e11s_l=cl_eem.*lx_m.^2.*sin2phi_m.^2;
e11s_l(reso/2+1,reso/2+1)=0;
e11s_l(reso/2+2:reso,:)=e11s_l(reso/2:-1:2,:);
e11s=ifft2(ifftshift(e11s_l))./(delta_x^2);
e11s(reso/2+1,reso/2+1)=0;
e11s(1,1)=real(e11s(1,1));
e11s(1,reso/2+1)=real(e11s(1,reso/2+1));
e11s(reso/2+1,1)=real(e11s(reso/2+1,1));
e11s(reso/2+2:reso,1)=conj(e11s(reso/2:-1:2,1));
e11s(1,reso/2+2:reso)=conj(e11s(1,reso/2:-1:2));
e11s(reso/2+1:reso,reso/2+1:reso)=conj(e11s(reso/2+1:-1:2,reso/2+1:-1:2));
e11s(reso/2+2:reso,2:reso/2)=conj(e11s(reso/2:-1:2,reso:-1:reso/2+2));
% %--------------------------------------------------------------------------
% es_l=
% 1/(cl_ee+cl_np)*sin(2phi)^2
es_l=cl_eem1.*sin2phi_m.^2;
es_l(reso/2+1,reso/2+1)=0;
es_l(reso/2+2:reso,:)=es_l(reso/2:-1:2,:);
es=ifft2(ifftshift(es_l))./(delta_x^2);
es(reso/2+1,reso/2+1)=0;
es(1,1)=real(es(1,1));
es(1,reso/2+1)=real(es(1,reso/2+1));
es(reso/2+1,1)=real(es(reso/2+1,1));
es(reso/2+2:reso,1)=conj(es(reso/2:-1:2,1));
es(1,reso/2+2:reso)=conj(es(1,reso/2:-1:2));
es(reso/2+1:reso,reso/2+1:reso)=conj(es(reso/2+1:-1:2,reso/2+1:-1:2));
es(reso/2+2:reso,2:reso/2)=conj(es(reso/2:-1:2,reso:-1:reso/2+2));
% %--------------------------------------------------------------------------
% e11c_l=
% cl_ee^2/(cl_ee+cl_np)*lx^2*cos(2phi)^2
e11c_l=cl_eem.*lx_m.^2.*cos2phi_m.^2;
e11c_l(reso/2+1,reso/2+1)=0;
e11c_l(reso/2+2:reso,:)=e11c_l(reso/2:-1:2,:);
e11c=ifft2(ifftshift(e11c_l))./(delta_x^2);
e11c(reso/2+1,reso/2+1)=0;
e11c(1,1)=real(e11c(1,1));
e11c(1,reso/2+1)=real(e11c(1,reso/2+1));
e11c(reso/2+1,1)=real(e11c(reso/2+1,1));
e11c(reso/2+2:reso,1)=conj(e11c(reso/2:-1:2,1));
e11c(1,reso/2+2:reso)=conj(e11c(1,reso/2:-1:2));
e11c(reso/2+1:reso,reso/2+1:reso)=conj(e11c(reso/2+1:-1:2,reso/2+1:-1:2));
e11c(reso/2+2:reso,2:reso/2)=conj(e11c(reso/2:-1:2,reso:-1:reso/2+2));
% %--------------------------------------------------------------------------
% ec_l=
% 1/(cl_ee+cl_np)*cos(2phi)^2
ec_l=cl_eem1.*cos2phi_m.^2;
ec_l(reso/2+1,reso/2+1)=0;
ec_l(reso/2+2:reso,:)=ec_l(reso/2:-1:2,:);
ec=ifft2(ifftshift(ec_l))./(delta_x^2);
ec(reso/2+1,reso/2+1)=0;
ec(1,1)=real(ec(1,1));
ec(1,reso/2+1)=real(ec(1,reso/2+1));
ec(reso/2+1,1)=real(ec(reso/2+1,1));
ec(reso/2+2:reso,1)=conj(ec(reso/2:-1:2,1));
ec(1,reso/2+2:reso)=conj(ec(1,reso/2:-1:2));
ec(reso/2+1:reso,reso/2+1:reso)=conj(ec(reso/2+1:-1:2,reso/2+1:-1:2));
ec(reso/2+2:reso,2:reso/2)=conj(ec(reso/2:-1:2,reso:-1:reso/2+2));
% %--------------------------------------------------------------------------
% e11sc_l=
% cl_ee^2/(cl_ee+cl_np)*lx^2*sin(2phi)*cos(2phi)
e11sc_l=cl_eem.*lx_m.^2.*sin2phi_m.*cos2phi_m;
e11sc_l(reso/2+1,reso/2+1)=0;
e11sc_l(reso/2+2:reso,:)=-e11sc_l(reso/2:-1:2,:);
e11sc=ifft2(ifftshift(e11sc_l))./(delta_x^2);
e11sc(reso/2+1,reso/2+1)=0;
e11sc(1,1)=real(e11sc(1,1));
e11sc(1,reso/2+1)=real(e11sc(1,reso/2+1));
e11sc(reso/2+1,1)=real(e11sc(reso/2+1,1));
e11sc(reso/2+2:reso,1)=conj(e11sc(reso/2:-1:2,1));
e11sc(1,reso/2+2:reso)=conj(e11sc(1,reso/2:-1:2));
e11sc(reso/2+1:reso,reso/2+1:reso)=conj(e11sc(reso/2+1:-1:2,reso/2+1:-1:2));
e11sc(reso/2+2:reso,2:reso/2)=conj(e11sc(reso/2:-1:2,reso:-1:reso/2+2));
%--------------------------------------------------------------------------
% esc_l=
% 1/(cl_ee+cl_np)*lx^2*sin(2phi)*cos(2phi)
esc_l=cl_eem1.*sin2phi_m.*cos2phi_m;
esc_l(reso/2+1,reso/2+1)=0;
esc_l(reso/2+2:reso,:)=-esc_l(reso/2:-1:2,:);
esc=ifft2(ifftshift(esc_l))./(delta_x^2);
esc(reso/2+1,reso/2+1)=0;
esc(1,1)=real(esc(1,1));
esc(1,reso/2+1)=real(esc(1,reso/2+1));
esc(reso/2+1,1)=real(esc(reso/2+1,1));
esc(reso/2+2:reso,1)=conj(esc(reso/2:-1:2,1));
esc(1,reso/2+2:reso)=conj(esc(1,reso/2:-1:2));
esc(reso/2+1:reso,reso/2+1:reso)=conj(esc(reso/2+1:-1:2,reso/2+1:-1:2));
esc(reso/2+2:reso,2:reso/2)=conj(esc(reso/2:-1:2,reso:-1:reso/2+2));
% %--------------------------------------------------------------------------
% e1sc_l=
% cl_eeg/(cl_ee+cl_np)*lx*sin(2phi)*cos(2phi)
e1sc_l=cl_eem2.*lx_m.*cos2phi_m.*sin2phi_m;
e1sc_l(reso/2+1,reso/2+1)=0;
e1sc_l(reso/2+2:reso,:)=-e1sc_l(reso/2:-1:2,:);
e1sc=ifft2(ifftshift(e1sc_l))./(delta_x^2);
e1sc(reso/2+1,reso/2+1)=0;
e1sc(1,1)=real(e1sc(1,1));
e1sc(1,reso/2+1)=real(e1sc(1,reso/2+1));
e1sc(reso/2+1,1)=real(e1sc(reso/2+1,1));
e1sc(reso/2+2:reso,1)=conj(e1sc(reso/2:-1:2,1));
e1sc(1,reso/2+2:reso)=conj(e1sc(1,reso/2:-1:2));
e1sc(reso/2+1:reso,reso/2+1:reso)=conj(e1sc(reso/2+1:-1:2,reso/2+1:-1:2));
e1sc(reso/2+2:reso,2:reso/2)=conj(e1sc(reso/2:-1:2,reso:-1:reso/2+2));
%--------------------------------------------------------------------------
% e1c_l=
% cl_eeg/(cl_ee+cl_np)*lx*cos(2phi)^2
e1c_l=cl_eem2.*lx_m.*cos2phi_m.^2;
e1c_l(reso/2+1,reso/2+1)=0;
e1c_l(reso/2+2:reso,:)=e1c_l(reso/2:-1:2,:);
e1c=ifft2(ifftshift(e1c_l))./(delta_x^2);
e1c(reso/2+1,reso/2+1)=0;
e1c(1,1)=real(e1c(1,1));
e1c(1,reso/2+1)=real(e1c(1,reso/2+1));
e1c(reso/2+1,1)=real(e1c(reso/2+1,1));
e1c(reso/2+2:reso,1)=conj(e1c(reso/2:-1:2,1));
e1c(1,reso/2+2:reso)=conj(e1c(1,reso/2:-1:2));
e1c(reso/2+1:reso,reso/2+1:reso)=conj(e1c(reso/2+1:-1:2,reso/2+1:-1:2));
e1c(reso/2+2:reso,2:reso/2)=conj(e1c(reso/2:-1:2,reso:-1:reso/2+2));
%--------------------------------------------------------------------------
% e1s_l=
% cl_eeg/(cl_ee+cl_np)*lx*sin(2phi)^2
e1s_l=cl_eem2.*lx_m.*sin2phi_m.^2;
e1s_l(reso/2+1,reso/2+1)=0;
e1s_l(reso/2+2:reso,:)=e1s_l(reso/2:-1:2,:);
e1s=ifft2(ifftshift(e1s_l))./(delta_x^2);
e1s(reso/2+1,reso/2+1)=0;
e1s(1,1)=real(e1s(1,1));
e1s(1,reso/2+1)=real(e1s(1,reso/2+1));
e1s(reso/2+1,1)=real(e1s(reso/2+1,1));
e1s(reso/2+2:reso,1)=conj(e1s(reso/2:-1:2,1));
e1s(1,reso/2+2:reso)=conj(e1s(1,reso/2:-1:2));
e1s(reso/2+1:reso,reso/2+1:reso)=conj(e1s(reso/2+1:-1:2,reso/2+1:-1:2));
e1s(reso/2+2:reso,2:reso/2)=conj(e1s(reso/2:-1:2,reso:-1:reso/2+2));
%--------------------------------------------------------------------------
%
%
%--------------------------------------------------------------------------
g11=e11s.*es+e11c.*ec+2*e11sc.*esc+e1s.^2+e1c.^2+2*e1sc.^2;
% g22=e22s.*b00c+e22c.*b00s-2*e22sc.*b00sc;
% g12=e12s.*b00c+e12c.*b00s-2*e12sc.*b00sc;
g11_l=real(fftshift(fft2(g11))*delta_x^2); 
% g22_l=real(fftshift(fft2(g22))*delta_x^2);
% g12_l=real(fftshift(fft2(g12))*delta_x^2);
al_eet= g11_l(1:reso/2+1,:).*lx_m.^2;% ...
%       +g22_l(1:reso/2+1,:).*ly_m.^2 ...
%       +g12_l(1:reso/2+1,:).*lx_m.*ly_m;


al_ee=l_m(reso/2+1,reso/2+1:reso).^2./al_eet(reso/2+1,reso/2+1:reso);
al_ee(1)=0;
end



