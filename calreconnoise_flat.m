function [cl_nn, l_out] = calreconnoise_flat(lmax,lreso,cl_ee,cl_bb,cl_enp,cl_bnp)

% generate reconstruction noise defined in
% arXiv:astro-ph/0111606 Eq.(11)
% ATTENSION: this code is only for EB estimator!!!
%            and the output is for lensing potential,
%            not deflection field!!!
%--------------------------------------------------------------------------
%
%  Wei-Hsiang Teng, NTU Dec 2010
%
%--------------------------------------------------------------------------
%
% Define the needed parameters:
% . the size of map (degree)  : ANGU_SIZE
% . the resolution            : RESO
%--------------------------------------------------------------------------
ANGU_SIZE=360/lreso;
reso=2*floor(lmax/lreso);
delta_x=ANGU_SIZE/reso*pi/180;
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
cl_bbm=zeros(reso/2+1,reso);
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
        l_m(j,i)=round(ll);
         if (ll<=lmax)
            cl_eem(j,i)=cl_ee(round(ll)+1)^2 ...
                /(cl_ee(round(ll)+1)+cl_enp(round(ll)+1));
            cl_bbm(j,i)=1/(cl_bb(round(ll)+1)+cl_bnp(round(ll)+1));    
         end
    end
end
%--------------------------------------------------------------------------
% calculate the needed matrix in Fourier space.
% there are totally 12 matries.
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
%--------------------------------------------------------------------------
% % e22s_l=
% % cl_ee^2/(cl_ee+cl_np)*ly^2*sin(2phi)^2
% e22s_l=cl_eem.*ly_m.^2.*sin2phi_m.^2;
% e22s_l(reso/2+1,reso/2+1)=0;
% e22s_l(reso/2+2:reso,:)=e22s_l(reso/2:-1:2,:);
% e22s=ifft2(ifftshift(e22s_l))./(delta_x^2);
% e22s(reso/2+1,reso/2+1)=0;
% e22s(1,1)=real(e22s(1,1));
% e22s(1,reso/2+1)=real(e22s(1,reso/2+1));
% e22s(reso/2+1,1)=real(e22s(reso/2+1,1));
% e22s(reso/2+2:reso,1)=conj(e22s(reso/2:-1:2,1));
% e22s(1,reso/2+2:reso)=conj(e22s(1,reso/2:-1:2));
% e22s(reso/2+1:reso,reso/2+1:reso)=conj(e22s(reso/2+1:-1:2,reso/2+1:-1:2));
% e22s(reso/2+2:reso,2:reso/2)=conj(e22s(reso/2:-1:2,reso:-1:reso/2+2));
%--------------------------------------------------------------------------
% % e12s_l=
% % cl_ee^2/(cl_ee+cl_np)*lx*ly*sin(2phi)^2
% e12s_l=cl_eem.*lx_m.*ly_m.*sin2phi_m.^2;
% e12s_l(reso/2+1,reso/2+1)=0;
% e12s_l(reso/2+2:reso,:)=-e12s_l(reso/2:-1:2,:);
% e12s=ifft2(ifftshift(e12s_l))./(delta_x^2);
% e12s(reso/2+1,reso/2+1)=0;
% e12s(1,1)=real(e12s(1,1));
% e12s(1,reso/2+1)=real(e12s(1,reso/2+1));
% e12s(reso/2+1,1)=real(e12s(reso/2+1,1));
% e12s(reso/2+2:reso,1)=conj(e12s(reso/2:-1:2,1));
% e12s(1,reso/2+2:reso)=conj(e12s(1,reso/2:-1:2));
% e12s(reso/2+1:reso,reso/2+1:reso)=conj(e12s(reso/2+1:-1:2,reso/2+1:-1:2));
% e12s(reso/2+2:reso,2:reso/2)=conj(e12s(reso/2:-1:2,reso:-1:reso/2+2));
%--------------------------------------------------------------------------
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
%--------------------------------------------------------------------------
% % e22c_l=
% % cl_ee^2/(cl_ee+cl_np)*ly^2*cos(2phi)^2
% e22c_l=cl_eem.*ly_m.^2.*cos2phi_m.^2;
% e22c_l(reso/2+1,reso/2+1)=0;
% e22c_l(reso/2+2:reso,:)=e22c_l(reso/2:-1:2,:);
% e22c=ifft2(ifftshift(e22c_l))./(delta_x^2);
% e22c(reso/2+1,reso/2+1)=0;
% e22c(1,1)=real(e22c(1,1));
% e22c(1,reso/2+1)=real(e22c(1,reso/2+1));
% e22c(reso/2+1,1)=real(e22c(reso/2+1,1));
% e22c(reso/2+2:reso,1)=conj(e22c(reso/2:-1:2,1));
% e22c(1,reso/2+2:reso)=conj(e22c(1,reso/2:-1:2));
% e22c(reso/2+1:reso,reso/2+1:reso)=conj(e22c(reso/2+1:-1:2,reso/2+1:-1:2));
% e22c(reso/2+2:reso,2:reso/2)=conj(e22c(reso/2:-1:2,reso:-1:reso/2+2));
% %--------------------------------------------------------------------------
% % e12c_l=
% % cl_ee^2/(cl_ee+cl_np)*lx*ly*cos(2phi)^2
% e12c_l=cl_eem.*ly_m.*lx_m.*cos2phi_m.^2;
% e12c_l(reso/2+1,reso/2+1)=0;
% e12c_l(reso/2+2:reso,:)=-e12c_l(reso/2:-1:2,:);
% e12c=ifft2(ifftshift(e12c_l))./(delta_x^2);
% e12c(reso/2+1,reso/2+1)=0;
% e12c(1,1)=real(e12c(1,1));
% e12c(1,reso/2+1)=real(e12c(1,reso/2+1));
% e12c(reso/2+1,1)=real(e12c(reso/2+1,1));
% e12c(reso/2+2:reso,1)=conj(e12c(reso/2:-1:2,1));
% e12c(1,reso/2+2:reso)=conj(e12c(1,reso/2:-1:2));
% e12c(reso/2+1:reso,reso/2+1:reso)=conj(e12c(reso/2+1:-1:2,reso/2+1:-1:2));
% e12c(reso/2+2:reso,2:reso/2)=conj(e12c(reso/2:-1:2,reso:-1:reso/2+2));
%--------------------------------------------------------------------------
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
% % e22sc_l=
% % cl_ee^2/(cl_ee+cl_np)*ly^2*sin(2phi)*cos(2phi)
% e22sc_l=cl_eem.*ly_m.^2.*sin2phi_m.*cos2phi_m;
% e22sc_l(reso/2+1,reso/2+1)=0;
% e22sc_l(reso/2+2:reso,:)=-e22sc_l(reso/2:-1:2,:);
% e22sc=ifft2(ifftshift(e22sc_l))./(delta_x^2);
% e22sc(reso/2+1,reso/2+1)=0;
% e22sc(1,1)=real(e22sc(1,1));
% e22sc(1,reso/2+1)=real(e22sc(1,reso/2+1));
% e22sc(reso/2+1,1)=real(e22sc(reso/2+1,1));
% e22sc(reso/2+2:reso,1)=conj(e22sc(reso/2:-1:2,1));
% e22sc(1,reso/2+2:reso)=conj(e22sc(1,reso/2:-1:2));
% e22sc(reso/2+1:reso,reso/2+1:reso)=conj(e22sc(reso/2+1:-1:2,reso/2+1:-1:2));
% e22sc(reso/2+2:reso,2:reso/2)=conj(e22sc(reso/2:-1:2,reso:-1:reso/2+2));
% %--------------------------------------------------------------------------
% % e12sc_l=
% % cl_ee^2/(cl_ee+cl_np)*lx*ly*sin(2phi)*cos(2phi)
% e12sc_l=cl_eem.*ly_m.*lx_m.*sin2phi_m.*cos2phi_m;
% e12sc_l(reso/2+1,reso/2+1)=0;
% e12sc_l(reso/2+2:reso,:)=e12sc_l(reso/2:-1:2,:);
% e12sc=ifft2(ifftshift(e12sc_l))./(delta_x^2);
% e12sc(reso/2+1,reso/2+1)=0;
% e12sc(1,1)=real(e12sc(1,1));
% e12sc(1,reso/2+1)=real(e12sc(1,reso/2+1));
% e12sc(reso/2+1,1)=real(e12sc(reso/2+1,1));
% e12sc(reso/2+2:reso,1)=conj(e12sc(reso/2:-1:2,1));
% e12sc(1,reso/2+2:reso)=conj(e12sc(1,reso/2:-1:2));
% e12sc(reso/2+1:reso,reso/2+1:reso)=conj(e12sc(reso/2+1:-1:2,reso/2+1:-1:2));
% e12sc(reso/2+2:reso,2:reso/2)=conj(e12sc(reso/2:-1:2,reso:-1:reso/2+2));
%--------------------------------------------------------------------------
% b00sc_l=
% 1/(cl_bb+cl_np)*sin(2phi)*cos(2phi)
b00sc_l=cl_bbm.*cos2phi_m.*sin2phi_m;
b00sc_l(reso/2+1,reso/2+1)=0;
b00sc_l(reso/2+2:reso,:)=-b00sc_l(reso/2:-1:2,:);
b00sc=ifft2(ifftshift(b00sc_l))./(delta_x^2);
b00sc(reso/2+1,reso/2+1)=0;
b00sc(1,1)=real(b00sc(1,1));
b00sc(1,reso/2+1)=real(b00sc(1,reso/2+1));
b00sc(reso/2+1,1)=real(b00sc(reso/2+1,1));
b00sc(reso/2+2:reso,1)=conj(b00sc(reso/2:-1:2,1));
b00sc(1,reso/2+2:reso)=conj(b00sc(1,reso/2:-1:2));
b00sc(reso/2+1:reso,reso/2+1:reso)=conj(b00sc(reso/2+1:-1:2,reso/2+1:-1:2));
b00sc(reso/2+2:reso,2:reso/2)=conj(b00sc(reso/2:-1:2,reso:-1:reso/2+2));
%--------------------------------------------------------------------------
% b00c_l=
% 1/(cl_bb+cl_np)*cos(2phi)^2
b00c_l=cl_bbm.*cos2phi_m.^2;
b00c_l(reso/2+1,reso/2+1)=0;
b00c_l(reso/2+2:reso,:)=b00c_l(reso/2:-1:2,:);
b00c=ifft2(ifftshift(b00c_l))./(delta_x^2);
b00c(reso/2+1,reso/2+1)=0;
b00c(1,1)=real(b00c(1,1));
b00c(1,reso/2+1)=real(b00c(1,reso/2+1));
b00c(reso/2+1,1)=real(b00c(reso/2+1,1));
b00c(reso/2+2:reso,1)=conj(b00c(reso/2:-1:2,1));
b00c(1,reso/2+2:reso)=conj(b00c(1,reso/2:-1:2));
b00c(reso/2+1:reso,reso/2+1:reso)=conj(b00c(reso/2+1:-1:2,reso/2+1:-1:2));
b00c(reso/2+2:reso,2:reso/2)=conj(b00c(reso/2:-1:2,reso:-1:reso/2+2));
%--------------------------------------------------------------------------
% b00s_l=
% 1/(cl_bb+cl_np)*sin(2phi)^2
b00s_l=cl_bbm.*sin2phi_m.^2;
b00s_l(reso/2+1,reso/2+1)=0;
b00s_l(reso/2+2:reso,:)=b00s_l(reso/2:-1:2,:);
b00s=ifft2(ifftshift(b00s_l))./(delta_x^2);
b00s(reso/2+1,reso/2+1)=0;
b00s(1,1)=real(b00s(1,1));
b00s(1,reso/2+1)=real(b00s(1,reso/2+1));
b00s(reso/2+1,1)=real(b00s(reso/2+1,1));
b00s(reso/2+2:reso,1)=conj(b00s(reso/2:-1:2,1));
b00s(1,reso/2+2:reso)=conj(b00s(1,reso/2:-1:2));
b00s(reso/2+1:reso,reso/2+1:reso)=conj(b00s(reso/2+1:-1:2,reso/2+1:-1:2));
b00s(reso/2+2:reso,2:reso/2)=conj(b00s(reso/2:-1:2,reso:-1:reso/2+2));
%--------------------------------------------------------------------------
%
%
%--------------------------------------------------------------------------
g11=e11s.*b00c+e11c.*b00s-2*e11sc.*b00sc;
% g22=e22s.*b00c+e22c.*b00s-2*e22sc.*b00sc;
% g12=e12s.*b00c+e12c.*b00s-2*e12sc.*b00sc;
g11_l=real(fftshift(fft2(g11))*delta_x^2); 
% g22_l=real(fftshift(fft2(g22))*delta_x^2);
% g12_l=real(fftshift(fft2(g12))*delta_x^2);

cl_nnm ...
    = g11_l(1:reso/2+1,:).*lx_m.^2;% ...
%      +g22_l(1:reso/2+1,:).*ly_m.^2 ...
%      +g12_l(1:reso/2+1,:).*lx_m.*ly_m;

cl_nnm=cl_nnm.^(-1);
cl_nnm(reso/2+1,reso/2+1)=0;
cl_nn=cl_nnm(reso/2+1,reso/2+1:reso);
l_out=0:lreso:lmax-lreso;
end

