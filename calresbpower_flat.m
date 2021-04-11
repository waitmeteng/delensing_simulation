function [ cl_res, l_out ] = calresbpower_flat(lmax,lreso,cl_eeg,cl_enp,cl_phi,cl_nn)
% -------------------------------------------------------------------------
% Generate the residual B power predicted by
% Smith et al. (2010), arXiv:1010.0048 using flat approximation.
% -------------------------------------------------------------------------
% INPUT:
%
% .The maximum L
% .The resolution in Fourier space
% .Unlensed EE power
% .Instrumental noise (micronK-arcmin)
% .Lensing potential power
% .Reconstruction noise for lensing potential
%
% OUTPUT:
%
% .Residual B power
%
% -------------------------------------------------------------------------
% Wei-Hsiang Teng, NTU, Dec 2010
% -------------------------------------------------------------------------
reso=2*floor(lmax/lreso);
angu_size=360/lreso;
delta_x=angu_size/reso*pi/180;
%--------------------------------------------------------------------------
% Define the needed matrix
%--------------------------------------------------------------------------
% l_m=zeros(reso/2+1,reso);
lx_m=zeros(reso/2+1,reso);
ly_m=zeros(reso/2+1,reso);
cl_eem=zeros(reso/2+1,reso);
cl_phim=zeros(reso/2+1,reso);
cl_eem1=zeros(reso/2+1,reso);
cl_phim1=zeros(reso/2+1,reso);
sin2phi_m2=zeros(reso/2+1,reso);
for i=1:reso
    for j=1:reso/2+1
        kx=i-(reso/2+1);
        ky=j-(reso/2+1);
        k=sqrt(kx^2+ky^2);
        ll=k*lreso;
        lx_m(j,i)=kx*lreso;
        ly_m(j,i)=ky*lreso;
        phi=atan2(ky,kx);
        sin2phi_m2(j,i)=(1-cos(4*phi))/2;
%         l_m(j,i)=round(ll); 
        if (ll<=lmax)
           cl_eem(j,i)=cl_eeg(round(ll)+1);  
           cl_phim(j,i)=cl_phi(round(ll)+1);  
           cl_eem1(j,i)=cl_eeg(round(ll)+1).^2 ...
               ./(cl_eeg(round(ll)+1)+cl_enp(round(ll)+1));  
           cl_phim1(j,i)=cl_phi(round(ll)+1).^2 ...
               ./(cl_phi(round(ll)+1)+cl_nn(round(ll)+1));
        end
    end
end

% cl_eem=cl_eeg(l_m+1);
% cl_phim=cl_phi(l_m+1);
%--------------------------------------------------------------------------
% e11_l=
% cl_ee*lx^2*sin(2phi)^2
e11_l=cl_eem.*lx_m.^2.*sin2phi_m2;
e11_l(reso/2+1,reso/2+1)=0;
e11_l(reso/2+2:reso,:)=e11_l(reso/2:-1:2,:);
e11=ifft2(ifftshift(e11_l))./(delta_x^2);
e11(reso/2+1,reso/2+1)=0;
e11(1,1)=real(e11(1,1));
e11(1,reso/2+1)=real(e11(1,reso/2+1));
e11(reso/2+1,1)=real(e11(reso/2+1,1));
e11(reso/2+2:reso,1)=conj(e11(reso/2:-1:2,1));
e11(1,reso/2+2:reso)=conj(e11(1,reso/2:-1:2));
e11(reso/2+1:reso,reso/2+1:reso)=conj(e11(reso/2+1:-1:2,reso/2+1:-1:2));
e11(reso/2+2:reso,2:reso/2)=conj(e11(reso/2:-1:2,reso:-1:reso/2+2));
%--------------------------------------------------------------------------
% % e22_l=
% % cl_ee*ly^2*sin(2phi)^2
e22_l=cl_eem.*ly_m.^2.*sin2phi_m2;
e22_l(reso/2+1,reso/2+1)=0;
e22_l(reso/2+2:reso,:)=e22_l(reso/2:-1:2,:);
e22=ifft2(ifftshift(e22_l))./(delta_x^2);
e22(reso/2+1,reso/2+1)=0;
e22(1,1)=real(e22(1,1));
e22(1,reso/2+1)=real(e22(1,reso/2+1));
e22(reso/2+1,1)=real(e22(reso/2+1,1));
e22(reso/2+2:reso,1)=conj(e22(reso/2:-1:2,1));
e22(1,reso/2+2:reso)=conj(e22(1,reso/2:-1:2));
e22(reso/2+1:reso,reso/2+1:reso)=conj(e22(reso/2+1:-1:2,reso/2+1:-1:2));
e22(reso/2+2:reso,2:reso/2)=conj(e22(reso/2:-1:2,reso:-1:reso/2+2));
%--------------------------------------------------------------------------
% % e12_l=
% % cl_ee*lx*ly*sin(2phi)^2
e12_l=cl_eem.*lx_m.*ly_m.*sin2phi_m2;
e12_l(reso/2+1,reso/2+1)=0;
e12_l(reso/2+2:reso,:)=-e12_l(reso/2:-1:2,:);
e12=ifft2(ifftshift(e12_l))./(delta_x^2);
e12(reso/2+1,reso/2+1)=0;
e12(1,1)=real(e12(1,1));
e12(1,reso/2+1)=real(e12(1,reso/2+1));
e12(reso/2+1,1)=real(e12(reso/2+1,1));
e12(reso/2+2:reso,1)=conj(e12(reso/2:-1:2,1));
e12(1,reso/2+2:reso)=conj(e12(1,reso/2:-1:2));
e12(reso/2+1:reso,reso/2+1:reso)=conj(e12(reso/2+1:-1:2,reso/2+1:-1:2));
e12(reso/2+2:reso,2:reso/2)=conj(e12(reso/2:-1:2,reso:-1:reso/2+2));
%--------------------------------------------------------------------------
% phi11_l=
% cl_phi*lx^2
phi11_l=cl_phim.*lx_m.^2;
phi11_l(reso/2+1,reso/2+1)=0;
phi11_l(reso/2+2:reso,:)=phi11_l(reso/2:-1:2,:);
phi11=ifft2(ifftshift(phi11_l))./(delta_x^2);
phi11(reso/2+1,reso/2+1)=0;
phi11(1,1)=real(phi11(1,1));
phi11(1,reso/2+1)=real(phi11(1,reso/2+1));
phi11(reso/2+1,1)=real(phi11(reso/2+1,1));
phi11(reso/2+2:reso,1)=conj(phi11(reso/2:-1:2,1));
phi11(1,reso/2+2:reso)=conj(phi11(1,reso/2:-1:2));
phi11(reso/2+1:reso,reso/2+1:reso)=conj(phi11(reso/2+1:-1:2,reso/2+1:-1:2));
phi11(reso/2+2:reso,2:reso/2)=conj(phi11(reso/2:-1:2,reso:-1:reso/2+2));
%--------------------------------------------------------------------------
% phi22_l=
% cl_phi*ly^2
phi22_l=cl_phim.*ly_m.^2;
phi22_l(reso/2+1,reso/2+1)=0;
phi22_l(reso/2+2:reso,:)=phi22_l(reso/2:-1:2,:);
phi22=ifft2(ifftshift(phi22_l))./(delta_x^2);
phi22(reso/2+1,reso/2+1)=0;
phi22(1,1)=real(phi22(1,1));
phi22(1,reso/2+1)=real(phi22(1,reso/2+1));
phi22(reso/2+1,1)=real(phi22(reso/2+1,1));
phi22(reso/2+2:reso,1)=conj(phi22(reso/2:-1:2,1));
phi22(1,reso/2+2:reso)=conj(phi22(1,reso/2:-1:2));
phi22(reso/2+1:reso,reso/2+1:reso)=conj(phi22(reso/2+1:-1:2,reso/2+1:-1:2));
phi22(reso/2+2:reso,2:reso/2)=conj(phi22(reso/2:-1:2,reso:-1:reso/2+2));
%--------------------------------------------------------------------------
% phi12_l=
% cl_phi*lx*ly
phi12_l=cl_phim.*lx_m.*ly_m;
phi12_l(reso/2+1,reso/2+1)=0;
phi12_l(reso/2+2:reso,:)=-phi12_l(reso/2:-1:2,:);
phi12=ifft2(ifftshift(phi12_l))./(delta_x^2);
phi12(reso/2+1,reso/2+1)=0;
phi12(1,1)=real(phi12(1,1));
phi12(1,reso/2+1)=real(phi12(1,reso/2+1));
phi12(reso/2+1,1)=real(phi12(reso/2+1,1));
phi12(reso/2+2:reso,1)=conj(phi12(reso/2:-1:2,1));
phi12(1,reso/2+2:reso)=conj(phi12(1,reso/2:-1:2));
phi12(reso/2+1:reso,reso/2+1:reso)=conj(phi12(reso/2+1:-1:2,reso/2+1:-1:2));
phi12(reso/2+2:reso,2:reso/2)=conj(phi12(reso/2:-1:2,reso:-1:reso/2+2));
%--------------------------------------------------------------------------
g=e11.*phi11+e22.*phi22+2*e12.*phi12;
g_l=real(fftshift(fft2(g))*delta_x^2);
% cl_res=g_l(reso/2+1,reso/2+1:reso);
% cl_res(1)=0;
%--------------------------------------------------------------------------
% cl_eem=cl_eeg(l_m+1).^2./(cl_eeg(l_m+1)+cl_enp(l_m+1));
% cl_phim=cl_phi(l_m+1).^2./(cl_phi(l_m+1)+cl_nn(l_m+1));
%--------------------------------------------------------------------------
% e11_l=
% cl_ee*lx^2*sin(2phi)^2
e11_l=cl_eem1.*lx_m.^2.*sin2phi_m2;
e11_l(reso/2+1,reso/2+1)=0;
e11_l(reso/2+2:reso,:)=e11_l(reso/2:-1:2,:);
e11=ifft2(ifftshift(e11_l))./(delta_x^2);
e11(reso/2+1,reso/2+1)=0;
e11(1,1)=real(e11(1,1));
e11(1,reso/2+1)=real(e11(1,reso/2+1));
e11(reso/2+1,1)=real(e11(reso/2+1,1));
e11(reso/2+2:reso,1)=conj(e11(reso/2:-1:2,1));
e11(1,reso/2+2:reso)=conj(e11(1,reso/2:-1:2));
e11(reso/2+1:reso,reso/2+1:reso)=conj(e11(reso/2+1:-1:2,reso/2+1:-1:2));
e11(reso/2+2:reso,2:reso/2)=conj(e11(reso/2:-1:2,reso:-1:reso/2+2));
%--------------------------------------------------------------------------
% % e22_l=
% % cl_ee*ly^2*sin(2phi)^2
e22_l=cl_eem1.*ly_m.^2.*sin2phi_m2;
e22_l(reso/2+1,reso/2+1)=0;
e22_l(reso/2+2:reso,:)=e22_l(reso/2:-1:2,:);
e22=ifft2(ifftshift(e22_l))./(delta_x^2);
e22(reso/2+1,reso/2+1)=0;
e22(1,1)=real(e22(1,1));
e22(1,reso/2+1)=real(e22(1,reso/2+1));
e22(reso/2+1,1)=real(e22(reso/2+1,1));
e22(reso/2+2:reso,1)=conj(e22(reso/2:-1:2,1));
e22(1,reso/2+2:reso)=conj(e22(1,reso/2:-1:2));
e22(reso/2+1:reso,reso/2+1:reso)=conj(e22(reso/2+1:-1:2,reso/2+1:-1:2));
e22(reso/2+2:reso,2:reso/2)=conj(e22(reso/2:-1:2,reso:-1:reso/2+2));
%--------------------------------------------------------------------------
% % e12_l=
% % cl_ee*lx*ly*sin(2phi)^2
e12_l=cl_eem1.*lx_m.*ly_m.*sin2phi_m2;
e12_l(reso/2+1,reso/2+1)=0;
e12_l(reso/2+2:reso,:)=-e12_l(reso/2:-1:2,:);
e12=ifft2(ifftshift(e12_l))./(delta_x^2);
e12(reso/2+1,reso/2+1)=0;
e12(1,1)=real(e12(1,1));
e12(1,reso/2+1)=real(e12(1,reso/2+1));
e12(reso/2+1,1)=real(e12(reso/2+1,1));
e12(reso/2+2:reso,1)=conj(e12(reso/2:-1:2,1));
e12(1,reso/2+2:reso)=conj(e12(1,reso/2:-1:2));
e12(reso/2+1:reso,reso/2+1:reso)=conj(e12(reso/2+1:-1:2,reso/2+1:-1:2));
e12(reso/2+2:reso,2:reso/2)=conj(e12(reso/2:-1:2,reso:-1:reso/2+2));
%--------------------------------------------------------------------------
% phi11_l=
% cl_phi*lx^2
phi11_l=cl_phim1.*lx_m.^2;
phi11_l(reso/2+1,reso/2+1)=0;
phi11_l(reso/2+2:reso,:)=phi11_l(reso/2:-1:2,:);
phi11=ifft2(ifftshift(phi11_l))./(delta_x^2);
phi11(reso/2+1,reso/2+1)=0;
phi11(1,1)=real(phi11(1,1));
phi11(1,reso/2+1)=real(phi11(1,reso/2+1));
phi11(reso/2+1,1)=real(phi11(reso/2+1,1));
phi11(reso/2+2:reso,1)=conj(phi11(reso/2:-1:2,1));
phi11(1,reso/2+2:reso)=conj(phi11(1,reso/2:-1:2));
phi11(reso/2+1:reso,reso/2+1:reso)=conj(phi11(reso/2+1:-1:2,reso/2+1:-1:2));
phi11(reso/2+2:reso,2:reso/2)=conj(phi11(reso/2:-1:2,reso:-1:reso/2+2));
%--------------------------------------------------------------------------
% phi22_l=
% cl_phi*ly^2
phi22_l=cl_phim1.*ly_m.^2;
phi22_l(reso/2+1,reso/2+1)=0;
phi22_l(reso/2+2:reso,:)=phi22_l(reso/2:-1:2,:);
phi22=ifft2(ifftshift(phi22_l))./(delta_x^2);
phi22(reso/2+1,reso/2+1)=0;
phi22(1,1)=real(phi22(1,1));
phi22(1,reso/2+1)=real(phi22(1,reso/2+1));
phi22(reso/2+1,1)=real(phi22(reso/2+1,1));
phi22(reso/2+2:reso,1)=conj(phi22(reso/2:-1:2,1));
phi22(1,reso/2+2:reso)=conj(phi22(1,reso/2:-1:2));
phi22(reso/2+1:reso,reso/2+1:reso)=conj(phi22(reso/2+1:-1:2,reso/2+1:-1:2));
phi22(reso/2+2:reso,2:reso/2)=conj(phi22(reso/2:-1:2,reso:-1:reso/2+2));
%--------------------------------------------------------------------------
% phi12_l=
% cl_phi*lx*ly
phi12_l=cl_phim1.*lx_m.*ly_m;
phi12_l(reso/2+1,reso/2+1)=0;
phi12_l(reso/2+2:reso,:)=-phi12_l(reso/2:-1:2,:);
phi12=ifft2(ifftshift(phi12_l))./(delta_x^2);
phi12(reso/2+1,reso/2+1)=0;
phi12(1,1)=real(phi12(1,1));
phi12(1,reso/2+1)=real(phi12(1,reso/2+1));
phi12(reso/2+1,1)=real(phi12(reso/2+1,1));
phi12(reso/2+2:reso,1)=conj(phi12(reso/2:-1:2,1));
phi12(1,reso/2+2:reso)=conj(phi12(1,reso/2:-1:2));
phi12(reso/2+1:reso,reso/2+1:reso)=conj(phi12(reso/2+1:-1:2,reso/2+1:-1:2));
phi12(reso/2+2:reso,2:reso/2)=conj(phi12(reso/2:-1:2,reso:-1:reso/2+2));
%--------------------------------------------------------------------------
g1=e11.*phi11+e22.*phi22+2*e12.*phi12;
g1_l=real(fftshift(fft2(g1))*delta_x^2);
cl_res=g_l(reso/2+1,reso/2+1:reso)-g1_l(reso/2+1,reso/2+1:reso);
cl_res(1)=0;
l_out=0:lreso:lmax-lreso;
end


