function cl_nn=cal_recon_noise(lmax,cl_ee,cl_bb,cl_enp,cl_bnp,d33,d3m3,d31,d3m1,d11,d1m1,d22,d2m2,weight)

% lmin=2;
% lsize=lmax-lmin+1;
% N=ceil((lmax+1)/2);
lvector=0:lmax;
lvector=lvector';
% the needed Wigner d function
% d33=zeros(N,lsize-1);
% d3m3=zeros(N,lsize-1);
% d31=zeros(N,lsize-1);
% d3m1=zeros(N,lsize-1);
% d11=zeros(N,lsize+1);
% d1m1=zeros(N,lsize+1);
% d22=zeros(N,lsize);
% d2m2=zeros(N,lsize);

% get the weights and abscissas 
% [weight,z,N]=gausslegendrecof(3*lmax,'jacobi',[-1 1]);
% 
% % let's cheat to obtain weight when lmax is large (for test)!
% 
% nanw=isnan(weight);
% infw=isinf(weight);
% ini=1;
% jni=1;
% for k=1:N
% if (nanw(k)==1 || infw(k)==1)
%     xi(ini)=k;
%     ini=ini+1;
% else
%     x(jni)=k;
%     y(jni)=weight(k);
%     jni=jni+1;
% end
% end
% 
% 
% weight(xi)=interp1(x,y,xi,'spline');
% 
% thetav=acos(z)*180/pi;

enn=cl_ee.^2./(cl_ee+cl_enp);
bnn=1./(cl_bb+cl_bnp);

lfactor=(2*lvector+1)/4/pi;
e33_l=lfactor(4:lmax+1).*enn(4:lmax+1) ...
    .*(lvector(4:lmax+1)-2).*(lvector(4:lmax+1)+3);
e31_l=lfactor(4:lmax+1).*enn(4:lmax+1) ...
    .*sqrt((lvector(4:lmax+1)-1).*(lvector(4:lmax+1)+2) ...
    .*(lvector(4:lmax+1)-2).*(lvector(4:lmax+1)+3));
e11_l=lfactor(2:lmax+1).*enn(2:lmax+1) ...
    .*(lvector(2:lmax+1)-1).*(lvector(2:lmax+1)+2);
b_l=lfactor(3:lmax+1).*bnn(3:lmax+1);

% for k=1:N
%    theta=thetav(k);
%    d33=blanco_m(lmax,3,3,thetav')';
    zeta_E33=d33*e33_l;
%    clear d33
%    d3m3=blanco_m(lmax,3,-3,thetav')';
    zeta_E3m3=d3m3*e33_l;
%    clear d3m3
%    d31=blanco_m(lmax,3,1,thetav')';
    zeta_E31=d31*e31_l;
%    clear d31
%    d3m1=blanco_m(lmax,3,-1,thetav')';
    zeta_E3m1=d3m1*e31_l;
%    clear d3m1
%    d11=blanco_m(lmax,1,1,thetav')';
    zeta_E11=d11*e11_l;
%    %clear d11
%    d1m1=blanco_m(lmax,1,-1,thetav')';
    zeta_E1m1=d1m1*e11_l;
%    %clear d1m1
%    d22=blanco_m(lmax,2,2,thetav')';
    zeta_Bp=d22*b_l;
%    clear d22
%    d2m2=blanco_m(lmax,2,-2,thetav')';
    zeta_Bm=d2m2*b_l;
%    clear d2m2
    fun1=d11'*((zeta_E33.*zeta_Bp-2*zeta_E3m1.*zeta_Bm+zeta_E11.*zeta_Bp).*weight);  
%    clear d11   
    fun2=d1m1'*((zeta_E3m3.*zeta_Bm-2*zeta_E31.*zeta_Bp+zeta_E1m1.*zeta_Bm).*weight);  
%    clear d1m1
% end

% compute Eq. (36~39)
    
cl_nn=pi/4*(fun1-fun2).*lvector(2:lmax+1).*(lvector(2:lmax+1)+1);
cl_nn=1./cl_nn;
