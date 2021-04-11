function d=blanco_m(L,m,mp,bta)

% store the Wigner d function for certain m & mp.
% the first element of d is for l=lmin
% the last element of d is for l=L.
lmin=max(abs(m),abs(mp));
lmin1=lmin+1;
btasize=size(bta);
N=btasize(2);

d=zeros(L-lmin+1,N);  
dcell=blanco_v(lmin1,bta);
lmin_i=lmin+1;
mi=lmin+m+1;
mpi=lmin+mp+1;
lmin1_i=lmin1+1;
mi1=lmin1+m+1;
mpi1=lmin1+mp+1;

if (L-lmin+1==1)
    d(1,:)=dcell{lmin_i}(mi,mpi,:);
elseif (L-lmin+1==1)
    d(1,:)=dcell{lmin_i}(mi,mpi,:);
    d(2,:)=dcell{lmin1_i}(mi1,mpi1,:);
else
d(1,:)=dcell{lmin_i}(mi,mpi,:);
d(2,:)=dcell{lmin1_i}(mi1,mpi1,:);

btar=bta/180*pi;
d1_00=cos(btar);

for k=1:N
if abs(btar(k)-pi/2)<eps
  d1_00(k)=0;
end
end

for l=lmin+2:L
   li=l-lmin+1; 
   fac1=l*(2*l-1)/sqrt((l^2-m^2)*(l^2-mp^2));
   fac2=d1_00-m*mp/l/(l-1);
   fac3=sqrt(((l-1)^2-m^2)*((l-1)^2-mp^2))/(l-1)/(2*l-1); 
   d(li,:)=fac1*(fac2.*d(li-1,:)...
	     -fac3*d(li-2,:));
   % avoid round-off error
      
   
end

   for l=lmin:L
       li=l-lmin+1;
   for k=1:N
   if (abs(d(li,k))<eps)
      d(li,k)=0; 
   end
   end
   end

end
                                           