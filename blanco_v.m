function d=blanco_v(L,bta)
% [DR,DC]=BLANCO(L,bta)
%
% Constructs spherical harmonics rotation matrices for use with spherical
% harmonics coefficients, using the algorithm of Blanco, Florez & Bermejo (1997).
%
% INPUT:
% 
% L        Maximum degree of the harmonics
% alp      Euler angle in degrees
% bta      Euler angle in degrees [default: 90]
% gam      Euler angle in degrees
%
% OUTPUT:
%
% DR       Rotation matrix for m>=0 only (lower right quarter)
% DC       Rotation matrix for m=-l:l (full matrix)
%
% EXAMPLE:
%
% Compare the difference; it works for the decomposition where it's just
% ninety 
%
% L=round(rand*180)
% d1=blanco(L,90); d2=dlmb(L); 
% for index=1:L+1; difer(d1{index}-d2{index}); end
%
% See also DLMB, PLM2ROT
%
% Last modified by fjsimons-at-alum.mit.edu, January 12th, 2003

%defval('bta',90)
bta=bta*pi/180;
btasize=size(bta);
N=btasize(2);

d=cell(1,L);
for l=1:L
    li=l+1;
    d{li}=zeros(2*l+1,2*l+1,N);
end
% if bta<0
%   warning('Only for positive rotations')
% end

% Eq. (48-52)
d{1}(1,1,1:N)=1;
d{2}(2,2,1:N)=cos(bta);
d{2}(3,1,1:N)=sin(bta/2).^2;
d{2}(3,2,1:N)=-1/sqrt(2)*sin(bta);
d{2}(3,3,1:N)=cos(bta/2).^2;
% Avoid round-off error
for k=1:N
if abs(bta(k)-pi/2)<eps
  d{2}(2,2,k)=0;
  d{2}(3,1,k)=1/2;
  d{2}(3,3,k)=1/2;
  
  if abs(bta(k)-pi)<eps
  d{2}(3,2,k)=0;
  d{2}(3,3,k)=0;
  end
end
end


% Use l, m, and mp as degrees and
% li, mi, mpi as indices of d{l}
% li1, mi1, mpi1 as indices of d{l-1}
% li2, mi2, mpi2 as indices of d{l-2}
% Loop over degrees l
for l=2:L
  % Start with zero to initialize
  li=l+1;
  li1=li-1;
  li2=li-2;
  % Use Eq. 64
  for m=0:l-2
    % Index of m for l
    mi =m+l+1;
    % Index of m for l-1
    mi1=m+l;
    % Index of m for l-2
    mi2=m+l-1;
    for mp=-m:m
      mpi =mp+l+1;
      mpi1=mp+l;
      mpi2=mp+l-1;
      fac1=l*(2*l-1)/sqrt((l^2-m^2)*(l^2-mp^2));
      fac2=d{2}(2,2,1:N)-m*mp/l/(l-1);
      fac3=sqrt(((l-1)^2-m^2)*((l-1)^2-mp^2))/(l-1)/(2*l-1);
      d{li}(mi,mpi,1:N)=fac1*(fac2.*d{li1}(mi1,mpi1,1:N)...
			  -fac3*d{li2}(mi2,mpi2,1:N));
    end
  end
  % Eq. 65
  % Last index of current l, double
  m=l;  mi=m+l+1;
  % Last index of previous l, double
  mi11=m+l-1; 
  d{li}(mi,mi,1:N)=d{2}(3,3,1:N).*d{li1}(mi11,mi11,1:N);
  % Eq. 66
  % One but last index of current l
  m=l-1; mi=m+l+1; mi11=m+l;
  d{li}(mi,mi,1:N)=(l*d{2}(2,2,1:N)-l+1).*d{li1}(mi11,mi11,1:N); 
  % Eq. 67
  % Last row of the matrix except its last column
  m=l;
  mi=m+l+1;
  for mp=l-1:-1:-l
    mpi=mp+l+1;
%     warning off
    fac1=-sqrt((l+mp+1)/(l-mp-1+1))...
	 *sqrt(d{2}(3,1,1:N)./d{2}(3,3,1:N));
%     warning on
    for k=1:N
    if isinf(sqrt(d{2}(3,1,k)/d{2}(3,3,k)))
      if mp==-l
	d{li}(mi,mpi,k)=1;
      else
	d{li}(mi,mpi,k)=0;
      end
    else
      d{li}(mi,mpi,k)=fac1(k)*d{li}(mi,mpi+1,k);
    end
    if isinf(fac1(k))
      d{li}(mi,mpi,k)=0;     
    end
    if isnan(fac1(k))
      d{li}(mi,mpi,k)=0;
    end
    end
  end
  % Eq. 68
  % One but last row
  m=l-1;
  mi=m+l+1;
  for mp=(l-2):-1:(1-l)
    mpi=mp+l+1;
%     warning off
    fac1=-(l*d{2}(2,2,1:N)-mp-1+1)./(l*d{2}(2,2,1:N)-mp-1)*...
	 sqrt((l+mp+1)/(l-mp-1+1))...
	 .*sqrt(d{2}(3,1,1:N)./d{2}(3,3,1:N));
%     warning on
    for k=1:N
    if isinf(sqrt(d{2}(3,1,k)/d{2}(3,3,k)))
      if mp==1-l
	d{li}(mi,mpi,k)=l-1;
      else
	d{li}(mi,mpi,k)=0;
      end
    else
      d{li}(mi,mpi,k)=fac1(k)*d{li}(mi,mpi+1,k);
    end
    if isinf(fac1(k))
      d{li}(mi,mpi,k)=0;     
    end
    if isnan(fac1(k))
      d{li}(mi,mpi,k)=0;
    end
    end
  end
end
% Something is wrong - either the l-1 elements are not good
% or we're not using enough of the symmetry properties.
% Symmetrize using mirror properties of Eq. 38
for l=1:L
  li=l+1;
  for k=1:N
  % First symmetry equation: flip across LLUR diagonal
  d{li}(:,:,k)=(d{li}(:,:,k)+rot90(d{li}(:,:,k)',2))./(flipud(eye(2*l+1))+1);
  % Second symmetry relation
  mmp=repmat(-l:l,2*l+1,1)+repmat([-l:l]',1,2*l+1);
  d{li}(:,:,k)=(d{li}(:,:,k)+(-1).^mmp.*d{li}(:,:,k)')./(eye(2*l+1)+1);
  end
end

% Collect only lower right quarter 
% for l=0:L
%   li=l+1;
%   dr{li}=d{li}(li:end,li:end);
% end
% 
% varn={'dr','d'};
% for index=1:nargout
%   varargout{index}=eval(varn{index});
% end
