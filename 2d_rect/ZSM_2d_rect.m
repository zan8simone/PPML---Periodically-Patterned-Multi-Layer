%%
function S = ZSM_2d_rect(ax,ay,L,...
   epssup,epssub,epsA,epsB,sigma,fx,fy,d,...
   halfnpw,k0,kparx,kpary)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is free software distributed under the BSD licence (see the 
%  containing folder).
% However, shall the results obtained through this code be included 
%  in an academic publication, we kindly ask you to cite the source 
%  website.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zero-order 4x4 scattering matrix   
%
% 2-d "rectangular" structure, rectangular lattice
%
% Simone Zanotto, Pisa, Jan-Feb 2013, Milano, June 2015, Florence, nov 2016
%                  and Pisa, 2017-2021
%                   
% -------------------------------------------------------------------------
%                              INPUTS
%
% ax             -> periodicity along x
% ay             -> periodicity along y
% L              -> number of layers excluded superstrate and substrate
% epssup, epssub -> super-and sub-strate permittivities 
% epsA           -> dielectric constant 
%                   for material A in the internal layers
%                                         [vector with L components]
% epsB           -> idem, for material B   
% sigma          -> conductivity of interfaces
%                                         [vector with (L+1) components]
% fx             -> x- duty cycle in the internal layers
%                                         [vector with L components]
% fy             -> y- duty cycle in the internal layers
%                                         [vector with L components]
% d              -> thicknesses of the layers, 
%                   including super- and sub-strate
%                                         [vector with (L+2) components]
% halfnpw        -> truncation order (square scheme adopted)
% k0             -> wavevector in vacuum  (= 2*pi/lambda0)
% kparx          -> x-projection of the incident wavevector 
% kpary          -> y-projection of the incident wavevector 
%
%%% ----------------------------------------------------------------------- 
%%%                             OUTPUTS
%
% S              -> zero-order scattering matrix
%                                              [4x4 complex]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% NOTES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% *the lengths (a,d) and inverse lengths (k0,kpar) must be set in the same
%  units (e.g. microns and inverse microns, respectively)
%
% *halfnpw = 0 sets the number of harmonic waves to 1, e.g. the scattering
%  matrix reduces to the ordinary 4x4 formalism for unpatterned multilayers
%
% *thickness of superstrate and substrate do not influence the result. 
%  The electromagnetic problem is solved de facto assuming semi-infinite 
%  super- and substrate boundary conditions.
%  You find them in view of further development (field profile calculation).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Building the square lattice

npw = (2*halfnpw+1)^2;
nx = zeros(npw,1); ny = zeros(npw,1);
i = 1;
for ix = -halfnpw:halfnpw
for iy = -halfnpw:halfnpw
nx(i) = ix; ny(i) = iy;
i = i+1;
end
end

normm = nx.^2+ny.^2;
[norm_sort,indices] = sort(normm);
nx_sort = nx(indices);
ny_sort = ny(indices);
nx = nx_sort;
ny = ny_sort;

clear nx_sort ny_sort indices norm_sort normm ix iy i 
%%% now   nx ny npw   are defined 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:npw 
for j = 1:npw
nx_ij(i,j) = (nx(i)-nx(j));
ny_ij(i,j) = (ny(i)-ny(j));
end
end

kx = (2*pi/ax)*nx + kparx;
ky = (2*pi/ay)*ny + kpary;

% Initializing the work variables
q   = zeros(2*npw,2*npw,L+2);
phi = zeros(2*npw,2*npw,L+2);
A   = zeros(2*npw,2*npw,L+2);

% Eigensolutions for super- and substrate
q(1:npw,1:npw,1)     = diag(sqrt_whittaker(epssup*k0^2 - kx.*kx - ky.*ky));
q(npw+1:2*npw,npw+1:2*npw,1) = q(1:npw,1:npw,1);
phi(:,:,1) = eye(2*npw);
eps = epssup;

kstorto = [(diag(ky)/eps)*diag(ky), -(diag(ky)/eps)*diag(kx);...
          -(diag(kx)/eps)*diag(ky), (diag(kx)/eps)*diag(kx)];
A(:,:,1)   = ((eye(2*npw)*k0^2 - kstorto)*phi(:,:,1))/q(:,:,1);

q(1:npw,1:npw,L+2)     = diag(sqrt_whittaker(epssub*k0^2 - kx.*kx - ky.*ky));
q(npw+1:2*npw,npw+1:2*npw,L+2) = q(1:npw,1:npw,L+2);
phi(:,:,L+2) = eye(2*npw);
eps = epssub;

kstorto = [(diag(ky)/eps)*diag(ky), -(diag(ky)/eps)*diag(kx);...
          -(diag(kx)/eps)*diag(ky), (diag(kx)/eps)*diag(kx)];
A(:,:,L+2)   = ((eye(2*npw)*k0^2 - kstorto)*phi(:,:,L+2))/q(:,:,L+2);


% Eigensolutions for internal layers
if L > 0
for l = 1:L

   F = fx(l)*fy(l)*sinc(fx(l)*nx_ij).*sinc(fy(l)*ny_ij);

   
   eps = (epsB(l) - epsA(l))*F + epsA(l)*eye(npw);
   
   
   Fx = zeros(2*halfnpw+1);
   Fy = zeros(2*halfnpw+1);
   for i = 1:2*halfnpw+1
   for j = 1:2*halfnpw+1
   Fx(i,j) = fx(l)*sinc(fx(l)*(i-j));
   Fy(i,j) = fy(l)*sinc(fy(l)*(i-j));
   end
   end
   
   cjan = inv((1/epsB(l) - 1/epsA(l))*Fx + eye(2*halfnpw+1)/epsA(l)) ...
           - epsA(l)*eye(2*halfnpw+1); 
   gjat = inv((1/epsB(l) - 1/epsA(l))*Fy + eye(2*halfnpw+1)/epsA(l)) ...
           - epsA(l)*eye(2*halfnpw+1);    
   
   eps_xy = zeros(npw);
   eps_yx = zeros(npw);
   for i = 1:npw 
   for j = 1:npw
   eps_xy(i,j) = cjan(nx(i)+halfnpw+1,nx(j)+halfnpw+1)*...
       Fy(ny(i)+halfnpw+1,ny(j)+halfnpw+1);
   eps_yx(i,j) = gjat(ny(i)+halfnpw+1,ny(j)+halfnpw+1)*...
       Fx(nx(i)+halfnpw+1,nx(j)+halfnpw+1);
   end
   end   
   
   eps_xy = eps_xy + eye(npw)*epsA(l);
   eps_yx = eps_yx + eye(npw)*epsA(l);
     
   % Eigensolution calculation
   kstorto = [(diag(ky)/eps)*diag(ky), -(diag(ky)/eps)*diag(kx);...
          -(diag(kx)/eps)*diag(ky), (diag(kx)/eps)*diag(kx)];
   kdritto = [diag(kx)*diag(kx), diag(kx)*diag(ky);...
              diag(ky)*diag(kx), diag(ky)*diag(ky)];
   epsilone = [eps_yx, zeros(npw); zeros(npw), eps_xy];
   [phhi,qq] = eig(epsilone*(eye(2*npw)*k0^2 - kstorto) ...
                    - kdritto);
   q(:,:,l+1) = diag(sqrt_whittaker(diag(qq)));
   phi(:,:,l+1) = phhi;
   A(:,:,l+1) = (eye(2*npw)*k0^2 - kstorto)*phhi/q(:,:,l+1);
   
end
end


% Forward propagation 
S1 = zeros(2*npw,2*npw,L+2); S2 = zeros(2*npw,2*npw,L+2);
S1(:,:,1) = eye(2*npw);
for l = 1:L+1
[S1(:,:,l+1),S2(:,:,l+1)] = smpropag_fw_cond(S1(:,:,l),S2(:,:,l),...
                                phi(:,:,l),phi(:,:,l+1),...
                                A(:,:,l),A(:,:,l+1),...
                                exp(1i*diag(q(:,:,l)*d(l))),exp(1i*diag(q(:,:,l+1)*d(l+1))),...
								sigma(l)*376.730/k0);
end

% Backward propagation
S3 = zeros(2*npw,2*npw,L+2); S4 = zeros(2*npw,2*npw,L+2);
S4(:,:,L+2) = eye(2*npw);
for ll = 1:L+1
    l = L+3-ll;
[S3(:,:,l-1),S4(:,:,l-1)] = smpropag_bw_cond(S3(:,:,l),S4(:,:,l),...
                                phi(:,:,l),phi(:,:,l-1),...
                                A(:,:,l),A(:,:,l-1),...
                                exp(1i*diag(q(:,:,l)*d(l))),exp(1i*diag(q(:,:,l-1)*d(l-1))),...
								sigma(l-1)*376.730/k0);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% calculating E-field S parameter.
% some working variables
kparvec  = [kparx, kpary, 0];
kpsup = [kparx, kpary,  q(1,1,1)];
kmsup = [kparx, kpary, -q(1,1,1)];

spsup = -cross(kparvec,[0,0,1]);
spsup = spsup/norm(spsup);
ppsup = cross(kpsup,spsup);
ppsup = ppsup/norm(ppsup);

smsup = spsup;
pmsup = cross(kmsup,smsup);
pmsup = pmsup/norm(pmsup);


kpsub = [kparx, kpary, -q(1,1,L+2)];
kmsub = [kparx, kpary,  q(1,1,L+2)];

spsub = spsup;
ppsub = cross(kpsub,spsub);
ppsub = ppsub/norm(ppsub);

smsub = spsup;
pmsub = cross(kmsub,smsub);
pmsub = pmsub/norm(pmsub);

X = zeros(1,2*npw); X(npw+1) = 1;
Y = zeros(1,2*npw); Y(1)     = 1;
Z = kpary*Y - kparx*X;

M12 = [ (-smsup(1)*X + smsup(2)*Y)*A(:,:,1)   + (smsup(3)/epssup)*Z ; ...
        (-pmsup(1)*X + pmsup(2)*Y)*A(:,:,1)   + (pmsup(3)/epssup)*Z];
		
M21 = [ ( smsub(1)*X - smsub(2)*Y)*A(:,:,L+2) + (smsub(3)/epssub)*Z ; ...
        ( pmsub(1)*X - pmsub(2)*Y)*A(:,:,L+2) + (pmsub(3)/epssub)*Z];
		
M = [ zeros(2,2*npw) , M12 ; M21 , zeros(2,2*npw)];

N11 = [ ( ( spsup(1)*X - spsup(2)*Y)*A(:,:,1)   + (spsup(3)/epssup)*Z )*exp(1i*q(1,1,1)*d(1)) ; ...
        ( ( ppsup(1)*X - ppsup(2)*Y)*A(:,:,1)   + (ppsup(3)/epssup)*Z )*exp(1i*q(1,1,1)*d(1))];

N22 = [ ( (-spsub(1)*X + spsub(2)*Y)*A(:,:,L+2) + (spsub(3)/epssub)*Z )*exp(1i*q(1,1,L+2)*d(L+2)) ; ...
        ( (-ppsub(1)*X + ppsub(2)*Y)*A(:,:,L+2) + (ppsub(3)/epssub)*Z )*exp(1i*q(1,1,L+2)*d(L+2))];


		
N = [ N11 , zeros(2,2*npw) ; zeros(2,2*npw) , N22 ];


SS = [S1(:,:,L+2), S2(:,:,L+2); ...
      S3(:,:,1)  , S4(:,:,1)  ];
	  
S = M*SS/N;
