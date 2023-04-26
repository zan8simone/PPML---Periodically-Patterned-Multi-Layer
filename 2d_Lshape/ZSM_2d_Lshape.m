%%
function S = ZSM_2d_Lshape(a,L,...
   epssup,epssub,epsA,epsB,sigma,f1,f2,f3,f4,d,...
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
% 2-d L-shaped structure
%
% Simone Zanotto, Pisa, Jan-Feb 2013; Milano, June 2015;
%                  and Pisa, 2017-2019
% -------------------------------------------------------------------------
%                              INPUTS
%
% a              -> periodicity 
% L              -> number of layers excluded superstrate and substrate
% epssup, epssub -> super-and sub-strate permittivities 
% epsA           -> dielectric constant 
%                   for material A in the internal layers
%                                         [vector with L components]
% epsB           -> idem, for material B   
% sigma          -> conductivity of interfaces
%                                         [vector with (L+1) components]
% f1, f2, f3 ,f4 -> define the shape of L inclusion (see notes)
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
% S           -> Zero-order scattering matrix. See manual for details.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% NOTES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% *the lengths (a,d) and inverse lengths (k0,kpar) must be set in the same
%  units (e.g. microns and inverse microns, respectively)
%
% *epssup must be real (otherwise the incident waves are ill-defined)
% 
% *epssub can be either real or complex. 
%
% *halfnpw = 0 sets the number of harmonic waves to 1, e.g. the scattering
%  matrix reduces to the ordinary 4x4 formalism for unpatterned multilayers
%
% *thickness of superstrate and substrate do not influence the result, apart for a 
% phase in reflection coefficients. The electromagnetic problem is solved assuming semi-infinite 
% super- and substrate boundary conditions.
%
% * Coefficients are determined assuming that the fields are evaluated at the interface between 
% superstrate and first internal layer. In the case there are no internal layers (L = 0), 
% the fields are evaluated at the interface between superstrate and substrate.
%
% *to run requires that the Fdata.mat file is in the folder of the calling
% script.  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if L>0
if (min(f1-f3) < -1e-3) || (min(f2-f4) < -1e-3)
warning('rxx_2d_Lshape:fconstr','Broken constraints on fs')
disp('error constraint replaced with warning');
end
end

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

kx = (2*pi/a)*nx + kparx;
ky = (2*pi/a)*ny + kpary;

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

[F, F1, F2, F3, F4, G1, G2, G13, G24] = fill_F(halfnpw,npw,nx,ny,f1(l),f2(l),f3(l),f4(l));

            
   eps = (epsB(l) - epsA(l))*F + epsA(l)*eye(npw);
 
 
   c0 = inv(eye(2*halfnpw+1)/epsA(l));
   c1 = inv(eye(2*halfnpw+1)/epsA(l) + (1/epsB(l) - 1/epsA(l))*F1);
   c2 = inv(eye(2*halfnpw+1)/epsA(l) + (1/epsB(l) - 1/epsA(l))*F2);
   c3 = inv(eye(2*halfnpw+1)/epsA(l) + (1/epsB(l) - 1/epsA(l))*F3);
   c4 = inv(eye(2*halfnpw+1)/epsA(l) + (1/epsB(l) - 1/epsA(l))*F4);
   
   
   eps_xy = zeros(npw);
   eps_yx = zeros(npw);
   for i = 1:npw 
   for j = 1:npw
   eps_xy(i,j) = c0(nx(i)+halfnpw+1,nx(j)+halfnpw+1) *G2(ny(i)+halfnpw+1,ny(j)+halfnpw+1) ...
               + c3(nx(i)+halfnpw+1,nx(j)+halfnpw+1)*G24(ny(i)+halfnpw+1,ny(j)+halfnpw+1) ...
               + c1(nx(i)+halfnpw+1,nx(j)+halfnpw+1) *F4(ny(i)+halfnpw+1,ny(j)+halfnpw+1);
           
   eps_yx(i,j) = c0(ny(i)+halfnpw+1,ny(j)+halfnpw+1) *G1(nx(i)+halfnpw+1,nx(j)+halfnpw+1) ...
               + c4(ny(i)+halfnpw+1,ny(j)+halfnpw+1)*G13(nx(i)+halfnpw+1,nx(j)+halfnpw+1) ...
               + c2(ny(i)+halfnpw+1,ny(j)+halfnpw+1) *F3(nx(i)+halfnpw+1,nx(j)+halfnpw+1);
   end
   end   
   
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

