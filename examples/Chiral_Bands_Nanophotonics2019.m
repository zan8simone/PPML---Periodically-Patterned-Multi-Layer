close all;
clearvars; clc
% Reproduces fig. 1c of "Photonic bands, superchirality, and inverse design
%         of a chiral minimal metasurface", Nanophotonics, 8(12), 2291-2301 (2019).  
%
%  Simone Zanotto, Pisa, oct. 2019


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is free software distributed under the BSD licence (see the 
%  containing folder).
% However, shall the results obtained through this code be included 
%  in an academic publication, we kindly ask you to cite the source 
%  website and, if applicable, the following paper:
%
% Simone Zanotto, Giacomo Mazzamuto, Francesco Riboli, Giorgio Biasiol, 
%         Giuseppe C. La Rocca, Alessandro Tredicucci, and Alessandro Pitanti, 
%         Photonic bands, superchirality, and inverse design
%         of a chiral minimal metasurface",  Nanophotonics 2019 
%                              DOI: https://doi.org/10.1515/nanoph-2019-0321
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%addpath('PPML_root')    % to be replaced by the proper path 
addpath(genpath('C:\Users\Simone\Documents\Z_Pisa2\Codici\PPML\PPML_v3.0'))

epsOx   = 1.5^2; epsGaAs = 3.37^2;
dOx = 5; dGaAs = 211;   

lambdav = [1200:20:1900];

a   = 1134;          
ff1 = 0.580;
ff2 = 0.760;
ff3 = 0.327;
ff4 = 0.327;
%%%%%%%%%%%

theta = 0.01;     % angle in degrees (never set 0)
phi   = 90;      % angle in degrees 

parfor i = 1:length(lambdav)
% Computation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i
lambda = lambdav(i);

halfnpw = 5;       

L = 3;     % number of internal layers 

%          air  |  oxide  |  patt. GaAs |  oxide |   air          
f1    = [          ff1       ff1           ff1            ]; % fraction of B in A
f2    = [          ff2       ff2           ff2            ]; % fraction of B in A
f3    = [          ff3       ff3           ff3            ]; % fraction of B in A
f4    = [          ff4       ff4           ff4            ]; % fraction of B in A
d     = [    0     dOx       dGaAs         dOx     0      ]; % nm
                                       
epsA  = [           epsOx    epsGaAs       epsOx          ]; % material A 
epsB  = [           1        1             1              ]; % material B
sigma = [       0         0            0           0      ];

epssup = 1;                                
epssub = 1;


k0   = 2*pi/lambda;         % wavevector in nm ^-1
kparx = k0*sin(theta*pi/180)*cos(phi*pi/180);
kpary = k0*sin(theta*pi/180)*sin(phi*pi/180);

% calling the scattering matrix subroutine for reflection
S = ZSM_2d_Lshape(a,L,...
   epssup,epssub,epsA,epsB,sigma,f1,f2,f3,f4,d,...
   halfnpw,k0,kparx,kpary);
T{i} = S(3:4,1:2);

end


clearvars a d dGaAs dOx epsA epsB epsGaAs epsOx epssup epssub f1 f2 f3 f4 ff1 ff2 ff3 ff4
clearvars halfnpw i j k0 k_GaAs kparx kpary L lambda_GaAs n_GaAs lambda S phi theta

save tmp
%%
load tmp
close all;


Lambda = (1/sqrt(2))*[1,  1; 1i, -1i];

for i = 1:length(lambdav)
Tcirc  = (Lambda')*T{i}*Lambda;
TR(i)  = sum( abs( Tcirc*[1; 0] ).^2 );
TL(i)  = sum( abs( Tcirc*[0; 1] ).^2 );
end

plot(lambdav,TR,'b'); hold on;
plot(lambdav,TL,'r'); hold on;

