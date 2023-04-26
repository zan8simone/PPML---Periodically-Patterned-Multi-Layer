clearvars; close all;
% Reproduces fig. 2a-b of "Multichannel remote polarization control 
%               enabled by nanostructured Liquid Crystalline Networks",
%               APL, May 2019.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is free software distributed under the BSD licence (see the 
%  containing folder).
% However, shall the results obtained through this code be included 
%  in an academic publication, we kindly ask you to cite the source 
%  website and, if applicable, the following paper:
%
% Simone Zanotto, Fabrizio Sgrignuoli, Sara Nocentini, Daniele Martella,
% Camilla Parmeggiani, Diederik S. Wiersma, "Multichannel remote polarization 
%         control enabled by nanostructured Liquid Crystalline Networks",
%               Applied Physics Letters, May 2019.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('PPML_root')    % to be replaced by the proper path 

dLCEv = [0:20:2500];
wstripe = 500; ne = 1.59; no = 1.47;

lambda = 633;      % HeNe
       
theta = 0.001;     % angle in degrees (never set 0 also for pol initialization)
phi   = 0;         % angle in degrees 

halfnpw = 5;          % good for convergence in this problem
                      % NOTE: this parameter is computation-dependent.
                      % Be careful in drawing conclusions!!!

epsLCEo = no^2;     % LCE permittivity tensor
epsLCEe = ne^2;
epsXLCE = epsLCEo;  % 
epsYLCE = epsLCEe;  % 
epszLCE = epsLCEo;

psi = 45;            % LC alignment rotation w.r. to pattern
epsxxLCE =  epsXLCE*cos(psi*pi/180)^2 + epsYLCE*sin(psi*pi/180)^2;
epsyyLCE =  epsXLCE*sin(psi*pi/180)^2 + epsYLCE*cos(psi*pi/180)^2;
epsxyLCE = (epsXLCE - epsYLCE)*cos(psi*pi/180)*sin(psi*pi/180);
                      
for i = 1:length(dLCEv)

dLCE = dLCEv(i);

a = 1500; % nm
L = 1;    % number of internal layers 

ff = wstripe/a; 
%           glass   | LCE     |    air          
f      = [           ff                   ]; % 
d      = [1500       dLCE         1500    ]; % nm
                   
epsA   = [          1                      ]; % material A 
epsxxB = [          epsxxLCE               ]; % material B
epsxyB = [          epsxyLCE               ]; % material B
epsyyB = [          epsyyLCE               ]; % material B
epszB  = [          epszLCE                ]; % material B
sigma  = [      0               0          ];

epssup = 1.44^2;  
epssub = 1;                                

k0    = 2*pi/lambda;         % wavevector in nm ^-1
kparx = k0*sin(theta*pi/180)*cos(phi*pi/180);
kpary = k0*sin(theta*pi/180)*sin(phi*pi/180);

% calling the PPML function
epar = epar_1d(a,L,...
   epssup,epssub,epsA,epsxxB,epsxyB,epsyyB,epszB,sigma,...
   f,d,halfnpw,k0,kparx,kpary,'p');

% extracting the field amplitudes for transmitted and diffracted waves.
npw = 2*halfnpw+1; 
ExI_T   = epar(npw+halfnpw+1);
ExI_Dp1 = epar(npw+halfnpw+1+1);
ExI_Dm1 = epar(npw+halfnpw+1-1);
EyI_T   = -epar(halfnpw+1);
EyI_Dp1 = -epar(halfnpw+1+1);
EyI_Dm1 = -epar(halfnpw+1-1);

% correction factor due to wave impedance difference between superstrate
%                                                       and substrate
ExIII_T = ExI_T/sqrt(sqrt(epssup));
EyIII_T = EyI_T/sqrt(sqrt(epssup));

ctI  = sqrt(1-(2*pi/(k0*a))^2); % diffraction angles 
ExIII_D = ExI_Dp1/ctI/sqrt(sqrt(epssup));
EyIII_D = EyI_Dp1/sqrt(sqrt(epssup));

% Stokes parameter for transmitted and diffracted waves
S0T(i) = ExIII_T*conj(ExIII_T) + EyIII_T*conj(EyIII_T); 
S1T(i) = ExIII_T*conj(ExIII_T) - EyIII_T*conj(EyIII_T); 
S2T(i) = ExIII_T*conj(EyIII_T) + EyIII_T*conj(ExIII_T); 
S0D(i) = ExIII_D*conj(ExIII_D) + EyIII_D*conj(EyIII_D); 
S1D(i) = ExIII_D*conj(ExIII_D) - EyIII_D*conj(EyIII_D); 
S2D(i) = ExIII_D*conj(EyIII_D) + EyIII_D*conj(ExIII_D); 
end

clearvars -except dLCEv S0T S1T S2T S0D S1D S2D

%%%%%%%%%%%%% plotting

cD = [100 0 0; 255 0 0; 255 160 160]/256;
cT = [70 1 164; 165 98 255; 215 186 255]/256;
cGrey = [1 1 1]*0.3;

hf = figure(1); set(hf,'Position',[200 200 450 360])
subplot(121)
plot([0 2.500],[0 0], 'Color',cGrey); hold on; grid on;
plot(dLCEv/1000,S0T,'Color',cT(1,:)); hold on;
plot(dLCEv/1000,S1T,'Color',cT(2,:)); hold on;
plot(dLCEv/1000,S2T,'Color',cT(3,:)); hold on;
axis([0 2.500 -1 1])
ax = gca; ax.XTick = [0 0.5 1 1.5 2 2.5];
xt = get(gca, 'XTick'); set(gca, 'FontSize', 9.5)
subplot(122)
plot([0 2.500],[0 0], 'Color',cGrey); hold on;  grid on;
plot(dLCEv/1000,S0D,'Color',cD(1,:)); hold on;
plot(dLCEv/1000,S1D,'Color',cD(2,:)); hold on;
plot(dLCEv/1000,S2D,'Color',cD(3,:)); hold on;
axis([0 2.500 -0.4 0.45])
ax = gca; ax.XTick = [0 0.5 1 1.5 2 2.5];
xt = get(gca, 'XTick'); set(gca, 'FontSize', 9.5)
