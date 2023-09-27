function optical_critical(a, h, theta, nu)
    % Calculates and plots the angularly-resolved reflectance
    % theoretical counterpart of Fig. 3a of "Optical critical coupling into
    %                highly confining metal-insulator-metal
    %                optical resonators" Appl. Phys. Lett. 103, 091110, (2013)
    %
    % Simone Zanotto, Firenze, feb. 2016
    % a: stripe spacing, aka lambda (microns)
    % h: etching depth
    % theta: angles in degrees
    % nu: frequencies in THz


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This is free software distributed under the BSD licence (see the
    %  containing folder).
    % However, shall the results obtained through this code be included
    %  in an academic publication, we kindly ask you to cite the source
    %  website and, if applicable, the following paper:
    %
    % J.-M- Manceau, S. Zanotto, I. Sagnes, G. Beaudoin, and R. Colombelli,
    %  "Optical critical coupling into highly confining metal-insulator-metal
    %   optical resonators" Appl. Phys. Lett. 103, 091110, (2015)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % Computation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numnu = length(nu);   % Number of frequency points
    numtheta = length(theta); % Number of angles of incidence

    halfnpw = 20;         % good for convergence in this problem
    % NOTE: this parameter is computation-dependent.
    % Be careful in drawing conclusions!!!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    L = 4;   % number of internal layers
    R = zeros(numnu, numtheta);

    %        superstr.  |  patt. Au  |  patt. GaAs  |  unp. GaAs  | dop. GaAs | subst.
    f     = [              .77          .77            1            1                       ]; % fraction of B in A
    d     = [    30        .2           h              (9-h)        1                 30    ]; % microns
    epsxA = [              1            1              1            1                       ]; % material A
    epszA = [              1            1              1            1                       ]; %
    sigma = [          0            0              0            0                0          ];

    epssup = 1;
    wb = waitbar(0, 'Starting', 'Name', 'Calculating reflectance');

    for i = 1:numnu
        % dispersive dielectric constants
        % mind the signs of imag parts (misprint in the APL)
        wn = nu(i)*100/3;        % wavenumber in cm^-1 !!!!!!
        om = 2*pi*nu(i)*1e12;       % ang freq in Hz
        epsAu = 1 - (7.25e4)^2/( wn^2 + 1i*216*wn ); % Ordal Au
        epsGaAs = 11*(1 + (292^2 - 268^2)/(268^2 - wn^2 - 1i*2.4*wn)); % GaAs with phonon
        epsGaAs_dop = 11 - (2.6e24*1.6e-19^2/(8.8e-12*0.067*9.1e-31*(om^2 + 1i*om/100e-15))); % doped GaAs
        %        superstr.  |  patt. Au  |  patt. GaAs  |  unp. GaAs  | dop. GaAs | subst.
        epsxB = [              epsAu        epsGaAs        epsGaAs      epsGaAs_dop             ]; % material B
        epszB = [              epsAu        epsGaAs        epsGaAs      epsGaAs_dop             ]; %
        epssub = epsAu;
        k0   = 2*pi*nu(i)*0.01/3;         % wavevector in micron ^-1

        for j = 1:numtheta
            % main cycles over energy and angles
            waitbar((i/numnu) + (j/(numnu*numtheta)), wb, sprintf("Frequency = %1.2f THz", nu(i)));

            kpar = k0*sin(theta(j)*pi/180);

            % calling the scattering matrix subroutine for reflection
            [RR,~,~] = RTA_1d_tm(a,L,...
                epssup,epssub,epsxA,epszA,epsxB,epszB,sigma,f,d,...
                halfnpw,k0,kpar);
            R(i,j) = RR;

        end
    end
    close(wb)
    %%
    if numtheta > 7
        clf
        s = pcolor(theta,nu,R);
        set(s,'edgecolor','none')
        xlabel('Incidence angle (deg)')
        ylabel('Frequency (THz)')
        colorbar
        title('Reflectance')
    else
        plot_tiled_reflectance(R, theta, nu)
    end
end
