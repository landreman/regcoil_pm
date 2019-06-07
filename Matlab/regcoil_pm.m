function regcoil_pm()

% This matlab script performs the same steps to the fortran program. The
% fortran and matlab versions are completely independent of each other. For
% identical inputs, they should give identical outputs to within roundoff
% error.

clear

symmetry_option = 3;
% 1 = Assume stellarator symmetry
% 2 = Test stellarator symmetry of inductance matrix
% 3 = Do not assume stellarator symmetry

load_bnorm = true;
%load_bnorm = false;

bnorm_filename = '/Users/mattland/Box Sync/MATLAB/bnorm.d23p4_tm';

% Resolution parameters:
% **********************************
%{
ntheta_plasma = 64;
ntheta_coil   = 64;
nzeta_plasma = 64;
nzeta_coil   = 64;
mpol_magnetization  = 12;
ntor_magnetization  = 12;
ns_magnetization = 1;
ns_integration = 2;
%}

ntheta_plasma = 30;
ntheta_coil   = 33;
nzeta_plasma = 36;
nzeta_coil   = 38;
mpol_magnetization  = 6;
ntor_magnetization  = 8;
ns_magnetization = 2;
ns_integration = 3;
d_initial = 0.01;
sign_normal = 1;

%{
ntheta_plasma = 38;
ntheta_coil   = 36;
nzeta_plasma = 33;
nzeta_coil   = 30;
mpol_magnetization  = 8;
ntor_magnetization  = 7;
ns_magnetization = 3;
ns_integration = 3;
d_initial = 0.005;
sign_normal = -1;
%}
%{
ntheta_plasma = 30;
ntheta_coil   = 33;
nzeta_plasma = 36;
nzeta_coil   = 38;
mpol_magnetization  = 6;
ntor_magnetization  = 8;
ns_magnetization = 1;
ns_integration = 5;
%}
%{
ntheta_plasma = 5;
ntheta_coil   = 4; %33;
nzeta_plasma = 3;
nzeta_coil   = 6; %38;
mpol_magnetization  = 2;
ntor_magnetization  = 3;
ns_magnetization = 1;
ns_integration = 5;
%}

% Options for the shape of the plasma surface:
% **********************************
geometry_option_plasma = 7;
R0_plasma = 3.0;
a_plasma = 1.0;
nfp_imposed = 1;
%woutFilename = 'C:\Users\landreman\Box Sync\MATLAB\20150601-01 Sfincs version 3\equilibria\wout_w7x_standardConfig.nc';
%woutFilename = '/Users/mattland/Box Sync/MATLAB/wout_d23p4_tm.nc';
%woutFilename = 'equilibria/wout_d23p4_tm.nc';
shape_filename_plasma = '../equilibria/tf_only_half_tesla.plasma';

% Options for the shape of the coil surface:
% **********************************
geometry_option_coil = 3;
R0_coil = 3.0;
a_coil = 1.7;
separation = 0.35;
%nescin_filename = 'nescin.w7x_standardConfig_separation0.3';
%nescin_filename = '/Users/mattland/Box Sync/MATLAB/nescin.w7x_winding_surface_from_Drevlak';
%nescin_filename = 'equilibria/nescin.w7x_winding_surface_from_Drevlak';
nescin_filename = '../equilibria/NCSX.vv';

%d_initial = 0.01;

s_integration_option = 'Gaussian';
%s_integration_option = 'Chebyshev';

% Options for the regularization parameter:
% **********************************
nlambda = 20;
lambda_min = 1e-20;
lambda_max = 1e-3;

% Plotting options:
% **********************************

plot_results = true;
%plot_results = false;

max_nlambda_for_contour_plots = 18;

%plot3DFigure = true;
plot3DFigure = false;

plotGrids = true;
%plotGrids = false;

plotVectors = true;
%plotVectors = false;

%stopAfterInitialPlots = true;
stopAfterInitialPlots = false;

figureOffset = 0;


% Options related to checking fortran version
% *******************************************

compareToFortran = true;
%compareToFortran = false;

%fortranNcFilename = 'C:\Users\landreman\Box Sync\MATLAB\bdistrib_out.compareToMatlab.nc';
%fortranNcFilename = '/Users/mattland/regcoil/examples/compareToMatlab1/regcoil_out.compareToMatlab1.nc';
%fortranNcFilename = '../examples/NCSX_low_resolution/regcoil_out.NCSX_low_resolution.nc';
%fortranNcFilename = '/Users/mattland/Box Sync/work19/20190526-01-testing_regcoil_pm/20190526-01-015-thetaZeta64_mpolNtor12_sMagnetization1_sIntegration2_Cheb_d0.01_1proc/regcoil_out.NCSX.nc';
%fortranNcFilename = '/Users/mattland/Box Sync/work19/20190526-01-testing_regcoil_pm/20190526-01-038-loRes_sMagnetization2_sIntegration3_Gauss_d0.01/regcoil_out.NCSX.nc';
%fortranNcFilename = '/Users/mattland/regcoil_pm/examples/compareToMatlab2/regcoil_out.compareToMatlab2.nc';
fortranNcFilename = '/Users/mattland/regcoil_pm/examples/compareToMatlab1_symmetry3/regcoil_out.compareToMatlab1_symmetry3.nc';

fortranComparisonThreshhold_abs = 1e-11;

% *************************************************************************
% *************************************************************************
% End of input parameters.
% *************************************************************************
% *************************************************************************

mu0 = 4*pi*(1e-7);


    function compareVariableToFortran(variableName, varargin)
        % Specify 'abs' as an argument to compare the absolute values.
        % This is useful for the singular vectors, which are only defined
        % up to a sign in practice.
        if ~ compareToFortran
            return
        end
        try
            fortranVariable = double(ncread(fortranNcFilename,variableName));
        catch
            fprintf(['*** Variable ',variableName,' does not exist in the fortran output file.\n'])
            return
        end
        matlabVariable = eval(variableName);
        assignin('base',[variableName,'_m'],matlabVariable)
        assignin('base',[variableName,'_f'],fortranVariable)
        if isrow(matlabVariable)
            matlabVariable = matlabVariable(:);
        end
        if isrow(fortranVariable)
            fortranVariable = fortranVariable(:);
        end
        try
            % Next lines will cause an exception if sizes are different:
            if nargin>1 && strcmp(varargin{1},'abs')
                differences = abs(abs(matlabVariable) - abs(fortranVariable)) > fortranComparisonThreshhold_abs;
            else
                differences = abs(matlabVariable - fortranVariable) > fortranComparisonThreshhold_abs;
                %differences = (abs(matlabVariable - fortranVariable) > fortranComparisonThreshhold_abs) && ;
            end
            if any(any(any(differences)))
                fprintf(['*** Variable ',variableName,' is the same size Matlab and fortran but differs in value. max|diff|=%g\n'],max(max(max(differences))))
            else
                fprintf(['    Variable ',variableName,' is the same in Matlab and fortran.\n'])
            end
        catch
            fprintf(['*** Variable ',variableName,' is a different size between Matlab and fortran.\n'])
        end
    end

compareVariableToFortran('ntheta_plasma')
compareVariableToFortran('ntheta_coil')
compareVariableToFortran('nzeta_plasma')
compareVariableToFortran('nzeta_coil')
compareVariableToFortran('geometry_option_plasma')
compareVariableToFortran('geometry_option_coil')
compareVariableToFortran('mpol_magnetization')
compareVariableToFortran('ntor_magnetization')
compareVariableToFortran('ns_magnetization')
compareVariableToFortran('ns_integration')

% *********************************************
% Set up range of lambda to try
% *********************************************

lambda = zeros(nlambda,1);
lambda(2:end) = lambda_min * exp((0:(nlambda-2))/(nlambda-2)*log(lambda_max/lambda_min));

compareVariableToFortran('nlambda')
compareVariableToFortran('lambda')

% *********************************************
% Initialize s grids:
% *********************************************

[s_magnetization, ~] = clencurt(ns_magnetization, 0, 1);
switch s_integration_option
    case 'Gaussian'
        [s_integration, s_weights] = glwt(ns_integration, 0, 1);
    case 'Chebyshev'
        [s_integration, s_weights] = clencurt(ns_integration, 0, 1);
    otherwise
        error('Unrecognized s_integration_option')
end
interpolate_magnetization_to_integration = m20130215_01_makeChebyshevInterpolationMatrix(ns_magnetization, 0, 1, s_integration);

compareVariableToFortran('s_magnetization')
compareVariableToFortran('s_integration')
compareVariableToFortran('s_weights')

% *********************************************
% Initialize the plasma surface:
% *********************************************

bfc = 0;
switch geometry_option_plasma
    case {0,1}
        % Plain axisymmetric circular torus
        nfp = nfp_imposed;
        mnmax = 2;
        xm = [0,1];
        xn = [0,0];
        rmnc = [R0_plasma; a_plasma];
        zmns = [0; a_plasma];
        whichSurface = 2;
        Rmajor_p = R0_plasma;
        
    case {2}
        % Load flux surface info from VMEC
        filename = woutFilename;
        ns = double(ncread(filename,'ns'));
        Rmajor_p = double(ncread(filename,'Rmajor_p'));
        nfp = double(ncread(filename,'nfp'));
        xm = double(ncread(filename,'xm'));
        xn = double(ncread(filename,'xn'));
        xm_nyq = double(ncread(filename,'xm_nyq'));
        xn_nyq = double(ncread(filename,'xn_nyq'));
        rmnc = double(ncread(filename,'rmnc'));
        zmns = double(ncread(filename,'zmns'));
        bmnc = double(ncread(filename,'bmnc'));
        mnmax = double(ncread(filename,'mnmax'));
        mnmax_nyq = double(ncread(filename,'mnmax_nyq'));
        
        whichSurface = ns; % Choose the outermost surface
        % Discard the other surfaces:
        rmnc = rmnc(:,whichSurface);
        zmns = zmns(:,whichSurface);
    case {7}
        % FOCUS format
        fid = fopen(shape_filename_plasma,'r');
        if fid<0
            fprintf('Error! Unable to open file %s\n',shape_filename_plasma)
            return
        end
        line = fgetl(fid);
        line = fgetl(fid);
        data = sscanf(line,'%d %d %d');
        mnmax_plasma = data(1);
        nfp = data(2);
        nbf = data(3);
        line = fgetl(fid);
        line = fgetl(fid);
        xm_plasma = zeros(mnmax_plasma,1);
        xn_plasma = zeros(mnmax_plasma,1);
        rmnc = zeros(mnmax_plasma,1);
        rmns = zeros(mnmax_plasma,1);
        zmnc = zeros(mnmax_plasma,1);
        zmns = zeros(mnmax_plasma,1);
        for j = 1:mnmax_plasma
            line = fgetl(fid);
            data = sscanf(line,'%d %d %g %g %g %g');
            xn_plasma(j) = data(1);
            xm_plasma(j) = data(2);
            rmnc(j) = data(3);
            rmns(j) = data(4);
            zmnc(j) = data(5);
            zmns(j) = data(6);
        end
        xn_plasma = xn_plasma * nfp;
        xm = xm_plasma;
        xn = xn_plasma;
        mnmax = mnmax_plasma;
        line = fgetl(fid);
        line = fgetl(fid);
        bfc = zeros(nbf,1);
        bfs = zeros(nbf,1);
        bfn = zeros(nbf,1);
        bfm = zeros(nbf,1);
        for j = 1:nbf
            line = fgetl(fid);
            data = sscanf(line,'%d %d %g %g');
            bfn(j) = data(1);
            bfm(j) = data(2);
            bfc(j) = data(3);
            bfs(j) = data(4);
        end
        bfn = bfn * nfp;
        fclose(fid);
    otherwise
        error('Invalid setting for geometry_option_plasma')
end

switch geometry_option_plasma
    case {2}
        % BNORM scales B_n by curpol=(2*pi/nfp)*bsubv(m=0,n=0)
        % where bsubv is the extrapolation to the last full mesh point of
        % bsubvmnc.  Let's undo this scaling now.
        bsubvmnc = ncread(woutFilename,'bsubvmnc');
        bsubv00 = 1.5*bsubvmnc(1,end) - 0.5*bsubvmnc(1,end-1);
        curpol = 2*pi/nfp*bsubv00;  % /1 since nfp=1.
        
        bvco = ncread(woutFilename,'bvco');
        net_poloidal_current_Amperes = (2*pi/mu0) * (1.5*bvco(end) - 0.5*bvco(end-1));
        fprintf('New value for net_poloidal_current_Amperes: %g\n',net_poloidal_current_Amperes)
    otherwise
        curpol = 1;
end


nzetal_plasma = nzeta_plasma * nfp;
nzetal_coil   = nzeta_coil   * nfp;

theta_plasma = linspace(0,2*pi,ntheta_plasma+1);
theta_plasma(end) = [];
zeta_plasma = linspace(0,2*pi/nfp,nzeta_plasma+1);
zeta_plasma(end) = [];
zetal_plasma = linspace(0,2*pi,nzetal_plasma+1);
zetal_plasma(end) = [];
[zetal_plasma_2D, theta_plasma_2D] = meshgrid(zetal_plasma, theta_plasma);

x = zeros(ntheta_plasma,nzetal_plasma);
y = zeros(ntheta_plasma,nzetal_plasma);
z = zeros(ntheta_plasma,nzetal_plasma);

dxdtheta = zeros(ntheta_plasma,nzetal_plasma);
dydtheta = zeros(ntheta_plasma,nzetal_plasma);
dzdtheta = zeros(ntheta_plasma,nzetal_plasma);

dxdzeta = zeros(ntheta_plasma,nzetal_plasma);
dydzeta = zeros(ntheta_plasma,nzetal_plasma);
dzdzeta = zeros(ntheta_plasma,nzetal_plasma);

for i=1:mnmax
    angle = xm(i)*theta_plasma_2D-xn(i)*zetal_plasma_2D;
    angle2 = zetal_plasma_2D;
    
    x = x + rmnc(i)*cos(angle).*cos(angle2);
    y = y + rmnc(i)*cos(angle).*sin(angle2);
    z = z + zmns(i)*sin(angle);
    
    dxdtheta = dxdtheta - xm(i)*rmnc(i)*sin(angle).*cos(angle2);
    dydtheta = dydtheta - xm(i)*rmnc(i)*sin(angle).*sin(angle2);
    dzdtheta = dzdtheta + xm(i)*zmns(i)*cos(angle);
    
    dxdzeta = dxdzeta + xn(i)*rmnc(i)*sin(angle).*cos(angle2) ...
        - rmnc(i)*cos(angle).*sin(angle2);
    dydzeta = dydzeta + xn(i)*rmnc(i)*sin(angle).*sin(angle2) ...
        + rmnc(i)*cos(angle).*cos(angle2);
    dzdzeta = dzdzeta - xn(i)*zmns(i)*cos(angle);
end

Nx = dydzeta .* dzdtheta - dzdzeta .* dydtheta;
Ny = dzdzeta .* dxdtheta - dxdzeta .* dzdtheta;
Nz = dxdzeta .* dydtheta - dydzeta .* dxdtheta;
norm_normal_plasma = sqrt(Nx.*Nx + Ny.*Ny + Nz.*Nz);
norm_normal_plasma = norm_normal_plasma(:,1:nzeta_plasma);
dtheta_plasma = theta_plasma(2)-theta_plasma(1);
dzeta_plasma = zeta_plasma(2)-zeta_plasma(1);
area_plasma = sum(sum(norm_normal_plasma)) * dtheta_plasma * dzeta_plasma * nfp;

r_plasma = zeros(3, ntheta_plasma, nzetal_plasma);
drdtheta_plasma = zeros(3, ntheta_plasma, nzetal_plasma);
drdzeta_plasma = zeros(3, ntheta_plasma, nzetal_plasma);
normal_plasma = zeros(3, ntheta_plasma, nzetal_plasma);

r_plasma(1,:,:) = x;
r_plasma(2,:,:) = y;
r_plasma(3,:,:) = z;
drdtheta_plasma(1,:,:) = dxdtheta;
drdtheta_plasma(2,:,:) = dydtheta;
drdtheta_plasma(3,:,:) = dzdtheta;
drdzeta_plasma(1,:,:) = dxdzeta;
drdzeta_plasma(2,:,:) = dydzeta;
drdzeta_plasma(3,:,:) = dzdzeta;
normal_plasma(1,:,:) = Nx;
normal_plasma(2,:,:) = Ny;
normal_plasma(3,:,:) = Nz;

compareVariableToFortran('nfp')
compareVariableToFortran('theta_plasma')
compareVariableToFortran('zeta_plasma')
compareVariableToFortran('zetal_plasma')

compareVariableToFortran('r_plasma')
compareVariableToFortran('drdtheta_plasma')
compareVariableToFortran('drdzeta_plasma')
compareVariableToFortran('normal_plasma')
compareVariableToFortran('norm_normal_plasma')
compareVariableToFortran('area_plasma')

% *********************************************
% Initialize the coil surface:
% *********************************************

    function [theta, zeta, zetal, theta_2D, zetal_2D, r, drdtheta, drdzeta, d2rdtheta2, d2rdthetadzeta, d2rdzeta2, ...
            normal, norm_normal, area, mean_curvature, Jacobian_coefficient, nX, nY, nZ] ...
            = initSurface(ntheta, nzeta, geometry_option, R0, a, separation, nescin_filename)
        
        nzetal = nzeta*nfp;
        theta = linspace(0,2*pi,ntheta+1);
        theta(end) = [];
        zeta = linspace(0,2*pi/nfp,nzeta+1);
        zeta(end) = [];
        zetal = linspace(0,2*pi,nzetal+1);
        zetal(end) = [];
        
        dtheta = theta(2)-theta(1);
        dzeta  = zeta(2)-zeta(1);
        [zetal_2D, theta_2D] = meshgrid(zetal, theta);
        
        x = zeros(size(theta_2D));
        y = zeros(size(theta_2D));
        z = zeros(size(theta_2D));
        dxdtheta = zeros(size(theta_2D));
        dydtheta = zeros(size(theta_2D));
        dzdtheta = zeros(size(theta_2D));
        dxdzeta = zeros(size(theta_2D));
        dydzeta = zeros(size(theta_2D));
        dzdzeta = zeros(size(theta_2D));
        d2xdtheta2 = zeros(size(theta_2D));
        d2ydtheta2 = zeros(size(theta_2D));
        d2zdtheta2 = zeros(size(theta_2D));
        d2xdthetadzeta = zeros(size(theta_2D));
        d2ydthetadzeta = zeros(size(theta_2D));
        d2zdthetadzeta = zeros(size(theta_2D));
        d2xdzeta2 = zeros(size(theta_2D));
        d2ydzeta2 = zeros(size(theta_2D));
        d2zdzeta2 = zeros(size(theta_2D));
        
        switch(geometry_option)
            case {0,1}
                if geometry_option == 0
                    R0_to_use = Rmajor_p;
                else
                    R0_to_use = R0;
                end
                
                x = (R0_to_use + a * cos(theta_2D)) .* cos(zetal_2D);
                y = (R0_to_use + a * cos(theta_2D)) .* sin(zetal_2D);
                z = a * sin(theta_2D);
                
                dxdtheta = -a * sin(theta_2D) .* cos(zetal_2D);
                dydtheta = -a * sin(theta_2D) .* sin(zetal_2D);
                dzdtheta = a * cos(theta_2D);
                
                dxdzeta = -(R0_to_use + a * cos(theta_2D)) .* sin(zetal_2D);
                dydzeta =  (R0_to_use + a * cos(theta_2D)) .* cos(zetal_2D);
                dzdzeta = zeros(size(theta_2D));
                
                if false
                    d2xdtheta2 = -a * cos(theta_2D) .* cos(zetal_2D);
                    d2ydtheta2 = -a * cos(theta_2D) .* sin(zetal_2D);
                    d2zdtheta2 = -a * sin(theta_2D);

                    d2xdthetadzeta =  a * sin(theta_2D) .* sin(zetal_2D);
                    d2ydthetadzeta = -a * sin(theta_2D) .* cos(zetal_2D);
                    d2zdthetadzeta = zeros(size(theta_2D));

                    d2xdzeta2 = -(R0_to_use + a * cos(theta_2D)) .* cos(zetal_2D);
                    d2ydzeta2 = -(R0_to_use + a * cos(theta_2D)) .* sin(zetal_2D);
                    d2zdzeta2 = zeros(size(theta_2D));
                end
            case 2
                error('geometry_option = 2 is not yet implemented for coil and outer surfaces.')
                
            case 3
                % Read coil surface from nescin file
                
                fid = fopen(nescin_filename,'r');
                search_string = '------ Current Surface';
                while true
                    line = fgetl(fid);
                    if strncmp(line,search_string,numel(search_string))
                        break
                    end
                end
                line = fgetl(fid); %Number of fourier modes in table
                line = fgetl(fid);
                mnmax_nescin = sscanf(line,'%d');
                fprintf('  Reading %d modes from nescin file %s\n',mnmax_nescin,nescin_filename)
                line = fgetl(fid); %Table of fourier coefficients
                line = fgetl(fid); %m,n,crc2,czs2,crs2,czc2
                xm_nescin = zeros(mnmax_nescin,1);
                xn_nescin = zeros(mnmax_nescin,1);
                rmnc_nescin = zeros(mnmax_nescin,1);
                zmns_nescin = zeros(mnmax_nescin,1);
                for i=1:mnmax_nescin
                    line = fgetl(fid);
                    data = sscanf(line,'%d %d %g %g %g %g %g %g');
                    xm_nescin(i) = data(1);
                    xn_nescin(i) = data(2);
                    rmnc_nescin(i) = data(3);
                    zmns_nescin(i) = data(4);
                end
                fclose(fid);
                % Done reading nescin file.
                                
                for i = 1:mnmax_nescin
                    angle = xm_nescin(i)*theta_2D + xn_nescin(i)*zetal_2D*nfp;
                    sinangle = sin(angle);
                    cosangle = cos(angle);
                    sinzeta = sin(zetal_2D);
                    coszeta = cos(zetal_2D);
                    
                    x = x + rmnc_nescin(i)*cosangle.*coszeta;
                    y = y + rmnc_nescin(i)*cosangle.*sinzeta;
                    z = z + zmns_nescin(i)*sinangle;
                    
                    dxdtheta = dxdtheta - xm_nescin(i)*rmnc_nescin(i)*sinangle.*coszeta;
                    dydtheta = dydtheta - xm_nescin(i)*rmnc_nescin(i)*sinangle.*sinzeta;
                    dzdtheta = dzdtheta + xm_nescin(i)*zmns_nescin(i)*cosangle;
                    
                    dxdzeta = dxdzeta - nfp*xn_nescin(i)*rmnc_nescin(i)*sinangle.*coszeta ...
                        - rmnc_nescin(i)*cosangle.*sinzeta;
                    dydzeta = dydzeta - nfp*xn_nescin(i)*rmnc_nescin(i)*sinangle.*sinzeta ...
                        + rmnc_nescin(i)*cosangle.*coszeta;
                    dzdzeta = dzdzeta + nfp*xn_nescin(i)*zmns_nescin(i)*cosangle;
                    
                    d2xdtheta2 = d2xdtheta2 - xm_nescin(i)*xm_nescin(i)*rmnc_nescin(i)*cosangle.*coszeta;
                    d2ydtheta2 = d2ydtheta2 - xm_nescin(i)*xm_nescin(i)*rmnc_nescin(i)*cosangle.*sinzeta;
                    d2zdtheta2 = d2zdtheta2 - xm_nescin(i)*xm_nescin(i)*zmns_nescin(i)*sinangle;

                    d2xdzeta2 = d2xdzeta2 - nfp*xn_nescin(i)*(nfp*xn_nescin(i))*rmnc_nescin(i)*cosangle.*coszeta ...
                        + 2*nfp*xn_nescin(i)*rmnc_nescin(i)*sinangle.*sinzeta ...
                        - rmnc_nescin(i)*cosangle.*coszeta;

                    d2xdthetadzeta = d2xdthetadzeta - nfp*xn_nescin(i)*xm_nescin(i)*rmnc_nescin(i)*cosangle.*coszeta ...
                         + xm_nescin(i)*rmnc_nescin(i)*sinangle.*sinzeta;
                    
                    d2ydzeta2 = d2ydzeta2 - nfp*xn_nescin(i)*nfp*xn_nescin(i)*rmnc_nescin(i)*cosangle.*sinzeta ...
                         - 2*nfp*xn_nescin(i)*rmnc_nescin(i)*sinangle.*coszeta ...
                         - rmnc_nescin(i)*cosangle.*sinzeta;
                     
                    d2ydthetadzeta = d2ydthetadzeta - nfp*xn_nescin(i)*xm_nescin(i)*rmnc_nescin(i)*cosangle.*sinzeta ...
                        - xm_nescin(i)*rmnc_nescin(i)*sinangle.*coszeta;
                    
                    d2zdzeta2 = d2zdzeta2 - nfp*xn_nescin(i)*nfp*xn_nescin(i)*zmns_nescin(i)*sinangle;
                    
                    d2zdthetadzeta = d2zdthetadzeta - nfp*xn_nescin(i)*xm_nescin(i)*zmns_nescin(i)*sinangle;
                end
                
            otherwise
                error('Invalid geometry_option')
        end
        
        Nx = dydzeta .* dzdtheta - dzdzeta .* dydtheta;
        Ny = dzdzeta .* dxdtheta - dxdzeta .* dzdtheta;
        Nz = dxdzeta .* dydtheta - dydzeta .* dxdtheta;
        
        r = zeros(3, ntheta, nzetal);
        drdtheta = zeros(3, ntheta, nzetal);
        drdzeta = zeros(3, ntheta, nzetal);
        normal = zeros(3, ntheta, nzetal);
        
        r(1,:,:) = x;
        r(2,:,:) = y;
        r(3,:,:) = z;
        drdtheta(1,:,:) = dxdtheta;
        drdtheta(2,:,:) = dydtheta;
        drdtheta(3,:,:) = dzdtheta;
        drdzeta(1,:,:) = dxdzeta;
        drdzeta(2,:,:) = dydzeta;
        drdzeta(3,:,:) = dzdzeta;
        normal(1,:,:) = Nx;
        normal(2,:,:) = Ny;
        normal(3,:,:) = Nz;
        
        if false
            d2rdtheta2 = zeros(3, ntheta, nzetal);
            d2rdthetadzeta = zeros(3, ntheta, nzetal);
            d2rdzeta2 = zeros(3, ntheta, nzetal);
            
            d2rdtheta2(1,:,:) = d2xdtheta2;
            d2rdtheta2(2,:,:) = d2ydtheta2;
            d2rdtheta2(3,:,:) = d2zdtheta2;
            
            d2rdthetadzeta(1,:,:) = d2xdthetadzeta;
            d2rdthetadzeta(2,:,:) = d2ydthetadzeta;
            d2rdthetadzeta(3,:,:) = d2zdthetadzeta;
            
            d2rdzeta2(1,:,:) = d2xdzeta2;
            d2rdzeta2(2,:,:) = d2ydzeta2;
            d2rdzeta2(3,:,:) = d2zdzeta2;
        else
            d2rdtheta2 = 0;
            d2rdthetadzeta = 0;
            d2rdzeta2 = 0;
        end
        
        norm_normal = sqrt(Nx.*Nx + Ny.*Ny + Nz.*Nz);
        dtheta = theta(2)-theta(1);
        dzeta  = zeta(2)-zeta(1);
        
        nX = Nx ./ norm_normal;
        nY = Ny ./ norm_normal;
        nZ = Nz ./ norm_normal;
        
        E = dxdtheta .* dxdtheta + dydtheta .* dydtheta + dzdtheta .* dzdtheta;
        F = dxdtheta .* dxdzeta  + dydtheta .* dydzeta  + dzdtheta .* dzdzeta;
        G = dxdzeta  .* dxdzeta  + dydzeta  .* dydzeta  + dzdzeta  .* dzdzeta;
        
        L = d2xdtheta2     .* nX + d2ydtheta2     .* nY + d2zdtheta2     .* nZ;
        M = d2xdthetadzeta .* nX + d2ydthetadzeta .* nY + d2zdthetadzeta .* nZ;
        P = d2xdzeta2      .* nX + d2ydzeta2      .* nY + d2zdzeta2      .* nZ;
        
        mean_curvature = (E.*P + G.*L - 2*F.*M) ./ (2*(E.*G - F.*F));
        Jacobian_coefficient = (M .* M - L .* P) ./ (norm_normal .* norm_normal);
        
        norm_normal = norm_normal(:,1:nzeta);
        mean_curvature = mean_curvature(:,1:nzeta);
        Jacobian_coefficient = Jacobian_coefficient(:,1:nzeta);
        area = sum(sum(norm_normal)) * dtheta * dzeta * nfp;

    end

tic
fprintf('Initializing coil surface.\n')
[theta_coil, zeta_coil, zetal_coil, theta_coil_2D, zetal_coil_2D, r_coil, drdtheta_coil, drdzeta_coil, ...
    d2rdtheta2_coil, d2rdthetadzeta_coil, d2rdzeta2_coil, normal_coil, norm_normal_coil, area_coil, ...
    mean_curvature_coil, Jacobian_coefficient, nX, nY, nZ] ...
    = initSurface(ntheta_coil, nzeta_coil, geometry_option_coil, R0_coil, a_coil, separation, nescin_filename);

dtheta_coil = theta_coil(2)-theta_coil(1);
dzeta_coil = zeta_coil(2)-zeta_coil(1);

d = d_initial * ones(ntheta_coil, nzeta_coil);
Jacobian_coil = zeros(ntheta_coil, nzeta_coil, ns_integration);
for js = 1:ns_integration
    Jacobian_coil(:,:,js) = d .* norm_normal_coil .* (-1 + sign_normal * 2 * s_integration(js) * d .* mean_curvature_coil ...
        + s_integration(js) * s_integration(js) * d .* d .* Jacobian_coefficient);
end
if any(any(Jacobian_coil >= 0)) 
    error('Jacobian_coil is not negative-definite!')
end
Jacobian_coil = abs(Jacobian_coil);

fprintf('Done. Took %g seconds.\n',toc)

compareVariableToFortran('sign_normal')
compareVariableToFortran('d_initial')
compareVariableToFortran('theta_coil')
compareVariableToFortran('zeta_coil')
compareVariableToFortran('zetal_coil')
compareVariableToFortran('r_coil')
compareVariableToFortran('drdtheta_coil')
compareVariableToFortran('drdzeta_coil')
compareVariableToFortran('normal_coil')
compareVariableToFortran('norm_normal_coil')
compareVariableToFortran('mean_curvature_coil')
compareVariableToFortran('Jacobian_coil')

% *********************************************
% Make 3D figure of surfaces
% *********************************************

%r_plasma = ncread(fortranNcFilename,'r_plasma');

if plot3DFigure
    r_plasma_toplot = r_plasma;
    r_coil_toplot = r_coil;
    
    % "Rotate" in theta so the seam in the plot is on the bottom
    nshift = round(ntheta_plasma*0.25);
    r_plasma_toplot = circshift(r_plasma_toplot, [0,nshift,0]);
    nshift = round(ntheta_coil*0.25);
    r_coil_toplot = circshift(r_coil_toplot, [0,nshift,0]);
    
    % Close surfaces for plotting:
    r_plasma_toplot(:,:,end+1) = r_plasma_toplot(:,:,1);
    r_plasma_toplot(:,end+1,:) = r_plasma_toplot(:,1,:);
    
    % For coil and outer surfaces, close in theta, but don't bother closing
    % in zeta:
    r_coil_toplot(:,end+1,:) = r_coil_toplot(:,1,:);
    
    
    
    mask = zetal_coil < 0.7*nfp;
    r_coil_toplot = r_coil_toplot(:,:,mask);
        
    figure(1 + figureOffset)
    clf
    set(gcf,'Color','w')
    faceColor = [1,0,0];
    surf(squeeze(r_plasma_toplot(1,:,:)),squeeze(r_plasma_toplot(2,:,:)),squeeze(r_plasma_toplot(3,:,:)),'FaceColor',faceColor,'EdgeColor','none')
    hold on
    if plotGrids
        plot3(squeeze(r_plasma_toplot(1,:,:)),squeeze(r_plasma_toplot(2,:,:)),squeeze(r_plasma_toplot(3,:,:)),'.r')
    end
    daspect([1,1,1])
    %shading interp
    axis vis3d
    hold on
    
    %faceColor = [1,0,1];
    faceColor = [0,1,0];
    surf(squeeze(r_coil_toplot(1,:,:)),squeeze(r_coil_toplot(2,:,:)),squeeze(r_coil_toplot(3,:,:)),'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.75)
    
    faceColor = [0,0,1];
    %surf(squeeze(r_outer_toplot(1,:,:)),squeeze(r_outer_toplot(2,:,:)),squeeze(r_outer_toplot(3,:,:)),'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.75)
    
    if plotGrids
        plot3(squeeze(r_coil_toplot(1,:,:)),squeeze(r_coil_toplot(2,:,:)),squeeze(r_coil_toplot(3,:,:)),'.m')
        %plot3(squeeze(r_outer_toplot(1,:,:)), squeeze(r_outer_toplot(2,:,:)), squeeze(r_outer_toplot(3,:,:)),'.b')
    end
    %surf(X_coil,Y_coil,Z_coil,'FaceColor',faceColor,'EdgeColor','none')
    %surf(X_coil,Y_coil,Z_coil,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.75)
    light
    lighting gouraud
    zoom(1.6)
    campos([  574.9370 -457.0244  424.3304])
    camva(1.0271)
    axis off
    
end

if stopAfterInitialPlots
    return
end

% *********************************************
% Set up Fourier arrays
% *********************************************

    function [mnmax, xm, xn] = setupFourierArrays(mpol,ntor)
        % xm is non-negative, while xn can be negative
        % xn is the rapidly increasing variable.
        
        % When xm=0, xn=1..ntor.
        % When xm>0, xn=-ntor..ntor.
        mnmax = ntor + mpol*(2*ntor+1);
        
        xm = zeros(mnmax,1);
        xn = zeros(mnmax,1);
        xn(1:ntor) = 1:ntor;
        nextIndex = ntor+1;
        for m = 1:mpol
            indices = nextIndex:(nextIndex+2*ntor);
            xm(indices) = m;
            xn(indices) = (-ntor):ntor;
            nextIndex = nextIndex + 2*ntor+1;
        end
        
    end

[mnmax_coil, xm_coil, xn_coil] = setupFourierArrays(mpol_magnetization, ntor_magnetization);
xn_coil = xn_coil * nfp;

compareVariableToFortran('symmetry_option')

% *********************************************
% Load BNORM data.
% *********************************************

Bnormal_from_TF_and_plasma_current = zeros(ntheta_plasma,nzeta_plasma);
[zeta_plasma_2D, theta_plasma_2D] = meshgrid(zeta_plasma, theta_plasma);
if numel(bfc)>1
    % Data must have been loaded from a FOCUS file.
    for j = 1:numel(bfc)
        angle = bfm(j) * theta_plasma_2D - bfn(j) * zeta_plasma_2D;
        Bnormal_from_TF_and_plasma_current = Bnormal_from_TF_and_plasma_current + bfc(j) * cos(angle) + bfs(j) * sin(angle);
    end
else
    if load_bnorm
        fid = fopen(bnorm_filename,'r');
        if fid<0
            error('Unable to open BNORM file %s.\n',bnorm_filename)
        end
        while ~ feof(fid);
            [fileline,count] = fscanf(fid,'%f %f %f\n',3);
            if count == 3
                mm  = fileline(1);
                nn  = fileline(2);
                amp = fileline(3);
                Bnormal_from_TF_and_plasma_current = Bnormal_from_TF_and_plasma_current + amp*sin(mm*theta_plasma_2D + nfp*nn*zeta_plasma_2D);
            end
        end
    else
    end
    Bnormal_from_TF_and_plasma_current = Bnormal_from_TF_and_plasma_current * curpol;
end
Bnormal_from_plasma_current_1D = reshape(Bnormal_from_TF_and_plasma_current, [ntheta_plasma*nzeta_plasma,1]);
compareVariableToFortran('Bnormal_from_TF_and_plasma_current')


% ***********************************************
% Compute the basis functions and f on the (theta,zeta) grids.
% ***********************************************

switch symmetry_option
    case {1,2}
        num_basis_functions = mnmax_coil + 1;
    case {3}
        num_basis_functions = mnmax_coil * 2 + 1;
    otherwise
        error('Invalid value for symmetry_option')
end
basis_functions_R = zeros(ntheta_coil*nzeta_coil, num_basis_functions);
basis_functions_zeta_Z = zeros(ntheta_coil*nzeta_coil, num_basis_functions);

% The first basis function is the constant function:
basis_functions_R(:,1) = 1;
basis_functions_zeta_Z(:,1) = 1;

fprintf('Computing Fourier functions.\n')
tic
[zeta_coil_2D, theta_coil_2D] = meshgrid(zeta_coil,theta_coil);
zeta_coil_indices = 1:nzeta_coil;
switch symmetry_option
    case {1,2}
        % sines only
        for imn = 1:mnmax_coil
            angle = xm_coil(imn)*theta_coil_2D - xn_coil(imn)*zeta_coil_2D;
            cosangle = cos(angle);
            sinangle = sin(angle);
            basis_functions_R(     :,imn+1) = reshape(sinangle, [ntheta_coil*nzeta_coil,1]);
            basis_functions_zeta_Z(:,imn+1) = reshape(cosangle, [ntheta_coil*nzeta_coil,1]);
        end
    case {3}
        % Both sines and cosines
        for imn = 1:mnmax_coil
            angle = xm_coil(imn)*theta_coil_2D - xn_coil(imn)*zeta_coil_2D;
            cosangle = cos(angle);
            sinangle = sin(angle);
            basis_functions_R(:,imn+1)            = reshape(sinangle, [ntheta_coil*nzeta_coil,1]);
            basis_functions_R(:,imn+1+mnmax_coil) = reshape(cosangle, [ntheta_coil*nzeta_coil,1]);
        end
        basis_functions_zeta_Z = basis_functions_R;
end
fprintf('Done. Took %g sec.\n',toc)
        
%compareVariableToFortran('basis_functions')


% *********************************************
% Compute g
% *********************************************
%return

fprintf('Computing inductance\n')
tic


sinzeta = sin(zetal_coil);
coszeta = cos(zetal_coil);
zeta_plasma_indices = 1:nzeta_plasma;
inductance = zeros(ntheta_plasma*nzeta_plasma, ntheta_coil*nzeta_coil, ns_magnetization, 3);
for ks = 1:ns_magnetization
    for js = 1:ns_integration
        for itheta_coil = 1:ntheta_coil
            for izeta_coil = 1:nzeta_coil
                index_coil = (izeta_coil-1)*ntheta_coil + itheta_coil;
                factor = Jacobian_coil(itheta_coil,izeta_coil,js) * s_weights(js) * interpolate_magnetization_to_integration(js,ks);
                for l_coil = 0:(nfp-1)
                    izetal_coil = izeta_coil + l_coil*nzeta_coil;
                    
                    dx = r_plasma(1,:,zeta_plasma_indices) - (r_coil(1,itheta_coil,izetal_coil) + sign_normal * s_integration(js) * d(itheta_coil,izeta_coil) * nX(itheta_coil,izetal_coil));
                    dy = r_plasma(2,:,zeta_plasma_indices) - (r_coil(2,itheta_coil,izetal_coil) + sign_normal * s_integration(js) * d(itheta_coil,izeta_coil) * nY(itheta_coil,izetal_coil));
                    dz = r_plasma(3,:,zeta_plasma_indices) - (r_coil(3,itheta_coil,izetal_coil) + sign_normal * s_integration(js) * d(itheta_coil,izeta_coil) * nZ(itheta_coil,izetal_coil));
                    dr2 = dx.*dx + dy.*dy + dz.*dz;
                    denominator = dr2 .* sqrt(dr2);
                    
                    normal_plasma_times_dr = (dx .* normal_plasma(1,:,zeta_plasma_indices) + dy .* normal_plasma(2,:,zeta_plasma_indices) + dz .* normal_plasma(3,:,zeta_plasma_indices));
                    
                    % e_R component
                    temp = (normal_plasma(1,:,zeta_plasma_indices)*coszeta(izetal_coil) ...
                        +   normal_plasma(2,:,zeta_plasma_indices)*sinzeta(izetal_coil) ...
                        - (3./dr2) .*  normal_plasma_times_dr ...
                        .* (dx * coszeta(izetal_coil) + dy * sinzeta(izetal_coil) )) ./ denominator;
                    inductance(:,index_coil,ks,1) = inductance(:,index_coil,ks,1) - factor * reshape(temp, [ntheta_plasma*nzeta_plasma,1]);
                    
                    % e_zeta component
                    % e_zeta = cos(zeta) * e_Y - sin(zeta) * e_X
                    temp = (normal_plasma(1,:,zeta_plasma_indices)*(-sinzeta(izetal_coil)) ...
                        +   normal_plasma(2,:,zeta_plasma_indices)*coszeta(izetal_coil) ...
                        - (3./dr2) .* normal_plasma_times_dr ...
                        .* (dx * (-sinzeta(izetal_coil)) + dy * coszeta(izetal_coil) )) ./ denominator;
                    inductance(:,index_coil,ks,2) = inductance(:,index_coil,ks,2) - factor * reshape(temp, [ntheta_plasma*nzeta_plasma,1]);
                    
                    % e_Z component
                    temp = (normal_plasma(3,:,zeta_plasma_indices) ...
                        - (3./dr2) .* normal_plasma_times_dr ...
                        .* ( dz )) ./ denominator;
                    inductance(:,index_coil,ks,3) = inductance(:,index_coil,ks,3) - factor * reshape(temp, [ntheta_plasma*nzeta_plasma,1]);
                end
            end
        end
    end
end
inductance = inductance * (mu0/(4*pi));
assignin('base','inductance_m',inductance)

if symmetry_option==2
    % Test symmetry of inductance matrix
    figure(1)
    clf
    signs = [1,-1,-1];
    for j_RZetaZ = 1:3
        itheta_coil = 2;
        izeta_coil = 3;
        index_coil = (izeta_coil-1)*ntheta_coil + itheta_coil;
        
        subplot(3,3,j_RZetaZ)
        d1 = reshape(inductance(:,index_coil,1,j_RZetaZ),[ntheta_plasma,nzeta_plasma]);
        imagesc(d1)
        colorbar
        
        itheta_coil = ntheta_coil + 2 - itheta_coil;
        izeta_coil = nzeta_coil + 2 - izeta_coil;
        index_coil = (izeta_coil-1)*ntheta_coil + itheta_coil;
        
        subplot(3,3,j_RZetaZ+3)
        d2 = rot90(signs(j_RZetaZ)*reshape(inductance(:,index_coil,1,j_RZetaZ),[ntheta_plasma,nzeta_plasma]),2);
        d2=circshift(d2,[1,1]);
        imagesc(d2)
        colorbar
        
        subplot(3,3,j_RZetaZ+6)
        imagesc(d1-d2)
        colorbar
    end
    return
end


fprintf('Done. Took %g sec.\n',toc)

%compareVariableToFortran('inductance')

tic1 = tic;
g = zeros(ntheta_plasma*nzeta_plasma, num_basis_functions, ns_magnetization, 3);
for js = 1:ns_magnetization
    g(:,:,js,1) = (dtheta_coil * dzeta_coil) * inductance(:,:,js,1) * basis_functions_R;
    g(:,:,js,2) = (dtheta_coil * dzeta_coil) * inductance(:,:,js,2) * basis_functions_zeta_Z;
    g(:,:,js,3) = (dtheta_coil * dzeta_coil) * inductance(:,:,js,3) * basis_functions_zeta_Z;
end
fprintf('inductance -> g: %g\n',toc(tic1))


compareVariableToFortran('g')
%return
% *********************************************
% Compute matrices and RHS for the normal equations:
% *********************************************

norm_normal_plasma_vec = reshape(norm_normal_plasma,[ntheta_plasma*nzeta_plasma,1]);
norm_normal_coil_vec   = reshape(norm_normal_coil,  [ntheta_coil*nzeta_coil,    1]);
diag_inv_norm_normal_plasma = diag(1./norm_normal_plasma_vec);
diag_inv_norm_normal_coil   = diag(1./norm_normal_coil_vec);
system_size = num_basis_functions * ns_magnetization * 3;
Bnormal_from_TF_and_plasma_current_1D = reshape(Bnormal_from_TF_and_plasma_current, [ntheta_plasma*nzeta_plasma,1]);

tic
fprintf('Computing RHS_B.\n')
RHS_B = zeros(system_size,1);
for j = 1:3
    for js = 1:ns_magnetization
        RHS_B((1:num_basis_functions) + (j-1)*ns_magnetization*num_basis_functions + (js-1)*num_basis_functions) ...
            = (-dtheta_plasma*dzeta_plasma)*(Bnormal_from_TF_and_plasma_current_1D' * g(:,:,js,j))';
    end
end
fprintf('Done. Took %g sec.\n',toc)

compareVariableToFortran('RHS_B')

tic
fprintf('Computing matrix_B.\n')
matrix_B = zeros(system_size);
for j = 1:3
    for k = 1:3
        for js = 1:ns_magnetization
            for ks = 1:ns_magnetization
                % Note here that the col/row indices must be consistent
                % with the indices in the left and right g matrices!
                col_indices = (1:num_basis_functions) + (j-1)*ns_magnetization*num_basis_functions + (js-1)*num_basis_functions;
                row_indices = (1:num_basis_functions) + (k-1)*ns_magnetization*num_basis_functions + (ks-1)*num_basis_functions;
                matrix_B(row_indices, col_indices) = (dtheta_plasma*dzeta_plasma)*( (g(:,:,ks,k)') * diag_inv_norm_normal_plasma * g(:,:,js,j) );
            end
        end
    end
end
fprintf('Done. Took %g sec.\n',toc)
compareVariableToFortran('matrix_B')

tic
fprintf('Computing matrix_regularization.\n')
matrix_regularization = zeros(system_size);
for j_RZetaZ = 1:3
    if j_RZetaZ==1
        basis_functions = basis_functions_R;
    else
        basis_functions = basis_functions_zeta_Z;
    end
    for js = 1:ns_magnetization
        for ks = 1:ns_magnetization
            block = zeros(num_basis_functions);
            for ps = 1:ns_integration
                block = block + s_weights(ps) * interpolate_magnetization_to_integration(ps,js) * interpolate_magnetization_to_integration(ps,ks) ...
                    * basis_functions' * (diag(reshape(Jacobian_coil(:,:,ps) .* d,[ntheta_coil*nzeta_coil,1])) * basis_functions);
            end
            row_indices = (j_RZetaZ-1) * ns_magnetization*num_basis_functions + (js-1)*num_basis_functions + (1:num_basis_functions);
            col_indices = (j_RZetaZ-1) * ns_magnetization*num_basis_functions + (ks-1)*num_basis_functions + (1:num_basis_functions);
            matrix_regularization(row_indices,col_indices) = block;
        end
    end
end
matrix_regularization = matrix_regularization * dtheta_coil * dzeta_coil;

fprintf('Done. Took %g sec.\n',toc)
compareVariableToFortran('matrix_regularization')


% *********************************************
% Solve the system for each lambda:
% *********************************************

single_valued_current_potential_mn = zeros(num_basis_functions, nlambda);
single_valued_current_potential_thetazeta = zeros(ntheta_coil, nzeta_coil, nlambda);
current_potential = zeros(ntheta_coil, nzeta_coil, nlambda);
[zeta_coil_2D, theta_coil_2D] = meshgrid(zeta_coil, theta_coil);
chi2_B = zeros(nlambda,1);
chi2_M = zeros(nlambda,1);
Bnormal_total = zeros(ntheta_plasma, nzeta_plasma, nlambda);
K2 = zeros(ntheta_coil, nzeta_coil, nlambda);

for ilambda=1:nlambda
    fprintf('Solving system for lambda = %g  (%d of %d)\n',lambda(ilambda), ilambda, nlambda)
    
    tic
    matrix = matrix_B + lambda(ilambda) * matrix_regularization;
    RHS    = RHS_B;
    %RHS    = RHS_B    + lambda(ilambda) * RHS_K;
    %matrix = matrix_B - lambda(ilambda) * matrix_K;
    %RHS    = RHS_B    - lambda(ilambda) * RHS_K;
    fprintf('  Summing matrices: %g sec.\n',toc)
    
    tic
    solution = matrix \ RHS;
    fprintf('  Solve: %g sec.\n',toc)
    
    tic
    temp = zeros(ntheta_plasma*nzeta_plasma,1);
    for j = 1:3
        for js = 1:ns_magnetization
            indices = (j-1)*ns_magnetization*num_basis_functions + (js-1)*num_basis_functions + (1:num_basis_functions);
            temp = temp + g(:,:,js,j) * solution(indices);
        end
    end
    this_Bnormal = reshape(temp,[ntheta_plasma,nzeta_plasma]) ./ norm_normal_plasma + Bnormal_from_TF_and_plasma_current;
    Bnormal_total(:,:,ilambda) = this_Bnormal;
    chi2_B(ilambda) = nfp*dtheta_plasma*dzeta_plasma*sum(sum(this_Bnormal .* this_Bnormal .* norm_normal_plasma));
    %{
    single_valued_current_potential_mn(:,ilambda) = solution;
    this_single_valued_current_potential_thetazeta = reshape(basis_functions*solution, [ntheta_coil,nzeta_coil]);
    single_valued_current_potential_thetazeta(:,:,ilambda) = this_single_valued_current_potential_thetazeta;
    current_potential(:,:,ilambda) = this_single_valued_current_potential_thetazeta ...
        + zeta_coil_2D * (net_poloidal_current_Amperes/(2*pi)) ...
        + theta_coil_2D * (net_toroidal_current_Amperes/(2*pi));
    
    this_Bnormal = Bnormal_from_TF_and_plasma_current + Bnormal_from_net_coil_currents ...
        + reshape(g*solution, [ntheta_plasma,nzeta_plasma]) ./ norm_normal_plasma;
    Bnormal_total(:,:,ilambda) = this_Bnormal;
    
    K_difference_x = d_x - f_x*solution;
    K_difference_y = d_y - f_y*solution;
    K_difference_z = d_z - f_z*solution;
    this_K2_over_N = reshape(K_difference_x.*K_difference_x + K_difference_y.*K_difference_y + K_difference_z.*K_difference_z, [ntheta_coil, nzeta_coil]) ./(norm_normal_coil);
    K2(:,:,ilambda) = this_K2_over_N ./ norm_normal_coil;
    chi2_K(ilambda) = nfp*dtheta_coil*dzeta_coil*sum(sum(this_K2_over_N));
    %}
    fprintf('  Diagnostics: %g sec.\n',toc)
    fprintf('  chi2_B: %g,   chi2_M: %g\n',chi2_B(ilambda),chi2_M(ilambda))
end

compareVariableToFortran('Bnormal_total')
compareVariableToFortran('chi2_B')
compareVariableToFortran('chi2_M')

if compareToFortran
    RHS_B_m = RHS_B;
    RHS_B_f = evalin('base','RHS_B_f');
    figure(100);
    clf;
    subplot(1,2,1);
    plot(RHS_B_m,'displayname','matlab');
    hold on;
    plot(RHS_B_f,':','displayname','fortran');
    legend show
    title('RHS')
    
    subplot(1,2,2);
    plot(RHS_B_m-RHS_B_f)
    title('RHS difference')

    matrix_B_m = matrix_B;
    matrix_B_f = evalin('base','matrix_B_f');
    figure(101);
    clf;
    subplot(1,3,1);
    imagesc(matrix_B_m);
    colorbar;
    title('matrix B matlab')
    subplot(1,3,2);
    imagesc(matrix_B_f);
    colorbar;
    title('matrix B fortran')
    subplot(1,3,3);
    imagesc(matrix_B_m-matrix_B_f);
    colorbar
    title('difference')
    
    
    matrix_regularization_m = matrix_regularization;
    matrix_regularization_f = evalin('base','matrix_regularization_f');
    figure(102);
    clf;
    subplot(1,3,1);
    imagesc(matrix_regularization_m);
    colorbar;
    title('matrix regularization matlab')
    subplot(1,3,2);
    imagesc(matrix_regularization_f);
    colorbar;
    title('matrix regularization fortran')
    subplot(1,3,3);
    imagesc(matrix_regularization_m-matrix_regularization_f);
    colorbar
    title('difference')
    
    
    chi2_B_m = chi2_B;
    chi2_B_f = evalin('base','chi2_B_f');
    figure(103);
    clf;
    subplot(1,2,1);
    semilogy(chi2_B_m,'displayname','matlab');
    hold on;
    semilogy(chi2_B_f,':','displayname','fortran');
    title('chi2 B')
    legend show
    
    subplot(1,2,2);
    semilogy(abs(chi2_B_m-chi2_B_f))
    title('difference')
    
end

return


% *********************************************
% Done with the main calculation.
% Now plot results.
% *********************************************

if ~ plot_results
    return
end

figure(2)
clf
numRows=2;
numCols=3;

subplot(numRows,numCols,1)
loglog(chi2_K, chi2_B,'o-')
xlabel('chi2 K')
ylabel('chi2 B')

subplot(numRows,numCols,2)
loglog(lambda, chi2_B,'o-')
xlabel('lambda')
ylabel('chi2 B')

subplot(numRows,numCols,3)
semilogy(lambda, chi2_B,'o-')
xlabel('lambda')
ylabel('chi2 B')

subplot(numRows,numCols,5)
loglog(lambda, chi2_K,'o-')
xlabel('lambda')
ylabel('chi2 K')

subplot(numRows,numCols,6)
semilogy(lambda, chi2_K,'o-')
xlabel('lambda')
ylabel('chi2 K')

% ***********************************************************************
% Plot single-valued part of the current potential

figure(3)
clf
numContours = 25;

numPlots = min([max_nlambda_for_contour_plots,nlambda]);
ilambda_to_plot = unique(round(linspace(1,nlambda,numPlots)));
numPlots = numel(ilambda_to_plot);

numCols=ceil(sqrt(numPlots));
numRows=ceil(numPlots / numCols);

for iplot = 1:numPlots
    subplot(numRows,numCols,iplot)
    contourf(zeta_coil_2D, theta_coil_2D, single_valued_current_potential_thetazeta(:,:,ilambda_to_plot(iplot)), numContours,'EdgeColor','none')
    colorbar
    xlabel('zeta')
    ylabel('theta')
    title(['Single valued current potential for lambda=',num2str(lambda(ilambda_to_plot(iplot)))])
end
    

% ***********************************************************************
% Plot full current potential

figure(4)
clf

for iplot = 1:numPlots
    subplot(numRows,numCols,iplot)
    contourf(zeta_coil_2D, theta_coil_2D, current_potential(:,:,ilambda_to_plot(iplot)), numContours,'EdgeColor','none')
    colorbar
    xlabel('zeta')
    ylabel('theta')
    title(['Current potential for lambda=',num2str(lambda(ilambda_to_plot(iplot)))])
end
    
% ***********************************************************************
% Plot K^2

figure(5)
clf

for iplot = 1:numPlots
    subplot(numRows,numCols,iplot)
    contourf(zeta_coil_2D, theta_coil_2D, K2(:,:,ilambda_to_plot(iplot)), numContours,'EdgeColor','none')
    colorbar
    xlabel('zeta')
    ylabel('theta')
    title(['K^2 for lambda=',num2str(lambda(ilambda_to_plot(iplot)))])
end

% ***********************************************************************
% Plot B_normal

figure(6)
clf

numCols=ceil(sqrt(numPlots+2));
numRows=ceil((numPlots+2) / numCols);

subplot(numRows,numCols,1)
contourf(zeta_plasma_2D, theta_plasma_2D, Bnormal_from_TF_and_plasma_current, numContours,'EdgeColor','none')
colorbar
xlabel('zeta')
ylabel('theta')
title('Bnormal from plasma current')

subplot(numRows,numCols,2)
contourf(zeta_plasma_2D, theta_plasma_2D, Bnormal_from_net_coil_currents, numContours,'EdgeColor','none')
colorbar
xlabel('zeta')
ylabel('theta')
title('Bnormal from net coil currents')

for iplot = 1:numPlots
    subplot(numRows,numCols,iplot+2)
    contourf(zeta_plasma_2D, theta_plasma_2D, Bnormal_total(:,:,ilambda_to_plot(iplot)), numContours,'EdgeColor','none')
    colorbar
    xlabel('zeta')
    ylabel('theta')
    title(['Total Bnormal for lambda=',num2str(lambda(ilambda_to_plot(iplot)))])
end

end

%    stringForTop = ['Singular vectors v of the transfer matrix: coil surface (threshold=',num2str(pseudoinverse_thresholds(whichThreshold)),')'];
%    annotation('textbox',[0 0.96 1 .04],'HorizontalAlignment','center',...
%        'Interpreter','none','VerticalAlignment','bottom',...
%        'FontSize',11,'LineStyle','none','String',stringForTop);
