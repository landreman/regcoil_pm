regcoil_out_filename = 'examples/compareToMatlab1/regcoil_out.compareToMatlab1.nc';
%regcoil_out_filename = 'examples/NCSX_vv_randomResolution1_iterate_d/regcoil_out.NCSX_vv_randomResolution1_iterate_d.nc';
%regcoil_out_filename = '/Users/mattland/Box Sync/work19/20190526-01-testing_regcoil_pm/20190526-01-044-vv_thetaZeta64_mpolNtor12_sMagnetization2_sIntegration2_d0.15/regcoil_out.NCSX.nc';

ilambda = 12;

nfp = ncread(regcoil_out_filename,'nfp');
sign_normal = double(ncread(regcoil_out_filename,'sign_normal'));
norm_normal_coil = ncread(regcoil_out_filename,'norm_normal_coil');
abs_M = ncread(regcoil_out_filename,'abs_M');
magnetization_vector = ncread(regcoil_out_filename,'magnetization_vector');
zetal_coil = ncread(regcoil_out_filename,'zetal_coil');
ns_magnetization = ncread(regcoil_out_filename,'ns_magnetization');
s_magnetization = ncread(regcoil_out_filename,'s_magnetization');
try
    normal_coil = ncread(regcoil_out_filename,'normal_coil');
catch
    error('Unable to read normal_coil from the output file. Probably save_level was set >1.')
end
r_plasma = ncread(regcoil_out_filename,'r_plasma');
r_coil = ncread(regcoil_out_filename,'r_coil');
nzetal_coil = ncread(regcoil_out_filename,'nzetal_coil');
ntheta_coil = ncread(regcoil_out_filename,'ntheta_coil');
d = ncread(regcoil_out_filename,'d');

norm_normal_coil = repmat(norm_normal_coil,[1,nfp]);
unit_normal_coil = normal_coil;
for j=1:3
    unit_normal_coil(j,:,:) = squeeze(unit_normal_coil(j,:,:)) ./ norm_normal_coil;
end
size(d)

r_coil_outer = r_coil;
for j = 1:3
    r_coil_outer(j,:,:) = squeeze(r_coil_outer(j,:,:)) + sign_normal * repmat(d(:,:,ilambda),[1,nfp]) .* squeeze(unit_normal_coil(j,:,:));
end

size(abs_M)

abs_M_inner = repmat(abs_M(:,:,  1,ilambda),[1,nfp]);
abs_M_outer = repmat(abs_M(:,:,end,ilambda),[1,nfp]);

figure(1)
clf

% Close plasma surface in theta and zeta:
r_plasma(:,end+1,:) = r_plasma(:,1,:);
r_plasma(:,:,end+1) = r_plasma(:,:,1);

surf(squeeze(r_plasma(1,:,:)), squeeze(r_plasma(2,:,:)), squeeze(r_plasma(3,:,:)),'facecolor','r','edgecolor','none')
hold on
light
daspect([1,1,1])
axis vis3d off
set(gca,'clipping','off')
rotate3d on
zoom(1.6)

big_d = repmat(d(:,:,ilambda),[1,nfp]);

% Close surface in theta:
r_coil(:,end+1,:) = r_coil(:,1,:);
r_coil_outer(:,end+1,:) = r_coil_outer(:,1,:);
abs_M_inner(end+1,:) = abs_M_inner(1,:);
abs_M_outer(end+1,:) = abs_M_outer(1,:);

% Only show a sector of the magnetization region:
max_zeta_index = round(nzetal_coil/2);
r_coil = r_coil(:,:,1:max_zeta_index);
r_coil_outer = r_coil_outer(:,:,1:max_zeta_index);
unit_normal_coil = unit_normal_coil(:,:,1:max_zeta_index);
abs_M_inner = abs_M_inner(:,1:max_zeta_index);
abs_M_outer = abs_M_outer(:,1:max_zeta_index);
big_d = big_d(:,1:max_zeta_index);

surf(squeeze(r_coil(1,:,:)), squeeze(r_coil(2,:,:)), squeeze(r_coil(3,:,:)),abs_M_inner,'edgecolor','none','facecolor','interp','facealpha',1)
surf(squeeze(r_coil_outer(1,:,:)), squeeze(r_coil_outer(2,:,:)), squeeze(r_coil_outer(3,:,:)),abs_M_outer,'edgecolor','none','facecolor','interp','facealpha',0.7)
colorbar
lighting gouraud

magnetization_vector = repmat(magnetization_vector(:,:,:,:,ilambda),[1,nfp,1,1]);
magnetization_vector = magnetization_vector(:,1:max_zeta_index,:,:);
MZ = magnetization_vector(:,:,:,3);
MX = MZ*0;
MY = MZ*0;
for izeta = 1:max_zeta_index
    coszeta = cos(zetal_coil(izeta));
    sinzeta = sin(zetal_coil(izeta));
    MX(:,izeta,:) = magnetization_vector(:,izeta,:,1) * coszeta + magnetization_vector(:,izeta,:,2) * (-sinzeta);
    MY(:,izeta,:) = magnetization_vector(:,izeta,:,1) * sinzeta + magnetization_vector(:,izeta,:,2) * coszeta;
end

X = MX * 0;
Y = MX * 0;
Z = MX * 0;
for js = 1:ns_magnetization
    X(:,:,js) = squeeze(r_coil(1,1:end-1,:)) + sign_normal * s_magnetization(js) * big_d .* squeeze(unit_normal_coil(1,:,:));
    Y(:,:,js) = squeeze(r_coil(2,1:end-1,:)) + sign_normal * s_magnetization(js) * big_d .* squeeze(unit_normal_coil(2,:,:));
    Z(:,:,js) = squeeze(r_coil(3,1:end-1,:)) + sign_normal * s_magnetization(js) * big_d .* squeeze(unit_normal_coil(3,:,:));
end
scale = 2;
quiver3(X,Y,Z,MX,MY,MZ,scale,'k')
