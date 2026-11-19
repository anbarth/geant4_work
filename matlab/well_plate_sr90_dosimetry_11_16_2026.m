%% Film-water correction factor -- open configuration
% e_low_MeV,e_high_MeV,e_mid_MeV,counts_per_decay,rel_error,counts_per_decay_per_MeV
spectrum_open = readmatrix("Sr90Film_spectrum\results\11.13.2026\spectrum.csv");

% plot the spectrum
e_mid  = spectrum_open(:,3);
counts = spectrum_open(:,4);
Ctot   = sum(counts);
figure; hold on; prettyplot; ax1=gca; ax1.YScale='log';
xlabel('Energy (MeV)'); ylabel('Normalized counts')
plot(e_mid, counts/Ctot,'-');

% calculate stopping powers
S_water_spectrum = read_stopping_powers("water_stopping_power.txt",e_mid,counts);
S_water = sum(counts .* S_water_spectrum) / Ctot;   % MeV cm2/g
S_ebtxd_spectrum = read_stopping_powers("ebt-xd_active_layer_stopping_power.txt",e_mid,counts);
S_ebtxd = sum(counts .* S_ebtxd_spectrum) / Ctot;   % MeV cm2/g

disp(['Open configuration - S_water: ' num2str(S_water) ' MeV cm2/g'])
disp(['Open configuration - S_film: ' num2str(S_ebtxd) ' MeV cm2/g'])
disp(['S_film/S_water: ' num2str(S_ebtxd/S_water)]);

%% Film-water correction factor -- well configuration
spectrum_plate_1 = readmatrix("Sr90Plate_spectrum\results\2026.11.13\moving_source_pos0\spectrum.csv");
spectrum_plate_2 = readmatrix("Sr90Plate_spectrum\results\2026.11.13\moving_source_pos1\spectrum.csv");
spectrum_plate_3 = readmatrix("Sr90Plate_spectrum\results\2026.11.13\moving_source_pos2\spectrum.csv");

% plot the spectra -- they are essentially identical
% some very small differences at lowest energies
spectra = {spectrum_plate_1,spectrum_plate_2,spectrum_plate_3};
figure; hold on; prettyplot;
xlabel('Energy (MeV)'); ylabel('Normalized counts');
ax1=gca; ax1.YScale='log';
for ii=1:length(spectra)
    spectrum = spectra{ii};
    e_mid = spectrum(:,3);
    counts = spectrum(:,4);
    Ctot = sum(counts);
    plot(e_mid, counts/Ctot,'-');
end

% average the counts and work with the average spectrum
e_mid = spectrum_plate_1(:,3);
spectrum_counts_avg = (spectrum_plate_1(:,4)+spectrum_plate_2(:,4)+spectrum_plate_3(:,4))/3;

% calculate stopping powers
S_water_spectrum = read_stopping_powers("water_stopping_power.txt",e_mid,counts);
S_water = sum(counts .* S_water_spectrum) / Ctot;   % MeV cm2/g
S_ebtxd_spectrum = read_stopping_powers("ebt-xd_active_layer_stopping_power.txt",e_mid,counts);
S_ebtxd = sum(counts .* S_ebtxd_spectrum) / Ctot;   % MeV cm2/g

disp(['Well configuration - S_water: ' num2str(S_water) ' MeV cm2/g'])
disp(['Well configuration - S_film: ' num2str(S_ebtxd) ' MeV cm2/g'])
disp(['S_film/S_water: ' num2str(S_ebtxd/S_water)]);

%% Well geometry correction factor
% t_mm,x_mm,y_mm,well,dose_Gy_per_h,uncertainty_Gy_per_h
spectrum_profile_1 = readmatrix("Sr90Plate_spectrum\results\2026.11.13\moving_source_pos0\line_x6.45_y6.45_0deg.csv");
spectrum_profile_2 = readmatrix("Sr90Plate_spectrum\results\2026.11.13\moving_source_pos1\line_x0_y6.45_0deg.csv");
spectrum_profile_3 = readmatrix("Sr90Plate_spectrum\results\2026.11.13\moving_source_pos2\line_x0_y0_45deg.csv");
% r_mm,dose_Gy_per_h,uncertainty_Gy_per_h,n_voxels
water_profile = readmatrix("Sr90Water\results\11.13.2026\radial_profile.csv");

convfac = 100/(60*60); % Gy/hr -> cGy/s
%convfac=1;
%kwell = 0.45;
kwell=1;

figure; hold on; prettyplot;
xlabel("Distance from source center (mm)")
ylabel("Dose rate (cGy/s)")
errorbar(spectrum_profile_1(:,1),spectrum_profile_1(:,5)*convfac,spectrum_profile_1(:,6)*convfac);
errorbar(spectrum_profile_2(:,1),spectrum_profile_2(:,5)*convfac,spectrum_profile_2(:,6)*convfac);
errorbar(spectrum_profile_3(:,1),spectrum_profile_3(:,5)*convfac,spectrum_profile_3(:,6)*convfac);
errorbar(water_profile(:,1),water_profile(:,2)*convfac*kwell,water_profile(:,3)*convfac*kwell,'k');
xlim([0 15])

