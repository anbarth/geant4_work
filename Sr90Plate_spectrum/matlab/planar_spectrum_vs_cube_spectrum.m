spectrum_plane = readmatrix("results\spectrum_above_well_floor.csv");
spectrum_cube = readmatrix("results\11.9.2026\spectrum_0.5_0.5.csv");

spectra = {spectrum_plane,spectrum_cube};

figure; hold on; prettyplot;
xlabel('Energy (MeV)')
ylabel('Normalized counts')
ax1=gca; ax1.YScale='log';
for ii=1:length(spectra)
    spectrum = spectra{ii};

    e_low = spectrum(:,1);
    e_high = spectrum(:,2);
    % for histogram plotting
    %energy_bin_edges = [e_low; e_high(end)];
    e_mid = spectrum(:,3);

    counts = spectrum(:,4);
    Ctot = sum(counts);

    plot(e_mid, counts/Ctot,'-');

end