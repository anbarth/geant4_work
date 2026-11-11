spectrum_1 = readmatrix("results\spectrum_exiting_steel.csv");
spectrum_2 = readmatrix("results\spectrum_above_well_floor.csv");
spectrum_3 = readmatrix("results\spectrum_2.5mm_above_well_floor.csv");

spectra = {spectrum_1,spectrum_2,spectrum_3};

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