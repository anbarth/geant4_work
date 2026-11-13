spectrum_1 = readmatrix("results\11.9.2026\spectrum_0.5_0.5.csv");
spectrum_2 = readmatrix("results\11.10.2026\0_0_extra_histories\spectrum.csv");
spectrum_3 = readmatrix("results\11.10.2026\0_0.5\spectrum.csv");

spectra = {spectrum_1,spectrum_2,spectrum_3};

figure; hold on; prettyplot;
xlabel('Energy (MeV)')
%ylabel('Counts per decay')
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

    %plot(e_mid, counts,'-');
    plot(e_mid, counts/Ctot,'-');

    disp(sum(e_mid.*counts) / Ctot);

end