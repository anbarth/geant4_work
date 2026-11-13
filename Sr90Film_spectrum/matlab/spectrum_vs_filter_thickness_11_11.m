spectrum_1 = readmatrix("results\11.11.2026\steel_filter_50um\spectrum.csv");
spectrum_2 = readmatrix("results\11.11.2026\steel_filter_100um\spectrum.csv");

spectra = {spectrum_1,spectrum_2};

figure; hold on; prettyplot;
xlabel('Energy (MeV)')
ylabel('Normalized counts')
%ylabel('Counts per decay')
%ax1=gca; ax1.YScale='log';
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