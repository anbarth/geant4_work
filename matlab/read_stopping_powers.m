% reads stopping powers from the given ESTAR file
% returns the stopping powers at the given energies (interpolates as
% needed)
% if necessary. we extrapolate and give a warning

function S_bin = read_stopping_powers(fname,energies,counts)


% Columns: kinetic energy [MeV], collision, radiative, total stopping power [MeV cm2/g]
% (the first 8 lines of the file are header text)
sp = readmatrix(fname, 'NumHeaderLines', 8);
E_tab   = sp(:,1);
S_coll  = sp(:,2);
S_rad   = sp(:,3);
S_tot   = sp(:,4);

% Choose which one to average: S_tot (total) or S_coll (collision, = what is deposited locally)
S_tab = S_tot;

% Interpolate onto the spectrum bins. Log-log linear interpolation suits stopping power well.
% Bins below the table's lowest energy (10 keV) are extrapolated -- see the warning below.
S_bin = exp(interp1(log(E_tab), log(S_tab), log(energies), 'linear', 'extrap'));

nBelow = energies < E_tab(1);
if any(nBelow)
    warning('%d bins (%.2g%% of the counts) lie below %.0f keV, the lowest ESTAR energy; their stopping power is extrapolated.', ...
        nnz(nBelow), 100*sum(counts(nBelow))/sum(counts), 1000*E_tab(1));
end

end
