%% compare the 3 different positions
pos1 = readmatrix("results\2026.10.29\cutoff_dist_2.5um\voxprofile_C_dose_rate.csv");
pos2 = readmatrix("results\2026.10.29\source_c4_c5_quarter\voxprofile_C_dose_rate.csv");
pos3 = readmatrix("results\2026.10.29\source_c4_c5\voxprofile_C_dose_rate.csv");

wellpitch = 13; % mm
figure; hold on; prettyplot;
xlabel('position along row C (mm)')
ylabel('dose rate (Gy/hr)')
errorbar(pos1(:,1),pos1(:,2),pos1(:,3),'-o');
errorbar(pos2(:,1)+wellpitch/4,pos2(:,2),pos2(:,3),'-o');
errorbar(pos3(:,1)+wellpitch/2,pos3(:,2),pos3(:,3),'-o');
legend('under C5','25% between C5 and C4','50% between C5 and C4')

%% compare the different energy cutoffs
cutoff2p5 = pos1;
cutoff5 = readmatrix("results\2026.10.29\cutoff_dist_5um\voxprofile_C_dose_rate.csv");
cutoff10 = readmatrix("results\2026.10.29\cutoff_dist_10um\voxprofile_C_dose_rate.csv");
figure; hold on; prettyplot;
xlabel('position along row C (mm)');
ylabel('dose rate (Gy/hr)');
errorbar(cutoff2p5(:,1),cutoff2p5(:,2),cutoff2p5(:,3),'-o');
errorbar(cutoff5(:,1),cutoff5(:,2),cutoff5(:,3),'-o');
errorbar(cutoff10(:,1),cutoff10(:,2),cutoff10(:,3),'-o');
legend('2.5 \mum','5 \mum','10 \mum');

%% compare the different physics lists
standard_opt4 = pos1;
penelope = readmatrix("results\2026.10.29\penelope\voxprofile_C_dose_rate.csv");
figure; hold on; prettyplot;
xlabel('position along row C (mm)');
ylabel('dose rate (Gy/hr)');
errorbar(standard_opt4(:,1),standard_opt4(:,2),standard_opt4(:,3),'-o');
errorbar(penelope(:,1),penelope(:,2),penelope(:,3),'-o');
legend('standard (opt 4)','penelope');

%% look at different cuts: fully centered
full_cut0 =  readmatrix("results\2026.10.29\cutoff_dist_2.5um\line_C5_0deg.csv");
full_cut23 =  readmatrix("results\2026.10.29\cutoff_dist_2.5um\line_C5_22.5deg.csv");
full_cut45 =  readmatrix("results\2026.10.29\cutoff_dist_2.5um\line_C5_45deg.csv");
full_cut90 =  readmatrix("results\2026.10.29\cutoff_dist_2.5um\line_C5_90deg.csv");
figure; hold on; prettyplot;
xlabel('distance from well center (mm)');
ylabel('dose rate (Gy/hr)');
errorbar(full_cut0(:,1),full_cut0(:,5),full_cut0(:,6),'-o');
errorbar(full_cut23(:,1),full_cut23(:,5),full_cut23(:,6),'-o');
errorbar(full_cut45(:,1),full_cut45(:,5),full_cut45(:,6),'-o');
errorbar(full_cut90(:,1),full_cut90(:,5),full_cut90(:,6),'-o');

%% look at different cuts: half centered
half_cut0 =  readmatrix("results\2026.10.29\source_c4_c5\line_x0_y6.5_0deg.csv");
half_cut23 =  readmatrix("results\2026.10.29\source_c4_c5\line_x0_y6.5_22.5deg.csv");
half_cut45 =  readmatrix("results\2026.10.29\source_c4_c5\line_x0_y6.5_45deg.csv");
figure; hold on; prettyplot;
xlabel('distance from well center (mm)');
ylabel('dose rate (Gy/hr)');
errorbar(half_cut0(:,1),half_cut0(:,5),half_cut0(:,6),'-o');
errorbar(half_cut23(:,1),half_cut23(:,5),half_cut23(:,6),'-o');
errorbar(half_cut45(:,1),half_cut45(:,5),half_cut45(:,6),'-o');