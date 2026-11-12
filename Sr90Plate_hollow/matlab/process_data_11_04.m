%% compare the 3 different positions
% columns: t_mm, x_mm, y_mm, well, dose, uncertainty  
pos0 = readmatrix("results\2026.11.04\moving_source_pos0\line_x6.5_y6.5_0deg.csv");
pos1 = readmatrix("results\2026.11.04\moving_source_pos1\line_x0_y6.5_0deg.csv");
pos2 = readmatrix("results\2026.11.04\moving_source_pos2\line_x0_y0_45deg.csv");

wellpitch = 13; % mm
figure; hold on; prettyplot;
xlabel('position along line (mm)')
ylabel('dose rate (Gy/hr)')
errorbar(pos0(:,1),pos0(:,5),pos0(:,6),'-');
errorbar(pos1(:,1),pos1(:,5),pos1(:,6),'-');
errorbar(pos2(:,1),pos2(:,5),pos2(:,6),'-');
xlim([-20 20])


%% at each position, look at different angles
pos0_angle0 = pos0;
pos0_angle45 = readmatrix("results\2026.11.04\moving_source_pos0\line_x6.5_y6.5_45deg.csv");
pos0_angle90 = readmatrix("results\2026.11.04\moving_source_pos0\line_x6.5_y6.5_90deg.csv");
figure; hold on; prettyplot;
xlabel('position along line (mm)')
ylabel('dose rate (Gy/hr)')
errorbar(pos0_angle0(:,1),pos0_angle0(:,5),pos0_angle0(:,6),'-');
errorbar(pos0_angle45(:,1),pos0_angle45(:,5),pos0_angle45(:,6),'-');
errorbar(pos0_angle90(:,1),pos0_angle90(:,5),pos0_angle90(:,6),'-');
legend('0','45','90')
xlim([-20 20])

%%
pos1_angle0 = pos1;
pos1_angle23 = readmatrix("results\2026.11.04\moving_source_pos1\line_x0_y6.5_22.5deg.csv");
pos1_angle45 = readmatrix("results\2026.11.04\moving_source_pos1\line_x0_y6.5_45deg.csv");

figure; hold on; prettyplot;
xlabel('position along line (mm)')
ylabel('dose rate (Gy/hr)')
errorbar(pos1_angle0(:,1),pos1_angle0(:,5),pos1_angle0(:,6),'-');
errorbar(pos1_angle23(:,1),pos1_angle23(:,5),pos1_angle23(:,6),'-');
errorbar(pos1_angle45(:,1),pos1_angle45(:,5),pos1_angle45(:,6),'-');
legend('0','22.5','45')
xlim([-20 20])

%%
pos2_angle23 = readmatrix("results\2026.11.04\moving_source_pos2\line_x0_y0_22.5deg.csv");
pos2_angle45 = pos2;

figure; hold on; prettyplot;
xlabel('position along line (mm)')
ylabel('dose rate (Gy/hr)')
errorbar(pos2_angle23(:,1),pos2_angle23(:,5),pos2_angle23(:,6),'-');
errorbar(pos2_angle45(:,1),pos2_angle45(:,5),pos2_angle45(:,6),'-');
legend('22.5','45')
xlim([-20 20])