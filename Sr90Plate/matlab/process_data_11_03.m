%% compare the 3 different positions
% columns: t_mm, x_mm, y_mm, well, dose, uncertainty  
pos1 = readmatrix("results\2026.11.02\moving_source_pos0\line_x0_y6.5_0deg.csv");
pos2 = readmatrix("results\2026.11.02\moving_source_pos0.25\line_x0_y6.5_0deg.csv");
pos3 = readmatrix("results\2026.11.02\moving_source_pos0.5\line_x0_y6.5_0deg.csv");

wellpitch = 13; % mm
figure; hold on; prettyplot;
xlabel('position along row C (mm)')
ylabel('dose rate (Gy/hr)')
errorbar(pos1(:,1),pos1(:,5),pos1(:,6),'-');
errorbar(pos2(:,1)-wellpitch/4,pos2(:,5),pos2(:,6),'-');
errorbar(pos3(:,1)-wellpitch/2,pos3(:,5),pos3(:,6),'-');


%% at each position, look at different angles
pos1_angle0 = pos1;
pos1_angle45 = readmatrix("results\2026.11.02\moving_source_pos0\line_x0_y6.5_45deg.csv");
figure; hold on; prettyplot;
xlabel('position along row C (mm)')
ylabel('dose rate (Gy/hr)')
errorbar(pos1_angle0(:,1),pos1_angle0(:,5),pos1_angle0(:,6),'-');
errorbar(pos1_angle45(:,1),pos1_angle45(:,5),pos1_angle45(:,6),'-');
legend('0','45','90')
xlim([-20 20])

%%
pos2_angle0 = pos2;
pos2_angle45 = readmatrix("results\2026.11.02\moving_source_pos0.25\line_x3.25_y6.5_45deg.csv");
pos2_angle90 = readmatrix("results\2026.11.02\moving_source_pos0.25\line_x3.25_y6.5_90deg.csv");
figure; hold on; prettyplot;
xlabel('position along row C (mm)')
ylabel('dose rate (Gy/hr)')
errorbar(pos2_angle0(:,1)-3.25,pos2_angle0(:,5),pos2_angle0(:,6),'-');
errorbar(pos2_angle45(:,1),pos2_angle45(:,5),pos2_angle45(:,6),'-');
errorbar(pos2_angle90(:,1),pos2_angle90(:,5),pos2_angle90(:,6),'-');
legend('0','45','90')
xlim([-20 20])

%%
pos3_angle0 = pos3;
pos3_angle45 = readmatrix("results\2026.11.02\moving_source_pos0.5\line_x6.5_y6.5_45deg.csv");
pos3_angle90 = readmatrix("results\2026.11.02\moving_source_pos0.5\line_x6.5_y6.5_90deg.csv");
figure; hold on; prettyplot;
xlabel('position along row C (mm)')
ylabel('dose rate (Gy/hr)')
errorbar(pos3_angle0(:,1)-6.5,pos3_angle0(:,5),pos3_angle0(:,6),'-');
errorbar(pos3_angle45(:,1),pos3_angle45(:,5),pos3_angle45(:,6),'-');
errorbar(pos3_angle90(:,1),pos3_angle90(:,5),pos3_angle90(:,6),'-');
legend('0','45','90')
xlim([-20 20])