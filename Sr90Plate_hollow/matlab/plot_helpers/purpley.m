function C = purpley(N)

% # Edit this gradient at https://eltos.github.io/gradient/#0:000000-19.3:2A0079-47.3:404AB3-100:B2DFFF

positions = [0 0.193 0.473 1.000];
colors = [0.000, 0.000, 0.000; 0.165, 0.000, 0.475; 0.251, 0.290, 0.702; 0.698, 0.875, 1.000];

C = customcolormap(positions,colors,N);

end
