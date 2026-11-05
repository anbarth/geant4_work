function C = warm_tones_to_blk(N)
% Edit this gradient at https://eltos.github.io/gradient/#0:FF00A9-17.3:A4004E-42.7:7E0010-69.9:5E0F00-100:000000

mymap = [[0.000, [1.000, 0.000, 0.663]],
    [0.173, [0.643, 0.000, 0.306]],
    [0.427, [0.494, 0.000, 0.063]],
    [0.699, [0.369, 0.059, 0.000]],
    [1.000, [0.000, 0.000, 0.000]]];

positions = mymap(:,1);
colors = mymap(:,2:end);

C = customcolormap(positions,colors,N);

end