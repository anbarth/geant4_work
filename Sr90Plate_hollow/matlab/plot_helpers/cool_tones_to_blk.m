function C = cool_tones_to_blk(N)

%  Edit this gradient at https://eltos.github.io/gradient/#0:E9FF00-28.6:4AA400-57.7:007E6B-84.5:002A5E-100:000000

mymap =     [[0.000, [0.914, 1.000, 0.000]],
    [0.189, [0.290, 0.643, 0.000]],
    [0.434, [0.000, 0.494, 0.420]],
    [0.691, [0.000, 0.165, 0.369]],
    [1.000, [0.000, 0.000, 0.000]]];

positions = mymap(:,1);
colors = mymap(:,2:end);

C = customcolormap(positions,colors,N);

end


