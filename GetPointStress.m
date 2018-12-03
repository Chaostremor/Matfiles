% 
lat = 9.029: 0.12: 21.029;
lon = -99.807: 0.12: -87.807;
data = load('G:\Coulomb3\coulomb34\coulomb34\cou2gmt_output\p2_horizontal.dat');
nlon = length(lon);
nlat= length(lat);
stress = zeros(nlat, nlon);
k=1;
for i = 1: nlon
    for j =1: nlat
    stress(j, i) = data(k, 3);
    k=k+1;
    end
end
[LAT, LON]=meshgrid(lat, lon);
B = [15.585, -94.369];
BS=interp2(LAT, LON, stress', B(1), B(2), 'spline');
aaa=interp2(LAT, LON, stress', 16, -94.25, 'spline');