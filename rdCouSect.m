% read the coulomb stress change result on a cross section and convert it
% to a compatible format

clear; close all;
% fid = fopen('G:\Coulomb3\coulomb34\coulomb34\output_files\coustress_on_p2_section.cou');
fid = fopen('G:\Coulomb3\coulomb34\coulomb34\output_files\p2_section.cou');
data = textscan(fid, '%f %f %f %f %f %f', 'HeaderLines', 3);
x = data{1,1};
y = data{1,2};
z = data{1,3};
stress = data{1,4};
nrow = 71;
ncol = length(x)/nrow;

for i = 1: ncol
    for j = 1: nrow
        dist(j, i) = sqrt( (x((i-1)*nrow+j)-x(j)).^2+(y((i-1)*nrow+j)-y(j)).^2 );
    end
end
dis = [];
for i = 1: ncol
    dis = [dis; dist(:, i)];
end
fdata=[dis z stress];
% save('G:\Coulomb3\coulomb34\coulomb34\cou2gmt_output\coustress_p2_section.dat', 'fdata','-ascii');
save('G:\Coulomb3\coulomb34\coulomb34\cou2gmt_output\p2_section.dat', 'fdata','-ascii');


fid = fopen('G:\Coulomb3\coulomb34\coulomb34\output_files\coustress_on_p2_section_depth.cou');
data = textscan(fid, '%f %f %f %f %f %f', 'HeaderLines', 3);
x = data{1,1};
y = data{1,2};
z = data{1,3};
stress = data{1,4};
nrow = 71;
ncol = length(x)/nrow;

for i = 1: ncol
    for j = 1: nrow
        dist(j, i) = sqrt( (x((i-1)*nrow+j)-x(j)).^2+(y((i-1)*nrow+j)-y(j)).^2 );
    end
end
dis = [];
for i = 1: ncol
    dis = [dis; dist(:, i)];
end
fdata=[dis z stress];
save('G:\Coulomb3\coulomb34\coulomb34\cou2gmt_output\coustress_p2_section_depth.dat', 'fdata','-ascii');
