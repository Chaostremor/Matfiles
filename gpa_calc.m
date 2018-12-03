% calculate the undergraduate GPA
% C. Song, 2017.11.22

clear;
load('C:\Users\Song Chao\Documents\MATLAB\myself\undergraduate_all.mat');
% fid = fopen('C:\Users\Song Chao\Documents\MATLAB\myself\1.txt');
% ss = textscan(fid, '%s');
% %
% fclose(fid) ;
% sort = char(ss{1});     % 台站名是第一列
grade_point = zeros(length(grade), 1);
for i = 1: length(grade)
    if grade(i) >=90 && grade(i) <=100
        grade_point(i) = 4.0;
    elseif grade(i) >=85 && grade(i) <90
        grade_point(i) = 3.7;
    elseif grade(i) >=82 && grade(i) <85
        grade_point(i) = 3.3;    
    elseif grade(i) >=78 && grade(i) <82
        grade_point(i) = 3.0;
    elseif grade(i) >=75 && grade(i) <78
        grade_point(i) = 2.7;
    elseif grade(i) >=72 && grade(i) <75
        grade_point(i) = 2.3;
    elseif grade(i) >=68 && grade(i) <72
        grade_point(i) = 2.0;
    elseif grade(i) >=64 && grade(i) <68
        grade_point(i) = 1.5;
    elseif grade(i) >=60 && grade(i) <64
        grade_point(i) = 1.0;
    else
        grade_point(i) = 0;
    end
end

credit_hour = credit.* grade_point;

% cumulative gpa
gpa_cumu = sum(credit_hour)/ sum(credit);

% advanced gpa, all courses completed after the second year
ind = 44;
gpa_ad = sum(credit_hour(1: ind))/ sum(credit(1: ind));

% major gpa, GPA for courses in the major field of study
index = find(sortt=='R');
gpa_major = sum(credit_hour(index))/ sum(credit(index));


