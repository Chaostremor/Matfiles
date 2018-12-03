% test regional maxium 

% A = 10*ones(10,10);
% A(2:4,2:4) = 22;    % maxima 12 higher than surrounding pixels
% A(6:8,6:8) = 33;    % maxima 23 higher than surrounding pixels
% A(2,7) = 44;
% A(3,8) = 45;
% A(4,9) = 44;
% regmax = imregionalmax(A, 8);
% maxval = [];
% [maxval(:, 1), maxval(:, 2)] = find(regmax);            %1st col = row indice of regional max. of data, 2nd col = col indice
% dim = length(maxval(:, 1));
% %maxval = zeros(dim,1);
% for i = 1:dim
%     maxval(i, 3) = A(maxval(i, 1), maxval(i, 2));         %3rd col = max. value of regional max. of data
% end
% 
% [locmax, indices] = localmax(A, [], false);
% %[irow, icol] = find(locmax);
% %[ar, ac] = ind2sub([10, 10], indices);
% maxval = [];
% [maxval(:, 1), maxval(:, 2)] = find(locmax);            %1st col = row indice of regional max. of data, 2nd col = col indice
% dim = length(maxval(:, 1));
% %maxval = zeros(dim,1);
% for i = 1:dim
%     maxval(i, 3) = A(maxval(i, 1), maxval(i, 2));         %3rd col = max. value of regional max. of data
% end


regmax = imregionalmax(pwssum1);
gmax = max(pwssum1(:));
[irow, icol] = find(regmax);            % 1st col = row indice of regional max. of data, 2nd col = col indice
dim = length(irow);
%maxval = zeros(dim,1);
maxset = [];
for i = 1:dim
    maxval = pwssum1(irow(i), icol(i));     % 3rd col = max. value of regional max. of data
    if (maxval >= 0.01*gmax)
        maxt = tuse1(irow(i));
        maxp = vel(icol(i));
        maxset = [maxset; irow(i) icol(i) maxt maxp maxval];
    end
end
figure
imagesc(tuse1(:, 1), vel, pwssum1'); hold on
colormap(jet);
colorbar;
xlabel('time');
ylabel('velocity');
plot(maxset(:, 3), maxset(:, 4), 'k.', 'MarkerSize', 8); hold on
% rectx = [85, 85, 90, 90, 85];
% recty = [0.017, -0.036, -0.036, 0.017, 0.017];
% plot(rectx, recty, 'linewidth', 1.5, 'color', 'w'); hold on
[~, indice] = max(pwssum1(:));
[row, col] = ind2sub([ntw1, nvel], indice);
plot(tuse1(row), vel(col), 'k*', 'MarkerSize', 10); 
[~, indice] = min(pwssum1(:));
[row, col] = ind2sub([ntw1, nvel], indice);
plot(tuse1(row), vel(col), 'kx', 'MarkerSize', 10);  


regmax = imregionalmax(pwssum2);
gmax = max(pwssum2(:));
[irow, icol] = find(regmax);            % 1st col = row indice of regional max. of data, 2nd col = col indice
dim = length(irow);
%maxval = zeros(dim,1);
maxset = [];
for i = 1:dim
    maxval = pwssum2(irow(i), icol(i));     % 3rd col = max. value of regional max. of data
    if (maxval >= 0.1*gmax)
        maxt = tuse2(irow(i));
        maxp = rayp(icol(i));
        maxset = [maxset; irow(i) icol(i) maxt maxp maxval];
    end
           
end
figure
imagesc(tuse2(:, 1), rayp, pwssum2'); hold on
colormap(jet);
colorbar;
xlabel('time');
ylabel('ray parameter');
plot(maxset(:, 3), maxset(:, 4), 'k.', 'MarkerSize', 8); hold on
% rectx = [85, 85, 90, 90, 85];
% recty = [0.017, -0.036, -0.036, 0.017, 0.017];
% plot(rectx, recty, 'linewidth', 1.5, 'color', 'w'); hold on
[~, indice] = max(pwssum2(:));
[row, col] = ind2sub([ntw2, nrayp], indice);
plot(tuse2(row), rayp(col), 'k*', 'MarkerSize', 10); 
[~, indice] = min(pwssum2(:));
[row, col] = ind2sub([ntw2, nrayp], indice);
plot(tuse2(row), rayp(col), 'kx', 'MarkerSize', 10); 



