clc;
clear all;

data.c = 299792458;
data.fc = 77e9;
data.B = 1e9;
data.PRF = 2000;
data.PRI = 1/data.PRF;
data.VPC_pos0 = [0,0];    % original VPC position x coordinate and y coordinate

% % define FoV used to test angular ambiguities 
% data.target = [15,7,1;12,10,1;18,13,1];    % 3 targets position
% data.FoV_min = [10,5];   % FoV closest point
% data.FoV_max = [20,15];  % FoV longest point

% define near field FoV
% data.target = [7,3,1;12,8,1;17,10,1];    % 3 targets position
% data.FoV_min = [5,0];   % FoV closest point
% data.FoV_max = [20,15];  % FoV longest point

% % define far field FoV
data.target = [37,33,5;42,38,5;47,40,5];   % 3 targets position
data.FoV_min = [35,30];   % FoV closest point
data.FoV_max = [50,45];  % FoV longest point

data.vego = [5,0]; % radar velocity
data.vtarget = [0,0]; % target velocity
data.N_pulse = 100;
data.Nch = 8;
data.dy = data.c/(4*data.fc);   % MIMO VPC spacing
% Forward projection
data = signalGen(data);
% Define Range-Doppler map
range_doppler_maps = cell(1, data.Nch);
for ch = 1:data.Nch
    range_compressed_data = data.sRC(:,:,ch);
    doppler_processed_data = fft(range_compressed_data, [], 1);
    doppler_shifted_data = fftshift(doppler_processed_data, 1);
    range_doppler_maps{ch} = abs(doppler_shifted_data);
end
combined_map = mean(cat(3, range_doppler_maps{:}), 3);

% doppler frequency axis
fdmax = data.PRF/2;
CPI = data.PRI*data.N_pulse;
delta_fd = 1/CPI;
fd = -fdmax:delta_fd:fdmax-delta_fd;
v_doppler = fd*data.c/(2*data.fc);
range_bin = linspace(round(mean(mean(data.Rmin))),round(mean(mean(data.Rmax))),size(combined_map,2));

% Display Range-Doppler map of phase history data
figure(1)
imagesc(fd, range_bin, fliplr(transpose(combined_map)));
% imagesc(v_doppler, range_bin, transpose(combined_map));
title_str = ['Vego: ', num2str(data.vego(1)), ' m/s,  Vtarget: ', num2str(data.vtarget(1)), ' m/s'];
title(title_str);
xlabel('doppler shift (Hz)')
% yticks([0 25 50 75 100 125 150 175 200 225 250]);
% yticklabels({'5','7','9','11','13','15','17','19','21','23','25'})
ylabel('range (m)')
colorbar;
axis xy;  % This makes sure the first row of the matrix is at the bottom
% colormap('gray');


% Back-projection process
% Scene size and pixel spacing
pixel_spacing = 0.02; % meters
% Define the pixel grid for the image
x_vec = data.FoV_min(1):pixel_spacing:data.FoV_max(1);
y_vec = data.FoV_min(2):pixel_spacing:data.FoV_max(2);
[data.x_mat, data.y_mat] = meshgrid(x_vec, y_vec);
% Operating back-projection function
data = BP(data);
image = data.image;
figure(2)
mesh(abs(image))
xlim([0 750])
ylim([0 750])
title_str = ['Vego: ', num2str(data.vego(1)), ' m/s,  Vtarget: ', num2str(data.vtarget(1)), ' m/s'];
title(title_str);
% xticks([0 50 100 150 200 250 300 350 400 450 500]);
% xticklabels({'10','11','12','13','14','15','16','17','18','19','20'})
% yticks([0 50 100 150 200 250 300 350 400 450 500]);
% yticklabels({'5','6','7','8','9','10','11','12','13','14','15'})
xticks([0 50 100 150 200 250 300 350 400 450 500 550 600 650 700 750]);
xticklabels({'35','36','37','38','39','40','41','42','43','44','45','46','47','48','49','50'})
yticks([0 50 100 150 200 250 300 350 400 450 500 550 600 650 700 750]);
yticklabels({'30','31','32','33','34','35','36','37','38','39','40','41','42','43','44','45'})
xlabel('x coordinate (m)')
ylabel('y coordinate (m)')
% colormap('gray');
