% Compute GRAPPA weights W with and without profile rotation
% This code performs MGRAPPA for several theta,x,y shifts and computes SNR.
% The SNR is computed w.r.t. a masked region -- brain
%
% Radu Balan, 8 March 2024 , based on Michael Rawson's code;
% Second change: 29 March 2024. Third change: 1 May 2024; 
% Fourth change: May 28, 2024; June 6, 2024.
% Fifth change: January 28, 2026; March 3, 2026
% Sixth change: October 22, 2026

close all
clear all

tic

% Data file
load_dicom = false;

if load_dicom
    calibration_offset = 24;
    partition_mode_squared = true;
    profile_smooth_level = 7/8;
    extrapolate_mode = 2;
    extrapolate_thresh = .1;
    
    topband_to_remove = 50;
    bottomband_to_remove = 35;
    pixelrange_x = 1+topband_to_remove:256-bottomband_to_remove;
    pixelrange_y = 1:160;

else
    % load *.nii file
    % image_file_name = 'sub-02_ses-01_T1w.nii'; -> inside setup_dataset
    calibration_offset = 24;
    partition_mode_squared = true;
    profile_smooth_level = 7/8;
    extrapolate_mode = 2;
    extrapolate_thresh = .1;

    topband_to_remove = 60; % 
    bottomband_to_remove = 50; %
    pixelrange_x = 1+topband_to_remove:256-bottomband_to_remove;
    pixelrange_y = 1:160; %

end

acceleration_factor = 2;

currdir = pwd;

%%%%%%%% Experiments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Vary theta
%x_eq_y = false;
%%theta_deg_range = [-30:5:30];
%%theta_deg_range = [-30:5:-5,-4:1:4,5:5:30];
%theta_deg_range = [-30:1:30];
%x_offset_range = 0;
%y_offset_range = 0;
%plotcomposite = true;

%%% Vary X
%x_eq_y = false;
%theta_deg_range = 0;
%%x_offset_range = [-30:5:30];
%%x_offset_range = [-30:5:-5,-4:1:4,5:5:30];
%x_offset_range = [-30:1:30];
%y_offset_range = 0;
%plotcomposite = true;

%%% Vary Y
x_eq_y = false;
theta_deg_range = 0;
x_offset_range = 0;
%y_offset_range = [-30:5:30];
%y_offset_range = [-30:5:-5,-4:1:4,5:5:30];
y_offset_range = [-30:1:30];
plotcomposite = true;

%%% Vary X and Y together 
%theta_deg_range = 9;
%x_offset_range = 5;
%x_eq_y = false;
%y_offset_range = 1;
%plotcomposite = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lf_thresh=.05;

beta = 1;
alpha = inf;

imrot_algo = 'bilinear'; % {'nearest','bicubic','bilinear'}

shrink_brain = false;

%slice_offsets = [-40, -30, -20, -10, 0, 10, 20, 30, 40];
%slice_offsets = [-40:40]; % Maximum number of slices -- dicom file

slice_offsets = [-68:12]; % For sub-02_ses-01_T1w.nii
%slice_offsets = [-1:1]; % For sub-02_ses-01_T1w.nii

err_l2 = zeros(size(x_offset_range,2)*size(y_offset_range,2),length(slice_offsets));
err_l1 = zeros(size(x_offset_range,2)*size(y_offset_range,2),length(slice_offsets));
err_l2_corrected = zeros(size(x_offset_range,2)*size(y_offset_range,2),length(slice_offsets));
err_l1_corrected = zeros(size(x_offset_range,2)*size(y_offset_range,2),length(slice_offsets));

for slice_offset_ind = 1 : length(slice_offsets)
slice_offset = slice_offsets(slice_offset_ind)

[ChannelImage, Image, N] = setup_dataset(shrink_brain,profile_smooth_level,slice_offset,load_dicom);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
resultsfoldername = 'results';

if x_eq_y
    plot_range = x_offset_range;
    plot_xlabel = 'XY translation';
    plot_ylabel = '||f-f_{recon}||_2';
elseif size(x_offset_range,2) > 1
    plot_range = x_offset_range;
    plot_xlabel = 'X translation (voxel)';
    plot_ylabel = '||f-f_{recon}||_2';
    filename_string = 'Xaxis';
elseif size(y_offset_range,2) > 1
    plot_range = y_offset_range;
    plot_xlabel = 'Y translation (voxel)';
    plot_ylabel = '||f-f_{recon}||_2';
    filename_string = 'Yaxis';
elseif size(theta_deg_range,2) > 1
    plot_range = theta_deg_range;
    plot_xlabel = '\theta rotation (degree)';
    plot_ylabel = '||f-f_{recon}||_2';
    filename_string = 'theta';
else
    plot_range = theta_deg_range;
    plot_xlabel = '\theta rotation (degree)';
    plot_ylabel = '||f-f_{recon}||_2';
    filename_string = 'theta';
end

err_ind = 1;

for theta_ind = 1:size(theta_deg_range,2) % rotation angle
    theta_deg = theta_deg_range(theta_ind);

    for x_offset_ind = 1:size(x_offset_range,2) % translation
        x_offset = x_offset_range(x_offset_ind);

        for y_offset_ind = 1:size(y_offset_range,2) % translation
            y_offset = y_offset_range(y_offset_ind);

            if x_eq_y
                y_offset = x_offset;
            end

            [ChannelImage_rot, channel_len] = getChannels(Image,0,theta_deg,x_offset,y_offset,profile_smooth_level);

            [rot_image, grappa_rot_image, grappa_rot_image_fixed_calib] ...
                = mgrappa(ChannelImage, ChannelImage_rot, partition_mode_squared, ...
                        calibration_offset, Lf_thresh, theta_deg, x_offset, ...
                        y_offset, imrot_algo, extrapolate_thresh, extrapolate_mode, beta, alpha, acceleration_factor);
% rot_image: target image that needs estimated
% grappa_rot_image : the image returned by grappa algorithm
% grappa_rot_image_fixed_calib: the image returned by mgrappa
            pmask = 2; % mask is (2pmax+1) x (2pmax+1) around detected brain voxel
            threshold = 0.1; % 0.0: no mask; 0.1: mask --for brain voxel decision
            image_mask = FindMask(rot_image,pmask,threshold); % binary mask that identifies the brain region



            grappa_rot_image_fixed_calib_resid = rot_image - grappa_rot_image_fixed_calib;

            grappa_rot_image_resid = rot_image - grappa_rot_image;
            
%            if abs(theta_deg - 9) < 10^-10 || abs(x_offset - 5) < 10^-10 ...
%                    || abs(y_offset - 1) < 10^-10
            if abs(theta_deg - 20) < 10^-10 || abs(x_offset - 20) < 10^-10 ...
                    || abs(y_offset - 20) < 10^-10
                rot_image_plot = rot_image;

                grappa_rot_image_fixed_calib_resid_plot = abs(grappa_rot_image_fixed_calib_resid)...
                        /max(max(abs(rot_image)));
                grappa_rot_image_fixed_calib_plot = grappa_rot_image_fixed_calib;

                grappa_rot_image_resid_plot = abs(grappa_rot_image_resid)/max(max(abs(rot_image)));
                grappa_rot_image_plot = grappa_rot_image;

            end



            err_l2(err_ind,slice_offset_ind) = vnorm(image_mask.*(rot_image-grappa_rot_image),2);
            err_l1(err_ind,slice_offset_ind) = vnorm(image_mask.*(rot_image-grappa_rot_image),1);

            err_l2_corrected(err_ind,slice_offset_ind) = vnorm(image_mask.*(rot_image-grappa_rot_image_fixed_calib),2);
            err_l1_corrected(err_ind,slice_offset_ind) = vnorm(image_mask.*(rot_image-grappa_rot_image_fixed_calib),1);
            err_ind = err_ind + 1;
        end
    end
end

if plotcomposite

figure
    subplot(2,3,1);
    rot_image_plot_stretch = imresize(rot_image_plot,[256,160]);
    imshowSc(rot_image_plot_stretch(pixelrange_x,pixelrange_y));
    if size(theta_deg_range,2) > 1
        title('Rotated Coil Profiles');
    else
        title('Translated Coil Profiles');
    end

    ax = gca;
    ax.FontSize = 8;

    subplot(2,3,2);
    grappa_rot_image_plot = abs(grappa_rot_image_plot)/max(max(abs(grappa_rot_image_plot)));
    grappa_rot_image_plot_stretch = imresize(grappa_rot_image_plot, [256,160]);
    imshow(sqrt(3*grappa_rot_image_plot_stretch(pixelrange_x,pixelrange_y)),[0,1]);
    title(sprintf('GRAPPA'));

    ax = gca;
    ax.FontSize = 8;

    subplot(2,3,3);
    grappa_rot_image_fixed_calib_plot = abs(grappa_rot_image_fixed_calib_plot)...
                        /max(max(abs(grappa_rot_image_fixed_calib_plot)));
    grappa_rot_image_fixed_calib_plot_stretch = imresize(grappa_rot_image_fixed_calib_plot, [256,160]);
    imshow(sqrt(3*grappa_rot_image_fixed_calib_plot_stretch(pixelrange_x,pixelrange_y)),[0,1]);
    title(sprintf('MGRAPPA'));

    ax = gca;
    ax.FontSize = 8;


    subplot(2,3,4);
    image_mask_stretch = imresize(image_mask, [256,160]);
    imshowSc(image_mask_stretch(pixelrange_x,pixelrange_y));
    title('Binary Mask');
    ax = gca;
    ax.FontSize = 8;
    
    subplot(2,3,5);
    grappa_rot_image_resid_plot_stretch = imresize(grappa_rot_image_resid_plot,[256,160]);
    imshow((grappa_rot_image_resid_plot_stretch(pixelrange_x,pixelrange_y)),...
        [0,max( ...
                max(max((grappa_rot_image_fixed_calib_resid_plot))),...
                max(max((grappa_rot_image_resid_plot))))]);
    title(sprintf('Residual with GRAPPA'));

    ax = gca;
    ax.FontSize = 8;

    subplot(2,3,6);
    grappa_rot_image_fixed_calib_resid_plot_strectch = imresize(grappa_rot_image_fixed_calib_resid_plot,[256,160]);
    imshow((grappa_rot_image_fixed_calib_resid_plot_strectch(pixelrange_x,pixelrange_y)),...
        [0,max( ...
                max(max((grappa_rot_image_fixed_calib_resid_plot))),...
                max(max((grappa_rot_image_resid_plot))))]);
    title(sprintf('Residual with MGRAPPA'));

    ax = gca;
    ax.FontSize = 8;
cd(resultsfoldername);
    exportgraphics(gcf,sprintf('%sresiduals%i.jpg',filename_string,slice_offset),'Resolution',300)
cd(currdir);

% Set of 3 pictures
figure
subplot(1,3,1);
rot_image_plot_stretch = imresize(rot_image_plot,[256,160]);
imshowSc(rot_image_plot_stretch(pixelrange_x,pixelrange_y));
%if size(theta_deg_range,2) > 1
%    title('Rotated Coil Profiles');
%else
%    title('Translated Coil Profiles');
%end

ax = gca;
ax.FontSize = 8;

subplot(1,3,2);
grappa_rot_image_resid_plot_stretch = imresize(grappa_rot_image_resid_plot,[256,160]);
imshow(grappa_rot_image_resid_plot_stretch(pixelrange_x,pixelrange_y));
% Scaled Images:
%imshow((grappa_rot_image_resid_plot_stretch(pixelrange_x,pixelrange_y)),...
%    [0,max( ...
%    max(max((grappa_rot_image_fixed_calib_resid_plot))),...
%    max(max((grappa_rot_image_resid_plot))))]);
%title(sprintf('Residual with GRAPPA'));

ax = gca;
ax.FontSize = 8;

subplot(1,3,3);
grappa_rot_image_fixed_calib_resid_plot_strectch = imresize(grappa_rot_image_fixed_calib_resid_plot,[256,160]);
imshow(grappa_rot_image_fixed_calib_resid_plot_strectch(pixelrange_x,pixelrange_y));
% Scaled Image
%imshow((grappa_rot_image_fixed_calib_resid_plot_strectch(pixelrange_x,pixelrange_y)),...
%    [0,max( ...
%    max(max((grappa_rot_image_fixed_calib_resid_plot))),...
%    max(max((grappa_rot_image_resid_plot))))]);
%title(sprintf('Residual with MGRAPPA'));

ax = gca;
ax.FontSize = 8;
cd(resultsfoldername);
exportgraphics(gcf,sprintf('%sresiduals%i_composite.jpg',filename_string,slice_offset),'Resolution',300)
cd(currdir);

% Uncropped
figure
    subplot(2,3,1);
    imshowSc(rot_image_plot);
    if size(theta_deg_range,2) > 1
        title('Rotated Coil Profiles');
    else
        title('Translated Coil Profiles');
    end

    ax = gca;
    ax.FontSize = 8;

    subplot(2,3,2);
    grappa_rot_image_plot = abs(grappa_rot_image_plot)/max(max(abs(grappa_rot_image_plot)));
    imshow((grappa_rot_image_plot),[0,1]);
    title(sprintf('GRAPPA'));

    ax = gca;
    ax.FontSize = 8;

    subplot(2,3,3);
    grappa_rot_image_fixed_calib_plot = abs(grappa_rot_image_fixed_calib_plot)...
                        /max(max(abs(grappa_rot_image_fixed_calib_plot)));
    imshow((grappa_rot_image_fixed_calib_plot),[0,1]);
    title(sprintf('MGRAPPA'));

    ax = gca;
    ax.FontSize = 8;
    
    subplot(2,3,5);
    imshow((grappa_rot_image_resid_plot),...
        [0,max( ...
                max(max((grappa_rot_image_fixed_calib_resid_plot))),...
                max(max((grappa_rot_image_resid_plot))))]);
    title(sprintf('Residual with GRAPPA'));

    ax = gca;
    ax.FontSize = 8;

    subplot(2,3,6);
    imshow((grappa_rot_image_fixed_calib_resid_plot),...
        [0,max( ...
                max(max((grappa_rot_image_fixed_calib_resid_plot))),...
                max(max((grappa_rot_image_resid_plot))))]);
    title(sprintf('Residual with MGRAPPA'));

    ax = gca;
    ax.FontSize = 8;

    cd(resultsfoldername);
    exportgraphics(gcf,sprintf('%sresiduals_uncropped%i.jpg',filename_string,slice_offset),'Resolution',300)
    cd(currdir);

    % set of 3 plots
    figure
    immax = max(max(max(abs(rot_image_plot))))/2;
    immaxres = immax/5;
    subplot(1,3,1);
    imshow(rot_image_plot,[0,immax]);
%    if size(theta_deg_range,2) > 1
%        title('Rotated Coil Profiles');
%    else
%        title('Translated Coil Profiles');
%    end

    ax = gca;
    ax.FontSize = 8;

    subplot(1,3,2);
    imshow(grappa_rot_image_resid_plot,[0,immaxres]);
    % Scaled Image
    %imshow((grappa_rot_image_resid_plot),...
    %    [0,max( ...
    %    max(max((grappa_rot_image_fixed_calib_resid_plot))),...
    %    max(max((grappa_rot_image_resid_plot))))]);
    %    title(sprintf('Residual with GRAPPA'));

    ax = gca;
    ax.FontSize = 8;

    subplot(1,3,3);
    imshow(grappa_rot_image_fixed_calib_resid_plot,[0,immaxres]);  
    % Scaled Image
    %imshow((grappa_rot_image_fixed_calib_resid_plot),...
    %    [0,max( ...
    %    max(max((grappa_rot_image_fixed_calib_resid_plot))),...
    %    max(max((grappa_rot_image_resid_plot))))]);
%    title(sprintf('Residual with MGRAPPA'));

    ax = gca;
    ax.FontSize = 8;

    cd(resultsfoldername);
    exportgraphics(gcf,sprintf('%sresiduals_uncropped%i_composite_NewNormalization.jpg',filename_string,slice_offset),'Resolution',600)
    cd(currdir);

% L2 errors
figure

%     subplot(1,2,1)
    plot(plot_range,err_l2_corrected(:,slice_offset_ind),...
        plot_range,err_l2(:,slice_offset_ind),'LineWidth',2);
    legend({'MGrappa','Grappa'});
    title(sprintf('L2 error'));
    xlabel(plot_xlabel);
    if (strcmp(filename_string , 'theta'))
        xticks([-30 -25 -20 -15 -10 -5 0 5 10 15 20 25 30]) % Same range for theta and xoffset
    else
        xticks([-30 -25 -20 -15 -10 -5 0 5 10 15 20 25 30])
    end
    ax = gca;
    ax.FontSize = 20;

    cd(resultsfoldername);
    exportgraphics(gcf,sprintf('%sL2error%i.jpg',filename_string,slice_offset),'Resolution',300)
    cd(currdir);

end

end

%%% Bar plot
figure

errorbar(plot_range,mean(err_l2_corrected,2),std(err_l2_corrected,1,2),'b');
hold on
errorbar(plot_range,mean(err_l2,2),std(err_l2,1,2),'r');
legend({'MGrappa','Grappa'},'Location','north');
%title(sprintf('L2 error'));
xlabel(plot_xlabel);
ylabel(plot_ylabel);
if (strcmp(filename_string , 'theta'))
    xticks([-30 -25 -20 -15 -10 -5 0 5 10 15 20 25 30])
else
    xticks([-30 -25 -20 -15 -10 -5 0 5 10 15 20 25 30])
end

ax = gca;
ax.FontSize = 20;

cd(resultsfoldername);
exportgraphics(gcf,sprintf('%sL2error_errorbars.jpg',filename_string),...
    'Resolution',600)
cd(currdir);

% Same figure, no bar plot
figure

plot(plot_range,mean(err_l2_corrected,2),'LineWidth',2);
hold on
plot(plot_range,mean(err_l2,2),'LineWidth',2);
grid on
grid minor
legend({'MGrappa','Grappa'},'Location','north');
%title(sprintf('L2 error'));
xlabel(plot_xlabel);
ylabel(plot_ylabel);
if (strcmp(filename_string , 'theta'))
    xticks([-30 -25 -20 -15 -10 -5 0 5 10 15 20 25 30])
else
    xticks([-30 -25 -20 -15 -10 -5 0 5 10 15 20 25 30])
end

ax = gca;
ax.FontSize = 20;

cd(resultsfoldername);
exportgraphics(gcf,sprintf('%sL2error_error_nobars.jpg',filename_string),...
    'Resolution',600)
cd(currdir);

%%%%% Save results, other than figures
cd(resultsfoldername);
save(sprintf('%sL2error.mat',filename_string),'plot_range','err_l2_corrected','err_l2')

writematrix(err_l2_corrected,sprintf('%s_err_l2_corrected.csv',filename_string))
writematrix(mean(err_l2_corrected,2),sprintf('%s_err_l2_corrected_mean.csv',filename_string))
writematrix(std(err_l2_corrected,1,2),sprintf('%s_err_l2_corrected_std.csv',filename_string))

writematrix(err_l2,sprintf('%s_err_l2.csv',filename_string))
writematrix(mean(err_l2,2),sprintf('%s_err_l2_mean.csv',filename_string))
writematrix(std(err_l2,1,2),sprintf('%s_err_l2_std.csv',filename_string))

writematrix(plot_range,sprintf('%s_plot_range.csv',filename_string))
cd(currdir);


toc

