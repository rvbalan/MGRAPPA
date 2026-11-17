function [ChannelImage, Image, N] = setup_dataset(shrink_brain,profile_smooth_level,slice_offset,load_dicom)

% image_file_names = {'dicomImage.mat',...
% '20201022_for_Michael\MPRageDicom3dicomImage.mat',   ...
% '20201022_for_Michael\MPRageDicom8dicomImage.mat',  ...
% '20201022_for_Michael\MPRageDicom4dicomImage.mat', ...
% '20201022_for_Michael\MPRageDicom9dicomImage.mat',...
% '20201022_for_Michael\MPRageDicom10dicomImage.mat',...
% '20201022_for_Michael\MPRageDicom5dicomImage.mat', ...
% '20201022_for_Michael\MPRageDicom1dicomImage.mat',...
% '20201022_for_Michael\MPRageDicom6dicomImage.mat',...
% '20201022_for_Michael\MPRageDicom2dicomImage.mat',...
% '20201022_for_Michael\MPRageDicom7dicomImage.mat'};
% 
% image_file_names = {'dicomImage.mat'};


if 0==1% Brain 
    
    load('ChannelImage_m0_pose0');
    load('ChannelImage_m0_pose1');
    load('ChannelImage_m1');

    N = size(ChannelImage_m0_pose0,1);
    channel_len = size(ChannelImage_m0_pose0,3);

elseif true
    % if load dicom:
    if load_dicom
        image_file_name = 'dicomImage.mat';
        dicomImage=0;
        image_file_name
        load(image_file_name, 'dicomImage');
        Image = squeeze(dicomImage(:,floor(size(dicomImage,2)/2)+slice_offset,:));
    else
        % If loading *.nii file:
        image_file_name = 'sub-02_ses-01_T1w.nii';
        info = niftiinfo(image_file_name);
        volume = niftiread(info);
        slice = squeeze(volume(:,floor(size(volume,2)/2)+slice_offset,:));
        slice = double(slice);
        mmm = min(min(slice));
        mmx = max(max(slice));
        Image = (slice - mmm)/(mmx - mmm); % Normalize to [0,1]
%        Image = slice/max(abs(mmx),abs(mmm)); % Normalize to [-1,1]
%        without changing its sign distribution --> bad idea!

        % Rotate if needed
        Image = rot90(Image);
    end

    figure
    imshowSc(Image)
    Image = imresize(Image, [size(Image,1), size(Image,1)], 'bilinear' );
    exportgraphics(gcf,sprintf('ScanImageStretch%i.jpg',slice_offset),'Resolution',300)

    figure
    imshowSc(Image)
    exportgraphics(gcf,sprintf('ScanImage%i.jpg',slice_offset),'Resolution',300)

%     Image = fft_smooth(Image, profile_smooth_level);
    
%     figure
%     mesh(abs(Image))
%     figure
%     mesh(angle(Image))
%     Image
% assert(false)

    N = size(Image,1);
    Image2 = imresize(Image, .75, 'bilinear' ); % squeeze the image to have space for rotations and translations
    Image = zeros(N,N);
    Image(N/2-size(Image2,1)/2:N/2+size(Image2,1)/2-1, N/2-size(Image2,2)/2:N/2+size(Image2,2)/2-1) = Image2;
    
%     Image = Image(size(Image,1)/2-N/2:size(Image,1)/2+N/2-1 , size(Image,1)/2-N/2:size(Image,1)/2+N/2-1);
    
%     Image = squeeze(dicomImage(:,:,floor(size(dicomImage,3)/2)));

    N = size(Image,1);

    if shrink_brain % shrink brain 
        % Down Sample
        Image_sig = fftshift(fft2(Image));
        Image = ifft2(ifftshift(Image_sig(N/2-N/4:N/2+N/4-1, N/2-N/4:N/2+N/4-1)));
        N = size(Image,1);

%         % Down Sample
%         Image_sig = fftshift(fft2(Image));
%         Image = ifft2(ifftshift(Image_sig(N/2-N/4:N/2+N/4-1, N/2-N/4:N/2+N/4-1)));
%         N = size(Image,1);

%         % Pad image
%         Image(N:2*N, N:2*N) = 0;
%         N = size(Image,1);
%         Image = Image([3*N/4:N 1:3*N/4-1], [3*N/4:N 1:3*N/4-1]);

        % Pad image
        Image(N:2*N, N:2*N) = 0;
        N = size(Image,1);
        Image = Image([3*N/4:N 1:3*N/4-1], [3*N/4:N 1:3*N/4-1]);
    else
%         % Pad image
%         Image(N:2*N, N:2*N) = 0;
%         N = size(Image,1);
%         Image = Image([3*N/4:N 1:3*N/4-1], [3*N/4:N 1:3*N/4-1]);
% 
%         % Down Sample
%         Image_sig = fftshift(fft2(Image));
%         Image = ifft2(ifftshift(Image_sig(N/2-N/4:N/2+N/4-1, N/2-N/4:N/2+N/4-1)));
%         N = size(Image,1);
        
    end
    
    
%     % Down Sample
%     Image_sig = fftshift(fft2(Image));
%     Image = ifft2(ifftshift(Image_sig(N/2-N/4:N/2+N/4-1, N/2-N/4:N/2+N/4-1)));
%     N = size(Image,1);

%     figure;imshowSc(Image);

    [ChannelImage, channel_len] = getChannels(Image,0,0,0,0,profile_smooth_level);
    
elseif 0==1 % Shepp-Logan
    N=2*64;
    Image = phantom('Modified Shepp-Logan',N);   
        
    % Pad image
    Image(N:2*N, N:2*N) = 0;
    N = size(Image,1);
    Image = Image([3*N/4:N 1:3*N/4-1], [3*N/4:N 1:3*N/4-1]);
    
    % Down Sample
    Image_sig = fftshift(fft2(Image));
    Image = ifft2(ifftshift(Image_sig(N/2-N/4:N/2+N/4-1, N/2-N/4:N/2+N/4-1)));
    N = size(Image,1);
    
    [ChannelImage, channel_len] = getChannels(Image,0,0);
    
else
    load('scan_sig_data.mat','ChannelImage','channel_len','N');

    ChannelImage = center_data(ChannelImage,N,channel_len);
    
    %Change data, put in notch 
    for icha = 1:channel_len
        ChannelImage(1:32,:,icha)=0;
        ChannelImage(64:70,:,icha)=0;
        ChannelImage(:,64:70,icha)=0;
        for i=1:128-6-1
            ChannelImage(i:i+6,i:i+6,icha)=0;
            ChannelImage(i:i+6,128-i-6:128-i,icha)=0;
        end
    end

    %get signal data
    s = zeros(N,N,channel_len);
    parfor icha = 1:channel_len
        s(:,:,icha)=fftshift(fft2(ChannelImage(:,:,icha)));
    end

    % Change N
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    s2 = zeros(64,64,20);
    parfor i = 1:channel_len
        s2(:,:,i) = s(32+1:32+64,32+1:32+64,i);
    end
    s = s2;
    N = 64;
    assert(N == 64);

    ChannelImage = zeros(size(s));
    parfor icha = 1:channel_len
        ChannelImage(:,:,icha)=ifft2(ifftshift(s(:,:,icha)));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DONE Change N

end
