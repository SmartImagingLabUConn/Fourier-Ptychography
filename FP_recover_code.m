close all;clear;clc;

%% Prepare the experimental data
% Add necessary folders into the current working directory
addpath(genpath(pwd));
% Load data file
data_name = 'MouseKidney_green';
data_dir = ['Data\' data_name '.mat']; 
load(data_dir); % refer to 'data_description.txt' for more details
% Display raw images
is_show = 'center'; % 'center' shows the first low-res raw image; 'all' dynamically shows all low-res images
if strcmp(is_show,'center')
    figure(1);
    set(gcf,'outerposition',get(0,'ScreenSize'))
    imshow(imlow_HDR(:,:,1),[]);
    title(['raw image ' num2str(1)]);
elseif strcmp(is_show,'all')
    for k = 1:size(imlow_HDR,3)
        figure(1);
        set(gcf,'outerposition',get(0,'ScreenSize'))
        imshow(imlow_HDR(:,:,k),[]);
        title(['raw image ' num2str(k)]); pause(0.1);
    end
end

%% Set up the experiment parameters
xstart = 18; ystart = 20; % absolute coordinate of initial LED
arraysize = 15; % side length of lit LED array
[xlocation, ylocation] = LED_location(xstart, ystart, arraysize);
H      = 90.88; % distance between LEDs and sample, in mm
LEDp   = 4;     % distance between adjacent LEDs, in mm
nglass = 1.52;  % refraction index of glass substrate
t      = 1;     % glass thickness, in mm
[kx, ky, NAt] = k_vector(xlocation-xstart, ylocation-ystart, H, LEDp, nglass, t, theta, xint, yint, arraysize^2);

%% Reconstruct by FP algorithm
NA          = 0.1;      % objective NA
spsize      = 1.845e-6; % pixel size of low-res image on sample plane, in m
upsmp_ratio = 4;        % upsampling ratio
psize       = spsize/upsmp_ratio; % pixel size of high-res image on sample plane, in m

opts.loopnum    = 10;   % iteration number
opts.alpha      = 1;    % '1' for ePIE, other value for rPIE
opts.beta       = 1;    % '1' for ePIE, other value for rPIE
opts.gamma_obj  = 1;    % the step size for object updating
opts.gamma_p    = 1;    % the step size for pupil updating
opts.eta_obj    = 0.2;  % the step size for adding momentum to object updating
opts.eta_p      = 0.2;  % the step size for adding momentum to pupil updating
opts.T          = 1;    % do momentum every T images. '0' for no momentum during the recovery; integer, generally (0, arraysize^2].
opts.aberration = aberration; % pre-calibrated aberration, if available

used_idx = 1:1:arraysize^2; % choose which raw image is used, for example, 1:2:arraysize^2 means do FPM recovery with No1 image, No3 image, No5 image......
imlow_used = imlow_HDR(:,:,used_idx);
kx_used = kx(used_idx);
ky_used = ky(used_idx);
[him, tt, fprobe, imlow_HDR1] = himrecover(imlow_used, kx_used, ky_used, NA, wlength, spsize, psize, z, opts);

figure;
set(gcf,'outerposition',get(0,'ScreenSize'))
subplot(121);imshow(abs(him(50:end-50,50:end-50)),[]);title('Amplitude');
subplot(122);imshow(angle(him(50:end-50,50:end-50)),[]);title('Phase');
disp(['Wavelength: ',num2str(wlength.*1e+9),' nm, Loop: ',num2str(opts.loopnum)]);
disp(['Maximum illumination NA = ',num2str(max(NAt(used_idx)))]);

%% Save the results
out_dir = 'Results';
mkdir(out_dir); addpath(out_dir);
out_name = [data_name '_result.mat'];
save([out_dir,'\',out_name],'him','fprobe','tt','imlow_HDR1'); 