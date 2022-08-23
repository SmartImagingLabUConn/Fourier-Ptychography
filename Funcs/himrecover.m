function [him,tt,fmaskpro,imseqlow]=himrecover(imseqlow,kx,ky,NA,wlength,spsize,psize,z,opts)
% FP algorithm to recover high-resolution image from low-resolution measured images
% Input:
%       imseqlow: low-res measurements, [m1 x n1 x numim] matrix
%       kx,ky: normalized wavevector of each LED, which should times k0 later
%       NA: objective NA
%       wlength: central peak wavelength of LED, in m
%       spsize: pixel size of low-res image on sample plane, in m
%       psize: pixel size of high-res image on sample plane, in m
%       z: known defocus distance, in m
%       opts: other optimization parameters
% Output:
%       him: recovered high-res image, [m x n] matrix
%       tt: recorded intensity ratio between estimated and measured low-res amplitude, used for intensity correction
%       fmaskpro: recovered pupil function
%       imseqlow: low-res amplitudes after intensity correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% input arguments check
if ~isfield(opts,'loopnum')
    opts.loopnum = 10;
end
if ~isfield(opts,'alpha')
    opts.alpha = 1;
end
if ~isfield(opts,'beta')
    opts.beta = 1;
end
if ~isfield(opts,'gamma_obj')
    opts.gamma_obj = 1;
end
if ~isfield(opts,'gamma_p')
    opts.gamma_p = 1;
end
if ~isfield(opts,'eta_obj')
    opts.eta_obj = 0;
end
if ~isfield(opts,'eta_p')
    opts.eta_p = 0;
end
if ~isfield(opts,'T')
    opts.eta_p = 0;
end
if ~isfield(opts,'aberration')
    opts.aberration = 0;
end
loopnum = opts.loopnum;
alpha = opts.alpha;
beta = opts.beta;
gamma_obj = opts.gamma_obj;
gamma_p = opts.gamma_p;
eta_obj = opts.eta_obj;
eta_p = opts.eta_p;
T = opts.T;
aberration = opts.aberration;

%% k-space parameterization
[m1, n1, numim] = size(imseqlow);
pratio = round(spsize/psize); % upsampling ratio
m = pratio*m1; n = pratio*n1;
k0 = 2*pi/wlength;
kx = k0*kx; ky = k0*ky;
NAfilx = NA*(1/wlength)*n*psize; NAfily = NA*(1/wlength)*m*psize; % m1*spsize = m*psize
kmax = pi/psize; % the max wave vector of the OTF
dkx = 2*pi/(psize*n); dky = 2*pi/(psize*m);
kx2 = -kmax:kmax/((n-1)/2):kmax; ky2 = -kmax:kmax/((m-1)/2):kmax; % odd N
[kxm, kym] = meshgrid(kx2,ky2); kzm = sqrt(k0^2-kxm.^2-kym.^2);

%% prior knowledge of aberration
H2 = exp(1j.*z.*real(kzm)).*exp(-abs(z).*abs(imag(kzm))); % define the defocus aberration if it is known or you want to test it
astigx = 0; astigy = 0; % define the astigmatism aberration if it is known or you want to test it
[M1, N1] = meshgrid(1:m1,1:n1);
zn = astigx*gzn(max(m1,n1),2*max(round(NAfily),round(NAfilx)),2,2)+...
     astigy*gzn(max(m1,n1),2*max(round(NAfily),round(NAfilx)),-2,2);
zn = imresize(zn,[m1,n1]);
if  any(aberration ~= 0,'all')
    fmaskpro = aberration; % pre-calibrated aberrations
else
    fmaskpro = 1.*double(((N1-(m1+1)/2)/NAfily).^2+((M1-(n1+1)/2)/NAfilx).^2<=1)... % low-pass filter
    .*H2(round((m+1)/2-(m1-1)/2):round((m+1)/2+(m1-1)/2),round((n+1)/2-(n1-1)/2):round((n+1)/2+(n1-1)/2))... % defocus aberration
    .*exp(pi*1j.*zn); % astigmatism aberration
    % In this example, we can test the effect of the defocus aberration (z) and astigmatism aberrations (astigx, astigy) 
    % If the aberration is unknown, one can test different z, astigx, astigy for the best result
    % Gradient descent can also be used to update z, astigx, astigy (did not implement here)
    % Higher-order Zernike modes can be tested in a similar manner
end

%% initialization
him = imresize(sum(imseqlow,3),[m,n]); 
himFT = fftshift(fft2(him));

%% main part to optimize estimate of high-res image
for i = 1:2 % 2 initial iterations to get a rough estimate of high-res image
    for i3 = 1:numim
        % when the image size is even, there will be a half pixel displacement for the cnter. 
        kxc=round((n+1)/2-kx(1,i3)/dkx);  
        kyc=round((m+1)/2-ky(1,i3)/dky);
        kyl=round(kyc-(m1-1)/2);kyh=round(kyc+(m1-1)/2);
        kxl=round(kxc-(n1-1)/2);kxh=round(kxc+(n1-1)/2);
        O_j=himFT(kyl:kyh,kxl:kxh);
        lowFT=O_j.*fmaskpro;
        im_lowFT=ifft2(ifftshift(lowFT));
        updatetemp=pratio^2.*imseqlow(:,:,i3);
        im_lowFT=updatetemp.*exp(1j.*angle(im_lowFT)); 
        lowFT_p=fftshift(fft2(im_lowFT));
        himFT(kyl:kyh,kxl:kxh)=himFT(kyl:kyh,kxl:kxh)+conj(fmaskpro)./(max(max((abs(fmaskpro)).^2))).*(lowFT_p-lowFT);
    end
end

countimg = 0;
tt = ones(1,loopnum*numim);

% for momentum method
vobj0 = zeros(m,n);
vp0 = zeros(m1, n1);
ObjT = himFT; 
PT = fmaskpro;

for i = 1:loopnum
    for i3 = 1:numim 
        countimg=countimg+1; 
        kxc=round((n+1)/2-kx(1,i3)/dkx);  
        kyc=round((m+1)/2-ky(1,i3)/dky);
        kyl=round(kyc-(m1-1)/2);kyh=round(kyc+(m1-1)/2);
        kxl=round(kxc-(n1-1)/2);kxh=round(kxc+(n1-1)/2);
        O_j=himFT(kyl:kyh,kxl:kxh); 
        lowFT=O_j.*fmaskpro;
        im_lowFT=ifft2(ifftshift(lowFT));
        tt(1,i3+(i-1)*numim)=(mean(mean(abs(im_lowFT)))/mean(mean(pratio^2*abs(imseqlow(:,:,i3))))); % LED intensity correctioin 
        if i>2
            imseqlow(:,:,i3)=imseqlow(:,:,i3).*tt(1,i3+(i-1)*numim); 
        end     

        updatetemp=pratio^2.*imseqlow(:,:,i3);
        im_lowFT=updatetemp.*exp(1j.*angle(im_lowFT)); 
        lowFT_p=fftshift(fft2(im_lowFT));

        himFT(kyl:kyh,kxl:kxh)=himFT(kyl:kyh,kxl:kxh)+...
        gamma_obj.*conj(fmaskpro).*((lowFT_p-lowFT))./((1-alpha).*abs(fmaskpro).^2 + alpha.*max(max(abs(fmaskpro).^2)));
        fmaskpro=fmaskpro+gamma_p.*conj(O_j).*((lowFT_p-lowFT))./((1-beta).*abs(O_j).^2 + beta.*max(max(abs(O_j).^2)));            

        if countimg == T % momentum method
            vobj = eta_obj.*vobj0 + (himFT - ObjT);
            himFT = ObjT + vobj;
            vobj0 = vobj;                  
            ObjT = himFT;

            vp = eta_p.*vp0 + (fmaskpro - PT);
            fmaskpro = PT + vp;
            vp0 = vp;
            PT = fmaskpro;

            countimg = 0;
        end
    end
end

him=ifft2(ifftshift(himFT));

end