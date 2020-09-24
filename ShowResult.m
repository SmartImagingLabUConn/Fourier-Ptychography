close all
clear
clc
addpath('Results')
%% blood
edge = 60;
figure;
load('bloodsmear_red_result.mat');
imBlood = zeros(size(him,1),size(him,1),3);
subplot(221),imshow(angle(him(edge+1:end-edge,edge+1:end-edge)),[-1 -1+pi]);title('red');
imBlood(:,:,1)=abs(him)/max(abs(him(:)))*1.1;
load('bloodsmear_green_result.mat');
subplot(222),imshow(angle(him(edge+1:end-edge,edge+1:end-edge)),[]);title('green')
imBlood(:,:,2)=abs(him)/max(abs(him(:)))*1.1;
load('bloodsmear_blue_result.mat');
subplot(223),imshow(angle(him(edge+1:end-edge,edge+1:end-edge)),[]);title('blue');
imBlood(:,:,3)=abs(him)/max(abs(him(:)))*1.05;
subplot(224),imshow(imBlood(edge+1:end-edge,edge+1:end-edge,:).^(1/0.75),[]);title('color');