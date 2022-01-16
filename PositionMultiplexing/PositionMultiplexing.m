%% Main.m
clear; close all; clc
%% load real positive images
Xsingle1 = imread('Message.png');
Xsingle1 = im2double(Xsingle1);
Xsingle2 = imread('Optical_256.png');
Xsingle2 = im2double(Xsingle2);
%% Parameters
[yN,xN]=size(Xsingle1);
lambda = 532e-6; % wavelength [mm]
%% Distances for position multiplexing
d1=100;
d2=160;
d3=180;
%% Construct random phase masks
RandomPhaseTEM1=exp(1j*2*pi*rand(yN,xN)); % Random phase mask 1
Diagonal1= RandomPhaseTEM1;
RandomPhaseTEM2=exp(1j*2*pi*rand(yN,xN)); % Random phase mask 2
Diagonal2= RandomPhaseTEM2;
X1 = Xsingle1;
L=5;
%% Y1=FRTλ{Y11}z in equation (31)
Y11=X1.*Diagonal1; % Y11=fn(a,b)*exp[i2πφ(a,b)] for n=1 in equation (31)
Y1=propTF(Y11,L,lambda,d1); % Y1=FRTλ{Y11}z in equation (31)
Y12=Y1.*Diagonal2; % Y12=Y1*exp[i2πψ(a',b')] in equation (31)
Y2=propTF(Y12,L,lambda,d2); % gn(a0,b0)=FRTλ{Y12}dn for n=1 in equation (31)
figure, imshow(Y2); % Show the first encrypted image
%% gn(a0,b0)=FRTλ{FRTλ{fn(a,b)exp[i2πφ(a,b)]}zexp[i2πψ(a',b')]}dn for n=2 equation (31)
X2 = Xsingle2;
Y21=X2.*Diagonal1; % Y21=fn(a,b)*exp[i2πφ(a,b)] for n=2 in equation (31)
Yx1=propTF(Y21,L,lambda,d1); % Yx1=FRTλ{Y21}z in equation (31)
Yx12=Yx1.*Diagonal2; % Yx12=Yx1*exp[i2πψ(a',b')] in equation (31)
Yx2=propTF(Yx12,L,lambda,d3); % gn(a0,b0)=FRTλ{Y12}dn for n=2 in equation (31)
figure, imshow(Yx2); % Show the second encrypted image
%% Ciphertext
%% g(a0,b0)=PNn=1gn (a0, b0)equation (32)
Cf = Y2+Yx2;
figure, imshow(Cf);
Cfx=real(Cf);
imwrite(Cfx,'CipherText.png'); % Show the single encrypted image
%% If Gaussian noise will be added, this part runs
% [s1, s2] = size(Cfx);
% sigma = 0.2; % σ
% Enoise = zeros(s1,s2);
% Enoise(:,:) = normrnd(0,sigma,s1,s2)+1i*normrnd(0,sigma,s1,s2);
% Cfx = Cfx+Enoise; % equation (48)
%% If occlusion attack is applied, this part runs
% Cfx(83:174,83:174);
%% Decryption stage
%%%%%%%-------------------------------------------------------------%%%%%%%
%% Xn(a,b)=(FRTλ{FRTλ{gk*(a0,b0)}dnexp[i2πψ(a',b')]}z)exp[i2πφ(a,b)] for n=1 and k=1 equation(34)
% Yd1=(FRTλ{gk*(a0,b0)}dn)*exp[i2πψ(a',b')] for n=1 and k=1 in equation (34)
Yd1=propTF(conj(Cf),L,lambda,d2).*Diagonal2;
% Xn(a,b)=(FRTλ{Yd1}z)*exp[i2πφ(a,b)] for n=1 in equation (34)
Yd2=propTF(Yd1,L,lambda,d1).*Diagonal1;
Ydd2=real(Yd2);
figure, imshow(Yd2); % Show the first decrypted image
imwrite(Ydd2,'DecryptedImage1.png');
%% If Gaussian noise or occlusion attack is applied, this part runs
% CC=corr2(Xsingle1, Ydd2); % Correlation coefficient
% MSE=MSE(Xsingle1, Ydd2); % Mean square error
% PSNR=psnr(Xsingle1, Ydd2); % Peak signal to noise ratio (PSNR)
% imwrite(Ydd2,'PosMultipNoise1.png')
% imwrite(Ydd2,'PosMultipOcclusion1.png')
%% Xn(a,b)=(FRTλ{FRTλ{gk*(a0,b0)}dnexp[i2πψ(a',b')]}z)exp[i2πφ(a,b)] for n=2 and k=2 equation(34)
% Ydx1=(FRTλ{gk*(a0,b0)}dn)*exp[i2πψ(a',b')] for n=2 and k=2 in equation (34)
Ydx1=propTF(conj(Cf),L,lambda,d3).*Diagonal2;
% Xn(a,b)=(FRTλ{Ydx1}z)*exp[i2πφ(a,b)] for n=2 in equation (34)
Ydx2=propTF(Ydx1,L,lambda,d1).*Diagonal1;
Ydxx2=real(Ydx2);
figure, imshow(Ydx2); % Show the second decrypted image
imwrite(Ydxx2,'DecryptedImage2.png'); %DecryptedImageWD1.png
%% Decryption with Wrong Wavelength
%% Xn(a,b)=(FRTλ{FRTλ{gk*(a0,b0)}dnexp[i2πψ(a',b')]}z)exp[i2πφ(a,b)] for n=1 and k=1 equation(34)
lambda = 500e-6;
Yd1=propTF(conj(Cf),L,lambda,d2).*Diagonal2;
Yd2=propTF(Yd1,L,lambda,d1).*Diagonal1;
Ydd2=real(Yd2);
figure, imshow(Yd2);
imwrite(Ydd2,'DecryptedImageW1.png');
%% Xn(a,b)=(FRTλ{FRTλ{gk*(a0,b0)}dnexp[i2πψ(a',b')]}z)exp[i2πφ(a,b)] for n=2 and k=2 equation(34)
Ydx1=propTF(conj(Cf),L,lambda,d3).*Diagonal2;
Ydx2=propTF(Ydx1,L,lambda,d1).*Diagonal1;
Ydxx2=real(Ydx2);
figure, imshow(Ydx2);
imwrite(Ydxx2,'DecryptedImageW2.png');
%% Decryption with Wrong z=(d1)
%% Xn(a,b)=(FRTλ{FRTλ{gk*(a0,b0)}dnexp[i2πψ(a',b')]}z)exp[i2πφ(a,b)] for n=1 and k=1 equation(34)
lambda=532e-6;
d1=140;
Yd1=propTF(conj(Cf),L,lambda,d2).*Diagonal2;
Yd2=propTF(Yd1,L,lambda,d1).*Diagonal1;
Ydd2=real(Yd2);
figure, imshow(Yd2);
imwrite(Ydd2,'DecryptedImageWZ1.png');
%% Xn(a,b)=(FRTλ{FRTλ{gk*(a0,b0)}dnexp[i2πψ(a',b')]}z)exp[i2πφ(a,b)] for n=2 and k=2 equation(34)
Ydx1=propTF(conj(Cf),L,lambda,d3).*Diagonal2;
Ydx2=propTF(Ydx1,L,lambda,d1).*Diagonal1;
Ydxx2=real(Ydx2);
figure, imshow(Ydx2);
imwrite(Ydxx2,'DecryptedImageWZ2.png');
%% Decryption with wrong d2 (for n=1)
%% Xn(a,b)=(FRTλ{FRTλ{gk*(a0,b0)}dnexp[i2πψ(a',b')]}z)exp[i2πφ(a,b)] for n=1 and k=1 equation(34)
d1=100;
d2=170;
Yd1=propTF(conj(Cf),L,lambda,d2).*Diagonal2;
Yd2=propTF(Yd1,L,lambda,d1).*Diagonal1;
Ydd2=real(Yd2);
figure, imshow(Yd2);
imwrite(Ydd2,'DecryptedImageWD1.png');
%% Xn(a,b)=(FRTλ{FRTλ{gk*(a0,b0)}dnexp[i2πψ(a',b')]}z)exp[i2πφ(a,b)] for n=2 and k=2 equation(34)
Ydx1=propTF(conj(Cf),L,lambda,d3).*Diagonal2;
Ydx2=propTF(Ydx1,L,lambda,d1).*Diagonal1;
Ydxx2=real(Ydx2);
figure, imshow(Ydx2);
imwrite(Ydxx2,'DecryptedImageWD2.png');
%% Decryption with wrong d3 (for n=2)
%% Xn(a,b)=(FRTλ{FRTλ{gk*(a0,b0)}dnexp[i2πψ(a',b')]}z)exp[i2πφ(a,b)] for n=1 and k=1 equation(34)
d2=160;
d3=200;
Yd1=propTF(conj(Cf),L,lambda,d2).*Diagonal2;
Yd2=propTF(Yd1,L,lambda,d1).*Diagonal1;
Ydd2=real(Yd2);
figure, imshow(Yd2);
imwrite(Ydd2,'DecryptedImageWDD1.png');
%% Xn(a,b)=(FRTλ{FRTλ{gk*(a0,b0)}dnexp[i2πψ(a',b')]}z)exp[i2πφ(a,b)] for n=2 and k=2 equation(34)
Ydx1=propTF(conj(Cf),L,lambda,d3).*Diagonal2;
Ydx2=propTF(Ydx1,L,lambda,d1).*Diagonal1;
Ydxx2=real(Ydx2);
figure, imshow(Ydx2);
imwrite(Ydxx2,'DecryptedImageWDD2.png');