%% Main.m
clear; close all; clc
%% Read real positive images
Xsingle1 = imread('Message.png');
Xsingle1 = im2double(Xsingle1);
Xsingle2 = imread('Optical_256.png');
Xsingle2 = im2double(Xsingle2);
[yN,xN] = size(Xsingle1);
%% Lambda (wavelength) parameters [mm]
lambda1 = 532e-6;
lambda2 = 132e-6;
%% Distance parameters
d1=100;
d2=160;
%% Construct random phase masks
RandomPhaseTEM1=exp(1j*2*pi*rand(yN,xN)); % Random phase mask 1
Diagonal1=RandomPhaseTEM1;
RandomPhaseTEM2=exp(1j*2*pi*rand(yN,xN)); % Random phase mask 2
Diagonal2=RandomPhaseTEM2;
X1=Xsingle1; % Take first image
L=5;
%% x1(a',b')=FRTλ1{ f(a,b)*exp[i2πφ(a,b)]}d1 equation (27)
Y11=X1.*Diagonal1; % Y11=f(a,b)*exp[i2πφ(a,b)] in equation (27)
Y1=propTF(Y11,L,lambda1,d1); % x1(a',b')=FRTλ1{Y11}d1 in equation (27)
%% c1(a0,b0)=FRTλ1{ x1(a',b')*exp[i2πψ(a',b')]}d2 in equation (28)
Y12=Y1.*Diagonal2; % Y12= x1(a',b')*exp[i2πψ(a',b')] in equation (28)
Y2=propTF(Y12,L,lambda1,d2); % c1(a0,b0)= FRTλ1{Y12}d2 in equation (28)
figure, imshow(Y2);
%% x2(a',b')=FRTλ2{ f(a,b)*exp[i2πφ(a,b)]}d1 equation (27)
X2= Xsingle2; % Take second image
Y21=X2.*Diagonal1; % Y21=f(a,b)*exp[i2πφ(a,b)] in equation (27)
Yx1=propTF(Y21,L,lambda2,d1); % x2(a',b')=FRTλ2{Y21}d1 in equation (27) with λ2.
%% c2(a0,b0)=FRTλ1{ x2(a',b')*exp[i2πψ(a',b')]}d2 in equation (28)
Yx12=Yx1.*Diagonal2; % Yx12=x2(a',b')*exp[i2πψ(a',b')] in equation (28)
Yx2=propTF(Yx12,L,lambda2,d2); % c2(a0,b0)=FRTλ2{Yx12}d2 in equation (28) with λ2.
figure, imshow(Yx2);
%% Encrypted image
Cf=Y2+Yx2; % Encrypted image
figure, imshow(Cf);
Cfx=real(Cf);
imwrite(Cfx,'CipherText.png');
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
Yd1=propTF(conj(Cf),L,lambda1,d2).*Diagonal2;
Yd2=propTF(Yd1,L,lambda1,d1).*Diagonal1;
Ydd2=real(Yd2); % First decrypted image
figure, imshow(Yd2);
imwrite(Ydd2,'DecryptedImage1.png');
%% If Gaussian noise or occlusion attack is applied, this part runs
% CC=corr2(Xsingle1, Ydd2); % Correlation coefficient
% MSE=MSE(Xsingle1, Ydd2); % Mean square error
% PSNR=psnr(Xsingle1, Ydd2); % Peak signal to noise ratio (PSNR)
% imwrite(Ydd2,'WaveMultipNoise1.png')
% imwrite(Ydd2,'WaveMultipOcclusion1.png')
%% Ydx1=FRTλ2{c*(a0,b0)}d2*exp[i2πψ(a',b')] in equation (29) with λ2.
Ydx1=propTF(conj(Cf),L,lambda2,d2).*Diagonal2;
% ḟ2(a,b)=FRTλ2{Ydx1}d1*exp[i2πφ(a,b)] in equation (29) with λ2.
Ydx2=propTF(Ydx1,L,lambda2,d1).*Diagonal1;
Ydxx2=real(Ydx2);
figure, imshow(Ydx2);
imwrite(Ydxx2,'DecryptedImage2.png');
%% Decryption with wrong distance, d1.
d1=180; % wrong distance, d1
% Yd1=FRTλ1{c*(a0,b0)}d2*exp[i2πψ(a',b')] in equation (29)
Yd1=propTF(conj(Cf),L,lambda1,d2).*Diagonal2;
Yd2=propTF(Yd1,L,lambda1,d1).*Diagonal1;
Ydd2=real(Yd2); % First decrypted image
figure, imshow(Yd2);
imwrite(Ydd2,'DecryptedImageW1d1.png');
%% Ydx1=FRTλ2{c*(a0,b0)}d2*exp[i2πψ(a',b')] in equation (29) with λ2.
Ydx1=propTF(conj(Cf),L,lambda2,d2).*Diagonal2;
Ydx2=propTF(Ydx1,L,lambda2,d1).*Diagonal1;
Ydxx2=real(Ydx2); % Second decrypted image
figure, imshow(Ydx2);
imwrite(Ydxx2,'DecryptedImageW2d1.png');
%% Decryption with wrong d2
d1=100; % True d1
d2=80; % Wrong d2
% Yd1=FRTλ1{c*(a0,b0)}d2*exp[i2πψ(a',b')] in equation (29)
Yd1=propTF(conj(Cf),L,lambda1,d2).*Diagonal2;
Yd2=propTF(Yd1,L,lambda1,d1).*Diagonal1;
Ydd2=real(Yd2); % First decrypted image
figure, imshow(Yd2);
imwrite(Ydd2,'DecryptedImageWr1d2.png');
%% Ydx1=FRTλ2{c*(a0,b0)}d2*exp[i2πψ(a',b')] in equation (29) with λ2.
Ydx1=propTF(conj(Cf),L,lambda2,d2).*Diagonal2;
Ydx2=propTF(Ydx1,L,lambda2,d1).*Diagonal1;
Ydxx2=real(Ydx2); % Second decrypted image
figure, imshow(Ydx2);
imwrite(Ydxx2,'DecryptedImageWr2d2.png');
%% Decryption with wrong wavelength1
d1=100; % True d1
d2=160; % True d2
lambda1 = 500e-6;
% Yd1=FRTλ1{c*(a0,b0)}d2*exp[i2πψ(a',b')] in equation (29)
Yd1=propTF(conj(Cf),L,lambda1,d2).*Diagonal2;
% ḟ1(a,b)=FRTλ1{Yd1}d1*exp[i2πφ(a,b)] in equation (29
Yd2=propTF(Yd1,L,lambda1,d1).*Diagonal1;
Ydd2=real(Yd2); % First decrypted image
figure, imshow(Yd2);
imwrite(Ydd2,'DecryptedImageWr1W1.png');
%% Ydx1=FRTλ2{c*(a0,b0)}d2*exp[i2πψ(a',b')] in equation (29) with λ2.
Ydx1=propTF(conj(Cf),L,lambda2,d2).*Diagonal2;
Ydx2=propTF(Ydx1,L,lambda2,d1).*Diagonal1;
Ydxx2=real(Ydx2); % Second decrypted image
figure, imshow(Ydx2);
imwrite(Ydxx2,'DecryptedImageWr2W1.png');
%% Decryption with wrong wavelength2
lambda1 = 525e-6; % True wavelength 1
lambda2 = 100e-6; % False wavelength 2
% Yd1=FRTλ1{c*(a0,b0)}d2*exp[i2πψ(a',b')] in equation (29)
Yd1=propTF(conj(Cf),L,lambda1,d2).*Diagonal2;
Yd2=propTF(Yd1,L,lambda1,d1).*Diagonal1;
Ydd2=real(Yd2); % First decrypted image
figure, imshow(Yd2);
imwrite(Ydd2,'DecryptedImageWr1W2.png');
%% Ydx1=FRTλ2{c*(a0,b0)}d2*exp[i2πψ(a',b')] in equation (29) with λ2.
Ydx1=propTF(conj(Cf),L,lambda2,d2).*Diagonal2;
Ydx2=propTF(Ydx1,L,lambda2,d1).*Diagonal1;
Ydxx2=real(Ydx2); % Second decrypted image
figure, imshow(Ydx2);
imwrite(Ydxx2,'DecryptedImageWr2W2.png');