%% propTF.m
function[u2]=propTF(u1,L,lambda,z)
% Fresnel propagation using the Transfer function
% Based on Computational Fourier Optics by Voelz 
% Assuming uniform sampling and presents reflections on the boundaries
%PARAMETERS
% 
% u1      - Complex Amplitude of the beam at the source plane
% L       - Sidelength of the simulation window of the source plane
% lambda  - Wavelength 
% z       - Propagation distance 
% u2      - Complex Amplitude of the beam at the observation plane
 
k=2*pi/lambda; %wavenumber
% Input array size
[M,~]=size(u1); 
% Sampling interval size
dx=L/M;   
% Frequency coordinates sampling 
fx=-1/(2*dx):1/L:1/(2*dx)-1/L; 
% Momentum/reciprocal Space 
[FX,FY]=meshgrid(fx,fx);
% Transfer function
H=exp(-1j*k*z).*exp(-1j*pi*lambda*z*(FX.^2+FY.^2)); 
H=fftshift(H);
% Fourier transform of the transfer fucntion
U1=fft2(fftshift(u1));
% Convolution of the system
U2=H.*U1;
% Fourier transform of the convolution to the observation plane
u2=ifftshift(ifft2(U2));  
end

