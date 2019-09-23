%% Script description
% Script to perform finite element elastic calculations 
% using elas3d code from NIST
% Performs Uniaxial strain calculations for X direction or Y direction
% The parameter file for emc3d code is located in the 
% Dependencies folder. 
% The parameter file contains the image size that can be changed manually. The 
% default image size is 480 x 512. This should be changed in the 
% elas3d-uniaxial_x.pam and elas3d-uniaxial_y.pam files in the
% Dependencies folder
% 
% Written by Vishal Das, September 2018


% Performing Digital Rock calculations and predicting compressibility
filename_out = 'test.raw';
fid = fopen(filename_out);
row = 480; col = 512; % Size of segmented image
temp_image = fread(fid, [row, col], 'uint8');

% Calculate mineral bulk and shear modulus 

addpath('Dependencies');
load('codeminerals.mat');

[unique_color_final, ~, n_final] = unique(temp_image);
color_counts_final = accumarray(n_final, 1);
frac_final = color_counts_final ./ sum(color_counts_final);

frac_final(unique_color_final == 1 | unique_color_final == 14 | ...
    unique_color_final == 17 | unique_color_final == 18 | ...
    unique_color_final == 20 | unique_color_final == 51 | ...
    unique_color_final == 53 | unique_color_final == 54 | ...
    unique_color_final == 55 | unique_color_final == 57 | ...
    unique_color_final == 58) = 0;

frac_mineral = frac_final./ sum(frac_final);

K_mineral = table2array(codeminerals(unique_color_final+1,3));
U_mineral = table2array(codeminerals(unique_color_final+1,4));

frac_mineral = frac_mineral.';
K_mineral = K_mineral.';
U_mineral = U_mineral.';

% Calculating Hashin-Strikman upper and lower bounds 
c=4/3;

kmx=max(K_mineral);
kmn=min(K_mineral);
umx=max(U_mineral);
umn=min(U_mineral);

k_u=1/sum(frac_mineral./(K_mineral+c*umx))-c*umx;	% HS upper bound
k_l=1/sum(frac_mineral./(K_mineral+c*umn))-c*umn;	% HS lower bound

etamx=umx*(9*kmx+8*umx)/(kmx+2*umx)/6;
etamn=umn*(9*kmn+8*umn)/(kmn+2*umn)/6;

u_u=1/sum(frac_mineral./(U_mineral+etamx))-etamx;	% HS upper bound
u_l=1/sum(frac_mineral./(U_mineral+etamn))-etamn;	% HS lower bound

K_HS=(k_u+k_l)/2;			% simple arithmetic average
U_HS=(u_u+u_l)/2;


% Compressibility calculations using Digital Rock simulation 
mode = 1;

% Convert binary image to ascii
filename_out_ascii = 'test.dat';
dlmwrite(filename_out_ascii, reshape(temp_image, [size(temp_image,1)*size(temp_image,2), 1]));
[ Cmin, Cmax ] = calcm( filename_out_ascii, row, col, mode, K_HS, U_HS );

