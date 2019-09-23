function [ Cmin, Cmax ] = calcm( image_name, row, col, mode, K_HS, U_HS )
%[ Cmin, Cmax ] = calccm( image_name, row, col, mode, K_HS, U_HS )
%   calcm takes a segmented image and uses Digital Rock to predict the
%   compressibility at 1500 psi depletion stress and maximum
%   compressibility
%
%   Inputs:
%   image_name = input segmented image in .raw format
%   row = length of the image (pixels)
%   col = width of the image (pixels)
%   mode = type of uniaxial experiment (1 = X direction, 2 = Y direction)
%   K_HS = Hashin-Shtrikman average mineral bulk modulus (GPa)
%   U_HS = Hashin-Shtrikman average mineral shear modulus (GPa)
%
%   Outputs:
%   Folder named <image_name>_simulations with all numerical results
%   Cmin = Predicted compressibility at 1500 psi (microsips)
%   Cmax = Predicted maximum compressibility (microsips)
%
%   Written by Vishal Das, September 2018

% Add path of Dependencies folder
addpath('Dependencies');


% Check if a folder exists with the same image_name
if (exist([image_name '_simulations'], 'dir') == 7)
    error(['Folder ' image_name '_simulations already exists. Delete folder and retry']);
end

% Remove temp directory
system('rmdir /s/q .\temp');
mkdir('temp');
cd('temp');
% Load and do calculations of stress and strains in NIST code
system('copy ..\Dependencies\elas3d.exe .\');
system('copy ..\Dependencies\elas3d-uniaxial_x.pam .\');
system('copy ..\Dependencies\elas3d-uniaxial_y.pam .\');
system(['copy ..\' image_name ' .\new_image.raw']);

if mode == 1
    
    % Uniaxial calculations
    mkdir('Uniaxial_X');
    system('copy elas3d.exe .\Uniaxial_X\');
    system('copy elas3d-uniaxial_x.pam .\Uniaxial_X\elas3d.pam');
    system('copy new_image.raw .\Uniaxial_X\');
    
    cd('.\Uniaxial_X');
    
    exec_status = system(['elas3d']);
    
    % Remove unnecessary files
    delete('elas3d.exe');
    
    % Loading results and calculating compressional modulus
    height = 1;
    stress_temp = dlmread('stressField.dat');
    stress11_temp = stress_temp(:,1);
    stress_11 = reshape(stress11_temp, row, col, height);
    
    strain_11 = 0.001;
    comp_modulus = mean(stress_11(:)./strain_11(:));
    
    % Mineral modulus calculation
    M_mineral = (K_HS+ 4*U_HS/3);
    
    Cm_DR = (1./(comp_modulus.*145037.73773)).*1e6; % From GPa to microsips
    Cm_mineral = (1./(M_mineral.*145037.73773)).*1e6; % From GPa to microsips
    n1 = 1.2557; % Exponent for Cmax
    n2 = 0.97932; % Exponent for Cmin

    % Empirical relationship for predicting Cmin and Cmax
    Cmin = (Cm_DR).^n2.*(Cm_mineral).^(1-n2);
    Cmax = (Cm_DR).^n1.*(Cm_mineral).^(1-n1);
    
    % Clean up
    cd ('..');
    delete('elas3d.exe');
    delete('elas3d-uniaxial_x.pam'); delete('elas3d-uniaxial_y.pam');
    delete('new_image.raw');
    
    % Rename temp folder to the sample name
    cd('..');
    dir_name = [image_name '_simulations'];
    fclose('all');
    movefile('.\temp', dir_name);
    
elseif (mode ==2)
    mkdir('Uniaxial_Y');
    system('copy elas3d.exe .\Uniaxial_Y\');
    system('copy elas3d-uniaxial_y.pam .\Uniaxial_Y\elas3d.pam');
    system('copy new_image.raw .\Uniaxial_Y\');
    
    cd('.\Uniaxial_Y');
    
    exec_status = system(['elas3d']);
    
    % Remove unnecessary files
    delete('elas3d.exe');
    
    % Loading results and calculating compressional modulus
    height = 1;
    stress_temp = dlmread('stressField.dat');
    stress22_temp = stress_temp(:,2);
    stress_22 = reshape(stress22_temp, row, col, height);

    
    strain_22 = 0.001;
    comp_modulus = mean(stress_22(:)./strain_22(:));
    
    % Mineral modulus calculation
    M_mineral = (K_HS+ 4*U_HS/3);
    
    Cm_DR = (1./(comp_modulus.*145037.73773)).*1e6; % From GPa to microsips
    Cm_mineral = (1./(M_mineral.*145037.73773)).*1e6; % From GPa to microsips
    n1 = 1.419; % Exponent for Cmax
    n2 = 1.1031; % Exponent for Cmin

    % Empirical relationship for predicting Cmin and Cmax
    Cmin = (Cm_DR).^n2.*(Cm_mineral).^(1-n2);
    Cmax = (Cm_DR).^n1.*(Cm_mineral).^(1-n1);       
    
    % Clean up
    cd ('..');
    delete('elas3d.exe');
    delete('elas3d-uniaxial_x.pam'); delete('elas3d-uniaxial_y.pam');
    delete('new_image.raw');
    
    % Rename temp folder to the sample name
    cd('..');
    dir_name = [image_name '_simulations'];
    fclose('all');
    movefile('.\temp', dir_name);
end


end

