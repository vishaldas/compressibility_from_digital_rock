function [ K_HS, U_HS ] = createsegimage(image_name, contact_image_name, summary_file_name, filename_out )
%[ K_HS, U_HS ] = createsegimage(image_name, contact_image_name, summary_file_name, filename_out )
%   Function to create segmented image from CAP image
%
%   createsegimage takes in CAP dataset and resegments them with codes
%   assigned as per codes in dictionary_elastic_properties.xslx file 
%   located in Dependencies folder. The outputs from the function are the 
%   resegmented image in .raw format and the Mineral bulk and shear modulus
%   
%
%   Inputs:
%   image_name = name of the image file
%   contact_image_name = name of the contact image file 
%   summary_filename = name of the summary file 
%   filename_out = name of the output raw file 
%
%   Outputs:
%   Resegmented image 
%   K_HS = Hashin-Shtrikman average mineral bulk modulus (GPa)
%   U_HS = Hashin-Shtrikman average mineral shear modulus (GPa)
%
%   Example:
%   image_name = 'sample1.tif';
%   summary_file_name = 'summary_sample1.txt';
%   contact_image_name = 'contact_sample1.tif';
%   filename_out = 'test.raw';
%   [ K_HS, U_HS ] = createsegimage(image_name, contact_image_name, summary_file_name, filename_out );
%
%   Written by Vishal Das, September 2018


% Add path of Dependencies folder
addpath('Dependencies');

% Image file name
filename = image_name;
rgb = imread(filename);
rgb_columns = reshape(rgb, [], 3);
[~, ~, n] = unique(rgb_columns, 'rows');
color_counts = accumarray(n, 1);
% Percent each color
percent_each_color = (color_counts./(size(rgb,1)*size(rgb,2))).*100;
[sort_rgb, sort_indx_rgb] = sort(percent_each_color, 'descend');
    
% Read summary text file and extract values for each row
    
% Determine the start and end of mineral summary
checkstart = 'MINERAL SUMMARY --------------------------------------------------------------------';
checkend = 'CONTACT SUMMARY --------------------------------------------------------------------';
filename_summary = summary_file_name;
fid  = fopen(filename_summary,'r');
i = 1;
tline = fgetl(fid);
while ischar(tline)
    i = i+1;
    if (strcmp(checkstart, tline) == 1)
        startRow = i+2; % 2 added since the actual data starts two lines below
    end
    if (strcmp(checkend, tline) == 1)
        endRow = i-3; % 3 subtracted since the actual data ends two lines before
    end
    tline = fgetl(fid);
end
fclose(fid);
% Read the data
startRow = startRow - 2;
endRow = endRow - 4;
formatSpec = '%s%s%s%s%s%[^\n\r]';
delimiter = ' ';
fileID = fopen(filename_summary,'r');
textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'ReturnOnError', false, 'EndOfLine', '\r\n');
fclose(fileID);
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[3,4,5]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end

rawNumericColumns = raw(:, [3,4,5]);
rawStringColumns = string(raw(:, [1,2]));

R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

summary1 = raw;

clearvars dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp rawNumericColumns rawStringColumns R;

% Reformating the summary1 file to two columns
summary_new = cell( size(summary1,1),2);

for i = 1:size(summary1,1)
    if (ischar(summary1{i,3}))
        summary_new{i,1} = strcat(summary1{i,1}, summary1{i,2});
        summary_new{i,2} = summary1{i,5};
    else
        summary_new{i,1} = strcat(summary1{i,1});
        summary_new{i,2} = summary1{i,4};
    end
end

% Comparing with output from resegmented image
summary_sort = sortrows(summary_new, -2);
    
% Assigning pixels with appropriate code numbers
length_new = min(length(sort_rgb), size(summary_sort,1));
new_image = zeros(size(rgb,1), size(rgb,2));
load('codeminerals.mat');
for i = 1: length_new
    % Get the code of the mineral from the code table
    name_mineral = summary_sort{i,1};
    code_mineral = codeminerals(strcmp(codeminerals.mineral_name, name_mineral), 2);
    
    % Replacing codes outside dictionary as pore
    if (isempty(code_mineral))
        code_mineral = codeminerals(strcmp(codeminerals.mineral_name, string(Intergranular)), 2);
    end
    idx = sort_indx_rgb(i);
    bw = n == idx;
    new_image(bw) = table2array(code_mineral);
    
end

% Reading in contact data
contact_data = imread(contact_image_name);
% Hard contact (quartz = code0)
new_image(contact_data == 255 | contact_data == 205) = 0;
% Soft contact (clay = code4)
new_image(contact_data == 128 | contact_data == 30) = 4;

% Calculate the mineral bulk and shear modulus
[unique_color_final, ~, n_final] = unique(new_image);
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

[~,~,~,~,K_HS,U_HS] = bound(1,frac_mineral.',K_mineral.',U_mineral.');
    
fid=fopen(filename_out,'w');
fwrite(fid,new_image,'uint8');
fclose(fid);

end

