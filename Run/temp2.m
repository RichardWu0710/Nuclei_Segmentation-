
% clear;
% clc;
% close all;

% celldata is a global variable from segmetation step. Contain binary
% information of the location of each single nuclei

ncell = length(celldata);


% Import Data
Data = RNA_Data;

% Extract X, Y, Gene names from the input
CenterX = Data.data(2:end,8);
CenterY = Data.data(2:end,7);
imgSize = Data.data(2,1);

if res == 1
    offset = 0.5*(2048-imgSize);
elseif res == 2
    offset = 0.5*(1024-imgSize);
end

Gene = Data.textdata(2:end,3);
% Gene = [];
% for i = 2:length(Data.textdata)
%     Genename = Data.textdata{i,2};
%     Gene = [Gene ; Genename];
% end


% get all the unique names of genes
UniqueGene = unique(Gene);

% Round x y pixel postions to integer 
for i = 1:length(CenterX)
    CenterX(i) = round(CenterX(i))+offset;
    CenterY(i) = round(CenterY(i))+offset;
    
end
% Plot the genes to check if distrubition makes sense 
a =(2048)/res;
img(1:a,1:a)=0;
img=logical(img);
pixel_val=[CenterX CenterY];
for i=1:size(pixel_val)
    rgb_img(pixel_val(i,1),pixel_val(i,2),:)=[255 0 0];
end
f6 = figure(6);
s = regionprops(ee, 'Centroid');
fuse = imfuse(rgb_img, cellpic);
imshow(fuse);
hold on
for k = 1:numel(s)
    c = s(k).Centroid;
    text(c(1), c(2), sprintf('%d', k), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle');
end
hold off

% Additional identical image for image output
fig = figure(7);
set(fig,'Visible', 'off'); 
imshow(fuse);
hold on
for k = 1:numel(s)
    c = s(k).Centroid;
    text(c(1), c(2), sprintf('%d', k), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle');
end
hold off

stop = 0;
    while stop == 0
        sign = input('Gene projection looks correct? 1 for yes, 2 to skip this position, 3 to repeat analysis on this position\n');
        if sign == 2
            err = 1;
            Output = 0;
            stop = 1;
        elseif sign == 3
            err = 2;
            Output = 0;
            stop = 1;
        elseif isempty(sign) || sign == 1
            stop = 1;
            err = 0;
            %Create a cell array named 'Output' for gene counts of all nuclei
            OutputSize = [length(UniqueGene)+2, ncell+2];
            Output = cell(OutputSize);
            Output{1,1} = 'Gene';
            for i = 1:length(UniqueGene)
                Output{i+1,1} = UniqueGene{i};
            end
            Output{end,1} = 'Total';
            Output{1,end} = 'Total';
            for i = 1:ncell
                temp = int2str(posnum);
                cellname = strcat(temp,'.' ,int2str(i));
                Output{1,i+1} = cellname;
            end

            % Count the number of each gene present in each nuclei
            for i = 1:length(UniqueGene)                    % loop through all genes
                PInd = strfind(Gene,UniqueGene(i));
                PInd = find(not(cellfun('isempty', PInd)));
                GeneX = CenterX(PInd);
                GeneY = CenterY(PInd);
                for j = 1:ncell                             % loop through all nuclei
                    mcell = celldata{j};
                    count = 0;
                    for k = 1:length(GeneX)                 % loop through particular gene
                        if mcell(GeneX(k),GeneY(k))
                            count = count+1;
                        end
                    end
                    Output{i+1,j+1} = count;                % Output total count of one gene in one cell
                end
            end

            %Compute average number of genes per cell
            % Sums of genes for each cell
            for i = 2:ncell+1
                sum = 0;
                for j = 2:length(UniqueGene)+1
                   sum = sum + Output{j,i};
                end
                Output{end,i} = sum;
            end

            % Total number of each type of genes
            for i = 2:length(UniqueGene)+2
                sum = 0;
                for j = 2:ncell +1
                   sum = sum + Output{i,j};
                end
                Output{i,end} = sum;
            end
        else
            disp('Error! Invalid input')
        end
    end
    set(f6,'Visible', 'off'); 



