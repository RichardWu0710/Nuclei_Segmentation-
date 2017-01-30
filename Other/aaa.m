

% celldata is a global variable from segmetation step. Contain binary
% information of the location of each single nuclei

ncell = length(celldata);

% Import Data
Data = importdata('SampleData.xlsx');

% Extract X, Y, Gene names from the input
CenterX = Data.data(2:end,7);
CenterY = Data.data(2:end,6);
Gene = Data.textdata(2:end,2);
% Gene = [];
% for i = 2:length(Data.textdata)
%     Genename = Data.textdata{i,2};
%     Gene = [Gene ; Genename];
% end


% get all the unique names of genes
UniqueGene = unique(Gene);

% Round x y pixel postions to integer 
for i = 1:length(CenterX)
    CenterX(i) = round(CenterX(i));
    CenterY(i) = round(CenterY(i));
    
end
% Plot the genes to check if distrubition makes sense 
a=1024;
img(1:a,1:a)=0;
img=logical(img);
pixel_val=[CenterX CenterY];
for i=1:size(pixel_val)
    rgb_img(pixel_val(i,1),pixel_val(i,2),:)=[255 0 0];
end
figure(101)
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
    cellname = ['Cell ' int2str(i)];
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

 Compute average number of genes per cell
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

xlswrite('Spacial_Analysis.xlsx',Output)
