% Segmentation 

clear;
clc;
close all;

I = imread('Position2002_z19_ch01.TIF');
figure (1)
imshow(I);
title('Original Image')

% Adjust Contrast
I3 = imadjust(I,[0.15,1]); % Adjust
figure (3)
imshow(I3);

% Turn image into binary format
% Adjust https://www.mathworks.com/help/images/ref/imbinarize.html
bw = imbinarize(I3);       

% Remove noise
bw = bwareaopen(bw, 50);
% figure (4)
% imshow(bw);

% cc contain pixels for all objects
cc = bwconncomp(bw, 18);
pic = false(size(bw));
pic(cc.PixelIdxList{10}) = true;
% figure (5)
% imshow(pic);

% allimage contains all matrix for each crude cell object
n = length(cc.PixelIdxList);
allimage = cell(n);
for i = 1:length(cc.PixelIdxList)
    pic = false(size(bw));
    pic(cc.PixelIdxList{i}) = true;
    allimage{i}= pic;
end

% smooth edges and fill hole for each crude cell, compute area
for i = 1:n
    bw = allimage{i};
    E = edge(bw);
    Ed = imdilate(E,strel('sphere',3));
    bw(Ed) = Ed(Ed);
    BW = imfill(bw,'holes');
    allimage{i} = BW;
    allimage{i,2} = regionprops(BW,'Area');
end

%%
% show entire picture
%
wholepic = zeros(size(bw));
for i = 1:n
    wholepic = wholepic+allimage{i};
end
figure (2)
imshow(wholepic);
title('All Objects After Initial Filtering')

%%
% plot  area value for each cell
% can be used to select threshold
%
areabar = zeros(n,1);
for i=1:n
    areabar(i)= allimage{i,2}.Area;
end

figure (3)
subplot(2,1,1)
bar(areabar);
grid on;
grid minor;
title('CellIndex vs. Area');
ylabel('Area');
xlabel('Index of the cell');
xlim([1 60]);
Aline = sort(areabar);

subplot(2,1,2)
bar(Aline);
grid on;
grid minor;
title('Sorted Area for selecting area threshold');
ylabel('Area');


%%
% seperate and plot cells based on area threshold
% multiObj, singleObj, debrisObj contain conditionaly sorted cell images
%

% Ask user input for debris area threshold
%debrislim define area limit bellow which is debri
debrislim = input('Area threshold for debris? hint:0.15*10^4 \n');
if isempty(debrislim)
     debrislim = 0.15*10^4;
end

% Ask user input for single cell area threshold
% singlelim define area limit above which is no long a single cell
singlelim = input('Area threshold for single cell? hint:1.1*10^4 \n');
if isempty(singlelim)
     singlelim = 1.1*10^4;
end

% find object indexes from the allimage array that
% satisfy each condition
singleObjInd = find(areabar<= singlelim ...
    & areabar>debrislim);
debrisObjInd = find(areabar<=debrislim);
multiObjInd = find(areabar > singlelim);

% number of cells for each condition
nmulti = sum(length(multiObjInd));
nsingle = sum(length(singleObjInd));
ndebris = sum(length(debrisObjInd));

% create structures contain cell images of each condition
multiObj = cell(nmulti);
singleObj = cell(nsingle);
debrisObj = cell(ndebris);

% Seperate cells into each condition structure,
% compute area for each cell
for i = 1:nmulti
    multiObj{i} = allimage{multiObjInd(i)};
    multiObj{i,2} = regionprops(multiObj{i},'Area');
end

for i = 1:nsingle
    singleObj{i} = allimage{singleObjInd(i)};
    singleObj{i,2} = regionprops(singleObj{i},'Area');
end

for i = 1:ndebris
    debrisObj{i} = allimage{debrisObjInd(i)};
    debrisObj{i,2} = regionprops(debrisObj{i},'Area');
end

%%
% plot result for above section: showing seperated cells 
% in each condition
figure (5)

multiObjpic = zeros(size(bw));
for i = 1:nmulti
    multiObjpic = multiObjpic+multiObj{i};
end
subplot(2,2,1)
imshow(multiObjpic);
title('All Clustersof Multiple Cells')

singleObjpic = zeros(size(bw));
for i = 1:nsingle
    singleObjpic = singleObjpic+singleObj{i};
end
subplot(2,2,2)
imshow(singleObjpic);
title('All Single Cells')

debrisObjpic = zeros(size(bw));
for i = 1:ndebris
    debrisObjpic = debrisObjpic+debrisObj{i};
end
subplot(2,2,3)
imshow(debrisObjpic);
title('All Debris')


%% 
% Modification by strel, gaussian smoothing and watershed 
% for multicell clusters
newMultiObj = cell(size(multiObj));
figure(11)
for i = 1:length(multiObj)
obj = multiObj{i};
BWnobord = imclearborder(obj, 1);
seD = strel('disk',1);
BWfinal = imerode(BWnobord,seD);
n=1;
for j = 1:n
BWfinal = imerode(BWfinal,seD);
end
    
obj = ~BWfinal;
D = bwdist(obj);
D = imgaussfilt(D,6);

D = -D;
D(obj) = -Inf;

conn = conndef(2,'maximal');
newMultiObj{i}= watershed(D);
rgb = label2rgb(newMultiObj{i},'jet',[0.5 0.5 0.5]);
subplot(2,length(multiObj),i)
imshow(rgb,'InitialMagnification','fit')
title('Watershed transform of multi cell')
end


%%
% Add newly seperated single cells to singleObj array and their areas
nnew = nsingle;

% nexted for loop goes through each array of newly seperated cells
for i = 1: length(newMultiObj)
    gg = bwconncomp(newMultiObj{i}, 4);
    oldSingle = singleObj;                  % temperory storage 
    cellsize = length(singleObj)+length(gg.PixelIdxList)-1;
    singleObj = cell(cellsize);             % create longer cell array
    % for loop record all old values back
    for j = 1: length(oldSingle)            
        singleObj{j} = oldSingle{j};
        singleObj{j,2} = oldSingle{j,2};
    end
    % for loop add new values
    for add = 1:length(gg.PixelIdxList)-1
        pic = false(size(newMultiObj{i}));
        pic(gg.PixelIdxList{add+1}) = true;
        newind = nnew+add;
        singleObj{newind}= pic;
        singleObj{newind,2} = regionprops(pic,'Area');  % Compute Area
    end
    nnew = nnew + length(gg.PixelIdxList)-1;         
end


%%
% Plot all single cells together after watershed
singleObjpic = zeros(size(bw));
for i = 1:length(singleObj)
    singleObjpic = singleObjpic+singleObj{i};
end
figure (15)
imshow(singleObjpic);
title('All Single Cells after Watershed')

%%
% % make single cell images
% for i = 1:length(singleObj)
%     str = strcat('Single_Cell_',num2str(i),'_tif');
%     imwrite(uint8(singleObj{i}),str,'tif');
% end

