% Segmentation 

clear;
clc;
close all;

% sample name
% MAX_Cycle7_Position1001_ch00.tif
% 'Sample.jpg'

I = imread('MAX_12-5-2016 Position2.tif');
figure (1)
imshow(I);
title('Original Image')

% Adjust Contrast
I3 = imadjust(I); % Adjust
% figure (3)
% imshow(I3);

% Turn image into binary format
% Adjust https://www.mathworks.com/help/images/ref/imbinarize.html
bw = imbinarize(I3);       

% Remove noise
bw = bwareaopen(bw, 50);
cc = bwconncomp(bw, 18);

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
    Ed = imdilate(E,strel('sphere',1));
    bw(Ed) = Ed(Ed);
    bw = imfill(bw,'holes');
    allimage{i} = bw;
    allimage{i,2} = regionprops(bw,'Area');
end

%%
% show entire picture
%
wholepic1 = zeros(size(bw));
for i = 1:n
    wholepic1 = wholepic1+allimage{i};
end
figure (2)
imshow(wholepic1);
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

% Ask user input for debris area threshold
%debrislim define area limit bellow which is debri
debrislim = input('Area threshold for debris? default:0.15*10^4 \n');
if isempty(debrislim)
     debrislim = 0.15*10^4;
end

% find object indexes from the allimage array that
debrisObjInd = find(areabar<=debrislim);

% number of cells for each condition
ndebris = sum(length(debrisObjInd));

% create structures contain cell images of each condition
debrisObj = cell(ndebris);

for i = 1:ndebris
    debrisObj{i} = allimage{debrisObjInd(i)};
    debrisObj{i,2} = regionprops(debrisObj{i},'Area');
end

%%
% select the multicell cluster from the rest
selectInd = 1:1:length(allimage);
for i = 1:length(debrisObjInd)
    badind = find(selectInd == debrisObjInd(i));
    selectInd(badind) = [];
end

wholepic = zeros(size(bw));
for i = 1:length(selectInd)
    wholepic = wholepic+allimage{selectInd(i)};
end
L = bwlabel(wholepic1);

figure (101)
s = regionprops(cc, 'Centroid');
imshow(wholepic)
hold on
for k = 1:numel(s)
    c = s(k).Centroid;
    text(c(1), c(2), sprintf('%d', k), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle');
end
hold off

% ask to input the index of multi cell
flag = 0;
multiObjInd = [];
while flag == 0
    num = input('Input the index of multicell cluster,one at a time. if no, input 0:\n');
    if num <= length(allimage) && num >= 1
        multiObjInd = [multiObjInd num];
    else
      flag = 1;
    end
end

nmulti = sum(length(multiObjInd));
multiObj = cell(nmulti);

for i = 1:nmulti
    multiObj{i} = allimage{multiObjInd(i)};
    multiObj{i,2} = regionprops(multiObj{i},'Area');
end

allInd = 1:1:length(allimage);
singleObjInd = setdiff(allInd,debrisObjInd);
singleObjInd = setdiff(singleObjInd,multiObjInd);
nsingle = sum(length(singleObjInd));
singleObj = cell(nsingle);
for i = 1:nsingle
    singleObj{i} = allimage{singleObjInd(i)};
    singleObj{i,2} = regionprops(singleObj{i},'Area');
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
newMultiObj = cell(nmulti);
figure(11)

for i = 1:nmulti
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
D = imgaussfilt(D,2);

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
    bw = bwareaopen(newMultiObj{i}, 50);
    gg = bwconncomp(bw, 26);
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

%% filter out small debris from newly seperated cells
n = length(singleObj);
Newdebrislim = input('Area threshold for newly seperated debris? default:100 \n');
if isempty(Newdebrislim)
     Newdebrislim = 100;
end
for i=1:n
    bw = singleObj{i};
    singleObj{i} = bw;
    singleObj{i,2} = regionprops(bw,'Area');
end
areabar = zeros(n,1);
for i=1:n
    areabar(i)= singleObj{i,2}.Area;
end
LargeInd = find(areabar>Newdebrislim);
numSingle = length(LargeInd);
finalSingleObj = cell(numSingle);
for i = 1:numSingle
    finalSingleObj{i} = singleObj{LargeInd(i)};
    finalSingleObj{i,2} = regionprops(finalSingleObj{i},'Area');
end

% Plot all single cells together after watershed
singleObjpic = zeros(size(bw));
for i = 1:length(finalSingleObj)
    singleObjpic = singleObjpic+finalSingleObj{i};
end
figure (15)
imshow(singleObjpic);
title('All Single Cells after Watershed')

global celldata 
celldata = finalSingleObj;
global cellpic
cellpic = singleObjpic;
