%Segmentation 

clear;
clc;
close all;



I = imread('Sample.jpg');
imshow(I);
background = imopen(I,strel('arbitrary',10));

% Display the Background Approximation as a Surface
% figure
% surf(double(background(1:8:end,1:8:end))),zlim([0 255]);
% ax = gca;
% ax.YDir = 'reverse';


I2 = background;
% figure (2)
% imshow(I2)

% Adjust Contrast
I3 = imadjust(I,[0.15,1]); % Adjust
% figure (3)
% imshow(I3);


bw = imbinarize(I3);       % Adjust https://www.mathworks.com/help/images/ref/imbinarize.html
bw = bwareaopen(bw, 50);
BW1 = bw;
% figure (4)
% imshow(bw);

% cc contain pixels for all objects
cc = bwconncomp(bw, 18);
pic = false(size(bw));
pic(cc.PixelIdxList{10}) = true;
% figure (5)
% imshow(pic);

% allimage contain matrix for each crude cell object
n = length(cc.PixelIdxList);
allimage = cell(n);
for i = 1:length(cc.PixelIdxList)
    pic = false(size(bw));
    pic(cc.PixelIdxList{i}) = true;
    allimage{i}= pic;
end

% smooth edges and fill hole for each crude cell
for i = 1:n
    bw = allimage{i};
    E = edge(bw);
    Ed = imdilate(E,strel('sphere',3));
    bw(Ed) = Ed(Ed);
    BW2 = imfill(bw,'holes');
    allimage{i} = BW2;
    allimage{i,2} = regionprops(BW2,'Area');
end

%%
% plot  area value for each cell
% can be used to select threshold
%
areabar = zeros(n,1);
for i=1:n
    areabar(i)= allimage{i,2}.Area;
end

figure (1)
bar(areabar);
grid on;
grid minor;
title('CellIndex vs. Area');
ylabel('Area');
xlabel('Index of the cell');

Aline = sort(areabar);
figure (2)
bar(Aline);
grid on;
grid minor;
title('Sorted Area for selecting area threshold');
ylabel('Area');


%%
% random ploting section for users
figure (3)
imshow(allimage{16});

%%
% show entire picture
%

wholepic = zeros(size(bw));
for i = 1:n
    wholepic = wholepic+allimage{i};
end
figure (4)
imshow(wholepic);

%%
% seperate and plot cells based on area threshold
% multiObj, singleObj, debrisObj contain conditionaly sorted cell images
%

% debrislim define area limit bellow which is debries
debrislim = 0.15*10^4;

% singlelim define area limit above which is not single cell
singlelim = 1.1*10^4;

% find object indexes of each condition
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

% Seperate cells into each condition structure
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

% plot whole picture for each condtion
multiObjpic = zeros(size(bw));
singleObjpic = zeros(size(bw));
debrisObjpic = zeros(size(bw));

for i = 1:nmulti
    multiObjpic = multiObjpic+multiObj{i};
end
figure (5)
imshow(multiObjpic);

for i = 1:nsingle
    singleObjpic = singleObjpic+singleObj{i};
end
figure (6)
imshow(singleObjpic);

for i = 1:ndebris
    debrisObjpic = debrisObjpic+debrisObj{i};
end
figure (7)
imshow(debrisObjpic);


%%
% random plotting section for users
figure (3)
imshow(multiObj{1});
%% 
%watershed for multi cells

newMultiObj = cell(size(multiObj));
for i = 1:length(multiObj)
obj = ~multiObj{i};
D = bwdist(obj);
% DD = imcrop(multiObj{1});
% figure
% imshow(D,[],'InitialMagnification','fit')

% title('Distance transform of ~bw');
D = -D;
D(obj) = -Inf;
conn = conndef(2,'maximal');
newMultiObj = watershed(D,26);
rgb = label2rgb(newMultiObj,'jet',[0.5 0.5 0.5]);
figure (11+i)
imshow(rgb,'InitialMagnification','fit')
title('Watershed transform of multi cell')
end
