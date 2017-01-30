I = imread('Sample.jpg');
% imshow(I);
background = imopen(I,strel('arbitrary',10));

% Display the Background Approximation as a Surface
% figure
% surf(double(background(1:8:end,1:8:end))),zlim([0 255]);
% ax = gca;
% ax.YDir = 'reverse';

I2 = background;
% figure (2)
% imshow(I2)

I3 = imadjust(I,[0.15,1]); % Adjust
% figure (3)
% imshow(I3);


bw = imbinarize(I3);       % Adjust https://www.mathworks.com/help/images/ref/imbinarize.html
bw = bwareaopen(bw, 50);
% figure (4)
% imshow(bw);

cc = bwconncomp(bw, 18);
pic = false(size(bw));
pic(cc.PixelIdxList{10}) = true;
% figure (5)
% imshow(pic);
n = length(cc.PixelIdxList);
allimage = cell(n);
for i = 1:length(cc.PixelIdxList)
    pic = false(size(bw));
    pic(cc.PixelIdxList{i}) = true;
    allimage{i}= pic;
end

for i = 1:n
    bw = allimage{i};
    E = edge(bw);
    Ed = imdilate(E,strel('sphere',3));
    bw(Ed) = Ed(Ed);
    BW2 = imfill(bw,'holes');
    allimage{i} = BW2;
end
%%
figure(17)
labeled = labelmatrix(picc);
whos labeled
RGB_label = label2rgb(labeled,@spring, 'c', 'shuffle');
imshow(RGB_label)
%%
%smooth image
E = edge(bw);
Ed = imdilate(E,strel('sphere',3));
% figure(7)
% imshow(E);
% figure (8)
% imshow(Ed);

bw(Ed) = Ed(Ed);

BW2 = imfill(bw,'holes');
figure (10)
imshow(BW2)
Ifilt = imfilter(BW2,fspecial('gaussian'));
figure(9)
imshow(Ifilt);

%% 
%watershed
D = bwdist(~BW2);
% figure
% imshow(D,[],'InitialMagnification','fit')

% title('Distance transform of ~bw');
D = -D;
D(~bw) = -Inf;
conn = conndef(2,'maximal');
L = watershed(D,conn);
rgb = label2rgb(L,'jet',[0.5 0.5 0.5]);
figure (12)
imshow(rgb,'InitialMagnification','fit')
title('Watershed transform of D')

