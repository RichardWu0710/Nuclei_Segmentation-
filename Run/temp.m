
% input 
% MIP = Max intensity image
% res = the resolution of MIP


% sample name
% MAX_12-5-2016 Position2.tif
% 'Sample.jpg'

I = MIP;
figure (1)
imshow(I);
title('Original Image')

% Ajust Brightness
I = I -25;
I = I*5;
bw = imbinarize(I);
for i = 1:2
bw = imfill(bw,'holes');
end
E = edge(bw);
    Ed = imdilate(E,strel('sphere',1));
    bw(Ed) = Ed(Ed);
    bw = imfill(bw,'holes');
% figure(2)
% imshow(bw)
%%
% Adjust Contrast
% I = imadjust(I); % Adjust
% figure (3)
% imshow(I3);

% Turn image into binary format
% Adjust https://www.mathworks.com/help/images/ref/imbinarize.html
% bw = imbinarize(I);       

% Remove noise
bw = bwareaopen(bw, 4);
cc = bwconncomp(bw, 18);

% allimage contains all matrix for each crude cell object
n = length(cc.PixelIdxList);
allimage = cell(n,2);
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
f2 = figure (2);
imshow(wholepic1);
s = regionprops(cc, 'Centroid');
imshow(wholepic1)
hold on
for k = 1:numel(s)
    c = s(k).Centroid;
    text(c(1), c(2), sprintf('%d', k), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle');
end
hold off
title('All Objects After Initial Filtering')


% plot  area value for each cell
% can be used to select threshold
%
areabar = zeros(n,1);
for i=1:n
    areabar(i)= allimage{i,2}.Area;
end

f3 = figure (3);
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


% seperate and plot cells based on area threshold
% multiObj, singleObj, debrisObj contain conditionaly sorted cell images

threshold = [0.15*10^4 700];
% Ask user input for debris area threshold
%debrislim define area limit bellow which is debri
threstr = int2str(threshold(res));



Question = ['Area threshold for debris? default:' threstr ',type 0 to skip entire position\n'];
debrislim = input(Question);

% minimize previous figure windows
set(f2,'Visible', 'off'); 
set(f3,'Visible', 'off'); 
%%
if debrislim == 0
    err = 1;
    celldata = 0;
    cellpic = 0;
    ee = 0;

else
    err = 0;
    if isempty(debrislim)
         debrislim = threshold(res);
    end

    % find object indexes from the allimage array that
    debrisObjInd = find(areabar<=debrislim);

    % number of cells for each condition
    ndebris = length(debrisObjInd);

        allimage(debrisObjInd,:)=[] ;
        cc.PixelIdxList(debrisObjInd)=[];
        cc.NumObjects=cc.NumObjects-ndebris;

    Func_PlotIndex(bw,allimage,cc,2)

    % ask to input the index of multi cell
    flag = 0;
    multiObjInd = [];
    while flag == 0
        num = input('Input the index of multicell cluster,one at a time. if no, input 0:\n');
        if isempty(num)  
             flag = 1;
        elseif num <= length(allimage) && num >= 1
            multiObjInd = [multiObjInd num];
        elseif length(num) > 1 || num > length(allimage)
            disp('Error! Invalid input')
        
        elseif num == 0
          flag = 1;
        else
            disp('Error! Invalid input')
        end
    end

    nmulti = length(multiObjInd);
    multiObj = cell(nmulti);

    for i = 1:nmulti
        multiObj{i} = allimage{multiObjInd(i)};
        multiObj{i,2} = regionprops(multiObj{i},'Area');
    end

    allInd = 1:1:length(allimage);
    singleObjInd = setdiff(allInd,multiObjInd);
    nsingle = length(singleObjInd);
    singleObj = cell(nsingle,2);
    for i = 1:nsingle
        singleObj{i} = allimage{singleObjInd(i)};
        singleObj{i,2} = regionprops(singleObj{i},'Area');
    end


    % Modification by strel, gaussian smoothing and watershed 
    % for multicell clusters
    newMultiObj = cell(nmulti);


    for i = 1:nmulti

        newMultiObj{i}= Func_Watershed(i,multiObj,1);

    end
    intermediat = singleObj;


    % Add newly seperated single cells to singleObj array and their areas
    nnew = nsingle;

    % reset single Obj incase this section false and have to rerun the section
    singleObj = intermediat;
    

    for i = 1: length(newMultiObj)
        %%
        bw = bwareaopen(newMultiObj{i}, 50);
        gg = bwconncomp(bw, 26);
        multiple = cell(length(gg.PixelIdxList)-1);
        tempObj = cell(length(gg.PixelIdxList)-1,1);
         for add = 1:length(gg.PixelIdxList)-1
            pic = false(size(newMultiObj{i}));
            pic(gg.PixelIdxList{add+1}) = true;
            tempObj{add}= pic;
         end

         segobj= tempObj;
         % image of watershed result with index
         Func_PlotIndex(bw,tempObj,gg,1);
       %%
        stop = 0;
        while stop ==0 
            good = input('Seperation looks correct? 1 for Yes. 2 for No\n');
            if isempty(good)
                freshsingle = tempObj;
                stop = 1;
            elseif good ==1
                freshsingle = tempObj;
                stop = 1;
            elseif good == 2
                stop = 1;
                flag = 0;
                count = 0;

                % keep segmenting user input

                while flag == 0
                Ind = input('Input index of cell to seperate (0 when done):\n');
                    if isempty(Ind)
                        flag = 1;
                    elseif length(Ind) >1
                        disp('Error! Seperation step requires one input at a time')  
                    elseif Ind == 0
                        flag = 1;
                    else
                        seg = Func_Watershed(Ind,segobj,2);
                        bw = bwareaopen(seg, 50);
                        ss = bwconncomp(bw, 4);
                        gg.PixelIdxList(Ind+1)=[];
                        gg.PixelIdxList=[gg.PixelIdxList ss.PixelIdxList(2:end)];
                        gg.NumObjects = length(gg.PixelIdxList);
                        tempObj = cell(length(gg.PixelIdxList)-1,1);
                         for add = 1:length(gg.PixelIdxList)-1
                            pic = false(size(newMultiObj{i}));
                            pic(gg.PixelIdxList{add+1}) = true;
                            tempObj{add}= pic;
                         end


                         % image of watershed result with index

                        f5 = Func_PlotIndex(bw,tempObj,gg,1);
                    end
                end

                % Stare Combining fragments
                flag = 0;
                allInd2 = 1:1:length(tempObj);
                while flag == 0
                disp(['There are ' num2str(length(allInd2)) ' fragments presented']);
                Ind = input('Input an array of Indexes of fragments that belong in one cell (0 when done):\n');
                    if isempty(Ind)
                        flag = 1;
                    elseif Ind == 0
                        flag = 1;
                    elseif length(Ind) == 1
                        disp('Error! Input more than one index in a array form to combine')
                    else
                        for i = 1:length(Ind)
                            num = find(allInd2 == Ind(i));
                            allInd2(num) = [];
                        end
                        count = count+1;
                        combine = false(size(bw));
                        for k = 1:length(Ind)
                            combine = combine + tempObj{Ind(k)};
                        end
                        combine = bwmorph(combine, 'bridge');
                        combine = imfill(combine,'holes');
                        multiple{count} = combine;
                    end
                end

                nrun = length(allInd2)+count;

                freshsingle = cell(nrun);
                for n = 1:count
                    freshsingle{n} = multiple{n};
                end

                for n = 1+count : nrun
                    freshsingle{n} = tempObj{allInd2(n-count)};
                end
            else
                disp('Error! Invalid input')
            end
        end


        oldSingle = singleObj;                  % temperory storage 
        cellsize = length(singleObj)+length(freshsingle);
        singleObj = cell(cellsize);             % create longer cell array
        % for loop record all old values back
        for j = 1: length(oldSingle)            
            singleObj{j} = oldSingle{j};
            singleObj{j,2} = oldSingle{j,2};
        end
        % for loop add new values
        for add = 1:length(freshsingle)
            newind = nnew+add;
            singleObj{newind}= freshsingle{add};
            singleObj{newind,2} = regionprops(pic,'Area');  % Compute Area
        end
        nnew = nnew + length(freshsingle);        

     end
%%
    finalSingleObj = singleObj;

    % Plot all single cells together after watershed
    singleObjpic = false(size(bw));
    for i = 1:length(finalSingleObj)
        singleObjpic = singleObjpic+finalSingleObj{i};
    end
    % figure (15)
    % imshow(singleObjpic);
    % title('All Single Cells after Watershed')


    % Restruct all single cells into a final array for analyzing use
    % Create global variables to use in different files
    dd = bwconncomp(singleObjpic,4);
    f5 = Func_PlotIndex(bw,finalSingleObj,dd,2);
    
    

%%   
    % adjustment option dilate/fill and combine single cell

        flag = 0;
    while flag == 0
        opt = input('Need further adjustment? Input 1 for deletion, 2 for dilation/fill holes, 3 for combining single nucleus\n, input 0 when done');
        if isempty(opt) 
            flag = 1;
        elseif opt == 0
            flag = 1;
        elseif opt == 1
            dd = Func_Deletion(dd,bw);
        elseif opt == 2
            dd = Func_DilateFill(dd,bw);
        elseif opt == 3
            dd = Func_CombineCell(dd,bw);
        else
            disp('Error! Invalid input')
        end
    end
   
    % allimage contains all matrix for each crude cell object

    finalSingleObj = Func_ObjStruct2Img(dd,bw);


    singleObjpic = false(size(bw));
    for i = 1:length(finalSingleObj)
        singleObjpic = singleObjpic+finalSingleObj{i};
    end

    celldata = finalSingleObj;

    cellpic = singleObjpic;

    ee = dd;
    set(f5,'Visible', 'off'); 
end



