function result = Func_CombineCell(dd,bw)
comObj = Func_ObjStruct2Img(dd,bw);
Func_PlotIndex(bw,comObj,dd,2);
flag = 0;
allInd3 = 1:1:length(dd.PixelIdxList)-1;
while flag == 0
    cInd = input('Input an array of Indexes of cells that belong in one cell (0 when done):\n');
        if isempty(cInd)
            flag = 1;
        elseif cInd == 0
            flag = 1;
        elseif length(cInd) == 1
            disp('Error! Input more than one index in a array form to combine')
        else
            for i = 1:length(cInd)
                num = find(allInd3 == cInd(i));
            end
            combine = false(size(bw));
            for k = 1:length(cInd)
                combine = combine + comObj{cInd(k)};
            end
            combine = bwmorph(combine, 'bridge');
            combine = imfill(combine,'holes');
            dd.PixelIdxList{cInd(1)} = find(combine);
            dd.NumObjects = dd.NumObjects -length(cInd)+1;
            for j = 2: length(cInd)
                dd.PixelIdxList(cInd(j)) = [];
            end
            comObj = Func_ObjStruct2Img(dd,bw);
            f5 = Func_PlotIndex(bw,comObj,dd,2);
        end
end
result = dd;
end