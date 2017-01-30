function Result = Func_ObjStruct2Img(dd,bw)


n = length(dd.PixelIdxList);
finalSingleObj = cell(n);
for i = 1:length(dd.PixelIdxList)
    pic = false(size(bw));
    pic(dd.PixelIdxList{i}) = true;
    finalSingleObj{i}= pic;
end
Result = finalSingleObj;
end