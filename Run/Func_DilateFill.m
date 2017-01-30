function result = Func_DilateFill(dd,bw)
comObj = Func_ObjStruct2Img(dd,bw);
    Func_PlotIndex(bw,comObj,dd,2);
    flag = 0;
    count = 0;
    oldInd = -1;
    while flag == 0
        fInd = input('Input the index of nucleous to dilate or fill in the holes, if none, input 0:, input [1 1] to go back to original shape\n');
        if isempty(fInd)
            flag = 1;
        elseif fInd == 0 
            flag = 1;
        elseif length(fInd) == 1
            if fInd > 0 && fInd < length(dd.PixelIdxList)-1
                count = count + 1;
                if fInd ~= oldInd
                    oldimg = comObj{fInd};
                    oldInd = fInd;
                end
                Ind = fInd;
                fillObj = comObj{fInd};
                E = edge(fillObj);
                Ed = imdilate(E,strel('sphere',1));
                fillObj(Ed) = Ed(Ed);
                fillObj = imfill(fillObj,'holes');
                dd.PixelIdxList{Ind} = find(fillObj);
                comObj = Func_ObjStruct2Img(dd,bw);
                Func_PlotIndex(bw,comObj,dd,2);
%                 E = edge(fillObj);
%                 Ed = imerode(E,strel('sphere',5));
%                 fillObj(Ed) = Ed(Ed);
%                 fillObj = imfill(fillObj,'holes');
            else
                disp('Error! Invalid input')
            end
        elseif length(fInd) == 2
            fillObj = oldimg;
            dd.PixelIdxList{Ind} = find(fillObj);
            comObj = Func_ObjStruct2Img(dd,bw);
            Func_PlotIndex(bw,comObj,dd,2);
        else
            disp('Error! Invalid input')
        end
    end
    result = dd;
end