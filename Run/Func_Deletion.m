function result = Func_Deletion(dd,bw)
comObj = Func_ObjStruct2Img(dd,bw);
Func_PlotIndex(bw,comObj,dd,2);
 flag = 0;
    while flag == 0
        delInd = input('Input the index of nucleous to delete, if none, input 0:\n');
        if isempty(delInd)
            flag = 1;
        elseif delInd == 0 
            flag = 1;
        elseif length(delInd) == 1
            if delInd > 0 || delInd < length(dd.PixelIdxList)-1
                dd.PixelIdxList(delInd)=[];
                comObj(delInd,:) = [];
                dd.NumObjects = dd.NumObjects -1;
                Func_PlotIndex(bw,comObj,dd,2);
            else
                disp('Error! Invalid input')
            end
        else
            disp('Error! Invalid input')
        end
    end
    result = dd;
end