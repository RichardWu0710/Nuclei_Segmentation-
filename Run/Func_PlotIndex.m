function f5 = Func_PlotIndex(bw,tempObj,gg,handle)
% this function is for plotting the index number on individual nuclei
% handle =  1 when the first obj of gg is the combination of all the rest
% objects
% handle = 2 when the first obj of gg is the actual first obj

tempic1 = zeros(size(bw));
        for j = 1:length(tempObj(:,1))
            tempic1 = tempic1+tempObj{j,1};
        end
        f5 = figure (5);
            imshow(tempic1);
        s = regionprops(gg, 'Centroid');
        imshow(tempic1)
        hold on
        for k = 1:numel(s)
            c = s(k).Centroid;
            text(c(1), c(2), sprintf('%d', k-2+handle), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle');
        end
        hold off


end