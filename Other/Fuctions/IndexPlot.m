function IndexPlot(struct, pic)
% Print Index in the middle of the cell

gg = struct;
tempic1 = pic;
    imshow(tempic1);
    s = regionprops(gg, 'Centroid');
    imshow(tempic1)
    hold on
    for k = 1:numel(s)
        c = s(k).Centroid;
        text(c(1), c(2), sprintf('%d', k-1), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle');
    end
    hold off

end