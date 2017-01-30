wholepic1 = zeros(size(bw));
for i = 1:n
    wholepic1 = wholepic1+allimage{i};
end
figure (2)
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