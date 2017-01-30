function Result = Func_Watershed(i,Obj,level)
% There are 2 levels of watershed designed
% level 1 is for first seperating multi cell clusters,
%   this level is gental(do no seperate in to large number of fragments
% level 2 is for for cutting the already seperated nucleus in case the 
%   level 1 result could not provide sufficient segmentation


    % Modification by strel, gaussian smoothing and watershed 
    % for multicell clusters
    
    % figure(11)

    if level == 1
          obj = Obj{i};

        BWnobord = imclearborder(obj, 8);
        seD = strel('disk',1);
        BWfinal = imerode(BWnobord,seD);
        n=1;
        for j = 1:n
        BWfinal = imerode(BWfinal,seD);
        end

        obj = ~BWfinal;
        D = bwdist(obj);
        D = imgaussfilt(D,2);

        D = -D;
        D(obj) = -Inf;
        Result = watershed(D);
    elseif level == 2
        obj = ~Obj{i};
        D = bwdist(obj);
        D = -D;
        D(obj) = -Inf;
        conn = conndef(2,'minimal');
        Result = watershed(D,conn);
    end

end