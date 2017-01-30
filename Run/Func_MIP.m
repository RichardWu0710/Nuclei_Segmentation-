function MIP = Func_MIP(ImgList,Image_File_dir,Imgdir,run_dir)
% Create Maximum Intensity Image from a input stack of images
    
    cd(Image_File_dir)
    cd(Imgdir)
    
    nImg = length(ImgList);
   
    I = imread(ImgList(end).name);
    SIZE = size(I);
    ComboImg = zeros(SIZE(1),SIZE(2),nImg);
    for j = 1:nImg
        ComboImg(:,:,j) = imread(ImgList(j).name);
    end
    
    MIP = max(ComboImg,[],3);
    MIP = uint8(MIP);
    
    cd(run_dir)

end