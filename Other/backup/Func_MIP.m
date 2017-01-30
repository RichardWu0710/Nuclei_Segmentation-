function MIP = Func_MIP(ImgList,Image_File_dir,Imgdir,Home_dir)
% Create Maximum Intensity Image from a input stack of images
    
    cd(Image_File_dir)
    cd(Imgdir)
    
    nImg = length(ImgList)-2;
    
    I = imread(ImgList(end).name);
    SIZE = size(I);
    ComboImg = zeros(SIZE(1),SIZE(2),nImg);
    for j = 1:nImg
        ComboImg(:,:,j) = imread(ImgList(j+2).name);
    end
    
    MIP = max(ComboImg,[],3);
    MIP = uint8(MIP);
    
    cd(Home_dir)

end