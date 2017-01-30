clear;
clc;
close all;

% set figure window position
p = [0 0 768 864];
set(0, 'DefaultFigurePosition', p)
%%
% Go to image dir 
Segmentation_Home_dir = 'C:\Users\weida\Desktop\Nuclei_Identification';
Image_File_dir = 'C:\Users\weida\Desktop\Nuclei_Identification\Image_File';
Output_dir = 'C:\Users\weida\Desktop\Nuclei_Identification\Output';

cd(Segmentation_Home_dir)

% make output dir


cd(Image_File_dir)
listing = dir;

nTotalMaxImg = length(listing) - 2;
ImgFolderList = cell(nTotalMaxImg,1);

for i = 1:nTotalMaxImg
    ImgFolderList{i} = listing(2+i).name;
end

for i = 1:nTotalMaxImg
    cd(ImgFolderList{i})
    ImgData = dir;
    I = imread(ImgData(end).name);
    
    cd(Segmentation_Home_dir)
    
    MIP = Func_MIP(ImgData,Image_File_dir,ImgFolderList{i},Segmentation_Home_dir);
   
    [celldata, cellpic, ee] = Func_Segmentation(MIP);
    
    Output = Func_Spacial_Analysis(celldata, cellpic, ee);
    
    cd(Output_dir)
    
    filename = ImgFolderList{i};
    Outputname = strcat(filename,'.xlsx');
    xlswrite(Outputname,Output)
    cd(Image_File_dir)
end
cd(Segmentation_Home_dir)
% 
% OutputFileName