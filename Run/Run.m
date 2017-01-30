clear;
clc;
close all;

% set figure window position
p = [0 0 768 864];
set(0, 'DefaultFigurePosition', p)
%%
% Go to image dir 
% Segmentation_Home_dir = 'C:\Users\weida\Desktop\Nuclei_Identification';
% Image_File_dir = 'C:\Users\weida\Desktop\DARTFISH_Example\DRAQ5_Restain';
% Output_dir = 'C:\Users\weida\Desktop\Nuclei_Identification\Output';
% RNA_Data_dir = 'C:\Users\weida\Desktop\Nuclei_Identification\RNA_Data';

Run_dir = pwd;
cd ../
Segmentation_Home_dir = pwd;
cd (Run_dir)
Image_File_dir = input('Input path to the directory containing images from all the positions\n','s');
RNA_Data_dir = input('Input path to the directory containing data of gene''s x and y positions from all the positions\n','s');

startpos = input('what image position to start? default input is 1');
if isempty(startpos)
    startpos = 1;
end
Output_dir = Segmentation_Home_dir;
cd(Output_dir)
mkdir('Output')
mkdir('Output','Gene_Projection_Data')
mkdir('Output','MIP')
mkdir('Output','Gene_Projection_Images')
cd Output
Output_dir = pwd;
cd(Segmentation_Home_dir)

cd(Image_File_dir)
listing = dir('*Pos*');

% Ask for user input on channel and resolution for the experiment
    flag = 0;
    while flag == 0
        Chan = input('what is the channel number (1 2 or 3) used to aquaire the image? default is 1\n');
        if isempty(Chan) 
            Channel = 'ch01';
            flag = 1;
        elseif Chan == 1 
            Channel = 'ch01';
            flag = 1;
        elseif Chan == 2
            Channel = 'ch02';
            flag = 1;
        elseif Chan == 3
            Channel = 'ch03';
            flag = 1;
        else
            disp('Error! Invalid input')
        end
    end
    
% Ask user to choose resolution of the image
disp('What is the resolution of the image? ')
flag = 0;
while flag == 0
    choice = input('input 1 for 2048 by 2048, input 2 for 1024 by 1024\n','s');
    if choice == '1'
        res = 1;
        flag = 1;
    elseif choice == '2'
        res = 2;
        flag = 1;
    else
        disp('Error! Invalid input')
    end
end

% Prep for the main loop 
nTotalMaxImg = length(listing);
ImgFolderList = cell(nTotalMaxImg,1);

for i = 1:nTotalMaxImg
    ImgFolderList{i} = listing(i).name;
end

ComboCount = 0;
for i = startpos:nTotalMaxImg
    cd(ImgFolderList{i})
    % use while loop so user has option to reapeat on a position  
    while true
        filename = ImgFolderList{i};

        disp(['You are on ' filename])

        cellpos = filename(end-2:end);
        posnum = str2num(cellpos);


        Imgname = strcat('*',Channel,'*');
        ImgData = dir(Imgname);
        I = imread(ImgData(end).name);

        cd(Run_dir)

        MIP = Func_MIP(ImgData,Image_File_dir,ImgFolderList{i},Run_dir);

        [celldata, cellpic, ee,err] = Func_Segmentation(MIP,res);
        
        if err == 1
            errorsign = strcat({'You have skipped'},{' '},{filename},{'\n'});
            disp(errorsign);
            break
        else
            cd(RNA_Data_dir)
            temp = int2str(posnum);
            RNA_Finder = strcat('*Pos',temp,'*');
            RNA_Data_Name = dir(RNA_Finder);
            RNA_Data = importdata(RNA_Data_Name(end).name);
            cd(Run_dir)
            [Output,err,fig] = Func_Spacial_Analysis(celldata,cellpic,ee,posnum,RNA_Data,res);
            if err ==1 
                errorsign = strcat({'You have skipped'},{' '},{filename});
                disp(errorsign);
                break
            elseif err ==2
                errorsign = strcat({'You are repeating'},{' '},{filename});
                disp(errorsign);
                cd(Image_File_dir)
                cd(ImgFolderList{i})
                % break the while loop to repeat on this position
            else
                if ComboCount == 0
                    Combo_Output = Output;
                else
                    Combo_Output = Combo_Cell(Combo_Output,Output);
                end
                ComboCount = ComboCount +1;
                cd([Output_dir '/Gene_Projection_Data'])
                Outputname = strcat(filename,'.xlsx');
                xlswrite(Outputname,Output)

                cd([Output_dir '/Gene_Projection_Images'])
                p = [0 0 768 864];
                set(0, 'DefaultFigurePosition', p)
                saveas(fig,[filename '.png'])
                break
            end
        end
    end
    cd([Output_dir '/MIP'])
        imwrite(MIP,[filename '.tif'])
    cd(Image_File_dir)
end
%%
cd([Output_dir '/Gene_Projection_Data'])
xlswrite('Combine_Cell_Stats.xlsx',Combo_Output)
cd(Segmentation_Home_dir)
