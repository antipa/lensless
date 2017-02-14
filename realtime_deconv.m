input_folder = 'Y:\Diffusers''nstuff\2d_images_to_process';
camera_type = 'pco';
colors = 'green';
%Scan input_folder for an image
im_found = 0;
k = 0;
files_to_delete = cell(1);
while ~im_found
    indir = dir(input_folder);
    for n = 3:length(indir)
        fname = [input_folder,'\',indir(n).name];
        try 
            im_data = imread(fname);   %Image data
            im_to_move = fname;    %Path to image before moving
            im_name = indir(n).name;  %Filename
            im_found = 1;   %Escape if an image is found.
        catch
            k = k+1;
            fprintf('non image file found. Deleting.\n')
            files_to_delete{k} = fname;

        end
    end
    
    pause(1)
end

%Erase non-image files
for n = 1:length(files_to_delete)
   delete(files_to_delete{n});
end

dots = findstr(im_name,'.');
im_base = im_name(1:dots-1);
% Get system time
time_stamp = datestr(now,'yyyymmdd_HHMMSS');
out_base = 'Y:\Diffusers''nstuff\Processing_results';
out_dir = [out_base,'\',time_stamp,'_',im_base,'_',colors];

res_dir = out_dir;
raw_dir = out_dir;

bin = (imread(im_to_move));
switch lower(camera_type)
    case('pco')
        b_dem =  demosaic(bin,'rggb');
        switch lower(colors)
            case('mono')
                bin = mean(double(b_dem),3);
            case('red')
                bin = double(b_dem(:,:,1));
            case('green')
                bin = double(b_dem(:,:,2));
            case('blue')
                bin = double(b_dem(:,:,3));
        end
    case('flea3')
        bin = double(bin);
end
        
%bin = mean(double(b_dem),3);


figure(2),clf
imagesc(bin)
axis image
colormap gray
%%
diffuser_2d_deconv_v2
%%
mkdir(out_dir)
movefile(im_to_move,raw_dir)
%%
xhat_save = (gather(xhat));
options.fighandle = [];
save([res_dir,'\',time_stamp,'_',im_base,'_processed.mat'],'xhat_save');
imwrite(uint8((max(xhat_save,0)/max(xhat_save(:))).^(1/1.6)*255),[res_dir,'\',time_stamp,'_',im_base,'_processed.png']);
