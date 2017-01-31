input_folder = 'Y:\Diffusers''nstuff\2d_images_to_process';

%Scan input_folder for an image
im_found = 0;
k = 0;
files_to_delete = cell(1);
while ~im_found
    indir = dir(input_folder);
    for n = 3:length(indir)
        fname = [input_folder,'\',indir(n).name];
        try 
            %im_data = imread(fname);   %Image data
            im_to_move = fname;    %Path to image before moving
            im_name = indir(n).name;  %Filename
            im_found = 1;   %Escape if an image is found.
        catch
            k = k+1;
            fprintf('non image file found. Deleting.')
            files_to_delete{k} = fname;
            im_found = 0;
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
out_dir = [out_base,'\',time_stamp,'_',im_base];

res_dir = out_dir;
raw_dir = out_dir;

bin = (imread(im_to_move));
b_dem =  demosaic(bin,'bggr');
bin = mean(double(b_dem),3);

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
xhat_save = crop(gather(xhat));

save([res_dir,'\',time_stamp,'_',im_base,'_processed.mat'],'xhat','options');
imwrite(uint8((max(xhat_save,0)/max(xhat_save(:))).^(1/2)*255),[res_dir,'\',time_stamp,'_',im_base,'_processed.png']);
