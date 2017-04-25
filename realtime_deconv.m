input_folder = 'Y:\Diffusers''nstuff\2d_images_to_process';
camera_type = 'pco';
process_color = 'rgb';
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
            im_found = 1;   %Escape if an image is founmd.
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
out_dir = [out_base,'\',time_stamp,'_',im_base,'_',process_color];

res_dir = out_dir;
raw_dir = out_dir;

bin = (imread(im_to_move));
switch lower(camera_type)
    case('pco')
        if numel(size(bin))==2
            b_dem =  demosaic(bin,'rggb');
        elseif numel(size(bin))==3
            b_dem = bin;
        end
        switch lower(process_color)
            case('mono')
                bstack = mean(double(b_dem),3);
            case('red')
                bstack = double(b_dem(:,:,1));
            case('green')
                bstack = double(b_dem(:,:,2));
            case('blue')
                bstack = double(b_dem(:,:,3));
            case('rgb')
                bstack = double(b_dem);
        end
    case('flea3')
        bstack = double(bin);
end
        
%bin = mean(double(b_dem),3);


% figure(2),clf
% imagesc(bin)
% axis image
% colormap gray
%%
%xhat_save = zeros(2*size(bstack,1),2*size(bstack,2),size(bstack,3));
for n = 1:size(bstack,3)
    if size(bstack,3)>1
        if n == 1
            colors = 'red';            
        elseif n == 2
            colors = 'green';
        elseif n == 3
            colors = 'blue';
        end
    else
        colors = process_color;
    end
           
        
    bin = bstack(:,:,n);
    diffuser_2d_deconv_v2
    if n == 1
        xhat_save = (gather(xhat));
    else
        xhat_save = cat(3,xhat_save,gather(xhat));
    end
end



mkdir(out_dir)
movefile(im_to_move,raw_dir)
options.fighandle = [];
save([res_dir,'\',time_stamp,'_',im_base,'_processed.mat'],'xhat_save');
imout = uint8(255*max(xhat_save/prctile(xhat_save(:),99.99),0).^(1/1.5));
imwrite(imout,[res_dir,'\',time_stamp,'_',im_base,'_processed.png']);
