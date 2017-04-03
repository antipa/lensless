imstack_dir = 'Y:\Diffusers''nstuff\video-grace';
process_dir = 'Y:\Diffusers''nstuff\2d_images_to_process';
out_base = 'Y:\Diffusers''nstuff\Processing_results';

im_name_base = 'scroll_video';

folders = dir(out_base)
c=0;
i = [];
dt = [];
for n = 1:numel(folders)
   if findstr(folders(n).name,im_name_base)
      
       
       folders(n).name;
       unders = strfind(folders(n).name,'_');
       num_str = [folders(n).name(1:unders(1)-1),folders(n).name(unders(1)+1:unders(2)-1)];
       dtp = str2double(num_str);
       if dtp >= 20170329124121
           c = c+1;
           i(c) = n;
           
           dt(c) = dtp;
       end
           
   end
end


movie_time_stamp = datestr(now,'yyyymmdd_HHMMSS');
movie_out_dir = [out_base,'\',movie_time_stamp,'_scrollMovieRemake'];
mkdir(movie_out_dir)
pause(1)
v = VideoWriter([movie_out_dir,'\',movie_time_stamp,'scroll.avi'],'Uncompressed AVI');
v.FrameRate = 15;
open(v)
for n = 1:length(dt)
    imname = folders(i(n)).name;
    unders = strfind(imname,'_');
    pngname = [imname(1:unders(end)-1),'_processed.png'];
    imin = imread([out_base,'\',imname,'\',pngname]);
     cm = round(.35*size(imin)+1);
    cmx = round(.65*size(imin));
    imcrop_rgb = imin(cm(1):cmx(1),cm(2):cmx(2),:);
    imcrop_tv = imcrop_rgb;%tvdenoise3d_diff(double(imcrop_rgb),.3,50,0);
     writeVideo(v,uint8(imcrop_tv))
    %imagesc(imcrop_rgb)
    n
end
close(v)