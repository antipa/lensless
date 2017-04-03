imstack_dir = 'Y:\Diffusers''nstuff\video-grace';
process_dir = 'Y:\Diffusers''nstuff\2d_images_to_process';
out_base = 'Y:\Diffusers''nstuff\Processing_results';


movie_time_stamp = datestr(now,'yyyymmdd_HHMMSS');
movie_out_dir = [out_base,'\',movie_time_stamp,'_scrollMovie'];
mkdir(movie_out_dir)
v = VideoWriter([movie_out_dir,'\',movie_time_stamp,'scroll.avi'],'Uncompressed AVI');
v.FrameRate = 15;

files = dir(imstack_dir);
fnum = [];
c = 0;
for n = 3:length(files)
    unders = findstr(files(n).name,'_');
    dots = findstr(files(n).name,'.');
    ext = files(n).name(dots(1)+1:end);
    if strcmpi(ext,'png')
        c = c+1;
        fnum(c) = str2double(files(n).name(unders(end)+1:(dots(1)-1)));  
    end
end
    
[~,I] = sort(fnum);
I = I+2;
open(v);

for m = 1:length(I)
    file_to_process = files(I(m)).name;
    copyfile([imstack_dir,'\',file_to_process],process_dir);
    pause(.5)  %Make sure all files are done moving
    realtime_deconv;   %Returns image imout rgb
    pause(.5)   
    cm = round(.35*size(imout)+1);
    cmx = round(.65*size(imout));
    writeVideo(v,imout(cm(1):cmx(1),cm(2):cmx(2),:))
end
close(v)
    