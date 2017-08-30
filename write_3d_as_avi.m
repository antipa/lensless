function write_3d_as_avi(arr,fname,cmap_in,frame_rate)
% write_3d_as_avi(arr,fname,cmap_in)
% inputs:
% arr : 3d array to be written. Dim 3 is always assumed to be time.
% fname : string name of the video file. It will be written to ~/[fname], so you
% can use a path here if you want, but it must be relative to your home
% directory.
% cmape_in : string specifying the colormap. Any MATLAB colormap should
% work.
% frame_rate (optional) : video framerate
%
% Outputs: 
% Saves video file, as well as a png of the color bar. The min and max for
% the colorbar are printed to the command prompt.
dots = strfind(fname,'.');
if ~isempty(dots)
    if strcmpi(fname(dots(end):end),'.avi')
        fpath = ['~/',fname];
        fname = fname(1:dots(end)-1);
    else
        error('File extension is not avi')
        return
    end
else
    fpath = ['~/',fname,'.avi'];
end

vout = VideoWriter(fpath);
vout.Quality = 100;
vout.FrameRate = frame_rate;
try
    vout.open;
catch
    error('cannot open file for writing')
    return
end

cmap = colormap(cmap_in);
current_frame = struct('cdata',[],'colormap',[]);


% Scale array to be from 0 to 1
arr_min = min(arr(:));
c_range = max(arr(:))-arr_min;
arr_sc = (arr-arr_min)/c_range;
for n = 1:size(arr,3)
    im_rgb = ind2rgb(floor(1+size(cmap,1)*arr_sc(:,:,n)),cmap);
    current_frame.cdata = uint8(255*im_rgb);
    writeVideo(vout,current_frame);
end
vout.close

c_bar = repmat(linspace(1,size(cmap,1),size(cmap,1))',[1,6]);
c_bar_rgb = ind2rgb(c_bar,cmap);
c_bar_fname = ['~/',fname,'colorbar.png'];
imwrite(uint8(255*c_bar_rgb),c_bar_fname)

fprintf(['file saved to ',fpath,'\n']);
fprintf(['color bar saved as ', c_bar_fname,' range: '])
fprintf('%.3e to %.3e\n',arr_min,arr_min+c_range)


