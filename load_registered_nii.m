% Load all Niftifiles in the given folder and save them in one variable
% called Im_all. 

d = dir(img_folder);
dsize = length(d);

% Select Nifti files in the folder
% filenames = [];
cnt = 0;
for i = 1 : dsize
    filename_tmp = d(i).name;
    isnii = strfind(filename_tmp, '.nii');
    if ~isempty(isnii)
        cnt = cnt + 1;
        filenames{cnt} = filename_tmp;
    end
end

info = niftiinfo([img_folder filenames{1}]);
Nim = size(filenames, 2);

% Load Nifti files
Im_all = [];
for i = 1 : Nim
    Im_all(:,:,:,i) = niftiread([img_folder filenames{i}]);
end 

% Check if they are registered
ctr_sl = round(size(Im_all,3)/2); % Define Central Slice to be visualized
figure
for i = 1 : Nim - 1
    imshowpair(Im_all(:,:,ctr_sl,i), Im_all(:,:,ctr_sl,i + 1))
    title(['Compare Im N° ' num2str(i) ' with Im N° ' num2str(i + 1) ])
    pause
end

