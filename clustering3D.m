% Run 3D clustering on previously loaded DCE images. 

global alea sigma;
alea = 1;

folder_out = [img_folder 'output/'];
if ~isfolder(folder_out)
    mkdir(folder_out)
end

sizeBorders = 4; % Size of borders to be ignored (in pixels)
v = - abs(2*t_low); % value applied to all pixels to be ignored
[nrows,ncols,nim,nph] = size(Im_all);

% show the central slices of the four phases of the raw data
if plotting
    figure
    for i = ctr_sl : ctr_sl + 1
        subplot(2,2,1);
        imagesc(Im_all(:,:,i,1)), colormap(gray), title(['Non injected phase - Slice N ' num2str(i)])
        subplot(2,2,2);
        imagesc(Im_all(:,:,i,2)), colormap(gray), title(['Arterial phase - Slice N ' num2str(i)]), 
        subplot(2,2,3);
        imagesc(Im_all(:,:,i,3)), colormap(gray), title(['Portal  phase - Slice N ' num2str(i)]), 
        subplot(2,2,4);
        imagesc(Im_all(:,:,i,4)), colormap(gray), title(['Late  phase - Slice N ' num2str(i)]), 
        sgtitle('Raw data') 
        pause
    end
end


%% Create matrix of x coordinates
for i = 1 : nrows
	for j = 1 : ncols
		for k = 1 : nim
			Im_coord_y(i,j,k) = i;
			Im_coord_x(i,j,k) = j;
			Im_coord_z(i,j,k) = k;
		end
	end
end

Im_coord = [Im_coord_x(:), Im_coord_y(:), Im_coord_z(:)];
Im_coord = norm_cols(Im_coord);

%% Images subsampling
subsampled = Im_all;
size_spx = 4; % size of 'superpixels'
for p = 1 : 4	
	for i = 0 : floor(nrows/size_spx) - 1
		for j = 0 : floor(ncols/size_spx) - 1
			for k = 0 : floor(nim/size_spx) - 1
				tmp = mean(Im_all((size_spx*i)+1 : size_spx*(i+1), (size_spx*j)+1 : size_spx*(j+1), (size_spx*k)+1 : size_spx*(k+1) , p)); 
				subsampled((size_spx*i)+1 : size_spx*(i+1), (size_spx*j)+1 : size_spx*(j+1), (size_spx*k)+1 : size_spx*(k+1), p) = mean(tmp(:));
			end
		end
	end
end

% Plot image subsampled in "superpixels"
if plotting
    figure
    for i = ctr_sl : ctr_sl + 1
        subplot(2,2,1);
        imagesc(subsampled(:,:,i,1)), colormap(gray), title(['Non injected phase - Slice N ' num2str(i)])
        subplot(2,2,2);
        imagesc(subsampled(:,:,i,2)), colormap(gray), title(['Arterial phase - Slice N ' num2str(i)])
        subplot(2,2,3);
        imagesc(subsampled(:,:,i,3)), colormap(gray), title(['Portal phase - Slice N ' num2str(i)])
        subplot(2,2,4);
        imagesc(subsampled(:,:,i,4)), colormap(gray), title(['Late phase - Slice N ' num2str(i)])
        sgtitle('Subsampled images')
        pause
    end
end


%% Compute matrix giving the average of the neighbours of every pixel, without the central pixel
if dif_avgnbr
	% Matrix with average of neighbours
	Im_all_avgnbr = Im_all; % not normalized yet
	n_nbr = 26;
	for p = 1 : 4
		for i = 2 : nrows - 1
			for j = 2 : ncols - 1
				for k = 2 : nim - 1
					sum_nb_k = sum(Im_all(i-1:i+1, j-1, k, p)) + sum(Im_all(i-1:i+1, j+1, k, p)) + Im_all(i-1, j, k, p) + Im_all(i+1, j, k, p);
					nb_kplus = Im_all(i-1:i+1, j-1:j+1, k+1, p);
					nb_kles = Im_all(i-1:i+1, j-1:j+1, k-1, p);
					Im_all_avgnbr(i, j, k, p) = (sum_nb_k + sum(nb_kplus(:)) + sum(nb_kles(:))) / n_nbr ;
				end
			end
		end
	end
	
	% Compute difference between current pixel and average of neighbours and then convert to vector
	Im_all_dif_avgnbr = Im_all - Im_all_avgnbr;
	for i = 1 : 4
		Im_all_dif_avgnbr_vect(:,i) = reshape(Im_all_dif_avgnbr(:,:,:,i), [nrows*ncols*nim,1]);
	end
	
	save([folder_out, 'Im_all_dif_avgnbr_vect'], 'Im_all_dif_avgnbr_vect');
else
	load ([folder_out, 'Im_all_dif_avgnbr_vect']);
end


%% Median filter (5 neighbours)
% takes the median value in a mf-by-mf neighbourhood
mf = 5;	% 3
for i = 1:4
	Im_all(:,:,:,i) = medfilt3(Im_all(:,:,:,i), [mf,mf,mf], 'symmetric'); %symmetric, otherwise it gives 0 on the corners
end

%% Compute time derivatives on enhancement curves 
% Image of Derivatives
for i = 1:3
	Im_Der(:,:,:,i) = (Im_all(:,:,:,i+1) - Im_all(:,:,:,i)) ./ (t(i+1) - t(i));
	Der_vect(:,i) = reshape(Im_Der(:,:,:,i), [nrows*ncols*nim,1]); % convert to vector
% 	figure 
% 	imagesc(Im_Der(:,:,ctr_sl,i)), title(['Derivative image at time ' num2str(i)])
end

%% Identify pixels to be removed using thresholds on intensities and derivatives
if Im_all_preproc
	Im_all = removeBones(Im_all, t_high, v);
	Im_all = removeNonEnhTissues(Im_all, Im_Der, v);
	
	for i = 1:nrows
		for j = 1:ncols
			for p = 1 : nim
				if max(Im_all(i,j,p,:)) > t_high
% 				if Im_all(i,j,p,4) > t_high % only on late phase
					Im_all(i,j,p,:) = v; % sature bones
				elseif min(Im_all(i,j,p,:)) < t_low && Im_all(i,j,p,1) ~= v
					Im_all(i,j,p,:) = v; % sature dark zones
				end
			end
		end
	end
		
	save([folder_out 'Im_all_pretr'], 'Im_all');
else
	load([folder_out 'Im_all_pretr']);
end

% Show all phases of preprocessed images
if plotting
    figure
%   for i = sizeBorders  : nim - sizeBorders
    for i = ctr_sl : ctr_sl + 1
        subplot(2,2,1);
        imagesc(Im_all(:,:,i,1)), colormap(gray), title(['Non injected phase - Slice N ' num2str(i)])
        subplot(2,2,2);
        imagesc(Im_all(:,:,i,2)), colormap(gray), title(['Arterial phase - Slice N ' num2str(i)])
        subplot(2,2,3);
        imagesc(Im_all(:,:,i,3)), colormap(gray), title(['Portal phase - Slice N ' num2str(i)])
        subplot(2,2,4);
        imagesc(Im_all(:,:,i,4)), colormap(gray), title(['Late phase - Slice N ' num2str(i)])
        impixelinfo
        sgtitle('Preprocessed images') 
        pause
    end
end


%% Transform images and subsampled images into vectors and remove borders
for p = 1:4
	Im_vect(:,p) = reshape(Im_all(:,:,:,p), [nrows*ncols*nim,1]);
	subsampled_vec(:,p) = reshape(subsampled(:,:,:,p), [nrows*ncols*nim,1]);
	% Remove borders to Im_vect and to subsampled_vec for pre-clustering
	Im_vect(:,p) = removeBorders3(Im_vect(:,p), sizeBorders, nrows, ncols, v);
	subsampled_vec(:,p) = removeBorders3(subsampled_vec(:,p), sizeBorders, nrows, ncols, v);
end

%% Indexes of pixels not to be taken for clustering
idx_nonv = find(Im_vect(:,1) ~= v);

% Remove pixels to be ignored (those equal to v)
Im_vec_cut = Im_vect(idx_nonv, :);
subsampled_vec_cut = subsampled_vec(idx_nonv, :);
Im_coord_cut = Im_coord(idx_nonv,:);


%% Normalize images and subsampled images to mu=0 and std=1
% for pre-clustering :
Im_vec_cut_norm = norm_cols(Im_vec_cut); 
subsampled_vec_cut = norm_cols(subsampled_vec_cut); 


%% preliminary k-means clustering on subsampled images
prelK = 1;
if prelK
	sigma = 6;
	n_clusters = 5; % 14
	[idx_p, C_p] = kmeans([Im_coord_cut, subsampled_vec_cut], n_clusters, 'MaxIter', 200);
% 	[idx_p, C_p] = kmeans3([Im_coord_cut, subsampled_vec_cut], n_clusters, 'MaxIter', 200); %'Start', 'cluster'
	
	for p = 1 : 4
		clusters_on_subsampled_vec(:,p) = Im_vec_cut_norm(:, p) - C_p(idx_p(:), p + 3); % + 3 if kmeans3 with Coord
	end
 	
	% To visualise clusters on subsampled images:
	idx_p_restore = zeros(size(Im_vect(:,1)));
	idx_p_restore(idx_nonv) = idx_p; 
	clusters_on_subsampled_vis = reshape(idx_p_restore,nrows,ncols,nim); % Show cluster Index
	figure
	for i = 1 : (nim)/size_spx
		imagesc(clusters_on_subsampled_vis(:,:,i*size_spx)), colormap(parula), colorbar % colorcube line prism
		title(['Preliminary K-means with N=' num2str(n_clusters) ' clusters - based on subsampled images']); 
		impixelinfo
		pause
	end
	save([folder_out 'clusters_on_subsampled_vec'], 'clusters_on_subsampled_vec');
else
	load([folder_out 'clusters_on_subsampled_vec']);
end


%% Distance to image center
Ctr = [round(ncols/2), round(nrows/2), round(nim/2)]; 

Im_distToCtr_vect = (Im_coord_x(:) - Ctr(1)).^2;
Im_distToCtr_vect = Im_distToCtr_vect + (Im_coord_y(:) - Ctr(2)).^2;
Im_distToCtr_vect = Im_distToCtr_vect + (Im_coord_z(:) - Ctr(3)).^2;


%% %% Remove pixels equal to v for all features  %%
 % (save original size for later)
idx_restore_size = zeros(size(Im_vect(:,1)));
 % ! size of vectors changes !
Im_vect = Im_vect(idx_nonv, :);
Der_vect = Der_vect(idx_nonv, :);
Im_all_dif_avgnbr_vect = Im_all_dif_avgnbr_vect(idx_nonv, :);
Im_distToCtr_vect = Im_distToCtr_vect(idx_nonv);

% Normalise remaining pixels (without all pixels equal to v)
Im_vect = 2*norm_cols(Im_vect); 
Der_vect = norm_cols(Der_vect);
Im_all_dif_avgnbr_vect = norm_cols(Im_all_dif_avgnbr_vect); 
Im_distToCtr_vect = norm_cols(Im_distToCtr_vect);



%% Put weigths on different phases
% portal and late phase are more important than not injected and arterial
% alpha = sqrt([0.1, 0.2, 0.35, 0.35]);
alpha = [0.1, 0.1, 0.4, 0.4];
% sum(alpha(:))

for i = 1:length(alpha)
	Im_vect_weighted(:,i) = Im_vect(:,i) * alpha(i); 
end

%% Put weigths on derivatives
% select derivatives that are more important in
% order to distinguish the tumor from healthy kidney
% beta = sqrt([0.4, 0.5, 0.1]);
beta = [0.1, 0.1, 0.8];
% sum(beta(:)) 

for i = 1:length(beta)
	Der_vect_weighted(:,i) = Der_vect(:,i) * beta(i);
end 


%% Choose clustering features

% M = [Im_vect(:,3:4), Der_vect(:,3), Im_all_dif_avgnbr_vect, Im_distToCtr_vect];
M = [Im_vect_weighted, Der_vect_weighted, Im_all_dif_avgnbr_vect, Im_distToCtr_vect];

M_all = [Im_coord_cut, M, clusters_on_subsampled_vec];


%% Plot couples of variables used for clustering
% figure
% % plot(Im_vect(:,1), Der_vect(:,1), '*')
% % plot(Im_vect(:,1), Im_coord_cut(:,1), '*')
% plot(Im_vect(:,3), clusters_on_subsampled_vec(:,1), '*')


%% Estimate number of clusters
if ~exist('n_clusters','var')
    eva = evalclusters(M,'kmeans','CalinskiHarabasz','KList',[4:30]);
    n_clusters = eva.OptimalK;
end

%% K-means
sigma = 7;
disp(['Launching k-Means with N=' num2str(n_clusters) ' clusters ...'  ]);
[idx, ~] = kmeans(M_all, n_clusters, 'MaxIter', 500, 'Start', 'uniform');
% [idx, ~] = kmeans3(M_all, n_clusters, 'MaxIter', 500, 'Start', 'uniform');

idx_restore_size(idx_nonv) = idx;
pixel_labels = reshape(idx_restore_size,nrows,ncols, nim); % Image showing the clusters index
figure
for i = sizeBorders+1 : nim-sizeBorders
 	imagesc(pixel_labels(:,:,i)) , colormap(parula) , colorbar %, 'AlphaData', 0.5
	title(['K-means with N.' num2str(n_clusters) ' clusters (tot ' num2str(size(M,2)+4) ' features)']);  
	impixelinfo
	pause
end

% save all clusters labels 
niftiwrite(single(pixel_labels), [folder_out 'all_clusters.nii'], info)

%% Plot Cluster curves
idx_clustokeep = [];
n_clusters = length(unique(pixel_labels)) - 1;
figure
map = parula(n_clusters);
for i = 1:n_clusters
	idxtoplot = find(idx==i); % find pixels with index i (= belonging to cluster i)
 	cl_sizes(i) = length(idxtoplot); % compute size of cluster i
	
%   save average of intensities for every phase, for cluster i
	for p = 1:4
		toplot(p) = mean(Im_vec_cut(idxtoplot, p)); 
    end
	% Compute time derivatives of cluster i
	for j = 1:3
		der(i,j) = (toplot(j+1) - toplot(j)) / (t(j+1) - t(j));
    end
    
	% Plot average curves in time of selected clusters
	if max(toplot) > 80 && toplot(4) < 130 && toplot(4) > 60 && der(i,1) > 0 && der(i,3) < 0 && toplot(4) > toplot(1)% && max(abs(der(i,:))) > 0.78 % && length(idxtoplot)>100 % Criteria on intensities
% 	if toplot(4) < 1150 && der(i,1) > 0 % && max(abs(der(i,:))) > 0.78 % && length(idxtoplot)>100 % Criteria on intensities
% 	if max(abs(der(i,:))) > 0.3 % Criterion on derivatives

        idx_clustokeep = [idx_clustokeep, i];
		plot(t, toplot, '*-', 'MarkerSize',10, 'LineWidth',2, 'DisplayName', ['Cluster ' num2str(i)], 'Color', map(i,:))
		legend('-DynamicLegend');
		legend('show');
		drawnow;
		hold on
        
	end
end
title(['Average curves in time of cluster(s) ' num2str(idx_clustokeep) ' over ' num2str(n_clusters) ' clusters']);
% cl_sizes % show cluster sizes to understand bad clusterings

% Create image with selected clusters
Iresult = ismember(pixel_labels, idx_clustokeep); % logic image with all selected clusters to 1
Iresult = Iresult.*pixel_labels; % Image labelled with indexes of selected clusters

% Visualise selected clusters
figure
for i = sizeBorders + 1 : nim-sizeBorders
	imagesc(Iresult(:,:,i)) , colormap([0,0,0; parula(n_clusters)]) , colorbar
	title(['Cluster(s) ' num2str(idx_clustokeep) ' over ' num2str(n_clusters) ' clusters - Slice N ' num2str(i) ]);
	impixelinfo
	pause
end


%% Extract the largest connected component (cc) (i.e. kidney) from THE selected cluster
I_cc = bwlabeln(Iresult); % label cc
[C,~,ic] = unique(I_cc); % get index of cc
a_counts = accumarray(ic,1); % count elements of every cc
value_counts = [C, a_counts]; % order cc by their size
[~,idx_kidney] = max(a_counts(2:end)); % get index of largest cc

% Visualize selected connected components of selected clusters
I_segmented = double(I_cc == idx_kidney);
figure
for i = sizeBorders +1 : nim-sizeBorders
	imagesc(I_segmented(:,:,i)) , colormap(gray)
	title(['Connected components - Clusters ' num2str(idx_clustokeep) ' over ' num2str(n_clusters) ' clusters - Slice N ' num2str(i) ]);
	pause
end

% Show the segmented four phases 
Im_all_segmented = Im_all .* I_segmented;
figure
for i = 1:nim
	subplot(2,2,1);
	imagesc(Im_all_segmented(:,:,i,1)), colormap(gray), title(['Non injected phase - Slice N  ' num2str(i)])
	subplot(2,2,2);
	imagesc(Im_all_segmented(:,:,i,2)), colormap(gray), title(['Arterial phase - Slice N  ' num2str(i)])
	subplot(2,2,3);
	imagesc(Im_all_segmented(:,:,i,3)), colormap(gray), title(['Portal phase - Slice N   ' num2str(i)])
	subplot(2,2,4);
	imagesc(Im_all_segmented(:,:,i,4)), colormap(gray), title(['Late phase - Slice N ' num2str(i)])
    sgtitle('Segmented images')
	pause
end

% save([folder_out 'Registration/Im_all_segmented'], 'Im_all_segmented'); % Save for further clustering


