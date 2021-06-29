function [ Im_all ] = removeBones( Im_all, t_high, v )

% Select and remove voxels belonging to the biggest connected component with values higher than t_high

[nrows, ncols, nim] = size(Im_all(:,:,:,1));
%  define I_max - Image with the maximum over the 4 times
for i = 1 : nrows
	for j = 1 : ncols
		for k = 1 : nim
			I_max(i,j,k) = max(Im_all(i,j,k,:));
		end
	end
end

% fig = figure;
% idx_central_im = round((size(I_max,3)/2));
% for i = idx_central_im : idx_central_im + 2
% 	imagesc(I_max(:,:,i)), colormap(gray), colorbar
% 	truesize(fig, [3*nrows 3*ncols])
% 	impixelinfo
%     title('Max value over the 4 times')
% 	pause
% end



I_binThrGr = double (I_max >= t_high);
% I_binThrGr = double (Im_all(:,:,:,4) >= 210); % On late phase to get bones and urine
% Compute connected components
I_cc = bwlabeln(I_binThrGr);

[C,~,ic] = unique(I_cc);
a_counts = accumarray(ic,1);
value_counts = [C, a_counts];

% Select only connected components with more than 1000 voxels
idx_cc_toremove = value_counts(value_counts(2:end,2)>1000); % 2:end to ignore the first row (background)
idx_cc_toremove = idx_cc_toremove + 1;  % to shift indexes due to line above
% [M,idx_cc_toremove] = max(a_counts(2:end));

% for visualisation
I_toremove = double(ismember(I_cc, idx_cc_toremove));
se = strel('sphere',1);
I_toremove_dil = imdilate(I_toremove,se);
% I_toremoveArt = I_toremove .* Im_all(:,:,:,2);


% figure 
% for i = round(nim/2) : round(nim/2) + 1
% 	subplot(2,2,1)
% 	imagesc(I_binThrGr(:,:,i)), colormap(gray), title(['Binary im thr=' num2str(t_high) ' - Slice ' num2str(i)])
% 	subplot(2,2,2)
% 	imagesc(I_cc(:,:,i)), colormap(gray), title(['Connected components - Slice  ' num2str(i)])
% 	subplot(2,2,3)
% 	imagesc(I_toremove(:,:,i).*Im_all(:,:,i,2)), colormap(gray), title(['Largest cc to be removed - Slice  ' num2str(i)])
% 	subplot(2,2,4)
% 	imagesc(I_toremove_dil(:,:,i).*Im_all(:,:,i,2)), colormap(gray), title(['Largest cc on Arterial phase - Slice ' num2str(i)])
% 	impixelinfo
%     sgtitle('Selection of connected components corresponding to bones')
% 	pause
% end


% Im_all(I_toremove_dil==1,:) = v;
for i = 1 : nrows
	for j = 1 : ncols
		for k = 1 : nim
			if I_toremove_dil(i,j,k) == 1
				Im_all(i,j,k,:) = v;
			end
		end
	end
end

% Show all phases
figure
for i = round(nim/2) : round(nim/2) + 1
	subplot(2,2,1);
	imagesc(Im_all(:,:,i,1)), colormap(gray), title(['Non injected - Slice ' num2str(i)])
	subplot(2,2,2);
	imagesc(Im_all(:,:,i,2)), colormap(gray), title(['Arterial phase - Slice ' num2str(i)])
	subplot(2,2,3);
	imagesc(Im_all(:,:,i,3)), colormap(gray), title(['Portal phase - Slice ' num2str(i)])
	subplot(2,2,4);
	imagesc(Im_all(:,:,i,4)), colormap(gray), title(['Late phase - Slice ' num2str(i)])
    sgtitle('Images without connected components corresponding to bones')
	pause
end

end