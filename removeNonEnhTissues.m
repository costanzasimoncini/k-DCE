function [ Im_all ] = removeNonEnhTissues( Im_all, Im_Der, v )

% Select and remove voxels belonging to the biggest connected component
% with low contrast enhancement derivative

I_binThrEnh = zeros(size(Im_all(:,:,:,1)));

[nrows, ncols, nim] = size(Im_all(:,:,:,1));
% Remove pixels with too small variation
for i = 1 : nrows
	for j = 1 : ncols
		for k = 1 : nim
			if max(abs(Im_Der(i,j,k,:))) < 0.15 % || max(abs(Im_Der(i,j,k,:))) > 7 % (to remove organs not registered)
				I_binThrEnh(i,j,k) = 1;
		end
		end
	end
end

I_cc = bwlabeln(I_binThrEnh);
[C,~,ic] = unique(I_cc);
a_counts = accumarray(ic,1);
value_counts = [C, a_counts]; 

% Select only biggest connected components
[~,idx_cc_toremove] = max(a_counts(2:end));
% idx_cc_toremove = value_counts(value_counts(2:end,2)>500000); % 2:end to ignore the first row (background)
% idx_cc_toremove = idx_cc_toremove + 1;  % to shift indexes due to line above


% For visualisation

% I_toremove = double(I_cc ~= idx_cc_toremove);
I_toremove = double(ismember(I_cc, idx_cc_toremove));
% I_toremoveArt = I_toremove .* Im_all(:,:,:,2);
se = strel('sphere',1);
I_toremove_dil = imdilate(I_toremove,se);


% figure 
% for i = round(nim/2) : round(nim/2) + 1
% 	subplot(2,2,1)
% 	imagesc(I_binThrEnh(:,:,i)), colormap(gray), title(['Binarized im of non enhancing tissue - Slice ' num2str(i)])
% 	subplot(2,2,2)
% 	imagesc(I_cc(:,:,i)), colormap(gray), title(['Connected components - Slice ' num2str(i)])
% 	subplot(2,2,3)
% 	imagesc(I_toremove(:,:,i)), colormap(gray), title(['Largest cc to be removed - Slice ' num2str(i)])
% 	subplot(2,2,4)
% 	imagesc(I_toremove_dil(:,:,i)), colormap(gray), title(['Largest cc on Arterial phase  - Slice ' num2str(i)])
% 	impixelinfo
%     sgtitle('Selection of connected components corresponding to non enhanced tissues')
% 	pause
% end


% Remove pixels with too small variation
% Im_all((I_cc==idx_cc_toremove), :) = v;	
for i = 1 : nrows
	for j = 1 : ncols
		for k = 1 : nim
% 			if I_cc(i,j,k) == idx_cc_toremove
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
	imagesc(Im_all(:,:,i,1)), colormap(gray), title(['Non injected phase - Slice ' num2str(i)])
	subplot(2,2,2);
	imagesc(Im_all(:,:,i,2)), colormap(gray), title(['Arterial phase - Slice ' num2str(i)])
	subplot(2,2,3);
	imagesc(Im_all(:,:,i,3)), colormap(gray), title(['Portal phase - Slice ' num2str(i)])
	subplot(2,2,4);
	imagesc(Im_all(:,:,i,4)), colormap(gray), title(['Late phase - Slice ' num2str(i)])
    sgtitle('Images without connected components corresponding to non enhancing tissues')
	pause
end

end

