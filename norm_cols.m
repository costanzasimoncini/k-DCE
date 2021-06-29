function [ Mnorm ] = norm_cols(M)

for i = 1 : size(M,2)
	Mnorm(:,i) = M(:,i) - mean(M(:,i)); 
	Mnorm(:,i) = Mnorm(:,i)/std(Mnorm(:,i));
end

end

