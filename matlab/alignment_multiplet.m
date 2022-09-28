function[X_temp,corr_al,sim_no,D,D_best] = alignment_multiplet(y,y_orig,s,range_C,set_comb,isfinal, minCorr, maxRange)
%%
% Performs to methods to align at best s to y_orig. The first method
% consider only the multiplets with good alignment (> minCorr). The second
% methods considers all the combinations of the multiplets and returns the
% combinations that lead to best fitting (LS method is used).
%%
s_temp = s;
y_temp_al = y;
corr_al = zeros(length(s(1,:)),1);
%% METHOD 1: Find the correlation score between each single component and the ICA's components and sum up the components with score > minCorr
for j = 1:length(s(1,:)) % for each multiplet component 
    [~,~,D] = alignsignals(y_temp_al,s_temp(:,j));
    s_temp(:,j) = circshift(s_temp(:,j),-D);
    len_temp = length(y_temp_al);
    len_s_temp = length(s_temp);
    if len_s_temp > len_temp
        s_temp_new(:,1) = s_temp(1:len_temp);
    else
        s_temp_new = zeros(len_temp,1);
        s_temp_new(1:len_s_temp,1) = s_temp(:,1);
    end
    corr_al(j) = corr(y_temp_al/norm(y_temp_al),s_temp_new/norm(s_temp_new));
end
% Sum up the simulated component which have corr score > minCorr
D_range = -maxRange:maxRange;
sim_sum = zeros(length(s(:,1)),1);
sim_no = zeros(length(corr_al),1);
s1 = zeros(size(s));
count_no = 1;
for j = 1:length(s(1,:))
    if corr_al(j) >= minCorr
       sim_sum(:,1) = sim_sum(:,1) + s(:,j);
       s1(:,j) = s(:,j);
       comb_final_corr(count_no) = j;
       count_no = count_no + 1;
    end    
end
s1 = circshift(s1,-D);
[X_temp1, D_best1] = improveAlignment(s1,y_orig,D_range);
[~,~,R2_1] = estimate_LS(X_temp1,y_orig);
%% METHOD 2: Find simulated components by alignment_comb
[comb_final_comb,~,~,D] = alignment_comb(y,s);
s2 = zeros(size(s));
for k = 1:length(comb_final_comb(:,1))
    sim_sum = sim_sum + s(:,comb_final_comb(k));
    s2(:,comb_final_comb(k)) = s(:,comb_final_comb(k));
end
s2 = circshift(s2,-D);
[X_temp2, D_best2] = improveAlignment(s2,y_orig,D_range);
[~,~,R2_2] = estimate_LS(X_temp2,y_orig);
if R2_1 > R2_2
    %disp("Method one gives better fitting!")
    X_temp = X_temp1;
    D_best = D_best1;
else
    X_temp = X_temp2;
    D_best = D_best2;
end
