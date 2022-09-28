function [index, vett_corr_max] = Find_match(template,data,range_C,range_H )
%% Computes crosscorrelation between argin(1) and the columns of data
% Inputs: template (i.i. ICA isgnal);data (2D NMD spectrum); range_C 
%(Carbon shifts range to restrict the spectrum area); range_H (
% Proton shifts range to restrict the spectrum area);
% Output: index (index of range_H where the correlation is the
% maximum);vett_corr_max (contains the corr valus between template and each
% column of data in range_H)
%%
range_indx = range_H(1):range_H(end);
vett_corr = zeros(length(range_indx),1);
for i = 1:length(range_indx)
    vett_corr(i,1) = corr(data(range_C(1):range_C(end),range_indx(i))/norm(data(range_C(1):range_C(end),range_indx(i))),template(range_C(1):range_C(end))); 
end
vett_corr_max = max(vett_corr);
if length(vett_corr_max) > 1
    keyboard;
end
index = range_indx(vett_corr == vett_corr_max(1));
if length(index) > 1
    keyboard;
end
end

