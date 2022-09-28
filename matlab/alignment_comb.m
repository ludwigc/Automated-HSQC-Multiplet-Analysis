function[comb_final,subset,DD,D] = alignment_comb(y,s)
    s(length(s(1,:)),1);
    D  = 0;
    for j = 1:length(s(1,:))
        combos = nchoosek(1:length(s(1,:)),j);
        corr_comb_sub = zeros(length(combos(:,1)),1);
        for k = 1:length(combos(:,1))
            sim_sum = zeros(size(s));
            for l = 1:length(combos(k,:))
                sim_sum(:,combos(k,l)) = s(:,combos(k,l));
            end
            if(j==1)
                [~,~,D] = alignsignals(y/norm(y),sum(sim_sum,2)/norm(sim_sum));  %
                DD(k)   = D;
            end
            s_temp = sim_sum;
            s_temp = circshift(s_temp,-D);
            len_temp = length(y);
            len_s_temp = length(s_temp);
            if len_s_temp > len_temp
                s_temp_new(:,1) = s_temp(1:end-1);
            end
            if len_s_temp < len_temp
                s_temp_new = zeros(len_temp,1);
                s_temp_new(1:len_s_temp,1) = s_temp(:,1);
            end
            if len_s_temp ~= len_temp
                [parEst,~,corr_comb_sub(k,1)] = estimate_LS(s_temp_new*max(y)/max(max(s_temp_new)),y/max(y));
                paramMat(k,:) = parEst;
            else
                [parEst,y_stim,corr_comb_sub(k,1)] = estimate_LS(s_temp*max(abs(y))/max(max(s_temp)),y/max(abs(y)));
                paramMat(k,:) = parEst;
            end
        end 
        index_max = find(max(corr_comb_sub) == corr_comb_sub);
        try
            corr_max_sub(j,1) = corr_comb_sub(index_max(1));
            corr_max_sub(j,2) = index_max(1);                        
            corr_max_sub(j,3) = j;                      
        catch
            keyboard
        end
        try
            corr_max_sub(j,3+(1:length(paramMat(index_max(1),:)))) = paramMat(index_max(1),:)';                      
        catch
            keyboard
        end
        if(j==1)
            index_first = find(max(corr_comb_sub) == corr_comb_sub);
            D = DD(index_first(1));
        end
    end

    corr_comb = find(max(corr_max_sub(:,1)) == corr_max_sub(:,1)); 
    corr_max = corr_comb(1);
    subset = corr_max_sub(corr_max,3); 
    comb_tab = nchoosek(1:length(s(1,:)),subset);
    comb_final(:,1) = (comb_tab(corr_max_sub(corr_max,2),:))';
end