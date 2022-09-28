function [s_final, D_best, scoreMax] = improveAlignment(s,y_orig,D_range)
    s_final = zeros(size(s));
    D_best = zeros(length(s(1,:)),1);
    for iter_mult = 1:length(s(1,:))
        score = zeros(length(D_range),2);
        for j = 1:length(D_range)
            X_temp = s;
            X_temp(:, iter_mult) = circshift(X_temp(:,iter_mult),-D_range(j));
            [~,~,R2] = estimate_LS(X_temp,y_orig);
            score(j,1) = R2;
            score(j,2) = D_range(j);
        end
        score_max = find(score(:,1) == max(score(:,1)));
        D_best(iter_mult) = score(score_max(1),2);
        scoreMax{iter_mult} = score;
        X_temp = s;
        s_final(:, iter_mult) = circshift(X_temp(:, iter_mult), -D_best(iter_mult));
    end
return