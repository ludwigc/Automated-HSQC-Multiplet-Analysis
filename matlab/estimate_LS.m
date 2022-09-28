function [p_LS,y_stim,index_good,index_good_component] = estimate_LS(X,y)
%%
% The contribution of each multiplet component of X is computed, fitting X
% to y by nonlinear least square method
%%
% Non-lin apprach
try
    options = optimset('MaxIter',1e8, 'TolX', 1e-2*max(max(X)));
catch
    options = optimset('MaxIter',1e8, 'TolX', 1e8);
end
[p_LS, resnorm,residual,exitflag,output,lambda] = lsqnonneg(X,y, options);
if exitflag == 0
    keyboard
    options.TolX = max(max(X));
    [p_LS, resnorm,residual,exitflag,output,lambda] = lsqnonneg(X,y, options);
    if(exitflag == 0)
        keyboard
    end
end
% estimated signal
y_stim = X*p_LS;
% compute residual
res = y - y_stim;
% R2 value (goodness of fit)
index_good = 1 - sum((res).^2)/sum((y - mean(y)).^2);
p_LS = p_LS/norm(p_LS);       % estimated parameters
% Compute the contribution of each component to the total R2 score
index_good_component = [];
if sum(isnan(p_LS))> 0
    p_LS = zeros(size(p_LS));
else
    index_good_component = zeros(length(p_LS),1);
    for c = 1:length(p_LS)
        p_LS_tmp = zeros(size(p_LS));
        p_LS_tmp(c) = p_LS(c);
        y_single = X*p_LS_tmp;            
        res = y - y_single;
        index_good_component(c) = sum((res).^2)/sum((y - mean(y)).^2); % R2 value 
    end
end
end