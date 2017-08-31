function t = tinv(alpha, v)
% % Variables: 
% % t: t-statistic
% % v: degrees of freedom
    if v <= 0
        t = nan;
        return;
    end
    persistent tinvCache975;
    persistent tinvCache950;
    if alpha == 0.5
        t = 0.5;
        return;    
    elseif alpha == 0.975 || 1-alpha == 0.975
        if numel(tinvCache975) >= v && tinvCache975(v) > 0
            t = tinvCache975(v);
            if alpha < 0.5
                t = -t;
            end
            return;
        end
    elseif alpha == 0.950 || 1-alpha == 0.950
        if numel(tinvCache950) >= v && tinvCache950(v) > 0
            t = tinvCache950(v);
            if alpha < 0.5
                t = -t;
            end
            return;
        end
    end
    tdist2T = @(t,v) (1-betainc(v/(v+t^2),v/2,0.5));                                % 2-tailed t-distribution
    tdist1T = @(t,v) 1-(1-tdist2T(t,v))/2;                                          % 1-tailed t-distribution
    t = fzero(@(tval) (max(alpha,(1-alpha)) - tdist1T(tval,v)), 5);  % T-Statistic Given Probability ‘alpha’ & Degrees-Of-Freedom ‘v’
    if alpha == 0.975 || 1-alpha == 0.975
        if alpha < 0.5
            t = -t;
        end
        tinvCache975(v) = t;
    elseif alpha == 0.950 || 1-alpha == 0.950
        if alpha < 0.5
            t = -t;
        end
        tinvCache950(v) = t;
    end
end