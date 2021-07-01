function [x_mean,y_mean,CI_up,CI_down] = plot_ci(x,y,color,varargin)
    %%plot the 95% confidence interval with x as vector, y the matrix of
    %%the data, with the same number of rows as x, and color
    % varargin{1} = omitnan; if 1 omitnan if 0 don't
    L = size(y,2);
    ts = tinv([0.025,.975],size(y,1)-1);
    if ~isempty(varargin)&&varargin{1}==1
        sigm = std(y,[],2,'omitnan');
    else
        sigm = std(y,[],2);
    end
    Conf_intervals(:,1) = ts(1)*sigm/sqrt(L);	% Confidence Intervals CI = mean(x)+- t * (s / square(n))
    Conf_intervals(:,2) = ts(2)*sigm/sqrt(L);
    
    if ~isempty(varargin)&&varargin{1}==1
        CI_up = mean(y,2,'omitnan') + Conf_intervals(:,1);
        CI_down = mean(y,2,'omitnan') - Conf_intervals(:,1);
    else
        CI_up = mean(y,2) + Conf_intervals(:,1);
        CI_down = mean(y,2) - Conf_intervals(:,1);        
    end
    
    CI_up (isnan(CI_up)) = 0;
    CI_down (isnan(CI_down)) = 0;
    
    XX= [x;flipud(x)];
    YY = [CI_up;flipud(CI_down)];
    hold on
    fill(XX,YY,1,...
        'facecolor',color, ...
        'edgecolor','none', ...
        'facealpha', 0.3);
    if ~isempty(varargin)&&varargin{1}==1  
        plot(x,mean(y,2,'omitnan'),'color',color)
        y_mean = mean(y,2,'omitnan');
    else
        plot(x,mean(y,2),'color',color)
        y_mean = mean(y,2);        
    end
    x_mean = x;
end