function detrended_data = detrend_timecourses(data)
     detrended_data = cell(size(data,1),1);
    for i=1:size(data,1)
        % Linear, first order sinusoidal and mean detrending
        nvols = size(data{i},2);
        linear_trend=1:nvols;
        sin_trend=sin((linear_trend)*2*pi/(nvols-1));
        cos_trend=cos((linear_trend)*2*pi/(nvols-1));
        mean_trend = ones(1,nvols);
        dt_design = [linear_trend;sin_trend;cos_trend;mean_trend];

        trend = data{i}/dt_design;
        est = trend*dt_design;
        detrended_data{i} = data{i} - est;
    end
end

