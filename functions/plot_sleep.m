function h = plot_sleep(x, y, sleepscore_time, wake_binary, sws_binary, REM_binary, MA_binary)
% PLOT_SLEEP plots a graph defined by X and Y with color coded zones
% marking wake, sws, REM, (and MA)
% Input arguments:
%   x: time vector for the graph to be plotted (chould match length of y)
%   y: values to be plotted (e.g. FibPho, EEG, or EMG trace)
%   sleepscore_time: time vector matching length of binary sleepscore vectors
%   wake_binary: binary vector where 1 indicates wakufulness (should be same length as sleepscore_time)
%   sws_binary: binary vector where 1 indicates slow-wave sleep (should be same length as sleepscore_time)
%   REM_binary: binary vector where 1 indicates REM sleep (should be same length as sleepscore_time)
%   MA_binary: binary vector where 1 indicates microarousals (shouold be
%   same length as sleepscore_time) - this input argument can be ommitted
% Output arguments:
%   h: the output is a plot, and no argument is needed

    ymin = min(y); ymax = max(y);
    time_matrix = [sleepscore_time fliplr(sleepscore_time)];
    y_height = ymax-ymin;
    h = patch(time_matrix, [ymin+zeros([1, length(wake_binary)]) ymin+y_height*fliplr(wake_binary)], 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    hold on
    patch(time_matrix, [ymin+zeros([1, length(sws_binary)]) ymin+y_height*fliplr(sws_binary)], 'r', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    patch(time_matrix, [ymin+zeros([1, length(REM_binary)]) ymin+y_height*fliplr(REM_binary)], 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    if nargin == 7
    patch(time_matrix, [ymin+zeros([1, length(MA_binary)]) ymin+y_height*fliplr(MA_binary)], 'y', 'EdgeColor', 'none');
    end
    plot(x, y, 'k')
    
end