clear;

%load data
ecg_signal = load_ecg();
idx_start = 5.3e4;
idx_end = 5.5e4-1;
excerpt = ecg_signal(idx_start:idx_end);

%slow version
wdw_length = 200;
plot_signal_and_anomaly(excerpt, wdw_length)
keyboard

%faster version after quantization
wdw_length_quantization = 10;
excerpt_quantized = apply_func_windows(excerpt, @(vals)(mean(vals)), wdw_length_quantization);
new_wdw_length = wdw_length / wdw_length_quantization;
plot_signal_and_anomaly(excerpt_quantized, new_wdw_length)


function [d_min, idx] = closest_subsequence(seq, subseq, func_dist, wdw_length, idcs_skip)
    disp(num2str(idcs_skip(end)))
    d = apply_func_windows(...
                    seq,...
                    @(subseq_)(func_dist(subseq_, subseq)),...
                    wdw_length,...
                    'dist_windows', 1);
    d(idcs_skip) = Inf;
    [d_min, idx] = min(d);
end

function [d_max, idx] = strangest_subsequence(seq, func_dist, wdw_length)
    f_idcs_skip = @(idcs)(max(idcs(1)-length(idcs)+1, 1) : idcs(end));
    d = apply_func_windows(...
                    [[1:length(seq)]', seq],...
                    @(subseq_)(closest_subsequence(seq, subseq_(:, 2), func_dist, wdw_length, f_idcs_skip(subseq_(:, 1)))),...
                    wdw_length,...
                    'dist_windows', 1);
    [d_max, idx] = max(d);
end

function plot_signal_and_anomaly(signal, wdw_length)
    figure; hold on; grid on; plot(signal, 'LineWidth',2);
    func_dist = @(x,y)(sum(abs(x-y)));
    [~, idx] = strangest_subsequence(signal, func_dist, wdw_length);
    idcs_anomaly = idx:idx+wdw_length-1;
    plot(idcs_anomaly, signal(idcs_anomaly), 'color', 'r', 'LineWidth', 2)
end


function signal_out = apply_func_windows(signal_in, func, samples_per_wdw, varargin)
    % define default behaviour
    dist_windows = samples_per_wdw;    %wdw_size can be diff from wdw_dist
    b_full_wdw_only = true;            %consider partial windows?
    
    nargin = length(varargin);
    idcs_vararg_names = 1:2:nargin;
    for idx_vararg_name = idcs_vararg_names
        name =  varargin{idx_vararg_name};
        switch name
            case 'dist_windows'
                dist_windows = varargin{idx_vararg_name+1};
            case 'b_full_wdw_only'
                b_full_wdw_only = varargin{idx_vararg_name+1};
            otherwise
                error(['Unknown input parameter: ' name])
        end
    end

    %input: column convention
    if (isvector(signal_in))
        signal_in = signal_in(:);
    end
    [num_samples, num_feats] = size(signal_in);
    assert(num_samples>num_feats); %defensive prog
    
    num_full_and_partial_windows = ceil(num_samples/dist_windows);
    num_full_windows = ceil((num_samples-samples_per_wdw+1)/dist_windows);
    if (b_full_wdw_only)
        num_windows = num_full_windows;
    else
        num_windows = num_full_and_partial_windows;
    end

    signal_out = [];
    for idx_wdw=1:num_windows
        %idcs are 1-based
        idx_start = (idx_wdw-1) * dist_windows +1;
        idx_end = idx_start + samples_per_wdw - 1;
        %last wdw could be shorter
        if (idx_end > num_samples)
            assert(~b_full_wdw_only);
            idx_end = num_samples;
        end
        windowed_signal_in = signal_in(idx_start:idx_end, :);
        windowed_result = func(windowed_signal_in);
        signal_out = [signal_out; windowed_result];
    end
end

function signal = load_ecg()
    file_name = 'sel102.dat';
    fid=fopen(file_name,'r');
    time=1000;
    f=fread(fid,2*360*time,'ubit12');
    signal=f(1:2:length(f));
end

function [] = draw_example_neighbors()
ys = [800, 1200];
x1 = 600; x2 = 800;
line([x1, x1],ys, 'LineWidth', 2, 'color', 'k')
line([x2, x2],ys, 'LineWidth', 2, 'color', 'k')

x1 = 200; x2 = 400;
line([x1, x1],ys, 'color','g')
line([x2, x2],ys,'color', 'g')

x1 =650; x2 = 850;
line([x1, x1],ys, 'color','r')
line([x2, x2],ys,'color', 'r')
end