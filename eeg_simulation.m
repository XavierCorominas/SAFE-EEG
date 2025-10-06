%% Simulate EEG

rng(45);
% Parameters
fs = 5000;                  % Sampling frequency
t  = 0:1/fs:2;              % ~2 seconds
N  = numel(t);              % length
T  = N/fs;                  % total duration (s)

% Frequency components
f_theta = 6; f_mu = 10; f_beta = 20;
a_theta = 1.5; a_mu = 3.0; a_beta = 2;

% Oscillatory components
theta_osc = a_theta * sin(2*pi*f_theta*t);
mu_osc    = a_mu    * sin(2*pi*f_mu*t);
beta_osc  = a_beta  * sin(2*pi*f_beta*t);

% Pink noise (1/f^alpha) via symmetric spectrum (no toolboxes)signal
alpha = 1;
k = 0:floor(N/2);                 % one-sided index
f = (k * fs / N);                 % Hz
f_safe = max(f, 1/T);             % avoid div-by-zero at DC with scalar 1/T
mag = 1 ./ (f_safe).^(alpha/2);   % amplitude ~ 1/f^(alpha/2)
mag(1) = 0;                       % remove DC

% Random phases
phi = 2*pi*rand(1, numel(k));
phi(1) = 0;
if mod(N,2)==0
    phi(end) = 0;                 % Nyquist bin must be real when N is even
end

% One-sided complex spectrum
Xpos = mag .* exp(1j*phi);

% Full spectrum with Hermitian symmetry
X = complex(zeros(1, N));
X(1:numel(k)) = Xpos;
if mod(N,2)==0
    % even N: mirror bins 2..(end-1)
    X(numel(k)+1:end) = conj(Xpos(end-1:-1:2));
else
    % odd N: mirror bins 2..end
    X(numel(k)+1:end) = conj(Xpos(end:-1:2));
end

pink_noise = real(ifft(X));
pink_noise = pink_noise / std(pink_noise);

% Combine all components
signal = theta_osc + mu_osc + beta_osc + 0.5 * pink_noise;
% Optional: Visualize first 2 seconds
figure;
plot(t(1:2000), signal(1:2000), 'k'); hold on;
plot(t(1:2000), theta_osc(1:2000), 'r--');
plot(t(1:2000), mu_osc(1:2000), 'b--');
plot(t(1:2000), beta_osc(1:2000), 'g--');
legend('Combined Signal', 'Theta', 'Mu', 'Beta');
xlabel('Time (s)'); ylabel('Amplitude');
title('Simulated EEG: Theta + Mu + Beta + 1/f Noise');
grid on;


% --- Plot 1: Time domain ---
figure('Position', [100, 100, 1000, 400]);
subplot(2,1,1);
plot(t(1:2000), signal(1:2000), 'k'); hold on;
plot(t(1:2000), mu_osc(1:2000), 'b--');
plot(t(1:2000), 0.5*pink_noise(1:2000), 'r:');
legend('Combined signal', '10 Hz mu', '1/f noise (scaled)');
xlabel('Time (s)');
ylabel('Amplitude');
title('Time-domain signal (first 2 seconds)');
grid on;


% --- Plot 2: Frequency domain (power spectrum) ---
[pxx, f] = pwelch(signal, hamming(10000), 5000, 0.5:0.5:100, fs);
subplot(2,1,2);
plot(f, 10*log10(pxx), 'k', 'LineWidth', 1.5);
xlim([2 50]);
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('Power Spectrum of Combined Signal');
grid on;


%%

%% Filtering and baseline correction for Original non-segmented data

signal_org = signal;
% Example parameters
Fs = 5000;            % Sampling frequency in Hz (example: 10 kHz)
x = signal; % Your input signal

% Bandpass filter limits
f_low = 0.1;   % Lower cutoff frequency in Hz
f_high = 2000; % Upper cutoff frequency in Hz

% Normalize cutoff frequencies (to Nyquist frequency)
Wn = [f_low f_high] / (Fs/2);

% Design a 4th-order Butterworth bandpass filter
[b, a] = butter(2, Wn, 'bandpass');

% Apply the filter
signal = filtfilt(b, a, x);   % zero-phase filtering (no phase distortion)



%% CUT the segments

% --- Parameters you can tweak ---
Fs        = 5000;         % Hz (sampling rate)
win_ms    = 3;            % segment length in milliseconds
start_ms  = 1000;         % first segment starts at 1000 ms (1.000 s)
end_ms    = 2000;         % last segment start no later than 2000 ms (2.000 s)
step_ms   = 100;          % base step between groups (1000,1100,1200,...)
offset_ms = [0, 5];       % segments at base and +5 ms (e.g., 1000 & 1005)

% --- Your data vector (replace with your EEG) ---
% x should be a vector of length N (e.g., N = 20001 as you mentioned)
% Example placeholder (delete next line and load your real data):
% x = randn(20001,1);  % <-- REMOVE this line and use your actual x

x = signal(:);               % ensure column vector
N = length(x);

% Build time vector (assuming t(1)=0)
t = (0:N-1).' / Fs;     % seconds

% Compute window length in samples
L = max(1, round(win_ms * Fs / 1000));  % 3 ms -> 3 samples at 1 kHz

% Compute all segment start times (ms): 1000,1005,1100,1105,...,2000
bases_ms  = start_ms:step_ms:end_ms;
starts_ms = sort([bases_ms, bases_ms + offset_ms(2)]); % add +5ms
% Keep only those <= end_ms
starts_ms = starts_ms(starts_ms <= end_ms);

% Convert each start time (ms) to a 1-based sample index
% idx = round(ms/1000*Fs) + 1 maps 0 ms -> 1
start_idx = round(starts_ms/1000 * Fs) + 1;

% Create a keep-mask and mark samples to remove
% Replace those windows with NaN
x_nan = x;
for k = 1:numel(start_idx)
    s = start_idx(k);
    e = min(s + L - 1, N);
    x_nan(s:e) = NaN;
end

% Interpolate the segments

% --- Fill NaN gaps via cubic interpolation ---
method = 'pchip';  % options: 'pchip' (shape-preserving cubic) or 'spline' (cubic spline)

valid = ~isnan(x_nan);
idx_nan = isnan(x_nan);

x_interp = x_nan;  % start from NaN-masked signal
if sum(valid) >= 2
    % Interpolate only at the NaN positions using surrounding valid samples
    x_interp(idx_nan) = interp1(t(valid), x_nan(valid), t(idx_nan), method);
else
    warning('Not enough valid points to interpolate.');
end


% --- Plot results ---
figure('Name','EEG Processing','Color','w');
subplot(3,1,1);
plot(t, x, 'k');
xlabel('Time (s)'); ylabel('Amplitude');
title('Original EEG Signal');
xlim([0 max(t)]);

subplot(3,1,2);
plot(t, x_nan, 'b');
xlabel('Time (s)'); ylabel('Amplitude');
title('EEG with 3 ms segments replaced by NaN (visual gaps)');
xlim([0 max(t)]);

subplot(3,1,3);
plot(t, x_interp, 'r');
xlabel('Time (s)'); ylabel('Amplitude');
title(sprintf('EEG after cubic interpolation (%s) over NaN gaps', method));
xlim([0 max(t)]);

signal_demeaned= x_interp;

%% Bandpas filter
%{
% Example parameters
Fs = 5000;            % Sampling frequency in Hz (example: 10 kHz)
x = x_interp; % Your input signal

% Bandpass filter limits
f_low = 0.1;   % Lower cutoff frequency in Hz
f_high = 2000; % Upper cutoff frequency in Hz

% Normalize cutoff frequencies (to Nyquist frequency)
Wn = [f_low f_high] / (Fs/2);

% Design a 4th-order Butterworth bandpass filter
[b, a] = butter(2, Wn, 'bandpass');

% Apply the filter
y = filtfilt(b, a, x);   % zero-phase filtering (no phase distortion)

% demean
signal_demeaned = detrend(y, 'linear');  % removes mean only




%% Filtering and baseline correction for Original non-segmented data


% Example parameters
Fs = 5000;            % Sampling frequency in Hz (example: 10 kHz)
x = signal; % Your input signal

% Bandpass filter limits
f_low = 0.1;   % Lower cutoff frequency in Hz
f_high = 2000; % Upper cutoff frequency in Hz

% Normalize cutoff frequencies (to Nyquist frequency)
Wn = [f_low f_high] / (Fs/2);

% Design a 4th-order Butterworth bandpass filter
[b, a] = butter(2, Wn, 'bandpass');

% Apply the filter
signal_2 = filtfilt(b, a, x);   % zero-phase filtering (no phase distortion)

% Demean
signal = detrend(signal_2, 'linear');  % removes mean only
%}


%% Plot power spectra of both signals.

figure;
plot(t(1:10001), signal(1:10001), 'k'); hold on;
plot(t(1:10001), theta_osc(1:10001), 'r--','LineWidth',1.5);
plot(t(1:10001), mu_osc(1:10001), 'b--', 'LineWidth',1.5);
plot(t(1:10001), beta_osc(1:10001), 'g--','LineWidth',1.5);
legend('Combined Signal', 'Theta', 'Mu', 'Beta');
xlabel('Time (s)'); ylabel('Amplitude');
title('Simulated EEG: Theta + Mu + Beta + 1/f Noise');
grid on;
set(gca, 'FontSize', 20, 'FontName', 'Arial');


% --- Plot 2: Frequency domain (power spectrum) ---
figure;
[pxx2, f2] = pwelch(signal, hamming(10000), 5000, 0.1:0.5:100, fs);
xlim([0 45]);
plot(f2, 10*log10(pxx2), 'k', 'LineWidth', 3);
xlim([0 45]);
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('Before Processing');
grid on
set(gca, 'FontSize', 20, 'FontName', 'Arial');


figure
[pxx, f] = pwelch(signal_demeaned, hamming(10000), 5000, 0.1:0.5:100, fs);
plot(f, 10*log10(pxx), 'r', 'LineWidth', 3);
xlim([0 45]);
xlim([0 45]);
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('After Processing');
grid on
set(gca, 'FontSize', 20, 'FontName', 'Arial');



figure
plot(f, 10*log10(pxx), 'r', 'LineWidth', 3);
hold on
plot(f2, 10*log10(pxx2), 'k', 'LineWidth', 3);
xlim([1 45]);
xlim([1 45]);
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('Processing effects');
grid on
set(gca, 'FontSize', 20, 'FontName', 'Arial');
legend({'Original Signal', 'Processed Signal'}, ...
       'Location', 'northeast', 'FontSize', 20);


%%
