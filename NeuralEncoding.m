load('spiketrain1.mat')
load('spiketrain2.mat')

%% Part I: Spike Train Statistics
%%
% *A) Binary Representations*

% GRAPH OF SPIKETRAIN RAW DATA
figure(1);
plot(timeaxis, Vm, 'LineWidth', 1.5);
xlabel('sec');
ylabel('mV');
title('Spiketrain 1 Raw Data');

timepts = length(Vm);
endtime = timeaxis(timepts);

% Q1: Sampling frequency
rate = timeaxis(1,2) - timeaxis(1,1);
samplingFrequency = ['Vm was sampled at a frequency of ', num2str(1/(1000*rate)), ' kHz'];
disp(samplingFrequency);

% Q2: Threshold
thresh = -40; % all action potentials surpass -40 mV
indices = {};

% Find all indices when threshold is surpassed
for j = 2:timepts
    results = find(Vm(j)>thresh & Vm(j-1)<thresh);
    if isempty(results) == 0
        indices = [indices, j];
    end
end

indices = cell2mat(indices);
spiketimes = indices * rate; % convert indices into spike times

% Q3: Binary representation

% Since the spiketrain is reported in seconds, the length of the binary
% vector is calculated by multiplying the endtime by the number of bins/sec
spiketrain_1ms = zeros(1, round(endtime * 1000));
spiketrain_05ms = zeros(1, round(endtime * 2000));
spiketrain_01ms = zeros(1, round(endtime * 10000));

% Spiketimes are converted to and placed in their corresponding bin number 
% by multiplying each time by bins/sec
for i = 1:length(indices)
    spiketrain_1ms(round(spiketimes(i) * 1000)) = 1000;
    spiketrain_05ms(round(spiketimes(i) * 2000)) = 2000;
    spiketrain_01ms(round(spiketimes(i) * 10000)) = 10000;
end

figure(2);
subplot(3,1,1);
plot(spiketrain_1ms, 'LineWidth', 1);
xlim([0, 10000]); ylim([-500, 1500]);
xlabel('Bin (1 msec)'); ylabel('Bins/sec');
title('Binary representations of spike trains at first 10000 bins');

subplot(3,1,2);
plot(spiketrain_05ms, 'LineWidth', 1);
xlim([0, 10000]); ylim([-1000, 3000]);
xlabel('Bin (0.5 msec)'); ylabel('Bins/sec');

subplot(3,1,3);
plot(spiketrain_01ms, 'LineWidth', 1);
xlim([0, 10000]); ylim([-5000, 15000]);
xlabel('Bin (0.1 msec)'); ylabel('Bins/sec');

%%
% *B) Interspike Interval Statistics*

% Q1: Interspike interval sequence
I = diff(spiketimes);

% Q2: Interspike interval histogram
figure(3);
I_hist = round(diff(spiketimes) * 1000); % convert to ms
I_hist = I_hist(I_hist <= 200); % remove all values above 200 ms
histogram(I_hist, 200);
title('Interspike interval histogram');
xlabel('Interval (msec)');
ylabel('Count');

% Q3: Variability
CV = std(I)/mean(I)
disp('This neuron is more variable than a Poisson process, which has CV = 1.')

% Q4: Correlation coefficients
lag = 200;
[P, lags] = xcorr(I-mean(I), lag, 'normalized');

figure(4);
plot(lags(1,201:251), P(1,201:251), '-x', 'LineWidth', 1);
xlabel('lag'); ylabel('Correlation coefficient (rho)');
title('Interspike interval correlation coefficients');

disp('None of the non-zero correlation coefficients are 0, so this spike train is NOT a renewal process.')
%%
% *C) Autocorrelation Function and Power Spectrum*
lag = 200;

% Q1: Autocorrelation functions
A_1ms = xcorr(spiketrain_1ms, lag, 'normalized');
A_05ms = xcorr(spiketrain_05ms, lag, 'normalized');
A_01ms = xcorr(spiketrain_01ms, lag, 'normalized');

figure(5);
subplot(3,1,1);
plot(lags, A_1ms, 'LineWidth', 1);
xlim([-30, 30]);
xlabel('lag');
title('Autocorrelation function for 1 msec binary representation');
subplot(3,1,2);
plot(lags, A_05ms, 'LineWidth', 1);
xlim([-30, 30]);
xlabel('lag'); ylabel('Correlation coefficient (rho)');
title('0.5 msec');
subplot(3,1,3);
plot(lags, A_01ms, 'LineWidth', 1);
xlim([-30, 30]);
xlabel('lag');
title('0.1 msec');

% Q2: Power spectra
[pxx1,~] = pwelch(spiketrain_1ms, bartlett(2048), 1024, 2048, 1000);
[pxx05,~] = pwelch(spiketrain_05ms, bartlett(2048), 1024, 2048, 2000);
[pxx01,f] = pwelch(spiketrain_01ms, bartlett(2048), 1024, 2048, 10000);

figure(6);
subplot(3,1,1);
plot(f, pxx1)
xlim([-20, 5000]); 
ylabel('Power (spk^2/s)'); xlabel('Freq (Hz)');
title('Power spectrum of 1 ms binary representation');
subplot(3,1,2);
plot(f, pxx05)
xlim([-20, 5000]);
ylabel('Power (spk^2/s)'); xlabel('Freq (Hz)');
title('0.5 ms');
subplot(3,1,3);
plot(f, pxx01);
xlim([-20, 5000]); 
ylabel('Power (spk^2/s)'); xlabel('Freq (Hz)');
title('0.1 ms');


%% Part II: Measures of Neural Encoding
%%
% *A) Raster Plots and PSTH*

% Q1: Raster plotting the data

figure(7);
subplot(2,1,1);
plot(data(:,2), data(:,1), '*');
xlabel('time (msec)'); ylabel('Trial number');
title('Spike times by trial');

subplot(2,1,2);
plot(data(:,2), data(:,1), '*');
xlabel('time (msec)'); ylabel('Trial number');
title('Sample of plot at first 50 msec');
xlim([0, 50]);
ylim([0, 100]);

% Q2: PSTH
endtime2 = data(end, 2);
PSTH = zeros(1, round(endtime2));

% bin the spike times and tally them in the PSTH
for i = 1:length(data)
   tt = round(data(i,2));
   PSTH(1,tt) = PSTH(1,tt) + 1;
end

binwidth = 0.001;
epochs = 100;
PSTH = PSTH / epochs / binwidth;

figure(8);
plot(PSTH);
xlabel('time (sec)');
ylabel('Firing rate (Hz)');
title('Sample of PSTH at first 100 sec');
xlim([0, 100]);
ylim([0, 250]);

%%
% *B) Cross-correlation function and transfer function*

% Q1: Cross-correlation functions

% Create a binary representation matrix for all trials
binaryrep = zeros(100, 2 * round(endtime2));
for i = 1:length(data)
    trial = data(i,1);
    timebin = round(data(i,2) * 2);
    binaryrep(trial, timebin) = 2;
end

trial1 = binaryrep(1,:);
trial20 = binaryrep(20,:);

[ccf1,lags] = xcorr(trial1, stim, 200, 'unbiased');
[ccf20,lags] = xcorr(trial20, stim, 200, 'unbiased');

% Compute average cross-correlation across trials
ccfavg = zeros(size(ccf1));
for i = 1:100
    ccfavg = ccfavg + xcorr(binaryrep(i,:), stim, 200, 'unbiased');
end
ccfavg = ccfavg / 100;

figure(9);
plot(lags, ccf1, 'LineWidth', 3);
hold on;
plot(lags, ccf20, 'LineWidth', 2);
plot(lags, ccfavg, 'LineWidth', 1.5);
xlabel('lag'); ylabel('cross-correlation');
title('Cross-correlation functions of trials and stimulus');
legend('Trial 1', 'Trial 20', 'Trial Average', 'Location', 'NE');
hold off;

% Q2: Cross-spectrum

[pxy, f]=cpsd(stim, trial1, bartlett(2048), 1024, 2048, 2000);

figure(10);
plot(f, pxy, 'LineWidth', 1.5);
xlim([0, 100]);
xlabel('frequency (Hz)'); ylabel('cross-spectrum');
title('Sample of cross-spectrum at 0-100 Hz');

gain = abs(pxy);
gain = gain / gain(1,1); % normalize
phase = angle(pxy);

figure(11);
subplot(2,1,1);
plot(f, gain, 'LineWidth', 1);
title('Gain and phase as a function of frequency at 0-100 Hz');
xlabel('frequency (Hz)'); ylabel('Gain (spk/s)/(deg/s)');
xlim([0, 200]);
subplot(2,1,2);
plot(f, phase, 'LineWidth', 1);
xlabel('frequency (Hz)'); ylabel('Phase (deg)');
xlim([0, 200]);


%%
% *C) Signal-to-Noise Ratio*

% Q1: Signal-to-Noise

trialavg = mean(binaryrep);
[Presp,f] = pwelch(trialavg, bartlett(2048), 1024, 2048, 2000);

% Compute power spectrum of the noise
Pnoise = zeros(size(Presp));
for i = 1:100
   trial = binaryrep(i,:) - trialavg;
   Pnoise = Pnoise + pwelch(trial, bartlett(2048), 1024, 2048, 2000);
end
Pnoise = Pnoise / 100;

figure(12);
hold on;
plot(f, Presp, 'LineWidth', 1.5);
plot(f, Pnoise, 'LineWidth', 1.5);
xlim([0, 300]);
xlabel('frequency (Hz)'); ylabel('Power (spk^2/s)');
title('Power spectra of response and noise at 0-300 Hz');
legend('Response', 'Noise', 'Location', 'NE');

SNR = Presp ./ Pnoise;

figure(13);
subplot(2,1,1);
hold on;
plot(f, gain, 'LineWidth', 1.5);
plot(f, SNR, 'LineWidth', 1.5);
set(gca, 'LineWidth', 1.5, 'yscale', 'log')
xlabel('Freq (Hz)')
ylabel('Power (spk^2/s)')
title('Gain and SNR as Function of Frequency')
legend('Gain', 'SNR', 'Location', 'NE');

subplot(2,1,2);
hold on;
plot(f, gain, 'LineWidth', 1.5);
plot(f, SNR, 'LineWidth', 1.5);
set(gca, 'LineWidth', 1.5, 'yscale', 'log')
xlim([0, 60])
xlabel('Freq (Hz)')
ylabel('Power (spk^2/s)')
title('Sample at 0-60 Hz')
legend('Gain', 'SNR', 'Location', 'NE');