%% Setup for Question 1: Temporal Waveforms

fs=10000; %10 kHz sampling rate
dt=1/fs; % 1 / frequency gives sampling rate
t=0:dt:10; % 10 seconds of time
colors = ['#0072BD'; '#D95319'; '#EDB120'; '#7E2F8E'; '#77AC30'; '#4DBEEE'; '#A2142F'; '#000000'; '#FF00FF'; '#964B00'];
% freq = [100, 150, 200, 250, 300, 400, 500, 600]; % Hz
% amp = [1, 2, 3, 4, 5, 10, 20, 100];
% phase = [0, 1/4, 1/2, 3/4, 1, 5/4, 3/2, 7/4, 2]; % pi

%%
% Varying frequency difference
freq = [100, 150, 200, 250, 300, 400, 500, 600];

for f = 1:8
    y1=sin(2*pi * t * 100); % a sinusoid going at 100 Hz
    y2=sin(2*pi * t * freq(f)); % a sinusoid going at variable frequency
    ysum = y1+y2;
    
    figure(1);
    subplot(4,2,f);
    plot(t, ysum, 'LineWidth', 1, 'Color', colors(f,:));
    set(gca, 'xlim', [0 .05], 'ylim', [-2.5, 2.5]); % zoom into x-axis so that you can see something
    xlabel('time (s)');
    ylabel('amp');
    titleline = ['100 Hz + ', num2str(freq(f)), ' Hz'];
    title(titleline);
    sgtitle('Varying frequency difference between two summed sinusoids', 'FontSize', 20)
end

% Varying amplitude
amp = [1, 2, 3, 4, 5, 10, 20, 100];

for a = 1:8
    y1= 1 * sin(2*pi*t*100); % a sinusoid with amplitude 1 (100 Hz)
    y2= amp(a) * sin(2*pi*t*200); % a sinusoid with variable amplitude (200 Hz)
    ysum = y1+y2;
    
    figure(2);
    subplot(4,2,a);
    plot(t, ysum, 'LineWidth', 1, 'Color', colors(a,:));
    set(gca, 'xlim', [0 .05]); % zoom into x-axis so that you can see something
    xlabel('time (s)');
    ylabel('amp');
    titleline = ['Amplitude: 1 + ', num2str(amp(a))];
    title(titleline);
    sgtitle('Varying amplitude of two summed sinusoids (100 & 200 Hz)', 'FontSize', 20)
end

% Swap static and variable frequencies
for a = 1:8
    y1= 1 * sin(2*pi*t*200); % a sinusoid with amplitude 1 (200 Hz)
    y2= amp(a) * sin(2*pi*t*100); % a sinusoid with variable amplitude (100 Hz)
    ysum = y1+y2;
    
    figure(3);
    subplot(4,2,a);
    plot(t, ysum, 'LineWidth', 1, 'Color', colors(a,:));
    set(gca, 'xlim', [0 .05]); % zoom into x-axis so that you can see something
    xlabel('time (s)');
    ylabel('amp');
    titleline = ['Amplitude: 1 + ', num2str(amp(a))];
    title(titleline);
    sgtitle('Varying amplitude of two summed sinusoids (200 & 100 Hz)', 'FontSize', 20)
end

phase = [0, 1/4, 1/2, 3/4, 1, 5/4, 3/2, 7/4, 2];

% Varying phase of two 100 Hz sinusoids
for p = 1:9
    y1= sin(2*pi*t*100); % sinusoid with zero phase (100 Hz)
    y2= sin(2*pi*t*100 - phase(p)*pi); % sinusoid with variable phase (100 Hz)
    ysum = y1+y2;
    
    figure(4);
    if p == 1; subplot(5,2,[1 2]); else subplot(5,2,p+1); end
    plot(t, ysum, 'LineWidth', 1, 'Color', colors(p,:));
    set(gca, 'xlim', [0 .1], 'ylim', [-2.25, 2.25]); % zoom into x-axis so that you can see something
    xlabel('time (s)');
    ylabel('amp');
    titleline = ['Shifted right by: ', num2str(phase(p)), ' pi'];
    title(titleline);
    sgtitle('Varying phase of two summed sinusoids (100 Hz each)', 'FontSize', 20)
end

%Varying phase of 200 Hz sinusoid (+ 100 Hz constant sinusoid)
for p = 1:9
    y1= sin(2*pi*t*100); % sinusoid with zero phase (100 Hz)
    y2= sin(2*pi*t*200 - phase(p)*pi); % sinusoid with variable phase (200 Hz)
    ysum = y1+y2;
    
    figure(5);
    if p == 1; subplot(5,2,[1 2]); else subplot(5,2,p+1); end
    plot(t, ysum, 'LineWidth', 1, 'Color', colors(p,:));
    set(gca, 'xlim', [0 .1], 'ylim', [-2.25, 2.25]); % zoom into x-axis so that you can see something
    xlabel('time (s)');
    ylabel('amp');
    titleline = ['200 Hz shifted right by: ', num2str(phase(p)), ' pi'];
    title(titleline);
    sgtitle('Varying phase of 200 Hz sinusoid + 100 Hz constant sinusoid', 'FontSize', 20)
end

%% Setup for Question 3

fs=10000; %10 kHz sampling rate
dt=1/fs; % 1 / frequency gives sampling rate
endtime = 100; % sec
t=0:dt:endtime; % 50 seconds of time
fc=1000; fm=100; % carrier and modulation frequency
m=0.5; %modulation depth
y=10*sin(2*pi*fc*t).*(1+m*sin(2*pi*fm*t)); %a sinusoid going at 1000 Hz, with amplitude modulated at 100 Hz.
%soundsc(y,fs); % signal sound

figure(6);
plot(t,y)
set(gca, 'xlim', [0 0.025]);
xlabel('time (s)'); ylabel('amplitude');
mstitle = ['Sample of modulated tone (', num2str(fc), ' Hz carrier + ', num2str(fm), ' Hz modulator)'];
title(mstitle);

fouriertransform = fft(y);
n = length(y); % number of samples
f = (0:n-1)*(fs/n); % frequency range
power = abs(fouriertransform).^2/n; % power of the Fourier transform

figure(7);
plot(f,power)
set(gca, 'LineWidth', 1.5, 'yscale', 'log');
xlabel('Frequency (Hz)')
ylabel('Power')
xlim([0 2*fc])
pstitle = ['Power spectrum of ', num2str(fc), ' Hz carrier + ', num2str(fm), ' Hz modulator'];
title(pstitle);

y = iosr.auditory.gammatoneFast(y,fc,fs);% gammatone filter centered at 1 kHz = carrier frequency
y(y<0) = 0; % half-wave rectification
output = lowpass(y,500,fs); %low-pass filter

figure(8);
plot(t,output);
set(gca, 'xlim', [0 .15]);
xlabel('time (sec)'); ylabel('amplitude');
filttitle = ['Gammatone filtered, half-wave rectified, and low-pass filtered tone (', num2str(fc), ' Hz carrier + ', num2str(fm), ' Hz modulator)'];
title(filttitle)

% Generate a poisson process
spiketimes=genPoisson(output,dt);

% Spiketimes are converted to and placed in their corresponding bin number 
% by multiplying each time by bins/sec
spiketrain = zeros(1, round(endtime * fs)); % sample according to carrier freq
for i = 1:length(spiketimes)
    spiketrain(round(spiketimes(i) * fs)) = fs;
end

figure(9);
plot(spiketrain, 'LineWidth', 1);
ylim([-0.5*fs, 1.5*fs]);
xlabel('Bin (0.1 msec)'); ylabel('Bins/sec');
brtitle = ['Binary representation of spike train (', num2str(fc), ' Hz carrier + ', num2str(fm), ' Hz modulator)'];
title(brtitle);
 
[ps,f] = pwelch(spiketrain,bartlett(2048),1024,2048,fs);
figure(10);
plot(f, ps)
xlim([0 2*fc])
xlabel('Frequency (Hz)'); ylabel('Power')
pstitle = ['Power spectrum of ', num2str(fc), ' Hz carrier + ', num2str(fm), ' Hz modulator-driven spike train'];
title(pstitle);

%% Calculate phase-locking without plotting figures
fs=10000;
dt=1/fs; 
endtime = 100;
t=0:dt:endtime;
fc=100; fm=10; 
m=1;
y=10*sin(2*pi*fc*t).*(1+m*sin(2*pi*fm*t));

y = iosr.auditory.gammatoneFast(y,fc,fs);
y(y<0) = 0;
output = lowpass(y,500,fs);

r_avg = 0; pval_avg = 0;

for i = 1:10
    spiketimes=genPoisson(output,dt);
    % Calculate phase-locking
    [r,r_pvalue]=calculateLocking(spiketimes,fc); %locking to carrier frequency
    r_avg = r_avg + r;
    pval_avg = pval_avg + r_pvalue;
end
r_avg = r_avg/10
pval_avg = pval_avg/10

function [r,r_pvalue]= calculateLocking(spiketimes,freq)
    for ind=1:length(freq)
        phasevals=spiketimes*freq*2*pi; %convert spike-time to phase
        r=circ_r(phasevals'); %calling toolbox function for vector strength
        r_pvalue=circ_rtest(phasevals'); %calculating pvalue from rayleigh test
    end 
end

function spiketimes=genPoisson(y,dt)
    %generates poisson spike-train
    %draw exponential random variable (rv), mean 1. integrate driving function till
    %crosses that value. generate spike at that time. reset integral to zero,
    %draw rv again and restart the process.
    %this can be efficiently achieved by using histc, which calculates the
    %number of values that lie between the two spiketimes for the mean 1
    %poisson process. multiplying that by dt, gives the duration/inter-spike
    %interval.
    
    driver=cumsum(y)*dt; %running integral
    exps=exprnd(1,1,round(20*dt*length(y)*mean(y))); %some number well above number of expected spikes, which can be calculated by mean(output)*total time duration of simulation 
    exptimes=cumsum(exps); %sampled spike-train for poisson process with mean 1
    isis=dt*histcounts(driver,[0 exptimes]);

    finalspike=find(isis>0, 1, 'last');
    isis=isis(1:finalspike); %removing exptimes that are too long, and so give counts of zero.

    spiketimes=cumsum(isis); % spike-train for 1 trial of stimulus presentation
end
