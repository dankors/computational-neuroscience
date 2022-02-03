% One dimensional filtering operation for simulating a reverse correlation
% experiment.  The input a is filtered by the simulated neuron and the
% response is returned in r.
% Christopher Pack, 11/07

impulse = [0 0 0 1 0 0 0];
output = one_d(impulse);
figure(1)
plot(output)
xlabel('time'); ylabel('response')
title('Temporal impulse response')

% Cross-correlation between impulse and temporal response
figure(2)
hold on
[xc,lags] = xcorr(output,impulse);
plot(lags, xc);
[val, idx] = max(xc);
plot(lags(idx), val, 'o', 'Color', 'r'); %we can see the lag is 5 timesteps
xlabel('lag'); ylabel('Correlation coefficient')
title('Cross-correlation between impulse and temporal impulse response') 

% create random white noise input
inputlength = 1000;
stim = (-2).*rand(1,inputlength) + 1; % values between -1 and 1
r = one_d(stim);
[cc,lags] = xcorr(r,stim); % compute cross-correlation
ntrials = 999;

% average cross-correlation across n trials to smooth it out
for i = 1:ntrials
    stim = (-2).*rand(1,inputlength) + 1;
    r = one_d(stim);
    [xctrial,~] = xcorr(r,stim);
    cc = cc + xctrial;
end
cc = cc / (ntrials + 1);

 % extract relevant part of cross-correlation
start = find(lags == 0);
cc = cc(:, start:start+49);
lags = lags(:, start:start+49);
figure(3);
plot(lags, cc);
xlabel('lag'); ylabel('Correlation coefficient')
title('Average cross-correlation between stimulus and neuron response across 1000 trials') 

% normalize filter
filt = cc / (length(stim) * std(stim)^2); 
figure(4);
plot(lags, filt);
xlabel('lag'); ylabel('Correlation coefficient')
title('Average cross-correlation NORMALIZED') 

% Use convolution to calculate predicted output of normalized filter
pred = conv(stim, filt);
pred(pred < 0) = 0;
pred = [pred, 0];

figure(5);
hold on
plot(r, 'LineWidth', 1);
plot(pred, 'LineWidth', 1);
xlim([0 length(r)])
legend('Real output', 'Predicted output', 'Location', 'NE');
xlabel('Time'); ylabel('Response')
title('Real and predicted neuron output calculated by convolution') 

% Visualize static nonlinearity
figure(6); hold on
scatter(pred,r,15);
xlabel('predicted'); ylabel('real');
title('Visualized static nonlinearity')

[~, indexOfLhalf] = min(abs(r-50)); % half of rmax
Lhalf = pred(indexOfLhalf);
F = staticnl(pred, [100 Lhalf 0.22]); %rmax, Lhalf, g1 (slope)

% Uncomment next line to superimpose sigmoid on static nonlinearity
scatter(pred,F,15);
legend('Real vs predicted plot', 'Static nonlinearity', 'Location', 'E');

figure(7);
hold on
plot(r, 'LineWidth', 1.5);
plot(pred, 'LineWidth', 1);
plot(F, 'LineWidth', 1);
xlim([0 length(r)])
legend('Real output', 'Predicted output (linear)', 'Predicted output (nonlinear)', 'Location', 'NE');
xlabel('Time'); ylabel('Response')
title('Real and predicted (linear and nonlinear) neuron output') 

%% Testing the filter

% MSE
inputlength = [100, 500, 1000, 2000, 5000];
ntrials = 1000;
linerr = zeros(size(inputlength));
nonlinerr = zeros(size(inputlength));

for len = 1:5
    avglin = 0;
    avgnonlin = 0;

    for i = 1:ntrials
        stim = (-2).*rand(1,inputlength(len)) + 1;
        r = one_d(stim); % actual response
        [cc,~] = xcorr(r,stim);
        start = find(lags == 0);
        cc = cc(:, start:start+49);
        filt = cc / (length(stim) * std(stim)^2);
        pred = conv(stim, filt); % predicted response
        pred(pred < 0) = 0;
        pred = [pred, 0];
        avglin = avglin + immse(pred,r);
        
        [~, indexOfLhalf] = min(abs(r-50)); % half of rmax
        Lhalf = pred(indexOfLhalf);
        F = staticnl(pred, [100 Lhalf 0.22]); %rmax, Lhalf, g1 (slope)
        avgnonlin = avgnonlin + immse(F,r);
    end
    
    avglin = avglin / ntrials;
    avgnonlin = avgnonlin / ntrials;
    linerr(len) = avglin;
    nonlinerr(len) = avgnonlin;
end

figure(8); clf; hold on
plot(inputlength, linerr, '-s', 'LineWidth', 1.5);
plot(inputlength, nonlinerr, '-s', 'LineWidth', 1.5);
xlabel('trial length'); ylabel('MSE')
title('Mean-squared error by input length (1000 trial average)');
legend('Linear MSE', 'Nonlinear MSE', 'Location', 'E');



function F=staticnl(L, params)
    rmax = params(1);
	Lhalf = params(2);
    g1 = params(3);
	
	F = rmax ./ (1 + exp(-g1 * (L - Lhalf)));
end

function r=one_d(a)

    k=25;  %Temporal scale factor 
    slow_t=temp_imp_resp(5,k,0:.02:1);
    fast_t=temp_imp_resp(3,k,0:.02:1);

    b=slow_t + fast_t; % linear filter
    c=conv(b,a); % convolve filter with stimulus
    c=c./max(c); % normalize

    d=100./(1.0 + exp(10*(0.5-c))); % static nonlinearity
    r=d; % for continuous output

end

function time_response=temp_imp_resp(n,k,t)
    %time_response=temp_imp_resp(n,k,t)
    %
    %Produces a temporal impulse response function using the from from
    %figure 1 in Adelson & Bergen (1985)
    %
    %It's pretty much a difference of Poisson functions with different
    %time constants.

    time_response=(k*t).^n .* exp(-k*t).*(1/factorial(n)-(k*t).^2/factorial(n+2));
end
