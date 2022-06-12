%% Load data
load DataPFC.mat

%% 1. Correlation matrices
Cpre = corrcoef(Qpre);
Crun = corrcoef(Qrun);
Cpost = corrcoef(Qpost);

f1 = figure(1);
f1.Position = [60 1000 300 600];
sgtitle('Correlation matrices')
subplot(3,1,1)
imagesc(Cpre); colorbar; axis square; title('PRE')
subplot(3,1,2)
imagesc(Crun); colorbar; axis square; title('RUN')
subplot(3,1,3)
imagesc(Cpost); colorbar; axis square; title('POST')

% The diagonal of the correlation matrix only has values of 1 because these 
% are the correlations between each neuron and itself, so they obviously
% have perfect correlation. The range of correlation values is between 0
% and 1 because each pair of neurons can correlate with each other well
% (closer to 1) or not well (closer to 0). There are no negative
% correlations here. 


%% 2. Cells 20 and 26

xc_pre = xcorr(Qpre(:,20), Qpre(:,26), 800);
xc_run = xcorr(Qrun(:,20), Qrun(:,26), 800);
[xc_post, lags] = xcorr(Qpost(:,20), Qpost(:,26), 800);
lags = lags / 10;

figure(2); hold on
plot(lags, xc_pre, 'LineWidth', 1)
plot(lags, xc_run, 'LineWidth', 1)
plot(lags, xc_post, 'LineWidth', 1)
title('Cross-correlation between cells 20 and 26 for each epoch') 
xlabel('lag (sec)'); ylabel('cross-correlation')
legend('PRE', 'RUN', 'POST')

% The cross-correlation between the two neurons in PRE (before exploration)
% is understandably close to zero, because nothing has been learned yet.
% During RUN (exploration), the cross-correlation becomes extremely high,
% as these two neurons learn to fire with respect to the same stimulus. 
% During wakefulness, neural activity is especially high and 
% well-correlated. Finally, during POST-exploration sleep, we see that the
% neurons are now correlated at zero lag. This is because during replay,
% the memory from wakefulness is replayed and consolidated; the neurons
% that learned to fire together during RUN are firing again during the POST
% replay.

%% 3. PCA weights

PCweights = pcacov(Crun);
PC1 = PCweights(:,1);
PCscore = zscore(Qrun)*PC1;
scoreRUN = PCscore;

%% 4. Firing rates

[len, num] = size(Qrun);
len = len / 10; % in seconds
runFR = zeros(1,num);

for n = 1:num
    runFR(n) = sum(Qrun(:,n)) / len;
end

figure(3)
scatter(runFR,abs(PC1))
ylabel('Absolute value of weights in PC1'); xlabel('Firing rates (Hz)')
title('Average firing rates against absolute value of PC1 weights')

% What can we conclude? The plot does not seem to demonstrate much of a
% correlation between the absolute value of the weights and neuron firing
% rates. There are a good number of points clustered at the origin (or with
% low firing rates in general), suggesting that most neurons were not very
% active during the exploration epoch. This would explain why there isn't 
% significant correlation between many of the neurons, as seen in the
% correlation matrices. Some neurons had high firing rates, which
% correlated with lower absolute value weights.

%% 5. Quantifying sleep reactivation

scorePRE = zscore(Qpre)*PC1;
scorePOST = zscore(Qpost)*PC1;

%% 6. Reactivation strength

singleNpre = (zscore(Qpre).^2)*(PC1.^2);
singleNpost = (zscore(Qpost).^2)*(PC1.^2);

reactPRE = (scorePRE.^2) - singleNpre;
reactPOST = (scorePOST.^2) - singleNpost;

f4 = figure(4);
f4.Position = [100 800 800 300];
tPRE = 1:length(reactPRE); tPRE = tPRE / 10; % convert to seconds
tPOST = 1:length(reactPOST); tPOST = tPOST / 10;
sgtitle('Reactivation strengths')
subplot(1,2,1)
plot(tPRE, reactPRE,'k'); axis tight
ylim([-30 160])
xlabel('Time (s)'); ylabel('Reactivation strength'); title('Pre')
subplot(1,2,2)
plot(tPOST, reactPOST,'k'); axis tight
xlabel('Time (s)'); ylabel('Reactivation strength'); title('Post')
ylim([-30 160])

% My plot is very similar to the white/black plot in figure in the
% assignment, with corresponding peaks, valleys, and reactivation patterns.
% The PRE reactivation strengths are all fairly weak, since the patterns
% have not yet been learned during exploration or reinforced during POST 
% sleep replay. Accordingly, the POST reactivation strengths are much
% stronger, especially visible in the 200-600 and 1000-1200 second
% intervals.

% As mentioned, the large positive peaks indicate stronger neuronal 
% reactivation POST exploration, meaning that the neurons are more likely
% to be reactivated after the learning process. This makes sense, because
% otherwise the animal does not learn anything. We can conclude that POST
% is characterized by stronger and more frequent reactivation strengths 
% than PRE. 

% We can also conclude from the figure in the assignment that the higher 
% reactivation strengths coincide with SWS (slow-wave sleep), especially in 
% the POST epoch (during replay). Non-SWS does not seem to significantly 
% increase reactivation strength, even in the POST epoch. 

