function [H, coeff] = build_cohereFilter(brainSignals, stimulusfeat, fs)

%Build a coherence based filter as decribed in https://www.biorxiv.org/content/early/2018/11/28/481572
% brainSignas is a matrix of time - electrodes - trials, represention
% several realizations of a random process
% stimulusfeat is a matrix of time - 1 - trials,
% representing teh stimulus that drives the process measured by the
% brainSignals
% fs is the sample rate
% H is the filter (in frequency domain) calculated using the traiing set.
% Coeff determines the weights given to each channel for the prediction of
% the stimulus.
% Jaime F. Delgado Saa, 2018

nfft = 4*fs;
brainSignals_ft = fft(brainSignals,nfft,1);
stimulusfeat_ft = fft(stimulusfeat,nfft,1);
    
Pyx = bsxfun(@times,stimulusfeat_ft,conj(brainSignals_ft));
Px = sum(brainSignals_ft.*conj(brainSignals_ft),3);
H = sum(Pyx,3)./Px;

x_train = ifft(bsxfun(@times,brainSignals_ft,H),nfft,1);
x_train = x_train(1:size(brainSignals,1),:,:);
x_train = squeeze(reshape(permute(x_train,[1 3 2]),[],1,size(x_train,2)));

y_train = stimulusfeat(:);
coeff = regress(y_train,[ones(length(x_train),1) x_train]);
return

