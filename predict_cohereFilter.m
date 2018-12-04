function predicted = predict_cohereFilter(brainSignals,H,coeff,fs)

% using the filter and coefficients learned using the function
% build_cohereFilter, this function allows to predict features of an
% stimulus given the brain signals
% H and coeff are calculated using build_cohereFilter
% brainSignas is a matrix of time-electrodes-trials and refers to the brain
% signals measured as response to the stimulus
% predicted is a matrix of time-1-trials with the predicted output
% Jaime F. Delgado Saa, 2018

nfft = 4*fs;
brainSignals_ft = fft(brainSignals,nfft,1);
inputData = ifft(bsxfun(@times,brainSignals_ft,H),nfft,1);
inputData = inputData(1:size(brainSignals,1),:,:);

predicted = zeros(size(brainSignals,1),1,size(brainSignals,3));
for i = 1:size(brainSignals,3)
    predicted(:,1,i) = [ones(size(inputData,1),1) inputData(:,:,i)]*coeff;
end

return

