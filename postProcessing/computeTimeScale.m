function Tu = computeTimeScale(time, vel, Nwindows)
    
% function that compute time scale of time series u:
% divide the time-series in N segments and average the N results

%Nwindows = 8;

W = floor(length(vel)/Nwindows);

for i = 1:Nwindows
    
    t = time((i-1)*W+1:i*W)-time((i-1)*W+1);
    u = vel((i-1)*W+1:i*W);
    N = length(t);
    
    %Ruu   = xcov(u);
    %RHOuu = Ruu(N:end)/Ruu(N);
    
    Nsamples = floor(1 * length(u) - 1);
    RHOuu = autocorr(u, Nsamples);
    
    %Tu_n_xcov(i) = t(min(find(RHOuu_xcov <= exp(-1))));
    Tu_n(i) = t(min(find(RHOuu <= exp(-1))));
    
%     % fit exponential
%     myFit   = fittype('exp(-x/a)', 'dependent',{'y'}, 'independent',{'x'}, 'coefficients',{'a'});
%     coeff   = fit(t(1:Nsamples+1), RHOuu, myFit, 'StartPoint', 1, 'Lower', 0);
%     Tu_n_fit(i) = coeff.a;

%     figure
%     set(gca, 'FontSize', 16);
%     hold on
%     plot(t, RHOuu, '.r')
%     plot(t, exp(-t/Tu_n(i)), 'b')
%     plot(t, exp(-t/Tu_n_fit(i)), 'k')
%     legend('data', sprintf('exp(-1): Tu = %6.3f', Tu_n(i)), sprintf('fit: Tu = %6.3f', Tu_n_fit(i)));
%     xlabel('time [s]'); ylabel('rho');
%     xlim([0 1])
end
Tu = mean(Tu_n);