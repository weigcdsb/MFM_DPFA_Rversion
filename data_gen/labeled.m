usr_dir = "D:\github";
addpath(genpath(usr_dir + "\MFM_DPFA"));
%%

close all
rng(2)
n = 5;
nClus = 2;
N = n*nClus;
p = 2;
T = 1000;

Lab = repelem(1:nClus, n);

delta = randn(n*nClus,1)*0.5;
C_trans = zeros(n*nClus, p*nClus);
for k = 1:nClus
    NtmpIdx = (n*(k-1)+1):(n*k);
    delta(NtmpIdx) = delta(NtmpIdx) - mean(delta(NtmpIdx));
    pIdx = id2id(k, p);
    C_trans(NtmpIdx, pIdx) = randn(n, p);
end


X = ones(p*nClus, T)*Inf;
x0 = zeros(p*nClus,1);
Q0 = eye(nClus*p)*1e-2;
X(:,1) = mvnrnd(x0, Q0)';

% Generate X offline (A unspecified)
for i=1:(nClus*p)
    while(sum(X(i,:) > 1.5) > 1)
        k = ceil(rand()*20)+10;
%         disp(k)
        X(i,:) = interp1(linspace(0,1,k),randn(k,1)*0.2,linspace(0,1,T),'spline');
    end
end

% generate mu
mu = ones(nClus, T)*Inf;
for i = 1:nClus
    while(sum(mu(i,:) > log(15)) > 1 || sum(mu(i,:) < log(2)) == T)
        k = ceil(rand()*25)+10;
        mu(i,:) = interp1(linspace(0,1,k),randn(k,1)*0.7,linspace(0,1,T),'spline');
    end
end

mu = mu - mean(mu,2);

% let's generate lambda
logLam = zeros(n*nClus, T);
logLam(:,1) = delta + C_trans*X(:,1) + repelem(mu(:,1), n);


for t=2:T
    logLam(:, t) = delta + C_trans*X(:,t) + repelem(mu(:,t), n);
end
%
Y = poissrnd(exp(logLam));

subplot(1,2,1)
imagesc(Y)
subplot(1,2,2)
imagesc(exp(logLam))

%% output simulated data
usr_dir = 'D:\github';
r_wd = [usr_dir '\MFM_DPFA_Rversion\data_gen'];

writematrix(Y, [r_wd '\Y.csv'])
writematrix(Lab, [r_wd '\lab.csv'])

writematrix(delta, [r_wd '\delta.csv'])
writematrix(mu, [r_wd '\mu.csv'])
writematrix(X, [r_wd '\X.csv'])
writematrix(logLam, [r_wd '\logLam.csv'])
