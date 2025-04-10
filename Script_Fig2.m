% =================== This Script reproduces Fig.2 in [1]=================%
% The scenario is a MIMO link with an unblocked direct link (NLoS, Rayleigh) and LoS channels
% through RIS/BD-RIS. The equivalent channel is H = Hd+(fa*fd')*Theta*(gd*ga').
% We compare the rate improvement provided by :
% i)  the optimal (proposed) BD-RIS [1] with optimal Tx covariance matrix; 
% ii) the optimal (proposed) BD-RIS [1] with isotropic Tx covariance matrix; 
% iii) the optimal single-strean transmission with optimized beamformers
% and BD-RIS;
% iv) BD-RIS with random coefficients
%
% The forward and backward channels , F and G, are rank-1 (pure LoS)
% channels
%
% [1] I. Santamaria, J. Gutierrez, M. Soleymani, E. Jorswieck 
% "Rate Analysis and Optimization of LoS Beyond Diagonal RIS-assisted MIMO Systems
% IEEE Comm Letters, 2025.

clear;
format compact;

M = 2:2:128;          % number of BD-RIS elements
Ntx = 4;              % number of transmit antennas
Nrx = 4;              % number of receive antennas (we asumme here Nr<=Nt)

PtdBm = 30;
Pt = 10.^(PtdBm/10); % total Tx power
numMC = 1e3;         % number of Monte Carlo simulations (increase this number for smoother curves)
B = 20;                       % Bandwidth MHz
NF = 10;                      % Noise Factor in dBs
noiseVariancedBm = -174 + 10*log10(B*10^6) + NF;
sigma2n = 10^(noiseVariancedBm/10);       % additive noise variance
Rxx_isotropic = (Pt/Ntx)*eye(Ntx);        % isotropic Tx covariance matrix
Rxxsqrt_isotropic = sqrt(Pt/(Ntx*sigma2n))*eye(Ntx); % Rxx^{-1/2} (note: we have included the noise)

%% Channel parameters
channelparams = struct;
channelparams.blocked = 0;         % Set to 1 if direct channel is blocked
channelparams.RiceRIS = Inf;       % Rician factor for the channel btw RIS and BS (if Inf-> pure LoS channels)
channelparams.RiceDirect = 0;      % Rician factor for the channel btw RIS and UEs (if 0 -> Rayleigh fading)
channelparams.pl_0 = -28;          % Path loss at a reference distance (d_0)
channelparams.alpha_RIS = 2;       % Path loss exponent for the RIS links
channelparams.alpha_direct = 3.75; % Path loss exponent for the direct links
channelparams.ray_fading = 0;      % Set to 1 if all channels Rayleigh

%% ------------ Scenario Definition -----------------%%
% --- Position of the users/RIS (units in meters)----%%
R = 200;         % Radius of the cell
% Tx position
x_tx = 0;
y_tx = 0;
z_tx = 3;

% BD-RIS position
x_ris = 20;
y_ris = 20;
z_ris = 20; %located higher than the Tx/Rx

%Rx position  
x_rx = R;
y_rx = R;
z_rx = 1.5;

PosRIS_XYZ = [x_ris, y_ris, z_ris];
PosTx_XYZ = [x_tx y_tx z_tx];
PosRx_XYZ = [x_rx y_rx z_rx];

% To store results
Rate_BDRIS = zeros(size(M));
Rate_BDRISrnd = zeros(size(M)); 
Rate_BDRIS1stream = zeros(size(M)); 
Rate_BDRIS_opt = zeros(size(M));

%% Loop starts here
for mm = 1:length(M)
    disp(['Number of BD-RIS elements = ' num2str(M(mm))])
    for tt = 1:numMC
       
        %% Generate channels for the direct MIMO link and for the LoS (rank-one) channels        
        [H,G,F] = ChannelsMIMO(M(mm),Nrx,Ntx,PosTx_XYZ, PosRx_XYZ,PosRIS_XYZ,channelparams);  % With this version the RiceanFactor also affects the direct links

        %% Extract LoS Channels 
        [UG,DG,VG] = svd(G);
        gd = sqrt(DG(1,1))*UG(:,1);      %we split the channel gain btw ga and gd
        ga = sqrt(DG(1,1))*VG(:,1);
        [UF,DF,VF] = svd(F);
        fa = sqrt(DF(1,1))*UF(:,1);      %we split the channel gain btw far fdr
        fd = sqrt(DF(1,1))*VF(:,1);

        %% Optimal BD-RIS with isotropic Rxx
        % Find Z and gamma3 (as in Proposition 1 of [1])
        A = H*Rxxsqrt_isotropic;
        g = Rxxsqrt_isotropic*gd;
        E = eye(Nrx) + A*A';
        gamma1 = fa'*(E\fa);
        gamma3 = (g'*A')*(E\fa);
        gamma2 = (g'*A')*(E\(A*g));
        alphaBDRIS = norm(fd)*norm(ga);
       
        Z = real(abs(gamma3).^2 + gamma1*(norm(g)^2-gamma2));
        Delta= real(Z*alphaBDRIS^2 + 2*alphaBDRIS*abs(gamma3));

        RatewoBDRIS = log2(real(det(eye(Nrx) + H*Rxx_isotropic*H'/sigma2n)));
        Rate_BDRIS(mm) = Rate_BDRIS(mm) + RatewoBDRIS+ log2(1 + Delta);
        
        %% Random BDRIS
        Qrnd = orth(randn(M(mm),M(mm)) + 1i*randn(M(mm),M(mm)));
        Thetarnd = (Qrnd*Qrnd.');  %symetric and unitary
        Heq = H + F*Thetarnd*G';
        RatewBDRISrnd = log2(real(det(eye(Nrx) + Heq*Rxx_isotropic*Heq'/sigma2n)));
        Rate_BDRISrnd(mm) = Rate_BDRISrnd(mm)+RatewBDRISrnd;

        %% Optimal single-stream solution 
        wtx = orth(randn(Ntx,1)+1i*randn(Ntx,1));
        wrx = orth(randn(Nrx,1)+1i*randn(Nrx,1));
        Q = orth(randn(M(mm),M(mm)) + 1i*randn(M(mm),M(mm)));
        BDRIS = (Q*Q.');  %symmetric and unitary
        Heq = H + F*BDRIS*G';
        heq = wrx'*Heq*wtx;
        Rateaux = log2(real(1 + Pt*abs(heq)^2/sigma2n));
        CMIMObeam = Rateaux;

        nitermax = 100;
        threshold = 1e-3;
        true = 1;
        iter = 0;
        while true ==1
            iter = iter +1;
            % Find optimal BD-RIS for fixed wtx wrx
            heq = wrx'*H*wtx;
            hr = F'*wrx;
            ht = G'*wtx;
            Ag = (hr*ht'+ conj(hr*ht')')/2;   % complex and symmetric
            [U,~,V] = svd(Ag);
            BDRIS = exp(1i*angle(heq))*[U(:,1:2) conj(V(:,3:end))]*V';
            % Find wtx, wrx for fixed BDRIS
            Heq = H + F*BDRIS*G';
            [U,D,V] = svd(Heq);
            wrx = U(:,1);
            wtx = V(:,1);
            heq = wrx'*Heq*wtx;
            Rateaux = log2(real(1 + Pt*abs(heq)^2/sigma2n));
            CMIMObeam = [CMIMObeam, Rateaux]; %#ok<AGROW>
            % Check convergence
            if (CMIMObeam(end)-CMIMObeam(end-1)<threshold)||(iter>=nitermax)
                true = 0;
            end

        end
        heq = wrx'*(H + F*BDRIS*G')*wtx;
        RatewBDRIS1stream = log2(real(1 + Pt*abs(heq)^2/sigma2n));
        Rate_BDRIS1stream(mm) = Rate_BDRIS1stream(mm)+RatewBDRIS1stream;

        %% Optimal BD-RIS and optimal Rxx 
        Q = orth(randn(M(mm),M(mm)) + 1i*randn(M(mm),M(mm)));
        BDRISopt = (Q*Q.');  %symmetric and unitary
        Heq = H + F*BDRISopt*G';
        Si = sigma2n*eye(Nrx);
        [Rxx,~,~] = OptTransmitCovMatrix(Heq,Si,Pt);
        Rateaux = log2(real(det(eye(Nrx) + Heq*Rxx*Heq'/sigma2n)));
        Copt = Rateaux;

        nitermax = 100;
        threshold = 1e-5;
        true = 1;
        iter = 0;
        while true ==1
            iter = iter +1;
            % Find optimal BD-RIS for fixed Rxx
            Rxxsqrt = sqrtm(Rxx);
            A = H*Rxxsqrt;
            g = Rxxsqrt*gd;
            E = eye(Nrx) + A*A';
            gamma3 = (g'*A')*(E\fa);
            Ag = (fd*ga'+ conj(fd*ga')')/2;   % complex and symmetric
            [U,~,V] = svd(Ag);
            BDRISopt = exp(-1i*angle(gamma3))*[U(:,1:2) conj(V(:,3:end))]*V';
            % Find Rxx for fixed BDRIS
            Heq = H + F*BDRISopt*G';
            [Rxx,~,~] = OptTransmitCovMatrix(Heq,Si,Pt);
            Rateaux = log2(real(det(eye(Nrx) + Heq*Rxx*Heq'/sigma2n)));
            Copt = [Copt, Rateaux]; %#ok<AGROW>
            % Check convergence
            if (Copt(end)-Copt(end-1)<threshold)||(iter>=nitermax)
                true = 0;
            end

        end
        [Rxx,Vxx,pxx] = OptTransmitCovMatrix(Heq,Si,Pt);
        Rate_BDRIS_opt(mm) = Rate_BDRIS_opt(mm)+log2(real(det(eye(Nrx) + Heq*Rxx*Heq'/sigma2n)));
        
    end
end

Rate_BDRIS =  Rate_BDRIS/numMC;
Rate_BDRISrnd =  Rate_BDRISrnd/numMC;
Rate_BDRIS1stream = Rate_BDRIS1stream/numMC; 
Rate_BDRIS_opt = Rate_BDRIS_opt/numMC;

%% Parameters for figures
fs = 12;   % fontsize
lw = 1.5;  % linewidth
ms = 8;    % markersize
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
figure(1);clf;
plot(M,Rate_BDRIS_opt,'k-o','MarkerSize',ms,'LineWidth',lw, 'MarkerIndices',[2,32,64]);
hold on
plot(M,Rate_BDRIS,'r-s','MarkerSize',ms,'LineWidth',lw,'MarkerIndices',[2,32,64]);
plot(M,Rate_BDRIS1stream,'m-*','MarkerSize',ms,'LineWidth',lw,'MarkerIndices',[2,32,64]);
plot(M,Rate_BDRISrnd,'b-d','MarkerSize',ms,'LineWidth',lw,'MarkerIndices',[2,32,64]);
legend('BD-RIS opt','BD-RIS isotropic','BD-RIS single-stream','BD-RIS random');
xlabel('M (number of elements)');
ylabel('Rate (bps/Hz)');
hold off

