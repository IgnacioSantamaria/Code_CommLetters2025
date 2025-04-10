% =================== This Script reproduces Fig.3 in [1]=================%
% The scenario is a MIMO link with an unblocked direct link (NLoS, Rayleigh) and LoS channels
% through RIS/BD-RIS. The equivalent channel is H = Hd+(fa*fd')*Theta*(gd*ga').
% We compare the rate improvement when the forward and backward channels
% are Ricean with Factor K.
%
% We compare the performance against several competing schemes described in
% [1]. We have not included in the simulation the algorithm of our SPACWC 2024
% paper to make the simulation faster. If needed, the code for this method
% can also be found at https://github.com/IgnacioSantamaria
%
% [1] I. Santamaria, J. Gutierrez, M. Soleymani, E. Jorswieck 
% "Rate Analysis and Optimization of LoS Beyond Diagonal RIS-assisted MIMO Systems
% IEEE Comm Letters, 2025.

clear;
format compact;

M = 64;               % number of BD-RIS elements
Ntx = 4;              % number of transmit antennas
Nrx = 4;              % number of receive antennas (we asumme here Nr<=Nt)

PtdBm = 10;
Pt = 10.^(PtdBm/10);  % total Tx power
B = 20;                       % Bandwidth MHz
NF = 10;                      % Noise Factor in dBs
noiseVariancedBm = -174 + 10*log10(B*10^6) + NF;
sigma2n = 10^(noiseVariancedBm/10);       % additive noise variance
numMC = 1e3;                  % number of Monte Carlo simulations (increase this number for smoother curves)
RiceFactor = 0:11;            % Ricean factors for the RIS channels
snr  = Pt/(Ntx*sigma2n);
%% Channel parameters
channelparams = struct;
channelparams.blocked = 0;         % Set to 1 if direct channel is blocked
channelparams.RiceRIS = Inf;       % Rician factor for the channel btw RIS and BS (if Inf-> pure LoS channels)
channelparams.RiceDirect = 0;      % Rician factor for the channel btw RIS and UEs (if 0 -> Rayleigh fading)
channelparams.pl_0 = -28;          % Path loss at a reference distance (d_0)
channelparams.alpha_RIS = 2;       % Path loss exponent for the RIS links
channelparams.alpha_direct = 3.75; % Path loss exponent for the direct links
channelparams.ray_fading = 0;      % Set to 1 if all channels Rayleigh

%% ------------ Scenario Definition ----------
% --- Position of the users/RIS (units in meters)
R = 200;         % Radius of the cell
% Tx position
x_tx = 0;
y_tx = 0;
z_tx = 3;

% BD-RIS position
x_ris = 20;
y_ris = 20;
z_ris = 20; % located higher than the Tx/Rx

%Rx position  
x_rx = R;
y_rx = R;
z_rx = 1.5;

PosRIS_XYZ = [x_ris, y_ris, z_ris];
PosTx_XYZ = [x_tx y_tx z_tx];
PosRx_XYZ = [x_rx y_rx z_rx];


%% To store the rates
Cave_OPA_BDRISnonrec = zeros(size(RiceFactor));
Cave_UPA_BDRIS = zeros(size(RiceFactor));
Cave_UPA_BDRISit = zeros(size(RiceFactor));
Cave_UPA_RIS = zeros(size(RiceFactor));
Cave_rnd_BDRIS = zeros(size(RiceFactor));
Cave_rnd_RIS = zeros(size(RiceFactor));
Cave_noRIS = zeros(size(RiceFactor));
Cave_optBDRIS = zeros(size(RiceFactor));

for mm = 1:length(RiceFactor)

    disp(['Rice Factor = ' num2str(RiceFactor(mm))])
    channelparams.RiceRIS =RiceFactor(mm);    % Channel Rice factor

    for tt = 1:numMC

        %% Generate channels for the direct MIMO link and for the LoS (rank-one) channels
        [Hd,G,F] = ChannelsMIMO(M,Nrx,Ntx,PosTx_XYZ, PosRx_XYZ,PosRIS_XYZ,channelparams);  % With this version the RiceanFactor also affects the direct links

        %% Extract the rank-one components (in case a Rician channel model is used, this is the dominant path)
        [UG,DG,VG] = svd(G);
        gd = sqrt(DG(1,1))*UG(:,1);      % we split the channel gain btw ga and gd
        ga = sqrt(DG(1,1))*VG(:,1);
        [UF,DF,VF] = svd(F);
        fa = sqrt(DF(1,1))*UF(:,1);       %we split the channel gain btw far fdr
        fd = sqrt(DF(1,1))*VF(:,1);

        %% Optimal BDRIS (LoS) + Waterfilling
        T = fd*ga' + (fd*ga').';    % Symmetric rank-two matrix
        %% Takagi
        [Uaux,~,Vaux] = svd(T);
        U1 = Uaux(:,1:2);
        V1 = Vaux(:,1:2);
        V2 = Vaux(:,3:end);
        Qrot = orth(randn(M-2,M-2) + 1i*randn(M-2,M-2));          % M-2 x M-2 rnd unitary matrix
        ThetaUPA_BDRIS = U1*V1' + conj(V2)*conj(Qrot)*Qrot'*V2';  % unitary + symmetric BDRIS
        
        alphaBDRIS = real(fd'*ThetaUPA_BDRIS*ga);
        thetadet = -angle(gd'*Hd'*((eye(Nrx)+snr*(Hd*Hd'))\fa));
        ThetaUPA_BDRIS = exp(1i*thetadet)*ThetaUPA_BDRIS;
        HUPA_BDRIS = Hd + F*ThetaUPA_BDRIS*G';
        [Qdet,~,~] = OptTransmitCovMatrix(HUPA_BDRIS,sigma2n*eye(Nrx),Pt);   % Optimal Tx cov. matrix
        Cave_UPA_BDRIS(mm) = Cave_UPA_BDRIS(mm) + log2(real(det(eye(Nrx)+(HUPA_BDRIS*Qdet*HUPA_BDRIS')/sigma2n)));

        %% Optimal RIS (LoS) + Waterfilling
        ThetaUPA_RIS = diag(exp(-1i*angle(conj(fd).*ga)));  % optimal DRIS
        alphaDRIS = real(fd'*ThetaUPA_RIS*ga);
        thetadet = -angle(gd'*Hd'*((eye(Nrx)+snr*(Hd*Hd'))\fa));
        ThetaUPA_RIS = exp(1i*thetadet)*ThetaUPA_RIS;
        HUPA_RIS = Hd + F*ThetaUPA_RIS*G';
        [Qdet,~,~] = OptTransmitCovMatrix(HUPA_RIS,sigma2n*eye(Nrx),Pt);  % Optimal Tx cov. matrix
        Cave_UPA_RIS(mm) = Cave_UPA_RIS(mm) + log2(real(det(eye(Nrx)+(HUPA_RIS*Qdet*HUPA_RIS')/sigma2n)));

        %% Random BDRIS + Waterfilling
        Qrnd = orth(randn(M,M) + 1i*randn(M,M));
        Thetarnd = (Qrnd*Qrnd.');   % symmetric and unitary
        Hrnd = Hd + F*Thetarnd*G';
        [Qdet,~,~] = OptTransmitCovMatrix(Hrnd,sigma2n*eye(Nrx),Pt); % Optimal Tx cov. matrix
        Cave_rnd_BDRIS(mm) = Cave_rnd_BDRIS(mm) + log2(real(det(eye(Nrx)+(Hrnd*Qdet*Hrnd')/sigma2n)));

        %% Random RIS + Waterfilling
        Thetarnd = diag(exp(1i*2*pi*rand(M,1)));
        Hrnd = Hd + F*Thetarnd*G';
        [Qdet,~,~] = OptTransmitCovMatrix(Hrnd,sigma2n*eye(Nrx),Pt); % Optimal Tx cov. matrix
        Cave_rnd_RIS(mm) = Cave_rnd_RIS(mm) + log2(real(det(eye(Nrx)+(Hrnd*Qdet*Hrnd')/sigma2n)));

        %% No RIS + Waterfilling
        [QnoRIS,~,~] = OptTransmitCovMatrix(Hd,sigma2n*eye(Nrx),Pt);
        Cave_noRIS(mm) = Cave_noRIS(mm) + log2(real(det(eye(Nrx)+(Hd*QnoRIS*Hd')/sigma2n)));

        %% Unitary BD-RIS (not symmetric) + Waterfilling (non-reciprocal channel!)
        [UF,DF,VF] = svd(F);
        [UG,DG,VG] = svd(G);
        Thetauni = VF*VG';
        Huni = Hd + F*Thetauni*G';
        [Quni,Vdet,pdet] = OptTransmitCovMatrix(Huni,sigma2n*eye(Nrx),Pt);  % Optimal Tx cov. matrix
        Cave_OPA_BDRISnonrec(mm) = Cave_OPA_BDRISnonrec(mm) + log2(real(det(eye(Nrx)+(Huni*Quni*Huni')/sigma2n)));

    end
end
Cave_OPA_BDRISnonrec =  Cave_OPA_BDRISnonrec/numMC;
Cave_UPA_BDRIS = Cave_UPA_BDRIS/numMC;
Cave_UPA_RIS = Cave_UPA_RIS/numMC;
Cave_rnd_BDRIS =  Cave_rnd_BDRIS/numMC;
Cave_rnd_RIS =  Cave_rnd_RIS/numMC;
Cave_noRIS =  Cave_noRIS/numMC;


%% Parameters for figures
fs = 12;   % fontsize
lw = 1.5;  % linewidth
ms = 8;    % markersize
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
figure(1);clf;
plot(RiceFactor,Cave_OPA_BDRISnonrec,'b:*','MarkerSize',ms,'LineWidth',lw);
hold on;
plot(RiceFactor,Cave_UPA_BDRIS,'m--d','MarkerSize',ms,'LineWidth',lw);
plot(RiceFactor,Cave_UPA_RIS,'k--+','MarkerSize',ms,'LineWidth',lw);
plot(RiceFactor,Cave_rnd_BDRIS,'m:hexagram','MarkerSize',ms,'LineWidth',lw);
plot(RiceFactor,Cave_rnd_RIS,'y:*','MarkerSize',ms,'LineWidth',lw);
plot(RiceFactor,Cave_noRIS,'g--v','MarkerSize',ms,'LineWidth',lw);
legend('BD-RIS non-rec. (NLoS)','BD-RIS (LoS)', 'RIS (LoS)','Random BD-RIS','Random RIS','No RIS','Location','best');
xlabel('K (Ricean Factor)');
ylabel('Achievable Rate (bps/Hz)');
hold off

