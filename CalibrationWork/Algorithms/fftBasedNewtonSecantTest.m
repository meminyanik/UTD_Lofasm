isTest = false;
signalModel = 2; % 1: Navigation paper, 2: Thomson Book

%% Test Algorithm
if isTest
    Ts = 1/1000; % Fs = 100;
    n = ((0:1000)*Ts).';
    
    modelledF = @(b,f,th) b*exp(1i*2*pi*f*n + th);
    modelledFdiff_f = @(b,f,th) b*(1i*2*pi*n).*exp(1i*2*pi*f*n + th);
    modelledFdiff2_f = @(b,f,th) b*(1i*2*pi*n).*(1i*2*pi*n).*exp(1i*2*pi*f*n + th);
    
    fIn = 34.7;
    bIn = 8;
    thIn = 0.2*pi;
    
    signal = modelledF(bIn,fIn,thIn);
    
    SNRdB = 40; % dB
    SNR = 10.^(SNRdB/10);
    varNoise = var(signal)/SNR;
    noise = sqrt(varNoise/2)*(randn(length(n),1)+1i*randn(length(n),1));
    
    measuredF = signal + noise; % 0.9*randn(length(n),1);
    nFFT = 2*1024; % 2^8;
    A = fft(measuredF,nFFT)/length(n);
    A = A(1:nFFT/2);
    [peak,k] = max(abs(A));
    % plot((0:nFFT/2-1)/nFFT/Ts,abs(A))
    
    fEst = (k-1)/nFFT/Ts;
    fEst1 = (k-2)/nFFT/Ts;
    fEst2 = (k)/nFFT/Ts;
    
    % thEst = angle(peak);
    % bEst = abs(peak);
    
    bIn = 1;
    thIn = 0;
    minFunc = @(f) (measuredF - modelledF(bIn,f,thIn)) .* conj(measuredF - modelledF(bIn,f,thIn));
    
    JacobianMinFunc_f = @(f) -modelledFdiff_f(bIn,f,thIn) .* conj(measuredF - modelledF(bIn,f,thIn)) + ...
        (measuredF - modelledF(bIn,f,thIn)) .* conj(-modelledFdiff_f(bIn,f,thIn));
    
    JacobianMinFunc2_f = @(f) -modelledFdiff2_f(bIn,f,thIn) .* conj(measuredF - modelledF(bIn,f,thIn)) + ...
        -modelledFdiff_f(bIn,f,thIn) .* conj(-modelledFdiff_f(bIn,f,thIn)) + ...
        -modelledFdiff_f(bIn,f,thIn) .* conj(-modelledFdiff_f(bIn,f,thIn)) + ...
        (measuredF - modelledF(bIn,f,thIn)) .* conj(-modelledFdiff2_f(bIn,f,thIn));
    
    %% Search Routine
    %     fSearch = linspace(fEst1,fEst2,200);
    %     % fSearch = linspace(18,23,200);
    %     mR = zeros(1,length(fSearch));
    %     for n = 1:length(fSearch)
    %         x_k_1 = fSearch(n);
    %         J = JacobianMinFunc_f(x_k_1);
    %         fm = minFunc(x_k_1);
    %         % mR(n) = abs(-(J.' * J) \ J.' * fm);
    %         mR(n) = -(J.' * J) \ J.' * fm;
    %     end
    %     plot(fSearch,mR)
    
    %% Iterative Algorithm
    % [x,k] = secant(JacobianMinFunc_f,fEst1,fEst2);
    % [x_k,k,all_x] = newton(minFunc,JacobianMinFunc_f,fEst1);
    [f_k,k,all_f] = newton(JacobianMinFunc_f,JacobianMinFunc2_f,fEst);
    % figure;
    % plot(all_f)
    plot(abs(fIn-all_f))
    
    %% Latitude and Longitude Estimation
else
    if (signalModel == 1)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Sky Source Coordinates
        %-------------------------------------------------------------------------%
        raCygnusA = 19 + 59/60 + 28.3/3600; % in Hours
        decCygnusA = 40 + 44/60 + 02/3600; % in Degrees
        % Cygnus A: 19h59m28.3s +40d44m02s
        raCasA = 23 + 23/60 + 27.9/3600; % in Hours
        decCasA = 58 + 48/60 + 42/3600; % in Degrees
        % Cas A: 23h23m27.9s +58d48m42s
        %-------------------------------------------------------------------------%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Define Total Scanned Hour Duration and start time
        %-------------------------------------------------------------------------%
        % Assume 6 hours of data that integration time is 10s
        % GST = 10; % in Hours
        h = (0:50/3600:6).';
        %-------------------------------------------------------------------------%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Define Baselines in meters
        %-------------------------------------------------------------------------%
        % x: towards the east, y: toward the north, z: towards zenith
        x = 200;
        y = 30;
        z = 10;
        %-------------------------------------------------------------------------%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Calculate Baselines and Theoric Delay
        %-------------------------------------------------------------------------%
        % l is unknown; % latitude
        % L is unknown % Longitude
        
        X = x;
        Y = @(l) z*sind(l) + y*cosd(l);
        Z = @(l) z*cosd(l) - y*sind(l);
        dY = @(l) z*cosd(l) - y*sind(l);
        dZ = @(l) -z*sind(l) - y*cosd(l);
        
        dataRecordTime = [5,45,11,21,7,2016];
        dataRecordTime = num2cell(dataRecordTime);
        GST = calculate_GST(dataRecordTime{:});
        LST = @(L) mod(L/15+GST,24);

        dCygA = decCygnusA;
        rCygA = raCygnusA;
        % HCygA = @(L) LST(L) + h - rCygA;
        % dHCygA = 1/15;
        HCygA = @(L) (LST(L) + h - rCygA)*pi/12; % L is in degrees
        dHCygA = 1;
        
        wCygA = @(l,L) -cosd(dCygA)*sin(HCygA(L))*X + sind(dCygA)*Y(l) + cosd(dCygA)*cos(HCygA(L))*Z(l);
        dwCygA_l = @(l,L) sind(dCygA)*dY(l) + cosd(dCygA)*cos(HCygA(L))*dZ(l);
        dwCygA_L = @(l,L) -cosd(dCygA)*cos(HCygA(L))*X - cosd(dCygA)*sin(HCygA(L))*Z(l);
        
        %     dCasA = decCasA;
        %     rCasA = raCasA;
        %     % HCasA = @(L) LST(L) + h - rCasA;
        %     % dHCasA = 1/15;
        %     HCasA = @(L) L + h - rCasA;
        %     dHCasA = 1;
        %     wCasA = @(l,L) -cosd(dCasA)*sin(HCasA(L)*pi/12)*X + sind(dCasA)*Y(l) + cosd(dCasA)*cos(HCasA(L)*pi/12)*Z(l);
        %     dwCasA_l = @(l,L) sind(dCasA)*dY(l) + cosd(dCasA)*cos(HCasA(L)*pi/12)*dZ(l);
        %     dwCasA_L = @(l,L) -cosd(dCasA)*cos(HCasA(L)*pi/12)*dHCasA*pi/12*X - cosd(dCasA)*sin(HCasA(L)*pi/12)*dHCasA*pi/12*Z(l);
        %-------------------------------------------------------------------------%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Define Frequency
        %-------------------------------------------------------------------------%
        % fc = linspace(0,100,1024)*1e6;
        f = 50e6; % 50 MHz
        c = physconst('LightSpeed');
        %-------------------------------------------------------------------------%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Define Signal Model
        %-------------------------------------------------------------------------%
        modelledR = @(l,L) exp(1i*2*pi*f*wCygA(l,L)/c);
        modelledRdiff_l = @(l,L) (1i*2*pi*f*dwCygA_l(l,L)/c) .* exp(1i*2*pi*f*wCygA(l,L)/c);
        modelledRdiff_L = @(l,L) (1i*2*pi*f*dwCygA_L(l,L)/c) .* exp(1i*2*pi*f*wCygA(l,L)/c);
        %-------------------------------------------------------------------------%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Define Measurement
        %-------------------------------------------------------------------------%
        lIn = 34.5; % Input latitude
        LIn = -116; % Input LST that is a function of Longitude
        measuredR = modelledR(lIn,LIn) + 0.9*randn(length(h),1);
        %-------------------------------------------------------------------------%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Define Min Function
        %-------------------------------------------------------------------------%
        minFunc = @(l,L) (measuredR - modelledR(l,L)) .* conj(measuredR - modelledR(l,L));
        JacobianMinFunc_l = @(l,L) -modelledRdiff_l(l,L) .* conj(measuredR - modelledR(l,L)) + ...
            (measuredR - modelledR(l,L)) .* conj(- modelledRdiff_l(l,L));
        JacobianMinFunc_L = @(l,L) -modelledRdiff_L(l,L) .* conj(measuredR - modelledR(l,L)) + ...
            (measuredR - modelledR(l,L)) .* conj(- modelledRdiff_L(l,L));
        JacobianMinFunc = @(l,L) [JacobianMinFunc_l(l,L), JacobianMinFunc_L(l,L)];
        %-------------------------------------------------------------------------%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Matched Filter Process
        %-------------------------------------------------------------------------%
        %% Iteration 1
        lSearch = linspace(0,90,10);
        LSearch = linspace(-180,180,100);
        % lSearch = lIn;
        % LSearch = LIn;
        mR = zeros(length(lSearch),length(LSearch));
        for nl = 1:length(lSearch)
            for nL = 1:length(LSearch)
                lS = lSearch(nl);
                LS = LSearch(nL);
                MF = modelledR(lS,LS)';
                mR(nl,nL) = MF * measuredR;
            end
        end
        figure; mesh(LSearch,lSearch,abs(mR));
        %         % plot(lSearch,abs(mR));
        %         % plot(LSearch*12/pi,abs(mR));
        %
        %     [~,k] = max(abs(mR(:)));
        %     [k_l, k_L] = ind2sub(size(mR),k);
        %
        % %     lEst1 = lSearch(k_l-1)
        % %     lEst2 = lSearch(k_l+1)
        %
        %     % LEst = LSearch(k_L)*12/pi;
        %     LEst1 = LSearch(k_L-1)*12/pi;
        %     LEst2 = LSearch(k_L+1)*12/pi;
        
        %% Iteration 2
        %     lSearch = linspace(0,90,110);
        %     LSearch = LIn;
        %     mR = zeros(1,length(lSearch));
        %     for nl = 1:length(lSearch)
        %         lS = lSearch(nl);
        %         LS = LSearch;
        %         MF = modelledR(lS,LS)';
        %         mR(nl) = MF * measuredR;
        %     end
        %     figure; plot(lSearch,abs(mR));
        %
        %     [~,k_l] = max(abs(mR));
        %     lEst1 = lSearch(k_l-1)
        %     lEst2 = lSearch(k_l+1)
        %-------------------------------------------------------------------------%
        
        
        %% Newton Search
        %     lEst1 = lIn-0.5; lEst2 = lIn+0.5;
        %     LEst1 = LIn-0.1; LEst2 = LIn+0.1;
        
        %     lSearch = linspace(36,40,100);
        %     LSearch = linspace(16*pi/12,18*pi/12,100);
        %     % lSearch = lIn;
        %     % LSearch = LIn;
        %     mR = zeros(length(lSearch),length(LSearch));
        %     for nl = 1:length(lSearch)
        %         for nL = 1:length(LSearch)
        %             x_k_1 = [lSearch(nl);LSearch(nL)];
        %             x_k_1_c = num2cell(x_k_1);
        %             J = JacobianMinFunc(x_k_1_c{:});
        %             fm = minFunc(x_k_1_c{:});
        %             mR(nl,nL) = norm(-(J.' * J) \ J.' * fm);
        %             % mR(n) = -(J.' * J) \ J.' * fm;
        %         end
        %     end
        %     mesh(lSearch,LSearch*12/pi,mR)
        %     % figure; plot(lSearch,mR)
        %     % figure; plot(LSearch*12/pi,mR)
        
        %% Newton and Secant Search for One Parameter
        %     minFunc_lat = @(l) minFunc(l,LIn);
        %     JacobianMinFunc_lat = @(l) JacobianMinFunc_l(l,LIn);
        %     % [l_k,~,all_l] = newton(minFunc_lat,JacobianMinFunc_lat,35);
        %     [l_k,~,all_l] = secant(minFunc_lat,39,37);
        %     figure; plot(all_l)
        %
        %     minFunc_lon = @(L) minFunc(lIn,L);
        %     JacobianMinFunc_lon = @(L) JacobianMinFunc_L(lIn,L);
        %     % [L_k,~,all_L] = newton(minFunc_lon,JacobianMinFunc_lon,16*pi/12);
        %     [L_k,~,all_L] = secant(minFunc_lon,16*pi/12,17*pi/12);
        %     figure; plot(all_L*12/pi)
        
        %% Newton and Secant Search for Two Parameter
        [lL_k,~,all_lL] = newton(minFunc,JacobianMinFunc,[38.8;-140]);
        % [lL_k,~,all_lL] = secant(minFunc,[38.8;17.3*pi/12],[39.4;18.5*pi/12]);
        figure; plot(all_lL(1,:))
        figure; plot(all_lL(2,:))
        
    elseif (signalModel == 2)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Sky Source Coordinates
        %-------------------------------------------------------------------------%
        raCygnusA = 19 + 59/60 + 28.3/3600; % in Hours
        decCygnusA = 40 + 44/60 + 02/3600; % in Degrees
        % Cygnus A: 19h59m28.3s +40d44m02s
        raCasA = 23 + 23/60 + 27.9/3600; % in Hours
        decCasA = 58 + 48/60 + 42/3600; % in Degrees
        % Cas A: 23h23m27.9s +58d48m42s
        %-------------------------------------------------------------------------%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Define Total Scanned Hour Duration and start time
        %-------------------------------------------------------------------------%
        % Assume 6 hours of data that integration time is 10s
        % GST = 10; % in Hours
        h = (0:50/3600:12).';
        %-------------------------------------------------------------------------%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Define Baselines in meters
        %-------------------------------------------------------------------------%
        % x: towards the east, y: toward the north, z: towards zenith
        D = 200; % Distance of baseline
        el = 10; % Relative elevation of baseline in degrees
        az = 20; % Relative azimuth of baseline in degrees
        %-------------------------------------------------------------------------%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Calculate Baselines and Theoric Delay
        %-------------------------------------------------------------------------%
        % l is unknown; % latitude
        % L is unknown % Longitude
        
        X = @(l) D * (cosd(l)*sind(el) - sind(l)*cosd(el)*cosd(az));
        Y = D * (cosd(el)*sind(az));
        Z = @(l) D * (sind(l)*sind(el) + cosd(l)*cosd(el)*cosd(az));
        
        dX = @(l) D * (-sind(l)*sind(el) - cosd(l)*cosd(el)*cosd(az));
        dZ = @(l) D * (cosd(l)*sind(el) - sind(l)*cosd(el)*cosd(az));
        
        d2X = @(l) D * (-cosd(l)*sind(el) + sind(l)*cosd(el)*cosd(az));
        d2Z = @(l) D * (-sind(l)*sind(el) - cosd(l)*cosd(el)*cosd(az));
        
        dataRecordTime = [5,45,11,21,7,2016];
        dataRecordTime = num2cell(dataRecordTime);
        GST = calculate_GST(dataRecordTime{:});
        LST = @(L) mod(L/15+GST,24);
        
        dCygA = decCygnusA;
        rCygA = raCygnusA;
        % HCygA = @(L) LST(L) + h - rCygA;
        % dHCygA = 1/15;
        HCygA = @(L) (LST(L) + h - rCygA)*pi/12; % L is in degrees
        % dHCygA = 1;
        wCygA = @(l,L) cosd(dCygA)*cos(HCygA(L))*X(l) - cosd(dCygA)*sin(HCygA(L))*Y + sind(dCygA)*Z(l);
        dwCygA_l = @(l,L) cosd(dCygA)*cos(HCygA(L))*dX(l) + sind(dCygA)*dZ(l);
        dwCygA_L = @(l,L) -cosd(dCygA)*sin(HCygA(L))*X(l) - cosd(dCygA)*cos(HCygA(L))*Y;
        
        dwCygA_l_l = @(l,L) cosd(dCygA)*cos(HCygA(L))*d2X(l) + sind(dCygA)*d2Z(l);
        dwCygA_L_L = @(l,L) -cosd(dCygA)*cos(HCygA(L))*X(l) + cosd(dCygA)*sin(HCygA(L))*Y;
        
        dwCygA_l_L = @(l,L) -cosd(dCygA)*sin(HCygA(L))*dX(l);
        dwCygA_L_l = @(l,L) -cosd(dCygA)*sin(HCygA(L))*dX(l);
        
        %     dCasA = decCasA;
        %     rCasA = raCasA;
        %     % HCasA = @(L) LST(L) + h - rCasA;
        %     % dHCasA = 1/15;
        %     HCasA = @(L) L + h - rCasA;
        %     dHCasA = 1;
        %     wCasA = @(l,L) -cosd(dCasA)*sin(HCasA(L)*pi/12)*X + sind(dCasA)*Y(l) + cosd(dCasA)*cos(HCasA(L)*pi/12)*Z(l);
        %     dwCasA_l = @(l,L) sind(dCasA)*dY(l) + cosd(dCasA)*cos(HCasA(L)*pi/12)*dZ(l);
        %     dwCasA_L = @(l,L) -cosd(dCasA)*cos(HCasA(L)*pi/12)*dHCasA*pi/12*X - cosd(dCasA)*sin(HCasA(L)*pi/12)*dHCasA*pi/12*Z(l);
        %-------------------------------------------------------------------------%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Define Frequency
        %-------------------------------------------------------------------------%
        % fc = linspace(0,100,1024)*1e6;
        f = (33:0.1:37)*1e6; % 50 MHz
        % f = 50*1e6; % 50 MHz
        c = physconst('LightSpeed');
        %-------------------------------------------------------------------------%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Define Signal Model
        %-------------------------------------------------------------------------%
        modelledR = @(l,L) reshape(exp(1i*2*pi*wCygA(l,L)/c*f),[],1);
        
        %% First Dervatives
        modelledRdiff_l = @(l,L) reshape((1i*2*pi*dwCygA_l(l,L)/c*f) .* exp(1i*2*pi*wCygA(l,L)/c*f),[],1);
        modelledRdiff_L = @(l,L) reshape((1i*2*pi*dwCygA_L(l,L)/c*f) .* exp(1i*2*pi*wCygA(l,L)/c*f),[],1);
        
        %% Second Derivatives
        modelledRdiff_l_l = @(l,L) reshape((1i*2*pi*dwCygA_l_l(l,L)/c*f) .* exp(1i*2*pi*wCygA(l,L)/c*f) + ...
            (1i*2*pi*dwCygA_l(l,L)/c*f) .* (1i*2*pi*dwCygA_l(l,L)/c*f) .* exp(1i*2*pi*wCygA(l,L)/c*f),[],1);
        
        modelledRdiff_L_L = @(l,L) reshape((1i*2*pi*dwCygA_L_L(l,L)/c*f) .* exp(1i*2*pi*wCygA(l,L)/c*f) + ...
            (1i*2*pi*dwCygA_L(l,L)/c*f) .* (1i*2*pi*dwCygA_L(l,L)/c*f) .* exp(1i*2*pi*wCygA(l,L)/c*f),[],1);
        
        modelledRdiff_l_L = @(l,L) reshape((1i*2*pi*dwCygA_l_L(l,L)/c*f) .* exp(1i*2*pi*wCygA(l,L)/c*f) + ...
            (1i*2*pi*dwCygA_l(l,L)/c*f) .* (1i*2*pi*dwCygA_L(l,L)/c*f) .* exp(1i*2*pi*wCygA(l,L)/c*f),[],1);
        
        modelledRdiff_L_l = @(l,L) reshape((1i*2*pi*dwCygA_L_l(l,L)/c*f) .* exp(1i*2*pi*wCygA(l,L)/c*f) + ...
            (1i*2*pi*dwCygA_L(l,L)/c*f) .* (1i*2*pi*dwCygA_l(l,L)/c*f) .* exp(1i*2*pi*wCygA(l,L)/c*f),[],1);
        %-------------------------------------------------------------------------%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Define Measurement
        %-------------------------------------------------------------------------%
        lIn = 35.4; % Input latitude
        LIn = -116; % Input LST that is a function of Longitude
        signal = modelledR(lIn,LIn);
        
        SNRdB = 0:10:50; % dB
        distanceError = zeros(1,length(SNRdB));
        
        for nS = 1:length(SNRdB)
            SNR = 10.^(SNRdB(nS)/10);
            varNoise = var(signal)/SNR;
            noise = sqrt(varNoise/2)*(randn(length(h),length(f))+1i*randn(length(h),length(f)));
            noise = reshape(noise,[],1);
            
            % SIRdB = 30; % dB
            SIR = 10.^(SNRdB(nS)/10);
            varRFI = var(signal)/SIR;
            rfiBands = [4,7,9,16];
            fRFI = f(rfiBands);
            % fRFI = f;
            %     binsRFI = rand(length(h),1);
            %     binsRFI(binsRFI>=0.5) = 1;
            %     binsRFI(binsRFI<0.5) = 0;
            %     RFI = binsRFI * sqrt(varRFI)*exp(1i*2*pi*fRFI*100/c);
            
            RFI = zeros(length(h),length(f));
            for nF = 1:length(rfiBands)
                RFI(:,rfiBands(nF)) = varRFI*exp(1i*2*pi*fRFI(nF)*100/c);
            end
            % RFI = sqrt(varRFI)*exp(1i*2*pi*fRFI*100/c);
            RFI = reshape(RFI,[],1);
            
            measuredR = signal + noise + RFI;
            %-------------------------------------------------------------------------%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Define Min Function
            %-------------------------------------------------------------------------%
            minFunc = @(l,L) (measuredR - modelledR(l,L)) .* conj(measuredR - modelledR(l,L));
            
            %% First Dervatives
            JacobianMinFunc_l = @(l,L) -modelledRdiff_l(l,L) .* conj(measuredR - modelledR(l,L)) + ...
                (measuredR - modelledR(l,L)) .* conj(-modelledRdiff_l(l,L));
            
            JacobianMinFunc_L = @(l,L) -modelledRdiff_L(l,L) .* conj(measuredR - modelledR(l,L)) + ...
                (measuredR - modelledR(l,L)) .* conj(-modelledRdiff_L(l,L));
            
            JacobianMinFunc = @(l,L) [JacobianMinFunc_l(l,L), JacobianMinFunc_L(l,L)];
            
            %% Second Derivatives
            JacobianMinFunc_l_l = @(l,L) -modelledRdiff_l_l(l,L) .* conj(measuredR - modelledR(l,L)) + ...
                -modelledRdiff_l(l,L) .* conj(-modelledRdiff_l(l,L)) + ...
                -modelledRdiff_l(l,L) .* conj(-modelledRdiff_l(l,L)) + ...
                (measuredR - modelledR(l,L)) .* conj(-modelledRdiff_l_l(l,L));
            
            JacobianMinFunc_L_L = @(l,L) -modelledRdiff_L_L(l,L) .* conj(measuredR - modelledR(l,L)) + ...
                -modelledRdiff_L(l,L) .* conj(-modelledRdiff_L(l,L)) + ...
                -modelledRdiff_L(l,L) .* conj(- modelledRdiff_L(l,L)) + ...
                (measuredR - modelledR(l,L)) .* conj(-modelledRdiff_L_L(l,L));
            
            JacobianMinFunc_l_L = @(l,L) -modelledRdiff_l_L(l,L) .* conj(measuredR - modelledR(l,L)) + ...
                -modelledRdiff_l(l,L) .* conj(-modelledRdiff_L(l,L)) + ...
                -modelledRdiff_L(l,L) .* conj(-modelledRdiff_l(l,L)) + ...
                (measuredR - modelledR(l,L)) .* conj(-modelledRdiff_l_L(l,L));
            
            JacobianMinFunc_L_l = @(l,L) -modelledRdiff_L_l(l,L) .* conj(measuredR - modelledR(l,L)) + ...
                -modelledRdiff_L(l,L) .* conj(-modelledRdiff_l(l,L)) + ...
                -modelledRdiff_l(l,L) .* conj(- modelledRdiff_L(l,L)) + ...
                (measuredR - modelledR(l,L)) .* conj(-modelledRdiff_L_l(l,L));
            
            JacobianMinFunc2 = @(l,L) [JacobianMinFunc_l_l(l,L), JacobianMinFunc_l_L(l,L);...
                JacobianMinFunc_L_l(l,L), JacobianMinFunc_L_L(l,L)];
            %-------------------------------------------------------------------------%
            
            lSearch = linspace(0,90,100);
            LSearch = linspace(-180,180,100);
            % lSearch = lIn;
            % LSearch = LIn;
            mR = zeros(length(lSearch),length(LSearch));
            for nl = 1:length(lSearch)
                for nL = 1:length(LSearch)
                    lS = lSearch(nl);
                    LS = LSearch(nL);
                    MF = modelledR(lS,LS)';
                    mR(nl,nL) = MF * measuredR(:);
                end
            end
            % figure;
            % mesh(LSearch,lSearch,abs(mR));
            % plot(lSearch,abs(mR));
            % plot(LSearch*12/pi,abs(mR));
            
            [~,k] = max(abs(mR(:)));
            [k_l, k_L] = ind2sub(size(mR),k);
            
            lEst = lSearch(k_l);
            lEst1 = lSearch(k_l-1);
            lEst2 = lSearch(k_l+1);
            
            LEst = LSearch(k_L);
            LEst1 = LSearch(k_L-1);
            LEst2 = LSearch(k_L+1);
            
            %% Newton and Secant Search for Two Parameter
            [lL_k,~,all_lL] = newton(minFunc,JacobianMinFunc,[lEst;LEst]);
            % [lL_k,~,all_lL] = newton(JacobianMinFunc,JacobianMinFunc2,[lEst;LEst]);
            
            % [lL_k,~,all_lL] = secant(JacobianMinFunc,[lEst1;LEst1],[lEst2;LEst2]);
            % [lL_k,~,all_lL] = secant(minFunc,[lEst1;LEst1],[lEst2;LEst2]);
            
            % figure;
            % plot(abs(lIn-all_lL(1,:)))
            % figure;
            % plot(all_lL(2,:))
            %         figure; plot(abs(lIn-all_lL(1,:)))
            %         figure; plot(abs(LIn-all_lL(2,:)))
            
            wgs84 = wgs84Ellipsoid('meters');
            error = distance(all_lL(1,:),all_lL(2,:),lIn,LIn,wgs84);
            distanceError(nS) = error(end);
            % figure; plot(error)
        end
        plot(SNRdB,distanceError)
    end
end


%% Secant Method
function [x_k,k,all_x] = secant(f,x_0,x_1)
% f = nonlinear function
% x_0 = first initial guess: row vector
% x_1 = second initial guess: row vector

% x = root of the function

% tol = 10^-5.5;   % tolerance value
nIter = 1e3; % iteration value

x_k_2 = x_0;    % x(k-2)
x_k_1 = x_1;    % x(k-1)

all_x = zeros(length(x_0),nIter);
k = 1;
while k <= nIter % abs(f(x_k)) > tol   % termination criterion
    
    x_k_2_c = num2cell(x_k_2);
    x_k_1_c = num2cell(x_k_1);
    
    f_x_k_2 = f(x_k_2_c{:});
    f_x_k_1 = f(x_k_1_c{:});
    J = f_x_k_1 - f_x_k_2;
    
    x_k = x_k_1 - ((J' * J) \ J') * f_x_k_1 * (x_k_1 - x_k_2);
    
    if length(x_0)>1
        all_x(:,k) = x_k;
    else
        all_x(k) = x_k;
    end
    
    x_k_2 = x_k_1;
    x_k_1 = x_k;
    
    k = k + 1;
end
end

%% Newton Method
function [x_k,k,all_x] = newton(f,df,x_0)
% f = nonlinear function
% df = derivative of the nonlinear function
% x_0 = initial guess: row vector

% x = root of the function
% k = number of iterations
% all_x = all roots at each iteraion

% tol = 10^-4;   % tolerance value
nIter = 1e3; % iteration value

x_k_1 = x_0;    % x(k-1)

all_x = zeros(length(x_0),nIter);
k = 1;
while k <= nIter % norm(x_k - x_k_1) > tol % termination criterion
    
    x_k_1_c = num2cell(x_k_1);
    J = df(x_k_1_c{:});
    f_k_1 = f(x_k_1_c{:});
    
    x_k = x_k_1 - 1*((J' * J) \ J') * f_k_1;
    
    if length(x_0)>1
        all_x(:,k) = x_k;
    else
        all_x(k) = x_k;
    end
    
    x_k_1 = x_k;
    
    k = k + 1;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the GST
%-------------------------------------------------------------------------%
function GST = calculate_GST(hour,minute,second,day,month,year)

% JD2000 = juliandate(2000,1,1,12,0,0);
JD2000 = juliandate(2000,1,1,11,58,55.816);

JD = juliandate(year,month,day);
UT = hour + minute/60 + second/3600;
T = (JD - JD2000)/36525.0;
GST = 6.697374558 + (2400.051336*T)+(0.000025862*T^2)+(UT*1.0027379093);
GST = mod(GST,24);
end
%-------------------------------------------------------------------------%