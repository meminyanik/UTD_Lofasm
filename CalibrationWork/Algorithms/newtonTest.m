%% LoFASM Test
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Calculate Total Scanned Hour Duration
% %-------------------------------------------------------------------------%
% h = (0 : 447*0.083886 : (504-1)*447*0.083886)/3600;
% %-------------------------------------------------------------------------%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Calculate Baselines and Theoric Delay
% %-------------------------------------------------------------------------%
% l = 35.247469;
% L = -116.791509;
% 
% X = 20;
% Y = 160;
% Z = 3;
% 
% % d = 40;
% r = 20;
% H = h - r;
% 
% w = @(d) -cosd(d)*sin(H*pi/12)*X + sind(d)*Y + cosd(d)*cos(H*pi/12)*Z;
% 
% dw = @(d) sind(d)*sin(H*pi/12)*X + cosd(d)*Y - sind(d)*cos(H*pi/12)*Z;
% %-------------------------------------------------------------------------%
% 
% fc = (0.09765625*100:0.09765625:0.09765625*100)*1e6;
% c = physconst('LightSpeed');
% 
% measuredR = reshape(exp(1i*2*pi*fc.'*(w(30)/c)),[],1);
% 
% modelledR = @(d) reshape(exp(1i*2*pi*fc.'*(w(d)/c)),[],1);
% 
% modelledRdiff = @(d) reshape(1i*2*pi*fc.'/c.*dw(d).*exp(1i*2*pi*fc.'*(w(d)/c)),[],1);
% 
% minFunc = @(d) (measuredR - modelledR(d)) .* conj(measuredR - modelledR(d));
% 
% JacobianMinFunc = @(d) (-modelledRdiff(d)) .* conj(measuredR - modelledR(d)) + ...
%     (measuredR - modelledR(d)) .* conj(-modelledRdiff(d));
% 
% [d_newton, k_newton, error] = newton(minFunc,JacobianMinFunc,20);


%% Simple Test 1 - Two Variable
% t = linspace(0,5*pi,100).';
% 
% modelledF = @(a,b) a^2*exp(1i*t) + b*cos(t).^2;
% modelledFdiff_a = @(a) 2*a*exp(1i*t);
% modelledFdiff_b = @(b) cos(t).^2;
% 
% % modelledFsecondDiff_a = @(a) 2*exp(1i*t);
% % modelledFsecondDiff_b = @(b) 0;
% 
% aDes = 135;
% bDes = 100;
% measuredF = aDes^2*exp(1i*t) + bDes*cos(t).^2 + 0.1*randn(length(t),1);
% 
% minFunc = @(a,b) (measuredF - modelledF(a,b)) .* conj(measuredF - modelledF(a,b));
% 
% JacobianMinFunc = @(a,b) [(-modelledFdiff_a(a) .* conj(measuredF - modelledF(a,b))) + ...
%     ((measuredF - modelledF(a,b)) .* conj(- modelledFdiff_a(a))) , ...
%     (-modelledFdiff_b(b) .* conj(measuredF - modelledF(a,b))) + ...
%     ((measuredF - modelledF(a,b)) .* conj(- modelledFdiff_b(b)))];
% 
% aT = linspace(10,300,150);
% bT = linspace(10,200,100);
% mR = zeros(length(aT),length(bT)); 
% for n1 = 1:length(aT)
%     for n2 = 1:length(bT)
%         x_k_1 = [aT(n1);bT(n2)];
%         x_k_1_c = num2cell(x_k_1);
%         J = JacobianMinFunc(x_k_1_c{:});
%         fm = minFunc(x_k_1_c{:});
%         mR(n1,n2) = norm(-(J.' * J) \ J.' * fm);
%     end
% end
% mesh(bT,aT,mag2db(mR))


%% Simple Test 2 - Two Variable
t = linspace(0,5*pi,100).';
modelledF = @(a,b) exp(1i*a*t) + b*cos(t).^2;
modelledFdiff_a = @(a) 1i*t.*exp(1i*a*t);
modelledFdiff_b = @(b) cos(t).^2;

aDes = 15;
bIn = 80;
measuredF = exp(1i*aDes*t) + bIn*cos(t).^2 + 0.1*randn(length(t),1);

minFunc = @(a,b) (measuredF - modelledF(a,b)) .* conj(measuredF - modelledF(a,b));

JacobianMinFunc = @(a,b) [(-modelledFdiff_a(a) .* conj(measuredF - modelledF(a,b))) + ...
    ((measuredF - modelledF(a,b)) .* conj(- modelledFdiff_a(a))) , ...
    (-modelledFdiff_b(b) .* conj(measuredF - modelledF(a,b))) + ...
    ((measuredF - modelledF(a,b)) .* conj(- modelledFdiff_b(b)))];

aT = linspace(13,16,150);
bT = linspace(10,200,100);
mR = zeros(length(aT),length(bT)); 
for n1 = 1:length(aT)
    for n2 = 1:length(bT)
        x_k_1 = [aT(n1);bT(n2)];
        x_k_1_c = num2cell(x_k_1);
        J = JacobianMinFunc(x_k_1_c{:});
        fm = minFunc(x_k_1_c{:});
        mR(n1,n2) = norm(-(J' * J) \ J' * fm);
        % mR(n1,n2) = norm(-real((J' * J) \ J' * fm));
    end
end
mesh(bT,aT,mag2db(mR))


%% Simple Test 3 - One Variable
% t = linspace(0,5*pi,100).';
% 
% % V1 -------------------------------------------------------------------------%
% % modelledR = @(a) exp(1i*a*t);
% % modelledRdiff = @(a) 1i*t.*exp(1i*a*t);
% % 
% % measuredR = modelledR(3) + 0.1*(randn(length(t),1)+1i*randn(length(t),1));
% % 
% % minFunc = @(a) (measuredR - modelledR(a)) .* conj(measuredR - modelledR(a));
% % 
% % JacobianMinFunc = @(a) (-modelledRdiff(a) .* conj(measuredR - modelledR(a))) + ...
% %     (measuredR - modelledR(a)) .* conj(-modelledRdiff(a));
% %-------------------------------------------------------------------------%
% 
% 
% % V2 -------------------------------------------------------------------------%
% modelledR = @(a) exp(1i*a*t);
% modelledRdiff = @(a) 1i*t.*exp(1i*a*t);
% 
% measuredR = modelledR(3.5) + 0.1*(randn(length(t),1)+1i*randn(length(t),1));
% 
% % minFunc = @(a) abs(measuredR - modelledR(a)).^2;
% 
% minFunc = @(a) (measuredR - modelledR(a));
% 
% % JacobianMinFunc = @(a) (-modelledRdiff(a) .* conj(measuredR - modelledR(a))) + ...
% %     (measuredR - modelledR(a)) .* conj(-modelledRdiff(a));
% 
% JacobianMinFunc = @(a) -modelledRdiff(a)/length(t); %   - conj(modelledRdiff(a));
%-------------------------------------------------------------------------%


% V3 -------------------------------------------------------------------------%
% realModelledR = @(a) cos(a*t);
% realModelledRdiff = @(a) -t.*sin(a*t);
% 
% imagModelledR = @(a) sin(a*t);
% imagModelledRdiff = @(a) t.*cos(a*t);
% 
% realMeasuredR = realModelledR(3) + 0.1*randn(length(t),1);
% imagMeasuredR = imagModelledR(3) + 0.1*randn(length(t),1);
% 
% minFunc = @(a) (realMeasuredR - realModelledR(a)).^2 + (imagMeasuredR - imagModelledR(a)).^2;
% JacobianMinFunc = @(a) -2*(realMeasuredR - realModelledR(a)).*realModelledRdiff(a) + ...
%     -2*(imagMeasuredR - imagModelledR(a)).*imagModelledRdiff(a);
%-------------------------------------------------------------------------%

% aT = linspace(0.01,5,100);
% mR = zeros(1,length(aT));
% for n = 1:length(aT)    
%     x_k_1 = aT(n);
%     J = JacobianMinFunc(x_k_1);
%     fm = minFunc(x_k_1);
%     mR(n) = - (J'*fm - J.' * conj(fm));
% end
% plot(aT,mR)
% [~,ind] = min(mR);
% aEstimated = aT(ind);

%% Find the minimum
% [a_newton, k_newton, allA, error] = newton(minFunc,JacobianMinFunc,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Newton Method
%-------------------------------------------------------------------------%
function [x_k,k,allx,error] = newton(f,df,x_0)
% f = nonlinear function
% df = derivative of the nonlinear function
% x_0 = initial guess

% x = root of the function
% k = number of iterations

% If root error is wanted to be minimized
% tol = 10^-5;   % tolerance value
nIter = 1e4; % iteration value

x_k_1 = x_0;

% x_k_1_c = num2cell(x_k_1);
% J = df(x_k_1_c{:});

J = df(x_k_1);
fm = f(x_k_1);
x_k = x_k_1 - (J.' * J) \ J.' * fm;
k = 1;

% If root error is wanted to be minimized
error = zeros(1,nIter);
allx = zeros(1,nIter);

while k < nIter % norm(x_k - x_k_1) > tol % termination criterion
    x_k_1 = x_k;
    
    % x_k_1_c = num2cell(x_k_1);
    % J = df(x_k_1_c{:});
    
    J = df(x_k_1);
    fm = f(x_k_1);
    
    x_k = x_k_1 - (J.' * J) \ J.' * fm;
    
    error(k) = norm(f(x_k));
    allx(k) = x_k;
    k = k + 1;
end

end
%-------------------------------------------------------------------------%