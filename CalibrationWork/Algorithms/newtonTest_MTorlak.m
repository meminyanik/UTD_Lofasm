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
% measuredF = modelledF(185,120) + 0.1*randn(1,length(t)).';
% 
% minFunc = @(a,b) (measuredF - modelledF(a,b)) .* conj(measuredF - modelledF(a,b));
% 
% JacobianMinFunc = @(a,b) [(-modelledFdiff_a(a) .* conj(measuredF - modelledF(a,b))) + ...
%     ((measuredF - modelledF(a,b)) .* conj(- modelledFdiff_a(a))) , ...
%     (-modelledFdiff_b(b) .* conj(measuredF - modelledF(a,b))) + ...
%     ((measuredF - modelledF(a,b)) .* conj(- modelledFdiff_b(b)))];
% 
% a1 = linspace(10,300,150);
% a2 = linspace(10,200,100);
% mR = zeros(length(a1),length(a2)); 
% for n1 = 1:length(a1)
%     for n2 = 1:length(a2)
%         x_k_1 = [a1(n1);a2(n2)];
%         x_k_1_c = num2cell(x_k_1);
%         J = JacobianMinFunc(x_k_1_c{:});
%         mR(n1,n2) = norm(-(J.' * J) \ J.' * minFunc(x_k_1_c{:}));
%     end
% end
% mesh(a2,a1,mag2db(mR))


%% Simple Test 2 - One Variable
t = linspace(0,50*pi,1000).';

modelledF = @(a,k) a*exp(1i*k*t);
modelledFdiff_a = @(a,k) exp(1i*k*t);
modelledFdiff_k = @(a,k) a*1i*t.*exp(1i*k*t);

% modelledF = @(a) cos(a*t) + 1i*sin(a*t);
% modelledFdiff = @(a) -t.*sin(a*t) + 1i*t.*cos(a*t);

measuredF = 5*exp(1i*2*t) + 0.1*(randn(length(t),1) + 1i*randn(length(t),1));

measuredFfft = abs(fft(measuredF,50));


% minFunc = @(a)( (measuredF - modelledF(a)) .* conj(measuredF - modelledF(a))/length(t));
minFunc = @(a,k)(measuredF - modelledF(a,k));
% JacobianMinFunc = @(a) (-modelledFdiff(a).' * conj(measuredF - modelledF(a)))/length(t) + ...
%     (measuredF - modelledF(a)).' * conj(-modelledFdiff(a))/length(t);
JacobianMinFunc = @(a,k) -1*[modelledFdiff_a(a,k) modelledFdiff_k(a,k)]; 

kT = linspace(0.01,5,100);
% aT = linspace(0.01,10,100);
aT = 5;
mF = zeros(length(aT),length(kT));
for nK = 1:length(kT)
    for nA = 1:length(aT)
        J = JacobianMinFunc(aT(nA),kT(nK));
        fm = minFunc(aT(nA),kT(nK));
        % mF(nA,nK) = norm(-imag((J' * J) \ J' * fm));
        TT = -imag((J' * J) \ J' * fm);
        mF(nA,nK) = TT(2);
    end
end
% mesh(aT,kT,mF)
plot(kT,mF)

%% Find the minimum
% [a_newton, k_newton, all_a, error] = newton(minFunc,JacobianMinFunc,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Newton Method
%-------------------------------------------------------------------------%
function [x_k,k,all_x,error] = newton(f,df,x_0)
% f = nonlinear function
% df = derivative of the nonlinear function
% x_0 = initial guess

% x = root of the function
% k = number of iterations

% If root error is wanted to be minimized
% tol = 10^-4;   % tolerance value
nIter = 1e4; % iteration value

x_k_1 = x_0;
%x_k_1_c = num2cell(x_k_1);
%J = df(x_k_1_c{:});
J = df(x_k_1);
fm = f(x_k_1);
x_k = x_k_1 - imag((J' * J) \ J' * fm);
k = 0;

% If root error is wanted to be minimized
error = zeros(1,nIter);
all_x = zeros(1,nIter);
while k < nIter % norm(x_k - x_k_1) > tol % termination criterion
    k = k + 1;
    error(k) = norm(f(x_k));
    all_x(k) = x_k;
    
    x_k_1 = x_k;
    J = df(x_k_1);
    fm = f(x_k_1);
    x_k = x_k_1 - imag((J' * J) \ J' * fm);
    
end

end
%-------------------------------------------------------------------------%