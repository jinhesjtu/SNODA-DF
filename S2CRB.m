clear all

syms theta

count = 0;
for sd =[1 2 4 8 16 25 32 49 64 81 100 200 400 1000]
    count = count + 1;
    disp(sprintf('sd: %.1f', sd));
    SnrdB = 20;
    Ms = 3;    % Inner sensors No.
    Ns = 3;     % Outer sensors No.
    Sensor = Ms + Ns;
    Snap = 1000;
    Lambda = 1;
    d = sd*Lambda;    D = sd*Lambda;
    Theta = pi/180*([40, 60]);
    SigAll = length(Theta);
    Sp1 = [0:Ms-1]*d;    % inner sensor location
    Sp2 = [Ms*d:(Ms+1)*d:(Ns*(Ms+1)-1)*d];  % outer sensor location
    Spd = [Sp1,Sp2];  Spl = Spd + D;
    Snr=sqrt(10.^(SnrdB/10));
    f = [0.2,0.4];
    for num = 1:SigAll
        Sig(num,:) = fmlin(Snap,f(num), f(num));
    end    
    
    
    b = [sin(theta); cos(theta)*exp(-j*2*pi*D*sin(theta))];
    q = exp(-j*2*pi/Lambda*Spd'*sin(theta));
    
    A = kron(b,q);
    
    t1 = Theta(1); t2 = Theta(2);
    p_theta1 = subs(diff(A,'theta'), {theta},{t1});
    p_theta2 = subs(diff(A,'theta'), {theta},{t2});
    
    J = [p_theta1'*p_theta1, p_theta1'*p_theta2 * Sig(1,:)*Sig(2,:)';
        p_theta2'*p_theta1 * Sig(2,:)*Sig(1,:)', p_theta2'*p_theta2 ];
    
    JJ = double(J);
    JJ = real(JJ);
    JJ = inv(JJ);
    
    Snr=sqrt(10.^(SnrdB/10));
    crb(:,count) = 180/pi*(sqrt(diag(JJ)/(2*Snap*(Snr^2))));
end
ErrS2Crb = crb;

save S2Crb ErrS2Crb





