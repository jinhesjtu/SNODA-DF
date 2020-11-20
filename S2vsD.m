%
clear
close all
%    Performance vs Sensor Spacing
count = 0;
TrialAll = 500;
Ms = 3;    % Inner sensors No.
Ns = 3;     % Outer sensors No.
Sensor = Ms + Ns;
Snap = 1000;        Lambda = 1;
Doa = [40 60];
SigAll = length(Doa);
SnrdB = [20];
for  sd = [1 2 4 8 16 25 32 49 64 81 100 200  400 1000]
    disp(sprintf('sd: %.1f', sd));
    count = count + 1;
    d = sd*Lambda;    D = sd*Lambda;
    Sp1 = [0:Ms-1]*d;    % inner sensor location
    Sp2 = [Ms*d:(Ms+1)*d:(Ns*(Ms+1)-1)*d];  % outer sensor location
    Spd = [Sp1,Sp2];
    Ad = exp(-j*2*pi/Lambda*Spd'*sind(Doa)); Al = Ad*diag(exp(j*2*pi*D*cosd(Doa)));
    SnrdB = 20;
    Snr=sqrt(10.^(SnrdB/10));        f = [0.2,0.4];
    for num = 1:SigAll
        Sig(num,:) = fmlin(Snap,f(num), f(num));
    end
    Sigh =  diag(sind(Doa)) *Sig;
    Sigv =  diag(cosd(Doa)) *Sig;
    Xdp = Ad*Sigh;
    Xlp =  Al*Sigv;
    X0 =  [Xdp;Xlp];
    nsucy = [0,0]; nsucz = [0,0]; nsucave = [0,0];
    Tfiny = zeros(TrialAll,SigAll); Tfinz = zeros(TrialAll,SigAll); Tfinave = zeros(TrialAll,SigAll);
    for trial = 1:TrialAll;
        X = awgn(X0, SnrdB,'measured');
        %%%%%%%%%%%%%%%%%%%%%%%%
        Xy = X(1:Sensor,:);  Xz = X(Sensor + 1:2 * Sensor,:);  Ryy = Xy*Xy'; Ryz = Xy*Xz'; Rzz = Xz*Xz'; Rzy = Xz*Xy';
        Zyy = NestR(Ryy,Ms,Ns).';  Zyz = NestR(Ryz,Ms,Ns).';  Zzz = NestR(Rzz,Ms,Ns).';  Zzy = NestR(Rzy,Ms,Ns).';
        L = (length(Zyy) + 1)/2;
        for num = 1:L
            bRyy(:,num) = Zyy(num:length(Zyy)-L+num);   bRyz(:,num) = Zyz(num:length(Zyz)-L+num);
            bRzz(:,num) = Zzz(num:length(Zzz)-L+num);   bRzy(:,num) = Zzy(num:length(Zzy)-L+num);
        end
        bRy = [bRyy;bRyz]; bRz = [bRzz;bRzy];
        [U S V] = svd(bRy); Esy = U(:,1:SigAll);   [U S V] = svd(bRz); Esz = U(:,1:SigAll);
        
        %%%  %%%%%%%%%%%%%%%%%%%%%    initial angle
        Ey1 = Esy(1:L,:);  Ey2 = Esy(L+1:2*L,:);   Ez1 = Esz(1:L,:);  Ez2 = Esz(L+1:2*L,:);
        [Tiy, Phiy] = eig(pinv(Ey1)*Ey2);  [Tiz, Phiz] = eig(pinv(Ez1)*Ez2);
        tiniy = acotd(abs((diag(Phiy))));  tiniz = atand(abs((diag(Phiz))));
        if abs(tiniy(1)) > abs(tiniy(2))
            Phiy =diag(rot90(eye(2),3)*diag(Phiy));
            tiniy = rot90(eye(2),3)*tiniy;
        end
        if abs(tiniz(1)) > abs(tiniz(2))
            Phiz =diag(rot90(eye(2),3)*diag(Phiz));
            tiniz = rot90(eye(2),3)*tiniz;
        end
        tiniyb = -acotd(abs((diag(Phiy)))); tinizb = -atand(abs((diag(Phiz))));
        %%%  %%%%%%%%%%%%%%%%%%%%%    initial angle
        
        
        nx = -2*D/Lambda:2*D/Lambda;
        %%%%%%%%%%%%   Fine angle
        for ns= 1:SigAll
            py = (angle(Phiy(ns,ns))+ 2*pi*nx)/(-2*pi*D);
            ind =  (py >=0) & (py <=1);
            phiy = py(find(ind == 1)); tally = acosd(phiy);
            thaty(ns) = tally(find(abs(tally - tiniy(ns)) == min(abs(tally - tiniy(ns)))));
            Refy(ns) = min(abs(tally - tiniy(ns)));
            pyb = (angle(Phiy(ns,ns)) + pi + 2*pi*nx)/(-2*pi*D);
            ind =  (pyb >=0) & (pyb <=1);
            phiyb = pyb(find(ind == 1)); tallyb = -acosd(phiyb);
            thatyb(ns) = tallyb(find(abs(tallyb - tiniyb(ns)) == min(abs(tallyb - tiniyb(ns)))));
            Refyb(ns) = min(abs(tallyb - tiniyb(ns)));
            Thaty(ns) = (Refy(ns) < Refyb(ns))*thaty(ns) + (Refy(ns) > Refyb(ns))*thatyb(ns);
            
            
            pz = (angle(Phiz(ns,ns))+ 2*pi*nx)/(2*pi*D);
            ind =  (pz >=0) & (pz <=1);
            phiz = pz(find(ind == 1)); tallz = acosd(phiz);
            thatz(ns) = tallz(find(abs(tallz - tiniz(ns)) == min(abs(tallz - tiniz(ns)))));
            Refz(ns) = min(abs(tallz - tiniz(ns)));
            pzb = (angle(Phiz(ns,ns)) + pi + 2*pi*nx)/(2*pi*D);
            ind =  (pzb >=0) & (pzb <=1);
            phizb = pzb(find(ind == 1)); tallzb = -acosd(phizb);
            thatzb(ns) = tallzb(find(abs(tallzb - tinizb(ns)) == min(abs(tallzb - tinizb(ns)))));
            Refzb(ns) = min(abs(tallzb - tinizb(ns)));
            Thatz(ns) = (Refz(ns) < Refzb(ns))*thatz(ns) + (Refz(ns) > Refzb(ns))*thatzb(ns);
            
            %%%% success rate for resolving +- 1 ambiguty
        ay = (sign(Thaty(ns)) == sign(Doa(ns))); az = (sign(Thatz(ns)) == sign(Doa(ns)));    a1= (sign(Thatz(ns)) == sign(Thaty(ns)));
        a2 = (a1 == a1.*ay.*az);
        
        if prod(ay) == 1
            nsucy(ns) = nsucy(ns) + 1;
            Tfiny(nsucy(ns),ns) = Thaty(ns);
        end
        if prod(az) == 1
            nsucz(ns) = nsucz(ns) + 1;
            Tfinz(nsucz(ns),ns) = Thatz(ns);
        end
        if prod(a1) * prod(ay) == 1 || sum(a1) == 0 ||  prod(a2) == 1
            nsucave(ns) = nsucave(ns) + 1;
        end
        end       
        %%%% success rate for resolving +- 1 ambiguty

        Trefy(trial,:) = sort(tiniy);   Trefz(trial,:) = sort(tiniz);
        %%%% success rate for resolving +- 1 ambiguty
    end
    sry(:,count) = nsucy'/TrialAll;  srz(:,count) = nsucz'/TrialAll;
    srave(:,count) = nsucave'/TrialAll;   
    
    %%%%%%%%%%%%  Fine angle
    Tref = (Trefy + Trefz)/2;  
     erefy = rmse(Trefy, Doa); erefz = rmse(Trefz, Doa);
     eref(:,count) = rmse(Tref, Doa);    
     
     for ns = 1: SigAll
         efiny(ns,count) = rmse(Tfiny(find(Tfiny(:,ns) ~= 0),ns), Doa(ns));
         efinz(ns,count) = rmse(Tfinz(find(Tfinz(:,ns) ~= 0),ns), Doa(ns));
         len = min(length(Tfiny(find(Tfiny(:,ns) ~= 0),ns)), length(Tfinz(find(Tfinz(:,ns) ~= 0),ns)));
         Tfine = 0.5*(Tfiny(1:len,ns) + Tfinz(1:len,ns));
         efine(ns,count) = rmse(Tfine, Doa(ns));
     end


end

ErrS2prop = [sry;srz;srave;eref;efiny;efinz;efine];

save S2prop ErrS2prop









