clear
close all
%% Coupling leakage of different array geometries.
count = 0;
for dd = [0.5:0.5:100];
    count = count + 1;    
    cc = 0.5*exp(j*pi/3);
    %%%%%%%%%%%%%%%%%%%%%%%% ULA L = 6
    L = 6;
    d = 0.5;
    Sp = [0:L-1]*d;
    for ni = 1:L
        for nj = 1:L
            ds = abs(Sp(ni) - Sp(nj))/d;
            if ds == 0
                Cula(ni,nj) = 1;
            elseif ds == 1
                Cula(ni,nj) = cc;
            else
                Cula(ni,nj) = cc/ds * exp(j*(ds-1)/8);
            end
        end
    end
    Lula(count) = norm(Cula- diag(diag(Cula)),'fro')/norm(Cula,'fro');    
     %%%%%%%%%%%%%%%%%%%%%%%% ULA L = 6 
     
     
         %%%%%%%%%%%%%%%%%%%%%%%% ULA L = 12
    L = 12;
    d = 0.5;
    Sp = [0:L-1]*d;
    for ni = 1:L
        for nj = 1:L
            ds = abs(Sp(ni) - Sp(nj))/d;
            if ds == 0
                Culab(ni,nj) = 1;
            elseif ds == 1
                Culab(ni,nj) = cc;
            else
                Culab(ni,nj) = cc/ds * exp(j*(ds-1)/8);
            end
        end
    end
    Lulab(count) = norm(Culab- diag(diag(Culab)),'fro')/norm(Culab,'fro');    
     %%%%%%%%%%%%%%%%%%%%%%%% ULA L = 12 
    
     %%%%%%%%%%%%%%%%%%%%%%%%  NA L = 6     
     Ms = 3;    % Inner sensors No.
     Ns = 3;     % Outer sensors No.     
     Sp1 = [0:Ms-1]*d;    % inner sensor location
     Sp2 = [Ms*d:(Ms+1)*d:(Ns*(Ms+1)-1)*d];  % outer sensor location
     Sp = [Sp1,Sp2];    
     for ni = 1:length(Sp)
         for nj = 1:length(Sp)
             ds = abs(Sp(ni) - Sp(nj))/d;
             if ds == 0
                 Cna(ni,nj) = 1;
             elseif ds == 1
                 Cna(ni,nj) = cc;
             else
                 Cna(ni,nj) = cc/ds * exp(j*(ds-1)/8);
             end
         end
     end     
    Lna(count) = norm(Cna - diag(diag(Cna)),'fro')/norm(Cna,'fro');
   %%%%%%%%%%%%%%%%%%%%%%%%  NA L = 6  
   
   
   
        %%%%%%%%%%%%%%%%%%%%%%%%  NA L = 12     
     Ms = 6;    % Inner sensors No.
     Ns = 6;     % Outer sensors No.     
     Sp1 = [0:Ms-1]*d;    % inner sensor location
     Sp2 = [Ms*d:(Ms+1)*d:(Ns*(Ms+1)-1)*d];  % outer sensor location
     Sp = [Sp1,Sp2];    
     for ni = 1:length(Sp)
         for nj = 1:length(Sp)
             ds = abs(Sp(ni) - Sp(nj))/d;
             if ds == 0
                 Cnab(ni,nj) = 1;
             elseif ds == 1
                 Cnab(ni,nj) = cc;
             else
                 Cnab(ni,nj) = cc/ds * exp(j*(ds-1)/8);
             end
         end
     end     
    Lnab(count) = norm(Cnab - diag(diag(Cnab)),'fro')/norm(Cnab,'fro');
   %%%%%%%%%%%%%%%%%%%%%%%%  NA L = 12  

   %%%%%%%%%%%%%%%%%%%%%%%%  SNA
   Ms = 3;    % Inner sensors No.
   Ns = 3;     % Outer sensors No.
   Sensor = Ms + Ns;
   DD = dd;   
   Sp1 = [0:Ms-1]*dd;    % inner sensor location
   Sp2 = [Ms*dd:(Ms+1)*dd:(Ns*(Ms+1)-1)*dd];  % outer sensor location
   Spd = [Sp1,Sp2];  Spl = Spd + DD;   
   Spa = [zeros(Sensor,1), Spd'];  Spb = [Spl', Spd'];
   Sp = [Spa;Spb];   
   for ni = 1:(2*Sensor)
       for nj = 1:(2*Sensor)
           ds = norm(Sp(ni,:) - Sp(nj,:))/d;
           if ds == 0
               Csna(ni,nj) = 1;
           elseif ds == 1
               Csna(ni,nj) = cc;
           else
               Csna(ni,nj) = cc/ds * exp(j*(ds-1)/8);
           end
       end
   end
   
   Lsna(count) = norm(Csna - diag(diag(Csna)),'fro')/norm(Csna,'fro');
end

%%%%%%%%%%%%%%%%%%%%%%%%  SNA 

Dp = [0.5:0.5:100];
figure;  set(gcf,'DefaultLineLineWidth',1.5)
semilogx(Dp, Lula, '-', Dp, Lulab, '-', Dp, Lna, '-', Dp, Lnab, '-', Dp, Lsna, '-k'); grid on;
xlabel('Sensor/element Spacing, in wavelengths')
ylabel('Coupling leakage')
title('Coupling leakage for SNODA with different sensor/element Spacing')
legend('ULA, 6 sensors', 'ULA, 12 sensors', 'Nested array, 6 sensors', 'Nested array, 12 sensors', 'SNODA')



            
