clear all;
close all;
clc;
%////////////////////////////////////////
%/////////////////erwthmata A////////////
%////////////////////////////////////////
T = 0.001;
over =10;
Ts = T/over;
A = 3;
a = 0.5;
Nf = 2048;
%Creation SRRC pulses
%/////////////////////A.1///////////////////////
Fs = 1/Ts;
[srrc, t] = srrc_pulse(T,Ts,A,a);
F = abs(fftshift(fft(srrc, Nf))*Ts).^2;
%t_F2 = -Fs/2:Fs/Nf:Fs/2 - Fs/Nf
t_F=linspace(-Fs/2,Fs/2,Nf);%Xwrizoume to Faxis se Φ f pou isapexoun
figure();
semilogy(t_F,F)
title('Φασματική πυκνότητα ενέργειας |Φ(F)|^2')
xlabel('συχνότητα(Hz)')
ylabel('|Φ(F)|^2')


%//////////////////A.2////////////////////
N_bits = 50;
%timeN = [0:Ts:N*T];

%Creation of series of symbols
b = (sign(randn(N_bits,1))+1)/2;
%Decoding to 2-PAM
Xn = bits_to_2PAM(b);
sum_upsample=(1/Ts)*upsample(Xn,over);
t_sum_upsample=0:Ts:N_bits*T-Ts;
X=conv(sum_upsample,srrc)*Ts;
t_X=t(1)+t_sum_upsample(1):Ts:t(end)+t_sum_upsample(end);
figure();
plot(t_X,X)
title('X(t),σε 2PAM')
xlabel('χρονος(s)')
ylabel('Χ(t)')

%var(X)επιστρέφει τη διακύμανση(διασπορα) των στοιχείων Aκατά μήκος της πρώτης διάταξης 
%του πίνακα των οποίων το μέγεθος δεν ισούται με 1.

% v=(var(X));
% Sx=(F/T).*v ;

%%%%%%%%
var1= 1/2*(1)^2  +  1/2*(-1)^2 %var=e(x^2)-[e(x)]^2=e(x^2)=pi8anotita tou 1*  to 1^2 + pi8anotita tou -1*  to (-1)^2 
Sx=(F/T).*var1 ;
%%%%%%%%%%

%////////////////////////A.3///////////////////
%//////////1o meros


t_Px=linspace(-Fs/2,Fs/2,Nf);
F_X_abs=abs(fftshift(fft(X,length(t_Px)))*Ts).^2;
Ttotal=length(X).*Ts;
Px=F_X_abs/Ttotal;

figure;
semilogy(t_Px,Px);
title('Px(F)με semilogy')
xlabel('συχνοτητα(Hz)')
ylabel('Px')

figure;
plot(t_Px,Px);
title('Px(F)με plot')
xlabel('συχνοτητα(Hz)')
ylabel('Px')


%//////////////2o meros
k=1000;
Pxnew=zeros(1,Nf);
for (i=1:k)

%Creation of series of symbols
b = (sign(randn(N_bits,1))+1)/2;
Xn = bits_to_2PAM(b);
sum_upsample=(1/Ts)*upsample(Xn,over);
t_sum_upsample=0:Ts:N_bits*T-Ts;
X=conv(sum_upsample,srrc)*Ts;
t_X=t(1)+t_sum_upsample(1):Ts:t(end)+t_sum_upsample(end);

% v=(var(X));
% Sx=(F/T).*v;

t_Px=linspace(-Fs/2,Fs/2,Nf);
F_X_abs=abs(fftshift(fft(X,length(t_Px)))*Ts).^2;
Ttotal=length(X).*Ts;
Px=F_X_abs/Ttotal;
Pxnew=Pxnew+Px;

end
Epxnew=Pxnew./k;
figure;
semilogy(t_Px,Sx,'b',t_Px,Epxnew,'r'); %thewritiki mple, mesi timi kokkini
title('Φασματική πυκνότητα ισχύος πάνω σε K=100 υλοποιήσεις περιοδογραμμάτων,σε 2PAM')
xlabel('συχνοτητα(Hz)')
ylabel('Px(meso)')

%///////////////////A.4/////////////////////////////

 Nf=2048;
 T = 0.001;
 over =10;
 Ts = T/over;
 N_bits=50;
 Fs = 1/Ts;


 b2 = (sign(randn(N_bits,1))+1)/2;
 Xn2 = bits_to_4PAM(b2);
 t2=0:N_bits/2-1;
%figure()
 %%stem(t2,Xn2);
 sum_upsample2=(1/Ts)*upsample(Xn2,over);
 t_sum_upsample2=0:Ts:((N_bits*T-Ts)/2);
 
 X2=conv(sum_upsample2,srrc)*Ts;
 t_X2=t(1)+t_sum_upsample2(1):Ts:t(end)+t_sum_upsample2(end);
 figure();
 plot(t_X2,X2);
 title('X(t),σε 4PAM')
 xlabel('χρονος(s)')
 ylabel('Χ(t)')

 t_Px2=linspace(-Fs/2,Fs/2,Nf);
 F_X_abs2=abs(fftshift(fft(X2,length(t_Px2)))*Ts).^2;
 Ttotal2=length(X2).*Ts;
 Px2=F_X_abs2/Ttotal2;
 
    

var2= 1/4*(1)^2  +  1/4*(-1)^2 +  1/4*(3)^2  +  1/4*(-3)^2 ;
Sx2=(F/T).*var2 ;
figure;
semilogy(t_Px2,Px2); %thewritiki mple mesi timi kokkini
title('Φασματική πυκνότητα ισχύος σε 4PAM')
xlabel('συχνοτητα(Hz)')
ylabel('Px2')
legend('simology')

figure;
plot(t_Px2,Px2); %thewritiki mple mesi timi kokkini
title('Φασματική πυκνότητα ισχύος σε 4PAM')
xlabel('συχνοτητα(Hz)')
ylabel('Px2')
legend('plot')


k=1000;
 Nf=2048;
 T = 0.001;
 over =10;
 %Ts = T/over
 N_bits=50;
Fs = 1/Ts;
Pxnew2=zeros(1,Nf);
for (i=1:k)
  

  b2 = (sign(randn(N_bits,1))+1)/2;
  Xn2 = bits_to_4PAM(b2);
  sum_upsample2=(1/Ts)*upsample(Xn2,over);
  t_sum_upsample2=0:Ts:((N_bits*T-Ts)/2);
   X2=conv(sum_upsample2,srrc)*Ts;
   t_X2=t(1)+t_sum_upsample2(1):Ts:t(end)+t_sum_upsample2(end);
   
    t_Px2=linspace(-Fs/2,Fs/2,Nf);
    F_X_abs2=abs(fftshift(fft(X2,length(t_Px2)))*Ts).^2;
    Ttotal2=length(X2).*Ts;
    Px2=F_X_abs2/Ttotal2;
    Pxnew2=Pxnew2+Px2;
    
end
%var2= 1/4*(1)^2  +  1/4*(-1)^2 +  1/4*(3)^2  +  1/4*(-3)^2 ;
%Sx2=(F/T).*var2 ;
Pxnew2=Pxnew2/k;
figure;
semilogy(t_Px2,Sx2,'b',t_Px2,Pxnew2,'r'); %thewritiki mple mesi timi kokkini
title('Φασματική πυκνότητα ισχύος πάνω σε K=1000 υλοποιήσεις περιοδογραμμάτων,σε 4PAM')
xlabel('συχνοτητα(Hz)')
ylabel('Px2(meso)')
legend('T= 0.001')
%////////////////A.5//////////////////////

 k=1000;
 Nf=2048;
 Tnew = 2*0.001;
 over =10;
 %Ts = T/over
 N_bits=50;
 Fs = 1/Ts;


 b2 = (sign(randn(N_bits,1))+1)/2;
 Xn2 = bits_to_4PAM(b2);
 t2=0:N_bits/2-1;
 sum_upsample2=(1/Ts)*upsample(Xn2,over);
 t_sum_upsample2=0:Ts:((N_bits*Tnew-Ts)/2);
 
 X2=conv(sum_upsample2,srrc)*Ts;
 t_X2=t(1)+t_sum_upsample2(1):Ts:t(end)+t_sum_upsample2(end);

 t_Px2=linspace(-Fs/2,Fs/2,Nf);
 F_X_abs2=abs(fftshift(fft(X2,length(t_Px2)))*Ts).^2;
 Ttotal2=length(X2).*Ts;
 Px2=F_X_abs2/Ttotal2;
 
    

var2= 1/4*(1)^2  +  1/4*(-1)^2 +  1/4*(3)^2  +  1/4*(-3)^2 ;
Sx2=(F/Tnew).*var2 ;
figure;
semilogy(t_Px2,Px2); %thewritiki mple mesi timi kokkini
title('Φασματική πυκνότητα ισχύος σε 4PAM σε simology')
xlabel('συχνοτητα(Hz)')
ylabel('Px2')
legend('T=0.002')

figure;
plot(t_Px2,Px2); %thewritiki mple mesi timi kokkini
title('Φασματική πυκνότητα ισχύος σε 4PAM σε plot')
xlabel('συχνοτητα(Hz)')
ylabel('Px2')
legend('T=0.002')




Pxnew2=zeros(1,Nf);
for (i=1:k)
  

  b2 = (sign(randn(N_bits,1))+1)/2;
  Xn2 = bits_to_4PAM(b2);
  sum_upsample2=(1/Ts)*upsample(Xn2,over);
  t_sum_upsample2=0:Ts:((N_bits*Tnew-Ts)/2);
   X2=conv(sum_upsample2,srrc)*Ts;
   t_X2=t(1)+t_sum_upsample2(1):Ts:t(end)+t_sum_upsample2(end);
   
    t_Px2=linspace(-Fs/2,Fs/2,Nf);
    F_X_abs2=abs(fftshift(fft(X2,length(t_Px2)))*Ts).^2;
    Ttotal2=length(X2).*Ts;
    Px2=F_X_abs2/Ttotal2;
    Pxnew2=Pxnew2+Px2;
    
end
%var2= 1/4*(1)^2  +  1/4*(-1)^2 +  1/4*(3)^2  +  1/4*(-3)^2 ;
%Sx2=(F/T).*var2 ;
Pxnew2=Pxnew2/k;
figure;
semilogy(t_Px2,Sx2,'b',t_Px2,Pxnew2,'r'); %thewritiki mple mesi timi kokkini
title('Φασματική πυκνότητα ισχύος πάνω σε K=1000 υλοποιήσεις περιοδογραμμάτων,σε 4PAM')
xlabel('συχνοτητα(Hz)')
ylabel('Px2(meso)')
legend('T=0.002')


