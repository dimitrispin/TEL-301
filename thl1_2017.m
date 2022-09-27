clear all;
close all;
clc;
% A.1)
T = 0.01;
over = 10;
Ts = T/over ;
A = 5;
a = [0 0.5 1]; %roll-off
%Square Root Raised Cosine
[srrc_blue,t]= srrc_pulse(T, Ts, A, a(1));
[srrc_black,t] = srrc_pulse(T, Ts, A, a(2));
[srrc_red,t] = srrc_pulse(T, Ts, A, a(3));

figure();
plot(t,srrc_blue,'b', t,srrc_black,'k', t,srrc_red,'r');
xlabel('Time(s)');
ylabel('srrcblue srrcblack  srrcred ');
title('diagram 1:SRRC roll-off a=0(blue),a=0.5(black),a=1(red)');





%% Erothma A.2)
Fs = 1/Ts ; %Periodos
Nf = 1024 ;
step = Fs/Nf ; %Vhma
%Dhmiourgoume tous treis Metasxhmatismous Fourier ton palmon me thn
%entolh fft kai tous kentraroume sto mhden
Fourier_blue = fftshift(fft(srrc_blue,Nf))*Ts;
Fourier_black= fftshift(fft(srrc_black,Nf))*Ts;
Fourier_red= fftshift(fft(srrc_red,Nf))*Ts;

%Fasmatiki puknotita energeias ton palmon
Fasmatimi_puknotita_blue = (abs(Fourier_blue)).^2 ;
Fasmatimi_puknotita_black = (abs(Fourier_black)).^2 ;
Fasmatimi_puknotita_red = (abs(Fourier_red)).^2 ;
%aksonas suxnothton
t_Suxn = [-Fs/2:step:Fs/2-step];

figure();
%Sxediazoume sto idio plot tis fasmatikes pyknothtes energeias ton trion palmon
plot(t_Suxn,Fasmatimi_puknotita_blue,'b',t_Suxn,Fasmatimi_puknotita_black,'k',t_Suxn,Fasmatimi_puknotita_red,'r' );
xlabel('Frequency(HZ)');
ylabel('|Ö(F)|^2 ton trion palmon');
title('diagram 2:Fasmatikh pyknothta energeias');
figure();
%Sxediazoume se koino semilogy tis fasmatikes pyknothtes energeias ton trion palmon
semilogy(t_Suxn,Fasmatimi_puknotita_blue,'b',t_Suxn,Fasmatimi_puknotita_black,'k',t_Suxn,Fasmatimi_puknotita_red,'r');
xlabel('Frequency(HZ))');
ylabel('|Ö(F)|^2 ton trion palmon');
title('diagram 3:Fasmatikh puknothta energeias sth logarithmikh klimaka(mesemilogy)');








%% erothma A.3)
%Theorhtiko eyros fasmatos kathenos apo tous parapano palmous(xoris na exoun apokopei sto xrono
% T = 0.01;
% a = [0 0.5 1]; %syntelesths roll-off
BW_blue = (1+a(1))/(2*T)
BW_black = (1+a(2))/(2*T)
BW_red = (1+a(3))/(2*T)
%Dhmiourgoume tis 2 orizonties grammes pou mas zhteitai (C0=T/100 kai C1=T/(10^6))
for k=1:length(t_Suxn)
 c3(k)=T/10^3;
 c5(k)=T/10^5;
end
figure();
%Sxediazoume sto koino semilogy tis fasmatikes pyknothtes energeias
%me thn orizontia grammh C0=T/1000
semilogy(t_Suxn,Fasmatimi_puknotita_blue,'b',t_Suxn,Fasmatimi_puknotita_black,'k',t_Suxn,Fasmatimi_puknotita_red,'r',t_Suxn,c3,'g');
xlabel('Frequency(HZ))');
ylabel('|Ö(F)|^2 ton trion palmon');
title('diagram 4:Fasmatikh puknothta energeias me semilogy gia apokomenes times kato apo c0=T/10^3');
figure();
%Sxediazoume sto koino semilogy tis fasmatikes pyknothtes energeias
%me thn orizontia grammh C5=T/(10^5)
semilogy(t_Suxn,Fasmatimi_puknotita_blue,'b',t_Suxn,Fasmatimi_puknotita_black,'k',t_Suxn,Fasmatimi_puknotita_red,'r',t_Suxn,c5,'g');
xlabel('Frequency(HZ))'); 
ylabel('|Ö(F)|^2 ton trion palmon');
title('diagram 5:Fasmatikh puknothta energeias me semilogy gia apokomenes times katw apo c5=T/10^5');






%% Erothma B
%%% erothma B.1)
for i=0:2*A

 f_t_kt_blue = [zeros(1,i*T/Ts) srrc_blue(1:end-i*T/Ts)];
 multi_blue = srrc_blue .* f_t_kt_blue;
 
 f_t_kt_black = [zeros(1,i*T/Ts) srrc_black(1:end-i*T/Ts)];
 multi_black =  srrc_black.*f_t_kt_black;
 
 f_t_kt_red = [zeros(1,i*T/Ts) srrc_red(1:end-i*T/Ts)];
 multi_red = srrc_red.*f_t_kt_red;

 %Sxediazoume ta shmata mas gia 2 times tou k,K=0,1
 if (i==0)

   figure();
   subplot(2,3,1);
   plot(multi_blue,'b');
   xlabel('Time(sec)');
   title('diagram 6:Ö(t)*Ö(t-kT) me ton srrcblue gia k=0');
 
   subplot(2,3,2);
   plot(multi_black,'k');
   xlabel('Time(sec)');
   title('diagram 7:Ö(t)*Ö(t-kT) me ton srrcblack gia k=0');
  
   subplot(2,3,3);
   plot(multi_red,'r');
   xlabel('Time(sec)');
   title('diagram 8:Ö(t)*Ö(t-kT) me ton srrcred gia k=0');
   sum_multi_blue0  = sum(multi_blue)*Ts
   sum_multi_black0 = sum(multi_black)*Ts
   sum_multi_red0  = sum(multi_red)*Ts
   

 elseif(i==1)
   
    subplot(2,3,4);
   plot(multi_blue,'b');
   xlabel('Time(sec)');
   title('diagram 9:Ö(t)*Ö(t-kT) me ton srrcblue gia k=1');

   
   subplot(2,3,5);
   plot(multi_black,'k');
   xlabel('Time(sec)');
   title('diagram 10:Ö(t)*Ö(t-kT) me ton srrcblack gia k=1');

   
   subplot(2,3,6);
   plot(multi_red,'r');
   xlabel('Time(sec)');
   title('diagram 11:Ö(t)*Ö(t-kT) me ton srrcred gia k=1');
   
   sum_multi_blue1  = sum(multi_blue)*Ts
   sum_multi_black1 = sum(multi_black)*Ts
   sum_multi_red1   = sum(multi_red)*Ts
 end
end
for i=0:2*A
    
     f_t_kt_blue = [zeros(1,i*T/Ts) srrc_blue(1:end-i*T/Ts)];
     f_t_kt_black = [zeros(1,i*T/Ts) srrc_black(1:end-i*T/Ts)];
     f_t_kt_red = [zeros(1,i*T/Ts) srrc_red(1:end-i*T/Ts)];
 
   if(i==0)
       
    figure();
    subplot(2,3,1);
    plot(t,srrc_blue,'b',t,f_t_kt_blue,'c');
    xlabel('Time(sec)');
    title('diagram 12:Ö(t),Ö(t-kT) me ton srrcblue gia k=0');
   
    subplot(2,3,2);
    plot(t,srrc_black,'k', t,f_t_kt_black,'c');
    xlabel('Time(sec)');
    title('diagram 13:Ö(t),Ö(t-kT) me ton srrcblack gia k=0');
   
    subplot(2,3,3);
    plot(t,srrc_red,'r', t,f_t_kt_red,'c');
    xlabel('Time(sec)');
    title('diagram 14:Ö(t),Ö(t-kT) me ton srrcred gia k=0');
 
   elseif(i==1)
     
   subplot(2,3,4);
   plot(t,srrc_blue,'b', t,f_t_kt_blue,'c');
   xlabel('Time(sec)');
   title('diagram 15:Ö(t),Ö(t-kT) me ton srrcblue gia k=1');
   
   subplot(2,3,5);
   plot(t,srrc_black,'k', t,f_t_kt_black,'c');
   xlabel('Time(sec)');
   title('diagram 16:Ö(t),Ö(t-kT) me ton srrcblack gia k=1');
   
   subplot(2,3,6);
   plot(t,srrc_red,'r', t,f_t_kt_red,'c');
   xlabel('Time(sec)');
   title('diagram 17:Ö(t),Ö(t-kT) me ton srrcred gia k=1');
   
 
   end
end


 
% Erothma C
% erothma C.1)
%Pairnoume gia N=50 bits
T=1;
a=0.5;
N = 50;
over=10;
Ts = T/over ;
A = 5;
%Square Root Raised Cosine
[srrc,t] = srrc_pulse(T, Ts, A, a);


% erothma C.2) 2PAM systhma
%a)
%kai ta dhmiourgoume me thn b pou mas dinetai
b = (sign(randn(N,1))+1)/2;
%Kaloume th synarthsh bits_to_2PAM() kai
%apothikeyoume ta stoixeia ths sto X_2PAM
X_2PAM =bits_to_2PAM(b);
%b)
X_delta= (1/Ts)*upsample(X_2PAM,over);
%Pairnoume ta oria tou xronou apo 0 mexri (N-Ts)*T
Oria_x_2PAM = [0:Ts:(N-Ts)*T];
%apikonizw to xdelta san diakrito shma me thn voh8eia ths entolis stem
figure();
stem(Oria_x_2PAM,X_delta);
xlabel('Time(sec)');
title('diagram 18:Xdelta');

%ã
%kanw sineli3h ton dia shmaton me thn voh8eia ths entolhs conv kai meta
%pollaplasiazw me to me to vhma Ts gt h conv den periexei plhroforia gia
%ton xrono
X = conv(X_delta,srrc)*Ts;
%ta oria  ths x einai ta a8roismata ths arxhs kai tou telous ton duo shmatwn
%pou ekana shnelu3h me vhma Ts
oria_X =[(Oria_x_2PAM(1)+t(1)):Ts:(Oria_x_2PAM(end)+t(end))];
figure();
plot(oria_X,X);
xlabel('Time(sec)');
title('suneli3h Xdelta me to srrc palmo');
%d)
 %kanw suneli3h tou X me to srrc(ö(t)) kai oxi me to ö(-t) giati to
 %srrc(ö(t)) einai artio 
  Z = conv(X,srrc)*Ts;
  %Pairnoume  oria tou xronou gia to Z opws kai parapanw
  oria_Z = [oria_X(1)+t(1):Ts:oria_X(end)+t(end)];
  %vriskw to Z(KT) arxizw apo 2*A*over+1 ews 2*A*over+1+over*(N-1) ,
  %2*A gt gia na dhmiourgi8ei to Z exw kanei 2 shneli3eis me to srrc  pou einai sto 
  %diashma -A ws A ,to vhma over(T/Ts) giati toso apexoun h korifes meta3u tous
  %(N-1) einai apo thn ekfwnisi
  Zkt=Z(2*A*over+1:over:2*A*over+1+over*(N-1));%over=T/Ts
  Zkt
  figure();
%Sxediazoume to shma Z
 plot(oria_Z,Z);
 %sxediazoume to Zkt me to stem mas dinete sthn ekfonish sto idio figure me
 %to Z
  hold on
  stem([0:N-1]*T,Zkt);
  xlabel('Time(sec)');
  title('diagram 19: Z(kt) se stem kai to Z(t) me plot');
  hold off




