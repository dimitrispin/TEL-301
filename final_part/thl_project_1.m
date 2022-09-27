clear ;
close all;
clc;

%%%%%%%%%% A %%%%%%%%%%

%%%%%%%%% A.2  %%%%%%%%%%
N=100;
b = (sign(randn(4*N,1))+1)/2;
b=b';
bI = b(1:2*N);
%plot(bI)
bQ = b(2*N+1:4*N);
A=1;
%%%%%%%%% A.3  %%%%%%%%%%%
XI= bits_to_4PAM(bI,A);
XQ= bits_to_4PAM(bQ,A);
%  scatterplot(XI +1i*XQ)
%  grid on;
%%%%%%%%% A.4  %%%%%%%%%%%
T=1;
over=10;
Ts=T/over;
Fs=1/Ts;
a=1;
Nf=2048;
[srrc, t] = srrc_pulse(T,Ts,A,a);

XI_up_t = linspace(0,N*T,N*over);
XI_up=(1/Ts)*upsample(XI,over);

XQ_up_t = linspace(0,N*T,N*over);
XQ_up=(1/Ts)*upsample(XQ,over);

tmin_XI=XI_up_t(1)+t(1);
tmax_XI=XI_up_t(end)+t(end);
XI_conv = conv(srrc,XI_up)*Ts;
XI_conv_t=linspace(tmin_XI,tmax_XI,length(XI_conv));

XQ_conv = conv(srrc,XQ_up)*Ts;
%XQ_conv_t= XQ_up_t(1)+t(1):Ts:XQ_up_t(end)+t(end-1);
tmin_XQ=XQ_up_t(1)+t(1);
tmax_XQ=XQ_up_t(end)+t(end);
XQ_conv_t=linspace(tmin_XQ,tmax_XQ,length(XQ_conv));

figure(1);
plot(XI_conv_t,XI_conv)
title(' convolution XI(t)')
ylabel('XI(t)');
xlabel('Time t(s)');

figure(2);
plot(XQ_conv_t,XQ_conv)
title(' convolution XQ(t)')
ylabel('XQ(t)');
xlabel('Time t(s)');


Px_XI_t=linspace(-Fs/2,Fs/2,Nf);
fft_XI_conv=fftshift(fft(XI_conv,Nf)*Ts);
spec_XI_conv = abs(fft_XI_conv).^2;
Ttotal_XI_up=length(XI_up_t)*Ts;
Px_XI=spec_XI_conv/Ttotal_XI_up;

figure(3);
plot(Px_XI_t,Px_XI)
title('periodogram of XI(t)')
xlabel('Frequency (Hz)');
ylabel(' XI(t)');

Px_XQ_t=linspace(-Fs/2,Fs/2,Nf);
fft_XQ_conv=fftshift(fft(XQ_conv,Nf)*Ts);
spec_XQ_conv = abs(fft_XQ_conv).^2;
Ttotal_XQ_up=length(XQ_up_t)*Ts;
Px_XQ=spec_XQ_conv/Ttotal_XQ_up;

figure(4);
plot(Px_XQ_t,Px_XQ)
title('periodogram of XQ(t)')
xlabel('Frequency (Hz)');
ylabel(' XQ(t)');
%%%%%%%%%%% A.5  %%%%%%%%%%%%%%%
Fo=2;
XI_mod=XI_conv.*(2*cos(2*pi*Fo*XI_conv_t));
XQ_mod=XQ_conv.*((-2)*sin(2*pi*Fo*XQ_conv_t));

figure(5);
plot(XI_conv_t,XI_mod);  %exei to idio a3ona xronou me XI
title('XImod signal for XI');
xlabel('Time t(s)');
ylabel('XImod(t)');

figure(6);
plot(XQ_conv_t,XQ_mod);  %exei to idio xrono a3ona xronou me XQ
title('XQmod signal for XQ');
xlabel('Time t(s)');
ylabel('XQmod(t)');

Px_XI_mod_t=linspace(-Fs/2,Fs/2,Nf);
fft_XI_mod=fftshift(fft(XI_mod,Nf)*Ts);
spec_XI_mod = abs(fft_XI_mod).^2;
Ttotal_XI_mod=length(XI_mod)*Ts;
Px_XI_mod=spec_XI_mod/Ttotal_XI_mod;

figure(7);
plot(Px_XI_mod_t,Px_XI_mod)
title('periodogram of XImod(t)')
xlabel('Frequency (Hz)');
ylabel(' XImod(t)');

Px_XQ_mod_t=linspace(-Fs/2,Fs/2,Nf);
fft_XQ_mod=fftshift(fft(XQ_mod,Nf)*Ts);
spec_XQ_mod = abs(fft_XQ_mod).^2;
Ttotal_XQ_mod=length(XQ_mod)*Ts;
Px_XQ_mod=spec_XQ_mod/Ttotal_XQ_mod;

figure(8);
plot(Px_XQ_mod_t,Px_XQ_mod)
title('periodogram of XQmod(t)')
xlabel('Frequency (Hz)');
ylabel(' XQmod(t)');

%%%%%%%%%%% A.6  %%%%%%%%%%%%%%%
X_mod=XI_mod+XQ_mod;
figure(9);
plot(XI_conv_t,X_mod);  %exei to idio a3ona xronou me XI_conv-XQ_conv
title('Xmod signal for XI,XQ');
xlabel('Time t(s)');
ylabel('XImod+XQmod(t)');

Px_X_mod_t=linspace(-Fs/2,Fs/2,Nf);
fft_X_mod=fftshift(fft(X_mod,Nf)*Ts);
spec_X_mod = abs(fft_X_mod).^2;
Ttotal_X_mod=length(X_mod)*Ts;
Px_X_mod=spec_X_mod/Ttotal_X_mod;

figure(10);
plot(Px_X_mod_t,Px_X_mod)
title('periodogram of Xmod(t)')
xlabel('Frequency (Hz)');
ylabel(' Xmod(t)');

%%%%%%%%%%%% A.7  %%%%%%%%%%%%%%%
%exoume idaniko kanali
%%%%%%%%%%%% A.8  %%%%%%%%%%%%%%%
SNR_dB=20;
s_w=(10*(A^2))/(Ts*( 10 ^(SNR_dB/10)));   %s(W) = σ_W ^ 2
%white Gaussian noise
W_noise=sqrt(s_w).*randn(1,length(X_mod)); %%%den eimai sgros gia thn riza sqrt
Y = X_mod + W_noise;

% Y = awgn(X_mod,SNR_dB);
%  figure();
%  plot(XI_conv_t,X_mod)
%  figure();
%  plot(XI_conv_t,Y,)
%  figure();
%  plot(W_noise);
%%%%%%%%%%%% A.9  %%%%%%%%%%%%%%%
Y_t=XI_conv_t;
YI=Y.*cos(2*pi*Fo*Y_t);  %Fo=2 apo parapanw
YQ=Y.*(-sin(2*pi*Fo*Y_t));

figure(11);
plot(Y_t,YI);
title('ενθόρυβη κυματομορφή YI');
xlabel('Time t(s)');
ylabel('YI(t)');

figure(12);
plot(Y_t,YQ);
title('ενθόρυβη κυματομορφή YQ');
xlabel('Time t(s)');
ylabel('YQ(t)');


Px_YI_t=linspace(-Fs/2,Fs/2,Nf);
fft_YI=fftshift(fft(YI,Nf)*Ts);
spec_YI = abs(fft_YI).^2;
Ttotal_YI=length(YI)*Ts;
Px_YI=spec_YI/Ttotal_YI;

figure(13);
plot(Px_YI_t,Px_YI)
title('periodogram of YI(t)')
xlabel('Frequency (Hz)');
ylabel(' YI(t)');

Px_YQ_t=linspace(-Fs/2,Fs/2,Nf);
fft_YQ=fftshift(fft(YQ,Nf)*Ts);
spec_YQ = abs(fft_YQ).^2;
Ttotal_YQ=length(YQ)*Ts;
Px_YQ=spec_YQ/Ttotal_YQ;

figure(14);
plot(Px_YQ_t,Px_YQ)
title('periodogram of YQ(t)')
xlabel('Frequency (Hz)');
ylabel(' YQ(t)');

%%%%%%%%%%%% A.10  %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% YI_down=(1/Ts)*downsample(YI,over);
% YQ_down=(1/Ts)*downsample(YQ,over);
% YI_down_t= linspace(0,N*T,N);
% YQ_down_t=linspace(0,N*T,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmin_YI=Y_t(1)+t(1);
tmax_YI=Y_t(end)+t(end);
YI_conv = conv(srrc,YI)*Ts;
YI_conv_t=linspace(tmin_YI,tmax_YI,length(YI_conv));

tmin_YQ=Y_t(1)+t(1);
tmax_YQ=Y_t(end)+t(end);
YQ_conv = conv(srrc,YQ)*Ts;
YQ_conv_t=linspace(tmin_YQ,tmax_YQ,length(YQ_conv));

figure(15);
plot(YI_conv_t,YI_conv)
title('convolution YIconv');
ylabel('YIconv(t)');
xlabel('Time t(s)');

figure(16);
plot(YQ_conv_t,YQ_conv)
title('convolution YQconv');
ylabel('YQconv(t)');
xlabel('Time t(s)');


Px_YI_t=linspace(-Fs/2,Fs/2,Nf);
fft_YI_conv=fftshift(fft(YI_conv,Nf)*Ts);
spec_YI_conv = abs(fft_YI_conv).^2;
Ttotal_YI=length(Y_t)*Ts;
Px_YI=spec_YI_conv/Ttotal_YI;

figure(17);
plot(Px_YI_t,Px_YI)
title('periodogram of YI(t)')
xlabel('Frequency (Hz)');
ylabel(' YI(t)');

Px_YQ_t=linspace(-Fs/2,Fs/2,Nf);
fft_YQ_conv=fftshift(fft(YQ_conv,Nf)*Ts);
spec_YQ_conv = abs(fft_YQ_conv).^2;
Ttotal_YQ=length(Y_t)*Ts;
Px_YQ=spec_YQ_conv/Ttotal_YQ;

figure(18);
plot(Px_YQ_t,Px_YQ)
title('periodogram of YQ(t)')
xlabel('Frequency (Hz)');
ylabel(' YQ(t)');

%%%%%%%%%%%% A.11  %%%%%%%%%%%%%%%
YI_conv_cut=YI_conv(2*Fs*A:(length(YI_conv)-1)-2*Fs*A );
YQ_conv_cut=YQ_conv(2*Fs*A:(length(YI_conv)-1)-2*Fs*A );
YI_cut=downsample(YI_conv_cut,over);
YQ_cut=downsample(YQ_conv_cut,over);
scatterplot(YI_cut+1i*YQ_cut);
grid on

%%%%%%%%%%%% A.12  %%%%%%%%%%%%%%%
YI_detect = detect_4_PAM(YI_cut,A);
YQ_detect= detect_4_PAM(YQ_cut,A);

%%%%%%%%%%%% A.13  %%%%%%%%%%%%%%%
num_of_symbol_errors=0;
for i=1:length(YI_detect)

     if( YI_detect(i)~=XI(i))
        num_of_symbol_errors = num_of_symbol_errors+ 1;
    end  
     if( YQ_detect(i)~=XQ(i))
        num_of_symbol_errors = num_of_symbol_errors+ 1;
    end 
end
num_of_symbol_errors

[numErrors,ser] = symerr([XI,XQ],[YI_detect,YQ_detect])
%%%%%%%%%%%% A.14  %%%%%%%%%%%%%%%
bI_new = PAM_4_to_bits(YI_detect,A);
bQ_new = PAM_4_to_bits(YQ_detect,A);
%    figure
%   plot(bI_new)
 

%%%%%%%%%%%% A.15  %%%%%%%%%%%%%%%

 num_of_bit_errors=0;
for i=1:1:length(bI)
    if (bI(i) ~= bI_new(i))
        num_of_bit_errors = num_of_bit_errors+1;
    end
    if ( bQ(i) ~= bQ_new(i))
        num_of_bit_errors = num_of_bit_errors+1;
    end
end
 num_of_bit_errors
 % close all;
 	[numErrors,ber]=biterr([bI,bQ],[bI_new,bQ_new]) 

