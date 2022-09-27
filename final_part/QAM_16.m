function [ num_of_symbol_errors, num_of_bit_errors] = QAM_16(N, SNR_db)
%%%%%%%%% A.2  %%%%%%%%%%
b = (sign(randn(4*N,1))+1)/2;
b=b';
bI = b(1:2*N);
bQ = b(2*N+1:4*N);
A=1;
%%%%%%%%% A.3  %%%%%%%%%%%
XI= bits_to_4PAM(bI,A);
XQ= bits_to_4PAM(bQ,A);
%%%%%%%%% A.4  %%%%%%%%%%%
T=1;
over=10;
Ts=T/over;
Fs=1/Ts;
a=1;
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
tmin_XQ=XQ_up_t(1)+t(1);
tmax_XQ=XQ_up_t(end)+t(end);
XQ_conv_t=linspace(tmin_XQ,tmax_XQ,length(XQ_conv));
%%%%%%%%%%% A.5  %%%%%%%%%%%%%%%
Fo=2;
XI_mod=XI_conv.*(2*cos(2*pi*Fo*XI_conv_t));
XQ_mod=XQ_conv.*(-2*sin(2*pi*Fo*XQ_conv_t));
%%%%%%%%%%% A.6  %%%%%%%%%%%%%%%
X_mod=XQ_mod+XI_mod;
%%%%%%%%%%%% A.7  %%%%%%%%%%%%%%%
%exoume idaniko kanali
%%%%%%%%%%%% A.8  %%%%%%%%%%%%%%%
s_w=(10*(A^2))/(Ts*( 10 ^(SNR_db/10)));   %s(W) = ó_W ^ 2
W_noise=sqrt(s_w)*randn(1,length(X_mod));
Y = X_mod + W_noise;
%%%%%%%%%%%% A.9  %%%%%%%%%%%%%%%
Y_t=XI_conv_t;
YI=Y.*cos(2*pi*Fo*Y_t);  %Fo=2 apo parapanw
YQ=Y.*(-sin(2*pi*Fo*Y_t));
%%%%%%%%%%%%%%%%%%%%%% A.10 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
YI_conv = conv(srrc,YI)*Ts;
YQ_conv = conv(srrc,YQ)*Ts;
%%%%%%%%%%%%%%%% A.11%%%%%%%%%%%%%%%
YI_conv_cut=YI_conv(2*Fs*A:1:(length(YI_conv)-1)-2*Fs*A );
YQ_conv_cut=YQ_conv(2*Fs*A:1:(length(YI_conv)-1)-2*Fs*A );
YI_cut=downsample(YI_conv_cut,over);
YQ_cut=downsample(YQ_conv_cut,over);
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
%%%%%%%%%%%% A.14  %%%%%%%%%%%%%%%
bI_new = PAM_4_to_bits(YI_detect,A);
bQ_new = PAM_4_to_bits(YQ_detect,A);
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
end