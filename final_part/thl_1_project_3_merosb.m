clear ;
close all;
clc;
%%%%%%%%%%%%%%%% B.1%%%%%%%%%%%

K=50;
N=100;
M=16;
SNR_db=0:2:M;
bps=log2(M); % bits per symbol
symbols=N*K;
bits=4*N*K;
P_symbol_error=ones(1,length(SNR_db));
P_bit_error=ones(1,length(SNR_db));
theory_symbol_error = ones(1,length(SNR_db));
theory_bit_error = ones(1,length(SNR_db));
i=0;
for SNRdb=0:2:M
    sum_symbol_e=0;
    sum_bit_e=0;
    i=i+1;
    for j=1:K
        [num_of_symbol_errors, num_of_bit_errors] = QAM_16(N, SNRdb);
        sum_symbol_e=sum_symbol_e+num_of_symbol_errors;
        sum_bit_e=sum_bit_e+num_of_bit_errors;
    end
    P_symbol_error(i)=sum_symbol_e/symbols; %pi8anothta la8os sumvolou
    P_bit_error(i)=sum_bit_e/bits;       %pi8anothta la8os bit
 
    s_w = (10*(1^2))/(0.1*10^(SNRdb/10));
    s_n = sqrt(s_w*0.1/2);  %T=0.1
    theory_symbol_error(i) = 3*Q(1/s_n);   %pi8anothta la8os sumvolou 8ewritika
    theory_bit_error(i)=theory_symbol_error(i)/bps;  %pi8anothta la8os bit 8ewritika
end

%%%%%%%%%%%%%%%  B.2 %%%%%%%%%%%
figure(1);
semilogy(SNR_db,P_symbol_error);
hold on;
semilogy(SNR_db,theory_symbol_error,'green');
title('πιθανότητα σΦάλματος συμβόλου για 16-QAM');
xlabel('SNR in dB');
ylabel('πιθανότητα σΦάλματος συμβόλου');
legend('πειραματικό SER','θεωρητικό SER');
%%%%%%%%%%%% B.3 %%%%%%%%%%%%%%%
figure(2);
semilogy(SNR_db,P_bit_error);
hold on;
semilogy(SNR_db,theory_bit_error,'green');
title('πιθανότητα σΦάλματος bit για 16-QAM');
xlabel('SNR in dB');
ylabel('πιθανότητα σΦάλματος bit');
legend('πειραματικό SER','θεωρητικό SER');