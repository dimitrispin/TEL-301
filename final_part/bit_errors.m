function [num_of_bit_errors] = bit_errors(b, b_new)
%h synarthsh ayth exei ws orismata ta bits kai thn ektimwmenh akolouthia
%bits pou vrhkame mesw ths detect 
num_of_bit_errors = 0;%counter 
for i=1:1:length(b)
    if b(i) ~= b_new(i) %oso diaferoun t bits dhladh exoume lathos sthn ektimhsh ayxanoume ton counter dhladh vgazoume posa lathh exoume kata mhkos ths akolouthias
        num_of_bit_errors = num_of_bit_errors+1;
    end
end
 
end