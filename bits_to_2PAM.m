function [Symbol] = bits_to_2PAM(b)
%Pairnoume tis times ton bits kai an einai 0 exoume symbolo +1 enw sthn
%antitheth periptosh (an einai 1) exoume symbolo -1
for i=1:length(b)
 if b(i)==1
  Symbol(i)=-1;
 else 
  Symbol(i)=+1;
 end
end
return;
end