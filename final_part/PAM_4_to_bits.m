function [est_bit] = PAM_4_to_bits(X, A)
  counter=1; 
  for i=1:length(X)  
        if (X(i)==A)% i decide that i send the symbol +A 
               est_bit(counter)=1;
               est_bit(counter+1)=1;    
        elseif (X(i)==-A )% i decide that i send the symbol -A 
               est_bit(counter)=0;
               est_bit(counter+1)=1;  
                
        elseif ( X(i)==3*A )% i decide that i send the symbol 3A
               est_bit(counter)=1;
               est_bit(counter+1)=0; 
               
        elseif (X(i)==-3*A )%i decide that i send the symbol -3A
               est_bit(counter)=0;
               est_bit(counter+1)=0; 
        end       
        counter=counter+2;     
   end
end