function [Xnew]=bits_to_4PAM(b)
m=length(b);
counter=1;
for i=1:2:(m)
    if (b(i)==0 & b(i+1)==0)
            Xnew(counter)=+3;
    end
    if (b(i)==0 & b(i+1)==1)
            Xnew(counter)=+1;
    end
    if (b(i)==1 & b(i+1)==1)
            Xnew(counter)=-1;
    end
    if (b(i)==1 & b(i+1)==0)
            Xnew(counter)=-3;
    end
    counter=counter+1;
end