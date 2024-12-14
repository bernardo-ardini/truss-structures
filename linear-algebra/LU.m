function [L,U]=LU(A)
    n=size(A,2);
    
    for i=1:n
        for j=2:i
            A(i,j-1)=A(i,j-1)/A(j-1,j-1);
            for k=1:j-1
                A(i,j)=A(i,j)-A(i,k)*A(k,j);
            end
        end
        for j=i+1:n
            for k=1:i-1
                A(i,j)=A(i,j)-A(i,k)*A(k,j);
            end
        end
    end

    U=triu(A);
    L=tril(A,-1);
    L=L+eye(n);
end