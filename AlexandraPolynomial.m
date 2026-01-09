X={[1,3,1,2],[3,1,3,4],[4,1,4,2],[2,3,2,4]};
function A = Alexandra(X)
n=length(X);
%Below is getting the matrix, this may be very large and sparse for high n
syms t
lst=[1,t,-t,-1];
AMatrix=repmat(t,n,n);
AMatrix(:)=AMatrix(:)-t;
for i = 1:n
    for j = 1:4
        AMatrix(i,X{i}(j))=AMatrix(i,X{i}(j)) + lst(j);
    end
end

AMatrix(:,1)=[];
AMatrix(1,:)=[];
% Below is getting the determinent of the matrix via Gaussian elimination
% detmult being the multiplier after the different operations.
detmult = sym(1);
for i = 1:n-1
    if AMatrix(i,i)==0
        swap=AMatrix(i,:);
        AMatrix(i,:)=AMatrix(i+1,:);
        AMatrix(i+1,:)=swap;
        detmult=detmult*-1;
    end
    if AMatrix(i,i)~=1
        mult=AMatrix(i,i);
        AMatrix(i,:)=AMatrix(i,:)/AMatrix(i,i);
        detmult=detmult*mult;
    end
    for j = 1:n-1-i
        mult=AMatrix(i+j,i);
        AMatrix(i+j,:)=AMatrix(i+j,:)-mult*AMatrix(i,:);
    end
end
A=simplify(detmult);
A=expand(A);
end
Alexandra(X)