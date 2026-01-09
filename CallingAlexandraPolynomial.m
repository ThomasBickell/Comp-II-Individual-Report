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


lst={{[1,1,1,1]},{[1,2,1,1],[1,2,1,1]},{[1,2,1,3],[2,3,2,1],[3,1,3,2]},{[1,3,1,2],[3,1,3,4],[4,1,4,2],[2,3,2,4]},{[3,1,3,5],[1,4,1,3],[4,2,4,1],[2,5,2,4],[5,3,5,2]},{[1,4,1,3],[3,1,3,5],[5,2,5,1],[4,3,4,2],[2,5,2,4]},{[1,4,1,3],[3,2,3,1],[2,5,2,6],[5,2,5,3],[4,1,4,6],[6,5,6,4]},{[1,4,1,5],[6,3,6,4],[5,3,5,2],[4,6,4,1],[3,1,3,2],[2,6,2,5]},{[1,3,1,4],[3,6,3,1],[4,6,4,5],[5,3,5,2],[6,1,6,2],[2,5,2,4]},{[1,5,1,4],[5,2,5,1],[2,6,2,5],[6,3,6,2],[3,7,3,6],[7,4,7,3],[4,1,4,7]},{[1,5,1,4],[4,2,4,1],[2,4,2,3],[3,7,3,6],[6,3,6,2],[7,6,7,5],[5,1,5,7]},{[3,7,3,1],[7,3,7,4],[4,6,4,7],[1,4,1,5],[5,1,5,2],[2,5,2,6],[6,2,6,3]},{[5,1,5,8],[1,5,1,4],[4,2,4,1],[2,4,2,3],[3,6,3,7],[7,2,7,3],[6,8,6,7],[8,6,8,5]},{[6,2,6,3],[7,6,7,5],[2,6,2,7],[5,8,5,7],[8,5,8,4],[3,1,3,8],[1,4,1,3],[4,2,4,1]},{[1,4,1,3],[4,2,4,1],[6,8,6,1],[3,6,3,5],[2,6,2,7],[8,2,8,3],[5,8,5,7],[7,5,7,4]},{[1,6,1,5],[5,1,5,9],[6,5,6,4],[4,7,4,6],[7,4,7,3],[3,8,3,7],[8,3,8,2],[2,9,2,8],[9,2,9,1]},{[1,6,1,5],[6,2,6,1],[2,1,2,9],[9,3,9,2],[3,9,3,8],[8,4,8,3],[5,8,5,7],[4,7,4,6],[7,5,7,4]},{[1,7,1,8],[7,1,7,2],[5,2,5,3],[7,4,7,5],[4,8,4,9],[1,5,1,6],[3,6,3,7],[6,9,6,1],[9,3,9,4]},{[1,6,1,7],[6,1,6,2],[9,6,9,5],[2,1,2,10],[7,3,7,2],[10,4,10,3],[3,8,3,7],[4,10,4,9],[5,9,5,8],[8,5,8,4]},{[7,10,7,1],[1,5,1,4],[8,7,8,6],[10,7,10,8],[3,9,3,10],[5,3,5,4],[4,2,4,1],[6,9,6,8],[2,5,2,6],[9,2,9,3]},{[5,2,5,1],[9,6,9,5],[2,8,2,9],[8,4,8,5],[6,2,6,3],[1,8,1,7],[3,10,3,9],[4,6,4,7],[7,1,7,10],[10,4,10,3]}}
time=zeros(2,length(lst));
for i = 1:length(lst)
    tic
    Alexandra(lst{i});
    t=toc;
    time(:,i)=[length(lst{i}),t];
end
scatter(time(1,:),time(2,:),100,'b','x','linewidth',2)
xlabel('Number of crossings')
ylabel('Time to complete (s)')
title('Calculating the Alexandra polynomial')
set(gca,'FontSize',25)
hold on
xFit=linspace(min(time(1,:)),max(time(1,:)),100);
p = polyfit(time(1,:),time(2,:), 3);
yFit = polyval(p, xFit);
plot(xFit, yFit, 'r-', 'LineWidth', 2)


legend('Data points','Best fit line','location','northwest')
