syms x y
a = [1,2,-1]
b = [0,2,3]
c = [a;b];
c
c(1,2)
indexs = zeros(1,size(a,2));
for i =1:size(a,2)
    if c(:,i) >= 0 == ones(1,3)
            indexs(i)=1;
    end
end
c(:,indexs==1)
[]==0