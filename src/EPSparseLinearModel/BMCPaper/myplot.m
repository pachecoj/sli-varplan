function myplot(x,y,lb,ub,skip,dskip,labels)

N = length(x);
if nargin < 5
    skip = round(N/20);
end
if nargin < 6
    dskip = skip * 3;
end

ii = 1:dskip:N-size(y,2)+1;
ii2 = 1:skip:N-size(y,2)+1;


for i = 1:size(y,2)
    iii = ii2 + i - 1;
    iii = [1 iii N];
    plot(x(iii),y(iii,i),[col(i) '-' typ(i)],'MarkerFaceColor',col(i),'Markersize',8); %3
%     plot(x,y(:,i),[col(i) '-' typ(i)],'MarkerFaceColor',col(i),'Markersize',8);
end
for i = 1:size(y,2)
    iii = ii + i - 1;
    errorbar(x(iii),y(iii,i),lb(iii,i),ub(iii,i),[col(i) '.']);
end
for i = size(y,2):-1:1
    iii = [1 iii N];
    iii = ii2 + i - 1;
    plot(x(iii),y(iii,i),[col(i) '-' typ(i)],'MarkerFaceColor',col(i),'Markersize',8);%3
end

function c = col(i)
base_col = {'b','g','r','c','m','k'};
c = base_col{mod(i-1,length(base_col)) + 1};

function t = typ(i)
base_typ = {'s','^','v','>','<','*','+','x','v','>'};
t = base_typ{mod(i,length(base_typ)) + 1 + floor(i/length(base_typ))};
