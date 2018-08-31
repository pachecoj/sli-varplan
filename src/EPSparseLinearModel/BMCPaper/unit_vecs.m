function X = unit_vecs(N,k,m);
% return m random N x 1 vectors with k ones
% does not !! cvheck whether two vectors are equal

if k == 1,
    if m <= N,
        X = eye(N);
    else
        fold = floor(m/N) + 1;
        X = repmat(eye(N),1,fold);
    end
    X = X(:,1:m);
else
    X = zeros(N,m);
    for i = 1:m
        ii = randperm(N);
        ii = ii(1:k);
        X(ii,i) = 1;
    end
end
    
    

% X = [];
% for i = 1:k
%     X = [X unit_vecsn(N,i)];
% end
% 
% 
% function X = unit_vecsn(N,k);
% if k == 1
%     X = eye(N);
% else
%     X = [];
%     for i = 1:N
%         XX = unit_vecsn(N-i,k-1);
%         L = size(XX,2);
%         X = [X [zeros(i-1,L);ones(1,L);XX]];
%     end
% end
