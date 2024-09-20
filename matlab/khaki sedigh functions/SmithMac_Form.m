function [G_S, SA, P, Z, Poles, Zeros] = SmithMac_Form(G)
s = sym('s');
switch class(G)
    case 'tf'
        G_TF = G; % Transfer Function Form
        G_TF_backup = G;
        [num, den] = tfdata(G_TF);
        for i = 1:size(num, 1) % number of Outputs
            for j = 1:size(num, 2) % number of Inputs
                num_S(i, j) = poly2sym(num{i, j}, s);
                den_S(i, j) = poly2sym(den{i, j}, s);
            end
        end
        G_S = simplify(num_S ./ den_S);
        G_S_backup = G_S;
    case 'sym'
        G_S = simplify(G); % Symbolic Form
        G_S_backup = G_S;
end
[num_S, den_S] = numden(G_S);
row = size(G_S,1);
col = size(G_S,2);
n = min(row,col);
minors = cell(1,n);
D0 = 1;
D0 = sym(D0);
D = sym(NaN(1,n));
invFact = sym(NaN(1,n));
for i = 1:n % Different Submatrices
    rowindex = false(1,row);
    rowindex(1:i) = true;
    rowperms = unique(perms(rowindex),'rows');
    colindex = false(1,col);
    colindex(1:i) = true;
    colperms = unique(perms(colindex),'rows');
    rownum = size(rowperms,1);
    colnum = size(colperms,1);
    minors{i} = sym(NaN(rownum,colnum));
    for j=1:rownum % Minors' Calculation
        for k=1:colnum
            Atmp = G_S;
            Atmp = Atmp(rowperms(j,:),:);
            Atmp = Atmp(:,colperms(k,:));
            minors{i}(j,k) = det(Atmp);
        end
    end
    rowlen = rownum*colnum; %(row - (i-1))*(col - (i-1));
    minors{i} = reshape(minors{i},1,rowlen);
    minors{i}(minors{i} == 0) = [];
    D(i) = gcd(minors{i});
    if i == 1 % Invarint factors' calculation
        invFact(i) = D(i)/D0;
    else
        invFact(i) = D(i)/D(i-1);
    end
end
SA = diag(invFact);
for i = 1:n
    f = factor(SA(i, i));
    SA(i, i) = prod(f);
    El(i) = SA(i, i);
end
[N, D] = numden(El);
Z = prod(N);
P = prod(D);
Zeros = (roots(sym2poly(Z)));
Poles = (roots(sym2poly(P)));
if row>col
    zerorows = zeros(row-col,col);
    SA = [SA;zerorows];
elseif col>row
    zerocols = zeros(row,col-row);
    SA = [SA, zerocols];
end
end