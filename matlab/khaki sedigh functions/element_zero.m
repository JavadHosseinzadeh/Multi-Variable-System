%% Element Zeros
row = size(G1, 1);
col = size(G1, 2);
EL_Z = cell(row,col);
for i = 1:row
for j = 1:col
[N, D] = numden(G_S(i, j));
if N==0 | size(sym2poly(N))==1 % Seperating the exceptions
continue
else
EL_Z{i, j}=solve(G_S(i, j)==0); % Calculating zero elements
end
end
end
EL_Z = EL_Z(~cellfun('isempty', EL_Z)); % Deleting empty cells
Element_Zeros = [];
for i = 1:size(EL_Z, 1)
Element_Zeros(i, :) = EL_Z{i};
end
Element_Zeros