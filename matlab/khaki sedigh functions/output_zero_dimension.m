%% Output Zero Direction
G_S_backup = G_S;
k = 1;
syms s
G_Z = subs(G_S_backup, s, -0.75); % Calculating G(z)
y_zH1 = null(G_Z')'; % Calculating Output zero direction
y_zH1 = y_zH1/sqrt(y_zH1*y_zH1'); % Normalized Output zero direction
for i = 1:size(y_zH1, 2)
y_zH(i) = sym2poly(y_zH1(i));
end
y_zH