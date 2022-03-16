% Generate two data sets which are rank deficient
Y1 = [1, 2, 3, 3, 5, 8, 13, 9;
      1, 2, 3, 3, 5, 8, 1, 9;
      1, 2, 3, 3, 5, 8, 13, 9;
      1, 2, 3, 3, 5, 8, 13, 9];


Y2 = [2, 1, 3, 4, 6, 7, 9, 10;
      2, 1, 3, 4, 6, 7, 9, 10;
      2, 1, 3, 4, 6, 8, 9, 10];

[m1, n] = size(Y1);
[m2, ~] = size(Y2);
%% Direct GSVD (CS-D)
[U1, U2, X, S1, S2] = gsvd(Y1, Y2);
plot(X)
%%
[U1, S1, V1] = svd(Y1);
plot(diag(S1), 'o')
grid on

[U2, S2, V2] = svd(Y2);
plot(diag(S2), 'o')
grid on

Y = [Y1; Y2];
[U, S, V] = svd(Y);
plot(diag(S), 'o')
grid on

% Degree of truncation
r = 2;
% Projection onto the TSVD Basis
Yt = Y*V(:,1:r);

Y1t = Yt(1:m1,:);

Y2t = Yt(m1+1:m1+m2,:);

[Ub1,Ub2,Vb,Cb,Sb] = gsvd(Y1t,Y2t);

Phi = V(:,1:r)*Vb;

Phi1 = Phi(1:m1,:);
Phi2 = Phi(m1+1:m1+m2,:);

plot(normalize(V(1:m1, 1:r)))
hold on
plot(normalize(Phi1))