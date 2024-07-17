syms t s Y

right_hand_side_eq = 4*exp(-6*t);
y_init = 6;
yprime_init = 0;
ydoubleprimecoeff = 0;
yprime_coeff = 1;
y_coeff = 7;


RHS = laplace(right_hand_side_eq);
Y1 = s * Y - y_init;
Y2 = s * Y1 - yprime_init;
Ys = solve(ydoubleprimecoeff * Y2 + yprime_coeff * Y1 + y_coeff * Y == RHS, Y);
Ys = simplify(Ys);
Pf = partfrac(Ys, 'FactorMode', 'full');
iL = ilaplace(Ys);

disp('Ys: ');
disp(Ys);
disp('Partial Fractions: ');
disp(Pf);
disp('Inverse Laplace: ');
disp(iL);
