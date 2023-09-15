clear all;
close all;
clc;

mu = 1;
a  = 1;
c = 1000;

S = [10000, 0, 0;
    0, 10000, 0;
    0, 0, 10000;
    10000, 10000, 10000];


S2 = [1000, 0, 0;
    0, 1000, 0;
    0, 0, 1000;
    1000, 1000, 1000;
    1000, 1000, 0;
    0, 1000, 1000;
    1000, 0, 1000;
    -1000, -1000, -1000];
xt = [1000,1000,1000];

t = vecnorm(S-xt,2,2)/c + [0;0.001;0;0];

% t = [ 9;
%     10;
%     10;
%     17.3205];

t2 = [ 1;
    1;
    1;
    1.73;
    1.4142;
    1.4142;
    1.4142;
    1.732];


x0 = [0;
    0;
    0];
x02 = [10;
    10;
    10];


rho = t*c;

rho2 = t2*c;

for i=1:1
S-x0';
dp = vecnorm(S - x0',2,2) - c*t 

% dp2 = vecnorm(S2 - x02',2,2) - c*t2; 

A = [(x0(1)-S(1,1))/rho(1),(x0(2)-S(1,2))/rho(1),(x0(3)-S(1,3))/rho(1),c;
     (x0(1)-S(2,1))/rho(2),(x0(2)-S(2,2))/rho(2),(x0(3)-S(2,3))/rho(2),c;
     (x0(1)-S(3,1))/rho(3),(x0(2)-S(3,2))/rho(3),(x0(3)-S(3,3))/rho(3),c;
     (x0(1)-S(4,1))/rho(4),(x0(2)-S(4,2))/rho(4),(x0(3)-S(4,3))/rho(4),c];

A2 = [(x0(1)-S2(1,1))/rho2(1),(x0(2)-S2(1,2))/rho2(1),(x0(3)-S2(1,3))/rho2(1),c;
     (x0(1)-S2(2,1))/rho2(2),(x0(2)-S2(2,2))/rho2(2),(x0(3)-S2(2,3))/rho2(2),c;
     (x0(1)-S2(3,1))/rho2(3),(x0(2)-S2(3,2))/rho2(3),(x0(3)-S2(3,3))/rho2(3),c;
     (x0(1)-S2(4,1))/rho2(4),(x0(2)-S2(4,2))/rho2(4),(x0(3)-S2(4,3))/rho2(4),c;
     (x0(1)-S2(5,1))/rho2(5),(x0(2)-S2(5,2))/rho2(5),(x0(3)-S2(5,3))/rho2(5),c;
     (x0(1)-S2(6,1))/rho2(6),(x0(2)-S2(6,2))/rho2(6),(x0(3)-S2(6,3))/rho2(6),c;
     (x0(1)-S2(7,1))/rho2(7),(x0(2)-S2(7,2))/rho2(7),(x0(3)-S2(8,3))/rho2(7),c;
     (x0(1)-S2(8,1))/rho2(8),(x0(2)-S2(8,2))/rho2(8),(x0(3)-S2(8,3))/rho2(8),c];
 
 
dx = inv(A'*A)*A'*dp
% dx2 = inv(A2'*A2)*A2'*dp2;

x0 = x0-dx(1:3);
% x02 = x02-dx2(1:3);

end
C = inv(A'*A);
eigs(C)
sum(eig(C))



C2 = inv(A2'*A2);
eigs(C2)
sum(eig(C2))
x0
x02

%%
% 
% randn('state',13);
% n = 6;
% P0 = randn(n); P0 = P0'*P0 + eps*eye(n);
% P1 = randn(n); P1 = P1'*P1;
% P2 = randn(n); P2 = P2'*P2;
% P3 = randn(n); P3 = P3'*P3;
% q0 = randn(n,1); q1 = randn(n,1); q2 = randn(n,1); q3 = randn(n,1);
% r0 = randn(1); r1 = randn(1); r2 = randn(1); r3 = randn(1);
% 
% fprintf(1,'Computing the optimal value of the QCQP and its dual... ');
% 
% cvx_begin
%     variable x(n)
%     dual variables lam1 lam2 lam3
%     minimize( 0.5*quad_form(x,P0) + q0'*x + r0 )
%     lam1: 0.5*quad_form(x,P1) + q1'*x + r1 <= 0;
%     lam2: 0.5*quad_form(x,P2) + q2'*x + r2 <= 0;
%     lam3: 0.5*quad_form(x,P3) + q3'*x + r3 <= 0;
% cvx_end
% 
% obj1 = cvx_optval;
% P_lam = P0 + lam1*P1 + lam2*P2 + lam3*P3;
% q_lam = q0 + lam1*q1 + lam2*q2 + lam3*q3;
% r_lam = r0 + lam1*r1 + lam2*r2 + lam3*r3;
% obj2 = -0.5*q_lam'*inv(P_lam)*q_lam + r_lam;
% 
% fprintf(1,'Done! \n');
% 
% % Displaying results
% disp('------------------------------------------------------------------------');
% disp('The duality gap is equal to ');
% disp(obj1-obj2)
%%

dp = vecnorm(S - x0',2,2) - c*t;
A = [(x0(1)-S(1,1))/rho(1),(x0(2)-S(1,2))/rho(1),(x0(3)-S(1,3))/rho(1),c;
     (x0(1)-S(2,1))/rho(2),(x0(2)-S(2,2))/rho(2),(x0(3)-S(2,3))/rho(2),c;
     (x0(1)-S(3,1))/rho(3),(x0(2)-S(3,2))/rho(3),(x0(3)-S(3,3))/rho(3),c;
     (x0(1)-S(4,1))/rho(4),(x0(2)-S(4,2))/rho(4),(x0(3)-S(4,3))/rho(4),c];

n = 4;
r = 0.1;
B = [eye(3),zeros(3,1)];
y = [100,100,100]';
z = [100,100,100]';
l = 1.2;
u = 0.8;
P0 = A'*A + r*B'*B;
P1 = l*u*B'*B;




q0 = -1*(dp'*A+r*(z-y)'*B);
q1 = -(l*(y-u*z)+u*(y-l*z))'*B;

r0 = dp'*dp + r*(z-y)'*(z-y);
r1 = (y-l*z)'*(y-u*z);


fprintf(1,'Computing the optimal value of the QCQP and its dual... ');

cvx_begin
    variable x(n)
    dual variables lam1 lam2 lam3
    minimize( 0.5*quad_form(x,P0) + q0*x + r0 )
    lam1: 0.5*quad_form(x,P1) + q1*x + r1 <= 0;
cvx_end

obj1 = cvx_optval;
P_lam = P0 + lam1*P1;
q_lam = q0 + lam1*q1;
r_lam = r0 + lam1*r1;
obj2 = -0.5*q_lam*inv(P_lam)*q_lam' + r_lam;

fprintf(1,'Done! \n');

% Displaying results
disp('------------------------------------------------------------------------');
disp('The duality gap is equal to ');
disp(obj1-obj2)


x