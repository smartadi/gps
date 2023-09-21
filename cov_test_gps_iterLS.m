clear all;
close all;
clc;

mu = 1;
a  = 1;
c = 1000;

% S = [10000, 0, 0;
%     0, 10000, 0;
%     0, 0, 10000;
%     10000, 10000, 10000];

% S = [500, 9000, 10000;
%     10000, 9000, 500;
%     9000, 500, 10000;
%     10000, 10000, 10000];

% S = [10000, 1000, 0;
%     1000, 10000, 1000;
%     11000, 500, 500;
%     500, 11000, 500;
%     10000, 500, 1500];
S = [10000, 1000, 0;
    1000, 10000, 1000;
    11000, 500, 500;
    500, 11000, 500];



xt = [1000,1000,1000];

t = vecnorm(S-xt,2,2)/c + [-0.01;0.01;-0.01;0.01];

x0 = [900;
    800;
    900];



rho = t*c;



for i=1:10
S-x0';
dp = vecnorm(S - x0',2,2) - c*t 

% dp2 = vecnorm(S2 - x02',2,2) - c*t2; 

% A = [(x0(1)-S(1,1))/rho(1),(x0(2)-S(1,2))/rho(1),(x0(3)-S(1,3))/rho(1),c;
%      (x0(1)-S(2,1))/rho(2),(x0(2)-S(2,2))/rho(2),(x0(3)-S(2,3))/rho(2),c;
%      (x0(1)-S(3,1))/rho(3),(x0(2)-S(3,2))/rho(3),(x0(3)-S(3,3))/rho(3),c;
%      (x0(1)-S(4,1))/rho(4),(x0(2)-S(4,2))/rho(4),(x0(3)-S(4,3))/rho(4),c;
%      (x0(1)-S(5,1))/rho(5),(x0(2)-S(5,2))/rho(5),(x0(3)-S(5,3))/rho(5),c];
 
 A = [(x0(1)-S(1,1))/rho(1),(x0(2)-S(1,2))/rho(1),(x0(3)-S(1,3))/rho(1),c;
     (x0(1)-S(2,1))/rho(2),(x0(2)-S(2,2))/rho(2),(x0(3)-S(2,3))/rho(2),c;
     (x0(1)-S(3,1))/rho(3),(x0(2)-S(3,2))/rho(3),(x0(3)-S(3,3))/rho(3),c;
     (x0(1)-S(4,1))/rho(4),(x0(2)-S(4,2))/rho(4),(x0(3)-S(4,3))/rho(4),c];

% A2 = [(x0(1)-S2(1,1))/rho2(1),(x0(2)-S2(1,2))/rho2(1),(x0(3)-S2(1,3))/rho2(1),c;
%      (x0(1)-S2(2,1))/rho2(2),(x0(2)-S2(2,2))/rho2(2),(x0(3)-S2(2,3))/rho2(2),c;
%      (x0(1)-S2(3,1))/rho2(3),(x0(2)-S2(3,2))/rho2(3),(x0(3)-S2(3,3))/rho2(3),c;
%      (x0(1)-S2(4,1))/rho2(4),(x0(2)-S2(4,2))/rho2(4),(x0(3)-S2(4,3))/rho2(4),c;
%      (x0(1)-S2(5,1))/rho2(5),(x0(2)-S2(5,2))/rho2(5),(x0(3)-S2(5,3))/rho2(5),c;
%      (x0(1)-S2(6,1))/rho2(6),(x0(2)-S2(6,2))/rho2(6),(x0(3)-S2(6,3))/rho2(6),c;
%      (x0(1)-S2(7,1))/rho2(7),(x0(2)-S2(7,2))/rho2(7),(x0(3)-S2(8,3))/rho2(7),c;
%      (x0(1)-S2(8,1))/rho2(8),(x0(2)-S2(8,2))/rho2(8),(x0(3)-S2(8,3))/rho2(8),c];
 
 
dx = inv(A'*A)*A'*dp
% dx2 = inv(A2'*A2)*A2'*dp2;

x0 = x0-dx(1:3);
% x02 = x02-dx2(1:3);

end
x0
C = inv(A'*A);
eig(C)
sum(eig(C))



% C2 = inv(A2'*A2);
% eigs(C2)
% sum(eig(C2))


p = 0.5 - rand(1000,3);
p = p./vecnorm(p,2,2);



m=[];
M=[];
d=[];
R=[];
Q = [eye(3),zeros(3,1);
    zeros(1,4)];
Q = eye(1);
for i=1:length(p)
    B = [A;p(i,:),0];
    m = [m;log(det(inv(Q*B'*B)))];
    M = [M;(det(inv(Q*B'*B)))];
    
    R = [R;(det(inv(A'*A + [p(i,:)';0]*[p(i,:),0])))];

    
end

[MM1,im] = min(m);

[Mm1,iM] = max(m);

[MM2,Im] = min(M);

[Mm2,IM] = max(M);

[Rm,iR] = min(M);

[Rm,IR] = max(M);



p(Im,:)
p(IM,:)

p(iR,:)
p(IR,:)



%%
ss  = 2*S./vecnorm(S,2,2);


%%
close all;
z = zeros(4,1);
zz = zeros(1000,1);
figure()
% plot3(ss(:,1),ss(:,2),ss(:,3),'bo','Linewidth');hold on;
quiver3(z,z,z,ss(:,1),ss(:,2),ss(:,3));hold on
quiver3(zz,zz,zz,p(:,1),p(:,2),p(:,3))
quiver3(0,0,0,2*p(Im,1),2*p(Im,2),2*p(Im,3))
quiver3(0,0,0,2*p(IM,1),2*p(IM,2),2*p(IM,3))

xlabel('x')
ylabel('y')
zlabel('z')


