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
S = [10000, 000, 0;
    1000, 10000, 1000;
    11000, 500, 500;
    500, 11000, 500];

S = [10000, -100, 100;
    10000, 100, -100;
    10000, -100, -100;
    10000, 100, 100];

xt = [1000,1000,1000];

t = vecnorm(S-xt,2,2)/c + [0.0;0.0;0.0;0.0];

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
 
 

 
 
dx = inv(A'*A)*A'*dp
% dx2 = inv(A2'*A2)*A2'*dp2;

x0 = x0-dx(1:3);
% x02 = x02-dx2(1:3);
end


x0
C = inv(A'*A);
eig(C)
sum(eig(C))

e = eig(C)

% C2 = inv(A2'*A2);
% eigs(C2)
% sum(eig(C2))

%%
p = 0.5 - rand(50000,3);
p = [p;0,0.5,0.5];
p = 1*p./vecnorm(p,2,2);


m=[];
M=[];
d=[];
R=[];
s=[];
s2=[];
o =[];
Q = [eye(3),zeros(3,1);
    zeros(1,3),1];
%  Q = eye(1);

T = [eye(3),zeros(3,1)];
P = T*inv(A'*A)*T';
for i=1:length(p)
     
     rh = norm(x0 - 10000*p(i,:)');
     
     A2 = [(x0(1)-S(1,1))/rho(1),(x0(2)-S(1,2))/rho(1),(x0(3)-S(1,3))/rho(1),c;
     (x0(1)-S(2,1))/rho(2),(x0(2)-S(2,2))/rho(2),(x0(3)-S(2,3))/rho(2),c;
     (x0(1)-S(3,1))/rho(3),(x0(2)-S(3,2))/rho(3),(x0(3)-S(3,3))/rho(3),c;
     (x0(1)-S(4,1))/rho(4),(x0(2)-S(4,2))/rho(4),(x0(3)-S(4,3))/rho(4),c
     (x0(1)-10000*p(i,1))/rh,(x0(2)-10000*p(i,2))/rh,(x0(3)-10000*p(i,3))/rh,c];

%     P = inv(A'*A);
%     Q = [p(i,:)';0]*[p(i,:),0];
%     K = P*Q'*inv(Q*P*Q' + eye(4));
%     CC = P-K*(Q*P*Q'+eye(4))*K';
        
    Q = p(i,:)'*p(i,:);
    K = P*Q'*inv(Q*P*Q' + eye(3));
    CC = P-K*(Q*P*Q'+eye(3))*K';
    
     m = [m;det(CC)];
     
     M = [M;trace(inv(T*A2'*A2*T'))];
     
     R = [R;det((T*A2'*A2*T'))];
     
     e = eigs(T*A2'*A2*T');
     
     s = [s;min(e)];
     % m = [m;trace(CC)];
%     m = [m;det(inv(CC))];
%     M = [M;(det((Q*B'*B)))];
%     
%     R = [R;(det(inv(T*A'*A*T' + p(i,:)'*p(i,:))))];
%     
%     ss = [ss;sum(svd(B'*B))];
%     s2 = [s2;sum(svd(T*A'*A*T' + p(i,:)'*p(i,:)))];


    
end

[MM1,im] = min(m);
[Mm1,iM] = max(m);

[MM2,Im] = min(M);
[Mm2,IM] = max(M);
% 
[Rm,iR] = min(R);
[Rm,IR] = max(R);
% 
[Sm,is] = min(s);
[SM,Is] = max(s);
% 
% [Sm2,is2] = min(s2);
% [SM2,Is2] = max(s2);
% 
% [om,io] = min(o);
% [oM,Io] = max(o);
% 

p(Im,:)
% p(IM,:)
% p(im,:)
p(iM,:)
% p(iR,:)
p(IR,:)
% 
% p(is,:)
p(Is,:)
% 
% p(is2,:)
% p(Is2,:)
% 
% p(io,:)
% p(Io,:)
% 
% Cc = (A'*A);
% E = eig(Cc)
%%
ss  = 2*S./vecnorm(S,2,2);


%%
close all;
z = zeros(4,1);
zz = zeros(1000,1);
figure()
quiver3(0,0,0,500*p(Im,1),500*p(Im,2),500*p(Im,3));hold on
quiver3(0,0,0,500*p(IM,1),500*p(IM,2),500*p(IM,3));hold on;
quiver3(0,0,0,500*p(im,1),500*p(im,2),500*p(im,3));hold on;
quiver3(0,0,0,500*p(iM,1),500*p(iM,2),500*p(iM,3));hold on;
quiver3(0,0,0,500*p(iR,1),500*p(iR,2),500*p(iR,3));hold on;
quiver3(0,0,0,500*p(IR,1),500*p(IR,2),500*p(IR,3));hold on;
quiver3(0,0,0,500*p(Is,1),500*p(Is,2),500*p(Is,3));hold on;
quiver3(0,0,0,500*p(is,1),500*p(is,2),500*p(is,3));hold on;
quiver3(z,z,z,200*ss(:,1),200*ss(:,2),200*ss(:,3));hold on;
ellipsoid(0,0,0,e(1),e(2),e(3));hold on;
% ellipsoid(0,0,0,500,500,500,'FaceAlpha',.3,'EdgeColor','none');hold on;

% quiver3(zz,zz,zz,0.5*p(:,1),0.5*p(:,2),0.5*p(:,3));hold on;
% legend('det max','det min''log det max','log det min','s')
axis([-500 500 -500 500 -500 500])
xlabel('x')
ylabel('y')
zlabel('z')


