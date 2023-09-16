clc;
close all;
clear all;
xt = [0,0,0]';
z = [100,100,100]';
y = [100,100,100]';

N=100000;
X = 50*randn(3,N)

x = [];
l = 1.01;
u = 0.0;
n = 0.99;
for i=1:length(X)
    xx = X(:,i);
%     if (y - l*(z-xx))'*(y - u*(z-xx))<0
%         x=[x,xx];
%     end
    if (z-xt)'*(z - xx) > n*sqrt((z-xt)'*(z-xt)*(z - xx)'*(z - xx))
        x=[x,xx];
    end
end

figure()
plot3(x(1,:),x(2,:),x(3,:),'ro');hold on
plot3(z(1),z(2),z(3),'go')
xlabel('x')
ylabel('y')
zlabel('z')


% xt = [100,0,0]'
% (y - l*(z-xt))'*(y - u*(z-xt))