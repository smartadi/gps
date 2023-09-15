function [errmagnification] = example2(delta_t, phis, thetas)

% 'true' position of observer
true_x = 0;
true_y = 0;
true_z = 6370;
true_d = 0.0001;

n = 20;             % number of iterations of multivariate newton's method
c = 299792.458;     % speed of light
rho = 26570;        % distance (in km) of satellites from center of earth
sat = zeros(4, 3);  % matrix of satellite positions (known)
times = zeros(1,4); % vector of observer time delays (fudged by backwards error)
for n=1:4
    sat(n,1) = rho * cos(phis(n)) * cos(thetas(n)); % sat #n’s x coord
    sat(n,2) = rho * cos(phis(n)) * sin(thetas(n)); % sat #n’s y coord
    sat(n,3) = rho * sin(phis(n));                  % sat #n’s z coord
    dist = sqrt(sat(n,1)^2 + sat(n,2)^2 + (sat(n,3) - 6370)^2);
    times(n) = dist/c + true_d;
end
%sat   % output calculated sat coordinates
%times % output calculated time delays

maxerr_disp = []; % input error vector in the maximum error case
maxerr_vect = []; % calculated observer position in the maximum error case
maxerr = 0;       % maximum position error (in individual coordinate)

%try all possible combinations of error, (-1, -1, -1, -1) to (1, 1, 1, 1)
for t1_err=-1:1
    for t2_err=-1:1
        for t3_err=-1:1
            for t4_err=-1:1
                %calculate fudged time delays recorded by observer
                t1 = times(1) + t1_err * delta_t;
                t2 = times(2) + t2_err * delta_t;
                t3 = times(3) + t3_err * delta_t;
                t4 = times(4) + t4_err * delta_t;
                
                % init_vec is close to actual value but off by 10km initially
                init_vec = [10; 10; 6380; 0];
                
                % multivariate newton's method copied from part 1
                for i = 1:n
                    x = init_vec(1,1);
                    y = init_vec(2,1);
                    z = init_vec(3,1);
                    d = init_vec(4,1);
                    
                    f1 = (x-sat(1,1))^2 + (y-sat(1,2))^2 + (z-sat(1,3))^2 - (c*(t1-d))^2;
                    f2 = (x-sat(2,1))^2 + (y-sat(2,2))^2 + (z-sat(2,3))^2 - (c*(t2-d))^2;
                    f3 = (x-sat(3,1))^2 + (y-sat(3,2))^2 + (z-sat(3,3))^2 - (c*(t3-d))^2;
                    f4 = (x-sat(4,1))^2 + (y-sat(4,2))^2 + (z-sat(4,3))^2 - (c*(t4-d))^2;
                    
                    f_sol = [f1; f2; f3; f4];
                    
                    Jacobian = 2*[x - sat(1,1), y - sat(1,2), z - sat(1,3), c^2*(t1 - d);...
                                  x - sat(2,1), y - sat(2,2), z - sat(2,3), c^2*(t2 - d);...
                                  x - sat(3,1), y - sat(3,2), z - sat(3,3), c^2*(t3 - d);...
                                  x - sat(4,1), y - sat(4,2), z - sat(4,3), c^2*(t4 - d);];
                    
                    s = Jacobian \ (-f_sol);
                    
                    init_vec = init_vec + s;
                end
                
                err = max(abs([init_vec(1)-true_x, init_vec(2)-true_y, init_vec(3)-true_z]));
                if err > maxerr
                    maxerr_disp = [t1_err, t2_err, t3_err, t4_err] * delta_t;
                    maxerr_vect = init_vec;
                    maxerr = err;
                end
            end
        end
    end
end

scatter3([sat(1,1); sat(2,1); sat(3,1); sat(4,1); true_x; maxerr_vect(1)], [sat(1,2); sat(2,2); sat(3,2); sat(4,2); true_y; maxerr_vect(2)], [sat(1,3); sat(2,3); sat(3,3); sat(4,3); true_z; maxerr_vect(3)], [20 ; 20 ; 20 ; 20; 80 ; 80])
xlabel('X');
ylabel('Y');
zlabel('Z');
[globex, globey, globez] = sphere(16);
figure
scatter3(sat(:,1), sat(:,2), sat(:,3))
hold on
surf(6370*globex, 6370*globey, 6370*globez, ones(size(globez)));
axis equal
maxerr_disp;
maxerr
errmagnification = maxerr / (c * delta_t);

