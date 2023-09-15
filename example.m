format long
c = 299792.458;

sat_1 = [15600; 7540; 20140];
sat_2 = [18760; 2750; 18610];
sat_3 = [17610; 14630; 13480];
sat_4 = [19170; 610; 18390];

t1 = 0.07074;
t2 = 0.07220;
t3 = 0.07690;
t4 = 0.07242;

init_vec_far = [1000; 1000; 50000; 0];

n = 1000;

for i = 1:n

    x = init_vec_far(1,1);
    y = init_vec_far(2,1);
    z = init_vec_far(3,1);
    d = init_vec_far(4,1);

    f1 = (x - sat_1(1,1))^2 + (y - sat_1(2,1))^2 + (z - sat_1(3,1))^2 - (c*(t1 - d))^2;
    f2 = (x - sat_2(1,1))^2 + (y - sat_2(2,1))^2 + (z - sat_2(3,1))^2 - (c*(t2 - d))^2;
    f3 = (x - sat_3(1,1))^2 + (y - sat_3(2,1))^2 + (z - sat_3(3,1))^2 - (c*(t3 - d))^2;
    f4 = (x - sat_4(1,1))^2 + (y - sat_4(2,1))^2 + (z - sat_4(3,1))^2 - (c*(t4 - d))^2;

    f_sol = [f1; f2; f3; f4];

    Jacobian = [2*(x - sat_1(1,1)), 2*(y - sat_1(2,1)), 2*(z - sat_1(3,1)), 2*c^2*(t1 - d);...
                2*(x - sat_2(1,1)), 2*(y - sat_2(2,1)), 2*(z - sat_2(3,1)), 2*c^2*(t2 - d);...
                2*(x - sat_3(1,1)), 2*(y - sat_3(2,1)), 2*(z - sat_3(3,1)), 2*c^2*(t3 - d);...
                2*(x - sat_4(1,1)), 2*(y - sat_4(2,1)), 2*(z - sat_4(3,1)), 2*c^2*(t4 - d)];

    s = Jacobian \ (-f_sol);

    init_vec_far = init_vec_far + s;
end

init_vec_far
dist_far = sqrt(init_vec_far(1,1)^2 + init_vec_far(2,1)^2 + init_vec_far(3,1)^2) - 6370

init_vec_near = [0; 0; 6370; 0];

for i = 1:n

    x = init_vec_near(1,1);
    y = init_vec_near(2,1);
    z = init_vec_near(3,1);
    d = init_vec_near(4,1);

    f1 = (x - sat_1(1,1))^2 + (y - sat_1(2,1))^2 + (z - sat_1(3,1))^2 - (c*(t1 - d))^2;
    f2 = (x - sat_2(1,1))^2 + (y - sat_2(2,1))^2 + (z - sat_2(3,1))^2 - (c*(t2 - d))^2;
    f3 = (x - sat_3(1,1))^2 + (y - sat_3(2,1))^2 + (z - sat_3(3,1))^2 - (c*(t3 - d))^2;
    f4 = (x - sat_4(1,1))^2 + (y - sat_4(2,1))^2 + (z - sat_4(3,1))^2 - (c*(t4 - d))^2;

    f_sol = [f1; f2; f3; f4];

    Jacobian = [2*(x - sat_1(1,1)), 2*(y - sat_1(2,1)), 2*(z - sat_1(3,1)), 2*c^2*(t1 - d);...
                2*(x - sat_2(1,1)), 2*(y - sat_2(2,1)), 2*(z - sat_2(3,1)), 2*c^2*(t2 - d);...
                2*(x - sat_3(1,1)), 2*(y - sat_3(2,1)), 2*(z - sat_3(3,1)), 2*c^2*(t3 - d);...
                2*(x - sat_4(1,1)), 2*(y - sat_4(2,1)), 2*(z - sat_4(3,1)), 2*c^2*(t4 - d)];

    s = Jacobian \ (-f_sol);

    init_vec_near = init_vec_near + s;
end

init_vec_near
dist_near = sqrt(init_vec_near(1,1)^2 + init_vec_near(2,1)^2 + init_vec_near(3,1)^2) - 6370

