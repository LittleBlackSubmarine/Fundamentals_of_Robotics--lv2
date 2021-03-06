function W = joint2cartesian(robot,q)

W = [];

[m, n] = size(q);

for k = 1:n
    robot = dirkin(robot,q(:,k));
    W = [W robot.L(robot.n).T(1:3,4)];
end
