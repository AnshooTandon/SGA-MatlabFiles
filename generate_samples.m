function [rho12, rho13, rho14, rho23, rho24, rho34] = ...
    generate_samples(num_samples, theta)

% Generate samples for a 4-node Markov chain

d = 4;
x = zeros(d, num_samples);

x(1,:) = (rand(1, num_samples) > 0.5);

tmp_0 = (rand(1, num_samples) > 1-theta);
tmp_1 = (rand(1, num_samples) > theta);
I_0 = find(x(1,:)==0);
x(2,I_0) = tmp_0(I_0);
I_1 = find(x(1,:)~=0);
x(2,I_1) = tmp_1(I_1);

tmp_0 = (rand(1, num_samples) > 1-theta);
tmp_1 = (rand(1, num_samples) > theta);
I_0 = find(x(2,:)==0);
x(3,I_0) = tmp_0(I_0);
I_1 = find(x(2,:)~=0);
x(3,I_1) = tmp_1(I_1);

tmp_0 = (rand(1, num_samples) > 1-theta);
tmp_1 = (rand(1, num_samples) > theta);
I_0 = find(x(3,:)==0);
x(4,I_0) = tmp_0(I_0);
I_1 = find(x(3,:)~=0);
x(4,I_1) = tmp_1(I_1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

jt_12 = [sum(((x(1,:) == 0) & (x(2,:) == 0))) sum(((x(1,:) == 0) & (x(2,:) == 1))) ; 
    sum(((x(1,:) == 1) & (x(2,:) == 0))) sum(((x(1,:) == 1) & (x(2,:) == 1)))]/num_samples;

jt_13 = [sum(((x(1,:) == 0) & (x(3,:) == 0))) sum(((x(1,:) == 0) & (x(3,:) == 1))) ; 
    sum(((x(1,:) == 1) & (x(3,:) == 0))) sum(((x(1,:) == 1) & (x(3,:) == 1)))]/num_samples;

jt_14 = [sum(((x(1,:) == 0) & (x(4,:) == 0))) sum(((x(1,:) == 0) & (x(4,:) == 1))) ; 
    sum(((x(1,:) == 1) & (x(4,:) == 0))) sum(((x(1,:) == 1) & (x(4,:) == 1)))]/num_samples;

jt_23 = [sum(((x(2,:) == 0) & (x(3,:) == 0))) sum(((x(2,:) == 0) & (x(3,:) == 1))) ; 
    sum(((x(2,:) == 1) & (x(3,:) == 0))) sum(((x(2,:) == 1) & (x(3,:) == 1)))]/num_samples;

jt_24 = [sum(((x(2,:) == 0) & (x(4,:) == 0))) sum(((x(2,:) == 0) & (x(4,:) == 1))) ; 
    sum(((x(2,:) == 1) & (x(4,:) == 0))) sum(((x(2,:) == 1) & (x(4,:) == 1)))]/num_samples;

jt_34 = [sum(((x(3,:) == 0) & (x(4,:) == 0))) sum(((x(3,:) == 0) & (x(4,:) == 1))) ; 
    sum(((x(3,:) == 1) & (x(4,:) == 0))) sum(((x(3,:) == 1) & (x(4,:) == 1)))]/num_samples;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho12 = jt_12(1,1) + jt_12(2,2) - jt_12(1,2) - jt_12(2,1);
rho13 = jt_13(1,1) + jt_13(2,2) - jt_13(1,2) - jt_13(2,1);
rho14 = jt_14(1,1) + jt_14(2,2) - jt_14(1,2) - jt_14(2,1);
rho23 = jt_23(1,1) + jt_23(2,2) - jt_23(1,2) - jt_23(2,1);
rho24 = jt_24(1,1) + jt_24(2,2) - jt_24(1,2) - jt_24(2,1);
rho34 = jt_34(1,1) + jt_34(2,2) - jt_34(1,2) - jt_34(2,1);

end


