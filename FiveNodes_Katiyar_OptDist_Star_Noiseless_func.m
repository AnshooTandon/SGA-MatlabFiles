function [sol, fval] = FiveNodes_Katiyar_OptDist_Star_Noiseless_func(theta,init_vec)

x = optimvar('x',31,'LowerBound',0);

y = zeros(32,1);
y(1) = 0.5*((1-theta)^3);
y(2) = 0.5*theta*((1-theta)^2);
y(3) = 0.5*theta*((1-theta)^2);
y(4) = 0.5*(theta^2)*(1-theta);
y(5) = 0.5*theta*((1-theta)^2);
y(6) = 0.5*(theta^2)*(1-theta);
y(7) = 0.5*(theta^2)*(1-theta);
y(8) = 0.5*(theta^3);
for ii = 9:16
    y(ii) = y(17-ii);
end
y = y./sum(y);

y_new = zeros(32,1);
y_new(1:8) = y(1:8)*(1-theta);
y_new(9:16) = y(1:8)*theta;
y_new(17:24) = y(9:16)*theta;
y_new(25:32) = y(9:16)*(1-theta);
y_new = y_new./sum(y_new);

obj = fcn2optimexpr(@MYobjfunx,x,y_new);
prob = optimproblem('Objective',obj);

prob.Constraints.cons1 = sum(x) - 1 <= 0;

rho = 1 - 2*theta;
thresh = (1+rho*rho)/2;

z2 = fcn2optimexpr(@MYconstrnt_func2,x);
prob.Constraints.cons2 = z2 <= thresh;
prob.Constraints.con3 = z2 >= 0.2;

z3 = fcn2optimexpr(@MYconstrnt_func3,x);
prob.Constraints.con4 = z3 >= thresh;
prob.Constraints.con5 = z3 <= 2;

x0.x = init_vec(1:31);

options = optimoptions('fmincon','ConstraintTolerance',1e-20,...
    'StepTolerance',1e-22,'OptimalityTolerance',1e-15,...
    'MaxIterations',1e6,'MaxFunctionEvaluations',2e7);
%show(prob);
[sol, fval] = solve(prob,x0,'Options',options);
end

function f = MYobjfunx(u,v)
u = [u; 1-sum(u)];
f = myKLdiv(u,v);
end

function d = myKLdiv(p,q)
    ind = find(p~=0);
    d = sum(p(ind).*log2(p(ind)./q(ind)));
end

function z2 = MYconstrnt_func2(x)
    x32 = 1 - sum(x);
    
    t12_pos = x(1)+x(2)+x(3)+x(4)+x(13)+x(14)+x(15)+x(16)+...
            x(17)+x(18)+x(19)+x(20)+x(29)+x(30)+x(31)+x32;
    t12_neg = x(5)+x(6)+x(7)+x(8)+x(9)+x(10)+x(11)+x(12)+...
            x(21)+x(22)+x(23)+x(24)+x(25)+x(26)+x(27)+x(28);
    rho12 = t12_pos - t12_neg;
    
    t34_pos = x(1)+x(4)+x(5)+x(8)+x(9)+x(12)+x(13)+x(16)+...
            x(17)+x(20)+x(21)+x(24)+x(25)+x(28)+x(29)+x32;
    t34_neg = x(2)+x(3)+x(6)+x(7)+x(10)+x(11)+x(14)+x(15)+...
            x(18)+x(19)+x(22)+x(23)+x(26)+x(27)+x(30)+x(31);
    rho34 = t34_pos - t34_neg;
    
    t13_pos = x(1)+x(2)+x(5)+x(6)+x(11)+x(12)+x(15)+x(16)+...
            x(17)+x(18)+x(21)+x(22)+x(27)+x(28)+x(31)+x32;
    t13_neg = x(3)+x(4)+x(7)+x(8)+x(9)+x(10)+x(13)+x(14)+...
            x(19)+x(20)+x(23)+x(24)+x(25)+x(26)+x(29)+x(30);
    rho13 = t13_pos - t13_neg;
    
    t24_pos = x(1)+x(3)+x(6)+x(8)+x(9)+x(11)+x(14)+x(16)+...
            x(17)+x(19)+x(22)+x(24)+x(25)+x(27)+x(30)+x32;
    t24_neg = x(2)+x(4)+x(5)+x(7)+x(10)+x(12)+x(13)+x(15)+...
            x(18)+x(20)+x(21)+x(23)+x(26)+x(28)+x(29)+x(31);
    rho24 = t24_pos - t24_neg;
    
    z2 = (rho13*rho24)/(rho12*rho34);
end


function z3 = MYconstrnt_func3(x)
    x32 = 1 - sum(x);
    
     t13_pos = x(1)+x(2)+x(5)+x(6)+x(11)+x(12)+x(15)+x(16)+...
            x(17)+x(18)+x(21)+x(22)+x(27)+x(28)+x(31)+x32;
    t13_neg = x(3)+x(4)+x(7)+x(8)+x(9)+x(10)+x(13)+x(14)+...
            x(19)+x(20)+x(23)+x(24)+x(25)+x(26)+x(29)+x(30);
    rho13 = t13_pos - t13_neg;
    
    t24_pos = x(1)+x(3)+x(6)+x(8)+x(9)+x(11)+x(14)+x(16)+...
            x(17)+x(19)+x(22)+x(24)+x(25)+x(27)+x(30)+x32;
    t24_neg = x(2)+x(4)+x(5)+x(7)+x(10)+x(12)+x(13)+x(15)+...
            x(18)+x(20)+x(21)+x(23)+x(26)+x(28)+x(29)+x(31);
    rho24 = t24_pos - t24_neg;
        
    t14_pos = x(1)+x(3)+x(5)+x(7)+x(10)+x(12)+x(14)+x(16)+...
            x(17)+x(19)+x(21)+x(23)+x(26)+x(28)+x(30)+x32;
    t14_neg = x(2)+x(4)+x(6)+x(8)+x(9)+x(11)+x(13)+x(15)+...
            x(18)+x(20)+x(22)+x(24)+x(25)+x(27)+x(29)+x(31);
    rho14 = t14_pos - t14_neg;
    
    t23_pos = x(1)+x(2)+x(7)+x(8)+x(9)+x(10)+x(15)+x(16)+...
            x(17)+x(18)+x(23)+x(24)+x(25)+x(26)+x(31)+x32;
    t23_neg = x(3)+x(4)+x(5)+x(6)+x(11)+x(12)+x(13)+x(14)+...
            x(19)+x(20)+x(21)+x(22)+x(27)+x(28)+x(29)+x(30);
    rho23 = t23_pos - t23_neg;

    z3 = (rho13*rho24)/(rho14*rho23);
end