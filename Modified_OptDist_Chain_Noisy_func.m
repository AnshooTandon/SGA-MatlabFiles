function [sol, fval] = Modified_OptDist_Chain_Noisy_func(theta,q,init_vec)

x = optimvar('x',15,'LowerBound',0);

y = zeros(16,1);
y(1) = 0.5*((1-theta)^3);
y(2) = 0.5*theta*((1-theta)^2);
y(3) = 0.5*(theta^2)*(1-theta);
y(4) = 0.5*theta*((1-theta)^2);
y(5) = 0.5*(theta^2)*(1-theta);
y(6) = 0.5*(theta^3);
y(7) = 0.5*(theta^2)*(1-theta);
y(8) = 0.5*theta*((1-theta)^2);
for ii = 9:16
    y(ii) = y(17-ii);
end
y_orig = y./sum(y);
y = zeros(16,1);
y(1:2:15) = (1-q)*y_orig(1:2:15) + q*y_orig(2:2:16);
y(2:2:16) = q*y_orig(1:2:15) + (1-q)*y_orig(2:2:16);
y = y./sum(y);

obj = fcn2optimexpr(@MYobjfunx,x,y);
prob = optimproblem('Objective',obj);

prob.Constraints.cons1 = sum(x) - 1 <= 0;
rho = 1 - 2*theta;
thresh = (1+rho*rho)/2;
z2 = fcn2optimexpr(@MYconstrnt_func2,x);
prob.Constraints.cons2 = z2 >= thresh;
prob.Constraints.con3 = z2 <= 2;

x0.x = init_vec(1:15);

options = optimoptions('fmincon','ConstraintTolerance',1e-20,...
    'StepTolerance',1e-22,'OptimalityTolerance',1e-15,...
    'MaxIterations',1e6,'MaxFunctionEvaluations',2e7);
%show(prob);
[sol, fval] = solve(prob,x0,'Options',options);
end

function f = MYobjfunx(x,y)
u = [x; 1-sum(x)];
v = y./sum(y);
f = myKLdiv(u,v);
end

function d = myKLdiv(p,q)
    ind = find(p~=0);
    d = sum(p(ind).*log2(p(ind)./q(ind)));
end

function z = MYconstrnt_func2(x)
    x16 = 1 - sum(x);
    
    t12_pos = x(1)+x(2)+x(3)+x(4)+x(13)+x(14)+x(15)+x16;
    t12_neg = x(5)+x(6)+x(7)+x(8)+x(9)+x(10)+x(11)+x(12);
    rho12 = t12_pos - t12_neg;
    
    t34_pos = x(1)+x(4)+x(5)+x(8)+x(9)+x(12)+x(13)+x16;
    t34_neg = x(2)+x(3)+x(6)+x(7)+x(10)+x(11)+x(14)+x(15);
    rho34 = t34_pos - t34_neg;
    
    t13_pos = x(1)+x(2)+x(5)+x(6)+x(11)+x(12)+x(15)+x16;
    t13_neg = x(3)+x(4)+x(7)+x(8)+x(9)+x(10)+x(13)+x(14);
    rho13 = t13_pos - t13_neg;
    
    t24_pos = x(1)+x(3)+x(6)+x(8)+x(9)+x(11)+x(14)+x16;
    t24_neg = x(2)+x(4)+x(5)+x(7)+x(10)+x(12)+x(13)+x(15);
    rho24 = t24_pos - t24_neg;
        
    t14_pos = x(1)+x(3)+x(5)+x(7)+x(10)+x(12)+x(14)+x16;
    t14_neg = x(2)+x(4)+x(6)+x(8)+x(9)+x(11)+x(13)+x(15);
    rho14 = t14_pos - t14_neg;
    
    t23_pos = x(1)+x(2)+x(7)+x(8)+x(9)+x(10)+x(15)+x16;
    t23_neg = x(3)+x(4)+x(5)+x(6)+x(11)+x(12)+x(13)+x(14);
    rho23 = t23_pos - t23_neg;
    
    z = sqrt(rho13*rho24*rho14*rho23)/(rho12*rho34);
end