% pde assignment 8
n = 81;
m = 9;

A = zeros(n);

for k=1:1:m
    for l=1:1:m
        A(k*l,(k-1)*m+ l) = -1/H;
        A(k*l,(k)*m+ l) = 4/H;
		A(k*l,(k)*m+ l-1) = -1/H;
		A(k*l,(k)*m+ l+1) = -1/H;
		A(k*l,(k+1)*m+ l) = 1/H;
    end
end

u_sol = zeros(1,n);  % solution

% array f

