% Assignment-2
% Matthieu Bourbeau
% 100975211

%% PART 1B

% Using the Finite Difference Method to solve for the electrostatic potential
% in a rectangular region of L × W.
% Solving for the case where V = Vo at x = 0, x = L and V = 0 at y = 0, y = W.

nx = 50;
ny = (3/2)*50;

% Created a G matrix and an operation matrix to solve for x (G/Op).
% Created a loop to apply G matrix's bulk nodes and boundary conditions.
% Solved for numerical solution, along with analytical solution while 
% iterating to create a summation of the infinite series.

G = sparse(nx*ny,nx*ny);
Op = sparse(nx*ny,1);

for x = 1: nx
    for y = 1: ny
        n = y + (x-1)*ny;       
        if x == 1
            G(n, :) = 0;
            G(n, n) = 1;
            Op(n) = 1;
        elseif x == nx
            G(n, :) = 0;
            G(n, n) = 1;
            Op(n) = 1;
        elseif y == 1
            G(n, :) = 0;
            G(n, n) = 1;
        elseif y == ny
            G(n, :) = 0;
            G(n, n) = 1;
        else
            G(n, n) = -4;
            G(n, n+1) = 1;
            G(n, n-1) = 1;
            G(n, n+ny) = 1;
            G(n, n-ny) = 1;
        end
    end
end

Voltage = G\Op;
sol = zeros(nx,ny,1);

for x = 1: nx
    for y = 1: ny       
        n = y+(x-1)*ny;
        sol(x,y) = Voltage(n);
    end  
end   

figure(1)
surf(sol)
axis tight
title("Plot of Surface Voltage Using Numerical Method in 2D")
xlabel("X position")
ylabel("Y position")
zlabel("Voltage")

a = ny;
b = nx/2;

x2 = linspace(-nx/2,nx/2, 50);
y2 = linspace(0,ny,ny);

[i,j] = meshgrid(x2,y2);
sol2 = sparse(ny,nx);

for n = 1:2:600  
    sol2 = (sol2+(cosh(n*pi*i/a).*sin(n*pi*j/a))./(n*cosh(n*pi*b/a)));
    figure(2)
    surf(x2,y2,(4/pi)*sol2)
    axis tight
    title("Plot of Surface Voltage Using Analytical Method in 2D")
    xlabel("X position")
    ylabel("Y position")
    zlabel("Voltage")
    view(-130,30);
    pause(0.001)
end