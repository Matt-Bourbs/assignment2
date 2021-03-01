% Assignment-2
% Matthieu Bourbeau
% 100975211

%% PART 1A
% Using the Finite Difference Method to solve for the electrostatic potential
% in a rectangular region of L × W.
% Solving for the simple case where V = Vo at x = 0 and V = 0 at x = L.

nx = 50;
ny = (3/2)*50;

% Created a G matrix and an operation matrix to solve for x (G/Op).
% Created a loop to apply G matrix's bulk nodes and boundary conditions.

G = sparse(nx*ny,nx*ny);
Op = zeros(nx*ny,1);

for x = 1: nx
    for y = 1: ny       
        n = y+(x-1)*ny;       
        if x == 1   
            G(n, :) = 0;
            G(n, n) = 1;
            Op(n) = 1;      
        elseif x == nx        
            G(n, :) = 0;
            G(n, n) = 1;
            Op(n) = 0;           
        elseif y == 1          
            G(n, :) = 0;
            G(n, n) = -3;
            G(n, n+1) = 1;
            G(n, n+ny) = 1;
            G(n, n-ny) = 1;     
        elseif y == ny        
            G(n, n) = -3;
            G(n, n-1) = 1;
            G(n, n+ny) = 1;
            G(n, n-ny) = 1;            
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
        n = y + (x-1)*ny;
        sol(x,y) = Voltage(n);
    end  
end   

figure(1)
surf(sol)
title("Plot of Voltage Using Finite Difference Method in 1D")
xlabel("X position")
ylabel("Y position")
zlabel("Voltage")
view(-130,30)