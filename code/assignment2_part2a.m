% Assignment-2
% Matthieu Bourbeau
% 100975211

%% PART 2A

% Using the Finite Difference Method to solve for the current flow in a 
% rectangular region of L × W.
% Calculating the current flow at the two contacts and generating plots 
% of σ(x,y), V(x,y), E(x,y), and J⃗(x,y).

nx = 50;                
ny = (3/2)*nx; 

% Created a G matrix and an operation matrix to solve for x (G/Op).
% Created sigma matrix, assigned values, and defined the box's dimensions.
% Created a loop to apply G matrix's bulk nodes, boundary conditions, and
% corresponding bottleneck conditions.

G = sparse(nx*ny);      
Op = zeros(1, nx*ny);    

Sigmatrix = zeros(nx, ny);    
Sig1 = 1;                
Sig2 = 10^-2;  

box = [nx*2/5 nx*3/5 ny*2/5 ny*3/5]; 

for i = 1: nx 
    for j = 1: ny
        if i > box(1) && i < box(2) && (j < box(3)||j > box(4))
            Sigmatrix(i, j) = Sig2;    
        else
            Sigmatrix(i, j) = Sig1;
        end
    end
end

for x = 1: nx
    for y = 1: ny
        n = y+(x-1)*ny;
        nposx = y+(x+1-1)*ny;
        nnegx = y+(x-1-1)*ny;
        nposy = y+1+(x-1)*ny;
        nnegy = y-1+(x-1)*ny;     
        if x == 1
            G(n, :) = 0;
            G(n, n) = 1;
            Op(n) = 1;         
        elseif x == nx 
            G(n, :) = 0;
            G(n, n) = 1;
            Op(n) = 0;            
        elseif y == 1
            G(n, nposx) = (Sigmatrix(x+1, y)+Sigmatrix(x,y))/2;
            G(n, nnegx) = (Sigmatrix(x-1, y)+Sigmatrix(x,y))/2;
            G(n, nposy) = (Sigmatrix(x, y+1)+Sigmatrix(x,y))/2;            
            G(n, n) = -(G(n,nposx)+G(n,nnegx)+G(n,nposy));         
        elseif y == ny         
            G(n, nposx) = (Sigmatrix(x+1, y)+Sigmatrix(x,y))/2;
            G(n, nnegx) = (Sigmatrix(x-1, y)+Sigmatrix(x,y))/2;
            G(n, nnegy) = (Sigmatrix(x, y-1)+Sigmatrix(x,y))/2;
            G(n, n) = -(G(n,nposx)+G(n,nnegx)+G(n,nnegy));         
        else           
            G(n, nposx) = (Sigmatrix(x+1, y)+Sigmatrix(x,y))/2;
            G(n, nnegx) = (Sigmatrix(x-1, y)+Sigmatrix(x,y))/2;
            G(n, nposy) = (Sigmatrix(x, y+1)+Sigmatrix(x,y))/2;
            G(n, nnegy) = (Sigmatrix(x, y-1)+Sigmatrix(x,y))/2;
            G(n, n) = -(G(n,nposx)+G(n,nnegx)+G(n,nposy)+G(n,nnegy));           
        end
    end
end

figure(1)
surf(Sigmatrix);
xlabel("X position")
ylabel("Y position")
zlabel("Sima")
axis tight
view([40 30]);
title("Plot of Surface Sigma")

Voltage = G\Op';
sol = zeros(ny, nx, 1);

for i = 1: nx
    for j = 1: ny
        n = j + (i-1)*ny;
        sol(j,i) = Voltage(n);
    end
end

figure(2)
surf(sol)
axis tight
xlabel("X position")
ylabel("Y position")
zlabel("Voltage")
view([40 30]);
title("Plot of Surface Voltage with Bottle-neck Conditions")

[elecx, elecy] = gradient(sol);

figure(3)
surf(-elecx)
axis tight
xlabel("X position")
ylabel("Y position")
zlabel("Electric Field")
view([40 30]);
title("Plot of X-component of Surface Electric Field")

figure(4)
surf(-elecy)
axis tight
xlabel("X position")
ylabel("Y position")
zlabel("Electric Field")
view([40 30]);
title("Plot of Y-component of Surface Electric Field")
J_x = Sigmatrix'.*elecx;
J_y = Sigmatrix'.*elecy;
J = sqrt(J_x.^2 + J_y.^2);

figure(5)
surf(J)
axis tight
xlabel("X position")
ylabel("Y position")
zlabel("Current Density")
view([40 30]);
title("Plot of Surface Current Density")