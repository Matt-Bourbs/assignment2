% Assignment-2
% Matthieu Bourbeau
% 100975211

%% PART 2B

% Using the Finite Difference Method to solve for the current flow in a 
% rectangular region of L × W.
% Investigating mesh density.

% Created a G matrix and an operation matrix to solve for x (G/Op).
% Created loop that incrementally increased mesh size and multiplied G and 
% operation matrix by the mesh size.
% Created sigma matrix, assigned values, and defined the box's dimensions
% with varying mesh size.

for meshsize = 10:10:100
    ny = (3/2)*meshsize;         
    G = sparse(meshsize*ny);   
    Op = zeros(1, meshsize*ny);    
    
    Sigmatrix = zeros(ny, meshsize);       
    Sig1 = 1;                              
    Sig2 = 10^-2;                               
   
    box = [meshsize*2/5 meshsize*3/5 ny*2/5 ny*3/5];  
     
    for x = 1: meshsize       
        for y = 1: ny            
            n = y+(x-1)*ny;            
            if x == 1
                G(n, :) = 0;
                G(n, n) = 1;
                Op(n) = 1;               
            elseif x == meshsize
                G(n, :) = 0;
                G(n, n) = 1;
                Op(n) = 0;               
            elseif y == 1                
                if x > box(1) && x < box(2)
                    G(n, n) = -3;
                    G(n, n+1) = Sig2;
                    G(n, n+ny) = Sig2;
                    G(n, n-ny) = Sig2;                   
                else                   
                    G(n, n) = -3;
                    G(n, n+1) = Sig1;
                    G(n, n+ny) = Sig1;
                    G(n, n-ny) = Sig1;                    
                end                
            elseif y == ny              
                if x > box(1) && x < box(2)                    
                    G(n, n) = -3;
                    G(n, n+1) = Sig2;
                    G(n, n+ny) = Sig2;
                    G(n, n-ny) = Sig2;                   
                else                    
                    G(n, n) = -3;
                    G(n, n+1) = Sig1;
                    G(n, n+ny) = Sig1;
                    G(n, n-ny) = Sig1;                   
                end               
            else               
                if x > box(1) && x < box(2) && (y < box(3)||y > box(4))               
                    G(n, n) = -4;
                    G(n, n+1) = Sig2;
                    G(n, n-1) = Sig2;
                    G(n, n+ny) = Sig2;
                    G(n, n-ny) = Sig2;                   
                else                    
                    G(n, n) = -4;
                    G(n, n+1) = Sig1;
                    G(n, n-1) = Sig1;
                    G(n, n+ny) = Sig1;
                    G(n, n-ny) = Sig1;                 
                end
            end
        end
    end
    
    for Length = 1: meshsize        
        for Width = 1: ny           
            if Length >= box(1) && Length <= box(2)
                Sigmatrix(Width, Length) = Sig2;              
            else               
                Sigmatrix(Width, Length) = Sig1;               
            end            
            if Length >= box(1) && Length <= box(2) && Width >= box(3) && Width <= box(4)               
                Sigmatrix(Width, Length) = Sig1;           
            end
        end
    end
       
    Voltage = G\Op';   
    sol = zeros(ny, meshsize, 1);
    
    for i = 1: meshsize       
        for j = 1: ny         
            n = j + (i-1)*ny;
            sol(j,i) = Voltage(n);           
        end
    end
    
    [elecx, elecy] = gradient(sol);
    
    J_x = Sigmatrix.*elecx;
    J_y = Sigmatrix.*elecy;
    J = sqrt(J_x.^2 + J_y.^2);
                                        
    figure(1)
    hold on    
    if meshsize == 10        
        Curr = sum(J, 1);                 
        Currtot = sum(Curr);
        Currold = Currtot;
        plot([meshsize, meshsize], [Currold, Currtot])      
    end   
    if meshsize > 10        
        Currold = Currtot;
        Curr = sum(J, 2);
        Currtot = sum(Curr);
        plot([meshsize-10, meshsize], [Currold, Currtot])
        xlabel("Mesh Size")
        ylabel("Current Density")       
    end     
    title("Plot of Effect of Mesh Size on Current Density")

end