% Assignment-2
% Matthieu Bourbeau
% 100975211

%% PART 2D

% Using the Finite Difference Method to solve for the current flow in a 
% rectangular region of L × W.
% Investigating varying the σ of the box.

% Created a G matrix and an operation matrix to solve for x (G/Op).
% Created loop that incrementally increases σ.
% Created sigma matrix, assigned values, and defined the box's dimensions
% with a varying σ.

for sigma = 1e-2:1e-2:0.9
    nx = 50;            
    ny = nx*3/2;        
    G = sparse(nx*ny);  
    Op = zeros(1, nx*ny);
    
    Sigmatrix = zeros(ny, nx);           
    Sig1 = 1;                              
    Sig2 = sigma;                         
   
    box = [nx*2/5 nx*3/5 ny*2/5 ny*3/5]; 
         
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
          
    for Length = 1: nx  
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
    sol = zeros(ny, nx, 1);
    
    for x = 1: nx       
        for y = 1: ny           
            n = y+(x-1)*ny;           
            sol(y,x) = Voltage(n);            
        end
    end
               
    [elecx, elecy] = gradient(sol);
                
    J_x = Sigmatrix.*elecx;
    J_y = Sigmatrix.*elecy;
    J = sqrt(J_x.^2 + J_y.^2);
             
    figure(1)
    hold on
    if sigma == 0.01
        Curr = sum(J, 2);
        Currtot = sum(Curr);
        Currold = Currtot;
        plot([sigma, sigma], [Currold, Currtot])
    end
    if sigma > 0.01
        Currold = Currtot;
        Curr = sum(J, 2);
        Currtot = sum(Curr);
        plot([sigma-0.01, sigma], [Currold, Currtot])
        xlabel("Sigma")
        ylabel("Current Density")
    end
    title("Plot of Effect of Varying Sigma on Current Density")
    
end