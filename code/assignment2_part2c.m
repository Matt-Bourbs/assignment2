% Assignment-2
% Matthieu Bourbeau
% 100975211

%% PART 2C

% Using the Finite Difference Method to solve for the current flow in a 
% rectangular region of L × W.
% Investigating narrowing the "bottle-neck".

% Created a G matrix and an operation matrix to solve for x (G/Op).
% Created loop that incrementally increases bottle-neck.
% Created sigma matrix, assigned values, and defined the box's dimensions
% with a narrowed bottle-neck.

for bottleneck = 0.1:0.01:0.9 
    nx = 50;            
    ny = nx*3/2;        
    G = sparse(nx*ny);  
    Op = zeros(1, nx*ny);
    
    Sigmatrix = zeros(ny, nx); 
    Sig1 = 10^-2;              
    Sig2 = 1;                   
    
    box = [nx*2/5 nx*3/5 ny*bottleneck ny*(1-bottleneck)]; 
    
    for i = 1: nx        
        for j = 1: ny           
            n = j + (i-1)*ny;           
            if i == 1               
                G(n, :) = 0;
                G(n, n) = 1;
                Op(n) = 1;                
            elseif i == nx               
                G(n, :) = 0;
                G(n, n) = 1;
                Op(n) = 0;               
            elseif j == 1                
                if i > box(1) && i < box(2)                  
                    G(n, n) = -3;
                    G(n, n+1) = Sig1;
                    G(n, n+ny) = Sig1;
                    G(n, n-ny) = Sig1;                   
                else                   
                    G(n, n) = -3;
                    G(n, n+1) = Sig2;
                    G(n, n+ny) = Sig2;
                    G(n, n-ny) = Sig2;                  
                end               
            elseif j == ny                
                if i > box(1) && i < box(2)
                    G(n, n) = -3;
                    G(n, n+1) = Sig1;
                    G(n, n+ny) = Sig1;
                    G(n, n-ny) = Sig1;                    
                else                    
                    G(n, n) = -3;
                    G(n, n+1) = Sig2;
                    G(n, n+ny) = Sig2;
                    G(n, n-ny) = Sig2;                    
                end                
            else                
                if i > box(1) && i < box(2) && (j < box(3)||j > box(4))                 
                    G(n, n) = -4;
                    G(n, n+1) = Sig1;
                    G(n, n-1) = Sig1;
                    G(n, n+ny) = Sig1;
                    G(n, n-ny) = Sig1;                
                else                  
                    G(n, n) = -4;
                    G(n, n+1) = Sig2;
                    G(n, n-1) = Sig2;
                    G(n, n+ny) = Sig2;
                    G(n, n-ny) = Sig2;                   
                end
            end
        end
    end
       
    for Length = 1: nx   
        for Width = 1: ny           
            if Length >= box(1) && Length <= box(2)
                Sigmatrix(Width, Length) = Sig1;              
            else
                Sigmatrix(Width, Length) = Sig2;                
            end           
            if Length >= box(1) && Length <= box(2) && Width >= box(3) && Width <= box(4)               
                Sigmatrix(Width, Length) = Sig2;
            end
        end
    end
          
    Voltage = G\Op';   
    sol = zeros(ny, nx, 1);
    
    for i = 1: nx       
        for j = 1: ny           
            n = j+(i-1)*ny;
            sol(j,i) = Voltage(n);            
        end
    end
        
    [elecx, elecy] = gradient(sol);

    J_x = Sigmatrix.*elecx;
    J_y = Sigmatrix.*elecy;
    J = sqrt(J_x.^2 + J_y.^2);
        
    figure(1)
    hold on   
    if bottleneck == 0.1      
        Curr = sum(J, 2);
        Currtot = sum(Curr);
        Currold = Currtot;
        plot([bottleneck, bottleneck], [Currold, Currtot])      
    end    
    if bottleneck > 0.1       
        Currold = Currtot;
        Curr = sum(J, 2);
        Currtot = sum(Curr);
        plot([bottleneck-0.01, bottleneck], [Currold, Currtot])
        xlabel("Bottleneck");
        ylabel("Current Density");     
    end 
    title("Plot of Effect of Narrowing Bottle-neck on Current Density")

end