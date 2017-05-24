%
%  ajay dsouza 
%
% Implementation of Newtons Method
% For solving a system of Non Linear Equations
% with Boundary Value Condition
%

% clear
clear;

% setting seed for any random number generation
rand('seed', 1234567);
 
% change current folder
cd('C:/wk/odrive/Amazon Cloud Drive/ajays_stuff/georgia_tech_ms/math_cse6644_interative_methods_for_system_of_equations/project');

% tolerance for terminating the iteration
tolerance = 1e-8;

%potential
global V
V=@(x,y) (1-x^2)^2/4+(y+x^2-1)^2/2;

%gradient
global gradV
gradV=@(x,y) [-(1-x^2)*x+2*x*(y+x^2-1);y+x^2-1];

%Hessian
global HessV
HessV=@(x,y) [-1+3*x^2+2*(y+x^2-1)+4*x^2 2*x;2*x 1];

%F_x
global F_x
F_x=@(x,y) (27*x^5 - 34*x^3 + 7*x  + ...
    24*x^3*y + 4*x*y^2 - 10*x*y);

%F_y
global F_y
F_y=@(x,y) ( 6*x^4 - 5*x^2 + ...
    4*x^2*y + y - 1);

%dF_1_dx
global dF_1_dx
dF_1_dx=@(x,y) (135*x^4 - 102*x^2 + 7 + 72*x^2*y + 4*y^2-10*y);

%dF_2_dx
global dF_2_dx
dF_2_dx=@(x,y) (24*x^3 - 10*x + 8*x*y);

%dF_1_dy
global dF_1_dy
dF_1_dy=@(x,y) (24*x^3 + 8*x*y - 10*x);

%dF_2_dy
global dF_2_dy
dF_2_dy=@(x,y) (4*x^2 + 1);

%initial value
global X0
X0=[-1;0];

%endvalue
global XT
XT=[1;0];

% The list of T to be used
% T= M*dh => T= n*dh
Tl = [1,5,10,20,100];


% time interval size for each T
% T = M*dh => T= n*dh
dhL = [0.1,0.1,0.1,0.2,0.2];

% color for plotting for each T
colr = ['r','g','b','k','m'];

% save x and y for plotting
Xsave = zeros((max(Tl)/min(dhL))+1,length(Tl));
Ysave = zeros((max(Tl)/min(dhL))+1,length(Tl));

% save performance statistics
perfMetrics = zeros(length(Tl),5);

% run newtons method for each T
for Ti = 1:length(Tl)
    
    % get the T
    T = Tl(Ti);
    
    % get the timestamp for this T
    dh = dhL(Ti);
    
    % no of intervals M, so points is M+1
    % T = M*dh => T= n*dh
    % No o points is M+1 from 0..M or in matlab 1...M+1
    M = T/dh;
  
    
    % Set up for newtons method for each T        
    %Ybar
    Ybar = zeros(2*(M-1),1);
    Ybar(1,1) = -1;
    Ybar(2*(M-2)+1,1) = 1;
    
    % A
    A = zeros(2*(M-1));
    % check if matrix is positive definite - P is 0
    [~,p] = chol(A);
    
    for i=1:(M-1)
        
        idx = (i-1)*2;
        
        A(idx+1,idx+1) = 2;
        A(idx+2,idx+2) = 2;
        
        if ( i > 1 )
            A(idx+1,idx-1) = -1;
            A(idx+2,idx) = -1;
        end
        
        if ( i < (M-1) )
            A(idx+1,idx+3) = -1;
            A(idx+2,idx+4) = -1;
        end
        
    end
    
    
    % X
    % initilized to 0 
    X =  zeros(2*(M+1),1);
    X(1,1) = -1;
    X((2*M)+1,1) = 1;
    %X(3:2*M,1) = rand(2*(M-1),1)-(1/2);
    
    
    %RHS F(X_i)
    F = zeros(2*(M-1),1);
    
    for i=1:(M-1)
        idx = (i-1)*2;
        F(idx+1:idx+2,1) = ...
            HessV(X(idx+2),X(idx+3))*gradV(X(idx+2),X(idx+3));
    end
    
    
    %Jacobian
    % J = A+h^2J(F(X))
    J = zeros(2*(M-1));
    
    for i=1:(M-1)
        
        idx = (i-1)*2;
        
        J(idx+1:idx+2,idx+1:idx+2) = ...
            dh^2*([dF_1_dx(X(idx+2),X(idx+3)),dF_1_dy(X(idx+2),X(idx+3));
            dF_2_dx(X(idx+2),X(idx+3)),dF_2_dy(X(idx+2),X(idx+3))]);
        
        
        J(idx+1,idx+1) = 2 + J(idx+1,idx+1);
        J(idx+2,idx+2) = 2 + J(idx+1,idx+1);
        
        if ( i > 1 )
            J(idx+1,idx-1) = J(idx+1,idx-1)-1;
            J(idx+2,idx) = J(idx+2,idx)-1;
        end
        
        if ( i < (M-1) )
            J(idx+1,idx+3) = J(idx+1,idx+3)-1;
            J(idx+2,idx+4) = J(idx+2,idx+4)-1;
        end
    end
        
    %start time in secs
    startTime = datevec(now);
    
    % Newton Method - Begin Iterations
    
    % get the initial residue to start the loop
    b = A*X(3:2*M) + dh^2*F - Ybar;
    
    iterations = 0;
    while(max(abs(b))>tolerance)
        
        % compute the new x from old x
        % as x(t+1) = x(t) - inv(Jacobian)*b
        X(3:2*M) = X(3:2*M) - inv(J)*b;

        % get the reside b with the new x
        b = A*X(3:2*M) + dh^2*F - Ybar;
        
        iterations = iterations + 1;
        
        % limit on iterations 
        if iterations > 20000
            fprintf('Reached Limit Breaking \n T=%d , Iteration:%d  - norm:%f\n',...
                 T,iterations,norm(b));
            break;
        end
        
        % print progress of the iterations
        if rem(iterations,50) == 0
             fprintf('T=%d , Iteration:%d  - norm:%f\n',...
                 T,iterations,norm(b));
        end
    end
    
    % clear what we do not bneed ad newtons method has exited
    % for this T
    clear J;
    clear A;
    clear F;
    clear Ybar;
    
    %time taken in secs
    endTime = datevec(now);
    
    fprintf('T= %d\n',T);
    fprintf('Time Taken(s) %f\n',etime(endTime,startTime));
    fprintf('Converging Norm %f\n\n\n',norm(b));
    
    % get the X and Y values into two different arrays
    Y=X(2:2:end);
    X(2:2:end)=[];
    
    % save them to a struct whcih keeps teach of X,Y values for
    % all T , will save them to file later for plotting
    Xsave(1:length(X),Ti) = X;
    Ysave(1:length(Y),Ti) = Y;
    
    % save the performance metric fro this T
    % will be saved to file for plotting
    perfMetrics(Ti,:) = [T,iterations,...
       etime(endTime,startTime),norm(b),dh];

    % A trivial plot for checking
    plot(X,Y,'color',colr(Ti));
    hold on;

    % clear the X,Y from newtons as this has been saved and iteration
    % is complete
    clear X;
    clear Y;
    
    %[X,Y] = meshgrid(-1:.1:1);
    %mesh(X,Y,V(X,Y))
      
end


% save the results to a csv file, can be read by R to plot in by ggplot
fname = sprintf('project_1_x.csv');
csvwrite(fname,Xsave);


% save the results to a csv file, can be read by R to plot in by ggplot
fname = sprintf('project_1_y.csv');
csvwrite(fname,Ysave);


% save the results to a csv file, can be read by R to plot in by ggplot
fname = sprintf('project_1_perf.csv');
csvwrite(fname,perfMetrics);


%--------------------------END ------------------------------------
