%% ASSIGNMENT-4  Advance CFD, Kunal Harinkhede (231290402)
%% Solution of Navier Stokes equation using SIMPLE algorithm to solve for flow over a rectangular obstacle.

clear
format long
clc

%% reading geometry and parameter files to store input data

fileID=fopen("history.txt",'w');
fileID=fopen("history.txt",'a');
fprintf(fileID,"solving the given problem using SIMPLE algorithm..\n\n")
fprintf(fileID,"reading geometric parameters of the problem from 'geometryInfo.txt' file..\n")

fID = fopen('problemParameters.txt', 'r');
data = textscan(fID, '%f %f %f %f %f', 'Delimiter', ',');
fclose(fID);

H=data{1};
L=data{2};
La=data{3};
B=data{4};
Re=data{5};


fprintf(fileID,"reading mesh parameters for the problem from 'meshParameters.txt' file..\n");
fID = fopen('solverParameters.txt', 'r');
data = textscan(fID, '%f %f %f %f %f %f %f %f', 'Delimiter', ',');
fclose(fID);

nxMax=data{1};
nyMax=data{2};
tolerance=data{3};
endTime=data{4};
maxTimeStep=data{5};
maxCourantNo=data{6};
alphaU=data{7};
alphaP=data{8};

% discretizing the domain and setting the boundary cells


% this section will create some files that store the grid data i.e. the
% data of cell centres and the vertices of the cells. 
% Vertices of the cells are termed as the "Nodes", whereas the cell 

%creating file for node points

delX=L/(nxMax);
delY=H/(nyMax); % these nxMax and nyMax represents number of cells I have in orthogonal directions

% to introduce ghost cells on all sides of the domain, we'll increase the
% domain by one series of cells on all sides

nyMax=nyMax+2;
nxMax=nxMax+2;

% to find the position of the cuboid obstale and iolate the cells that are
% covered by this obstcle,

nxBodyMin=floor((La-B/2)/delX);
nxBodyMax=ceil((La+B/2)/delX);
nyBodyMin=floor((H-B)/(2*delY));
nyBodyMax=floor((H+B)/(2*delY));

% INITIAlIZATION: let's create the primitive variable matrices

time=0; %starting time

P=zeros(nyMax,nxMax);   % zero gauge pressure
U=ones(nyMax,nxMax);
V=zeros(nyMax,nxMax);

% % removing all the cells that cover the area of the blockage
% for j=nyBodyMin:nyBodyMax
%     for i=nxBodyMin:nxBodyMax
%         P(j,i)=NaN;
%         U(j,i)=NaN;
%         V(j,i)=NaN;
%     end
% end

% creating x and y axes for cell centre values

x=delX*(0.5:1:(nxMax-2.5));
y=delY*(0.5:1:(nyMax-2.5));

fprintf(fileID,"time: %f \n",time);
fprintf(fileID,"standard INITIALIZATION complete !!\n\n");

%% Time loop


while time<endTime
    
    % first calculate time step
    % time step is calculated on the basis of the maximum local courant
    % number 
    maxLocalTimeStepInv=0;
    for j=2:nyMax-1
        for i=2:nxMax-1
            if ((abs(U(j,i))/delX+abs(V(j,i))/delY)/maxCourantNo)>=maxLocalTimeStepInv
                maxLocalTimeStepInv=(abs(U(j,i))/delX+abs(V(j,i))/delY)/maxCourantNo;
            end
        end
    end
    timeStepInv=max(1/maxTimeStep,maxLocalTimeStepInv);
    deltaT=1/timeStepInv;
    time=time+deltaT
    
    % finding maximum courant number based on set current time step:
    maxLocalCourantNo=deltaT*maxLocalTimeStepInv*maxCourantNo;
    
    fprintf(fileID,"time: %f \n",time);
    fprintf(fileID,"deltaT: %f, maximum local Courant number: %f \n\n",[deltaT,maxLocalCourantNo]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % iterative loop
    rmsError=100;
    iter=0;
    % calculation for the coefficiets that are common for one time step and
    % not necessary to be calculated at all iterations
    
    % there will be eight matrices for solving ADI for u equation and same
    % eight matrices for solving v equations for 
    
    Utilda=U;   %taking the values from the previous time-step
    Vtilda=V;
    Ptilda=P;
    
    while rmsError>tolerance    % the iterative loop will keep correcting the field varibles
        % let's call the provisional velocity and pressure fields as Utilda
        % etc.
                
        
        % ALL THE COEFFICIENTS ARE CALUCLATED USING MESH PARAMTERS AND
        % U,V,P DATA FIELDS i.e. THE FIELD FROM PREVIOUS TIME STEP
        

        
        % Correcting boundary condition:
        
        % Top and Bottom wall boundary (No Slip)
        for i=2:nxMax-1
            % for U
            Utilda(1,i)=-Utilda(2,i); % no slip boundary at bottom boundary
            Utilda(nyMax,i)=-Utilda(nyMax-1,i); % no slip boundary at top boundary
            % for V
            Vtilda(1,i)=0; % bottom wall: no penetration formula
            Vtilda(nyMax-1,i)=0; % top wall: no penetration condition
            % for P
            Ptilda(1,i)=Ptilda(2,i); % zero gradient pressure
            Ptilda(nyMax,i)=Ptilda(nyMax-1,i);
        end
        
        % Inlet Boundary
        for j=1:nyMax
            % for U
            Utilda(j,1)=1;  % inlet velocity
            % for V
            Vtilda(j,1)=-Vtilda(j,2);   % zero v component of velocity
            % for P
            Ptilda(j,1)=Ptilda(j,2); % zero gradient pressure
        end
        
        % Outlet boundary
        for j=1:nyMax
            % for U
            Utilda(j,nxMax)=Utilda(j,nxMax-1);  % zero gradient boundary
            % for V
            Vtilda(j,nxMax)=Vtilda(j,nxMax-1);
            % for P
            Ptilda(j,nxMax)=0;  % constant dump pressure
        end
        
        
        % solving Utilda: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        % first the row sweep:
        for j=2:nyMax-1
            
                A=[];   % the coefficient matrix
                for m1=1:nxMax-2
                    for n1=1:nxMax-2
                        m=m1+1;
                        n=n1+1;
                        if m1==n1                            
                            A(m1,n1)=delX*delY/deltaT+0.25*(U(j,m)+U(j,m+1))*delY+...
                                2*(delX/Re/delY+delY/Re/delX)-0.25*(U(j,m-1)+U(j,m))*delY+...
                                0.25*(V(j,m)+V(j,m+1))*delX-0.25*(V(j-1,m)+V(j-1,m+1))*delX;
                        elseif n1==m1+1
                            A(m1,n1)=0.25*delY*(U(j,m)+U(j,m+1))-delY/Re/delX;
                        elseif n1==m1-1
                            A(m1,n1)=-0.25*delY*(U(j,m)+U(j,m-1))-delY/Re/delX;
                        else 
                            A(m1,n1)=0;
                        end
                    end
                end
                
                b=[];   % RHS vector                
                for m1=1:nxMax-2
                    m=m1+1;
                    b(m1)=-Utilda(j+1,m)*(0.25*delX*(V(j,m)+V(j,m+1))-delX/Re/delY)...
                        -Utilda(j-1,m)*(-0.25*delX*(V(j-1,m)+V(j-1,m+1))-delX/Re/delY)...
                        +delX*delY/deltaT*U(j,m)...
                        +(Ptilda(j,m)-Ptilda(j,m+1))*delY;
                    if m1==1
                        b(m1)=b(m1)-Utilda(j,m-1)*(-0.25*delY*(U(j,m-1)+U(j,m))-delY/Re/delX);
                    elseif m1==nxMax-2
                        b(m1)=b(m1)-Utilda(j,m+1)*(0.25*delY*(U(j,m)+U(j,m+1))-delY/Re/delX);
                    end
                end
                
                Utilda(j,2:nxMax-1)=thomasAlgorithm(A,b);  
                
        end
        
        % then column sweeep
        for i=2:nxMax-1
            
            A=[];
            for m1=1:nyMax-2
                    for n1=1:nyMax-2
                        m=m1+1;
                        n=n1+1;
                        if m1==n1                            
                            A(m1,n1)=delX*delY/deltaT+0.25*(U(m,i)+U(m,i+1))*delY+...
                                    2*(delX/Re/delY+delY/Re/delX)-0.25*(U(m,i-1)+U(m,i))*delY+...
                                    0.25*(V(m,i)+V(m,i+1))*delX-0.25*(V(m-1,i)+V(m-1,i+1))*delX;
                        elseif n1==m1+1
                            A(m1,n1)=0.25*delX*(V(m,i)+V(m,i+1))-delX/Re/delY;
                        elseif n1==m1-1
                            A(m1,n1)=-0.25*delX*(V(m-1,i)+V(m-1,i+1))-delX/Re/delY;
                        else 
                            A(m1,n1)=0;
                        end
                    end
            end
             
            b=[];
            for m1=1:nyMax-2
                    m=m1+1;
                    b(m1)=-Utilda(m,i+1)*(0.25*delY*(U(m,i)+U(m,i+1))-delY/Re/delX)...
                        -Utilda(m,i-1)*(-0.25*delY*(U(m,i-1)+U(m,i))-delY/Re/delX)...
                        +delX*delY/deltaT*U(m,i)...
                        +(Ptilda(m,i)-Ptilda(m,i+1))*delY;
                    if m1==1
                        b(m1)=b(m1)-Utilda(m-1,i)*(-0.25*delX*(V(m-1,i)+V(m-1,i+1))-delX/Re/delY);
                    elseif m1==nxMax-2
                        b(m1)=b(m1)-Utilda(m+1,i)*(0.25*delX*(V(m,i)+V(m,i+1))-delX/Re/delY);
                    end
             end
                
             Utilda(2:nyMax-1,i)=thomasAlgorithm(A,b);                            
            
        end
        
        % solving for Vtilda: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % first row sweep:
        for j=2:nyMax-2 %(last rows are not required)
            
            A=[];
            for m1=1:nxMax-2
                    for n1=1:nxMax-2
                        m=m1+1;
                        n=n1+1;
                        if m1==n1                            
                            A(m1,n1)=delX*delY/deltaT+0.25*(U(j,m)+U(j+1,m))*delY+...
                                    2*(delX/Re/delY+delY/Re/delX)-0.25*(U(j,m-1)+U(j+1,m-1))*delY...
                                    +0.25*(V(j,m)+V(j+1,m))*delX-0.25*(V(j-1,m)+V(j,m))*delX;
                        elseif n1==m1+1
                            A(m1,n1)=0.25*delY*(U(j,m)+U(j+1,m))-delY/Re/delX;
                        elseif n1==m1-1
                            A(m1,n1)=-0.25*delY*(U(j,m-1)+U(j+1,m-1))-delY/Re/delX;
                        else 
                            A(m1,n1)=0;
                        end
                    end
            end
            
            b=[];
            for m1=1:nxMax-2
                    m=m1+1;
                    b(m1)=-Vtilda(j+1,m)*(0.25*delX*(V(j,m)+V(j+1,m))-delX/Re/delY)...
                        -Vtilda(j-1,m)*(-0.25*delX*(V(j-1,m)+V(j,m))-delX/Re/delY)...
                        +delX*delY/deltaT*V(j,m)...
                        +(Ptilda(j,m)-Ptilda(j+1,m))*delX;
                    if m1==1
                        b(m1)=b(m1)-Vtilda(j,m-1)*(-0.25*delY*(U(j,m-1)+U(j+1,m-1))-delY/Re/delX);
                    elseif m1==nxMax-2
                        b(m1)=b(m1)-Vtilda(j,m+1)*(0.25*delY*(U(j,m)+U(j+1,m))-delY/Re/delX);
                    end
             end
                
             Vtilda(j,2:nxMax-1)=thomasAlgorithm(A,b);   
                    
        end
        
        % now consider column sweep
        for i=2:nxMax-1
            
            A=[];
            for m1=1:nyMax-3
                    for n1=1:nyMax-3
                        m=m1+1;
                        n=n1+1;
                        if m1==n1                            
                            A(m1,n1)=delX*delY/deltaT+0.25*(U(m,i)+U(m+1,i))*delY+...
                                    2*(delX/Re/delY+delY/Re/delX)-0.25*(U(m,i-1)+U(m+1,i-1))*delY...
                                    +0.25*(V(m,i)+V(m+1,i))*delX-0.25*(V(m-1,i)+V(m,i))*delX;
                        elseif n1==m1+1
                            A(m1,n1)=0.25*delX*(V(m,i)+V(m+1,i))-delX/Re/delY;
                        elseif n1==m1-1
                            A(m1,n1)=-0.25*delX*(V(m-1,i)+V(m,i))-delX/Re/delY;
                        else 
                            A(m1,n1)=0;
                        end
                    end
            end
             
            b=[];
            for m1=1:nyMax-3
                    m=m1+1;
                    b(m1)=-Vtilda(m,i+1)*(0.25*delY*(U(m,i)+U(m+1,i))-delY/Re/delX)...
                        -Vtilda(m,i-1)*(-0.25*delY*(U(m,i-1)+U(m+1,i-1))-delY/Re/delX)...
                        +delX*delY/deltaT*V(m,i)...
                        +(Ptilda(m,i)-Ptilda(m+1,i))*delX;
                    if m1==1
                        b(m1)=b(m1)-Vtilda(m-1,i)*(-0.25*delX*(V(m-1,i)+V(m,i))-delX/Re/delY);
                    elseif m1==nxMax-2
                        b(m1)=b(m1)-Vtilda(m+1,i)*(0.25*delX*(V(m,i)+V(m+1,i))-delX/Re/delY);
                    end
             end
                
             Vtilda(2:nyMax-2,i)=thomasAlgorithm(A,b);                            
            
        end
        
        % solving for pressure correction (Poisson Solver) %%%%%%%%%%%%%%%
        Pcorr=zeros(nyMax,nxMax);
        
        % first row sweep
        
        for j=2:nyMax-1
            
            A=[];   % coefficient matrix
            for m1=1:nxMax-2
                for n1=1:nxMax-2
                       m=m1+1;     
                       a1=delX*delY/deltaT+0.25*(U(j,m)+U(j,m+1))*delY+...
                                2*(delX/Re/delY+delY/Re/delX)-0.25*(U(j,m-1)+U(j,m))*delY+...
                                0.25*(V(j,m)+V(j,m+1))*delX-0.25*(V(j-1,m)+V(j-1,m+1))*delX;
                       E1=deltaT*a1/(delX*delY);
                        
                        
                        if m==2    % inlet boundary
                            a2=1;
                            E2=0;
                        else
                            a2=delX*delY/deltaT+0.25*(U(j,m-1)+U(j,m))*delY+...
                                2*(delX/Re/delY+delY/Re/delX)-0.25*(U(j,m-2)+U(j,m-1))*delY+...
                                0.25*(V(j,m-1)+V(j,m))*delX-0.25*(V(j-1,m-1)+V(j-1,m))*delX;
                            E2=deltaT*a2/(delX*delY);
                        end
                        
                        if j==nyMax-1 % top wall boundary
                            a3=1;
                            E3=0;
                        else
                            a3=delX*delY/deltaT+0.25*(U(j,m)+U(j+1,m))*delY+...
                                    2*(delX/Re/delY+delY/Re/delX)-0.25*(U(j,m-1)+U(j+1,m-1))*delY...
                                    +0.25*(V(j,m)+V(j+1,m))*delX-0.25*(V(j-1,m)+V(j,m))*delX;

                            E3=deltaT*a3/(delX*delY);
                        end
                            
                        if j==2 % bottom wall boundary
                            a4=1;
                            E4=0;
                        else
                            a4=delX*delY/deltaT+0.25*(U(j-1,m)+U(j,m))*delY+...
                                2*(delX/Re/delY+delY/Re/delX)-0.25*(U(j-1,m-1)+U(j,m-1))*delY...
                                +0.25*(V(j-1,m)+V(j,m))*delX-0.25*(V(j-2,m)+V(j-1,m))*delX;
                            E4=deltaT*a4/(delX*delY);
                        end
                        
                    if m1==n1
                        A(m1,n1)=E1*delY^2/(1+E1)/a1+E2*delY^2/(1+E2)/a2 +E3*delY^2/(1+E3)/a3 +E4*delY^2/(1+E4)/a4;                        
                    elseif n1==m1+1
                        A(m1,n1)=-E1*delY^2/(1+E1)/a1;
                    elseif n1==m1-1
                        A(m1,n1)=-E2*delY^2/(1+E2)/a2;                        
                    end

                end
            end
            
            b=[];   % RHS vector
            
            for m1=1:nxMax-2
                m=m1+1;
                b(m1)=(Utilda(j,m-1)-Utilda(j,m))*delY+(Vtilda(j-1,m)-Vtilda(j,m))*delX...
                    -Pcorr(j+1,m)*(-E3*delY^2/(1+E3)/a3)-Pcorr(j-1,m)*(-E4*delY^2/(1+E4)/a4);
            end
            
            Pcorr(j,2:nxMax-1)=thomasAlgorithm(A,b);
            
        end
        
        % now the column sweep
        
        for i=2:nxMax-1
            
            A=[];
            for m1=1:nyMax-2
                for n1=1:nyMax-2
                       m=m1+1;     
                       a1=delX*delY/deltaT+0.25*(U(m,i)+U(m,i+1))*delY+...
                                    2*(delX/Re/delY+delY/Re/delX)-0.25*(U(m,i-1)+U(m,i))*delY+...
                                    0.25*(V(m,i)+V(m,i+1))*delX-0.25*(V(m-1,i)+V(m-1,i+1))*delX;
                       E1=deltaT*a1/(delX*delY);
                        
                        
                        if i==2    % inlet boundary
                            a2=1;
                            E2=0;
                        else
                            a2=delX*delY/deltaT+0.25*(U(m,i-1)+U(m,i))*delY+...
                                    2*(delX/Re/delY+delY/Re/delX)-0.25*(U(m,i-2)+U(m,i-1))*delY+...
                                    0.25*(V(m,i-1)+V(m,i))*delX-0.25*(V(m-1,i-1)+V(m-1,i))*delX;
                            E2=deltaT*a2/(delX*delY);
                        end
                        
                        if m==nyMax-1 % top wall boundary
                            a3=1;
                            E3=0;
                        else
                            a3=delX*delY/deltaT+0.25*(U(m,i)+U(m+1,i))*delY+...
                                    2*(delX/Re/delY+delY/Re/delX)-0.25*(U(m,i-1)+U(m+1,i-1))*delY...
                                    +0.25*(V(m,i)+V(m+1,i))*delX-0.25*(V(m-1,i)+V(m,i))*delX;
                            E3=deltaT*a3/(delX*delY);
                        end
                            
                        if m==2 % bottom wall boundary
                            a4=1;
                            E4=0;
                        else
                            a4=delX*delY/deltaT+0.25*(U(m-1,i)+U(m,i))*delY+...
                                    2*(delX/Re/delY+delY/Re/delX)-0.25*(U(m-1,i-1)+U(m,i-1))*delY...
                                    +0.25*(V(m-1,i)+V(m,i))*delX-0.25*(V(m-2,i)+V(m-1,i))*delX;
                            E4=deltaT*a4/(delX*delY);
                        end
                        
                    if m1==n1
                        A(m1,n1)=E1*delY^2/(1+E1)/a1+E2*delY^2/(1+E2)/a2 +E3*delY^2/(1+E3)/a3 +E4*delY^2/(1+E4)/a4;                        
                    elseif n1==m1+1
                        A(m1,n1)=-E3*delY^2/(1+E3)/a3;
                    elseif n1==m1-1
                        A(m1,n1)=-E4*delY^2/(1+E4)/a4;                        
                    end

                end
            end
            
            b=[];   % RHS vector
            
            for m1=1:nyMax-2
                m=m1+1;
                b(m1)=(Utilda(m,i-1)-Utilda(m,i))*delY+(Vtilda(m-1,i)-Vtilda(m,i))*delX...
                    -Pcorr(m,i+1)*(-E1*delY^2/(1+E1)/a1)-Pcorr(m,i-1)*(-E2*delY^2/(1+E2)/a2);
            end
            
            Pcorr(2:nyMax-1,i)=thomasAlgorithm(A,b);
                        
            
        end
        
        % CORRECTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % now that Pcorr is calculated,
        % we'll find the velocity corrections using the pressure
        % corrections
        Ucorr=zeros(nyMax,nxMax);
        Vcorr=zeros(nyMax,nxMax);
        
        for i=2:nxMax-1
            for j=2:nyMax-1
                a1=delX*delY/deltaT+0.25*(U(j,i)+U(j,i+1))*delY+...
                                    2*(delX/Re/delY+delY/Re/delX)-0.25*(U(j,i-1)+U(j,i))*delY+...
                                    0.25*(V(j,i)+V(j,i+1))*delX-0.25*(V(j-1,i)+V(j-1,i+1))*delX;
                E1=deltaT*a1/(delX*delY);
                Ucorr(j,i)=E1*delY*(Pcorr(j,i)-Pcorr(j,i+1))/(1+E1)/a1;
                
                a2=delX*delY/deltaT+0.25*(U(j,i)+U(j+1,i))*delY+...
                                    2*(delX/Re/delY+delY/Re/delX)-0.25*(U(j,i-1)+U(j+1,i-1))*delY...
                                    +0.25*(V(j,i)+V(j+1,i))*delX-0.25*(V(j-1,i)+V(j,i))*delX;
                E2=deltaT*a2/(delX*delY);
                Vcorr(j,i)=E2*delY*(Pcorr(j,i)-Pcorr(j+1,i))/(1+E2)/a2;
            end
        end
                               
        
        % finally calculating the error in continuity equation
        bp=0;
        noOfTerms=0;
        for i=2:nxMax-1
            for j=2:nyMax-1
                bp=bp+((Utilda(j,i-1)-Utilda(j,i))*delY+(Vtilda(j-1,i)-Vtilda(j,i))*delX)^2;
                noOfTerms=noOfTerms+1;
            end
        end
        rmsError=sqrt(bp/noOfTerms)
        
        % correcting fields if error is high
        
        if rmsError>tolerance
            Utilda=Utilda+alphaU*Ucorr;
            Vtilda=Vtilda+alphaU*Vcorr;
            Ptilda=Ptilda+alphaP*Pcorr;
        end
            
        iter=iter+1
        
        %%%%%%%%%%%%%%%%%%%%%%%% end of Iterative loop %%%%%%%%%%%%
    end
    
    % saving present time step data
    
    U=Utilda;
    V=Vtilda;
    P=Ptilda;
    
    % Writing Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % writing fields into plt files
     % before writing files, we need to interpolate the surface scalars from th eface centres to
     % the internal cell centres.
     % for interpolation function IDWinterp is used which is defined below
    
     % Ucc means U interpolated back to cell centre points
     Pcc=P(2:nyMax-1,2:nxMax-1);
     Ucc=[];
     Vcc=[];
     
     for j=1:nyMax-2
         for i=1:nxMax-2
            % collect the neighbouring points into 
            m=j+1;
            n=i+1;
            u2=U(m,n-1);
            u4=U(m,n);
            u1=0.25*(U(m,n)+U(m,n-1)+U(m-1,n-1)+U(m-1,n));
            u3=0.25*(U(m,n)+U(m,n-1)+U(m+1,n)+U(m+1,n-1));
            
            v1=V(m-1,n);
            v2=0.25*(V(m-1,n)+V(m,n)+V(m-1,n-1)+V(m,n-1));
            v3=V(m,n);
            v4=0.25*(V(m-1,n)+V(m,n)+V(m-1,n+1)+V(m,n+1));
            
            Ucc(j,i)=((u2+u4)/delX+(u1+u3)/delY)/(2/delX+2/delY);
            Vcc(j,i)=((v2+v4)/delX+(v1+v3)/delY)/(2/delX+2/delY);
         end
     end
     
h=pcolor(x,y,Ucc);
set(h, 'EdgeColor', 'none');
%colormap hot
r=colorbar('fontweight','bold','fontname','times','fontsize',15);
ylabel(r,'U');
xlabel('{\it x}','fontweight','bold','FontName','times','FontSize',25)
ylabel('{\it y}','fontweight','bold','FontName','times','FontSize',25)
pause(0.0005);

end
    
%% for plotting purpose

figure;
plot(y,Ucc(:,89))

figure;
imagesc(U)
%% functions

% Thomas Algorithm 

function result=thomasAlgorithm(A,b)
% here, A= the coefficient matrix (LHS) and b=constant column vector (RHS)
result=[];
m=length(b);

% transforming the vector and matrix

% modification:
for i=1:m-1
    if i==1
        A(i,i+1)=A(i,i+1)/A(i,i);        
    else
        A(i,i+1)=A(i,i+1)/(A(i,i)-A(i,i-1)*A(i-1,i));        
    end    
end

for i=1:m
    if i==1        
        b(i)=b(i)/A(i,i);
    else
        b(i)=(b(i)-A(i,i-1)*b(i-1))/(A(i,i)-A(i,i-1)*A(i-1,i));
    end    
end


% back substitution:
for i=m:(-1):1
    if i==m
        result(i,1)=b(m);
    else
        result(i,1)=b(i)-A(i,i+1)*result(i+1,1);
    end
end

end


% Inverse Distance Weighting (IDW) interpolation (interp)

function result=IDWinterp(A,p)

% A is an input matrix that has different inputs in each rows
% the first three column represents the position vector and the forth
% column represents the variable to be interpolated

% p is the 3d position vector (row vector) where we want the interpolated
% value

num=0;
denom=0;

for i=1:length(A(:,1))
    r=A(i,1:3)-p;
    dist=norm(r);
    
    num=num+A(i,4)/dist;
    denom=denom+1/dist;
end

% result is the value of the variable at point 'p'
result = num/denom;
end