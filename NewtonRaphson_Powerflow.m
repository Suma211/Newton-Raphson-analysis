 clear all;
Com_data ='ieee14cdf.txt';
cdf=fopen(Com_data); % open data file
%% Bus data
k_str=fgetl(cdf);  %first line of data
S=str2num(k_str(32:37));
k_str=fgetl(cdf);  %second line of data
Busdata=[]; % empty bus data matrrix
nb=0;   %no of buses
while ischar(k_str)  %loop till the end 
      k_str=fgetl(cdf);
      if(strcmp(k_str(1:4),'-999') ==1);
          break;
      end
      indx=19; k_str_nm=k_str(indx:end);
      k_nm=str2num(k_str_nm);  % K_nm is value of bus data
      nb=nb+1;    %incrementing bus number
      Busdata=[Busdata; [nb k_nm]];
end         
%% Line Data
k_str=fgetl(cdf);
Linedata=[];  % empty linedata matrix
nm=0; % number of lines
while ischar(k_str)
      k_str=fgetl(cdf); % first line of linedata
      if(strcmp(k_str(1:4),'-999') ==1);  %k_str(1:4) reads from 1 to 4
          break;
      end
      indx=1; k_str_nm=k_str(indx:end);
      k_nm=str2num(k_str_nm);   % K_nm is variable to store values
      nm=nm+1;   % no. of lines
      Linedata=[Linedata; [nm k_nm]];
end
% Y bus matrix computation
Y_m=zeros(nb);   %Y matrix
B_m=zeros(nb);   %B matrix
Z_fb=Linedata(:,2);  %  from bus
Z_tb=Linedata(:,3);  %  to bus
R=Linedata(:,8);  % Resistance 
X=Linedata(:,9);   % Reactance
B=Linedata(:,10).*1j;  %susceptance
tran_rt=Linedata(:,16);   %tansformer turns ratio
shun_cap=Busdata(:,15)+(Busdata(:,16)*1i);  %Shunt capacitor at bus 9
Z=R+1j*X;  % line Impedance
Y=1./Z;  % Line Admittance 
for n=1:nm  %nm is no of lines 
    if tran_rt~=0  % finding off diagonal elements of trans adm matrix
        Y_m(Z_fb(n),Z_tb(n))=-1/tran_rt(n)*Y(n);
        Y_m(Z_tb(n),Z_fb(n))=Y_m(Z_fb(n),Z_tb(n));
        Y_m(Z_fb(n),Z_fb(n))=Y_m(Z_fb(n),Z_fb(n))+Y(n)*(1/tran_rt(n))^2;
        Y_m(Z_tb(n),Z_tb(n))=Y_m(Z_tb(n),Z_tb(n))+Y(n);
    else
        Y_m(Z_fb(n),Z_tb(n))=-Y(n);  
        Y_m(Z_tb(n),Z_fb(n))=Y_m(Z_fb(n),Z_tb(n));
        Y_m(Z_fb(n),Z_fb(n))=Y_m(Z_fb(n),Z_fb(n))+Y(n)+B(n)/2;
        Y_m(Z_tb(n),Z_tb(n))=Y_m(Z_tb(n),Z_tb(n))+Y(n)+B(n)/2;
    end 
end 
for r=1:nb  % YBus plus diagonal valid only for 9 th bus
    Y_m(r,r)=Y_m(r,r)+shun_cap(r);
end

%%Reading the busdata from file to calculate Jacobian Matrix

V_mag = Busdata(:,5); % magnitude of voltage magnitude in P.U.
V_ang = Busdata(:,6);% voltage angle
P_load= Busdata(:,7)/S;  % Active power of load in P.U.
Q_load= Busdata(:,8)/S;  % Reactive power of load in p.u.
P_gen = Busdata(:,9)/S;   % Active power of the generation in p.u.
Q_gen= Busdata (:,10)/S;  % Reactive power of the generation in P.U.
Q_max= Busdata(:,13)/S;   % Maximum Reactive power limit in P.U
Q_min = Busdata(:,14)/S;   %MInimum Reactive power limit in P.U

% Types of buses
Slack = find(Busdata(:,4)==3);  % Slack Bus
PQ= find(Busdata(:,4)==0|Busdata(:,4)==1);  % PQ bus
PV= find (Busdata(:,4)==2);  %PV Bus
nSlack=length(Slack);
nPQ=length(PQ);
nPV=length(PV);
% Changing the Ybus matrix from rectangular to polar form
[Theta, Y_mag] = cart2pol(real(Y_m), imag(Y_m)); % get angle and magnitude of admittance matrix
 B = imag(Y_m); % susceptance matrix
 G = real(Y_m); % conductance matrix
% Initializing the parameters
V_mag= Busdata(:,12); %Initial Voltage
V_mag(~V_mag)=1;   %replacing V_mag with value 1
V_delta=zeros(nb,1); %Initial angle
P_sch=P_gen-P_load;  % net power scheduled at bus
Q_sch=Q_gen-Q_load;
dif_Voltage= zeros(nb-1+nPQ,1);

% Power calculation formulas
P_cal=zeros(nb,1);  %initialize active power vector
Q_cal=zeros(nb,1);  %initialize reactive power vector
     for i=1:nb
         for j=1:nb
             P_cal(i)=P_cal(i)+V_mag(i)*V_mag(j)*Y_mag(i,j)*cos(Theta(i,j)-V_delta(i)+V_delta(j));
             Q_cal(i)=Q_cal(i)-V_mag(i)*V_mag(j)*Y_mag(i,j)*sin(Theta(i,j)-V_delta(i)+V_delta(j));
         end 
     end


 %% Newton raphson power flow calculation
 tol_max=0.01; % iteration tolerance of the solution
 iter=0; %times of iteration
 tol=1; %initial value of tol
 
 while (tol>tol_max)
    [P_cal,Q_cal]=cal_PQ(V_mag,Y_mag,Theta,V_delta,nb); %calculate the power at buses 
    [dif_PQ]=difference_PQ(P_sch,Q_sch,P_cal,Q_cal,PQ,nPQ); % mismatches vector
    [J]=Jacobian_Matrix(V_mag,P_cal,Q_cal,Y_mag,Theta, V_delta,nb,PQ,nPQ,B,G); %call the Jacobian matrix
    [dif_Voltage]=LU_factor_PQ(J,dif_PQ); %dif_Voltage=inv(J)*dif_PQ; %get correction vector
    dif_D=dif_Voltage(1:nb-1); % angle correction vector
    dif_V=dif_Voltage(nb:end); % magnitude correction vector
    V_delta(2:end)=V_delta(2:end)+dif_D; %correct the results, angle and voltage
    V_mag(PQ)=V_mag(PQ)+dif_V;
    tol=max(abs(dif_PQ));
    iter=iter+1;
 end
    if iter>5
        disp('bad, not converge');
    else
        V_result=[V_mag,rad2deg(V_delta)];        
    end  
    %% verify the calculation results
if max(V_result-[V_mag,V_delta])<=0.1
    fprintf('Congratulation,converge!, times of iteration=%d.\n',iter);
else
    disp('converge, but results are not correct, please go back to check!');
end
    
 
 %% Jacobian matrix calculation  
function [J]=Jacobian_Matrix(V_mag,P_cal,Q_cal,Y_mag,Theta, V_delta,nb,PQ,nPQ,B,G)
J1=zeros(nb-1);
J2=zeros(nb-1,nPQ);
J3=zeros(nPQ,nb-1);
J4=zeros(nPQ,nPQ);

%J1=dP/dDelta
for i=2:nb %row position of a bus
    for j=2:nb %column position of bus
        if(j==i) 
                       
            J1(i-1,j-1)=-Q_cal(i)-V_mag(i)^2*B(i,j); %diagonal elements  
        else
            J1(i-1,j-1)= -V_mag(i)*V_mag(j)*Y_mag(i,j)*sin(Theta(i,j)-V_delta(i)+V_delta(j));
        end
    end
end

 %J2=dP/dV
 for i=2:nb %position of row
     for j=1:nPQ  %position of column
         jj=PQ(j);         
         if(jj==i) %diagonal elements     
             J2(i-1,j)=P_cal(i)/V_mag(i)+V_mag(i)*G(i,i);
            
         else %off-diagonal elements
            J2(i-1,j)= V_mag(i)*Y_mag(i,jj)*cos(Theta(i,jj)-V_delta(i)+V_delta(jj));
                 
         end
     end
 end
             
 %J3=dQ/dDelta
 for i=1:nPQ %position of row
     ii=PQ(i);
     for j=2:nb %position of column
         if(j==ii) 
             J3(i,j-1)=P_cal(ii)-G(ii,ii)*(V_mag(ii)^2); %diagonal elements
            
         else  %off-diagonal elements
             J3(i,j-1)=-V_mag(ii)*V_mag(j)*Y_mag(ii,j)*cos(Theta(ii,j)-V_delta(ii)+V_delta(j));
         end
     end
 end
 
 %J4=dQ/dV
 for i=1:nPQ %position of row
     ii=PQ(i);
     for j=1:nPQ %position of colume
         jj=PQ(j);
         if(jj==ii) %diagonal elements 
             J4(i,j)=Q_cal(ii)/V_mag(ii)-V_mag(ii)*B(ii,ii);
             
         else %off-diagonal elements
             J4(i,j)=-V_mag(ii)*Y_mag(ii,jj)*sin(Theta(ii,jj)-V_delta(ii)+V_delta(jj));
                
         end
     
 end
   J=[J1 J2;J3 J4];  
 end
 end

        



      
      
      
    
      
 
      
      
          
      



