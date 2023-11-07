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


        



      
      
      
    
      
 
      
      
          
      



