function [S,nb,nl,Busdata,Linedata]=read_data(Com_data)
clear all;
Com_data ='ieee14cdf.txt';
cdf=fopen(Com_data); % open data file

%% Bus data
k_str=fgetl(cdf);  %first line of data
S=str2num(k_str(32:37));
k_str=fgetl(cdf);  %second line of data
Busdata=[];
nb=0;   %no of buses
while ischar(k_str)  %loop till the end 
      k_str=fgetl(cdf);
      if(strcmp(k_str(1:4),'-999') ==1);
          break;
      end
      indx=19; k_str_nm=k_str(indx:end);
      k_nm=str2num(k_str_nm);
      nb=nb+1;
      Busdata=[Busdata; [nb k_nm]];
end         
%% Line Data
k_str=fgetl(cdf);
Linedata=[];
nm=0; % number of lines
while ischar(k_str)
      k_str=fgetl(cdf);
      if(strcmp(k_str(1:4),'-999') ==1);
          break;
      end
      indx=1; k_str_nm=k_str(indx:end);
      k_nm=str2num(k_str_nm);
      nm=nm+1;
      Linedata=[Linedata; [nm k_nm]];
end

      
      
      
    
      
 
      
      
          
      



