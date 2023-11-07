
function [P_cal,Q_cal]=cal_PQ(V_mag,Y_mag,Theta,V_delta,nb) %calculate nodes'power flow

P_cal=zeros(nb,1);  %initialize active power vector
Q_cal=zeros(nb,1);  %initialize reactive power vector
     for i=1:nb
         for j=1:nb
             P_cal(i)=P_cal(i)+V_mag(i)*V_mag(j)*Y_mag(i,j)*cos(Theta(i,j)-V_delta(i)+V_delta(j));
             Q_cal(i)=Q_cal(i)-V_mag(i)*V_mag(j)*Y_mag(i,j)*sin(Theta(i,j)-V_delta(i)+V_delta(j));
         end
     end
end


