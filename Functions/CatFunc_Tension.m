function [Function] = CatFunc_Tension(Angle_A,VerForce_B,Tension_A,SpWeight,Stiffness,HorDist,VerDist)
HorForce = Tension_A*cos(Angle_A);
VerForce_A = Tension_A*sin(Angle_A);
Function = [
 (HorForce/SpWeight)*((VerForce_B - VerForce_A)/Stiffness + asinh(VerForce_B/HorForce) - asinh(VerForce_A/HorForce)) - HorDist;
 (1/SpWeight)*((VerForce_B^2 - VerForce_A^2)/(2*Stiffness) + HorForce*(sqrt(1 + (VerForce_B/HorForce)^2) - sqrt(1 + (VerForce_A/HorForce)^2))) - VerDist];
 
 
 