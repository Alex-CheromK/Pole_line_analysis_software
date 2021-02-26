function [Function,Jacobian] = CatFunc_Length(HorForce,VerForce_B,Weight,SpWeight,Stiffness,HorDist,VerDist)
VerForce_A = VerForce_B - Weight;
Function = [
 (HorForce/SpWeight)*((VerForce_B - VerForce_A)/Stiffness + asinh(VerForce_B/HorForce) - asinh(VerForce_A/HorForce)) - HorDist;
 (1/SpWeight)*((VerForce_B^2 - VerForce_A^2)/(2*Stiffness) + HorForce*(sqrt(1 + (VerForce_B/HorForce)^2) - sqrt(1 + (VerForce_A/HorForce)^2))) - VerDist];
J11 = (asinh(VerForce_B/HorForce) + asinh((Weight - VerForce_B)/HorForce) + Weight/Stiffness)/SpWeight - (HorForce*((Weight - VerForce_B)/(HorForce^2*((Weight - VerForce_B)^2/HorForce^2 + 1)^(1/2)) + VerForce_B/(HorForce^2*(VerForce_B^2/HorForce^2 + 1)^(1/2))))/SpWeight;
J12 = (HorForce*(1/(HorForce*(VerForce_B^2/HorForce^2 + 1)^(1/2)) - 1/(HorForce*((Weight - VerForce_B)^2/HorForce^2 + 1)^(1/2))))/SpWeight;
J21 = -(((Weight - VerForce_B)^2/HorForce^2 + 1)^(1/2) - (VerForce_B^2/HorForce^2 + 1)^(1/2) + HorForce*(VerForce_B^2/(HorForce^3*(VerForce_B^2/HorForce^2 + 1)^(1/2)) - (Weight - VerForce_B)^2/(HorForce^3*((Weight - VerForce_B)^2/HorForce^2 + 1)^(1/2))))/SpWeight;
J22 = (HorForce*((2*Weight - 2*VerForce_B)/(2*HorForce^2*((Weight - VerForce_B)^2/HorForce^2 + 1)^(1/2)) + VerForce_B/(HorForce^2*(VerForce_B^2/HorForce^2 + 1)^(1/2))) + Weight/Stiffness)/SpWeight;
Jacobian = [
 J11 J12;
 J21 J22];