function [Output,Pole,Guy,SideBar] = GuyAssemblyFunc(Vars,Pole,Guy,Anch,SideBar,Solver)
OriginalVars = Vars;
for PoleSideBarNum = 1:Pole.NumSideBars
 SideBarNum = Pole.SideBar(PoleSideBarNum).SideBarNum;
 SideBar(SideBarNum).RelTipDef = Vars(PoleSideBarNum);
end
Vars(1:Pole.NumSideBars) = [];
Vars = reshape(Vars,4,[]);
Pole.PointDispVec(:,2:Pole.NumPointsTotal) = Vars(1:2,:);
Pole.PointRotVec(:,2:Pole.NumPointsTotal) = Vars(3:4,:);
Pole.PointPosVec(1:2,2:Pole.NumPointsTotal) = [
 Pole.PosX + Pole.PointDispVec(1,2:Pole.NumPointsTotal);
 Pole.PosY + Pole.PointDispVec(2,2:Pole.NumPointsTotal)];
Pole.PointForceVec = zeros(size(Pole.PointForceVec));
Pole.PointMomVec = zeros(size(Pole.PointMomVec));

%% Compute loads from guy wires
for PoleGuyNum = 1:Pole.NumGuys
 GuyNum = Pole.Guy(PoleGuyNum).GuyNum;
 AnchNum = Guy(GuyNum).AnchNum;
 SideBarNum = Guy(GuyNum).SideBarNum;
 for GuySpanNum = 1:Guy(GuyNum).NumSpans
  if (Guy(GuyNum).SideBarNum > 0) && (GuySpanNum == 1)
   % Update sidebar end point position vectors
   SideBar(SideBarNum).BasePosVec = [
    interp1(Pole.PointPosVec(3,1:Pole.NumPointsTotal),Pole.PointPosVec(1,1:Pole.NumPointsTotal),SideBar(SideBarNum).MountHeight);
    interp1(Pole.PointPosVec(3,1:Pole.NumPointsTotal),Pole.PointPosVec(2,1:Pole.NumPointsTotal),SideBar(SideBarNum).MountHeight);
    Pole.PosZ + SideBar(SideBarNum).MountHeight];
   SideBar(SideBarNum).TipDisp =...
    -SideBar(SideBarNum).Length*cosd(SideBar(SideBarNum).Angle)*interp1(Pole.PointPosVec(3,1:Pole.NumPointsTotal),Pole.PointRotVec(1,1:Pole.NumPointsTotal),SideBar(SideBarNum).MountHeight)...
    - SideBar(SideBarNum).Length*sind(SideBar(SideBarNum).Angle)*interp1(Pole.PointPosVec(3,1:Pole.NumPointsTotal),Pole.PointRotVec(2,1:Pole.NumPointsTotal),SideBar(SideBarNum).MountHeight)...
    + SideBar(SideBarNum).RelTipDef;
   SideBar(SideBarNum).TipPosVec = [
    SideBar(SideBarNum).BasePosVec(1) + SideBar(SideBarNum).Length*cosd(SideBar(SideBarNum).Angle);
    SideBar(SideBarNum).BasePosVec(2) + SideBar(SideBarNum).Length*sind(SideBar(SideBarNum).Angle);
    SideBar(SideBarNum).BasePosVec(3) + SideBar(SideBarNum).TipDisp];
   
   % Update span end point position vectors
   Guy(GuyNum).Span(GuySpanNum).PointA_PosVec = [
    Anch(AnchNum).PosX;
    Anch(AnchNum).PosY;
    Anch(AnchNum).PosZ];
   Guy(GuyNum).Span(GuySpanNum).PointB_PosVec = SideBar(SideBarNum).TipPosVec;
   
   % Solve linear cable problem
   Guy(GuyNum).Span(GuySpanNum).BaseLength = norm(Guy(GuyNum).Span(GuySpanNum).PointB_PosVec - Guy(GuyNum).Span(GuySpanNum).PointA_PosVec)/(1 + Guy(GuyNum).Tension/Guy(GuyNum).Stiffness);
   Guy(GuyNum).Span(GuySpanNum).Weight = Guy(GuyNum).Span(GuySpanNum).BaseLength*Guy(GuyNum).BaseSpWeight;
   
   % Discretize span
   Guy(GuyNum).Span(GuySpanNum).NumElm = ceil(Guy(GuyNum).Span(GuySpanNum).BaseLength/Solver.CableElmLength);
   Guy(GuyNum).Span(GuySpanNum).NumPointsTotal = Guy(GuyNum).Span(GuySpanNum).NumElm + 1;
   Guy(GuyNum).Span(GuySpanNum).NumPoints = Guy(GuyNum).Span(GuySpanNum).NumPointsTotal - 2;
   Guy(GuyNum).Span(GuySpanNum).PointWeight = Guy(GuyNum).Span(GuySpanNum).Weight/Guy(GuyNum).Span(GuySpanNum).NumPointsTotal;
   Guy(GuyNum).Span(GuySpanNum).ElmLength = Guy(GuyNum).Span(GuySpanNum).BaseLength/Guy(GuyNum).Span(GuySpanNum).NumElm;
   Guy(GuyNum).Span(GuySpanNum).ElmStiff = Guy(GuyNum).Stiffness/Guy(GuyNum).Span(GuySpanNum).ElmLength;
   
   % Build initial span point position vector
   GuySpanPointPosVecInit = [
    linspace(Guy(GuyNum).Span(GuySpanNum).PointA_PosVec(1),Guy(GuyNum).Span(GuySpanNum).PointB_PosVec(1),Guy(GuyNum).Span(GuySpanNum).NumPointsTotal);
    linspace(Guy(GuyNum).Span(GuySpanNum).PointA_PosVec(2),Guy(GuyNum).Span(GuySpanNum).PointB_PosVec(2),Guy(GuyNum).Span(GuySpanNum).NumPointsTotal);
    linspace(Guy(GuyNum).Span(GuySpanNum).PointA_PosVec(3),Guy(GuyNum).Span(GuySpanNum).PointB_PosVec(3),Guy(GuyNum).Span(GuySpanNum).NumPointsTotal)];
   GuySpanPointPosVecInit(:,1) = []; GuySpanPointPosVecInit(:,end) = [];
   InitVars = reshape(GuySpanPointPosVecInit,[],1);
   
   % Solve span LP problem
   Func = @(Vars) LP_Cable(Vars,Guy(GuyNum).Span(GuySpanNum));
   Vars = fsolve(Func,InitVars,Solver.FsolveOptions);
   [~,PointPosVec,SpringForceA,SpringForceB] = LP_Cable(Vars,Guy(GuyNum).Span(GuySpanNum));
   Guy(GuyNum).Span(GuySpanNum).PointPosVec(:,1:Guy(GuyNum).Span(GuySpanNum).NumPointsTotal) = [Guy(GuyNum).Span(GuySpanNum).PointA_PosVec,PointPosVec,Guy(GuyNum).Span(GuySpanNum).PointB_PosVec];
   Guy(GuyNum).Span(GuySpanNum).SpringForceA(:,1:Guy(GuyNum).Span(GuySpanNum).NumPoints) = SpringForceA;
   Guy(GuyNum).Span(GuySpanNum).SpringForceB(:,1:Guy(GuyNum).Span(GuySpanNum).NumPoints) = SpringForceB;
  
   % Apply loads to pole
   [DistRatio,PointNum_A,PointNum_B] = PoleInterFunc(Pole,SideBar(SideBarNum).MountHeight);
   Pole.PointForceVec(:,PointNum_A) = Pole.PointForceVec(:,PointNum_A) - (1 - DistRatio)*Guy(GuyNum).Span(GuySpanNum).SpringForceB(1:2,Guy(GuyNum).Span(GuySpanNum).NumPoints);
   Pole.PointForceVec(:,PointNum_B) = Pole.PointForceVec(:,PointNum_B) - DistRatio*Guy(GuyNum).Span(GuySpanNum).SpringForceB(1:2,Guy(GuyNum).Span(GuySpanNum).NumPoints);
   Pole.PointMomVec(:,PointNum_A) = Pole.PointMomVec(:,PointNum_A) + [
    (1 - DistRatio)*Guy(GuyNum).Span(GuySpanNum).SpringForceB(3,Guy(GuyNum).Span(GuySpanNum).NumPoints)*SideBar(SideBarNum).Length*cosd(SideBar(SideBarNum).Angle);
    (1 - DistRatio)*Guy(GuyNum).Span(GuySpanNum).SpringForceB(3,Guy(GuyNum).Span(GuySpanNum).NumPoints)*SideBar(SideBarNum).Length*sind(SideBar(SideBarNum).Angle)];
   Pole.PointMomVec(:,PointNum_B) = Pole.PointMomVec(:,PointNum_B) + [
    DistRatio*Guy(GuyNum).Span(GuySpanNum).SpringForceB(3,Guy(GuyNum).Span(GuySpanNum).NumPoints)*SideBar(SideBarNum).Length*cosd(SideBar(SideBarNum).Angle);
    DistRatio*Guy(GuyNum).Span(GuySpanNum).SpringForceB(3,Guy(GuyNum).Span(GuySpanNum).NumPoints)*SideBar(SideBarNum).Length*sind(SideBar(SideBarNum).Angle)];
  
  elseif (Guy(GuyNum).SideBarNum > 0) && (GuySpanNum == 2)
   % Update sidebar end point position vectors
   SideBar(SideBarNum).BasePosVec = [
    interp1(Pole.PointPosVec(3,1:Pole.NumPointsTotal),Pole.PointPosVec(1,1:Pole.NumPointsTotal),SideBar(SideBarNum).MountHeight);
    interp1(Pole.PointPosVec(3,1:Pole.NumPointsTotal),Pole.PointPosVec(2,1:Pole.NumPointsTotal),SideBar(SideBarNum).MountHeight);
    Pole.PosZ + SideBar(SideBarNum).MountHeight];
   SideBar(SideBarNum).TipDisp =...
    -SideBar(SideBarNum).Length*cosd(SideBar(SideBarNum).Angle)*interp1(Pole.PointPosVec(3,1:Pole.NumPointsTotal),Pole.PointRotVec(1,1:Pole.NumPointsTotal),SideBar(SideBarNum).MountHeight)...
    - SideBar(SideBarNum).Length*sind(SideBar(SideBarNum).Angle)*interp1(Pole.PointPosVec(3,1:Pole.NumPointsTotal),Pole.PointRotVec(2,1:Pole.NumPointsTotal),SideBar(SideBarNum).MountHeight)...
    + SideBar(SideBarNum).RelTipDef;
   SideBar(SideBarNum).TipPosVec = [
    SideBar(SideBarNum).BasePosVec(1) + SideBar(SideBarNum).Length*cosd(SideBar(SideBarNum).Angle);
    SideBar(SideBarNum).BasePosVec(2) + SideBar(SideBarNum).Length*sind(SideBar(SideBarNum).Angle);
    SideBar(SideBarNum).BasePosVec(3) + SideBar(SideBarNum).TipDisp];
   
   % Update span end point position vectors
   Guy(GuyNum).Span(GuySpanNum).PointA_PosVec = SideBar(SideBarNum).TipPosVec;
   Guy(GuyNum).Span(GuySpanNum).PointB_PosVec = [
    interp1(Pole.PointPosVec(3,1:Pole.NumPointsTotal),Pole.PointPosVec(1,1:Pole.NumPointsTotal),Guy(GuyNum).AttachHeight);
    interp1(Pole.PointPosVec(3,1:Pole.NumPointsTotal),Pole.PointPosVec(2,1:Pole.NumPointsTotal),Guy(GuyNum).AttachHeight);
    Pole.PosZ + Guy(GuyNum).AttachHeight];
   
   % Solve catenary problem
   SpanHorPosVec = Guy(GuyNum).Span(GuySpanNum).PointB_PosVec(1:2) - Guy(GuyNum).Span(GuySpanNum).PointA_PosVec(1:2);
   SpanHorDist = norm(SpanHorPosVec);
   SpanVerDist = Guy(GuyNum).Span(GuySpanNum).PointB_PosVec(3) - Guy(GuyNum).Span(GuySpanNum).PointA_PosVec(3);
   SpanTension_A = Guy(GuyNum).Tension;
   SpanAngle_A_Init = atan(SpanVerDist/SpanHorDist);
   SpanVerForce_B_Init = SpanTension_A*sin(SpanAngle_A_Init);
   Func = @(Vars) CatFunc_Tension(Vars(1),Vars(2),SpanTension_A,Guy(GuyNum).BaseSpWeight,Guy(GuyNum).Stiffness,SpanHorDist,SpanVerDist);
   Vars = fsolve(Func,[SpanAngle_A_Init SpanVerForce_B_Init]',Solver.FsolveOptions);
   Guy(GuyNum).Span(GuySpanNum).Angle_A = Vars(1);
   Guy(GuyNum).Span(GuySpanNum).HorForce = SpanTension_A*cos(Guy(GuyNum).Span(GuySpanNum).Angle_A);
   Guy(GuyNum).Span(GuySpanNum).VerForce_A = SpanTension_A*sin(Guy(GuyNum).Span(GuySpanNum).Angle_A);
   Guy(GuyNum).Span(GuySpanNum).VerForce_B = Vars(2);
   
   % Find length-constrained span properties
   Guy(GuyNum).Span(GuySpanNum).Weight = Guy(GuyNum).Span(GuySpanNum).VerForce_B - Guy(GuyNum).Span(GuySpanNum).VerForce_A;
   Guy(GuyNum).Span(GuySpanNum).BaseLength = Guy(GuyNum).Span(GuySpanNum).Weight/Guy(GuyNum).BaseSpWeight;
   
   % Discretize span
   Guy(GuyNum).Span(GuySpanNum).NumElm = ceil(Guy(GuyNum).Span(GuySpanNum).BaseLength/Solver.CableElmLength);
   Guy(GuyNum).Span(GuySpanNum).NumPointsTotal = Guy(GuyNum).Span(GuySpanNum).NumElm + 1;
   Guy(GuyNum).Span(GuySpanNum).NumPoints = Guy(GuyNum).Span(GuySpanNum).NumPointsTotal - 2;
   Guy(GuyNum).Span(GuySpanNum).PointWeight = Guy(GuyNum).Span(GuySpanNum).Weight/Guy(GuyNum).Span(GuySpanNum).NumPointsTotal;
   Guy(GuyNum).Span(GuySpanNum).ElmLength = Guy(GuyNum).Span(GuySpanNum).BaseLength/Guy(GuyNum).Span(GuySpanNum).NumElm;
   Guy(GuyNum).Span(GuySpanNum).ElmStiff = Guy(GuyNum).Stiffness/Guy(GuyNum).Span(GuySpanNum).ElmLength;
  
   % Build initial span point position vector
   GuySpanHorVec = Guy(GuyNum).Span(GuySpanNum).PointB_PosVec(1:2) - Guy(GuyNum).Span(GuySpanNum).PointA_PosVec(1:2);
   GuySpanHorDist = norm(GuySpanHorVec);
   HeightA = Guy(GuyNum).Span(GuySpanNum).PointA_PosVec(3);
   HeightB = Guy(GuyNum).Span(GuySpanNum).PointB_PosVec(3);
   a = Guy(GuyNum).Span(GuySpanNum).HorForce/Guy(GuyNum).SpWeight;
   A = exp(-GuySpanHorDist/a) - 1;
   B = -2*(HeightB - HeightA)/a;
   C = exp(GuySpanHorDist/a) - 1;
   Lambda = (-B - sqrt(B^2 - 4*A*C))/(2*A);
   if Lambda <= 0
    error('Solution to catenary profile is a complex number!');
   end
   k1 = a*log(Lambda);
   k2 = HeightA - a*cosh(k1/a);
   GuySpanHorAxisVec = linspace(0,GuySpanHorDist,Guy(GuyNum).Span(GuySpanNum).NumPointsTotal)';
   GuySpanHorUnitVec = GuySpanHorVec/GuySpanHorDist;
   GuySpanHorPointA_PosVec = Guy(GuyNum).Span(GuySpanNum).PointA_PosVec(1:2);
   GuySpanPointPosVecInit = [
    GuySpanHorPointA_PosVec + GuySpanHorAxisVec'.*GuySpanHorUnitVec;
    a*cosh((GuySpanHorAxisVec' - k1)/a) + k2];
   GuySpanPointPosVecInit(:,1) = []; GuySpanPointPosVecInit(:,end) = [];
   InitVars = reshape(GuySpanPointPosVecInit,[],1);
   
   % Solve span LP problem
   Func = @(Vars) LP_Cable(Vars,Guy(GuyNum).Span(GuySpanNum));
   Vars = fsolve(Func,InitVars,Solver.FsolveOptions);
   [~,PointPosVec,SpringForceA,SpringForceB] = LP_Cable(Vars,Guy(GuyNum).Span(GuySpanNum));
   Guy(GuyNum).Span(GuySpanNum).PointPosVec(:,1:Guy(GuyNum).Span(GuySpanNum).NumPointsTotal) = [Guy(GuyNum).Span(GuySpanNum).PointA_PosVec,PointPosVec,Guy(GuyNum).Span(GuySpanNum).PointB_PosVec];
   Guy(GuyNum).Span(GuySpanNum).SpringForceA(:,1:Guy(GuyNum).Span(GuySpanNum).NumPoints) = SpringForceA;
   Guy(GuyNum).Span(GuySpanNum).SpringForceB(:,1:Guy(GuyNum).Span(GuySpanNum).NumPoints) = SpringForceB;
   
   % Apply loads to pole
   [DistRatio,PointNum_A,PointNum_B] = PoleInterFunc(Pole,SideBar(SideBarNum).MountHeight);
   Pole.PointForceVec(:,PointNum_A) = Pole.PointForceVec(:,PointNum_A) - (1 - DistRatio)*Guy(GuyNum).Span(GuySpanNum).SpringForceA(1:2,1);
   Pole.PointForceVec(:,PointNum_B) = Pole.PointForceVec(:,PointNum_B) - DistRatio*Guy(GuyNum).Span(GuySpanNum).SpringForceA(1:2,1);
   Pole.PointMomVec(:,PointNum_A) = Pole.PointMomVec(:,PointNum_A) + [
    (1 - DistRatio)*Guy(GuyNum).Span(GuySpanNum).SpringForceA(3,1)*SideBar(SideBarNum).Length*cosd(SideBar(SideBarNum).Angle);
    (1 - DistRatio)*Guy(GuyNum).Span(GuySpanNum).SpringForceA(3,1)*SideBar(SideBarNum).Length*sind(SideBar(SideBarNum).Angle)];
   Pole.PointMomVec(:,PointNum_B) = Pole.PointMomVec(:,PointNum_B) + [
    DistRatio*Guy(GuyNum).Span(GuySpanNum).SpringForceA(3,1)*SideBar(SideBarNum).Length*cosd(SideBar(SideBarNum).Angle);
    DistRatio*Guy(GuyNum).Span(GuySpanNum).SpringForceA(3,1)*SideBar(SideBarNum).Length*sind(SideBar(SideBarNum).Angle)];
   [DistRatio,PointNum_A,PointNum_B] = PoleInterFunc(Pole,Guy(GuyNum).AttachHeight);
   Pole.PointForceVec(:,PointNum_A) = Pole.PointForceVec(:,PointNum_A) - (1 - DistRatio)*Guy(GuyNum).Span(GuySpanNum).SpringForceB(1:2,Guy(GuyNum).Span(GuySpanNum).NumPoints);
   Pole.PointForceVec(:,PointNum_B) = Pole.PointForceVec(:,PointNum_B) - DistRatio*Guy(GuyNum).Span(GuySpanNum).SpringForceB(1:2,Guy(GuyNum).Span(GuySpanNum).NumPoints);
   
  else
   % Update span end point position vectors
   Guy(GuyNum).Span(GuySpanNum).PointA_PosVec = [
    Anch(AnchNum).PosX;
    Anch(AnchNum).PosY;
    Anch(AnchNum).PosZ];
   Guy(GuyNum).Span(GuySpanNum).PointB_PosVec = [
    interp1(Pole.PointPosVec(3,1:Pole.NumPointsTotal),Pole.PointPosVec(1,1:Pole.NumPointsTotal),Guy(GuyNum).AttachHeight);
    interp1(Pole.PointPosVec(3,1:Pole.NumPointsTotal),Pole.PointPosVec(2,1:Pole.NumPointsTotal),Guy(GuyNum).AttachHeight);
    Pole.PosZ + Guy(GuyNum).AttachHeight];
   
   % Solve catenary problem
   SpanHorPosVec = Guy(GuyNum).Span(GuySpanNum).PointB_PosVec(1:2) - Guy(GuyNum).Span(GuySpanNum).PointA_PosVec(1:2);
   SpanHorDist = norm(SpanHorPosVec);
   SpanVerDist = Guy(GuyNum).Span(GuySpanNum).PointB_PosVec(3) - Guy(GuyNum).Span(GuySpanNum).PointA_PosVec(3);
   SpanTension_A = Guy(GuyNum).Tension;
   SpanAngle_A_Init = atan(SpanVerDist/SpanHorDist);
   SpanVerForce_B_Init = SpanTension_A*sin(SpanAngle_A_Init);
   Func = @(Vars) CatFunc_Tension(Vars(1),Vars(2),SpanTension_A,Guy(GuyNum).BaseSpWeight,Guy(GuyNum).Stiffness,SpanHorDist,SpanVerDist);
   Vars = fsolve(Func,[SpanAngle_A_Init SpanVerForce_B_Init]',Solver.FsolveOptions);
   Guy(GuyNum).Span(GuySpanNum).Angle_A = Vars(1);
   Guy(GuyNum).Span(GuySpanNum).HorForce = SpanTension_A*cos(Guy(GuyNum).Span(GuySpanNum).Angle_A);
   Guy(GuyNum).Span(GuySpanNum).VerForce_A = SpanTension_A*sin(Guy(GuyNum).Span(GuySpanNum).Angle_A);
   Guy(GuyNum).Span(GuySpanNum).VerForce_B = Vars(2);
   
   % Find length-constrained span properties
   Guy(GuyNum).Span(GuySpanNum).Weight = Guy(GuyNum).Span(GuySpanNum).VerForce_B - Guy(GuyNum).Span(GuySpanNum).VerForce_A;
   Guy(GuyNum).Span(GuySpanNum).BaseLength = Guy(GuyNum).Span(GuySpanNum).Weight/Guy(GuyNum).BaseSpWeight;
   
   % Discretize span
   Guy(GuyNum).Span(GuySpanNum).NumElm = ceil(Guy(GuyNum).Span(GuySpanNum).BaseLength/Solver.CableElmLength);
   Guy(GuyNum).Span(GuySpanNum).NumPointsTotal = Guy(GuyNum).Span(GuySpanNum).NumElm + 1;
   Guy(GuyNum).Span(GuySpanNum).NumPoints = Guy(GuyNum).Span(GuySpanNum).NumPointsTotal - 2;
   Guy(GuyNum).Span(GuySpanNum).PointWeight = Guy(GuyNum).Span(GuySpanNum).Weight/Guy(GuyNum).Span(GuySpanNum).NumPointsTotal;
   Guy(GuyNum).Span(GuySpanNum).ElmLength = Guy(GuyNum).Span(GuySpanNum).BaseLength/Guy(GuyNum).Span(GuySpanNum).NumElm;
   Guy(GuyNum).Span(GuySpanNum).ElmStiff = Guy(GuyNum).Stiffness/Guy(GuyNum).Span(GuySpanNum).ElmLength;
  
   % Build initial span point position vector
   GuySpanHorVec = Guy(GuyNum).Span(GuySpanNum).PointB_PosVec(1:2) - Guy(GuyNum).Span(GuySpanNum).PointA_PosVec(1:2);
   GuySpanHorDist = norm(GuySpanHorVec);
   HeightA = Guy(GuyNum).Span(GuySpanNum).PointA_PosVec(3);
   HeightB = Guy(GuyNum).Span(GuySpanNum).PointB_PosVec(3);
   a = Guy(GuyNum).Span(GuySpanNum).HorForce/Guy(GuyNum).SpWeight;
   A = exp(-GuySpanHorDist/a) - 1;
   B = -2*(HeightB - HeightA)/a;
   C = exp(GuySpanHorDist/a) - 1;
   Lambda = (-B - sqrt(B^2 - 4*A*C))/(2*A);
   if Lambda <= 0
    error('Solution to catenary profile is a complex number!');
   end
   k1 = a*log(Lambda);
   k2 = HeightA - a*cosh(k1/a);
   GuySpanHorAxisVec = linspace(0,GuySpanHorDist,Guy(GuyNum).Span(GuySpanNum).NumPointsTotal)';
   GuySpanHorUnitVec = GuySpanHorVec/GuySpanHorDist;
   GuySpanHorPointA_PosVec = Guy(GuyNum).Span(GuySpanNum).PointA_PosVec(1:2);
   GuySpanPointPosVecInit = [
    GuySpanHorPointA_PosVec + GuySpanHorAxisVec'.*GuySpanHorUnitVec;
    a*cosh((GuySpanHorAxisVec' - k1)/a) + k2];
   GuySpanPointPosVecInit(:,1) = []; GuySpanPointPosVecInit(:,end) = [];
   InitVars = reshape(GuySpanPointPosVecInit,[],1);
   
   % Solve span LP problem
   Func = @(Vars) LP_Cable(Vars,Guy(GuyNum).Span(GuySpanNum));
   Vars = fsolve(Func,InitVars,Solver.FsolveOptions);
   [~,PointPosVec,SpringForceA,SpringForceB] = LP_Cable(Vars,Guy(GuyNum).Span(GuySpanNum));
   Guy(GuyNum).Span(GuySpanNum).PointPosVec(:,1:Guy(GuyNum).Span(GuySpanNum).NumPointsTotal) = [Guy(GuyNum).Span(GuySpanNum).PointA_PosVec,PointPosVec,Guy(GuyNum).Span(GuySpanNum).PointB_PosVec];
   Guy(GuyNum).Span(GuySpanNum).SpringForceA(:,1:Guy(GuyNum).Span(GuySpanNum).NumPoints) = SpringForceA;
   Guy(GuyNum).Span(GuySpanNum).SpringForceB(:,1:Guy(GuyNum).Span(GuySpanNum).NumPoints) = SpringForceB;
   
   % Apply loads to pole
   [DistRatio,PointNum_A,PointNum_B] = PoleInterFunc(Pole,Guy(GuyNum).AttachHeight);
   Pole.PointForceVec(:,PointNum_A) = Pole.PointForceVec(:,PointNum_A) - (1 - DistRatio)*Guy(GuyNum).Span(GuySpanNum).SpringForceB(1:2,Guy(GuyNum).Span(GuySpanNum).NumPoints);
   Pole.PointForceVec(:,PointNum_B) = Pole.PointForceVec(:,PointNum_B) - DistRatio*Guy(GuyNum).Span(GuySpanNum).SpringForceB(1:2,Guy(GuyNum).Span(GuySpanNum).NumPoints);
  end
 end
end

%% Compute pole displacements
% Stiffness matrix
StiffMat = Pole.StiffMat(1:2*Pole.NumPointsTotal,1:2*Pole.NumPointsTotal);
StiffMat(1:2,:) = [];StiffMat(:,1:2) = [];

% x direction
LoadVec = [
 Pole.PointForceVec(1,1:Pole.NumPointsTotal);
 Pole.PointMomVec(1,1:Pole.NumPointsTotal)];
LoadVec(:,1) = [];
LoadVec = reshape(LoadVec,[],1);
Tmp = reshape(StiffMat\LoadVec,2,[]);
Pole.PointDispVec(1,2:Pole.NumPointsTotal) = Tmp(1,:);
Pole.PointRotVec(1,2:Pole.NumPointsTotal) = Tmp(2,:);

% y direction
LoadVec = [
 Pole.PointForceVec(2,1:Pole.NumPointsTotal);
 Pole.PointMomVec(2,1:Pole.NumPointsTotal)];
LoadVec(:,1) = [];
LoadVec = reshape(LoadVec,[],1);
Tmp = reshape(StiffMat\LoadVec,2,[]);
Pole.PointDispVec(2,2:Pole.NumPointsTotal) = Tmp(1,:);
Pole.PointRotVec(2,2:Pole.NumPointsTotal) = Tmp(2,:);

%% Compute sidebar displacements
for PoleSideBarNum = 1:Pole.NumSideBars
 SideBarNum = Pole.SideBar(PoleSideBarNum).SideBarNum;
 SideBarTipForce = 0;
 for SideBarGuyNum = 1:SideBar(SideBarNum).NumGuys
  GuyNum = SideBar(SideBarNum).Guy(SideBarGuyNum).GuyNum;
  SideBarTipForce = SideBarTipForce - Guy(GuyNum).Span(1).SpringForceB(3,Guy(GuyNum).Span(1).NumPoints);
  SideBarTipForce = SideBarTipForce - Guy(GuyNum).Span(2).SpringForceA(3,1);
 end
 SideBar(SideBarNum).RelTipDef = SideBarTipForce*(SideBar(SideBarNum).Length^3)/(3*SideBar(SideBarNum).ElasticMod*SideBar(SideBarNum).AMOI);
end

%% Build output vector
SideBarOutput = zeros(Pole.NumSideBars,1);
for PoleSideBarNum = 1:Pole.NumSideBars
 SideBarNum = Pole.SideBar(PoleSideBarNum).SideBarNum;
 SideBarOutput(PoleSideBarNum) = SideBar(SideBarNum).RelTipDef;
end
PoleOutput = [
 Pole.PointDispVec(:,2:Pole.NumPointsTotal);
 Pole.PointRotVec(:,2:Pole.NumPointsTotal)];
Output = [
 SideBarOutput;
 reshape(PoleOutput,[],1)];
Output = Output - OriginalVars;