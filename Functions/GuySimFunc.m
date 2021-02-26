function [Output,Pole,Guy,SideBar] = GuySimFunc(Vars,Pole,Line,Guy,Anch,SideBar,Solver)
OriginalVars = Vars;
for PoleSideBarNum = 1:Pole.NumSideBars
 SideBarNum = Pole.SideBar(PoleSideBarNum).SideBarNum;
 SideBar(SideBarNum).RelTipDef = Vars(PoleSideBarNum);
end
Vars(1:Pole.NumSideBars) = [];
TmpVars = reshape(Vars,4,[]);
Pole.PointDispVec(:,2:Pole.NumPointsTotal) = TmpVars(1:2,:);
Pole.PointRotVec(:,2:Pole.NumPointsTotal) = TmpVars(3:4,:);
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
   
   % Define initial conditions
   InitVars = Guy(GuyNum).Span(GuySpanNum).PointPosVec(:,2:Guy(GuyNum).Span(GuySpanNum).NumPointsTotal - 1);
   InitVars = reshape(InitVars,[],1);
   
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
   
   % Define initial conditions
   InitVars = Guy(GuyNum).Span(GuySpanNum).PointPosVec(:,2:Guy(GuyNum).Span(GuySpanNum).NumPointsTotal - 1);
   InitVars = reshape(InitVars,[],1);
   
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
   
   % Define initial conditions
   InitVars = Guy(GuyNum).Span(GuySpanNum).PointPosVec(:,2:Guy(GuyNum).Span(GuySpanNum).NumPointsTotal - 1);
   InitVars = reshape(InitVars,[],1);
   
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

%% Compute loads from support wires
for PoleLineAttachNum = 1:Pole.NumLineAttach
 LineNum = Pole.Attach(PoleLineAttachNum).LineNum;
 LineAttachNum = Pole.Attach(PoleLineAttachNum).LineAttachNum;
 LineSpanNum = LineAttachNum;
 if LineAttachNum == 1
  % Span B
  [DistRatio,PointNum_A,PointNum_B] = PoleInterFunc(Pole,Line(LineNum).Attach(LineAttachNum).Height);
   Pole.PointForceVec(:,PointNum_A) = Pole.PointForceVec(:,PointNum_A) - (1 - DistRatio)*Line(LineNum).Span(LineSpanNum).SpringForceA(1:2,1);
   Pole.PointForceVec(:,PointNum_B) = Pole.PointForceVec(:,PointNum_B) - DistRatio*Line(LineNum).Span(LineSpanNum).SpringForceA(1:2,1);
 elseif LineAttachNum == Line(LineNum).NumAttach
  % Span A
  [DistRatio,PointNum_A,PointNum_B] = PoleInterFunc(Pole,Line(LineNum).Attach(LineAttachNum).Height);
   Pole.PointForceVec(:,PointNum_A) = Pole.PointForceVec(:,PointNum_A) - (1 - DistRatio)*Line(LineNum).Span(LineSpanNum - 1).SpringForceB(1:2,Line(LineNum).NumSpans);
   Pole.PointForceVec(:,PointNum_B) = Pole.PointForceVec(:,PointNum_B) - DistRatio*Line(LineNum).Span(LineSpanNum - 1).SpringForceB(1:2,Line(LineNum).NumSpans);
 else
  % Span A
  [DistRatio,PointNum_A,PointNum_B] = PoleInterFunc(Pole,Line(LineNum).Attach(LineAttachNum).Height);
   Pole.PointForceVec(:,PointNum_A) = Pole.PointForceVec(:,PointNum_A) - (1 - DistRatio)*Line(LineNum).Span(LineSpanNum - 1).SpringForceB(1:2,Line(LineNum).NumSpans);
   Pole.PointForceVec(:,PointNum_B) = Pole.PointForceVec(:,PointNum_B) - DistRatio*Line(LineNum).Span(LineSpanNum - 1).SpringForceB(1:2,Line(LineNum).NumSpans);
   % Span B
  [DistRatio,PointNum_A,PointNum_B] = PoleInterFunc(Pole,Line(LineNum).Attach(LineAttachNum).Height);
   Pole.PointForceVec(:,PointNum_A) = Pole.PointForceVec(:,PointNum_A) - (1 - DistRatio)*Line(LineNum).Span(LineSpanNum).SpringForceA(1:2,1);
   Pole.PointForceVec(:,PointNum_B) = Pole.PointForceVec(:,PointNum_B) - DistRatio*Line(LineNum).Span(LineSpanNum).SpringForceA(1:2,1);
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
LoadVec = LoadVec + reshape(Pole.WindLoadVecX(1:2*Pole.NumPointsTotal),2,[]);
LoadVec(:,1) = [];
LoadVec = reshape(LoadVec,[],1);
Tmp = reshape(StiffMat\LoadVec,2,[]);
Pole.PointDispVec(1,2:Pole.NumPointsTotal) = Tmp(1,:);
Pole.PointRotVec(1,2:Pole.NumPointsTotal) = Tmp(2,:);

% y direction
LoadVec = [
 Pole.PointForceVec(2,1:Pole.NumPointsTotal);
 Pole.PointMomVec(2,1:Pole.NumPointsTotal)];
LoadVec = LoadVec + reshape(Pole.WindLoadVecY(1:2*Pole.NumPointsTotal),2,[]);
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