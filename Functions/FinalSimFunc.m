function [Output,Pole,Line,Guy,SideBar] = FinalSimFunc(Global,Vars,Pole,Line,Guy,Anch,SideBar,Solver,PoleNum,System)
%% Update pole structures
OriginalVars = Vars;
for PoleSideBarNum = 1:Pole(PoleNum).NumSideBars
 SideBarNum = Pole(PoleNum).SideBar(PoleSideBarNum).SideBarNum;
 SideBar(SideBarNum).RelTipDef = Vars(PoleSideBarNum);
end
Vars(1:Pole(PoleNum).NumSideBars) = [];
Vars = reshape(Vars,4,[]);
Pole(PoleNum).PointDispVec(:,2:Pole(PoleNum).NumPointsTotal) = Vars(1:2,:);
Pole(PoleNum).PointRotVec(:,2:Pole(PoleNum).NumPointsTotal) = Vars(3:4,:);
Pole(PoleNum).PointPosVec(1:2,2:Pole(PoleNum).NumPointsTotal) = [
 Pole(PoleNum).PosX + Pole(PoleNum).PointDispVec(1,2:Pole(PoleNum).NumPointsTotal);
 Pole(PoleNum).PosY + Pole(PoleNum).PointDispVec(2,2:Pole(PoleNum).NumPointsTotal)];
Pole(PoleNum).PointForceVec = zeros(size(Pole(PoleNum).PointForceVec));
Pole(PoleNum).PointMomVec = zeros(size(Pole(PoleNum).PointMomVec));

%% Update line tensions
Line = LineSimFunc(Solver,System,Pole,Line,Global);

% Apply loads to (selected) pole
for PoleLineAttachNum = 1:Pole(PoleNum).NumLineAttach
 LineNum = Pole(PoleNum).Attach(PoleLineAttachNum).LineNum;
 LineAttachNum = Pole(PoleNum).Attach(PoleLineAttachNum).LineAttachNum;
 LineSpanNum = LineAttachNum;
 if LineAttachNum == 1
  % Span B
  [DistRatio,PointNum_A,PointNum_B] = PoleInterFunc(Pole(PoleNum),Line(LineNum).Attach(LineAttachNum).Height);
   Pole(PoleNum).PointForceVec(:,PointNum_A) = Pole(PoleNum).PointForceVec(:,PointNum_A) - (1 - DistRatio)*Line(LineNum).Span(LineSpanNum).SpringForceA(1:2,1);
   Pole(PoleNum).PointForceVec(:,PointNum_B) = Pole(PoleNum).PointForceVec(:,PointNum_B) - DistRatio*Line(LineNum).Span(LineSpanNum).SpringForceA(1:2,1);
 elseif LineAttachNum == Line(LineNum).NumAttach
  % Span A
  [DistRatio,PointNum_A,PointNum_B] = PoleInterFunc(Pole(PoleNum),Line(LineNum).Attach(LineAttachNum).Height);
   Pole(PoleNum).PointForceVec(:,PointNum_A) = Pole(PoleNum).PointForceVec(:,PointNum_A) - (1 - DistRatio)*Line(LineNum).Span(LineSpanNum - 1).SpringForceB(1:2,Line(LineNum).NumSpans);
   Pole(PoleNum).PointForceVec(:,PointNum_B) = Pole(PoleNum).PointForceVec(:,PointNum_B) - DistRatio*Line(LineNum).Span(LineSpanNum - 1).SpringForceB(1:2,Line(LineNum).NumSpans);
 else
  % Span A
  [DistRatio,PointNum_A,PointNum_B] = PoleInterFunc(Pole(PoleNum),Line(LineNum).Attach(LineAttachNum).Height);
   Pole(PoleNum).PointForceVec(:,PointNum_A) = Pole(PoleNum).PointForceVec(:,PointNum_A) - (1 - DistRatio)*Line(LineNum).Span(LineSpanNum - 1).SpringForceB(1:2,Line(LineNum).NumSpans);
   Pole(PoleNum).PointForceVec(:,PointNum_B) = Pole(PoleNum).PointForceVec(:,PointNum_B) - DistRatio*Line(LineNum).Span(LineSpanNum - 1).SpringForceB(1:2,Line(LineNum).NumSpans);
   % Span B
  [DistRatio,PointNum_A,PointNum_B] = PoleInterFunc(Pole(PoleNum),Line(LineNum).Attach(LineAttachNum).Height);
   Pole(PoleNum).PointForceVec(:,PointNum_A) = Pole(PoleNum).PointForceVec(:,PointNum_A) - (1 - DistRatio)*Line(LineNum).Span(LineSpanNum).SpringForceA(1:2,1);
   Pole(PoleNum).PointForceVec(:,PointNum_B) = Pole(PoleNum).PointForceVec(:,PointNum_B) - DistRatio*Line(LineNum).Span(LineSpanNum).SpringForceA(1:2,1);
 end
end

%% Compute loads from guy wires
for PoleGuyNum = 1:Pole(PoleNum).NumGuys
 GuyNum = Pole(PoleNum).Guy(PoleGuyNum).GuyNum;
 AnchNum = Guy(GuyNum).AnchNum;
 SideBarNum = Guy(GuyNum).SideBarNum;
 for GuySpanNum = 1:Guy(GuyNum).NumSpans
  if (Guy(GuyNum).SideBarNum > 0) && (GuySpanNum == 1)
   % Update sidebar end point position vectors
   SideBar(SideBarNum).BasePosVec = [
    interp1(Pole(PoleNum).PointPosVec(3,1:Pole(PoleNum).NumPointsTotal),Pole(PoleNum).PointPosVec(1,1:Pole(PoleNum).NumPointsTotal),SideBar(SideBarNum).MountHeight);
    interp1(Pole(PoleNum).PointPosVec(3,1:Pole(PoleNum).NumPointsTotal),Pole(PoleNum).PointPosVec(2,1:Pole(PoleNum).NumPointsTotal),SideBar(SideBarNum).MountHeight);
    Pole(PoleNum).PosZ + SideBar(SideBarNum).MountHeight];
   SideBar(SideBarNum).TipDisp =...
    -SideBar(SideBarNum).Length*cosd(SideBar(SideBarNum).Angle)*interp1(Pole(PoleNum).PointPosVec(3,1:Pole(PoleNum).NumPointsTotal),Pole(PoleNum).PointRotVec(1,1:Pole(PoleNum).NumPointsTotal),SideBar(SideBarNum).MountHeight)...
    - SideBar(SideBarNum).Length*sind(SideBar(SideBarNum).Angle)*interp1(Pole(PoleNum).PointPosVec(3,1:Pole(PoleNum).NumPointsTotal),Pole(PoleNum).PointRotVec(2,1:Pole(PoleNum).NumPointsTotal),SideBar(SideBarNum).MountHeight)...
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
   [DistRatio,PointNum_A,PointNum_B] = PoleInterFunc(Pole(PoleNum),SideBar(SideBarNum).MountHeight);
   Pole(PoleNum).PointForceVec(:,PointNum_A) = Pole(PoleNum).PointForceVec(:,PointNum_A) - (1 - DistRatio)*Guy(GuyNum).Span(GuySpanNum).SpringForceB(1:2,Guy(GuyNum).Span(GuySpanNum).NumPoints);
   Pole(PoleNum).PointForceVec(:,PointNum_B) = Pole(PoleNum).PointForceVec(:,PointNum_B) - DistRatio*Guy(GuyNum).Span(GuySpanNum).SpringForceB(1:2,Guy(GuyNum).Span(GuySpanNum).NumPoints);
   Pole(PoleNum).PointMomVec(:,PointNum_A) = Pole(PoleNum).PointMomVec(:,PointNum_A) + [
    (1 - DistRatio)*Guy(GuyNum).Span(GuySpanNum).SpringForceB(3,Guy(GuyNum).Span(GuySpanNum).NumPoints)*SideBar(SideBarNum).Length*cosd(SideBar(SideBarNum).Angle);
    (1 - DistRatio)*Guy(GuyNum).Span(GuySpanNum).SpringForceB(3,Guy(GuyNum).Span(GuySpanNum).NumPoints)*SideBar(SideBarNum).Length*sind(SideBar(SideBarNum).Angle)];
   Pole(PoleNum).PointMomVec(:,PointNum_B) = Pole(PoleNum).PointMomVec(:,PointNum_B) + [
    DistRatio*Guy(GuyNum).Span(GuySpanNum).SpringForceB(3,Guy(GuyNum).Span(GuySpanNum).NumPoints)*SideBar(SideBarNum).Length*cosd(SideBar(SideBarNum).Angle);
    DistRatio*Guy(GuyNum).Span(GuySpanNum).SpringForceB(3,Guy(GuyNum).Span(GuySpanNum).NumPoints)*SideBar(SideBarNum).Length*sind(SideBar(SideBarNum).Angle)];
   
  elseif (Guy(GuyNum).SideBarNum > 0) && (GuySpanNum == 2)
   % Update sidebar end point position vectors
   SideBar(SideBarNum).BasePosVec = [
    interp1(Pole(PoleNum).PointPosVec(3,1:Pole(PoleNum).NumPointsTotal),Pole(PoleNum).PointPosVec(1,1:Pole(PoleNum).NumPointsTotal),SideBar(SideBarNum).MountHeight);
    interp1(Pole(PoleNum).PointPosVec(3,1:Pole(PoleNum).NumPointsTotal),Pole(PoleNum).PointPosVec(2,1:Pole(PoleNum).NumPointsTotal),SideBar(SideBarNum).MountHeight);
    Pole(PoleNum).PosZ + SideBar(SideBarNum).MountHeight];
   SideBar(SideBarNum).TipDisp =...
    -SideBar(SideBarNum).Length*cosd(SideBar(SideBarNum).Angle)*interp1(Pole(PoleNum).PointPosVec(3,1:Pole(PoleNum).NumPointsTotal),Pole(PoleNum).PointRotVec(1,1:Pole(PoleNum).NumPointsTotal),SideBar(SideBarNum).MountHeight)...
    - SideBar(SideBarNum).Length*sind(SideBar(SideBarNum).Angle)*interp1(Pole(PoleNum).PointPosVec(3,1:Pole(PoleNum).NumPointsTotal),Pole(PoleNum).PointRotVec(2,1:Pole(PoleNum).NumPointsTotal),SideBar(SideBarNum).MountHeight)...
    + SideBar(SideBarNum).RelTipDef;
   SideBar(SideBarNum).TipPosVec = [
    SideBar(SideBarNum).BasePosVec(1) + SideBar(SideBarNum).Length*cosd(SideBar(SideBarNum).Angle);
    SideBar(SideBarNum).BasePosVec(2) + SideBar(SideBarNum).Length*sind(SideBar(SideBarNum).Angle);
    SideBar(SideBarNum).BasePosVec(3) + SideBar(SideBarNum).TipDisp];
   
   % Update span end point position vectors
   Guy(GuyNum).Span(GuySpanNum).PointA_PosVec = SideBar(SideBarNum).TipPosVec;
   Guy(GuyNum).Span(GuySpanNum).PointB_PosVec = [
    interp1(Pole(PoleNum).PointPosVec(3,1:Pole(PoleNum).NumPointsTotal),Pole(PoleNum).PointPosVec(1,1:Pole(PoleNum).NumPointsTotal),Guy(GuyNum).AttachHeight);
    interp1(Pole(PoleNum).PointPosVec(3,1:Pole(PoleNum).NumPointsTotal),Pole(PoleNum).PointPosVec(2,1:Pole(PoleNum).NumPointsTotal),Guy(GuyNum).AttachHeight);
    Pole(PoleNum).PosZ + Guy(GuyNum).AttachHeight];
   
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
   [DistRatio,PointNum_A,PointNum_B] = PoleInterFunc(Pole(PoleNum),SideBar(SideBarNum).MountHeight);
   Pole(PoleNum).PointForceVec(:,PointNum_A) = Pole(PoleNum).PointForceVec(:,PointNum_A) - (1 - DistRatio)*Guy(GuyNum).Span(GuySpanNum).SpringForceA(1:2,1);
   Pole(PoleNum).PointForceVec(:,PointNum_B) = Pole(PoleNum).PointForceVec(:,PointNum_B) - DistRatio*Guy(GuyNum).Span(GuySpanNum).SpringForceA(1:2,1);
   Pole(PoleNum).PointMomVec(:,PointNum_A) = Pole(PoleNum).PointMomVec(:,PointNum_A) + [
    (1 - DistRatio)*Guy(GuyNum).Span(GuySpanNum).SpringForceA(3,1)*SideBar(SideBarNum).Length*cosd(SideBar(SideBarNum).Angle);
    (1 - DistRatio)*Guy(GuyNum).Span(GuySpanNum).SpringForceA(3,1)*SideBar(SideBarNum).Length*sind(SideBar(SideBarNum).Angle)];
   Pole(PoleNum).PointMomVec(:,PointNum_B) = Pole(PoleNum).PointMomVec(:,PointNum_B) + [
    DistRatio*Guy(GuyNum).Span(GuySpanNum).SpringForceA(3,1)*SideBar(SideBarNum).Length*cosd(SideBar(SideBarNum).Angle);
    DistRatio*Guy(GuyNum).Span(GuySpanNum).SpringForceA(3,1)*SideBar(SideBarNum).Length*sind(SideBar(SideBarNum).Angle)];
   [DistRatio,PointNum_A,PointNum_B] = PoleInterFunc(Pole(PoleNum),Guy(GuyNum).AttachHeight);
   Pole(PoleNum).PointForceVec(:,PointNum_A) = Pole(PoleNum).PointForceVec(:,PointNum_A) - (1 - DistRatio)*Guy(GuyNum).Span(GuySpanNum).SpringForceB(1:2,Guy(GuyNum).Span(GuySpanNum).NumPoints);
   Pole(PoleNum).PointForceVec(:,PointNum_B) = Pole(PoleNum).PointForceVec(:,PointNum_B) - DistRatio*Guy(GuyNum).Span(GuySpanNum).SpringForceB(1:2,Guy(GuyNum).Span(GuySpanNum).NumPoints);
   
  else
   % Update span end point position vectors
   Guy(GuyNum).Span(GuySpanNum).PointA_PosVec = [
    Anch(AnchNum).PosX;
    Anch(AnchNum).PosY;
    Anch(AnchNum).PosZ];
   Guy(GuyNum).Span(GuySpanNum).PointB_PosVec = [
    interp1(Pole(PoleNum).PointPosVec(3,1:Pole(PoleNum).NumPointsTotal),Pole(PoleNum).PointPosVec(1,1:Pole(PoleNum).NumPointsTotal),Guy(GuyNum).AttachHeight);
    interp1(Pole(PoleNum).PointPosVec(3,1:Pole(PoleNum).NumPointsTotal),Pole(PoleNum).PointPosVec(2,1:Pole(PoleNum).NumPointsTotal),Guy(GuyNum).AttachHeight);
    Pole(PoleNum).PosZ + Guy(GuyNum).AttachHeight];
   
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
   [DistRatio,PointNum_A,PointNum_B] = PoleInterFunc(Pole(PoleNum),Guy(GuyNum).AttachHeight);
   Pole(PoleNum).PointForceVec(:,PointNum_A) = Pole(PoleNum).PointForceVec(:,PointNum_A) - (1 - DistRatio)*Guy(GuyNum).Span(GuySpanNum).SpringForceB(1:2,Guy(GuyNum).Span(GuySpanNum).NumPoints);
   Pole(PoleNum).PointForceVec(:,PointNum_B) = Pole(PoleNum).PointForceVec(:,PointNum_B) - DistRatio*Guy(GuyNum).Span(GuySpanNum).SpringForceB(1:2,Guy(GuyNum).Span(GuySpanNum).NumPoints);
  end
 end
end

%% Compute pole displacements
% Stiffness matrix
StiffMat = Pole(PoleNum).StiffMat(1:2*Pole(PoleNum).NumPointsTotal,1:2*Pole(PoleNum).NumPointsTotal);
StiffMat(1:2,:) = [];StiffMat(:,1:2) = [];

% x direction
LoadVec = [
 Pole(PoleNum).PointForceVec(1,1:Pole(PoleNum).NumPointsTotal);
 Pole(PoleNum).PointMomVec(1,1:Pole(PoleNum).NumPointsTotal)];
LoadVec = LoadVec + reshape(Pole(PoleNum).WindLoadVecX(1:2*Pole(PoleNum).NumPointsTotal),2,[]);
LoadVec(:,1) = [];
LoadVec = reshape(LoadVec,[],1);
Tmp = reshape(StiffMat\LoadVec,2,[]);
Pole(PoleNum).PointDispVec(1,2:Pole(PoleNum).NumPointsTotal) = Tmp(1,:);
Pole(PoleNum).PointRotVec(1,2:Pole(PoleNum).NumPointsTotal) = Tmp(2,:);

% y direction
LoadVec = [
 Pole(PoleNum).PointForceVec(2,1:Pole(PoleNum).NumPointsTotal);
 Pole(PoleNum).PointMomVec(2,1:Pole(PoleNum).NumPointsTotal)];
LoadVec = LoadVec + reshape(Pole(PoleNum).WindLoadVecY(1:2*Pole(PoleNum).NumPointsTotal),2,[]);
LoadVec(:,1) = [];
LoadVec = reshape(LoadVec,[],1);
Tmp = reshape(StiffMat\LoadVec,2,[]);
Pole(PoleNum).PointDispVec(2,2:Pole(PoleNum).NumPointsTotal) = Tmp(1,:);
Pole(PoleNum).PointRotVec(2,2:Pole(PoleNum).NumPointsTotal) = Tmp(2,:);

%% Compute sidebar displacements
for PoleSideBarNum = 1:Pole(PoleNum).NumSideBars
 SideBarNum = Pole(PoleNum).SideBar(PoleSideBarNum).SideBarNum;
 SideBarTipForce = 0;
 for SideBarGuyNum = 1:SideBar(SideBarNum).NumGuys
  GuyNum = SideBar(SideBarNum).Guy(SideBarGuyNum).GuyNum;
  SideBarTipForce = SideBarTipForce - Guy(GuyNum).Span(1).SpringForceB(3,Guy(GuyNum).Span(1).NumPoints);
  SideBarTipForce = SideBarTipForce - Guy(GuyNum).Span(2).SpringForceA(3,1);
 end
 SideBar(SideBarNum).RelTipDef = SideBarTipForce*(SideBar(SideBarNum).Length^3)/(3*SideBar(SideBarNum).ElasticMod*SideBar(SideBarNum).AMOI);
end

%% Build output vector
SideBarOutput = zeros(Pole(PoleNum).NumSideBars,1);
for PoleSideBarNum = 1:Pole(PoleNum).NumSideBars
 SideBarNum = Pole(PoleNum).SideBar(PoleSideBarNum).SideBarNum;
 SideBarOutput(PoleSideBarNum) = SideBar(SideBarNum).RelTipDef;
end
PoleOutput = [
 Pole(PoleNum).PointDispVec(:,2:Pole(PoleNum).NumPointsTotal);
 Pole(PoleNum).PointRotVec(:,2:Pole(PoleNum).NumPointsTotal)];
Output = [
 SideBarOutput;
 reshape(PoleOutput,[],1)];
Output = Output - OriginalVars;