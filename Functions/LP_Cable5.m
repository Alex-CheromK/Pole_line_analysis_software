function [NetForce,PointPosVec,SpringForceA,SpringForceB] = LP_Cable5(Vars,Span,Global)
%% Extract variables
PointPosVec = reshape(Vars,3,[]);

%% Initialize arrays
SpringForceA = zeros(3,Span.NumPoints);
SpringForceB = zeros(3,Span.NumPoints);
GravForce = zeros(3,Span.NumPoints);
DragForce = zeros(3,Span.NumPoints);
NetForce = zeros(3,Span.NumPoints);

% Build updated main span structure
TmpSpan.NumPoints = Span.NumPoints;
TmpSpan.PointA_PosVec = Span.PointA_PosVec;
TmpSpan.PointB_PosVec = Span.PointB_PosVec;
TmpSpan.PointPosVec = [Span.PointA_PosVec,PointPosVec,Span.PointB_PosVec];

%% Compute total load on each point
for PointNum = 1:Span.NumPoints
 % Segment position vectors
 if PointNum == 1
  SegPosVecA = Span.PointA_PosVec - PointPosVec(:,PointNum);
  SegPosVecB = PointPosVec(:,PointNum + 1) - PointPosVec(:,PointNum);
  ParPosVec = PointPosVec(:,PointNum + 1) - Span.PointA_PosVec;
 elseif PointNum == Span.NumPoints
  SegPosVecA = PointPosVec(:,PointNum - 1) - PointPosVec(:,PointNum);
  SegPosVecB = Span.PointB_PosVec - PointPosVec(:,PointNum);
  ParPosVec = Span.PointB_PosVec - PointPosVec(:,PointNum - 1);
 else
  SegPosVecA = PointPosVec(:,PointNum - 1) - PointPosVec(:,PointNum);
  SegPosVecB = PointPosVec(:,PointNum + 1) - PointPosVec(:,PointNum);
  ParPosVec = PointPosVec(:,PointNum + 1) - PointPosVec(:,PointNum - 1);
 end
 
 % Spring force
 SegLengthA = norm(SegPosVecA);
 SegLengthB = norm(SegPosVecB);
 SegUnitVecA = SegPosVecA/SegLengthA;
 SegUnitVecB = SegPosVecB/SegLengthB;
 SpringForceA(:,PointNum) = Span.ElmStiff*(SegLengthA - Span.ElmLength)*SegUnitVecA;
 SpringForceB(:,PointNum) = Span.ElmStiff*(SegLengthB - Span.ElmLength)*SegUnitVecB;
  
 % Gravitational force
 GravForce(:,PointNum) = [0 0 -Span.PointWeight]';
 
 % Drag force
 ParUnitVec = ParPosVec/norm(ParPosVec);
 ParWindVelVec = dot(Global.WindVelVec,ParUnitVec)*ParUnitVec;
 PerpWindVelVec = Global.WindVelVec - ParWindVelVec;
 % ParDragForce = 0.5*Global.ParDragCoeff*Global.AirDensity*pi*Span.Dia*Span.ElmLength*norm(ParWindVelVec)*ParWindVelVec;
 ParDragForce = 0;
 PerpDragForce = 0.5*Global.PerpDragCoeff*Global.AirDensity*Span.Dia*Span.ElmLength*norm(PerpWindVelVec)*PerpWindVelVec;
 DragForce(:,PointNum) = ParDragForce + PerpDragForce;
 
 % Numerical force
 NumForce = (1e-6)*ones(3,1);
 
 % Total force
 NetForce(:,PointNum) = SpringForceA(:,PointNum) + SpringForceB(:,PointNum) + GravForce(:,PointNum) + DragForce(:,PointNum) + NumForce;
end

% Apply weight from equipment
for SpanEquipNum = 1:Span.NumEquip
 [EquipDistRatio,EquipPointNum_A,EquipPointNum_B,~,~,~] = InterFunc(TmpSpan,Span.Equip(SpanEquipNum).MountDist);
 if EquipPointNum_A == 0
  NetForce(3,EquipPointNum_B) = NetForce(3,EquipPointNum_B) - EquipDistRatio*Span.Equip(SpanEquipNum).Weight;
 elseif EquipPointNum_B == 0
  NetForce(3,EquipPointNum_A) = NetForce(3,EquipPointNum_A) - (1 - EquipDistRatio)*Span.Equip(SpanEquipNum).Weight;
 else
  NetForce(3,EquipPointNum_A) = NetForce(3,EquipPointNum_A) - (1 - EquipDistRatio)*Span.Equip(SpanEquipNum).Weight;
  NetForce(3,EquipPointNum_B) = NetForce(3,EquipPointNum_B) - EquipDistRatio*Span.Equip(SpanEquipNum).Weight;
 end
end