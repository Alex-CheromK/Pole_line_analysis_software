function [NetForce,PointPosVec,SpringForceA,SpringForceB] = LP_Cable(Vars,Span)
%% Extract variables
PointPosVec = reshape(Vars,3,[]);

%% Initialize arrays
SpringForceA = zeros(3,Span.NumPoints);
SpringForceB = zeros(3,Span.NumPoints);
GravForce = zeros(3,Span.NumPoints);
NetForce = zeros(3,Span.NumPoints);

%% Compute total load on each point
for PointNum = 1:Span.NumPoints
 % Segment position vectors
 if PointNum == 1
  SegPosVecA = Span.PointA_PosVec - PointPosVec(:,PointNum);
  SegPosVecB = PointPosVec(:,PointNum + 1) - PointPosVec(:,PointNum);
 elseif PointNum == Span.NumPoints
  SegPosVecA = PointPosVec(:,PointNum - 1) - PointPosVec(:,PointNum);
  SegPosVecB = Span.PointB_PosVec - PointPosVec(:,PointNum);
 else
  SegPosVecA = PointPosVec(:,PointNum - 1) - PointPosVec(:,PointNum);
  SegPosVecB = PointPosVec(:,PointNum + 1) - PointPosVec(:,PointNum);
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
 
 % Numerical force
 NumForce = (1e-6)*ones(3,1);
 
 % Total force
 NetForce(:,PointNum) = SpringForceA(:,PointNum) + SpringForceB(:,PointNum) + GravForce(:,PointNum) + NumForce;
end