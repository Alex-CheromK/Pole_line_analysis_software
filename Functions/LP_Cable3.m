function [Output,MainPointPosVec,MainSpringForceA,MainSpringForceB,InterPointPosVec,InterSpringForceA,InterSpringForceB] = LP_Cable3(Vars,MainSpan,InterSpan,NumInterSpan_Active,TotalNumInterSpanPoints,Solver)
%% Extract variables
Vars = reshape(Vars,3,[]);
MainPointPosVec = Vars(:,1:MainSpan.NumPoints);
Vars(:,1:MainSpan.NumPoints) = [];
InterPointPosVec = zeros(3,Solver.MaxNumCablePoints,NumInterSpan_Active);
for InterSpanNum = 1:NumInterSpan_Active
 InterPointPosVec(:,1:InterSpan(InterSpanNum).NumPoints,InterSpanNum) = Vars(:,1:InterSpan(InterSpanNum).NumPoints);
 Vars(:,1:InterSpan(InterSpanNum).NumPoints) = [];
end

% Build updated main span structure
TmpMainSpan.NumPoints = MainSpan.NumPoints;
TmpMainSpan.PointA_PosVec = MainSpan.PointA_PosVec;
TmpMainSpan.PointB_PosVec = MainSpan.PointB_PosVec;
TmpMainSpan.PointPosVec = [MainSpan.PointA_PosVec,MainPointPosVec,MainSpan.PointB_PosVec];

%% Initialize arrays
MainSpringForceA = zeros(3,Solver.MaxNumCablePoints);
MainSpringForceB = zeros(3,Solver.MaxNumCablePoints);
MainGravForce = zeros(3,Solver.MaxNumCablePoints);
MainNetForce = zeros(3,Solver.MaxNumCablePoints);
InterSpringForceA = zeros(3,Solver.MaxNumCablePoints,NumInterSpan_Active);
InterSpringForceB = zeros(3,Solver.MaxNumCablePoints,NumInterSpan_Active);
InterGravForce = zeros(3,Solver.MaxNumCablePoints,NumInterSpan_Active);
InterNetForce = zeros(3,Solver.MaxNumCablePoints,NumInterSpan_Active);

%% Compute loads for all LP intersecting cable points
DistRatio = zeros(1,NumInterSpan_Active);
PointNum_A = zeros(1,NumInterSpan_Active);
PointNum_B = zeros(1,NumInterSpan_Active);
Point_PosVec = zeros(3,NumInterSpan_Active);
for InterSpanNum = 1:NumInterSpan_Active
 % Find intersection point along main span
 [DistRatio(InterSpanNum),PointNum_A(InterSpanNum),PointNum_B(InterSpanNum),~,~,Point_PosVec(:,InterSpanNum)] = InterFunc(TmpMainSpan,InterSpan(InterSpanNum).InterDist);
 
 % Compute internal cable loads
 for InterPointNum = 1:InterSpan(InterSpanNum).NumPoints
  % Segment position vectors
  if InterPointNum == 1
   SegPosVecA = InterSpan(InterSpanNum).PointA_PosVec - InterPointPosVec(:,InterPointNum,InterSpanNum);
   SegPosVecB = InterPointPosVec(:,InterPointNum + 1,InterSpanNum) - InterPointPosVec(:,InterPointNum,InterSpanNum);
  elseif InterPointNum == InterSpan(InterSpanNum).NumPoints
   SegPosVecA = InterPointPosVec(:,InterPointNum - 1,InterSpanNum) - InterPointPosVec(:,InterPointNum,InterSpanNum);
   SegPosVecB = Point_PosVec(:,InterSpanNum) - InterPointPosVec(:,InterPointNum,InterSpanNum);
  else
   SegPosVecA = InterPointPosVec(:,InterPointNum - 1,InterSpanNum) - InterPointPosVec(:,InterPointNum,InterSpanNum);
   SegPosVecB = InterPointPosVec(:,InterPointNum + 1,InterSpanNum) - InterPointPosVec(:,InterPointNum,InterSpanNum);
  end
  
  % Spring force
  SegLengthA = norm(SegPosVecA);
  SegLengthB = norm(SegPosVecB);
  SegUnitVecA = SegPosVecA/SegLengthA;
  SegUnitVecB = SegPosVecB/SegLengthB;
  InterSpringForceA(:,InterPointNum,InterSpanNum) = InterSpan(InterSpanNum).ElmStiff*(SegLengthA - InterSpan(InterSpanNum).ElmLength)*SegUnitVecA;
  InterSpringForceB(:,InterPointNum,InterSpanNum) = InterSpan(InterSpanNum).ElmStiff*(SegLengthB - InterSpan(InterSpanNum).ElmLength)*SegUnitVecB;
  
  % Gravitational force
  InterGravForce(:,InterPointNum,InterSpanNum) = [0 0 -InterSpan(InterSpanNum).PointWeight]';
  
  % Numerical force
  NumForce = (1e-6)*ones(3,1);
  
  % Total force
  InterNetForce(:,InterPointNum,InterSpanNum) = InterSpringForceA(:,InterPointNum,InterSpanNum) + InterSpringForceB(:,InterPointNum,InterSpanNum) + InterGravForce(:,InterPointNum,InterSpanNum) + NumForce;
 end
end

%% Compute loads for LP main cable points
for MainPointNum = 1:MainSpan.NumPoints
 % Segment position vectors
 if MainPointNum == 1
  SegPosVecA = MainSpan.PointA_PosVec - MainPointPosVec(:,MainPointNum);
  SegPosVecB = MainPointPosVec(:,MainPointNum + 1) - MainPointPosVec(:,MainPointNum);
 elseif MainPointNum == MainSpan.NumPoints
  SegPosVecA = MainPointPosVec(:,MainPointNum - 1) - MainPointPosVec(:,MainPointNum);
  SegPosVecB = MainSpan.PointB_PosVec - MainPointPosVec(:,MainPointNum);
 else
  SegPosVecA = MainPointPosVec(:,MainPointNum - 1) - MainPointPosVec(:,MainPointNum);
  SegPosVecB = MainPointPosVec(:,MainPointNum + 1) - MainPointPosVec(:,MainPointNum);
 end
 
 % Spring force
 SegLengthA = norm(SegPosVecA);
 SegLengthB = norm(SegPosVecB);
 SegUnitVecA = SegPosVecA/SegLengthA;
 SegUnitVecB = SegPosVecB/SegLengthB;
 MainSpringForceA(:,MainPointNum) = MainSpan.ElmStiff*(SegLengthA - MainSpan.ElmLength)*SegUnitVecA;
 MainSpringForceB(:,MainPointNum) = MainSpan.ElmStiff*(SegLengthB - MainSpan.ElmLength)*SegUnitVecB;
 
 % Gravitational force
 MainGravForce(:,MainPointNum) = [0 0 -MainSpan.PointWeight]';
  
 % Numerical force
 NumForce = (1e-6)*ones(3,1);
 
 % Total force
 MainNetForce(:,MainPointNum) = MainSpringForceA(:,MainPointNum) + MainSpringForceB(:,MainPointNum) + MainGravForce(:,MainPointNum) + NumForce;
end

% Add loads from LP intersecting cables
for InterSpanNum = 1:NumInterSpan_Active
 if PointNum_A(InterSpanNum) == 0
  MainNetForce(:,PointNum_B(InterSpanNum)) = MainNetForce(:,PointNum_B(InterSpanNum)) - DistRatio(InterSpanNum)*InterSpringForceB(:,InterSpan(InterSpanNum).NumPoints,InterSpanNum);
 elseif PointNum_B(InterSpanNum) == 0
  MainNetForce(:,PointNum_A(InterSpanNum)) = MainNetForce(:,PointNum_A(InterSpanNum)) - (1 - DistRatio(InterSpanNum))*InterSpringForceB(:,InterSpan(InterSpanNum).NumPoints,InterSpanNum);
 else
  MainNetForce(:,PointNum_A(InterSpanNum)) = MainNetForce(:,PointNum_A(InterSpanNum)) - (1 - DistRatio(InterSpanNum))*InterSpringForceB(:,InterSpan(InterSpanNum).NumPoints,InterSpanNum);
  MainNetForce(:,PointNum_B(InterSpanNum)) = MainNetForce(:,PointNum_B(InterSpanNum)) - DistRatio(InterSpanNum)*InterSpringForceB(:,InterSpan(InterSpanNum).NumPoints,InterSpanNum);
 end
end

%% Collect final output
Tmp = zeros(3,TotalNumInterSpanPoints);
RefNumPoints = 0;
for InterSpanNum = 1:NumInterSpan_Active
 Tmp(:,RefNumPoints + 1:RefNumPoints + InterSpan(InterSpanNum).NumPoints) = InterNetForce(:,1:InterSpan(InterSpanNum).NumPoints,InterSpanNum);
 RefNumPoints = RefNumPoints + InterSpan(InterSpanNum).NumPoints;
end
Output = [
 reshape(MainNetForce(:,1:MainSpan.NumPoints),[],1);
 reshape(Tmp,[],1)];