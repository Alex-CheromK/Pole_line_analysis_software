function [DistRatio,PointNum_A,PointNum_B] = PoleInterFunc(Pole,InterDist)
CumDist = 0;
for PointNum = 2:Pole.NumPointsTotal
 CumDist = CumDist + Pole.ElmLength;
 if InterDist <= CumDist
  PointNum_A = PointNum - 1;
  PointNum_B = PointNum;
  break;
 end
end
DistRatio = (InterDist - Pole.PointPosVec(3,PointNum_A))/(Pole.PointPosVec(3,PointNum_B) - Pole.PointPosVec(3,PointNum_A));
