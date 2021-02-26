function [DistRatio,PointNum_A,PointNum_B,PointA_PosVec,PointB_PosVec,Point_PosVec] = InterFunc(Span,InterDist)
MainSpanVec = Span.PointB_PosVec - Span.PointA_PosVec;
MainSpanHorVec = [
 MainSpanVec(1:2);
 0];
MainSpanHorUnitVec = MainSpanHorVec/norm(MainSpanHorVec);
MainPointPosVec = Span.PointPosVec(:,2:Span.NumPoints + 1);
CumDist = 0;
Flag = 0;
for MainPointNum = 1:Span.NumPoints
 if MainPointNum == 1
  CumDist = CumDist + dot(MainPointPosVec(:,MainPointNum) - Span.PointA_PosVec,MainSpanHorUnitVec);
 else
  CumDist = CumDist + dot(MainPointPosVec(:,MainPointNum) - MainPointPosVec(:,MainPointNum - 1),MainSpanHorUnitVec);
 end
 if InterDist <= CumDist
  if MainPointNum == 1
   PointNum_A = 0;
   PointNum_B = MainPointNum;
   PointA_PosVec = Span.PointA_PosVec;
   PointB_PosVec = MainPointPosVec(:,MainPointNum);
  else
   PointNum_A = MainPointNum - 1;
   PointNum_B = MainPointNum;
   PointA_PosVec = MainPointPosVec(:,MainPointNum - 1);
   PointB_PosVec = MainPointPosVec(:,MainPointNum);
  end
  Flag = 1;
  break;
 end
end
if Flag < 1
 PointNum_A = Span.NumPoints;
 PointNum_B = 0;
 PointA_PosVec = MainPointPosVec(:,Span.NumPoints);
 PointB_PosVec = Span.PointB_PosVec;
end
DistRatio = (InterDist - dot(PointA_PosVec - Span.PointA_PosVec,MainSpanHorUnitVec))/dot(PointB_PosVec - PointA_PosVec,MainSpanHorUnitVec);
Point_PosVec = PointA_PosVec + DistRatio*(PointB_PosVec - PointA_PosVec);
