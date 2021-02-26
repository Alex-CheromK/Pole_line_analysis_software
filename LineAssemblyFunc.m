function Line = LineAssemblyFunc(Solver,System,Pole,Line)
%% Solve catenary problems for support wires
LineFlag = ones(System.NumLines,1);
while any(LineFlag) > 0
 % Line span tensions
 for LineNum = 1:System.NumLines
  LineSpanNum = 1;
  if (LineFlag(LineNum) > 0) && ((Line(LineNum).InterLineNum < 1) || ((Line(LineNum).InterLineNum > 0) && (LineFlag(Line(LineNum).InterLineNum) < 1)))
   for LineAttachNum = 1:Line(LineNum).NumAttach - 1
    if (LineSpanNum == Line(LineNum).NumSpans) && (Line(LineNum).InterLineNum > 0)
     % Update span endpoint positions
    PoleNum_A = Line(LineNum).Attach(LineAttachNum).PoleNum;
    Line(LineNum).Span(LineSpanNum).PointA_PosVec = [
     interp1(Pole(PoleNum_A).PointPosVec(3,1:Pole(PoleNum_A).NumPointsTotal),Pole(PoleNum_A).PointPosVec(1,1:Pole(PoleNum_A).NumPointsTotal),Line(LineNum).Attach(LineAttachNum).Height);
     interp1(Pole(PoleNum_A).PointPosVec(3,1:Pole(PoleNum_A).NumPointsTotal),Pole(PoleNum_A).PointPosVec(2,1:Pole(PoleNum_A).NumPointsTotal),Line(LineNum).Attach(LineAttachNum).Height);
     Pole(PoleNum_A).PosZ + Line(LineNum).Attach(LineAttachNum).Height];
     
     % Build main span structure
     MainLineNum = Line(LineNum).InterLineNum;
     MainLineSpanNum = Line(LineNum).InterLineSpanNum;
     MainSpan = Line(MainLineNum).Span(MainLineSpanNum);
     
     % Build other pre-solved intersecting span structures
     clear InterSpan;
     NumInterSpan_Active = 0;
     for InterByNum = 1:Line(MainLineNum).Span(MainLineSpanNum).NumInterBy
      InterLineNum = Line(MainLineNum).Span(MainLineSpanNum).InterByLineNum(InterByNum);
      if (InterLineNum ~= LineNum) && (LineFlag(InterLineNum) < 1)
       NumInterSpan_Active = NumInterSpan_Active + 1;
      end
     end
     InterLineNum_Active = zeros(1,NumInterSpan_Active);
     InterSpanNum = 1;
     for InterByNum = 1:Line(MainLineNum).Span(MainLineSpanNum).NumInterBy
      InterLineNum = Line(MainLineNum).Span(MainLineSpanNum).InterByLineNum(InterByNum);
      if (InterLineNum ~= LineNum) && (LineFlag(InterLineNum) < 1)
       InterLineNum_Active(InterSpanNum) = InterLineNum;
       InterSpanNum = InterSpanNum + 1;
      end
     end
     if NumInterSpan_Active > 0
      for InterSpanNum = NumInterSpan_Active:-1:1
       InterLineNum = InterLineNum_Active(InterSpanNum);
       InterLineSpanNum = Line(InterLineNum).NumSpans;
       InterSpan(InterSpanNum) = Line(InterLineNum).Span(InterLineSpanNum);
       InterSpan(InterSpanNum).InterDist = Line(InterLineNum).InterDist;
      end
     else
      InterSpan = struct([]);
     end
     
     % Build catenary span structure
     CatSpan.InterDist = Line(LineNum).InterDist;
     PoleNum_A = Line(LineNum).Attach(LineAttachNum).PoleNum;
     CatSpan.PointA_PosVec = Line(LineNum).Span(LineSpanNum).PointA_PosVec;
     if LineAttachNum == 1
      CatSpan.Tension_A = Line(LineNum).Tension;
     else
      CatSpan.Tension_A = sqrt(Line(LineNum).Span(LineSpanNum - 1).HorForce^2 + Line(LineNum).Span(LineSpanNum - 1).VerForce_B^2);
     end
     CatSpan.SpWeight = Line(LineNum).BaseSpWeight;
     CatSpan.Stiffness = Line(LineNum).Stiffness;
     
     % Define initial guesses for combined catenary and LP problem
     CatSpanAngle_A_Init = -10*pi/180;
     CatSpanVerForce_B_Init = Line(LineNum).BaseSpWeight*norm(MainSpan.PointA_PosVec - CatSpan.PointA_PosVec)/2;
     MainSpanPointPosVec_Init = MainSpan.PointPosVec(:,2:MainSpan.NumPointsTotal - 1);
     TotalNumInterSpanPoints = 0;
     if NumInterSpan_Active > 0
      for InterSpanNum = 1:NumInterSpan_Active
       TotalNumInterSpanPoints = TotalNumInterSpanPoints + InterSpan(InterSpanNum).NumPoints;
      end
      InterSpanPointPosVec_Init = zeros(3,TotalNumInterSpanPoints);
      RefNumPoints = 0;
      for InterSpanNum = 1:NumInterSpan_Active
       InterSpanPointPosVec_Init(:,RefNumPoints + 1:RefNumPoints + InterSpan(InterSpanNum).NumPoints) = InterSpan(InterSpanNum).PointPosVec(:,2:InterSpan(InterSpanNum).NumPointsTotal - 1);
       RefNumPoints = RefNumPoints + InterSpan(InterSpanNum).NumPoints;
      end
     end
     if NumInterSpan_Active > 0
      Vars_Init = [
       CatSpanAngle_A_Init;
       CatSpanVerForce_B_Init;
       reshape(MainSpanPointPosVec_Init,[],1)
       reshape(InterSpanPointPosVec_Init,[],1)];
     else
      Vars_Init = [
       CatSpanAngle_A_Init;
       CatSpanVerForce_B_Init;
       reshape(MainSpanPointPosVec_Init,[],1)];
     end
     
     % Solve combined catenary and LP problem
     Func = @(Vars) LP_Cable2(Vars,MainSpan,InterSpan,CatSpan,NumInterSpan_Active,TotalNumInterSpanPoints,Solver);
     Vars = fsolve(Func,Vars_Init,Solver.FsolveOptions);
     Line(LineNum).Span(LineSpanNum).Angle_A = Vars(1);
     Line(LineNum).Span(LineSpanNum).HorForce = CatSpan.Tension_A*cos(Line(LineNum).Span(LineSpanNum).Angle_A);
     Line(LineNum).Span(LineSpanNum).VerForce_A = CatSpan.Tension_A*sin(Line(LineNum).Span(LineSpanNum).Angle_A);
     Line(LineNum).Span(LineSpanNum).VerForce_B = Vars(2);
     Vars(1:2) = [];
     Vars = reshape(Vars,3,[]);
     Line(MainLineNum).Span(MainLineSpanNum).PointPosVec(:,1:Line(MainLineNum).Span(MainLineSpanNum).NumPointsTotal) = [Line(MainLineNum).Span(MainLineSpanNum).PointA_PosVec,Vars(:,1:Line(MainLineNum).Span(MainLineSpanNum).NumPoints),Line(MainLineNum).Span(MainLineSpanNum).PointB_PosVec];
     MainSpan = Line(MainLineNum).Span(MainLineSpanNum); % Update main span information
     Vars(:,1:Line(MainLineNum).Span(MainLineSpanNum).NumPoints) = [];
     if NumInterSpan_Active > 0
      for InterSpanNum = 1:NumInterSpan_Active
       InterLineNum = InterLineNum_Active(InterSpanNum);
       InterLineSpanNum = Line(InterLineNum).NumSpans;
       [~,~,~,~,~,Point_PosVec] = InterFunc(MainSpan,Line(InterLineNum).InterDist);
       Line(InterLineNum).Span(InterLineSpanNum).PointPosVec(:,1:Line(InterLineNum).Span(InterLineSpanNum).NumPointsTotal) = [Line(InterLineNum).Span(InterLineSpanNum).PointA_PosVec,Vars(:,1:Line(InterLineNum).Span(InterLineSpanNum).NumPoints),Point_PosVec];
       Vars(:,1:Line(InterLineNum).Span(InterLineSpanNum).NumPoints) = [];
      end
     end
     
     % Find length-constrained properties for catenary span
     Line(LineNum).Span(LineSpanNum).Weight = Line(LineNum).Span(LineSpanNum).VerForce_B - Line(LineNum).Span(LineSpanNum).VerForce_A;
     Line(LineNum).Span(LineSpanNum).BaseLength = Line(LineNum).Span(LineSpanNum).Weight/Line(LineNum).BaseSpWeight;
     
     % Discretize catenary span
     Line(LineNum).Span(LineSpanNum).NumElm = ceil(Line(LineNum).Span(LineSpanNum).BaseLength/Solver.CableElmLength);
     Line(LineNum).Span(LineSpanNum).NumPointsTotal = Line(LineNum).Span(LineSpanNum).NumElm + 1;
     Line(LineNum).Span(LineSpanNum).NumPoints = Line(LineNum).Span(LineSpanNum).NumPointsTotal - 2;
     Line(LineNum).Span(LineSpanNum).PointWeight = Line(LineNum).Span(LineSpanNum).Weight/Line(LineNum).Span(LineSpanNum).NumPointsTotal;
     Line(LineNum).Span(LineSpanNum).ElmLength = Line(LineNum).Span(LineSpanNum).BaseLength/Line(LineNum).Span(LineSpanNum).NumElm;
     Line(LineNum).Span(LineSpanNum).ElmStiff = Line(LineNum).Stiffness/Line(LineNum).Span(LineSpanNum).ElmLength;
     
     % Catenary span endpoint locations
     Line(LineNum).Span(LineSpanNum).PointA_PosVec = Line(LineNum).Span(LineSpanNum).PointA_PosVec;
     [~,~,~,~,~,Point_PosVec] = InterFunc(MainSpan,Line(LineNum).InterDist);
     Line(LineNum).Span(LineSpanNum).PointB_PosVec = Point_PosVec;
     
     % Build initial catenary span point position vector
     LineSpanHorVec = Line(LineNum).Span(LineSpanNum).PointB_PosVec(1:2) - Line(LineNum).Span(LineSpanNum).PointA_PosVec(1:2);
     LineSpanHorDist = norm(LineSpanHorVec);
     HeightA = Line(LineNum).Span(LineSpanNum).PointA_PosVec(3);
     HeightB = Line(LineNum).Span(LineSpanNum).PointB_PosVec(3);
     a = Line(LineNum).Span(LineSpanNum).HorForce/Line(LineNum).SpWeight;
     A = exp(-LineSpanHorDist/a) - 1;
     B = -2*(HeightB - HeightA)/a;
     C = exp(LineSpanHorDist/a) - 1;
     Lambda = (-B - sqrt(B^2 - 4*A*C))/(2*A);
     if Lambda <= 0
      error('Solution to catenary profile is a complex number!');
     end
     k1 = a*log(Lambda);
     k2 = HeightA - a*cosh(k1/a);
     LineSpanHorAxisVec = linspace(0,LineSpanHorDist,Line(LineNum).Span(LineSpanNum).NumPointsTotal)';
     LineSpanHorUnitVec = LineSpanHorVec/LineSpanHorDist;
     LineSpanHorPointA_PosVec = Line(LineNum).Span(LineSpanNum).PointA_PosVec(1:2);
     LineSpanPointPosVecInit = [
      LineSpanHorPointA_PosVec + LineSpanHorAxisVec'.*LineSpanHorUnitVec;
      a*cosh((LineSpanHorAxisVec' - k1)/a) + k2];
     Line(LineNum).Span(LineSpanNum).PointPosVec(:,1:Line(LineNum).Span(LineSpanNum).NumPointsTotal) = LineSpanPointPosVecInit;
     
     % Update intersecting span structure
     clear InterSpan;
     NumInterSpan_Active = 0;
     for InterByNum = 1:Line(MainLineNum).Span(MainLineSpanNum).NumInterBy
      InterLineNum = Line(MainLineNum).Span(MainLineSpanNum).InterByLineNum(InterByNum);
      if (InterLineNum == LineNum) || (LineFlag(InterLineNum) < 1)
       NumInterSpan_Active = NumInterSpan_Active + 1;
      end
     end
     InterLineNum_Active = zeros(1,NumInterSpan_Active);
     InterSpanNum = 1;
     for InterByNum = 1:Line(MainLineNum).Span(MainLineSpanNum).NumInterBy
      InterLineNum = Line(MainLineNum).Span(MainLineSpanNum).InterByLineNum(InterByNum);
      if (InterLineNum == LineNum) || (LineFlag(InterLineNum) < 1)
       InterLineNum_Active(InterSpanNum) = InterLineNum;
       InterSpanNum = InterSpanNum + 1;
      end
     end
     for InterSpanNum = NumInterSpan_Active:-1:1
      InterLineNum = InterLineNum_Active(InterSpanNum);
      InterLineSpanNum = Line(InterLineNum).NumSpans;
      InterSpan(InterSpanNum) = Line(InterLineNum).Span(InterLineSpanNum);
      InterSpan(InterSpanNum).InterDist = Line(InterLineNum).InterDist;
     end
     
     % Define initial guesses for full LP problem
     MainSpanPointPosVec_Init = MainSpan.PointPosVec(:,2:MainSpan.NumPointsTotal - 1);
     TotalNumInterSpanPoints = 0;
     for InterSpanNum = 1:NumInterSpan_Active
      TotalNumInterSpanPoints = TotalNumInterSpanPoints + InterSpan(InterSpanNum).NumPoints;
     end
     InterSpanPointPosVec_Init = zeros(3,TotalNumInterSpanPoints);
     RefNumPoints = 0;
     for InterSpanNum = 1:NumInterSpan_Active
      InterSpanPointPosVec_Init(:,RefNumPoints + 1:RefNumPoints + InterSpan(InterSpanNum).NumPoints) = InterSpan(InterSpanNum).PointPosVec(:,2:InterSpan(InterSpanNum).NumPointsTotal - 1);
      RefNumPoints = RefNumPoints + InterSpan(InterSpanNum).NumPoints;
     end
     Vars_Init = [
      reshape(MainSpanPointPosVec_Init,[],1)
      reshape(InterSpanPointPosVec_Init,[],1)];
     
     % Solve full LP problem
     Func = @(Vars) LP_Cable3(Vars,MainSpan,InterSpan,NumInterSpan_Active,TotalNumInterSpanPoints,Solver);
     Vars = fsolve(Func,Vars_Init,Solver.FsolveOptions);
     Vars = reshape(Vars,3,[]);
     Line(MainLineNum).Span(MainLineSpanNum).PointPosVec(:,1:Line(MainLineNum).Span(MainLineSpanNum).NumPointsTotal) = [Line(MainLineNum).Span(MainLineSpanNum).PointA_PosVec,Vars(:,1:Line(MainLineNum).Span(MainLineSpanNum).NumPoints),Line(MainLineNum).Span(MainLineSpanNum).PointB_PosVec];
     MainSpan = Line(MainLineNum).Span(MainLineSpanNum); % Update main span information
     Vars(:,1:Line(MainLineNum).Span(MainLineSpanNum).NumPoints) = [];
     for InterSpanNum = 1:NumInterSpan_Active
      InterLineNum = InterLineNum_Active(InterSpanNum);
      InterLineSpanNum = Line(InterLineNum).NumSpans;
      [~,~,~,~,~,Point_PosVec] = InterFunc(MainSpan,Line(InterLineNum).InterDist);
      Line(InterLineNum).Span(InterLineSpanNum).PointPosVec(:,1:Line(InterLineNum).Span(InterLineSpanNum).NumPointsTotal) = [Line(InterLineNum).Span(InterLineSpanNum).PointA_PosVec,Vars(:,1:Line(InterLineNum).Span(InterLineSpanNum).NumPoints),Point_PosVec];
      Vars(:,1:Line(InterLineNum).Span(InterLineSpanNum).NumPoints) = [];
     end
    else
     % Update span endpoint positions
    PoleNum_A = Line(LineNum).Attach(LineAttachNum).PoleNum;
    PoleNum_B = Line(LineNum).Attach(LineAttachNum + 1).PoleNum;
    Line(LineNum).Span(LineSpanNum).PointA_PosVec = [
     interp1(Pole(PoleNum_A).PointPosVec(3,1:Pole(PoleNum_A).NumPointsTotal),Pole(PoleNum_A).PointPosVec(1,1:Pole(PoleNum_A).NumPointsTotal),Line(LineNum).Attach(LineAttachNum).Height);
     interp1(Pole(PoleNum_A).PointPosVec(3,1:Pole(PoleNum_A).NumPointsTotal),Pole(PoleNum_A).PointPosVec(2,1:Pole(PoleNum_A).NumPointsTotal),Line(LineNum).Attach(LineAttachNum).Height);
     Pole(PoleNum_A).PosZ + Line(LineNum).Attach(LineAttachNum).Height];
    Line(LineNum).Span(LineSpanNum).PointB_PosVec = [
     interp1(Pole(PoleNum_B).PointPosVec(3,1:Pole(PoleNum_B).NumPointsTotal),Pole(PoleNum_B).PointPosVec(1,1:Pole(PoleNum_A).NumPointsTotal),Line(LineNum).Attach(LineAttachNum + 1).Height);
     interp1(Pole(PoleNum_B).PointPosVec(3,1:Pole(PoleNum_B).NumPointsTotal),Pole(PoleNum_B).PointPosVec(2,1:Pole(PoleNum_A).NumPointsTotal),Line(LineNum).Attach(LineAttachNum + 1).Height);
     Pole(PoleNum_B).PosZ + Line(LineNum).Attach(LineAttachNum + 1).Height];
     
     % Solve catenary problem
     LineSpanHorVec = [
      Line(LineNum).Span(LineSpanNum).PointB_PosVec(1) - Line(LineNum).Span(LineSpanNum).PointA_PosVec(1);
      Line(LineNum).Span(LineSpanNum).PointB_PosVec(2) - Line(LineNum).Span(LineSpanNum).PointA_PosVec(2)];
     LineSpanHorDist = norm(LineSpanHorVec);
     LineSpanVerDist = Line(LineNum).Span(LineSpanNum).PointB_PosVec(3) - Line(LineNum).Span(LineSpanNum).PointA_PosVec(3);
     if LineAttachNum == 1
      SpanTension_A = Line(LineNum).Tension;
     else
      SpanTension_A = sqrt(Line(LineNum).Span(LineSpanNum - 1).HorForce^2 + Line(LineNum).Span(LineSpanNum - 1).VerForce_B^2);
     end
     SpanAngle_A_Init = -10*pi/180;
     SpanVerForce_B_Init = Line(LineNum).BaseSpWeight*norm([LineSpanHorDist LineSpanVerDist]')/2;
     Func = @(Vars) CatFunc_Tension(Vars(1),Vars(2),SpanTension_A,Line(LineNum).BaseSpWeight,Line(LineNum).Stiffness,LineSpanHorDist,LineSpanVerDist);
     Vars = fsolve(Func,[SpanAngle_A_Init SpanVerForce_B_Init]',Solver.FsolveOptions);
     Line(LineNum).Span(LineSpanNum).Angle_A = Vars(1);
     Line(LineNum).Span(LineSpanNum).HorForce = SpanTension_A*cos(Line(LineNum).Span(LineSpanNum).Angle_A);
     Line(LineNum).Span(LineSpanNum).VerForce_A = SpanTension_A*sin(Line(LineNum).Span(LineSpanNum).Angle_A);
     Line(LineNum).Span(LineSpanNum).VerForce_B = Vars(2);
     
     % Find length-constrained span properties
     Line(LineNum).Span(LineSpanNum).Weight = Line(LineNum).Span(LineSpanNum).VerForce_B - Line(LineNum).Span(LineSpanNum).VerForce_A;
     Line(LineNum).Span(LineSpanNum).BaseLength = Line(LineNum).Span(LineSpanNum).Weight/Line(LineNum).BaseSpWeight;
     
     % Discretize span
     Line(LineNum).Span(LineSpanNum).NumElm = ceil(Line(LineNum).Span(LineSpanNum).BaseLength/Solver.CableElmLength);
     Line(LineNum).Span(LineSpanNum).NumPointsTotal = Line(LineNum).Span(LineSpanNum).NumElm + 1;
     Line(LineNum).Span(LineSpanNum).NumPoints = Line(LineNum).Span(LineSpanNum).NumPointsTotal - 2;
     Line(LineNum).Span(LineSpanNum).PointWeight = Line(LineNum).Span(LineSpanNum).Weight/Line(LineNum).Span(LineSpanNum).NumPointsTotal;
     Line(LineNum).Span(LineSpanNum).ElmLength = Line(LineNum).Span(LineSpanNum).BaseLength/Line(LineNum).Span(LineSpanNum).NumElm;
     Line(LineNum).Span(LineSpanNum).ElmStiff = Line(LineNum).Stiffness/Line(LineNum).Span(LineSpanNum).ElmLength;
     
     % Build initial span point position vector
     LineSpanHorVec = Line(LineNum).Span(LineSpanNum).PointB_PosVec(1:2) - Line(LineNum).Span(LineSpanNum).PointA_PosVec(1:2);
     LineSpanHorDist = norm(LineSpanHorVec);
     HeightA = Line(LineNum).Span(LineSpanNum).PointA_PosVec(3);
     HeightB = Line(LineNum).Span(LineSpanNum).PointB_PosVec(3);
     a = Line(LineNum).Span(LineSpanNum).HorForce/Line(LineNum).SpWeight;
     A = exp(-LineSpanHorDist/a) - 1;
     B = -2*(HeightB - HeightA)/a;
     C = exp(LineSpanHorDist/a) - 1;
     Lambda = (-B - sqrt(B^2 - 4*A*C))/(2*A);
     if Lambda <= 0
      error('Solution to catenary profile is a complex number!');
     end
     k1 = a*log(Lambda);
     k2 = HeightA - a*cosh(k1/a);
     LineSpanHorAxisVec = linspace(0,LineSpanHorDist,Line(LineNum).Span(LineSpanNum).NumPointsTotal)';
     LineSpanHorUnitVec = LineSpanHorVec/LineSpanHorDist;
     LineSpanHorPointA_PosVec = Line(LineNum).Span(LineSpanNum).PointA_PosVec(1:2);
     LineSpanPointPosVecInit = [
      LineSpanHorPointA_PosVec + LineSpanHorAxisVec'.*LineSpanHorUnitVec;
      a*cosh((LineSpanHorAxisVec' - k1)/a) + k2];
     LineSpanPointPosVecInit(:,1) = []; LineSpanPointPosVecInit(:,end) = [];
     InitVars = reshape(LineSpanPointPosVecInit,[],1);
     
     % Solve span LP problem
     Func = @(Vars) LP_Cable(Vars,Line(LineNum).Span(LineSpanNum));
     Vars = fsolve(Func,InitVars,Solver.FsolveOptions);
     Line(LineNum).Span(LineSpanNum).PointPosVec(:,1:Line(LineNum).Span(LineSpanNum).NumPointsTotal) = [Line(LineNum).Span(LineSpanNum).PointA_PosVec,reshape(Vars,3,[]),Line(LineNum).Span(LineSpanNum).PointB_PosVec];
    end
    
    % Go to next line span
    LineSpanNum = LineSpanNum + 1;
   end
   LineFlag(LineNum) = 0;
  end
 end
end