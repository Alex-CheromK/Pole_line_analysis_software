function Line = LineSimFunc(Solver,System,Pole,Line,Global)
%% 
for LineNum = 1:System.NumLines
 LineSpanNum = 1;
 for LineAttachNum = 1:Line(LineNum).NumAttach - 1
  if (LineSpanNum == Line(LineNum).NumSpans) && (Line(LineNum).InterLineNum > 0) % The current line span intersects another main span
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
   
   % Build intersecting span structures
   clear InterSpan;
   NumInterSpan_Active = Line(MainLineNum).Span(MainLineSpanNum).NumInterBy;
   InterLineNum_Active = zeros(1,NumInterSpan_Active);
   InterSpanNum = 1;
   for InterByNum = 1:Line(MainLineNum).Span(MainLineSpanNum).NumInterBy
    InterLineNum = Line(MainLineNum).Span(MainLineSpanNum).InterByLineNum(InterByNum);
    InterLineNum_Active(InterSpanNum) = InterLineNum;
    InterSpanNum = InterSpanNum + 1;
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
   Func = @(Vars) LP_Cable4(Vars,MainSpan,InterSpan,NumInterSpan_Active,TotalNumInterSpanPoints,Solver,Global);
   Vars = fsolve(Func,Vars_Init,Solver.FsolveOptions);
   [~,MainPointPosVec,MainSpringForceA,MainSpringForceB,InterPointPosVec,InterSpringForceA,InterSpringForceB] = LP_Cable4(Vars,MainSpan,InterSpan,NumInterSpan_Active,TotalNumInterSpanPoints,Solver,Global);
   Line(MainLineNum).Span(MainLineSpanNum).PointPosVec(:,1:Line(MainLineNum).Span(MainLineSpanNum).NumPointsTotal) = [Line(MainLineNum).Span(MainLineSpanNum).PointA_PosVec,MainPointPosVec,Line(MainLineNum).Span(MainLineSpanNum).PointB_PosVec];
   Line(MainLineNum).Span(MainLineSpanNum).SpringForceA(:,1:Line(MainLineNum).Span(MainLineSpanNum).NumPoints) = MainSpringForceA(:,1:Line(MainLineNum).Span(MainLineSpanNum).NumPoints);
   Line(MainLineNum).Span(MainLineSpanNum).SpringForceB(:,1:Line(MainLineNum).Span(MainLineSpanNum).NumPoints) = MainSpringForceB(:,1:Line(MainLineNum).Span(MainLineSpanNum).NumPoints);
   MainSpan = Line(MainLineNum).Span(MainLineSpanNum); % Update main span information
   for InterSpanNum = 1:NumInterSpan_Active
    InterLineNum = InterLineNum_Active(InterSpanNum);
    InterLineSpanNum = Line(InterLineNum).NumSpans;
    [~,~,~,~,~,Point_PosVec] = InterFunc(MainSpan,Line(InterLineNum).InterDist);
    Line(InterLineNum).Span(InterLineSpanNum).PointPosVec(:,1:Line(InterLineNum).Span(InterLineSpanNum).NumPointsTotal) = [Line(InterLineNum).Span(InterLineSpanNum).PointA_PosVec,InterPointPosVec(:,1:Line(InterLineNum).Span(InterLineSpanNum).NumPoints,InterSpanNum),Point_PosVec];
    Line(InterLineNum).Span(InterLineSpanNum).SpringForceA(:,1:Line(InterLineNum).Span(InterLineSpanNum).NumPoints) = InterSpringForceA(:,1:Line(InterLineNum).Span(InterLineSpanNum).NumPoints,InterLineSpanNum);
    Line(InterLineNum).Span(InterLineSpanNum).SpringForceB(:,1:Line(InterLineNum).Span(InterLineSpanNum).NumPoints) = InterSpringForceB(:,1:Line(InterLineNum).Span(InterLineSpanNum).NumPoints,InterLineSpanNum);
   end
  elseif Line(LineNum).Span(LineSpanNum).NumInterBy > 0 % The current line span is intersected by another span
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
   
   % Build main span structure
   MainLineNum = LineNum;
   MainLineSpanNum = LineSpanNum;
   MainSpan = Line(MainLineNum).Span(MainLineSpanNum);
   
   % Build intersecting span structures
   clear InterSpan;
   NumInterSpan_Active = Line(MainLineNum).Span(MainLineSpanNum).NumInterBy;
   InterLineNum_Active = zeros(1,NumInterSpan_Active);
   InterSpanNum = 1;
   for InterByNum = 1:Line(MainLineNum).Span(MainLineSpanNum).NumInterBy
    InterLineNum = Line(MainLineNum).Span(MainLineSpanNum).InterByLineNum(InterByNum);
    InterLineNum_Active(InterSpanNum) = InterLineNum;
    InterSpanNum = InterSpanNum + 1;
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
   Func = @(Vars) LP_Cable4(Vars,MainSpan,InterSpan,NumInterSpan_Active,TotalNumInterSpanPoints,Solver,Global);
   Vars = fsolve(Func,Vars_Init,Solver.FsolveOptions);
   [~,MainPointPosVec,MainSpringForceA,MainSpringForceB,InterPointPosVec,InterSpringForceA,InterSpringForceB] = LP_Cable4(Vars,MainSpan,InterSpan,NumInterSpan_Active,TotalNumInterSpanPoints,Solver,Global);
   Line(MainLineNum).Span(MainLineSpanNum).PointPosVec(:,1:Line(MainLineNum).Span(MainLineSpanNum).NumPointsTotal) = [Line(MainLineNum).Span(MainLineSpanNum).PointA_PosVec,MainPointPosVec,Line(MainLineNum).Span(MainLineSpanNum).PointB_PosVec];
   Line(MainLineNum).Span(MainLineSpanNum).SpringForceA(:,1:Line(MainLineNum).Span(MainLineSpanNum).NumPoints) = MainSpringForceA(:,1:Line(MainLineNum).Span(MainLineSpanNum).NumPoints);
   Line(MainLineNum).Span(MainLineSpanNum).SpringForceB(:,1:Line(MainLineNum).Span(MainLineSpanNum).NumPoints) = MainSpringForceB(:,1:Line(MainLineNum).Span(MainLineSpanNum).NumPoints);
   MainSpan = Line(MainLineNum).Span(MainLineSpanNum); % Update main span information
   for InterSpanNum = 1:NumInterSpan_Active
    InterLineNum = InterLineNum_Active(InterSpanNum);
    InterLineSpanNum = Line(InterLineNum).NumSpans;
    [~,~,~,~,~,Point_PosVec] = InterFunc(MainSpan,Line(InterLineNum).InterDist);
    Line(InterLineNum).Span(InterLineSpanNum).PointPosVec(:,1:Line(InterLineNum).Span(InterLineSpanNum).NumPointsTotal) = [Line(InterLineNum).Span(InterLineSpanNum).PointA_PosVec,InterPointPosVec(:,1:Line(InterLineNum).Span(InterLineSpanNum).NumPoints,InterSpanNum),Point_PosVec];
    Line(InterLineNum).Span(InterLineSpanNum).SpringForceA(:,1:Line(InterLineNum).Span(InterLineSpanNum).NumPoints) = InterSpringForceA(:,1:Line(InterLineNum).Span(InterLineSpanNum).NumPoints,InterLineSpanNum);
    Line(InterLineNum).Span(InterLineSpanNum).SpringForceB(:,1:Line(InterLineNum).Span(InterLineSpanNum).NumPoints) = InterSpringForceB(:,1:Line(InterLineNum).Span(InterLineSpanNum).NumPoints,InterLineSpanNum);
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
   
   % Define initial conditions
   InitVars = Line(LineNum).Span(LineSpanNum).PointPosVec(:,2:Line(LineNum).Span(LineSpanNum).NumPointsTotal - 1);
   InitVars = reshape(InitVars,[],1);
   
   % Solve span LP problem
   Func = @(Vars) LP_Cable5(Vars,Line(LineNum).Span(LineSpanNum),Global);
   Vars = fsolve(Func,InitVars,Solver.FsolveOptions);
   [~,PointPosVec,SpringForceA,SpringForceB] = LP_Cable5(Vars,Line(LineNum).Span(LineSpanNum),Global);
   Line(LineNum).Span(LineSpanNum).PointPosVec(:,1:Line(LineNum).Span(LineSpanNum).NumPointsTotal) = [Line(LineNum).Span(LineSpanNum).PointA_PosVec,PointPosVec,Line(LineNum).Span(LineSpanNum).PointB_PosVec];
   Line(LineNum).Span(LineSpanNum).SpringForceA(:,1:Line(LineNum).Span(LineSpanNum).NumPoints) = SpringForceA;
   Line(LineNum).Span(LineSpanNum).SpringForceB(:,1:Line(LineNum).Span(LineSpanNum).NumPoints) = SpringForceB;
  end
  
  % Move to next span
  LineSpanNum = LineSpanNum + 1;
 end
end