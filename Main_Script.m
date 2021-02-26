%% Remaining tasks
% Full system coupling
% GUI
% Database
% Effective diameter calculation
% Parallel drag component
% Verify drag moment direction
% Specify anchor positions relative to poles (using distance and angle)

%% Setup workspace
clear;
clc;
close all;
format long g;

%% Read input parameters from text file
% Open input file
[FileName,FilePath] = uigetfile('.txt','Select input parameters file...');
InFile = fopen(fullfile(FilePath,FileName),'r');

% Solver properties
Solver.CableElmLength = 1;
Solver.PoleElmLength = 2;
Solver.MaxNumCablePoints = 100;
Solver.MaxNumPolePoints = 20;
Solver.FsolveOptions = optimoptions('fsolve',...
 'Display','final',...
 'FunctionTolerance',1e-12,...
 'MaxFunctionEvaluations',10000,...
 'MaxIterations',10000,...
 'OptimalityTolerance',1e-12,...
 'UseParallel',false);

% Global properties
fgetl(InFile);fgetl(InFile);
GlobalPars = fscanf(InFile,'%f');
Global.IceThickness = GlobalPars(1);
Global.WindSpeed = GlobalPars(2);
Global.WindDirection = GlobalPars(3);
Global.AmbientTemp = GlobalPars(4);
Global.GravAccel = 9.81;
Global.AirDensity = 1.2;
Global.IceDensity = 1000;
Global.ParDragCoeff = 1;
Global.PerpDragCoeff = 0.2;
Global.CylDragCoeff = 1.2;
Global.BaseTemp = 20;
Global.WindSpeedX = Global.WindSpeed*cosd(Global.WindDirection);
Global.WindSpeedY = Global.WindSpeed*sind(Global.WindDirection);
Global.WindVelVec = [
 Global.WindSpeedX;
 Global.WindSpeedY;
 0];
clear GlobalPars;

% Pole properties
fgetl(InFile);fgetl(InFile);
PolePars = fscanf(InFile,'%f');
PolePars = reshape(PolePars,8,[])';
System.NumPoles = size(PolePars,1);
for PoleNum = System.NumPoles:-1:1
 Pole(PoleNum).PosX = PolePars(PoleNum,2);
 Pole(PoleNum).PosY = PolePars(PoleNum,3);
 Pole(PoleNum).PosZ = PolePars(PoleNum,4);
 Pole(PoleNum).Height = PolePars(PoleNum,5);
 Pole(PoleNum).BaseDia = PolePars(PoleNum,6);
 Pole(PoleNum).TopDia = PolePars(PoleNum,7);
 Pole(PoleNum).ElasticMod = PolePars(PoleNum,8);
 Pole(PoleNum).BaseRad = Pole(PoleNum).BaseDia/2;
 Pole(PoleNum).TopRad = Pole(PoleNum).TopDia/2;
end
clear PolePars;

% Line properties
fgetl(InFile);fgetl(InFile);
LinePars = fscanf(InFile,'%f');
if isempty(LinePars)
 System.NumLines = 0;
else
 LinePars = reshape(LinePars,10,[])';
 System.NumLines = size(LinePars,1);
 for LineNum = System.NumLines:-1:1
  Line(LineNum).Tension = LinePars(LineNum,2);
  Line(LineNum).BaseSpWeight = LinePars(LineNum,3);
  Line(LineNum).Dia = LinePars(LineNum,4);
  Line(LineNum).Stiffness = LinePars(LineNum,5);
  Line(LineNum).ThermExpCoeff = LinePars(LineNum,6);
  Line(LineNum).InterLineNum = LinePars(LineNum,7);
  Line(LineNum).InterPoleNum_A = LinePars(LineNum,8);
  Line(LineNum).InterPoleNum_B = LinePars(LineNum,9);
  Line(LineNum).InterDist = LinePars(LineNum,10);
  Line(LineNum).SpWeight = Line(LineNum).BaseSpWeight./(1 + Line(LineNum).ThermExpCoeff*(Global.AmbientTemp - Global.BaseTemp));
 end
end
clear LinePars;

% Line-pole attachments
fgetl(InFile);fgetl(InFile);
if System.NumLines > 0
 LinePoleAttach = fscanf(InFile,'%f');
 LinePoleAttach = reshape(LinePoleAttach,System.NumLines + 1,[])';
 System.MaxNumAttachPerLine =  size(LinePoleAttach,1);
 for LineNum = System.NumLines:-1:1
  for LineAttachNum = System.MaxNumAttachPerLine:-1:1
   Line(LineNum).Attach(LineAttachNum).PoleNum = LinePoleAttach(LineAttachNum,LineNum + 1);
  end
 end
end
clear LinePoleAttach;

% Line attachment heights
fgetl(InFile);fgetl(InFile);
if System.NumLines > 0
 LineAttachHeight = fscanf(InFile,'%f');
 LineAttachHeight = reshape(LineAttachHeight,System.NumLines + 1,[])';
 for LineNum = System.NumLines:-1:1
  for LineAttachNum = System.MaxNumAttachPerLine:-1:1
   Line(LineNum).Attach(LineAttachNum).Height = LineAttachHeight(LineAttachNum,LineNum + 1);
  end
 end
end
clear LineAttachHeight;

% Lashing information
fgetl(InFile);fgetl(InFile);
LashPars = fscanf(InFile,'%f');
if isempty(LashPars)
 System.NumLash = 0;
else
 LashPars = reshape(LashPars,6,[])';
 System.NumLash = size(LashPars,1);
 for LashNum = System.NumLash:-1:1
  Lash(LashNum).SpWeight = LashPars(LashNum,2);
  Lash(LashNum).Dia = LashPars(LashNum,3);
  Lash(LashNum).LineNum = LashPars(LashNum,4);
  Lash(LashNum).PoleNum_A = LashPars(LashNum,5);
  Lash(LashNum).PoleNum_B = LashPars(LashNum,6);
  end
end
clear LashPars;

% Heavy equipment information
fgetl(InFile);fgetl(InFile);
EquipPars = fscanf(InFile,'%f');
if isempty(EquipPars)
 System.NumEquip = 0;
else
 EquipPars = reshape(EquipPars,6,[])';
 System.NumEquip = size(EquipPars,1);
 for EquipNum = System.NumEquip:-1:1
  Equip(EquipNum).Mass = EquipPars(EquipNum,2);
  Equip(EquipNum).LineNum = EquipPars(EquipNum,3);
  Equip(EquipNum).PoleNum_A = EquipPars(EquipNum,4);
  Equip(EquipNum).PoleNum_B = EquipPars(EquipNum,5);
  Equip(EquipNum).MountDist = EquipPars(EquipNum,6);
  end
end
clear EquipPars;

% Guy properties
fgetl(InFile);fgetl(InFile);
GuyPars = fscanf(InFile,'%f');
if isempty(GuyPars)
 System.NumGuys = 0;
else
 GuyPars = reshape(GuyPars,10,[])';
 System.NumGuys = size(GuyPars,1);
 for GuyNum = System.NumGuys:-1:1
  Guy(GuyNum).Tension = GuyPars(GuyNum,2);
  Guy(GuyNum).BaseSpWeight = GuyPars(GuyNum,3);
  Guy(GuyNum).Dia = GuyPars(GuyNum,4);
  Guy(GuyNum).Stiffness = GuyPars(GuyNum,5);
  Guy(GuyNum).ThermExpCoeff = GuyPars(GuyNum,6);
  Guy(GuyNum).PoleNum = GuyPars(GuyNum,7);
  Guy(GuyNum).AttachHeight = GuyPars(GuyNum,8);
  Guy(GuyNum).AnchNum = GuyPars(GuyNum,9);
  Guy(GuyNum).SideBarNum = GuyPars(GuyNum,10);
  Guy(GuyNum).SpWeight = Guy(GuyNum).BaseSpWeight./(1 + Guy(GuyNum).ThermExpCoeff*(Global.AmbientTemp - Global.BaseTemp));
 end
end
clear GuyPars;

% Anchor properties
fgetl(InFile);fgetl(InFile);
AnchPars = fscanf(InFile,'%f');
if isempty(AnchPars)
 System.NumAnch = 0;
else
 AnchPars = reshape(AnchPars,4,[])';
 System.NumAnch = size(AnchPars,1);
 for AnchNum = System.NumAnch:-1:1
  Anch(AnchNum).PosX = AnchPars(AnchNum,2);
  Anch(AnchNum).PosY = AnchPars(AnchNum,3);
  Anch(AnchNum).PosZ = AnchPars(AnchNum,4);
  end
end
clear AnchPars;

% Sidebar properties
fgetl(InFile);fgetl(InFile);
SideBarPars = fscanf(InFile,'%f');
if isempty(SideBarPars)
 System.NumSideBars = 0;
 SideBar = struct([]);
else
 SideBarPars = reshape(SideBarPars,7,[])';
 System.NumSideBars = size(SideBarPars,1);
 for SideBarNum = System.NumSideBars:-1:1
  SideBar(SideBarNum).PoleNum = SideBarPars(SideBarNum,2);
  SideBar(SideBarNum).MountHeight = SideBarPars(SideBarNum,3);
  SideBar(SideBarNum).Length = SideBarPars(SideBarNum,4);
  SideBar(SideBarNum).Angle = SideBarPars(SideBarNum,5);
  SideBar(SideBarNum).Dia = SideBarPars(SideBarNum,6);
  SideBar(SideBarNum).ElasticMod = SideBarPars(SideBarNum,7);
  SideBar(SideBarNum).AMOI = pi/4*((SideBar(SideBarNum).Dia/2)^4);
 end
end
clear SideBarPars;

% Close text file
fclose(InFile);

%% Pole-line information arrays
% Number of attachments per line
for LineNum = System.NumLines:-1:1
 Line(LineNum).NumAttach = 0;
 for LineAttachNum = 1:System.MaxNumAttachPerLine
  if Line(LineNum).Attach(LineAttachNum).PoleNum > 0
   Line(LineNum).NumAttach = Line(LineNum).NumAttach + 1;
  end
 end
 if Line(LineNum).InterLineNum > 0
  Line(LineNum).NumAttach = Line(LineNum).NumAttach + 1;
 end
end

% Number of spans per line
Tmp = zeros(System.NumLines,1);
for LineNum = System.NumLines:-1:1
 Line(LineNum).NumSpans = Line(LineNum).NumAttach - 1;
 Tmp = Line(LineNum).NumSpans;
end
System.MaxNumSpansPerLine = max(Tmp);
clear Tmp;

% Pole attachment line information
for PoleNum = System.NumPoles:-1:1
 Pole(PoleNum).NumLineAttach = 0;
end
Tmp = zeros(System.NumPoles,1);
for PoleNum = 1:System.NumPoles
 for LineNum = 1:System.NumLines
  for LineAttachNum = 1:Line(LineNum).NumAttach
   if Line(LineNum).Attach(LineAttachNum).PoleNum == PoleNum
    Pole(PoleNum).NumLineAttach = Pole(PoleNum).NumLineAttach + 1;
    Tmp(PoleNum) = Pole(PoleNum).NumLineAttach;
   end
  end
 end
end
System.MaxNumLineAttachPerPole = max(Tmp);
for PoleNum = System.NumPoles:-1:1
 for PoleLineAttachNum = System.MaxNumLineAttachPerPole:-1:1
  Pole(PoleNum).Attach(PoleLineAttachNum).LineNum = 0;
  Pole(PoleNum).Attach(PoleLineAttachNum).LineAttachNum = 0;
 end
end
for PoleNum = 1:System.NumPoles
 PoleLineAttachNum = 1;
 for LineNum = 1:System.NumLines
  for LineAttachNum = 1:Line(LineNum).NumAttach
   if Line(LineNum).Attach(LineAttachNum).PoleNum == PoleNum
    Pole(PoleNum).Attach(PoleLineAttachNum).LineNum = LineNum;
    Pole(PoleNum).Attach(PoleLineAttachNum).LineAttachNum = LineAttachNum;
    PoleLineAttachNum = PoleLineAttachNum + 1;
   end
  end
 end
end

%% Line intersection information arrays
% Span number of intersected line
for LineNum = System.NumLines:-1:1
 Line(LineNum).InterLineSpanNum = 0;
end
for LineNum = 1:System.NumLines
 if Line(LineNum).InterLineNum > 0
  MainLineNum = Line(LineNum).InterLineNum;
  MainLineSpanNum = 1;
  for MainLineAttachNum = 1:Line(MainLineNum).NumAttach - 1
   MainPoleNum_A = Line(MainLineNum).Attach(MainLineAttachNum).PoleNum;
   MainPoleNum_B = Line(MainLineNum).Attach(MainLineAttachNum + 1).PoleNum;
   if (MainPoleNum_A == Line(LineNum).InterPoleNum_A) && (MainPoleNum_B == Line(LineNum).InterPoleNum_B)
    Line(LineNum).InterLineSpanNum = MainLineSpanNum;
   end
   MainLineSpanNum = MainLineSpanNum + 1;
  end
 end
end

% Lines and spans that intersect any given span
for LineNum = System.NumLines:-1:1
 for LineSpanNum = Line(LineNum).NumSpans:-1:1
  Line(LineNum).Span(LineSpanNum).NumInterBy = 0;
 end
end
Tmp = zeros(System.NumLines,System.MaxNumSpansPerLine);
for LineNum = 1:System.NumLines
 LineSpanNum = 1;
 for LineAttachNum = 1:Line(LineNum).NumAttach - 1
  for InterLineNum = 1:System.NumLines
   if (Line(InterLineNum).InterLineNum == LineNum) && (Line(InterLineNum).InterLineSpanNum == LineSpanNum)
    Line(LineNum).Span(LineSpanNum).NumInterBy = Line(LineNum).Span(LineSpanNum).NumInterBy + 1;
   end
  end
  Tmp(LineNum,LineSpanNum) = Line(LineNum).Span(LineSpanNum).NumInterBy;
  LineSpanNum = LineSpanNum + 1;
 end
end
System.MaxNumInterBy = max(max(Tmp));
clear Tmp;
for LineNum = System.NumLines:-1:1
 for LineSpanNum = Line(LineNum).NumSpans:-1:1
  Line(LineNum).Span(LineSpanNum).InterByLineNum = zeros(1,System.MaxNumInterBy);
 end
end
for LineNum = 1:System.NumLines
 LineSpanNum = 1;
 for LineAttachNum = 1:Line(LineNum).NumAttach - 1
  InterByNum = 1;
  for InterLineNum = 1:System.NumLines
   if (Line(InterLineNum).InterLineNum == LineNum) && (Line(InterLineNum).InterLineSpanNum == LineSpanNum)
    Line(LineNum).Span(LineSpanNum).InterByLineNum(InterByNum) = InterLineNum;
    InterByNum = InterByNum + 1;
   end
  end
  LineSpanNum = LineSpanNum + 1;
 end
end

%% Line-Lashing information arrays
for LashNum = System.NumLash:-1:1
 Lash(LashNum).LineSpanNum_A = 0;
 Lash(LashNum).LineSpanNum_B = 0;
end
for LashNum = 1:System.NumLash
 LineNum = Lash(LashNum).LineNum;
 LineSpanNum = 1;
 for LineAttachNum = 1:Line(LineNum).NumAttach
  PoleNum = Line(LineNum).Attach(LineAttachNum).PoleNum;
  if (LineAttachNum == Line(LineNum).NumAttach) && (Line(LineNum).InterLineNum > 0)
   Lash(LashNum).LineSpanNum_B = LineSpanNum - 1;
  elseif PoleNum == Lash(LashNum).PoleNum_B
   Lash(LashNum).LineSpanNum_B = LineSpanNum - 1;
  elseif PoleNum == Lash(LashNum).PoleNum_A
   Lash(LashNum).LineSpanNum_A = LineSpanNum;
  end
  LineSpanNum = LineSpanNum + 1;
 end
end

%% Equipment information arrays
% Equipment span number
for EquipNum = System.NumEquip:-1:1
 Equip(EquipNum).LineSpanNum = 0;
end
for EquipNum = 1:System.NumEquip
 LineNum = Equip(EquipNum).LineNum;
 LineSpanNum = 1;
 for LineAttachNum = 1:Line(LineNum).NumAttach - 1
  PoleNum = Line(LineNum).Attach(LineAttachNum).PoleNum;
  if PoleNum == Equip(EquipNum).PoleNum_A
   Equip(EquipNum).LineSpanNum = LineSpanNum;
  end
  LineSpanNum = LineSpanNum + 1;
 end
end

% Span equipment information
for LineNum = System.NumLines:-1:1
 for LineSpanNum = System.MaxNumSpansPerLine:-1:1
  Line(LineNum).Span(LineSpanNum).NumEquip = 0;
 end
end
Tmp = zeros(System.NumLines,System.MaxNumSpansPerLine);
for LineNum = 1:System.NumLines
 for LineSpanNum = 1:Line(LineNum).NumSpans
  for EquipNum = 1:System.NumEquip
   if (Equip(EquipNum).LineNum == LineNum) && (Equip(EquipNum).LineSpanNum == LineSpanNum)
    Line(LineNum).Span(LineSpanNum).NumEquip = Line(LineNum).Span(LineSpanNum).NumEquip + 1;
   end
  end
  Tmp(LineNum,LineSpanNum) = Line(LineNum).Span(LineSpanNum).NumEquip;
 end
end
System.MaxNumEquipPerLineSpan = max(max(Tmp));
for LineNum = System.NumLines:-1:1
 for LineSpanNum = System.MaxNumSpansPerLine:-1:1
  for LineSpanEquipNum = System.MaxNumEquipPerLineSpan:-1:1
   Line(LineNum).Span(LineSpanNum).Equip(LineSpanEquipNum).EquipNum = 0;
   Line(LineNum).Span(LineSpanNum).Equip(LineSpanEquipNum).Weight = 0;
   Line(LineNum).Span(LineSpanNum).Equip(LineSpanEquipNum).MountDist = 0;
  end
 end
end
for LineNum = 1:System.NumLines
 for LineSpanNum = 1:Line(LineNum).NumSpans
  LineSpanEquipNum = 1;
  for EquipNum = 1:System.NumEquip
   if (Equip(EquipNum).LineNum == LineNum) && (Equip(EquipNum).LineSpanNum == LineSpanNum)
    Line(LineNum).Span(LineSpanNum).Equip(LineSpanEquipNum).EquipNum = EquipNum;
    Line(LineNum).Span(LineSpanNum).Equip(LineSpanEquipNum).Weight = Global.GravAccel*Equip(EquipNum).Mass;
    Line(LineNum).Span(LineSpanNum).Equip(LineSpanEquipNum).MountDist = Equip(EquipNum).MountDist;
    LineSpanEquipNum = LineSpanEquipNum + 1;
   end
  end
 end
end

%% Guy information arrays
% Number of spans per guy
Tmp = ones(System.NumGuys,1);
for GuyNum = System.NumGuys:-1:1
 if Guy(GuyNum).SideBarNum > 0
  Tmp(GuyNum) = 2;
 end
 Guy(GuyNum).NumSpans = Tmp(GuyNum);
end
System.MaxNumSpansPerGuy = max(Tmp);
clear Tmp;

%% Pole-guy information arrays
for PoleNum = System.NumPoles:-1:1
 Pole(PoleNum).NumGuys = 0;
end
Tmp = zeros(System.NumPoles,1);
for PoleNum = 1:System.NumPoles
 for GuyNum = 1:System.NumGuys
  if Guy(GuyNum).PoleNum == PoleNum
   Pole(PoleNum).NumGuys = Pole(PoleNum).NumGuys + 1;
   Tmp(PoleNum) = Pole(PoleNum).NumGuys;
  end
 end
end
System.MaxNumGuysPerPole = max(Tmp);
for PoleNum = System.NumPoles:-1:1
 for PoleGuyNum = System.MaxNumGuysPerPole:-1:1
  Pole(PoleNum).Guy(PoleGuyNum).GuyNum = 0;
 end
end
for PoleNum = 1:System.NumPoles
 PoleGuyNum = 1;
 for GuyNum = 1:System.NumGuys
  if Guy(GuyNum).PoleNum == PoleNum
   Pole(PoleNum).Guy(PoleGuyNum).GuyNum = GuyNum;
   PoleGuyNum = PoleGuyNum + 1;
  end
 end
end

%% Pole-sidebar information arrays
for PoleNum = System.NumPoles:-1:1
 Pole(PoleNum).NumSideBars = 0;
end
Tmp = zeros(System.NumPoles,1);
for PoleNum = 1:System.NumPoles
 for SideBarNum = 1:System.NumSideBars
  if SideBar(SideBarNum).PoleNum == PoleNum
   Pole(PoleNum).NumSideBars = Pole(PoleNum).NumSideBars + 1;
   Tmp(PoleNum) = Pole(PoleNum).NumSideBars;
  end
 end
end
System.MaxNumSideBarsPerPole = max(Tmp);
for PoleNum = System.NumPoles:-1:1
 for PoleSideBarNum = System.MaxNumSideBarsPerPole:-1:1
  Pole(PoleNum).SideBar(PoleSideBarNum).SideBarNum = 0;
 end
end
for PoleNum = 1:System.NumPoles
 PoleSideBarNum = 1;
 for SideBarNum = 1:System.NumSideBars
  if SideBar(SideBarNum).PoleNum == PoleNum
   Pole(PoleNum).SideBar(PoleSideBarNum).SideBarNum = SideBarNum;
   PoleSideBarNum = PoleSideBarNum + 1;
  end
 end
end

%% Guy-sidebar information arrays
for SideBarNum = System.NumSideBars:-1:1
 SideBar(SideBarNum).NumGuys = 0;
end
Tmp = zeros(System.NumSideBars,1);
for SideBarNum = 1:System.NumSideBars
 for GuyNum = 1:System.NumGuys
  if Guy(GuyNum).SideBarNum == SideBarNum
   SideBar(SideBarNum).NumGuys = SideBar(SideBarNum).NumGuys + 1;
   Tmp(SideBarNum) = SideBar(SideBarNum).NumGuys;
  end
 end
end
System.MaxNumGuysPerSideBar = max(Tmp);
for SideBarNum = System.NumSideBars:-1:1
 for SideBarGuyNum = System.MaxNumGuysPerSideBar:-1:1
  SideBar(SideBarNum).Guy(SideBarGuyNum).GuyNum = 0;
 end
end
for SideBarNum = 1:System.NumSideBars
 SideBarGuyNum = 1;
 for GuyNum = 1:System.NumGuys
  if Guy(GuyNum).SideBarNum == SideBarNum
   SideBar(SideBarNum).Guy(SideBarGuyNum).GuyNum = GuyNum;
   SideBarGuyNum = SideBarGuyNum + 1;
  end
 end
end

%% Initialize arrays
% Line properties
for LineNum = System.NumLines:-1:1
 for LineSpanNum = System.MaxNumSpansPerLine:-1:1
  Line(LineNum).Span(LineSpanNum).Angle_A = 0;
  Line(LineNum).Span(LineSpanNum).HorForce = 0;
  Line(LineNum).Span(LineSpanNum).VerForce_A = 0;
  Line(LineNum).Span(LineSpanNum).VerForce_B = 0;
  Line(LineNum).Span(LineSpanNum).Weight = 0;
  Line(LineNum).Span(LineSpanNum).BaseLength = 0;
  Line(LineNum).Span(LineSpanNum).NumElm = 0;
  Line(LineNum).Span(LineSpanNum).NumPointsTotal = 0;
  Line(LineNum).Span(LineSpanNum).NumPoints = 0;
  Line(LineNum).Span(LineSpanNum).PointWeight = 0;
  Line(LineNum).Span(LineSpanNum).ElmLength = 0;
  Line(LineNum).Span(LineSpanNum).ElmStiff = 0;
  Line(LineNum).Span(LineSpanNum).Dia = Line(LineNum).Dia + 2*Global.IceThickness;
  Line(LineNum).Span(LineSpanNum).PointA_PosVec = zeros(3,1);
  Line(LineNum).Span(LineSpanNum).PointB_PosVec = zeros(3,1);
  Line(LineNum).Span(LineSpanNum).PointPosVec = zeros(3,Solver.MaxNumCablePoints);
  Line(LineNum).Span(LineSpanNum).InterDist = 0;
  Line(LineNum).Span(LineSpanNum).SpringForceA = zeros(3,Solver.MaxNumCablePoints - 2);
  Line(LineNum).Span(LineSpanNum).SpringForceB = zeros(3,Solver.MaxNumCablePoints - 2);
 end
end

% Guy properties
for GuyNum = System.NumGuys:-1:1
 for GuySpanNum = System.MaxNumSpansPerGuy:-1:1
  Guy(GuyNum).Span(GuySpanNum).Angle_A = 0;
  Guy(GuyNum).Span(GuySpanNum).HorForce = 0;
  Guy(GuyNum).Span(GuySpanNum).VerForce_A = 0;
  Guy(GuyNum).Span(GuySpanNum).VerForce_B = 0;
  Guy(GuyNum).Span(GuySpanNum).Weight = 0;
  Guy(GuyNum).Span(GuySpanNum).BaseLength = 0;
  Guy(GuyNum).Span(GuySpanNum).NumElm = 0;
  Guy(GuyNum).Span(GuySpanNum).NumPointsTotal = 0;
  Guy(GuyNum).Span(GuySpanNum).NumPoints = 0;
  Guy(GuyNum).Span(GuySpanNum).PointWeight = 0;
  Guy(GuyNum).Span(GuySpanNum).ElmLength = 0;
  Guy(GuyNum).Span(GuySpanNum).ElmStiff = 0;
  Guy(GuyNum).Span(GuySpanNum).Dia = Guy(GuyNum).Dia;
  Guy(GuyNum).Span(GuySpanNum).PointA_PosVec = zeros(3,1);
  Guy(GuyNum).Span(GuySpanNum).PointB_PosVec = zeros(3,1);
  Guy(GuyNum).Span(GuySpanNum).PointPosVec = zeros(3,Solver.MaxNumCablePoints);
  Guy(GuyNum).Span(GuySpanNum).SpringForceA = zeros(3,Solver.MaxNumCablePoints - 2);
  Guy(GuyNum).Span(GuySpanNum).SpringForceB = zeros(3,Solver.MaxNumCablePoints - 2);
 end
end

% Sidebar arrays
for SideBarNum = System.NumSideBars:-1:1
 SideBar(SideBarNum).BasePosVec = zeros(3,1);
 SideBar(SideBarNum).TipPosVec = zeros(3,1);
 SideBar(SideBarNum).TipDisp = 0;
 SideBar(SideBarNum).RelTipDef = 0;
end

%% Mesh poles
% Initialize arrays
for PoleNum = System.NumPoles:-1:1
 Pole(PoleNum).NumElm = 0;
 Pole(PoleNum).ElmLength = 0;
 Pole(PoleNum).NumPointsTotal = 0;
 Pole(PoleNum).NumPoints = 0;
 Pole(PoleNum).PointDispVec = zeros(2,Solver.MaxNumPolePoints);
 Pole(PoleNum).PointRotVec = zeros(2,Solver.MaxNumPolePoints);
 Pole(PoleNum).RefPointDispVec = zeros(2,Solver.MaxNumPolePoints);
 Pole(PoleNum).RefPointRotVec = zeros(2,Solver.MaxNumPolePoints);
 Pole(PoleNum).PointPosVec = zeros(3,Solver.MaxNumPolePoints);
 Pole(PoleNum).PointForceVec = zeros(2,Solver.MaxNumPolePoints);
 Pole(PoleNum).PointMomVec = zeros(2,Solver.MaxNumPolePoints);
 Pole(PoleNum).StiffMat = zeros(2*Solver.MaxNumPolePoints,2*Solver.MaxNumPolePoints);
end

% FEM setup
for PoleNum = 1:System.NumPoles
 % Pole discretization
 Pole(PoleNum).NumElm = ceil(Pole(PoleNum).Height/Solver.PoleElmLength);
 Pole(PoleNum).ElmLength = Pole(PoleNum).Height/Pole(PoleNum).NumElm;
 Pole(PoleNum).NumPointsTotal = Pole(PoleNum).NumElm + 1;
 Pole(PoleNum).NumPoints = Pole(PoleNum).NumPointsTotal - 1;
 Pole(PoleNum).PointPosVec(:,1:Pole(PoleNum).NumPointsTotal) = [
  Pole(PoleNum).PosX*ones(1,Pole(PoleNum).NumPointsTotal);
  Pole(PoleNum).PosY*ones(1,Pole(PoleNum).NumPointsTotal);
  (0:Pole(PoleNum).ElmLength:Pole(PoleNum).Height)];
 
 % FEM stiffnes matrix
 for PoleElmNum = 1:Pole(PoleNum).NumElm
  ElmStiffMat = ElmStiffMatFunc(Pole(PoleNum).ElmLength,Pole(PoleNum).ElasticMod,Pole(PoleNum).BaseRad,Pole(PoleNum).TopRad,Pole(PoleNum).PointPosVec(3,PoleElmNum),Pole(PoleNum).Height);
  Pole(PoleNum).StiffMat(2*PoleElmNum - 1:2*PoleElmNum + 2,2*PoleElmNum - 1:2*PoleElmNum + 2) = Pole(PoleNum).StiffMat(2*PoleElmNum - 1:2*PoleElmNum + 2,2*PoleElmNum - 1:2*PoleElmNum + 2) + ElmStiffMat;
 end
end

% Total number of pole points for simulation
TotalNumPolePoints = 0;
for PoleNum = 1:System.NumPoles
 TotalNumPolePoints = TotalNumPolePoints + Pole(PoleNum).NumPoints;
end

%% Wind load vectors on poles
for PoleNum = System.NumPoles:-1:1
 Pole(PoleNum).WindLoadVecX = zeros(2*Solver.MaxNumPolePoints,1);
 Pole(PoleNum).WindLoadVecY = zeros(2*Solver.MaxNumPolePoints,1);
end
for PoleNum = 1:System.NumPoles
 for PoleElmNum = 1:Pole(PoleNum).NumElm
  ElmWindLoadVecX = [
   (Global.CylDragCoeff*Pole(PoleNum).ElmLength*Global.WindSpeedX*Global.WindSpeed*Pole(PoleNum).BaseRad*Global.AirDensity)/2 - (Global.CylDragCoeff*Pole(PoleNum).ElmLength*Global.WindSpeedX*Global.WindSpeed*Global.AirDensity*(Pole(PoleNum).BaseRad - Pole(PoleNum).TopRad)*(3*Pole(PoleNum).ElmLength + 10*Pole(PoleNum).PointPosVec(3,PoleElmNum)))/(20*Pole(PoleNum).Height)
   (Global.CylDragCoeff*Pole(PoleNum).ElmLength^2*Global.WindSpeedX*Global.WindSpeed*Global.AirDensity*(5*Pole(PoleNum).Height*Pole(PoleNum).BaseRad - 2*Pole(PoleNum).ElmLength*Pole(PoleNum).BaseRad + 2*Pole(PoleNum).ElmLength*Pole(PoleNum).TopRad - 5*Pole(PoleNum).BaseRad*Pole(PoleNum).PointPosVec(3,PoleElmNum) + 5*Pole(PoleNum).TopRad*Pole(PoleNum).PointPosVec(3,PoleElmNum)))/(60*Pole(PoleNum).Height)
   (Global.CylDragCoeff*Pole(PoleNum).ElmLength*Global.WindSpeedX*Global.WindSpeed*Global.AirDensity*(10*Pole(PoleNum).Height*Pole(PoleNum).BaseRad - 7*Pole(PoleNum).ElmLength*Pole(PoleNum).BaseRad + 7*Pole(PoleNum).ElmLength*Pole(PoleNum).TopRad - 10*Pole(PoleNum).BaseRad*Pole(PoleNum).PointPosVec(3,PoleElmNum) + 10*Pole(PoleNum).TopRad*Pole(PoleNum).PointPosVec(3,PoleElmNum)))/(20*Pole(PoleNum).Height)
   -(Global.CylDragCoeff*Pole(PoleNum).ElmLength^2*Global.WindSpeedX*Global.WindSpeed*Global.AirDensity*(5*Pole(PoleNum).Height*Pole(PoleNum).BaseRad - 3*Pole(PoleNum).ElmLength*Pole(PoleNum).BaseRad + 3*Pole(PoleNum).ElmLength*Pole(PoleNum).TopRad - 5*Pole(PoleNum).BaseRad*Pole(PoleNum).PointPosVec(3,PoleElmNum) + 5*Pole(PoleNum).TopRad*Pole(PoleNum).PointPosVec(3,PoleElmNum)))/(60*Pole(PoleNum).Height)];
  ElmWindLoadVecY = [
   (Global.CylDragCoeff*Pole(PoleNum).ElmLength*Global.WindSpeedY*Global.WindSpeed*Pole(PoleNum).BaseRad*Global.AirDensity)/2 - (Global.CylDragCoeff*Pole(PoleNum).ElmLength*Global.WindSpeedY*Global.WindSpeed*Global.AirDensity*(Pole(PoleNum).BaseRad - Pole(PoleNum).TopRad)*(3*Pole(PoleNum).ElmLength + 10*Pole(PoleNum).PointPosVec(3,PoleElmNum)))/(20*Pole(PoleNum).Height)
   (Global.CylDragCoeff*Pole(PoleNum).ElmLength^2*Global.WindSpeedY*Global.WindSpeed*Global.AirDensity*(5*Pole(PoleNum).Height*Pole(PoleNum).BaseRad - 2*Pole(PoleNum).ElmLength*Pole(PoleNum).BaseRad + 2*Pole(PoleNum).ElmLength*Pole(PoleNum).TopRad - 5*Pole(PoleNum).BaseRad*Pole(PoleNum).PointPosVec(3,PoleElmNum) + 5*Pole(PoleNum).TopRad*Pole(PoleNum).PointPosVec(3,PoleElmNum)))/(60*Pole(PoleNum).Height)
   (Global.CylDragCoeff*Pole(PoleNum).ElmLength*Global.WindSpeedY*Global.WindSpeed*Global.AirDensity*(10*Pole(PoleNum).Height*Pole(PoleNum).BaseRad - 7*Pole(PoleNum).ElmLength*Pole(PoleNum).BaseRad + 7*Pole(PoleNum).ElmLength*Pole(PoleNum).TopRad - 10*Pole(PoleNum).BaseRad*Pole(PoleNum).PointPosVec(3,PoleElmNum) + 10*Pole(PoleNum).TopRad*Pole(PoleNum).PointPosVec(3,PoleElmNum)))/(20*Pole(PoleNum).Height)
   -(Global.CylDragCoeff*Pole(PoleNum).ElmLength^2*Global.WindSpeedY*Global.WindSpeed*Global.AirDensity*(5*Pole(PoleNum).Height*Pole(PoleNum).BaseRad - 3*Pole(PoleNum).ElmLength*Pole(PoleNum).BaseRad + 3*Pole(PoleNum).ElmLength*Pole(PoleNum).TopRad - 5*Pole(PoleNum).BaseRad*Pole(PoleNum).PointPosVec(3,PoleElmNum) + 5*Pole(PoleNum).TopRad*Pole(PoleNum).PointPosVec(3,PoleElmNum)))/(60*Pole(PoleNum).Height)];
  Pole(PoleNum).WindLoadVecX(2*PoleElmNum - 1:2*PoleElmNum + 2) = Pole(PoleNum).WindLoadVecX(2*PoleElmNum - 1:2*PoleElmNum + 2) + ElmWindLoadVecX;
  Pole(PoleNum).WindLoadVecY(2*PoleElmNum - 1:2*PoleElmNum + 2) = Pole(PoleNum).WindLoadVecY(2*PoleElmNum - 1:2*PoleElmNum + 2) + ElmWindLoadVecY;
 end
end

%% Solve catenary problems for guy wires
for PoleNum = 1:System.NumPoles
 InitSideBarVars = zeros(Pole(PoleNum).NumSideBars,1);
 for PoleSideBarNum = 1:Pole(PoleNum).NumSideBars
  SideBarNum = Pole(PoleNum).SideBar(PoleSideBarNum).SideBarNum;
  InitSideBarVars(PoleSideBarNum) = SideBar(SideBarNum).RelTipDef;
 end
 InitPoleVars = [
  Pole(PoleNum).PointDispVec(:,2:Pole(PoleNum).NumPointsTotal);
  Pole(PoleNum).PointRotVec(:,2:Pole(PoleNum).NumPointsTotal)];
 InitVars = [
  InitSideBarVars;
  reshape(InitPoleVars,[],1)];
 Func = @(Vars) GuyAssemblyFunc(Vars,Pole(PoleNum),Guy,Anch,SideBar,Solver);
 Vars = fsolve(Func,InitVars,Solver.FsolveOptions);
 [~,Pole(PoleNum),Guy,SideBar] = GuyAssemblyFunc(Vars,Pole(PoleNum),Guy,Anch,SideBar,Solver);
end

%% Solve catenary problems for support wires
while 1
 for PoleNum = 1:System.NumPoles
  InitSideBarVars = zeros(Pole(PoleNum).NumSideBars,1);
  for PoleSideBarNum = 1:Pole(PoleNum).NumSideBars
   SideBarNum = Pole(PoleNum).SideBar(PoleSideBarNum).SideBarNum;
   InitSideBarVars(PoleSideBarNum) = SideBar(SideBarNum).RelTipDef;
  end
  InitPoleVars = [
   Pole(PoleNum).PointDispVec(:,2:Pole(PoleNum).NumPointsTotal);
   Pole(PoleNum).PointRotVec(:,2:Pole(PoleNum).NumPointsTotal)];
  InitVars = [
   InitSideBarVars;
   reshape(InitPoleVars,[],1)];
  Func = @(Vars) AssemblyFunc(Vars,Pole,Line,Guy,Anch,SideBar,Solver,PoleNum,System);
  Vars = fsolve(Func,InitVars,Solver.FsolveOptions);
  [~,Pole,Line,Guy,SideBar] = AssemblyFunc(Vars,Pole,Line,Guy,Anch,SideBar,Solver,PoleNum,System);
 end
%  if (all(abs(Pole(PoleNum).PointDispVec(:,2:Pole(PoleNum).NumPointsTotal) - Pole(PoleNum).RefPointDispVec(:,2:Pole(PoleNum).NumPointsTotal)) < 1e-6,'all')) && (all(abs(Pole(PoleNum).PointRotVec(:,2:Pole(PoleNum).NumPointsTotal) - Pole(PoleNum).RefPointRotVec(:,2:Pole(PoleNum).NumPointsTotal)) < 1e-6,'all'))
 if (norm(abs(Pole(PoleNum).PointDispVec(:,2:Pole(PoleNum).NumPointsTotal) - Pole(PoleNum).RefPointDispVec(:,2:Pole(PoleNum).NumPointsTotal))) < 1e-6) && (norm(abs(Pole(PoleNum).PointRotVec(:,2:Pole(PoleNum).NumPointsTotal) - Pole(PoleNum).RefPointRotVec(:,2:Pole(PoleNum).NumPointsTotal))) < 1e-6)
  break;
 end
 Pole(PoleNum).RefPointDispVec(:,2:Pole(PoleNum).NumPointsTotal) = Pole(PoleNum).PointDispVec(:,2:Pole(PoleNum).NumPointsTotal);
 Pole(PoleNum).RefPointRotVec(:,2:Pole(PoleNum).NumPointsTotal) = Pole(PoleNum).PointRotVec(:,2:Pole(PoleNum).NumPointsTotal);
end

%% Add lashing, icing, heavy equipment, and temperature effects to lines
% Expand cable segment lengths due to temperature rise
for LineNum = 1:System.NumLines
 for LineSpanNum = 1:Line(LineNum).NumSpans
  Line(LineNum).Span(LineSpanNum).ElmLength = (1 + Line(LineNum).ThermExpCoeff*(Global.AmbientTemp - Global.BaseTemp))*Line(LineNum).Span(LineSpanNum).ElmLength;
 end
end

% Add weight and diameter from lashings
for LashNum = 1:System.NumLash
 LineNum = Lash(LashNum).LineNum;
 for LineSpanNum = Lash(LashNum).LineSpanNum_A:Lash(LashNum).LineSpanNum_B
  Line(LineNum).Span(LineSpanNum).Weight = Line(LineNum).Span(LineSpanNum).Weight + Lash(LashNum).SpWeight*Line(LineNum).Span(LineSpanNum).BaseLength;
  Line(LineNum).Span(LineSpanNum).Dia = Line(LineNum).Span(LineSpanNum).Dia + Lash(LashNum).Dia + 2*Global.IceThickness;
 end
end

% Add weight from icing
for LineNum = System.NumLines:-1:1
 for LineSpanNum = System.MaxNumSpansPerLine:-1:1
  Line(LineNum).Span(LineSpanNum).IceArea = 0;
 end
end
for LineNum = 1:System.NumLines
 for LineSpanNum = 1:Line(LineNum).NumSpans
  Line(LineNum).Span(LineSpanNum).IceArea = Line(LineNum).Span(LineSpanNum).IceArea + pi/4*((Line(LineNum).Dia + 2*Global.IceThickness)^2 - Line(LineNum).Dia^2);
 end
end
for LashNum = 1:System.NumLash
 LineNum = Lash(LashNum).LineNum;
 for LineSpanNum = Lash(LashNum).LineSpanNum_A:Lash(LashNum).LineSpanNum_B
  Line(LineNum).Span(LineSpanNum).IceArea = Line(LineNum).Span(LineSpanNum).IceArea + pi/4*((Lash(LashNum).Dia + 2*Global.IceThickness)^2 - Lash(LashNum).Dia^2);
 end
end
for LineNum = 1:System.NumLines
 for LineSpanNum = 1:Line(LineNum).NumSpans
  Line(LineNum).Span(LineSpanNum).Weight = Line(LineNum).Span(LineSpanNum).Weight + Global.IceDensity*Global.GravAccel*Line(LineNum).Span(LineSpanNum).IceArea*Line(LineNum).Span(LineSpanNum).BaseLength;
  Line(LineNum).Span(LineSpanNum).PointWeight = Line(LineNum).Span(LineSpanNum).Weight/Line(LineNum).Span(LineSpanNum).NumPointsTotal;
 end
end

% Add point weight due to heavy equipment
for LineNum = System.NumLines:-1:1
 for LineSpanNum = System.MaxNumSpansPerLine:-1:1
  Line(LineNum).Span(LineSpanNum).EquipPointWeightVec = zeros(3,Solver.MaxNumCablePoints - 2);
 end
end
for EquipNum = 1:System.NumEquip
 LineNum = Equip(EquipNum).LineNum;
 LineSpanNum = Equip(EquipNum).LineSpanNum;
 [DistRatio,PointNum_A,PointNum_B,~,~,~] = InterFunc(Line(LineNum).Span(LineSpanNum),Equip(EquipNum).MountDist);
 if PointNum_A == 0
  Line(LineNum).Span(LineSpanNum).EquipPointWeightVec(3,PointNum_B) = -DistRatio*Global.GravAccel*Equip(EquipNum).Mass;
 elseif PointNum_B == 0
  Line(LineNum).Span(LineSpanNum).EquipPointWeightVec(3,PointNum_A) = -(1 - DistRatio)*Global.GravAccel*Equip(EquipNum).Mass;
 else
  Line(LineNum).Span(LineSpanNum).EquipPointWeightVec(3,PointNum_A) = -(1 - DistRatio)*Global.GravAccel*Equip(EquipNum).Mass;
  Line(LineNum).Span(LineSpanNum).EquipPointWeightVec(3,PointNum_B) = -DistRatio*Global.GravAccel*Equip(EquipNum).Mass;
 end
end

%% Add icing and temperature effects to guys
% Expand cable segment lengths due to temperature rise
for GuyNum = 1:System.NumGuys
 for GuySpanNum = 1:Guy(GuyNum).NumSpans
  Guy(GuyNum).Span(GuySpanNum).ElmLength = (1 + Guy(GuyNum).ThermExpCoeff*(Global.AmbientTemp - Global.BaseTemp))*Guy(GuyNum).Span(GuySpanNum).ElmLength;
 end
end

% Add weight from icing
for GuyNum = System.NumGuys:-1:1
 for GuySpanNum = System.MaxNumSpansPerGuy:-1:1
  Guy(GuyNum).Span(GuySpanNum).IceArea = 0;
 end
end
for GuyNum = 1:System.NumGuys
 for GuySpanNum = 1:Guy(GuyNum).NumSpans
  Guy(GuyNum).Span(GuySpanNum).IceArea = Guy(GuyNum).Span(GuySpanNum).IceArea + pi/4*((Guy(GuyNum).Dia + 2*Global.IceThickness)^2 - Guy(GuyNum).Dia^2);
 end
end
for GuyNum = 1:System.NumGuys
 for GuySpanNum = 1:Guy(GuyNum).NumSpans
  Guy(GuyNum).Span(GuySpanNum).Weight = Guy(GuyNum).Span(GuySpanNum).Weight + Global.IceDensity*Global.GravAccel*Guy(GuyNum).Span(GuySpanNum).IceArea*Guy(GuyNum).Span(GuySpanNum).BaseLength;
  Guy(GuyNum).Span(GuySpanNum).PointWeight = Guy(GuyNum).Span(GuySpanNum).Weight/Guy(GuyNum).Span(GuySpanNum).NumPointsTotal;
 end
end

%% Solve final problems for all wires
while 1
 for PoleNum = 1:System.NumPoles
  InitSideBarVars = zeros(Pole(PoleNum).NumSideBars,1);
  for PoleSideBarNum = 1:Pole(PoleNum).NumSideBars
   SideBarNum = Pole(PoleNum).SideBar(PoleSideBarNum).SideBarNum;
   InitSideBarVars(PoleSideBarNum) = SideBar(SideBarNum).RelTipDef;
  end
  InitPoleVars = [
   Pole(PoleNum).PointDispVec(:,2:Pole(PoleNum).NumPointsTotal);
   Pole(PoleNum).PointRotVec(:,2:Pole(PoleNum).NumPointsTotal)];
  InitVars = [
   InitSideBarVars;
   reshape(InitPoleVars,[],1)];
  Func = @(Vars) FinalSimFunc(Global,Vars,Pole,Line,Guy,Anch,SideBar,Solver,PoleNum,System);
  Vars = fsolve(Func,InitVars,Solver.FsolveOptions);
  [~,Pole,Line,Guy,SideBar] = FinalSimFunc(Global,Vars,Pole,Line,Guy,Anch,SideBar,Solver,PoleNum,System);
 end
%  if (all(abs(Pole(PoleNum).PointDispVec(:,2:Pole(PoleNum).NumPointsTotal) - Pole(PoleNum).RefPointDispVec(:,2:Pole(PoleNum).NumPointsTotal)) < 1e-6,'all')) && (all(abs(Pole(PoleNum).PointRotVec(:,2:Pole(PoleNum).NumPointsTotal) - Pole(PoleNum).RefPointRotVec(:,2:Pole(PoleNum).NumPointsTotal)) < 1e-6,'all'))
 if (norm(abs(Pole(PoleNum).PointDispVec(:,2:Pole(PoleNum).NumPointsTotal) - Pole(PoleNum).RefPointDispVec(:,2:Pole(PoleNum).NumPointsTotal))) < 1e-6) && (norm(abs(Pole(PoleNum).PointRotVec(:,2:Pole(PoleNum).NumPointsTotal) - Pole(PoleNum).RefPointRotVec(:,2:Pole(PoleNum).NumPointsTotal))) < 1e-6)
  break;
 end
 Pole(PoleNum).RefPointDispVec(:,2:Pole(PoleNum).NumPointsTotal) = Pole(PoleNum).PointDispVec(:,2:Pole(PoleNum).NumPointsTotal);
 Pole(PoleNum).RefPointRotVec(:,2:Pole(PoleNum).NumPointsTotal) = Pole(PoleNum).PointRotVec(:,2:Pole(PoleNum).NumPointsTotal);
end

%% Plot
figure;
hold on;
LineSpanNum = 1;
for LineNum = 1:System.NumLines
 LineSpanNum = 1;
 for LineAttachNum = 1:Line(LineNum).NumAttach - 1
  plot3(...
   Line(LineNum).Span(LineSpanNum).PointPosVec(1,1:Line(LineNum).Span(LineSpanNum).NumPointsTotal),...
   Line(LineNum).Span(LineSpanNum).PointPosVec(2,1:Line(LineNum).Span(LineSpanNum).NumPointsTotal),...
   Line(LineNum).Span(LineSpanNum).PointPosVec(3,1:Line(LineNum).Span(LineSpanNum).NumPointsTotal),...
   'k','LineWidth',1.5);
  LineSpanNum = LineSpanNum + 1;
 end
end
for GuyNum = 1:System.NumGuys
 for GuySpanNum = 1:Guy(GuyNum).NumSpans
  plot3(...
   Guy(GuyNum).Span(GuySpanNum).PointPosVec(1,1:Guy(GuyNum).Span(GuySpanNum).NumPointsTotal),...
   Guy(GuyNum).Span(GuySpanNum).PointPosVec(2,1:Guy(GuyNum).Span(GuySpanNum).NumPointsTotal),...
   Guy(GuyNum).Span(GuySpanNum).PointPosVec(3,1:Guy(GuyNum).Span(GuySpanNum).NumPointsTotal),...
   'k','LineWidth',1.5);
 end
end
for AnchNum = 1:System.NumAnch
 plot3(Anch(AnchNum).PosX,Anch(AnchNum).PosY,Anch(AnchNum).PosZ,'kx','MarkerSize',6,'LineWidth',1);
end
for PoleNum = 1:System.NumPoles
 plot3(...
  Pole(PoleNum).PointPosVec(1,1:Pole(PoleNum).NumPointsTotal),...
  Pole(PoleNum).PointPosVec(2,1:Pole(PoleNum).NumPointsTotal),...
  Pole(PoleNum).PosZ + Pole(PoleNum).PointPosVec(3,1:Pole(PoleNum).NumPointsTotal),...
  'k','LineWidth',3.5);
end
Tmp = zeros(3,System.NumPoles);
for PoleNum = 1:System.NumPoles
 Tmp(:,PoleNum) = [
  Pole(PoleNum).PosX;
  Pole(PoleNum).PosY;
  Pole(PoleNum).PosZ];
end
for SideBarNum = 1:System.NumSideBars
 PoleNum = SideBar(SideBarNum).PoleNum;
 plot3(...
  [SideBar(SideBarNum).BasePosVec(1);SideBar(SideBarNum).TipPosVec(1)],...
  [SideBar(SideBarNum).BasePosVec(2);SideBar(SideBarNum).TipPosVec(2)],...
  [SideBar(SideBarNum).BasePosVec(3);SideBar(SideBarNum).TipPosVec(3)],...
  'k','LineWidth',2);
end
for EquipNum = 1:System.NumEquip
 LineNum = Equip(EquipNum).LineNum;
 LineSpanNum = Equip(EquipNum).LineSpanNum;
 [~,~,~,~,~,Point_PosVec] = InterFunc(Line(LineNum).Span(LineSpanNum),Equip(EquipNum).MountDist);
 plot3(Point_PosVec(1),Point_PosVec(2),Point_PosVec(3) + 0.5,'kv','MarkerSize',6,'LineWidth',1.5);
end
axis([...
 min(Tmp(1,:)) - 10, max(Tmp(1,:)) + 10,...
 min(Tmp(2,:)) - 10, max(Tmp(2,:)) + 10,...
 -inf, inf]);
daspect([1 1 1]);
view(3);
box on;
hold off;