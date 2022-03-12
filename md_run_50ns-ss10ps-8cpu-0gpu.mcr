


ph=7.4

ions='Na,Cl,0.9'

temperature='310K'

density=0.99342

pressurectrl='SolventProbe,Name=HOH,Density=(density)'




format='sim'

extension=10

cellshape='Cube'

if !count speed
  speed='normal'

if !count duration
  duration=50000

if speed=='fast'
  saveinterval=250000
else  
  saveinterval=10000

ForceField AMBER14

Cutoff 8

Boundary periodic

Longrange Coulomb

CorrectDrift On

RandomSeed 1


RequireVersion 15.1.1

WarnIsError On

if MacroTarget==''
  RaiseError "This macro requires a target. Either edit the macro file or click Options > Macro > Set target to choose a target structure"

if runWithMacro and ConsoleMode
  Processors CPUThreads=8,GPU=None

Clear
Console off
waterscene = FileSize (MacroTarget)_water.sce
solventscene = FileSize (MacroTarget)_solvent.sce

if waterscene
  LoadSce (MacroTarget)_water
elif solventscene 
  LoadSce (MacroTarget)_solvent
else
  scene = FileSize (MacroTarget).sce
  if scene
    LoadSce (MacroTarget)
    simcell = CountObj SimCell
    if !simcell
      RaiseError 'If you provide a scene, it must contain a simulation cell, but none was found in (MacroTarget).sce'
  else
    for type in 'yob','pdb'
      size = FileSize (MacroTarget).(type)
      if size
        break
    if !size
      RaiseError 'Initial structure not found, expected (MacroTarget).pdb or .yob. Make sure to create a project directory and place the structure there'
    Load(type) (MacroTarget)
    Unselect
    NiceOriAll
    DelBond N,C,LenMin=5
    DelRes Water with 0 arrows to all
    CleanAll
    pH (ph)
    if Structure
      OptHydAll
    Cell Auto,Extension=(extension),Shape=(cellshape)
    SaveSce (MacroTarget)
  bnd = Boundary
  if bnd=='Wall'
    ShowMessage "The simulation cell you created has wall boundaries, which reduces the simulation accuracy due to boundary effects..."
    Wait ContinueButton
    ShowMessage "You can click 'Simulation > Cell boundaries > Periodic' now to correct the problem, or 'Continue' immediately..."
    Wait ContinueButton
  if ph!='None'
    bnd = Boundary
    if bnd=='Wall'
      FillCellWater
    else
      Experiment Neutralization
        WaterDensity (density)
        pH (ph)
        Ions (ions)
        pKaFile (MacroTarget).pka
        Speed Fast
      Experiment On
      Wait ExpEnd
  filename = '(MacroTarget)_solvent.yob'
  solventfound = FileSize (filename)
  if solventfound
    obj2 = ListObj Water
    DelRes Water
    obj1 = LoadYOb (MacroTarget)_solvent
    dens = PropObj (obj1)
    if !dens
      RaiseError 'Please load the solvent molecule (filename), click Edit > Number > Property value > Obj X, choose the solvent density, then save the file again'
    CleanObj (obj1)
    FillCellObj (obj1),Density=(dens),BumpSum=4,RandomOri=Yes
    if obj1!=obj2
      JoinObj (obj)
    NameObj (obj2),Solvent
    SaveSce (MacroTarget)_solvent
  else
    SaveSce (MacroTarget)_water
Unselect

if speed=='fast'
  FixBond all,Element H
  FixHydAngle all
  tslist=2,2.5
else
  FreeBond all,all
  FreeAngle all,all,all
  if speed=='slow'
    tslist=2,1.0
  else
    tslist=2,1.25
    Brake 13000
_,_,gpu = Processors
if gpu
  SimSteps Screen=25,Pairlist=25
else    
  SimSteps Screen=10,Pairlist=10
ts=tslist2*tslist1
savesteps=saveinterval/ts
TimeStep (tslist)
Temp (temperature)
fixedlist() = ListAtom fixed
if count fixedlist
  if ConsoleMode
    FreeAll
  else
    MarkAtom (fixedlist1)
    ShowMessage '(count fixedlist) atoms are currently fixed. This will yield unrealistic trajectories, normally distances should be restrained instead, see user manual at Essentials > The 10 magic words > Bond. Click Simulation > Free > All if you agree...' 
    Wait ContinueButton




i=00000
if format=='sim'
  trajectfilename='(MacroTarget)(i).sim'
else  
  restartfilename='(MacroTarget).sim'
  trajectfilename='(MacroTarget).(format)'
  old = FileSize (MacroTarget)(i).xtc
  if old
    RenameFile (MacroTarget)(i).xtc,(trajectfilename)
running = FileSize (trajectfilename)
if not running
  Experiment Minimization
  Experiment On
  Wait ExpEnd
  Sim On
else
  ShowMessage "Simulation has been running before, loading last snapshot..."
  if format=='sim'
    do
      i=i+1
      found = FileSize (MacroTarget)(i).sim
    while found
    i=i-1
    LoadSim (MacroTarget)(i)
    if i>0
      t = Time
      savesteps=0+t/(ts*i)
  else
    found = FileSize (restartfilename)
    if found
      last,t = Load(format) (trajectfilename),1
      if !last
        last,t = Load(format) (trajectfilename),2
        savesteps=0+t/ts
      LoadSim (restartfilename)
    else
      do
        i=i+1
        last,t = Load(format) (trajectfilename),(i)
        ShowMessage 'Searching (format) trajectory for last snapshot, showing snapshot (i) at (0+t) fs'
        Sim Pause
        Wait 1
      while !last
      savesteps=0+t/(ts*(i-1))
      Sim Continue
HideMessage
  
TempCtrl Rescale
PressureCtrl (pressurectrl)





Save(format) (trajectfilename),(savesteps)
if format!='sim'
  SaveSim (restartfilename),(savesteps),Number=no

if duration=='forever'
  Console On
  if ConsoleMode
    Wait forever
else
  measurements=0
  do
    if ConsoleMode
      Console On
    Wait 10
    Console Off
    measurements=measurements+1
    t = Time
  while t<1000.*duration+1
  vallist() = Tab Default
  if count vallist
    SaveTab default,(MacroTarget)_duringsim,Format=Text,Columns=(count vallist/measurements),Header='Insert your own header here'
  Sim Off
if runWithMacro and ConsoleMode and !IndentationLevel
  Exit
  