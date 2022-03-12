
RequireVersion 20.1.1

MacroTarget='2ONC-cac-r2md'



soluteobj=1

pdbsaved=0

bfactorscale=1.0

if not count block
  firstsnapshot=0
  snapshots='all'
  snapshotstep=1
  
refsnapshot=0

central=0

rmsdmin=2.0

figurewidth=512

hiresplotted=0


casel='CA Protein or C1* NucAcid'




resnummin=36
resnummax=766





rdfsellist()=''


def AnalyzeInsideCell
  global ligandsel,resnummin,resnummax
  
  celllist1,celllist2,celllist3=Cell
  Plot celllist,'Simulation cell lengths','Length in Angstrom','CellLengthX CellLengthY CellLengthZ'
  
  elist()=EnergyAll All
  Plot sum elist,'Total potential energy of the system','Energy in (EnergyUnit)','TotalEnergy'
  Plot elist,'Potential energy components','Energy in (EnergyUnit)','Bond Angle Dihedral Planarity Coulomb VdW Packing1D Packing3D'
  
  Plot "SurfObj Solute",'Surface areas of the solute','Surface areas in Angstrom^2','SurfVdW SurfMol SurfAcc'

  hbolist()=ListHBoAtom Obj Solute,Obj Solute
  Plot count hbolist/2,'Number of hydrogen bonds in the solute','Hydrogen bonds','SoluteHBonds'

  hbolist()=ListHBoAtom Obj Solute,Obj Solvent
  Plot count hbolist,'Number of hydrogen bonds between solute and solvent','Hydrogen bonds','SltSlvHBonds'
  
  global hbosel,hbondsmax
  if hbosel==''
    hbosel=ligandsel
  if hbosel!='' and hbondsmax!=0
    acceptors=CountAtom Element N O S P and (hbosel)
    donors=CountAtom Element H with bond to Element N O S P and (hbosel)
    bonds=acceptors*2+donors
    if hbondsmax=='all' or hbondsmax>bonds
      hbondsmax=bonds
    hbolist()=ListHBoAtom (hbosel),!Water,Results=6
    hbonds=count hbolist/6
    hboname=''
    for j=1 to hbondsmax
      hboname=hboname+'HB(j)Atm1 HB(j)Atm2 HB(j)E HB(j)D '
      if j<=hbonds
        for k=0 to 1
          hboresultlist(4*j-3+k)=ListAtom (hbolist(j*6-5+k)),Format='ATOMNAME.RESNAME1RESNUM.MOLNAME'
        hboresultlist(4*j-1)=hbolist(j*6-3)
        hboresultlist(4*j)=hbolist(j*6)
      else
        for k=1 to 4
          hboresultlist(4*j-4+k)='-'
    if hbondsmax>0
      WriteTable hboresultlist,'Hydrogen bonds made by (hbosel)','(hboname)'
    else
      WriteTable 'None','Hydrogen bonds made by (hbosel)','HB1Atm1'
  
  proteins=CountMol Protein
  if proteins
    Plot "SecStr",'Protein secondary structure content','Secondary structure fractions in percent','Helix Sheet Turn Coil Helix310 HelixPi'
  
  sel='Protein Obj Solute and Res (resnummin)-(resnummax)'
  secstrnamelist()=  'Helix','Sheet','Turn','Coil','Helix310','HelixPi'
  secstrletterlist()='H',    'E',    'T',   'C',   'G',       'I'
  for i=1 to count secstrletterlist
    secstrnum(secstrletterlist(i))=i
    secstrcolorlist(i)=ColorPar SecStr,(secstrnamelist(i))
  secstrreslist()=ListRes (sel) and SecStr (join secstrnamelist)
  if count secstrreslist
    secstrlist()=SecStrRes (join secstrreslist())
    for i=1 to count secstrlist
      secstrlist(i)=secstrnum(secstrlist(i))
    PlotRes 'secstr',secstrreslist,secstrlist,1,6,'Per-residue protein secondary structure','(join secstrnamelist)','(join secstrcolorlist)'
  
  conreslist()=ListRes Protein or NucAcid and Obj Solute Res (resnummin)-(resnummax)
  if count conreslist
    conlist()=CountConRes (join conreslist),Obj Solute,Cutoff=0.5,Subtract=HBoRadii,Exclude=5,Occluded=No
    PlotRes 'con',conreslist,conlist,1,15,'Per-residue number of contacts','','Blue Yellow'
  
  if ligandsel!=''
    sel='Protein or NucAcid and Obj Solute Res (resnummin)-(resnummax) and not (ligandsel)'
    ligconreslist()=ListRes (sel)
    if count ligconreslist
      for i=1 to count ligconreslist
        vallist(i)=0
      intlist='HBonds','Hydrophobic','Ionic'
      occlist='',      'No',         'Yes'
      for i=1 to count intlist
        if i==1
          intreslist() = ListHBoRes (sel),(ligandsel),Results=1
        else
          intreslist() = ListIntRes (sel),(ligandsel),Type=(intlist(i)),Exclude=5,Occluded=(occlist(i)),Results=1
        idx=1
        intresidues=count intreslist
        for j=1 to count ligconreslist
          if idx>intresidues
            break
          if ligconreslist(j)==intreslist(idx)
            vallist(j)=vallist(j)|(1<<(i-1))
            idx=idx+1
      vallist()=vallist+1
      PlotRes 'ligcon',ligconreslist,vallist,1,8,'Per-residue contacts with ligand','None HB Hyd Hyd+HB Ion Ion+HB Ion+Hyd Ion+Hyd+HB','White Red Green Yellow Blue Magenta Cyan Gray'
  

  
  
  

  


  
  
  
  
  
 
 
  
  
  
def AnalyzeOutsideCell
  global fof
  
  if fof=='YASARA2'
    for checktype in 'dihedrals','packing1d','packing3d'
      zscore(checktype)=CheckObj Solute,(checktype)
    zscore=zscoredihedrals*0.145+zscorepacking1d*0.390+zscorepacking3d*0.465 
    Plot zscore,'Structure quality of the solute','Quality Z-score','Quality'
  
  Plot "RadiusObj Solute,Center=Mass,Type=Gyration",'Radius of gyration of the solute','Radius in Angstrom','RadGyration'
  
  

  
  

def AnalyzeChange
  global casel,caselref,ligandsel,ligandselref
  
  resultlist1=SupAtom (casel),(caselref)
  resultlist2=SupAtom Backbone Obj Solute,Backbone Obj SoluteRef
  resultlist3=SupAtom Element !H Obj Solute,Element !H Obj SoluteRef
  Plot resultlist,'Solute RMSD from the starting structure','RMSD in Angstrom','RMSDCa RMSDBb RMSDAll'
  
  if ligandsel!=''
    if casel!='None'
      SupAtom (casel) and not (ligandsel),(caselref) and not (ligandselref)
    else
      SupAtom Element !H Obj Solute and not (ligandsel),Element !H Obj SoluteRef and not (ligandselref)
    result=RMSDAtom (ligandsel),(ligandselref)
    Plot result,'Ligand movement RMSD after superposing on the receptor','RMSD in Angstrom','RMSDLigMove'
  
  if ligandsel!=''
    result=SupAtom (ligandsel),(ligandselref)
    Plot result,'Ligand conformation RMSD after superposing on the ligand','RMSD in Angstrom','RMSDLigConf'
  





fof=ForceField


CellLengthXText='Conformational changes of the simulated solute molecules lead to fluctuations in density. '
                'If the simulation box has a constant size, changes in density lead to changes in pressure. '
                'This is not realistic, because molecules normally "live" in a constant pressure environment. '
                'During the simulation the cell is therefore rescaled to maintain a constant cell '
                'pressure. Depending on the chosen pressure control mode, '
                'the three cell axes are either rescaled together [Manometer1D], partly together '
                '[X- and Z-axes, Manometer2D, used for membrane simulations], independently [Manometer3D], '
                'or not at all [Off]. You can deduce the pressure control mode from the plot below.'

TotalEnergyText='The total potential energy of the system is plotted, '
                'according to the (fof) force field. If you ran the simulation '
                'with a different force field, you need to adapt the `ForceField` command '
                'at the top of this macro accordingly.\n\n'
                'When the simulation is started from an energy-minimized "frozen" conformation, '
                'there is usually a sharp increase in energy during the first picoseconds, '
                'since the added kinetic energy is partly stored as potential energy. '
                'Also on a larger time-scale, the potential energy will often not decrease. '
                'A common reason are counter ions. These are initially placed at the '
                'positions with the lowest potential energy, usually close to charged solute groups, '
                'from where they detach to gain entropy, but also potential energy. '

def PrintBondInfo
  global fof
  text='The following individual components of the total potential energy are plotted: '
       'bond energies [Bond], bond angle energies [Angle], dihedral angle energies [Dihedral], '
       'planarity or improper dihedral energies [Planarity], Van der Waals energies [VdW]'
  if fof=='YASARA2'
    text=text+', electrostatic energies [Coulomb], 1D packing energies [Packing1D] and 3D packing energies [Packing3D]. '
  else
    text=text+' and electrostatic energies [Coulomb]. '
  text=text+'Force field energies help to judge the structural quality of a protein: '
            'distortions of local covalent geometry can be found by looking at the bond, angle and planarity energies. '
            'Unrealistically close contacts [bumps] lead to a high Van der Waals energy, '
            'just like a large number of hydrogen bonds [since they pull the atoms closer than '
            'their normal Van der Waals contact distance]. The Coulomb energy is the least '
            'informative, because it strongly depends on the amino acid composition '
            '[e.g. proteins with a net charge have a higher Coulomb energy].'
  WriteReport Paragraph,'(text)'

HelixText='The total percentages of alpha helices, beta sheets, turns, coils, 3-10 helices and pi helices are '
          'calculated and plotted. For clarification, a turn is simply a stretch of four residues that are not '
          'part of other secondary structure elements and form a hydrogen bond between the O of the first and '
          'the NH of the last residue. A coil is anything that does not fit into the other categories. '
          'Note that pi-helices [helices with hydrogen bonds between residues N and N+5] are rather unstable and '
          'thus do not normally occur in proteins, except for short bulges in alpha helices '
          '[which are often the result of single residue insertions and prolines]. '

SecStrText='The following plots show the protein secondary structure per residue as a function of '
           'simulation time. They are helpful to monitor protein folding and all other kinds of structural '
           'changes. The default secondary structure colors are used, you can change them at '
           'View > Color > Parameters > Secondary structure colors. One plot per protein molecule is shown.'

ConText='The number of contacts per residue as a function of simulation time is shown in the following plots. '
        'There is one plot for each protein or nucleic acid molecule. Even though contacts between atoms separated '
        'by up to four chemical bonds are excluded, neighboring residues in the molecule usually have enough close '
        'atoms to be counted as a contact. Consequently residues with zero contacts are very rare and often glycines. '
        'The number of contacts tells you how densely a certain residue range is packed and allows to identify structurally '
        'very important residues, e.g. a phenylalanine in the hydrophobic core can contact 15 or more other residues. '

LigConText='The following plots show the types of contacts made with the ligand `(ligandsel)` as a function of simulation time. '
            'There is one plot for each protein or nucleic acid molecule. Three types of contacts are shown: Hydrogen bonds [red] ',
            'hydrophobic contacts [green] and ionic interactions [blue]. Also mixtures of these three colors can show up if a '
            'certain residue is involved in more than one type of contact with the ligand [see plot legend]. '

SurfVdWText='The Van der Waals [SurfVdW], molecular [SurfMol] and solvent accessible [SurfAcc] surface areas of the '
            'solute in A^2 are plotted. The difference between these surface types can be summarized as follows:\n\n'
            '`Van der Waals surface`: if you think of atoms as spheres with a given Van der Waals radius, '
            'then the Van der Waals surface consists of all the points on these spheres that are not inside another sphere. '
            'In practice, the Van der Waals surface is of limited use, because it can be found throughout a protein and '
            'does not tell much about the interaction with the solvent.\n\n'
            '`Molecular surface`: this is the Van der Waals surface from the viewpoint of a solvent molecule, '
            'which is a much more useful concept. The water is assumed to be a sphere of a given radius '
            '[also called the water probe], that rolls over the solute. '
            'Those parts of the Van der Waals surface that the water probe can touch are simply copied to the molecular surface '
            '[and called the contact surface]. Clefts in the Van der Waals surface that are too narrow for the water probe to enter '
            'are replaced by the Van der Waals surface of the water probe itself [and called the reentrant surface]. '
            'So the molecular surface is a smooth composition of two Van der Waals surfaces: '
            'the one of the solute and the one of the solvent molecule while it traces the contours of the solute. '
            'Other common names for the molecular surface are solvent excluded surface or Connolly surface.\n\n'
            '`Solvent accessible surface`: this surface consists of all the points that the center of the water probe '
            '[i.e. the nucleus of the oxygen atom in the water molecule] can reach while rolling over the solute. '
            'The shortest possible distance between the water oxygen nucleus and a solute atom is simply '
            'the sum of the Van der Waals radii of the solute atom and the water probe.'

def PrintSoluteHBondsInfo
  WriteReport Paragraph,
    'The number of hydrogen bonds inside the solute is plotted below. '
    'One hydrogen bond per hydrogen atom is assigned at most, '
    'picking the better one if two acceptors are available.'
    'The following formula yields the bond energy in [kJ/mol] '
    'as a function of the Hydrogen-Acceptor distance in [A] and two scaling factors:'
  WriteReport Image,Filename=(YASARADir)/doc/ListHBoAtomResMolObj1.png,Style=Center,Name="formula_energyhbo0"
  WriteReport Paragraph,
    'The first scaling factor depends on the angle formed by Donor-Hydrogen-Acceptor:'
  WriteReport Image,Filename=(YASARADir)/doc/ListHBoAtomResMolObj2.png,Style=Center,Name="formula_energyhbo1"
  WriteReport Paragraph,
    'The second scaling factor is derived from the angle formed by Hydrogen-Acceptor-X, '
    'where X is the atom covalently bound to the acceptor. If X is a heavy atom:'
  WriteReport Image,Filename=(YASARADir)/doc/ListHBoAtomResMolObj3.png,Style=Center,Name="formula_energyhbo2"
  WriteReport Paragraph,
    'If X is a hydrogen, slightly smaller angles are allowed:'
  WriteReport Image,Filename=(YASARADir)/doc/ListHBoAtomResMolObj4.png,Style=Center,Name="formula_energyhbo3"
  WriteReport Paragraph,
    'A hydrogen bond is counted if the hydrogen bond energy obtained with this formula '
    'is better than 6.25 kJ/mol [or 1.5 kcal/mol], which is 25% of the optimum value 25 kJ/mol. '

SltSlvHBondsText='The plot shows the number of hydrogen bonds between solute and solvent. '
                 'Together with the plot above, it is a good indicator for successful protein folding, '
                 'indicated by a decreasing number of bonds with the solvent and a growing '
                 'number of bonds within the solute.'

def PrintHB1Atm1Info
  global hbosel,hbondsmax,EnergyUnit
  caller3 acceptors,donors,bonds
  WriteReport Paragraph,
    'The following table shows all hydrogen bonds made by (hbosel). '
    'With (acceptors) acceptors and (donors) donors, '
    'a total number of (bonds) hydrogen bonds are possible. '
    '(hbondsmax) hydrogen bonds are listed - labeled HB1 to HB(hbondsmax). '
    'The first atom of the bonding pair is labeled Atm1 and the second Atm2, respectively. '
    'The atom ID separates atom name, residue ID and molecule name with dots. A lower-case "h" '
    'indicates hetgroups. E and D are short for the hydrogen bonding energy in [(EnergyUnit)] '
    'and the distance between the bonding partners in [A]. '
    'To list other hydrogen bonds, edit the `hbosel` variable at the beginning of this macro. '
    'To list more or fewer hydrogen bonds, edit the `hbondsmax` variable at the beginning of this macro.'


def PrintQualityInfo
  WriteReport Paragraph,
    'When validating a structure, one can use a gold standard of highest resolution reference '
    'structures to obtain estimates for the expected average energy and its standard deviation, '
    'and then calculate how many standard deviations the actual energy is away from the average, '
    'thereby obtaining a Z-score:'
  WriteReport Image,Filename=(YASARADir)/doc/CheckAtomResObjAll1.png,Style=Center,Name="formula_zscore"
  WriteReport Paragraph,
    'In the formula above, `x` is the energy of the current structure, '
    '`my` and `sigma` are the average value and the standard deviation '
    'of the energies in the gold standard population. '
    'Assuming that the gold standard population has a standard distribution, '
    '95.4% of all proteins have Z-scores between -2 to +2. '
    'The others can be called outliers. '
    'Positive outliers are usually small perfect proteins like a single alpha helix, '
    'and negative outliers are proteins with serious errors.'
  MakeTab zscores,2,2
  descriptionlist()='disgusting','terrible','bad','poor','satisfactory','good'
  for i=-5 to 0
    Tabulate '< (i)','(descriptionlist(i+6))'
  Tabulate '> 0','optimal'
  WriteReport Table,zscores,"Mapping between Z-score and human language.",.0f,InfoColumnName='Z-score',DataColumnName='Description'

def PrintRadGyrationInfo()
  WriteReport Paragraph,
    'After determining the center of mass of the solute, the radius '
    'of gyration is calculated and plotted according to this formula:'
  WriteReport Image,Filename=(YASARADir)/doc/RadiusAtomResMolObjAll2.png,Style=Center,Name="formula_gyrrad"
  WriteReport Paragraph,
    'In this formula, `C` is the center of mass, and `Ri` is the position of atom `i` of `N`.'


def PrintRMSDCaInfo
  global casel,calphas
  WriteReport Paragraph,
    'The plot shows Calpha [RMSDCa], backbone [RMSDBb] and all-heavy atom [RMSDAll] RMSDs calculated '
    'according to this formula, where `Ri` is the vector linking the positions of atom `i` [of `N` atoms] '
    'in the reference snapshot and the current snapshot after optimal superposition: '
  WriteReport Image,Filename=(YASARADir)/doc/RMSDAtomResMolObj1.png,Style=Center,Name="formula_rmsd"
  if casel=='None'
    text='Less than three atoms matched the Calpha selection `(casel)`, therefore the Calpha RMSD plot '
         'graph is set to flat zero because at least three atoms are needed for structure superposition. '
  else
    text='The selection for the Calpha RMSD calculation is `(casel)`, this matched (calphas) atoms. '
    if casel=='CA Protein or C1* NucAcid and Obj Solute'
      text=text+'The Calpha selection thus includes the main backbone carbon C1* of nucleic acids, so '
                'the plot also shows a Calpha RMSD if you simulate just nucleic acids. In simulations '
                'of protein-DNA complexes, the Calpha RMSD therefore considers the DNA too. '
  text=text+'To change the Calpha selection, edit the `casel` variable at the beginning of this macro.'
  WriteReport Paragraph,'(text)'

RMSDLigMoveText='The following plot shows the RMSD of the ligand heavy atoms '
                'over time, measured after superposing the receptor on its reference structure. '
                'This procedure delivers information about the movement of the ligand in its binding pocket.'

RMSDLigConfText='This plot displays the RMSD of the ligand atoms '
                'over time, measured after superposing on the reference structure of the ligand. '
                'The gained data summarize the conformational changes of the ligand. '

UnselectAll

if MacroTarget==''
  RaiseError "This macro requires a target. Either edit the macro file or click Options > Macro > Set target to choose a target structure"

if resnummax<resnummin
   RaiseError 'The maximum residue number (resnummax) is smaller than the minimum (resnummin), please swap resnummin and resnummax'

Clear
Console Off
SurfPar Molecular=Numeric

if not count block
  id=''

waterscene=FileSize (MacroTarget)_water.sce
solventscene=FileSize (MacroTarget)_solvent.sce
if waterscene
  LoadSce (MacroTarget)_water
elif solventscene 
  LoadSce (MacroTarget)_solvent
else
  RaiseError 'Could not find initial scene file (MacroTarget)_water.sce. You must run a simulation with the macro md_run first'

objs1=CountObj Solute and not (soluteobj)
objs2=CountObj SoluteRef
if sum objs
  RaiseError "The object names 'Solute' and 'SoluteRef' are reserved for use by this macro, please rename your objects and try again"
solutename = NameObj (soluteobj)
NameObj (soluteobj),Solute
NameObj Water,Solvent

NameMol ' ','_' 

if casel==''
  calphas=0
else
  calphas=CountAtom (casel)
if calphas<3
  casel='None'

ShowMessage "Preparing analysis, please wait..."
Wait 1

old=FileSize (MacroTarget)00000.xtc
if old
  RenameFile (MacroTarget)00000.xtc,(MacroTarget).xtc

for format in 'xtc','mdcrd','sim'
  found=FileSize (MacroTarget).(format)
  if found
    break
  
LigandObjWarning=''
if ligandsel==''
  mols=CountMol Obj Solute
  if mols>1
    ligandreportstring='automatically by YASARA'
    reslist()=ListRes Hetgroup !Water Obj Solute with >0 bonds to all
    reslenmax=6
    ligandname=''
    carbohydlist()=''
    for res in reslist
      reslen=CountAtom (res)
      resname=NameRes (res)
      carbons=CountAtom Element C (res) with bond to Element O
      if carbons>4 and resname not in carbohydlist
        carbohydlist(1+count carbohydlist)=resname
      if reslen>reslenmax
        reslenmax=reslen
        ligandname=resname
    if ligandname!=''
      if ligandname in carbohydlist
        ligandsel='Res'
        for i=2 to count carbohydlist
          ligandsel=ligandsel+' "(carbohydlist(i))"'
      else
        ligandsel='Res "(ligandname)"'
else
  ligandobjlist()=ListObj (ligandsel) and not Membrane Solvent
  if !count ligandobjlist
    RaiseError 'Your ligand selection "(ligandsel)" did not match any atoms. If you want to select a residue add "Res" at the beginning of your selection.'
  if count ligandobjlist>1
    RaiseError 'Your ligand "(ligandsel)" is present in objects (ligandobjlist), but only one ligand object is allowed'
  solfirstatm,sollastatm=SpanAtom Obj Solute
  ligfirstatm,liglastatm=SpanAtom Obj (ligandobjlist1)
  if ligfirstatm==sollastatm+1
    ligatoms=CountAtom (ligandsel)
    ligandobjname=NameObj (ligandobjlist1)
    SegAtom (ligandsel),LigO
    ligandsel='Segment LigO'
    JoinObj (ligandobjlist),Solute
    if ligatoms!=liglastatm-ligfirstatm+1
      LigandObjWarning='`WARNING:` The selected ligand is only part of a larger object. The remaining atoms in '
                        'object (ligandobjlist1) with name `(ligandobjname)` were treated as solute atoms.'
  elif ligfirstatm!=solfirstatm 
    RaiseError 'Ligand (ligandsel) is not part of the solute, no ligand RMSD values can be calculated. '
               'You could click Edit > Join > Objects to join the ligand object (ligandobjlist1) with the solute object (soluteobj) '
               'in the original scene, but this would change the order of atoms and require to run the MD once again.'
  ligandreportstring='by the user'

if refsnapshot=='average'
  filename='(resultbase)_average(id).pdb'
  exists=FileSize (filename)
  if !exists
    RaiseError "No time-average structure has been calculated yet, cannot superpose onto it. Run the macro once with refsnapshot=0 (or any other number), then run again with refsnapshot='average'"
  refobj=LoadPDB (filename)
  SupObj (refobj),Solute
else
  if refsnapshot
    if format=='sim'
      LoadSim (MacroTarget)(00000+refsnapshot)
    else
      Load(format) (MacroTarget),(refsnapshot+1)
  refobj=DuplicateObj Solute
NameObj (refobj),SoluteRef
ligandselref=''
caselref='None'
if ligandsel!=''
  ligandselref=ligandsel+' Obj SoluteRef'
  ligandsel=ligandsel+' Obj Solute'
if casel!='None'
  caselref=casel+' and Obj SoluteRef'
  casel=casel+' and Obj Solute'
RemoveObj SoluteRef

if rmsdmin
  clusterobjlist()=refobj

if dccmsel!=''
  dccmunits=Count(dccmsel) and Obj Solute
  if !dccmunits and dccmsel!='Atom CA Protein or C1* NucAcid'
    RaiseError 'The DCCM selection (dccmsel) did not match any atoms'

task='AddHeader'
tables=0
plotrestabs=0

Sim On
AnalyzeInsideCell
Sim Off
AnalyzeOutsideCell
AddObj SoluteRef
AnalyzeChange
RemoveObj SoluteRef
if !count introwritten
  ShowSystem
  SaveScreenshot 'sim','simulated system'
  ShowSolute
  Style Ribbon,Stick
  StickRadius 20
  SaveScreenshot 'solute','solute object'
  if ligandsel!=''
    ligandobj=DuplicateRes (ligandsel)
    SwitchAll off
    SwitchObj (ligandobj),on
    NiceOriObj (ligandobj)
    LabelAtom Obj (ligandobj) Element !H,Format=ATOMNAME,Height=.3
    LabelAtom Obj (ligandobj) Element H,Format=ATOMNAME,Height=.2
    ZoomObj (ligandobj),Steps=0
    Style BallStick
    SaveScreenshot 'ligand','ligand'
    DelObj (ligandobj)
  ShowSystem

for i=1 to tables
  MakeTab (tabnamelist(i)),Dimensions=2,Columns=(tabcolslist(i))

ycolumn=3

if tabrowsmax=='all'
  tabrowsmax=0

task='Tabulate'
i=00000+firstsnapshot
emin=1e99
last=0
while !last and snapshots!=0
  lastsnapshot=i
  if format=='sim'
    sim=FileSize (MacroTarget)(i+1).sim
    if not sim
      last=1
    LoadSim (MacroTarget)(i)
  else
    last=Load(format) (MacroTarget),(i+1)
  Sim Pause
  if central
    _,_,_,_,_,_,cen()=Cell
    pos()=PosAtom (central)
    MoveAtom all,(cen-pos)
  simtime=Time
  ShowMessage 'Analyzing snapshot (0+i) at (0+(simtime/1000)) ps'
  Wait 1
  for j=1 to tables
    SelectTab (tabnamelist(j))
    Tabulate (simtime/1000),(simtime/1000000)
  AnalyzeInsideCell
  if count rdfsellist==4
    BinDistance (rdfsellist)
  e=EnergyObj Solute
  Sim Off
  if e<emin
    emin=e
    SaveSce (resultbase)_energymin(id)
    SavePDB Solute,(resultbase)_energymin(id)
  if last
    SaveSce (resultbase)_last
    SavePDB Solute,(resultbase)_last
  if pdbsaved
    SavePDB Solute,(resultbase)_(i)
  AnalyzeOutsideCell
  if espmethod!=''
    RemoveObj Solvent
    SaveESP (resultbase)(i).phi.gz,(espmethod)
    AddObj Solvent
  AddObj SoluteRef
  AnalyzeChange
  SupAtom Element !H Obj Solute,Element !H Obj SoluteRef
  AddPosAtom Obj Solute
  if rmsdmin
    AddObj (join clusterobjlist) and not SoluteRef
    for j=1 to count clusterobjlist
      aarmsd = SupAtom Element !H Obj Solute,Element !H Obj (clusterobjlist(j))
      if aarmsd<=rmsdmin
        break
    if aarmsd>rmsdmin
      obj = DuplicateObj Solute
      members=count clusterobjlist+1
      clusterobjlist(members)=obj
      NameObj (obj),Cluster(members)
      PropObj (obj),(i)
    RemoveObj (join clusterobjlist)
  else
    RemoveObj SoluteRef
  if !last
    for j=2 to snapshotstep
      if format=='sim'
        sim=FileSize (MacroTarget)(i+j).sim
        if not sim
          LoadSim (MacroTarget)(i+j-1)
          last=1
      else
        last=Load(format) (MacroTarget),(i+j)
      if last
        Sim Off
        SaveSce (resultbase)_last
        SavePDB Solute,(resultbase)_last
        break
  i=i+snapshotstep
  if snapshots!='all'
    snapshots=snapshots-1
if i==firstsnapshot
  RaiseError "This macro is meant to analyze a molecular dynamics trajectory created with md_run, but none was found in this directory"
snapshots=(-firstsnapshot+i)/snapshotstep

if rmsdmin
  if count block and block>1
    DelObj (clusterobjlist1)
    DelVar clusterobjlist1
  clustermembers=count clusterobjlist
  AddObj (join clusterobjlist)
  for j=1 to clustermembers
    obj=clusterobjlist(j)
    snum=PropObj (obj)
    clusterfilenamelist(j)='(resultbase)_cluster(id)_(zeroed clustermembers+j)_snapshot(0+snum)'
    SaveYOb (obj),(clusterfilenamelist(j))
  DelObj (join clusterobjlist) and not SoluteRef
  RemoveObj SoluteRef

simtime=Time
if simtime<1000000
  xcolumn=1
  plottimestring='picoseconds'
  tabtimestring='Time [ps]'
  simtime=0.00+simtime/1000
else
  xcolumn=2
  plottimestring='nanoseconds'
  tabtimestring='Time [ns]'
  simtime=0.00+simtime/1000000

AveragePosAtom Obj Solute
_ = RMSFAtom Obj Solute,Unit=BFactor
if bfactorscale!=1.0
  firstatm,lastatm=SpanAtom Obj Solute
  for i=firstatm to lastatm
    bf=BFactorAtom (i)
    BFactorAtom (i),(bf*bfactorscale)    
if refsnapshot!='average'
  SavePDB Solute,(resultbase)_average(id)
MakeTab RMSF
firstatm,lastatm=SpanAtom Obj Solute
rmsflist()=RMSFAtom Obj Solute
for i=firstatm to lastatm
  res=ListAtom (i),Format="'ATOMNAME','RESNAME','RESNUM','MOLNAME'"
  Tabulate '(i)',(res),(rmsflist(i-firstatm+1))
SaveTab RMSF,(resultbase)_rmsf(id),Format=Text,Columns=6,NumFormat=8.2f,"Table of atomic Root Mean Square Fluctuations in [A]"
MakeTab RESRMSF
reslist()=ListRes Obj Solute,Format="'RESNAME','RESNUM','MOLNAME'"
resrmsflist()=RMSFRes Obj Solute
for i=1 to count reslist
  Tabulate (reslist(i)),(resrmsflist(i))
SaveTab RESRMSF,(resultbase)_rmsfres(id),Format=Text,Columns=4,NumFormat=8.2f,"Table of residue Root Mean Square Fluctuations in [A]" 

bound=Boundary
ycolumn=3
if !count introwritten
  introwritten=1
  MakeTab Info,2,2
  infolist()='Protein molecules','Mol Protein','Protein residues','Res Protein','Protein atoms','Atom Protein',
             'Nucleic acid molecules','Mol NucAcid','Nucleic acid residues','Res NucAcid','Nucleic acid atoms','Atom NucAcid'
  for i=1 to count infolist step 2
    Tabulate '(infolist(i))'
    Tabulate Count(infolist(i+1))
  resnamelist()=ListRes !Protein and !NucAcid and !Water,Format='RESNAME'
  if count resnamelist
    unilist()=resnamelist1
    for resname in resnamelist
      if resname not in unilist
        unilist(count unilist+1)=resname
    for resname in unilist
      residues=CountRes '(resname)'
      resid=ListRes '(resname)'
      resatoms=CountAtom Res (resid)
      if resatoms==1
        infoatomlist(count infoatomlist+1)=ListAtom (resid),Format='Element ATOMElement'
        infoatomslist(count infoatomslist+1)=residues
      else 
        Tabulate 'Residue (resname) with (resatoms) atoms',(residues)
  for i=1 to count infoatomlist
    Tabulate '(infoatomlist(i))',(infoatomslist(i))
  Tabulate 'Water residues'
  Tabulate CountRes Water
  Tabulate 'Total number of atoms'
  Tabulate (Atoms)
  
  WriteReport Title,Filename=(resultbase)_report,Text='YASARA Molecular Dynamics Trajectory Analysis for (basename MacroTarget)'
  WriteReport Heading,2,'About the simulation'
  if not count block
    WriteReport Paragraph,
      'The trajectory `(MacroTarget)` has been analyzed with YASARA version (Version) over a period of '
      '(simtime) (plottimestring) with (snapshots) snapshots and the (fof) force field. Note that the MD '
      'simulation may have been run with a different force field, but (fof) was used to calculate the '
      'energies in this report. To change this, edit the ForceField setting at the start of this macro. '
  else
    WriteReport Paragraph,
      'The trajectory has been analyzed with YASARA version (Version) in blocks of (blocksnapshots) snapshots and the (fof) force field.'
  WriteReport Paragraph,
    'All plots and pictures in this report [like the simulated system below] are (figurewidth) pixels wide, you '
    'can change the `figurewidth` variable in this macro as needed.'
  caption='A ray-traced picture of the simulated system. The simulation cell boundary is set to (bound). '
          'Atoms that stick out of the simulation cell will '
  if bound=='periodic'
    caption=caption+'be wrapped to the opposite side of the cell during the simulation.'
  else
    caption=caption+'not be included in the simulation.'
  WriteReport Image,Filename=(resultbase)_sim.png,Style=Figure,Caption=(caption),Delete=Yes
  WriteReport Heading,3,'Composition of the system'
  WriteReport Paragraph,'The components of the system are shown in the table below. '
  WriteReport Table,Info,'Composition of the simulated system',.0f,InfoColumnName="Type",DataColumnName="Number"
  WriteReport Paragraph,
    'Object (soluteobj) with name `(solutename)` has been identified as the solute and is shown below. '
    'If this is not the intended solute, please change the `soluteobj` variable in this macro. (LigandObjWarning)'
  WriteReport Image,Filename=(resultbase)_solute.png,Style=Figure,Caption='The solute oriented along the major axes.',Delete=Yes
  if ligandsel!=''
    ligandresidues=CountRes (ligandsel)
    ligandatoms=CountAtom (ligandsel)
    WriteReport Heading,3,'The ligand'
    WriteReport Paragraph,
      'A special analysis has been performed for the ligand, chosen (ligandreportstring) with the selection `(ligandsel)`. '
      'The number of residues matching the ligand selection is (ligandresidues), with (ligandatoms) atoms. '
      'To change the ligand selection, edit the `ligandsel` variable at the beginning of this macro.'
    WriteReport Image,Filename=(resultbase)_ligand.png,Style=Figure,
      Caption='A ray-traced picture of the ligand (ligandsel). Bonds are colored by '
              'their order: Gray = 1, blue = 1.25, magenta = 1.33, red = 1.5, orange = 1.66, '
              'bright orange = 1.75, yellow = 2, lime green = 2.5, green = 3 and cyan = 4.',Delete=Yes
WriteReport Heading,2,'Analyses inside the simulation cell'
WriteReport Paragraph,
  'This section shows all analyses that have been performed inside the simulation cell, '
  'when all atoms share the common coordinate system of the simulation cell. '
if bound=='periodic'  
  WriteReport Paragraph,
    'Periodic boundaries are active and considered for distance measurements. Calculations '
    'that involve groups of atoms [center of mass, regression lines, enclosing spheres..] '
    'are ambiguous and should be placed in the next section, unless it is known that the '
    'atom group does not drift through a periodic boundary.'
task='WriteReport'
Sim on
AnalyzeInsideCell
Sim off
WriteReport Heading,2,'Analyses outside the simulation cell'
WriteReport Paragraph,
  'The following section presents data gathered outside the simulation cell, where each object '
  'has its own local coordinate system and no periodic boundaries are present. Calculations '
  'that involve the interaction between objects [common surface areas, contacts between objects..] '
  'must be placed in the previous section.'
AnalyzeOutsideCell
AddObj SoluteRef
if refsnapshot=='average'
  text='Analyses performed with respect to the time average structure'
else
  if !refsnapshot
    text='Analyses performed with respect to the starting structure'
  else
    text='Analyses performed with respect to snapshot (refsnapshot)'
WriteReport Heading,2,'(text)'
WriteReport Paragraph,
  '(text) are shown in this section. '
  'These are also done outside the simulation cell, where each object has its own local '
  'coordinate systems and no periodic boundaries are present. To choose another reference '
  'snapshot than (refsnapshot), edit the `refsnapshot` variable at the beginning of this macro.'
AnalyzeChange
RemoveObj SoluteRef
for func in 'Mean','Min','Max'
  SelectTab Main
  Tabulate "_"
  Tabulate '(func)'
  for i=3 to tabcolslist1
    vallist()=Tab Main,Column=(i)
    Tabulate ((func) vallist)
tabname='(resultbase)_analysis'
for i=1 to tables
  if i>1
    tabname=tabname+'_(tabnamelist(i))'
  SaveTab (tabnamelist(i)),(tabname)(id),Format=Text,Columns=(tabcolslist(i)),NumFormat=12.3f,(tabheaderlist(i))

rmsfreslist()=RMSFRes Obj Solute
residlist()=ListRes Obj Solute
resnumlist()=ListRes Obj Solute,'RESNUMWIC'
sel='Protein or NucAcid and Obj Solute'
molnamelist()=ListMol (sel),Format='Mol MOLNAME'
molnametitlelist()=ListMol (sel),Format='MOLNAME'
mols=count molnamelist
mmresidlist()=ListRes (sel)
mmresidues=count mmresidlist
if mmresidues
  mmrmsflist()=ShortList rmsfreslist,residlist,mmresidlist
  mmmolnamelist()=ListRes (sel),'Mol MOLNAME'
  mmresnumlist()=ListRes (sel),'RESNUMWIC'
  mmresnumlist()=0+mmresnumlist
  minlist()=min mmresnumlist,0+resnummin
  maxlist()=max mmresnumlist,0+resnummax
  rmsfresnummin=max minlist
  rmsfresnummax=min maxlist
  MakeTab MMRMSF,2,(mols+1),(rmsfresnummax-rmsfresnummin+1)
  Tab MMRMSF,Set=0
  for i=rmsfresnummin to rmsfresnummax
    Tab MMRMSF,1,(i-rmsfresnummin+1),Set=(i)
  for i=1 to mmresidues
    if mmresnumlist(i)>=rmsfresnummin and mmresnumlist(i)<=rmsfresnummax
      for j=1 to mols
        if mmmolnamelist(i)==molnamelist(j)
          Tab MMRMSF,(1+j),(mmresnumlist(i)-rmsfresnummin+1),Set=(mmrmsflist(i))

sel='!Protein !NucAcid !Water Obj Solute'
hetresinfolist()=ListRes (sel),"'MOLNAME','RESName``RESNUM','ATOMNUM'"
hetresidues=count hetresinfolist
if hetresidues
  hetresidlist()=ListRes (sel)
  hetrmsflist()=ShortList rmsflist,residlist,hetresidlist
  MakeTab HetRMSF,2,4
  for i=1 to hetresidues
    Tabulate (hetresinfolist(i)),(rmsfreslist(i))

if mmresidues or hetresidues
  WriteReport Heading,2,'Solute residue RMSF'
  WriteReport Paragraph,
    'The Root Mean Square Fluctuation [RMSF] per solute residue is calculated from the average RMSF of its constituting atoms. '
    'The RMSF of atom i with j runing from 1 to 3 for the x, y, and z coordinate of the atom position vector P and '
    'k runing over the set of N evaluated snapshots is given by following formula:'
  WriteReport Image,Filename=(YASARADir)/doc/RMSFAtomResMol1.png,Style=Center,Name="formula_rmsf"
  if mmresidues
    WriteReport Paragraph,'Each graph in the following plot represents one molecule, '
                          'so that you can easily see differences between molecules. '
                          'Note: Residue numbers are not unique, so graphs can overlap. '
    WriteReport Plot,'The Root Mean Square Fluctuation [vertical axis] per solute protein/nucleic acid residue [horizontal axis] '
                      'calculated from the average RMSF of the atoms constituting the residue. '
                      'A RMSF of exactly zero means that that residue number is not present in the molecule. '
                      'Atom RMSF table: "file://(basename resultbase)_rmsf(id).tab", '
                      'residue RMSF table: "file://(basename resultbase)_rmsfres(id).tab"',
                MMRMSF,Width=(figurewidth),Height=480,Title='Solute protein/nucleic acid residue RMSF',
                XColumn=1,YColumn=2,YColumns=(mols),XLabel='Residue number',
                YLabel="RMSF in Angstrom",LegendPos='Outside',Graphname=(molnamelist)
    if hiresplotted
      SavePlot Filename="LastReportPlot_hires",MMRMSF,Width=1600,Height=1200,Title='Solute protein/nucleic acid residue RMSF',
               XColumn=1,YColumn=2,YColumns=(mols),XLabel='Residue number',YLabel="RMSF in Angstrom",Graphname=(molnamelist)
    if mols>1
      WriteReport Paragraph,'In case the plot above is too crowded, '
                            'the per-residue RMSF values are shown separately for all (mols) molecules in the following plots:'
      for i=1 to mols
        WriteReport Plot,'The Root Mean Square Fluctuation [vertical axis] per solute protein/nucleic acid residue [horizontal axis] '
                          'calculated from the average RMSF of the atoms constituting the residue. '
                          'A RMSF of exactly zero means that that residue number is not present in the molecule. '
                          'Atom RMSF table: "file://(basename resultbase)_rmsf(id).tab", '
                          'residue RMSF table: "file://(basename resultbase)_rmsfres(id).tab"',
                    MMRMSF,Width=(figurewidth),Height=480,
                    Title='Solute protein/nucleic acid residue RMSF of molecule (molnametitlelist(i))',
                    XColumn=1,YColumn=(1+i),YColumns=1,XLabel='Residue number',
                    YLabel="RMSF in Angstrom",LegendPos='Outside',Graphname=(molnamelist(i))
        if hiresplotted
          SavePlot Filename="LastReportPlot_hires",MMRMSF,Width=1600,Height=1200,
                   Title='Solute protein/nucleic acid residue RMSF of molecule (molnametitlelist(i))',
                   XColumn=1,YColumn=(1+i),YColumns=1,XLabel='Residue number',YLabel="RMSF in Angstrom",Graphname=(molnamelist(i))
  if hetresidues
    WriteReport Table,HetRMSF,'RMSF in Angstrom for non-protein/nucleic acid residues in the solute.',InfoColumnName='Mol',DataColumnName='Residue','First atom','RMSF[A]'

if count rdfsellist==4
  rdflist()=RDF
  MakeTab RDF,2,2
  for i=1 to count rdflist
    Tabulate (rdfsellist4*i),(rdflist(i))
  SaveTab RDF,(resultbase)_rdf(id),Format=Text,Columns=2,NumFormat=6.3f,
          'Radial distribution function with parameters (rdfsellist) as a function of the radial distance in A'
  WriteReport Heading,2,'Radial Distribution Function'
  WriteReport Paragraph,
    'The Radial Distribution Function is calculated by first by determining distances '
    'between all (rdfsellist1) - (rdfsellist2) pairs and sorting them in (rdfsellist3) bins '
    'with a bin width of (rdfsellist4) A. The RDF is then computed with the following formula:'
  WriteReport Image,Filename=(YASARADir)/doc/RDF1.png,Style=Center,Name="formula_rdf"
  WriteReport Paragraph,
    'The RDF in bin `i` is thus calculated from the number of `CountsInBin i` divided by `Atoms1` '
    '[the number of atoms matching the first selection (rdfsellist1)] times the volume of the '
    'shell corresponding to bin i. To change the selection for the RDF edit the variable "rdfsel" '
    'at the beginning of this macro.'
  WriteReport Plot,'Radial Distribution Function [vertical axis] for (rdfsellist1) - (rdfsellist2) pairs '
                   'as a function of the radial distance r in units of A [horizontal axis] '
                   'calculated from (rdfsellist3) bins with (rdfsellist4) A bin width.'
                   'A table with the raw data is available here: "file://(basename resultbase)_rdf(id).tab"',
              RDF,Width=(figurewidth),Height=480,Title='RDF for (rdfsellist1) - (rdfsellist2) pairs',
              XColumn=1,YColumn=2,YColumns=1,XLabel='Radial distance r in A',YLabel="RDF(r)",LegendPos='Outside',Graphname='RDF'
  if hiresplotted
    SavePlot Filename="LastReportPlot_hires",RDF,Width=1600,Height=1200,Title='RDF for (rdfsellist1) - (rdfsellist2) pairs',
             XColumn=1,YColumn=2,YColumns=1,XLabel='Radial distance r in A',YLabel="RDF(r)"

if espmethod!=''
  SwitchObj not SimCell Solute,off
  esp=LoadESP (MacroTarget)(lastsnapshot)
  PointPar Radius=32.00
  SaveScreenshot 'esp','electrostatic potential'
  DelObj (esp)
  Switch On
  WriteReport Heading,2,'Electrostatic potential'
  WriteReport Paragraph,
    'The electrostatic potential along a grid inside the simulation cell has been calculated '
    'for each snapshot and saved in Delphi text format in "file://(basename resultbase)(00000+firstsnapshot).phi.gz". '
    'to "file://(basename resultbase)(lastsnapshot).phi.gz". '
  WriteReport Image,(resultbase)_esp.png,Style=Figure,
    'Visualization of the electrostatic potential at each point on a grid throughout the simulation cell with transparent dots. '
    'Red dots indicate a negative, blue dots a positive potential. '
  if (espmethod=='PBS')
    WriteReport Paragraph,
      'The chosen method "PBS" uses APBS, the Adaptive Poisson-Boltzmann Solver [Baker et al., 2001], and '
      'thus the Poisson-Boltzmann equation to include solvent and counter ion effects implicitly. '
      'This method treats the cell as non-periodic. '
  else
    WriteReport Paragraph,
      'The chosen method "PME" is based on the Particle Mesh Ewald approach [Essman et al., 1995]: '
      'the atom charges are distributed on a grid, then a fast Fourier transform is made to obtain '
      'the reciprocal space portion of the Coulomb potential. Compared to the real potential, '
      'this is a smoothed representation without short-range noise and singularities [Krieger, Nielsen et al., 2006]. '
      'The resulting ESP considers all atoms in the cell and nothing else, so it is a vacuum potential without '
      'implicit solvent. If the cell is not neutral, the PME method adds a uniform "neutralizing plasma" '
      'throughout the cell to avoid artifacts. This method requires that the cell is periodic. '

if dccmsel!='' and casel!='None'
  if !dccmunits
    WriteReport Heading,2,'Dynamic Cross-Correlation Matrix'
    WriteReport Paragraph,
      'The DCCM selection (dccmsel) did not match any atoms and therefore the DCCM could not be calculated. '
      'You can change the selection by editing the `dccmsel` variable at the beginning of the macro.'
  else
    DelObj SoluteRef
    refobj=DuplicateObj Solute
    NameObj (refobj),SoluteRef
    RemoveObj SoluteRef
    i=00000+firstsnapshot
    last=0
    if count block
      snapshots=blocksnapshots
    while !last and snapshots!=0
      if format=='sim'
        sim=FileSize (MacroTarget)(i+1).sim
        if not sim
          last=1
        LoadSim (MacroTarget)(i)
      else
        last=Load(format) (MacroTarget),(i+1)
      Sim Pause
      ShowMessage 'Calculating dynamic cross-correlation matrix, analyzing snapshot (0+i)...'
      Wait 1
      Sim Off
      AddObj SoluteRef
      SupAtom (casel),(caselref)
      AddDisp(dccmsel) and Obj Solute,(dccmsel) and Obj SoluteRef
      RemoveObj SoluteRef
      if format!='sim' and !last
        for j=2 to snapshotstep
          last=Load(format) (MacroTarget),(i+j)
          if last
            break
      if snapshots!='all'
        snapshots=snapshots-1
      i=i+snapshotstep
    HideMessage
    MakeTab DCCM,Dimensions=2,Columns=(dccmunits)
    Tabulate DCCM
    coordsys=CoordSys
    pointwidth=1.
    height=5.
    dccmobj1=ShowTab DCCM,Width=(pointwidth),Range=(height),Min=-1,MinCol=(dccmcol(1)),Max=1.0,MaxCol=(dccmcol(2))
    MoveMesh (dccmobj1),Z=(-height*0.5)
    dccmobj2=ShowTab DCCM,Width=(pointwidth),Range=0,Min=-1,Max=1.0
    dccmobj3=ShowWireObj (dccmobj2),Static,Mesh=Solid
    DelObj (dccmobj2)
    PointPar Radius=0.5,Plastic=No
    NameObj (dccmobj3),ZeroLevel
    ScaleObj (dccmobj1) (dccmobj3),X=(coordsys)
    RotateObj (dccmobj1) (dccmobj3),X=180
    textwidth=pointwidth*dccmunits*2
    idlist()=List(dccmsel) and Obj Solute,Format='MOLNAME RESName RESNUM' 
    if textwidth<4500
      textobj1=MakeTextObj Units,Width=(textwidth),Height=(textwidth)
      Font Arial,Height=(pointwidth*0.6),Color=Yellow,Depth=0.5,DepthCol=Red
      for i=1 to dccmunits
        PosText X=(textwidth*0.5+pointwidth*-0.5*(dccmunits+8)),
                Y=(textwidth*0.5+pointwidth*(0.5*dccmunits-i)),justify=left
        Print (idlist(i))
      textobj2=DuplicateObj (textobj1)
      RotateObj (textobj2),Z=(90*coordsys)
    r=RadiusObj (dccmobj1)
    r=r+5
    s=PixToA*(ScrSizeY-50)*0.5/sqrt (r*r/2)
    PosObj all,Z=(EyeDis/s-EyeDis)
    DelObj not (dccmobj1) (dccmobj3) Units Solute
    if not count block
      SwitchObj Solute,off
      SwitchObj ZeroLevel,off
      SaveScreenshot 'dccm1','DCCM'
      SwitchObj Solute,on
      SwitchObj ZeroLevel,on
    SaveTab DCCM,(resultbase)_dccm(id),Format=Text,Columns=(dccmunits),NumFormat=6.3f,
            'Dynamic Cross-Correlation Matrix for (dccmunits) selected units'
    MakeTab DCCM2,Dimensions=2,Columns=(dccmunits+1)
    Tabulate "DCCM"
    Tabulate List(dccmsel) and Obj Solute,Format='MOLNAME\\n RESName\\n RESNUM' 
    idlist()=List(dccmsel) and Obj Solute,Format='MOLNAME``RESName``RESNUM'
    for i=1 to dccmunits
      Tabulate '(idlist(i))'
      Tabulate Tab DCCM,Row=(i)
    Style Ribbon
    ShowAtom (dccmsel)
    idlist()=List(dccmsel) and Obj Solute
    dccmlist()=DCCM
    ShowSolute
    Style Ribbon,Stick
    HideRes Protein
    Show(dccmsel)
    for i=1 to dccmunits
      atm1=ListAtom (idlist(i)) Visible
      for j=1 to dccmunits
        corr=dccmlist((i-1)*dccmunits+j)
        if corr<-dccmcut
          col='Blue'
        elif corr>dccmcut
          col='Red'
        else
          continue
        atm2=ListAtom (idlist(j)) Visible
        ShowArrow Start=AtAtom,(atm1),End=AtAtom,(atm2),Radius=0.1,Heads=0,Color=(col)
    SaveYOb Solute,(resultbase)_dccm(id)
    if not count block
      SaveScreenshot 'dccm2','motion correlations'
    SwitchAll On
    DelObj Solute
    NumberObj all
    PosObj all,Z=(EyeDis/s-EyeDis)
    SaveSce (resultbase)_dccm(id)
    WriteReport Heading,2,'Dynamic Cross-Correlation Matrix'
    WriteReport Paragraph,
      'The dynamic cross-correlation matrix [DCCM] is a square matrix, whose rows and columns '
      'match the selected units `(dccmsel)`. To change this selection, edit the `dccmsel` variable '
      'at the beginning of this macro. The DCCM shows how the movements of all selected pairs correlate. '
      'The values in the DCCM range from -1 [perfectly anti-correlated] to +1 [perfectly correlated]. '
      'The values along the diagonal are always +1 [because the motion of an atom is perfectly correlated to itself]. '
      'The DCCM element for units i and j is obtained with the following formula: '
    Writereport Image,Filename=(YASARADir)/doc/DCCM1.png,Style=Center,Name="formula_dccm"
    WriteReport Paragraph,
      'Here `d` is the displacement between the current position and the average position of the '
      'selected unit, and the angle brackets indicate the average over all samples. '
      'The highest correlations off the diagonal can often be found for bridged cysteines.'
    if not count block
      WriteReport Paragraph,'The image below shows the correlation directly in the solute object:'
      WriteReport Image,(resultbase)_dccm2.png,Style=Figure,
                  'Blue and red lines are shown between strongly anti- and correlated residue pairs. '
                  'To change the threshold value for the correlation lines edit the `dccmcut` variable at the beginning of this macro. '
                  'To look at this structure interactively, open the file "file://(basename resultbase)_dccm.yob" in YASARA.',Delete=yes
      WriteReport Paragraph,
        'In the image below, the DCCM is visualized with colors ranging from '
        '(dccmcol(1)) [-1, fully anti-correlated] to (dccmcol(2)) [+1, fully correlated]. '
      WriteReport Image,(resultbase)_dccm1.png,Style=Figure,
        'Visualization of the dynamic cross-correlation matrix. Open the file "file://(basename resultbase)_dccm.sce" in YASARA '
        'to look at this matrix visualization interactively. '
        'In the scene file, the zero level [0, not correlated] is indicated with a wire-frame grid. ',Delete=Yes
    else
      WriteReport Paragraph,
      'Open the file "file://(basename resultbase)_dccm(id).sce" in YASARA to look at the dynamic cross-correlation matrix visualization interactively. '
      'In the scene file, the DCCM is visualized with colors ranging from (dccmcol(1)) [-1, fully anti-correlated] to (dccmcol(2)) [+1, fully correlated] and '
      'the zero level [0, not correlated] is indicated with a wire-frame grid. '
      'To view the correlation directly in the solute object open the file "file://(basename resultbase)_dccm(id).yob" in YASARA. '
      'Blue and red lines are shown between strongly anti- and correlated residue pairs. '
      'To change the threshold value for the correlation lines edit the `dccmcut` variable at the beginning of this macro.'
    captiontext='Dynamic cross-correlation matrix. '
                'The full table is also available in text format, you need a proper text editor '
                'without line wrapping to look at this file: "file://(basename resultbase)_dccm(id).tab". '
    if tabrowsmax
      captiontext=captiontext+'Note: At most (tabrowsmax) rows of the DCCM are shown above. '
                              'Change the `tabrowsmax` variable in the macro to adjust this number. '
    WriteReport Table,DCCM2,'(captiontext)',RowsMax=(tabrowsmax)
    
Writereport Heading,2,'Additional files'
WriteReport Paragraph,'The following additional files have been created:'
Writereport Heading,3,'The main data table'
WriteReport Paragraph,
  'The main table contains all collected data in a single file. The column names match the '
  'names used above for graphs in plots and columns in tables. You can find a more detailed '
  'explanation of this table in the user manual at Recipes > Run a molecular dynamics simulation > '
  'Analyzing a trajectory. If you parse this file automatically, keep in mind that the number of '
  'columns can change any time, so you have to use the names in the first table row to find the '
  'columns of interest: "file://(basename resultbase)_analysis(id).tab"'
if (tables>1)
  Writereport Heading,3,'Extra data tables'
  WriteReport Paragraph,'User defined extra data tables:'
  for i=2 to tables
    WriteReport Paragraph,'"file://(basename resultbase)_analysis_(tabnamelist(i))(id).tab"'
if (plotrestabs>1)
  Writereport Heading,3,'Per-residue data tables'
  WriteReport Paragraph,'Data of the per-residue plots:'
  for i=1 to plotrestabs
    WriteReport Paragraph,'"file://(basename plotrestablist(i)).tab"'
Writereport Heading,3,'The structures'
paragraphtext=''
if count block
  paragraphtext=' of the current block'
WriteReport Paragraph,
  'The `time averaged structure`(paragraphtext) in PDB format: "file://(basename resultbase)_average(id).pdb"'
WriteReport Paragraph,
  'The `snapshot with the minimum solute energy`. Either just the solute in PDB format '
  '"file://(basename resultbase)_energymin(id).pdb", or the complete system including '
  'solvent as a YASARA scene "file://(basename resultbase)_energymin(id).sce".'
if last
  WriteReport Paragraph,
    'The `last snapshot` of the simulation. Either just the solute in PDB format '
    '"file://(basename resultbase)_last.pdb", or the complete system including solvent '
    'as a YASARA scene "file://(basename resultbase)_last.sce"'
if rmsdmin
  WriteReport Paragraph,'One `representative structure for each cluster` in YOb format:'
  for i=1 to clustermembers
    WriteReport Paragraph,'"file://(basename clusterfilenamelist(i)).yob"'
Writereport Heading,3,'The RMSF tables'
WriteReport Paragraph,
  'A table that lists the Root Mean Square Fluctuations [RMSFs] of all atoms in [A] is available '
  'here: "file://(basename resultbase)_rmsf(id).tab". The RMSFs have also been converted '
  'to B-factors and stored in the B-factor field of the time-average structure above. '
WriteReport Paragraph,
  'A table with average atom RMSFs per residue can be found here: "file://(basename resultbase)_rmsfres(id).tab".'
if hiresplotted and (last or not count block)
  Writereport Heading,3,'High resolution plots'
  WriteReport Paragraph,
    'To facilitate publication, high resolution versions of the plots above have been '
    'created with a 4:3 aspect ratio suited for printing in a single column of a typical '
    'journal article. Just look at the figure number above to find the right file:'
  for i=1 to 9999
    size=FileSize (resultbase)_report_figure(i)_hires.png
    if size
      WriteReport Paragraph,'"file://(basename resultbase)_report_figure(i)_hires.png"'

if last or not count block
  WriteReport End
  HideMessage
  if runWithMacro and ConsoleMode
    Exit
  elif not ConsoleMode
    ShowURL file://(resultbase)_report.html

def Plot valuelist(),title,ylabel,graphnamestr
  ShowData 'Plot',0,valuelist,title,"None",ylabel,graphnamestr

def ScatterPlot valuelist(),xvalue,title,xlabel,ylabel,graphnamestr
  ShowData 'Plot',xvalue,valuelist,title,xlabel,ylabel,graphnamestr

def WriteTable valuelist(),title,datacolumnstr
  ShowData 'Main',0,valuelist,title,"None","None",datacolumnstr

def WriteExtraTable tablename,valuelist(),datacolumnstr
  ShowData tablename,0,valuelist,"None","None","None",datacolumnstr

def ShowData resulttype,xvalue,valuelist(),title,xlabel,ylabel,graphnamestr
  global resultbase,task,hiresplotted,header,xcolumn,ycolumn,plottimestring,tabtimestring,figurewidth,tabrowsmax,snapshots,tabnamelist,tabheaderlist,tabcolslist,tables,id
  graphsmax=30
  graphcolorlist()=''
  command=''
  if type valuelist1=='StrongString'
    if count valuelist!=1
      RaiseError 'When plotting the output of a YASARA command as "(title)", exactly one command must be provided, not (count valuelist)'
    if task=='WriteReport' and valuelist1=="SecStr"
      secstrnamelist()='Helix','Sheet','Turn','Coil','Helix310','HelixPi'
      for i=1 to count secstrnamelist
        graphcolorlist(i)=ColorPar SecStr,(secstrnamelist(i))
    command=''+valuelist1
    valuelist()=(command)
  values=count valuelist
  if !values
    RaiseError 'The list of values for (resulttype) "(title)" is empty. Please check the selection'
  if values>graphsmax and resulttype=='Plot'
    RaiseError 'Plot "(title)" has (values) graphs. The maximum allowed number of graphs per plot is (graphsmax)'
  graphnamelist()=split graphnamestr
  graphnames=count graphnamelist
  if values<graphnames
    for i=graphnames to values+1 step -1
      DelVar graphnamelist(i)
  elif graphnamelist(graphnames)=='Auto'
    for i=graphnames to values
      graphnamelist(i)='Graph(i)'
  elif values>graphnames
    if resulttype=='Plot'
      valuename='graph'
    else
      valuename='data column'
    RaiseError 'You tried to add (values) (valuename)s to the (resulttype) "(title)", but provided only (graphnames) (valuename) names "(graphnamestr)". '
                'Choose "Auto" as the last or only one (valuename) name for automatic naming'
  if task=='Tabulate'
    SelectTab Main
    if resulttype!='Plot' and resulttype!='Main'
      SelectTab '(resulttype)'
    for i=1 to count valuelist
      if type valuelist(i)=='WeakString'
        Tabulate '(valuelist(i))'
      else
        Tabulate (valuelist(i))
  elif task=='AddHeader'
    table=1
    if !tables or (resulttype!='Plot' and resulttype!='Main')
      tables=tables+1
      table=tables
      tabnamelist(table)=resulttype
      if resulttype=='Plot'
        tabnamelist(table)='Main'
      tabcolslist(table)=2
      tabheaderlist(table)="    Time[ps]     Time[ns]"
    for i=1 to values
      graphname=graphnamelist(i)
      graphnamelen=strlen graphname
      if graphnamelen>12
        RaiseError 'Graph or table column name "(graphname)" is (graphnamelen) characters long, but at most 12 characters are allowed'
      tabheaderlist(table)=tabheaderlist(table)+" "*(13-graphnamelen)+graphname
    tabcolslist(table)=tabcolslist(table)+values
  elif task=='WriteReport' and (resulttype=='Plot' or resulttype=='Main')
    WriteReport Heading,3,'(title)'
    charlist()=crack graphnamelist1
    if ' ' not in charlist and '.' not in charlist
      global (graphnamelist1)Text
      if isfunction Print(graphnamelist1)Info
        Print(graphnamelist1)Info
      elif count (graphnamelist1)Text
        WriteReport Paragraph,'((graphnamelist1)Text)'
    if resulttype=='Plot'
      if xvalue
        caption=title+' [vertical axis] as a function of (graphnamelist(xvalue)) [horizontal axis]'
        DelVar graphnamelist(xvalue)
        plotvalues=values-1
        plottype='Scatter'
        plotxcolumn=ycolumn-1+xvalue
      else
        caption=title+' [vertical axis] as a function of simulation time [horizontal axis]'
        xlabel='Simulation time in (plottimestring)'
        plotvalues=values
        plottype='Line'
        plotxcolumn=xcolumn
      if command!=''
        caption=caption+', obtained with the command "(command)"'
      caption=caption+'.'
      if values>1 and !xvalue
        text=''
        textlist()=' Graph',', graph',' graph',' and graph',' has',' have'
        idx=0
        for i=1 to values
          donelist(i)=0
          columnlist()=Tab Main,(ycolumn+i-1)
          sumlist(i)=Sum columnlist
          if !sumlist(i)
            if idx<2
              idx=idx+1
            text=text+textlist(idx)+' `(graphnamelist(i))`'
        if idx
          text=text+textlist(idx+4)+' all zero values.'
        idx=0
        for i=values to 2 step -1
          if !donelist(i)
            for j=i-1 to 1 step -1
              if sumlist(j)==sumlist(i) and !donelist(j)
                donelist(j)=1
                if idx<2
                  text=text+textlist(idx+1)+' `(graphnamelist(i))` completely covers '
                  idx=3
                text=text+textlist(idx)+' `(graphnamelist(j))`'
                idx=4
            if idx
              idx=1
        if idx
          text=text+', they share the same values.'
        if text!=''
          caption=caption+' Note:'+text
      if graphnamelist1=='TotalEnergy' and snapshots>1
        tempval1=Tab Main,(ycolumn),1
        tempval2=Tab Main,(ycolumn),2
        if tempval1<tempval2
          caption=caption+' Note: The first value of the plot [(0.00+tempval1)], coming from the energy minimized starting structure, '
                          'has been replaced with the second value of the plot [(0.00+tempval2)] to show this plot with a smaller energy range '
                          'and thus a higher resolution. '
          Tab Main,(ycolumn),1,Set=(tempval2)
      WriteReport Plot,Caption='(caption)',Main,Width=(figurewidth),Height=480,Title='(title)',
                  Type=(plottype),XColumn=(plotxcolumn),YColumn=(ycolumn),YColumns=(plotvalues),XLabel=(xlabel),
                  YLabel=(ylabel),LegendPos='Outside',Graphname=(graphnamelist),Graphcol=(graphcolorlist)
      if hiresplotted
        SavePlot Filename="LastReportPlot_hires",Main,Width=1600,Height=1200,Title='(title)',
                 Type=(plottype),XColumn=(plotxcolumn),YColumn=(ycolumn),YColumns=(plotvalues),XLabel=(xlabel),
                 YLabel=(ylabel),Graphname=(graphnamelist),Graphcol=(graphcolorlist)
      if graphnamelist1=='TotalEnergy' and snapshots>1
        Tab Main,(ycolumn),1,Set=(tempval1)
    else
      caption=title+' as a function of simulation time [first column]'
      if command!=''
        caption=caption+', obtained with the command "(command)"'
      caption=caption+'.'
      if tabrowsmax
        caption=caption+' Note: At most (tabrowsmax) table rows are shown. '
                        'Change the `tabrowsmax` variable in the macro to adjust the number of shown table rows. '
                        'The full table can be found in "file://(basename resultbase)_analysis(id).tab".'
      WriteReport Table,Main,Caption='(caption)',NumFormat=.2f,RowsMax=(tabrowsmax),
                  InfoColumn=(xcolumn),DataColumn=(ycolumn),DataColumns=(values),
                  InfoColumnName='(tabtimestring)',DataColumnName=(graphnamelist)
    ycolumn=ycolumn+values

def PlotRes name,reslist(),valuelist(),valmin,valmax,title,graphnamestr,graphcolstr
  global resultbase,task,figurewidth,hiresplotted,plottimestring,simtime,xcolumn,plotrestabs,plotrestablist(),(name)Text,id
  
  resnumlist()=ListRes (join reslist),RESNUMWIC
  mollist()=ListMol (join reslist)
  start=1
  for i=1 to count mollist
    residues=CountRes (mollist(i)) (join reslist)
    if residues
      if task=='AddHeader'
        MakeTab (name)Mol(i),2,(residues+2)
        Tabulate "Time_ps","ns\ResNum"
        for j=start to start+residues-1
          Tabulate (resnumlist(j))
        start=start+residues
      elif task=='Tabulate'
        SelectTab (name)Mol(i)
        Tabulate (simtime/1000),(simtime/1000000)
        for j=start to (start+residues-1)
          Tabulate (valuelist(j))
        start=start+residues
        SelectTab Main
      elif task=='WriteReport'
        if i==1
          WriteReport Heading,3,'(title)'
          if isfunction Print(name)Info
            Print(name)Info
          elif count (name)Text
            WriteReport Paragraph,'((name)Text)'
        molname=ListMol (mollist(i)),Format='MOLNAME'
        moltitle='(title) of molecule (molname)'
        plotrestabs=plotrestabs+1
        plotrestablist(plotrestabs)='(resultbase)_plotres_(name)Mol(molname)(id)'
        if graphnamestr==''
          vals=valmax-valmin+1
          for j=1 to vals
            graphnamelist(j)='(j+valmin-1)'
          graphnames=count graphnamelist
          legendcaption=''
        else
          graphnamelist()=split graphnamestr
          graphnames=count graphnamelist
          legendcaption='. Values 1-(graphnames) in the table correspond to the (graphnames) labels in the plot legend.'
        if (count graphnamelist==2 and graphnamelist1=='Min' and graphnamelist2=='Max')
          used=1
          graphnamelist1='Min (valmin)'
          graphnamelist2='Max (valmax)'
          caption='data'
          SaveTab (name)Mol(i),(plotrestablist(plotrestabs)),NumFormat=12.3f,(moltitle)
        else
          used=0
          cap()='',''
          shift=1-valmin
          tablist()=Tab (name)Mol(i)
          cols=residues+2
          rows=count tablist/cols
          for j=1 to graphnames*cols
            percentlist(j)=0
          for j=2 to rows
            for k=1 to cols
              if k<3
                time=Tab (name)Mol(i),(k),(j)
                time=time//1
                Tab (name)Mol(i),(k),(j),Set=(time)
              else
                idx=(j-1)*cols+k
                if tablist(idx)<valmin
                  tablist(idx)=valmin
                  cap1='<='
                elif tablist(idx)>valmax
                  tablist(idx)=valmax
                  cap2='>='
                tablist(idx)=tablist(idx)+shift
                if tablist(idx)>1
                  used=1
                num=k+(tablist(idx)-1)*cols
                percentlist(num)=1.+percentlist(num)
          if graphnamelist1=='(valmin)'
            graphnamelist1='(cap1)(valmin)'
          if graphnamelist(graphnames)=='(valmax)'
            graphnamelist(graphnames)='(cap2)(valmax)'
          for j=1 to graphnames
            idx=(j-1)*cols
            percentlist(idx+1)=""+'(graphnamelist(j))'
            percentlist(idx+2)="%"
            for k=3 to cols
              percentlist(idx+k)=round (percentlist(idx+k)*100./(rows-1))
          SelectTab (name)Mol(i)
          Tabulate (percentlist)
          caption='data including percentages'
          SaveTab (name)Mol(i),(plotrestablist(plotrestabs)),NumFormat=10.0f,(moltitle)
          DelTab (name)Mol(i)
          MakeTab (name)Mol(i),2,(cols)
          tablist1=""+'(tablist1)'
          tablist2=""+'(tablist2)'
          Tabulate (tablist)
          SelectTab Main
        if used
          WriteReport Plot,'(title) as a function of simulation time [horizontal axis] for each residue number [vertical axis]. '
                           'A table with the raw (caption) is available here: "file://(basename plotrestablist(plotrestabs)).tab"(legendcaption)',
                      (name)Mol(i),Width=(figurewidth),Height=480,Title=(moltitle),Type=Heatmap,
                      XColumn=(xcolumn),YColumn=3,YColumns=(residues),XLabel='Simulation time in (plottimestring)',
                      YLabel='Residue number',LegendPos='Outside',Graphname=(graphnamelist),Graphcol=(split graphcolstr)
          if hiresplotted
            SavePlot Filename="LastReportPlot_hires",(name)Mol(i),Width=1600,Height=1200,Title=(moltitle),Type=Heatmap,
                     XColumn=(xcolumn),YColumn=3,YColumns=(residues),XLabel='Simulation time in (plottimestring)',
                     YLabel='Residue number',Graphname=(graphnamelist),Graphcol=(split graphcolstr)

def ShowSystem
  global ligandsel 
  SwitchAll On
  alpha,beta,gamma = OriObj SimCell
  RotateAll Y=(-beta)
  RotateAll Z=(-gamma)
  RotateAll X=(-alpha)
  Style BallStick
  Style Ribbon,Stick
  ColorBonds Order
  if ligandsel!=''
    BallRes (ligandsel)
  ZoomAll Steps=0

def ShowSolute
  global ligandsel 
  SwitchAll Off
  SwitchObj Solute,On
  NiceOriObj Solute
  Style BallStick
  BallAtom all with 0 bonds to all
  if ligandsel!=''
    BallRes (ligandsel)
    CenterObj Solute
    TransformObj Solute
    _,_,cenz = GroupCenter (ligandsel) Obj Solute
    if cenz>0
      RotateObj Solute,Y=180
  ZoomObj Solute,Steps=0
  
def SaveScreenshot fileid,description
  global figurewidth,resultbase
  ShowMessage 'Creating ray-traced picture of the (description)...'
  Wait 1
  RayTrace Filename=(resultbase)_(fileid).png,X=(figurewidth),Zoom=1,LabelShadow=No,Display=Off,Outline=0,Background=On

def ShortList longdatalist(),longidlist(),shortidlist()
  ids=count longidlist
  idx=1
  for i=1 to count shortidlist
    while longidlist(idx)!=shortidlist(i)
      idx=idx+1
    shortdatalist(i)=longdatalist(idx)
  return (shortdatalist)
