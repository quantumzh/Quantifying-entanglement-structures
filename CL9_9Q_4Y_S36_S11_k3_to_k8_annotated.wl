(* ::Package:: *)

(* CL9 9-qubit Mathematica script *)
(* Converted from CL9_9Q_4Y_S36_S11_k3_to_k8.nb *)
(* Includes S36 d=1..4, S11 d=5..9, k-separability k=3..8. *)
(* File paths for real experimental data must be configured locally. *)
(* Run with Get["path/to/CL9_9Q_4Y_S36_S11_k3_to_k8.wl"]. *)
(* NOTE: Mathematical implementations retained; corrected two S36 output field labels to Upsilon. Not kernel-validated. *)

(* ==== Input cell 1 ==== *)
(*
SECTION 1 | Quantum state, measurement settings, and Omega
Paper correspondence:
  N = 9                    -> n
  j = 4                    -> j  (exactly four Pauli Y operators)
  h                        -> h  (maximum number of swapped blocks)
  Omega                    -> omega, four ordered original pairs (1,8),(2,7),(3,6),(4,5)
  |m> of each (|m>,|n>)    -> ms
  X/Y/Z encoding           -> 1/2/3
  rowsByY[j]               -> original 1-based measurement row numbers for j Ys
  pauliCodes[j]            -> nine-entry local Pauli strings, lexicographic order
  rowCodes[row]            -> observed Pauli string on that row

IMPORTANT: the originally named first group contains EIGHT Ys and is
stored as rowsByY[8]. Original row r corresponds to file
basis_measurement(r-1).txt. row 340 is Z^tensor9.
*)

ClearAll["Global`*"];
n=9; h=1; j=4;
rowsByY=<|2->{2,3,4,6,7,8,10,11,12,14,19,20,21,23,28,40,41,42,44,49,58,79,80,81,83,88,97,118,157,158,159,161,166,175,196,235},4->{9,13,15,16,17,22,24,25,26,29,30,31,33,34,35,43,45,46,47,50,51,52,54,55,56,59,60,61,63,64,65,67,68,69,71,82,84,85,86,89,90,91,93,94,95,98,99,100,102,103,104,106,107,108,110,119,120,121,123,124,125,127,128,129,131,136,137,138,140,145,160,162,163,164,167,168,169,171,172,173,176,177,178,180,181,182,184,185,186,188,197,198,199,201,202,203,205,206,207,209,214,215,216,218,223,236,237,238,240,241,242,244,245,246,248,253,254,255,257,262,274,275,276,278,283,292},6->{36,57,66,70,72,73,74,96,105,109,111,112,113,126,130,132,133,134,139,141,142,143,146,147,148,150,151,152,174,183,187,189,190,191,204,208,210,211,212,217,219,220,221,224,225,226,228,229,230,243,247,249,250,251,256,258,259,260,263,264,265,267,268,269,277,279,280,281,284,285,286,288,289,290,293,294,295,297,298,299,301,302,303,305},8->{153,231,270,291,300,304,306,307,308}|>;
pauliCodes=Association@Table[
  yy->Sort[ReplacePart[ConstantArray[1,n],
     Thread[#->ConstantArray[2,yy]]]&/@Subsets[Range[n],{yy}]],
  {yy,{2,4,6,8}}];
rowCodes=Association[Flatten[Table[
 Thread[rowsByY[yy]->pauliCodes[yy]],{yy,{2,4,6,8}}],1]];
rowCodes[1]=ConstantArray[1,9];
rowCodes[340]=ConstantArray[3,9];
If[!And@@Table[Length[rowsByY[k]]==Binomial[9,k],{k,{2,4,6,8}}]||
 !DuplicateFreeQ[Flatten[Values[rowsByY]]]||
 !DuplicateFreeQ[Values[rowCodes]],Print["Pauli mapping error"];Abort[]];
stateStrings={"000000000","000000111","000111011","000111100",
 "111000011","111000100","111111000","111111111"};
toBits[s_String]:=ToExpression/@Characters[s];
states=toBits/@stateStrings;
omega=states[[#]]&/@{{1,8},{2,7},{3,6},{4,5}};
ms=omega[[All,1]];
Print["Pauli row counts: ",Length/@rowsByY];
Print["Four pairs: ",omega//InputForm];


(* ==== Input cell 2 ==== *)
(*
SECTION 2 | Import local 341 measurement files / repetition
  folders                  -> 10 separately acquired experimental runs
  loadRep[folder]          -> 341 settings x 512 outcomes (or Failure)
  probabilities[measure]   -> normalized 512-outcome distribution of one setting
  correlators[block]       -> corr[[row]] = <P_row> = parity-weighted sum of outcomes
  computationalProbabilities -> p_b = <b|rho|b> from all-Z setting

The experiment's outcome index is x=1+Sum[(iq-1)2^(q-1)], q=1..9:
qubit 1 is the least significant bit, so binIndex reverses the
written bit-string before FromDigits. WARNING: paths rep0..rep9 are
PLACEHOLDERS; change them to actual folders before real measurements.
*)

(* IMPORTANT: replace folders with the ten ACTUAL separate folders.
   Originally all ten paths were rep0, which would duplicate data. *)
baseDir="E:\\Detection of k partite entanglement\\Experimental_special_bases\\9qcluster_special_bases\\special_bases";
folders=Table[FileNameJoin[{baseDir,"rep"<>ToString[k],
 "9Q_Cluster_1500ns","readout_cali"}],{k,0,9}];
(* Edit folders[[k]] if the actual repetition folder differs. *)
ClearAll[loadRep,probabilities,correlators];
loadRep[folder_String]:=Module[{files,missing},
 files=Table[FileNameJoin[{folder,"basis_measurement"<>ToString[t]<>".txt"}],
   {t,0,340}];
 missing=Select[files,!FileExistsQ[#]&];
 If[missing=!={},Return[Failure["MissingFiles",
  <|"Folder"->folder,"MissingCount"->Length[missing],
    "FirstMissing"->First[missing]|>]]];
 Import[#,"Table"]&/@files];
(* Qubit 1 is LSB of the outcome number x-1, as in user's nested loops. *)
paritySigns=Table[(-1)^DigitCount[x,2,1],{x,0,511}];
probabilities[measure_List]:=Module[{v},
 If[Length[measure]!=512||!And@@(Length[#]>=1&/@measure),
  Return[Failure["WrongOutcomeCount",<||>]]];
 v=N[measure[[All,1]]];
 If[!VectorQ[v,NumericQ]||Min[v]<0||Total[v]<=0,
  Return[Failure["BadProbabilities",<||>]]];
 v/Total[v]];
correlators[block_List]:=Module[{p},
 If[Length[block]!=341,Return[Failure["WrongSettingCount",<||>]]];
 Table[p=probabilities[block[[r]]];
  If[Head[p]===Failure,Return[p]];
  paritySigns.p,{r,1,341}]];
(* Row 340 corresponds to basis_measurement339.txt and Z^9.
   This distribution is indexed Q1 as the least significant bit. *)
computationalProbabilities[block_List]:=probabilities[block[[340]]];


(* ==== Input cell 3 ==== *)
(*
SECTION 3 | Witness observables and paper's f_j
  signForCode[m,code]       -> s_S(m)=(-1)^(Sum_{i in S} m_i)
  wFromCorrelations        -> <W_j^9(m)> = Sum_{|S|=j}s_S(m)<P_S>
  allW                     -> one signed <W_j^9(m)> for EACH of four Omega pairs
  Total[Abs[allW[...]]]    -> Sum_{(m,n) in Omega}|<W_j^9(m)>|
  pairPopulation[p]        -> Sum_Omega sqrt(p_m p_n)
  fPaper[j,N]              -> paper's combinatorial f_j(N)

Take absolute values AFTER each signed sum, not term-by-term.
p and corr are from different measured bases and must correspond
to the SAME experimental repetition.
*)

ClearAll[signForCode,wFromCorrelations,allW,binIndex,pairPopulation,fPaper];
signForCode[m_List,code_List]:=(-1)^Total[Pick[m,code,2]];
wFromCorrelations[m_List,jj_Integer,observations_List]:=
 Total[MapThread[
   signForCode[m,#2]*observations[[#1]]&,
   {rowsByY[jj],pauliCodes[jj]}]];
allW[jj_Integer,observations_List]:=
 wFromCorrelations[#,jj,observations]&/@ms;
(* x=1+(i1-1)+2(i2-1)+...+256(i9-1) *)
binIndex[b_List]:=FromDigits[Reverse[b],2]+1;
pairPopulation[p_List]:=Total[
  Sqrt[p[[binIndex[#[[1]]]]]*p[[binIndex[#[[2]]]]]]&/@omega];
fPaper[2,3]=1;fPaper[2,4]=2;
fPaper[2,k_Integer]/;k>=5:=(k-1)(k-4)/2;
fPaper[3,5]=2;fPaper[3,6]=4;
fPaper[3,k_Integer]/;k>=7:=Sum[fPaper[2,u],{u,6,k-1}];
fPaper[4,5]=3;fPaper[4,6]=5;fPaper[4,7]=10;fPaper[4,8]=14;
fPaper[4,k_Integer]/;k>=9:=Sum[fPaper[3,u],{u,8,k-1}];
Print["f4(9) = ",fPaper[4,9],"; Binomial[9,4] = ",Binomial[9,4]];


(* ==== Input cell 4 ==== *)
(*
SECTION 4 | Weighted swap sets: Notation S3, S4, S5, S6
Paper symbols and names in THIS implementation:
  gamma (partition)                -> part
  D_d^* (allowed d-depth partitions)-> minPartitions[d]
  I_gamma^h (selected unions of blocks) -> subsetChoices[part,hh]
  subsystem exchange              -> partialSwap[pair,sites]
  Gamma_gamma^(Omega,h)           -> gammaUpsilonForPartition[...,hh][[1]]
  Upsilon_gamma^(Omega,h)         -> gammaUpsilonForPartition[...,hh][[2]]
  Gamma_dD^(Omega,h)              -> gammaD[d,hh]
  Upsilon_dD^(Omega,h)            -> upsilonD[d,hh]
  Upsilon_1D^(Omega,Floor[N/2])   -> refUpsilon[]
  overbar(Upsilon_dD^(Omega,h))   -> hatUpsilon[d,hh]
  C_d^h                            -> cDh[d,hh]

Single swapped pair (a,b) contributes:
  Gamma entries (a,1/2) and (b,1/2) [single-state weights];
  Upsilon entry (a,b,1) [joint-pair weight].
"Addition union" (S5,S6) -> Merge[...,Total] within each partition,
across selected subsystems and Omega.
"Maximum union" (S3,S4) -> associationMax across allowed partitions.
The set difference in hatUpsilon removes every overlapping KEY,
not merely subtracts its coefficient.

KEY VALIDATION CAVEAT: minPartitions[d] chooses partitions having
Ceiling[9/d] blocks. Check that this is exactly D_d^* in the paper.
The half-subset convention for h=|gamma|/2 and any special symmetry
normalization should also be verified against an author test case.
S3-S6 arithmetic has been encoded but not Mathematica-kernel validated.

  epsGamma[weightedSingles,p] -> Sum_(a,x) x p_a
  epsUpsilon[weightedPairs,p] -> Sum_(a,b,x) x sqrt(p_a p_b)
For S36 (as pictured by the user), BOTH RHS epsilon terms use Upsilon.
The distinct Gamma terms are needed in the high-depth S11 branch.
*)

ClearAll[boundedPartitions,minPartitions,subsetChoices,partialSwap,
 stateKey,jointKey,decodeState,decodeJoint,
 gammaUpsilonForPartition,associationAdd,associationMax,
 gammaD,upsilonD,refUpsilon,hatUpsilon,
 cDh,epsGamma,epsUpsilon];
boundedPartitions[0,d_Integer?Positive]:={{}};
boundedPartitions[nn_Integer?Positive,d_Integer?Positive]:=
 boundedPartitions[nn,d]=Module[{prev=boundedPartitions[nn-1,d]},
 Join[(Append[#,{nn}]&/@prev),
 Flatten[Table[If[Length[q[[k]]]<d,
   {ReplacePart[q,k->Append[q[[k]],nn]]},{}],
   {q,prev},{k,Length[q]}],2]]];
(* D*_d: partitions with the minimal number of blocks. *)
minPartitions[d_Integer?Positive]:=
 minPartitions[d]=Select[boundedPartitions[9,d],
    Length[#]==Ceiling[9/d]&];
(* I_gamma^h: choose u blocks, for 1<=u<=h; each subset is ONE B. *)
subsetChoices[part_List,hh_Integer?Positive]:=
 Flatten[Table[Flatten/@Subsets[part,{u}],
   {u,1,Min[hh,Length[part]]}],1];
partialSwap[pair_List,sites_List]:=Module[
 {a=pair[[1]],b=pair[[2]],temp},
 temp=a[[sites]];a[[sites]]=b[[sites]];
 b[[sites]]=temp;{a,b}];
stateKey[s_List]:=StringJoin[ToString/@s];
jointKey[ab_List]:=StringRiffle[stateKey/@Sort[ab],"|"];
decodeState[k_String]:=toBits[k];
decodeJoint[k_String]:=toBits/@StringSplit[k,"|"];
associationAdd[arr_List]:=Module[{keys=Union[Flatten[Keys/@arr]]},
 AssociationThread[keys,
  Table[Total[Lookup[arr,key,0]],{key,keys}]]];
associationMax[arr_List]:=Module[{keys=Union[Flatten[Keys/@arr]]},
 AssociationThread[keys,
  Table[Max[Lookup[arr,key,0]],{key,keys}]]];
gammaUpsilonForPartition[part_List,hh_Integer?Positive]:=
 Module[{choices=subsetChoices[part,hh],pairs,gEntries,uEntries},
 pairs=Flatten[Table[partialSwap[ab,b],
     {ab,omega},{b,choices}],1];
 gEntries=Flatten[({stateKey[#[[1]]]->1/2,
                   stateKey[#[[2]]]->1/2}&/@pairs),1];
 uEntries=(jointKey[#]->1&/@pairs);
 {
  Merge[gEntries,Total],
  Merge[uEntries,Total]
 }];
gammaD[d_Integer?Positive,hh_Integer?Positive]:=
 gammaD[d,hh]=associationMax[
  (gammaUpsilonForPartition[#,hh][[1]]&/@minPartitions[d])];
upsilonD[d_Integer?Positive,hh_Integer?Positive]:=
 upsilonD[d,hh]=associationMax[
  (gammaUpsilonForPartition[#,hh][[2]]&/@minPartitions[d])];
(* reference Upsilon(1D, Floor[N/2]): singleton blocks, up to four chosen *)
refUpsilon[]:=refUpsilon[]=gammaUpsilonForPartition[
 ({#}&/@Range[9]),Floor[9/2]][[2]];
(* Set difference: fully remove shared pair keys; do not subtract weights *)
hatUpsilon[d_Integer?Positive,hh_Integer?Positive]:=
 hatUpsilon[d,hh]=KeySelect[refUpsilon[],
   !KeyExistsQ[upsilonD[d,hh],#]&];
cDh[d_Integer?Positive,hh_Integer?Positive]:=
 cDh[d,hh]=Min[
  (Sum[Binomial[Length[#],u],{u,1,Min[hh,Length[#]]}]&
    /@minPartitions[d])];
epsGamma[c_Association,p_List]:=Total[KeyValueMap[
  Function[{k,x},x p[[binIndex[decodeState[k]]]]],c]];
epsUpsilon[c_Association,p_List]:=Total[KeyValueMap[
  Function[{k,x},With[{ab=decodeJoint[k]},
    x Sqrt[p[[binIndex[ab[[1]]]]]*p[[binIndex[ab[[2]]]]]]]],c]];


(* ==== Input cell 5 ==== *)
(*
SECTION 5 | Proposition S2 Eq. (S36), N=9, j=4, d=1..4
Paper: C_d^h c_L <= c_R.
  W values  -> wv
  f_4(9)    -> f (=14)
  C_d^h     -> C
  Binomial[9,4]+f_4(9)(C_d^h-1) -> den
  Sum_Omega sqrt(p_m p_n)        -> pairPop
  epsilon_D^(Omega,h)           -> e  (IMPORTANT: Upsilon-based)
  epsilonHat_D^(Omega,h)        -> eh (complement Upsilon-based)
  LEFT = C_d^h c_L              -> left
  RIGHT = c_R                   -> right
  LEFT-RIGHT                    -> Delta
  positive Delta                -> violation of the tested d-producible bound

LEFT = C*(Sum_Omega Abs[W] - 2*f*pairPop)/(2*den).
RIGHT = e + (f*C/den)*eh.
If the code uses epsGamma in e, it must be corrected to epsUpsilon;
this annotated edition retains only the S36 RHS correction already
made in the previous edition.
*)

ClearAll[oneS36,allS36];
oneS36[d_Integer,hh_Integer,corr_List,p_List]:=
 Module[{wv,f,C,den,left,right,e,eh,pairPop},
 If[!(1<=d<=4 && hh>=1 && Length[corr]>=340 &&
   Length[p]==512 && VectorQ[corr,NumericQ]&&VectorQ[p,NumericQ]),
   Return[Failure["BadInput",<||>]]];
 wv=allW[4,corr];f=fPaper[4,9];C=cDh[d,hh];
 pairPop=pairPopulation[p];den=Binomial[9,4]+f(C-1);
 left=N[C (Total[Abs[wv]]-2 f pairPop)/(den)];
 e=N[epsUpsilon[upsilonD[d,hh],p]];
 eh=N[epsUpsilon[hatUpsilon[d,hh],p]];
 right=N[2e+2f*C/den*eh];
 <|"d"->d,"h"->hh,
   "W(1,8)"->N[wv[[1]]],"W(2,7)"->N[wv[[2]]],
   "W(3,6)"->N[wv[[3]]],"W(4,5)"->N[wv[[4]]],
   "SumAbsW"->N[Total[Abs[wv]]],"Cdh"->C,
   "PairPopulationSum"->N[pairPop],
   "EpsilonUpsilon"->e,"EpsilonHatUpsilon"->eh,
   "LEFT"->left,"RIGHT"->right,"Delta"->N[left-right],
   "ViolationQ"->TrueQ[left>right]|>
 ];
allS36[hh_Integer,corr_List,p_List]:=
 Dataset[Table[oneS36[d,hh,corr,p],{d,1,4}]];


(* ==== Input cell 6 ==== *)
(*
SECTION 6 | S36 repeated experiment driver
  allRepetitions[folders,h] -> for each repetition and d=1..4:
  oneS36[d,h,corr,p] -> LEFT, RIGHT, Delta, ViolationQ.
Not a k-separability or S11 calculation; those appear below.
*)

ClearAll[allRepetitions];
allRepetitions[folderList_List,hh_Integer:1]:=
 Module[{results={},raw,corr,p,ans},
 Do[
  raw=loadRep[folderList[[rep]]];
  If[Head[raw]===Failure,
   AppendTo[results,<|"Repeat"->rep,"Error"->raw|>];
   Continue[]];
  corr=correlators[raw];p=computationalProbabilities[raw];
  If[Head[corr]===Failure||Head[p]===Failure,
   AppendTo[results,<|"Repeat"->rep,
     "Error"->"Invalid measurement file format"|>];
   Continue[]];
  Do[
   ans=oneS36[d,hh,corr,p];
   AppendTo[results,Join[<|"Repeat"->rep|>,ans]],
   {d,1,4}],
  {rep,1,Length[folderList]}];
 Dataset[results]];
(* Actual experiment:
results4Y=allRepetitions[folders,1];
results4Y
*)


(* ==== Input cell 7 ==== *)
(*
SECTION 7 | Ideal-state diagnostic example (NOT experiment)
  idealKet -> |CL_9> with eight amplitudes 1/(2 sqrt(2)).
  idealP -> p_b=|<b|CL_9>|^2.
  idealPauli[code] -> calculated <CL_9|P_code|CL_9>.
  idealCorr -> 341-element vector (only specified 4Y rows set);
  allW[4,idealCorr] -> four signed joint-observable expectations.
The result is a theoretical consistency test, not experimental evidence.
*)

idealKet=ConstantArray[0,512];
Do[idealKet[[binIndex[toBits[s]]]]=1/(2 Sqrt[2]),{s,stateStrings}];
idealP=Abs[idealKet]^2;
idealPauli[code_List]:=Module[{flipSites=Flatten[Position[code,1|2]],
 yzSites=Flatten[Position[code,2|3]],yy=Count[code,2]},
 FullSimplify[Sum[Module[{b=IntegerDigits[x,2,9],nb,phase},
  nb=b;nb[[flipSites]]=1-nb[[flipSites]];
  phase=I^yy (-1)^Total[b[[yzSites]]];
  Conjugate[idealKet[[binIndex[nb]]]]*phase*
      idealKet[[binIndex[b]]]],
 {x,0,511}]]];
idealCorr=ConstantArray[0,341];
Do[idealCorr[[rowsByY[4][[k]]]]=idealPauli[pauliCodes[4][[k]]],
 {k,126}];
Print["CL9 norm squared = ",Conjugate[idealKet].idealKet];
Print["four signed W values = ",allW[4,idealCorr]];
Print["four pair-population root terms = ",pairPopulation[idealP]];
Print["IDEAL 4Y S36 LEFT / RIGHT (h=1):"];
allS36[1,idealCorr,idealP]


(* ==== Input cell 8 ==== *)
(*
SECTION 8 | Eq. (S11): high depth d=5..9
Paper's notation (for ellOmega=9, j=4, f_4(9)=14):
  2c_L = (Sum_Omega |<V_4|J^9>| - 28 Sum_Omega sqrt(p_m p_n))/126
  2c_R = epsilon_D^Omega + (14/126) epsilonHat_D^Omega
      = epsilon_D^Omega + (1/9) epsilonHat_D^Omega.
  referenceGamma1[]       -> Gamma_(1D)^(Omega,1)
  gammaHatS11[d]          -> Gamma_(N-1,D)^(Omega,1) \ Gamma_dD^(Omega,1)
  epsGamma[... ,p]        -> sum_(a,x) x p_a
  oneS11[d,corr,p]        -> S11 LEFT, RIGHT, Delta
  compareAllDepths[...]   -> S36(d=1..4) + S11(d=5..9)
  allRepetitionsBoth[...] -> repeated experimental outputs for both

CAUTION: d=9 and the exact Gamma membership/normalization require
full-paper verification; this is inherited from the preceding file.
Do not interpret unverified numerical violation as certification.
*)

ClearAll[referenceGamma1, gammaHatS11, oneS11, allS11,
         compareAllDepths, allRepetitionsBoth];
(* Reuse gammaD and associated operations from Section 4. *)
referenceGamma1[] := referenceGamma1[] = gammaD[1,1];
(* Paper: overlined Gamma_dD = Gamma_(N-1,D) \ Gamma_dD.
   The set difference removes keys in Gamma_dD entirely. *)
gammaHatS11[d_Integer] :=
 gammaHatS11[d] =
   KeySelect[gammaD[8,1],!KeyExistsQ[gammaD[d,1],#]&];

oneS11[d_Integer,corr_List,p_List] :=
 Module[{wv,L,e,eHat,R},
  If[!(5<=d<=9 && Length[corr]>=340 && Length[p]==512
      &&VectorQ[corr,NumericQ]&&VectorQ[p,NumericQ]),
   Return[Failure["BadS11Input",<||>]]];
  wv=allW[4,corr];
  L=N[(Total[Abs[wv]]-28 pairPopulation[p])/126];
  e=N[epsGamma[referenceGamma1[],p]];
  eHat=N[epsGamma[gammaHatS11[d],p]];
  R=N[e+(1/9)eHat];
  <|"d"->d,"Equation"->"S11","W(1,8)"->N[wv[[1]]],
   "W(2,7)"->N[wv[[2]]],"W(3,6)"->N[wv[[3]]],
   "W(4,5)"->N[wv[[4]]],
   "LEFT"->L,"EpsilonGamma1"->e,"EpsilonGammaHat"->eHat,
   "RIGHT"->R,"Delta"->N[L-R],
   "ViolationQ"->TrueQ[L>R],
   "GammaReferenceKeys"->Length[referenceGamma1[]],
   "GammaHatKeys"->Length[gammaHatS11[d]]|>
 ];
allS11[corr_List,p_List] :=
 Dataset[Table[oneS11[d,corr,p],{d,5,9}]];

compareAllDepths[hh_Integer,corr_List,p_List] :=
 Dataset[Join[
  Table[Join[<|"Equation"->"S36"|>,
    oneS36[d,hh,corr,p]],{d,1,4}],
  Table[oneS11[d,corr,p],{d,5,9}]]];

allRepetitionsBoth[folderList_List,hh_Integer:1] :=
 Module[{records={},raw,corr,p},
 Do[
  raw=loadRep[folderList[[rep]]];
  If[Head[raw]===Failure,
    AppendTo[records,<|"Repeat"->rep,"Error"->raw|>];
    Continue[]];
  corr=correlators[raw];p=computationalProbabilities[raw];
  If[Head[corr]===Failure||Head[p]===Failure,
    AppendTo[records,<|"Repeat"->rep,
      "Error"->"Bad experiment file"|>];
    Continue[]];
  records=Join[records,Table[
     Join[<|"Repeat"->rep,"Equation"->"S36"|>,
       oneS36[d,hh,corr,p]],{d,1,4}],
    Table[Join[<|"Repeat"->rep|>,
       oneS11[d,corr,p]],{d,5,9}]],
  {rep,Length[folderList]}];
 Dataset[records]];


(* ==== Input cell 9 ==== *)
(* SECTION 9 | Execute ideal-state comparisons for S11 and S36.
These cells output numbers when evaluated with Mathematica; the script
does not contain any real experimental observation values. *)

(* Run Section 7 once to define idealCorr and idealP. *)
Print["S11 high-depth comparison (ideal CL9):"];
allS11[idealCorr,idealP]
Print["S36 + S11 comparison (ideal CL9):"];
compareAllDepths[1,idealCorr,idealP]
(* To read 10 independent real experimental repetitions:
  resultsBoth=allRepetitionsBoth[folders,1];
  resultsBoth
*)


(* ==== Input cell 10 ==== *)
(*
SECTION 10 | 3...8 separability: k-partition version of Eq. S36
Paper: replace D_d^* by S_k^* and 1D reference by 1S reference.
  S_k^*                     -> kPartitions9[k] (exactly k nonempty blocks)
  Upsilon_kS^(Omega,h)      -> kUpsilon9[k,h]
  Upsilon_1S^(Omega,Floor[N/2]) -> kReference9[]
  overbar(Upsilon_kS)       -> kHat9[k,h]
  C_k^h                     -> kC9[k,h] = Sum_{u=1}^h Binomial[k,u]
  epsilon_kS                -> epsUpsilon[kUpsilon9[k,h],p]
  epsilonHat_kS             -> epsUpsilon[kHat9[k,h],p]
  oneKS9                    -> LEFT,RIGHT,Delta for a fixed k
  allKS9                    -> k=3,...,8
  allRepeatedKS9            -> all experimental repetitions

LEFT=C_k^h*(SumOmega|W|-28*PairPopulation)/
                  (252+28(C_k^h-1)).
RIGHT=epsilon_kS + 14*C_k^h/(126+14(C_k^h-1))*epsilonHat_kS.
For k-separability this is a DIFFERENT question than entanglement
depth d. k=3..8 denotes number of tensor-product blocks.
*)

ClearAll[kPartitions9,kUpsilon9,kReference9,kHat9,kC9,
  oneKS9,allKS9,allRepeatedKS9];
(* Exactly k blocks, unconstrained block sizes. *)
kPartitions9[k_Integer?Positive] /; 1<=k<=9 :=
 kPartitions9[k] =
   Select[boundedPartitions[9,9],Length[#]==k&];

(* Each partition uses the previously defined Gamma/Upsilon generator.
   Only joint-pair Upsilon component [[2]] is used. *)
kUpsilon9[k_Integer,hh_Integer] :=
 kUpsilon9[k,hh] =
 associationMax[
   (gammaUpsilonForPartition[#,hh][[2]] & /@ kPartitions9[k])];

(* Special reference 1-separable construction: nine singleton blocks,
   allow up to Floor[9/2]=4 swapped blocks. *)
kReference9[] :=
 kReference9[] =
 gammaUpsilonForPartition[({#}&/@Range[9]),Floor[9/2]][[2]];

(* Set difference: delete ALL shared pair keys, never subtract weights. *)
kHat9[k_Integer,hh_Integer] :=
 kHat9[k,hh] =
 KeySelect[kReference9[],!KeyExistsQ[kUpsilon9[k,hh],#]&];

kC9[k_Integer,hh_Integer] :=
 Sum[Binomial[k,u],{u,1,Min[k,hh]}];

oneKS9[k_Integer,hh_Integer,corr_List,p_List] :=
 Module[{wv,C,f=14,den,left,right,e,eHat,pairPop},
  If[!(3<=k<=8&&hh>=1&&Length[corr]>=340&&
        Length[p]==512&&VectorQ[corr,NumericQ]&&
        VectorQ[p,NumericQ]&&Min[p]>=0),
   Return[Failure["BadInput",<||>]]];
  wv=allW[4,corr];
  C=kC9[k,hh];den=126+f(C-1);
  pairPop=pairPopulation[p];
  left=N[C*(Total[Abs[wv]]-2 f pairPop)/(den)];
  e=N[epsUpsilon[kUpsilon9[k,hh],p]];
  eHat=N[epsUpsilon[kHat9[k,hh],p]];
  right=N[2e+(2f C/den)*eHat];
  <|"k"->k,"h"->hh,"C_k_h"->C,"SumAbsW"->N[Total[Abs[wv]]],
    "Epsilon"->e,"EpsilonHat"->eHat,
    "LEFT"->left,"RIGHT"->right,
    "Delta"->left-right,"ViolationQ"->TrueQ[left>right],
    "PartitionCount"->Length[kPartitions9[k]]|>
 ];
allKS9[hh_Integer,corr_List,p_List] :=
 Dataset[Table[oneKS9[k,hh,corr,p],{k,3,8}]];

allRepeatedKS9[folderList_List,hh_Integer:1] :=
 Module[{rows={},raw,corr,p},
  Do[
   raw=loadRep[folderList[[rep]]];
   If[Head[raw]===Failure,
    AppendTo[rows,<|"Repeat"->rep,"Error"->raw|>];Continue[]];
   corr=correlators[raw];
   p=computationalProbabilities[raw];
   If[Head[corr]===Failure||Head[p]===Failure,
    AppendTo[rows,<|"Repeat"->rep,
      "Error"->"Invalid measured rows"|>];Continue[]];
   rows=Join[rows,Table[
     Join[<|"Repeat"->rep|>,oneKS9[k,hh,corr,p]],
     {k,3,8}]],
   {rep,Length[folderList]}];
 Dataset[rows]
 ];


(* ==== Input cell 11 ==== *)
(* SECTION 11 | Examples / experiment calls
  allKS9[1,idealCorr,idealP] and allKS9[2,...] = ideal k-tests.
  allRepeatedKS9[folders,1|2] = measured k-tests.
No claim that the ideal example or code was executed in a Wolfram kernel.
*)

(* First evaluate the earlier notebook cells to define idealCorr and idealP. *)
Print["k-separability (ideal CL9, h=1):"];
allKS9[1,idealCorr,idealP]
Print["k-separability (ideal CL9, h=2):"];
allKS9[2,idealCorr,idealP]

(* Experiment: supply ten real, DISTINCT folders before running. *)
(* measuredKS_h1=allRepeatedKS9[folders,1]; *)
(* measuredKS_h2=allRepeatedKS9[folders,2]; *)
(* measuredKS_h2 *)

(* Show the corrected S36 d-depth rows, both epsilon terms now from Upsilon *)
(* allS36[2,idealCorr,idealP] *)

