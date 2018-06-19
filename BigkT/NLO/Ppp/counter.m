(* ::Package:: *)

If[ $FrontEnd === Null,
		$FeynCalcStartupMessages = False;
		Print["Computation of the matrix element squared for the q_i qbar_i -> g g scattering in QCD at tree level"];
];
If[$Notebooks === False, $FeynCalcStartupMessages = False];
$LoadFeynArts= True;
<<FeynCalc`
$FAVerbose = 0;
subcnt={DiracGamma[6]->1/2,DiracGamma[7]->1/2,dMf1[3,1]->0,dZAA1->0,dZe1->0,dZZA1->0,dZg1->-((gs^2 (1+\[Epsilon] Subscript[C, MSB]) (11 N-4 Subscript[n, f] Subscript[T, F]))/(96 \[Pi]^2 \[Epsilon])),dZGG1->gs^2/(12*Pi^2)*(5N/4-Subscript[T, F]*Subscript[n, f])(1/\[Epsilon]+Subscript[C, MSB]),Conjugate[dZfR1[3,1,1]]->dZfR1[3,1,1],dZfR1[3,1,1]->-(gs^2/(1*16*Pi^2))*C2(1/\[Epsilon]+Subscript[C, MSB]),Conjugate[dZfL1[3,1,1]]->dZfL1[3,1,1],dZfL1[3,1,1]->-(gs^2/(1*16*Pi^2))*C2(1/\[Epsilon]+Subscript[C, MSB]),gs->1,Subscript[C, MSB]->Log[4*Pi]-EulerGamma};
sublorind={Lor1->Lorcon1,Lor2->Lorcon2,Lor3->Lorcon3,Lor4->Lorcon4,Lor5->Lorcon5,Lor6->Lorcon6,Lor7->Lorcon7};



cntamps[processname_]:=Module[
{process,instree,ins1loopcnt,amps},

process=Which[
  processname == "Agtoqqb" , "{V[1],V[5]}-> {F[3, {1}], -F[3, {1}]}",
  processname == "Aqtogq" , "{V[1],F[3, {1}]}-> {V[5],F[3, {1}]}",
  processname == "Aqtoqg" , "{V[1],F[3, {1}]}-> {F[3, {1}],V[5]}"
  ]//ToExpression;
instree = InsertFields[CreateTopologies[0, 2 -> 2,ExcludeTopologies->{Tadpoles,WFCorrections}],  process,
		InsertionLevel -> {Classes},Model -> SMQCD, ExcludeParticles -> {S[1], S[2],S[3], V[1],V[2],V[3],V[4],U[5]}];
(*Paint[instree, ColumnsXRows -> {2, 1}, Numbering\[Rule]Simple,SheetHeader->None,ImageSize->{768,256}];
*)
ins1loopcnt = InsertFields[CreateCTTopologies[1, 2 -> 2,ExcludeTopologies->{Tadpoles,WFCorrectionCTs}], process,
		InsertionLevel -> {Classes},Model -> SMQCD, ExcludeParticles -> {S[1], S[2],S[3], V[1],V[2],V[3],V[4],U[5]}];
(*Paint[ins1loopcnt, ColumnsXRows -> {2, 1}, Numbering\[Rule]Simple,SheetHeader->None,ImageSize->{768,256}];*)
amps={instree,ins1loopcnt}
]



CNT[processname_,diags1_,diags2_,fac_]:=
Module[
{ins1,amp1,ins,amp2,s,t,u,p1,p2,q1,q2,VBmom,Int},
amp1=FCFAConvert[CreateFeynAmp[diags1,Truncated -> False],IncomingMomenta->{p1,p2},OutgoingMomenta->{q1,q2},
DropSumOver->True,ChangeDimension->D,UndoChiralSplittings->True,List->False,SMP->True]/.Pair[LorentzIndex[Lor_, D], Momentum[Polarization[p1, I], D]]->FVD[p2,Lor](*Pair[LorentzIndex[Lor, D], Momentum[p1, D]]*);
(*Paint[diags2, ColumnsXRows -> {2, 1}, Numbering\[Rule]Simple,SheetHeader->None,ImageSize->{768,256}];*)
(*Print[amp1];*)
amp2=FCFAConvert[CreateFeynAmp[diags2,Truncated -> False],IncomingMomenta->{p1,p2},OutgoingMomenta->{q1,q2},
DropSumOver->True,ChangeDimension->D,UndoChiralSplittings->True,List->False,SMP->True]*fac/.Pair[LorentzIndex[Lor_, D], Momentum[Polarization[p1, I], D]]->FVD[p2,Lor](*Pair[LorentzIndex[Lor, D], Momentum[p1, D]]*)//.subcnt;
(*Print[Style[amp2,FontColor\[Rule]Red]];*)
ClearScalarProducts;
(*Print[amp];*)
SetMandelstam[s,t,u,p1,p2,-q1,-q2,I*Q,0,0,0];

VBmom=Which[
  processname == "Agtoqqb" , {p1,p2},
  processname == "Aqtogq" , {p1,q1},
  processname == "Aqtoqg" , {p1,q2}
  ];
Int=(amp1*
		(ComplexConjugate[amp2]/.sublorind//FCRenameDummyIndices))//FeynAmpDenominatorSplit//PropagatorDenominatorExplicit//
		SUNSimplify[#,Explicit->True,SUNNToCACF->False]&//FermionSpinSum[#, ExtraFactor -> 1]&//
		Contract//ReplaceAll[#,{DiracTrace->Tr,SUNN->3}]&(*//DoPolarizationSums[#,VBmom[[1]],0]&*)//
		DoPolarizationSums[#,VBmom[[2]],0]&//Factor;
Int=Int//.{u -> U,t->T,s->S,D->4-2\[Epsilon],SMP["e"]->g,SMP["g_s"]->1,SMP["m_u"]->0}(*{u -> -s v,t -> -s (1-v),D\[Rule]4-2\[Epsilon]}*)//Simplify
](*1/3 for color average,  ExtraFactor ->  1/2/2 for spin average. note  the gluon spin average factor 1/2/(1-\[Epsilon]). For virtual photon spin average put 1/2 just to be consistent with other parts. In actual calculation there is no need to average over intermediate particle spin*)




counter[processname_]:=Module[
{process,instree,ins1loopcnt,t1,t2,amps},
amps=cntamps[processname];
renfac=-dZGG1/2-dZfL1[3,1,1];
t1=CNT[processname,amps[[1]],amps[[2]],1]/.{C2->4/3,N->3}/.{S->s,T->t,U->-Q^2-s-t}//Simplify;
t2=CNT[processname,amps[[1]],amps[[1]],renfac]/.{C2->4/3,N->3}/.{S->s,T->t,U->-Q^2-s-t}//Simplify;
(*spin and color average already included from the above*)
cnttot=
Series[avefac[processname]*prefac222virt*Simplify[2*(t1+t2)9/4/.Subscript[t, F]->1/2],{\[Epsilon],0,0}]//Normal(*factor 2 is to include complex conjugate contribution, factor 9/4 is to set u quark charge to 1 instead of 2/3, as is done in real and virtual part*)
]



countersinglepole[xsc_]:=Module[
{xsc1,xsc2,xsc3},
If[xsc=!=0,Expand[Select[Expand[xsc],MatchQ[#,\[Epsilon]^-1__]&]]//Simplify,0]
]
