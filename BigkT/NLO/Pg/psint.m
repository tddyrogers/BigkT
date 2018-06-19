(* ::Package:: *)

(* ::Text::RGBColor[1, 0, 0]:: *)
(*integrals from Appdx. C of "QCD corrections to heavy-quark production in pp collisions" by Van Neerven at el.*)


(*IntnoabueqABC[i_,j_,AA_,BB_,CC_]:=Which[
i==0&&j==2,2Pi*1/(AA^2-BB^2-CC^2),
i==0&&j==1,Pi/Sqrt[BB^2+CC^2] Log[(AA+Sqrt[BB^2+CC^2])/(AA-Sqrt[BB^2+CC^2])],
i==0&&j==-1,2Pi*AA,
i==0&&j==-2,2Pi*(AA^2+1/3(BB^2+CC^2))
]*)


(*InteqabueqABC[i_,j_,aa_,AA_,BB_,CC_]:=Which[
(*i\[Equal]0&&j\[Equal]2,2Pi*1/(AA^2-BB^2-CC^2),
i\[Equal]0&&j==1,Pi/Sqrt[BB^2+CC^2]Log[(AA+Sqrt[BB^2+CC^2])/(AA-Sqrt[BB^2+CC^2])],
i==0&&j==-1,2Pi*AA,
i\[Equal]0&&j==-2,2Pi*(AA^2+1/3(BB^2+CC^2)),*)
i==1&&j==1,Pi*1/(aa(AA+BB)) (2/(n-4)+1*Log[(AA+BB)^2/(AA^2-BB^2-CC^2)]+1*(n-4)/2 (Log[(AA-Sqrt[BB^2+CC^2])/(AA+BB)]^2-1/2 Log[(AA+Sqrt[BB^2+CC^2])/(AA-Sqrt[BB^2+CC^2])]^2+2PolyLog[2,-((BB+Sqrt[BB^2+CC^2])/(AA-Sqrt[BB^2+CC^2]))]-2PolyLog[2,(BB-Sqrt[BB^2+CC^2])/(AA+BB)]))/.n->4-2\[Epsilon],
i==1&&j==2,Pi*1/(aa (AA+BB)^2) (2/(n-4)+Log[(AA+BB)^2/(AA^2-BB^2-CC^2)]+(2(BB^2+CC^2+AA*BB))/(AA^2-BB^2-CC^2))/.n->4-2\[Epsilon],
i==1&&j==-1,Pi*(AA+BB)/aa (2/(n-4)-(2BB)/(AA+BB))/.n->4-2\[Epsilon],
i==1&&j==-2,Pi*(AA+BB)^2/aa (2/(n-4)+(CC^2-4AA*BB-2BB^2)/(AA+BB)^2)/.n->4-2\[Epsilon],
i==2&&j==1,Pi*1/(aa^2(AA+BB)) ((BB^2+AA*BB+CC^2)/(AA+BB)^2*(2/(n-4)+Log[(AA+BB)^2/(AA^2-BB^2-CC^2)])-(2*CC^2)/(AA+BB)^2-1)/.n->4-2\[Epsilon],
i==2&&j==2,Pi*1/(aa^2(AA+BB)^2) (((3CC^2)/(AA+BB)^2+(2BB)/(AA+BB))*(2/(n-4)+Log[(AA+BB)^2/(AA^2-BB^2-CC^2)])-(8*CC^2)/(AA+BB)^2+(2(BB^2+CC^2))/(AA^2-BB^2-CC^2)-1)/.n->4-2\[Epsilon]
]*)


(* ::Text::RGBColor[1, 0, 0]:: *)
(*integrals computed according to Eq.26 in the draft.*)


IntueqabeqABC[i_,j_,D_,c_,Iep_]:=Module[
{out},
out=Which[
i==2&&j==0,(2 \[Pi])/(-1+D^2)+(D \[Pi] \[Epsilon] Log[(-1+D)/(1+D)])/(-1+D^2)/.\[Epsilon]->-2\[Epsilon],
i==1&&j==0,\[Pi] Log[(1+D)/(-1+D)]+1/12 \[Pi] \[Epsilon] (-2 \[Pi]^2+3(Log[16]-Log[(-1+D)(1+D)^3])*Log[(-1+D)/(1+D)]+12 PolyLog[2,(-1+D)/(1+D)])/.\[Epsilon]->-2\[Epsilon],
i==-1&&j==0,2 D \[Pi]-2 D \[Pi] \[Epsilon]/.\[Epsilon]->-2\[Epsilon],
i==-2&&j==0,-(2/9) (-3-9 D^2) \[Pi]-2/9 (4+9 D^2) \[Pi] \[Epsilon]/.\[Epsilon]->-2\[Epsilon],
i==2&&j==-1,1/(12 (-1+D^2)) \[Pi] (24-24 c D+12 c Log[(-1+D)/(1+D)]-12 c D^2 Log[(-1+D)/(1+D)]+\[Epsilon] (-2 c (-1+D^2) \[Pi]^2+3 Log[(-1+D)/(1+D)] (4 D-4 c Log[2]+c D^2 (-4+Log[16])-c (-1+D^2) Log[(-1+D)(1+D)^3])+12 c (-1+D^2) PolyLog[2,(-1+D)/(1+D)]))/.\[Epsilon]->-2\[Epsilon],
i==1&&j==-1,1/12 \[Pi] (24 c+12 (-1+c D) Log[(-1+D)/(1+D)]+\[Epsilon] (-2 \[Pi]^2+2 c (-12+D \[Pi]^2)+3 (-1+c D) Log[(-1+D)/(1+D)] Log[(-1+D)(1+D)^3/16]-12 (-1+c D) PolyLog[2,(-1+D)/(1+D)]))/.\[Epsilon]->-2\[Epsilon],
i==-1&&j==-1,-(2/9) (-3 c-9 D) \[Pi]-2/9 (4 c+9 D) \[Pi] \[Epsilon]/.\[Epsilon]->-2\[Epsilon],
i==-2&&j==-1,-(2/9) (-3-6 c D-9 D^2) \[Pi]-2/9 (4+8 c D+9 D^2) \[Pi] \[Epsilon]/.\[Epsilon]->-2\[Epsilon],
i==2&&j==-2,\[Pi] (4-4 c D-2 D^2+c^2 (-4+6 D^2)+(-1+D) (1+D) (-D+c (-2+3 c D)) Log[(-1+D)/(1+D)])/(-1+D^2)/.\[Epsilon]->-2\[Epsilon],
i==1&&j==-2,-(1/24) \[Pi] (-96 c-24 D+72 c^2 D+12 (3-4 c D-D^2+c^2 (-1+3 D^2)) Log[(-1+D)/(1+D)])-1*1/24 \[Pi] \[Epsilon] (96 c+36 D-84 c^2 D-6 (1-c^2-D^2+c^2 D^2) Log[(-1+D)/(1+D)]-3(-3+4 c D+D^2+c^2 (1-3 D^2)) *(2/3 \[Pi]^2+Log[(-1+D)(1+D)^3/16] Log[(-1+D)/(1+D)]-4  PolyLog[2,(-1+D)/(1+D)]))/.\[Epsilon]->-2\[Epsilon],
i==-1&&j==-2,-(2/9) (-6 c-12 D) \[Pi]-2/9 (8 c+13 D) \[Pi] \[Epsilon]/.\[Epsilon]->-2\[Epsilon],
i==-2&&j==-2,4/15 (3+c^2+10 c D+10 D^2) \[Pi]-2/225 (123+46 c^2+400 c D+325 D^2) \[Pi] \[Epsilon]/.\[Epsilon]->-2\[Epsilon],
i==1&&j==2,((-2+2 c D) \[Pi])/((c-D)^3 \[Epsilon])+(\[Pi] (2-c^2-2 c D+D^2-(-1+c D) Log[(D^2-1)/(D-c)^2]))/(c-D)^3/.\[Epsilon]->-2\[Epsilon],
i==2&&j==2,-((2 (-3+c^2+2 c D) \[Pi])/((c-D)^4 \[Epsilon]))+1/((c-D)^4 (-1+D^2)) \[Pi] (8+2 c D (-3+D^2)-D^2 (5+D^2)+c^2 (-5+7 D^2)-(-3+c^2+2 c D) (-1+D^2) Log[(-c+D)^2/(-1+D^2)])/.\[Epsilon]->-2\[Epsilon],
i==-1&&j==2,(2 c \[Pi])/\[Epsilon]+(c \[Pi]-D \[Pi])+1/4 (\[Pi]+4 c \[Pi]-3 D \[Pi]) \[Epsilon]/.\[Epsilon]->-2\[Epsilon],
i==-2&&j==2,(2 \[Pi]-6 c^2 \[Pi]+4 c D \[Pi])/\[Epsilon]+(-2 \[Pi]+3 c^2 \[Pi]+2 c D \[Pi]-D^2 \[Pi])+1/4 (8 \[Pi]-2 c \[Pi]-21 c^2 \[Pi]+2 D \[Pi]+8 c D \[Pi]-3 D^2 \[Pi]) \[Epsilon]/.\[Epsilon]->-2\[Epsilon],
i==1&&j==1,-((2 \[Pi])/((c-D) \[Epsilon]))+1*(\[Pi]  Log[(D^2-1)/(D-c)^2])/(c-D)+1/(12 (c-D)) \[Pi] \[Epsilon] (-2 \[Pi]^2+12(Log[1+c]*Log[(D-c)/(D+1)]+Log[D-c]*Log[(D-1)/(D-c)])-3Log[(D-1)/(D+1)]*Log[(D-1)(D+1)^3]-12 PolyLog[2,(-1+c)/(-1+D)]+12 PolyLog[2,(-c+D)/(1+D)])/.\[Epsilon]->-2\[Epsilon],
i==2&&j==1,(2 \[Pi])/((c-D)^2 \[Epsilon])+(\[Pi] (2-2 c D+(-1+D^2) Log[(c-D)^2/(-1+D^2)]))/((c-D)^2 (-1+D^2))+1/(12 (c-D)^2 (-1+D^2)) \[Pi] \[Epsilon] (12 (-1+c+(-1+D) D+(-1+D^2) Log[1+c]) Log[1+D]+3 (-1+D^2) Log[(-1+D)/(1+D)] Log[(-1+D) (1+D)^3]-12 Log[-1+D] (1+c-D (1+D)+(-1+D^2) Log[-c+D])+2 (-1+D^2) (\[Pi]^2+6 Log[-c+D] (-2+Log[(-c+D)/(1+c)]))+12 (-1+D^2) (PolyLog[2,(-1+c)/(-1+D)]-PolyLog[2,(-c+D)/(1+D)]))/.\[Epsilon]->-2\[Epsilon],
i==-2&&j==1,(2 (c-D)^2 \[Pi])/\[Epsilon]+(1-3 c^2+4 c D) \[Pi]+1/2 (-3+7 c^2-8 c D) \[Pi] \[Epsilon]/.\[Epsilon]->-2\[Epsilon],
i==-1&&j==1,(2 (-c+D) \[Pi])/\[Epsilon]+2 c \[Pi]-2 (c \[Pi]) \[Epsilon]/.\[Epsilon]->-2\[Epsilon],
True,Print[{i,j},"not implemented"];0
];(*need (2,1),(1,1) and (-1,2)*)
(*If[Iep==1&&Coefficient[out,\[Epsilon],1]===0,Print["insufficient accuracy for case ",{i,j}],Print[{i,j}]];*)
If[Iep==1,Print[{i,j}," O(\[Epsilon]) term needed"],Print[{i,j}]];
(*If[j!=2,out=0];*)(*check if contribution of j=2 terms cancel*)
(*Print[out];*)
out
]


Intnomass[i_,j_,c_]:=2Pi*Gamma[1-2\[Epsilon]]/Gamma[1-\[Epsilon]] *2^(i+j)*Beta[1-\[Epsilon]+i,1-\[Epsilon]+j]/Gamma[1-\[Epsilon]]*Hypergeometric2F1[-i,-j,1-\[Epsilon],c]


(*for the case i=j=1, use following expansion*)


Intnomassexp[s_,c_]:=(-Pi/\[Epsilon])*(s)^(-1-\[Epsilon]) (1+\[Epsilon]^2PolyLog[2,c])


subtheta2={Sin[\[Theta]2/2]->Sqrt[-((s23 u)/((s-s23) (s23-t)))],Cos[\[Theta]2/2]->Sqrt[-((Q^2 s23+s t)/((s-s23) (s23-t)))]}


subut={sQ2->s+Q^2,su1->s+Q^2+u,p10->(s+t)/(2*Sqrt[s23]),p10z->Sqrt[p10^2+Q^2],p20->(s23-t)/(2*Sqrt[s23]),k10->-((s23-s)/(2*Sqrt[s23])),k20->Sqrt[s23]/2,k30->Sqrt[s23]/2,Cos[\[Theta]1]->(-s^2+s (s23-t)-s23 (2 Q^2+t))/((-s+s23) Sqrt[4 Q^2 s23+(s+t)^2]),Sin[\[Theta]1]->(2Sqrt[s23*u*(Q^2 s23+s t)])/((s-s23) Sqrt[4 Q^2 s23+(s+t)^2]),Sin[\[Alpha]1/2]->Sqrt[-((s23 u)/((s-s23) (s23-t)))],Cos[\[Alpha]1/2]->Sqrt[-((Q^2 s23+s t)/((s-s23) (s23-t)))],Cos[\[Alpha]2]->(-2 Q^2 s23+(s23-t) t-s (s23+t))/( (s23-t) Sqrt[4 Q^2 s23+(s+t)^2]),Sin[\[Alpha]2]->(2Sqrt[s23*u*(Q^2 s23+s t)])/((s23-t) Sqrt[4 Q^2 s23+(s+t)^2]),Sin[\[Theta]2/2]->Sqrt[-((s23 u)/((s-s23) (s23-t)))],Cos[\[Theta]2/2]->Sqrt[-((Q^2 s23+s t)/((s-s23) (s23-t)))],Cos[\[Alpha]3]->(-2 Q^2 s23+(s23-t) t-s (s23+t))/( (s23-t) Sqrt[4 Q^2 s23+(s+t)^2]),Sin[\[Alpha]3]->(2Sqrt[s23*u*(Q^2 s23+s t)])/((s23-t) Sqrt[4 Q^2 s23+(s+t)^2]),Cos[\[Theta]3]->(-s^2+s (s23-t)-s23 (2 Q^2+t))/((-s+s23) Sqrt[4 Q^2 s23+(s+t)^2]),Sin[\[Theta]3]->(2Sqrt[s23*u*(Q^2 s23+s t)])/((s-s23) Sqrt[4 Q^2 s23+(s+t)^2]),t1->t,u1->u,D->4-2*\[Epsilon],u->s23-Q^2-s-t}


psint[term_]:=Module[
{l2,l3,j2,j3,k2,k3,m,n,c,Iep,out1,out2},
l2=Exponent[term,t2];
l3=Exponent[term,t3];
j2=Exponent[term,u2];
j3=Exponent[term,u3];
k2=Exponent[term,s12];
k3=Exponent[term,s13];
Iep=If[Exponent[term,s23]==-1,1,0];
out1=Which[
(*two ADMVs*)
l2!=0&&j2!=0,
m=l2;
n=j2;
term/.t2^m*u2^n->Simplify[(-2*k20*p10z)^m*(-2p20*k20)^n*IntueqabeqABC[-m,-n,(2*p10*k20+Q^2)/(2*k20*p10z),Cos[\[Alpha]3],Iep]//.subut,Assumptions->{s^2+4 Q^2 s23+2 s t+t^2>0,s23>0}],
l2!=0&&j3!=0,
m=l2;
n=j3;
term/.t2^m*u3^n->Simplify[(-2*k20*p10z)^m*(-2p20*k30)^n*IntueqabeqABC[-m,-n,(2*p10*k20+Q^2)/(2*k20*p10z),-Cos[\[Alpha]3],Iep]//.subut,Assumptions->{s^2+4 Q^2 s23+2 s t+t^2>0,s23>0}],
l3!=0&&j2!=0,
m=l3;n=j2;term/.t3^m*u2^n->Simplify[(-2*k30*p10z)^m*(-2p20*k20)^n*IntueqabeqABC[-m,-n,(2*p10*k30+Q^2)/(2*k30*p10z),-Cos[\[Alpha]3],Iep]//.subut,Assumptions->{s^2+4 Q^2 s23+2 s t+t^2>0,s23>0}],
l3!=0&&j3!=0,
m=l3;n=j3;term/.t3^m*u3^n->Simplify[(-2*k30*p10z)^m*(-2p20*k30)^n*IntueqabeqABC[-m,-n,(2*p10*k30+Q^2)/(2*k30*p10z),Cos[\[Alpha]3],Iep]//.subut,Assumptions->{s^2+4 Q^2 s23+2 s t+t^2>0,s23>0}],

l2!=0&&k2!=0,
m=l2;n=k2;term/.t2^m*s12^n->Simplify[(-2*k20*p10z)^m*(2k10*k20)^n*IntueqabeqABC[-m,-n,(2*p10*k20+Q^2)/(2*k20*p10z),Cos[\[Theta]3],Iep]//.subut,Assumptions->{s^2+4 Q^2 s23+2 s t+t^2>0,s23>0}],
l2!=0&&k3!=0,
m=l2;n=k3;term/.t2^m*s13^n->Simplify[(-2*k20*p10z)^m*(2k10*k30)^n*IntueqabeqABC[-m,-n,(2*p10*k20+Q^2)/(2*k20*p10z),-Cos[\[Theta]3],Iep]//.subut,Assumptions->{s^2+4 Q^2 s23+2 s t+t^2>0,s23>0}],
l3!=0&&k2!=0,
m=l3;n=k2;term/.t3^m*s12^n->Simplify[(-2*k30*p10z)^m*(2k10*k20)^n*IntueqabeqABC[-m,-n,(2*p10*k20+Q^2)/(2*k20*p10z),-Cos[\[Theta]3],Iep]//.subut,Assumptions->{s^2+4 Q^2 s23+2 s t+t^2>0,s23>0}],
l3!=0&&k3!=0,
m=l3;n=k3;term/.t3^m*s13^n->Simplify[(-2*k30*p10z)^m*(2k10*k30)^n*IntueqabeqABC[-m,-n,(2*p10*k20+Q^2)/(2*k20*p10z),Cos[\[Theta]3],Iep]//.subut,Assumptions->{s^2+4 Q^2 s23+2 s t+t^2>0,s23>0}],

j2!=0&&k2!=0&&(j2>=-1||k2>=-1),
m=j2;
n=k2;
term/.u2^m*s12^n->Simplify[(-2*k20*p20)^m*(2k10*k20)^n*If[m==-1&&n==-1,Intnomassexp[Sin[\[Theta]2/2]^2,Cos[\[Theta]2/2]^2],Intnomass[m,n,Cos[\[Theta]2/2]^2]]//.subut,Assumptions->{s^2+4 Q^2 s23+2 s t+t^2>0,s23>0}],
j2!=0&&k3!=0&&(j2>=-1||k3>=-1),
m=j2;
n=k3;
term/.u2^m*s13^n->Simplify[(-2*k20*p20)^m*(2k10*k30)^n*If[m==-1&&n==-1,Intnomassexp[Cos[\[Theta]2/2]^2,Sin[\[Theta]2/2]^2],Intnomass[m,n,Sin[\[Theta]2/2]^2]]//.subut,Assumptions->{s^2+4 Q^2 s23+2 s t+t^2>0,s23>0}],
j3!=0&&k2!=0&&(j3>=-1||k2>=-1),
m=j3;
n=k2;
term/.u3^m*s12^n->Simplify[(-2*k30*p20)^m*(2k10*k20)^n*If[m==-1&&n==-1,Intnomassexp[Cos[\[Theta]2/2]^2,Sin[\[Theta]2/2]^2],Intnomass[m,n,Sin[\[Theta]2/2]^2]]//.subut,Assumptions->{s^2+4 Q^2 s23+2 s t+t^2>0,s23>0}],
j3!=0&&k3!=0&&(j3>=-1||k3>=-1),
m=j3;
n=k3;
term/.u3^m*s13^n->Simplify[(-2*k30*p20)^m*(2k10*k30)^n*If[m==-1&&n==-1,Intnomassexp[Sin[\[Theta]2/2]^2,Cos[\[Theta]2/2]^2],Intnomass[m,n,Cos[\[Theta]2/2]^2]]//.subut,Assumptions->{s^2+4 Q^2 s23+2 s t+t^2>0,s23>0}],

(*one ADMV*)
l2!=0,
m=l2;
n=0;
term/.t2^m->Simplify[(-2*k20*p10z)^m*IntueqabeqABC[-m,-n,(2*p10*k20+Q^2)/(2*k20*p10z),Cos[\[Alpha]3],Iep]//.subut,Assumptions->{s^2+4 Q^2 s23+2 s t+t^2>0,s23>0}],
l3!=0,
m=l3;
n=0;
term/.t3^m->Simplify[(-2*k30*p10z)^m*IntueqabeqABC[-m,-n,(2*p10*k20+Q^2)/(2*k20*p10z),Cos[\[Alpha]3],Iep]//.subut,Assumptions->{s^2+4 Q^2 s23+2 s t+t^2>0,s23>0}],
j2!=0,
m=j2;
n=0;
term/.u2^m->Simplify[(-2*k20*p20)^m*Intnomass[m,n,Cos[\[Theta]2/2]^2]//.subut,Assumptions->{s^2+4 Q^2 s23+2 s t+t^2>0,s23>0}],
j3!=0,
m=j3;
n=0;
term/.u3^m->Simplify[(-2*k30*p20)^m*Intnomass[m,n,Cos[\[Theta]2/2]^2]//.subut,Assumptions->{s^2+4 Q^2 s23+2 s t+t^2>0,s23>0}],
k2!=0,
m=0;
n=k2;
term/.s12^n->Simplify[(2k10*k20)^n*Intnomass[m,n,Cos[\[Theta]2/2]^2]//.subut,Assumptions->{s^2+4 Q^2 s23+2 s t+t^2>0,s23>0}],
k3!=0,
m=0;
n=k3;
term/.s13^n->Simplify[(2k10*k30)^n*Intnomass[m,n,Cos[\[Theta]2/2]^2]//.subut,Assumptions->{s^2+4 Q^2 s23+2 s t+t^2>0,s23>0}],

(*no ADMV*)
l2==0&&l3==0&&j2==0&&j3==0&&k2==0&&k3==0,
term*Intnomass[0,0,Cos[\[Theta]2/2]^2],
True,
Print["integral not implemented ",{l2,l3,j2,j3,k2,k3}];Abort[]
]
]


SimpNum[terms_]:=Module[
{list,num,den,i},
list=Level[Expand[terms],1];
Do[
num=Numerator[list[[i]]//.{t1->t,u1->u,u->s23-Q^2-s-t}];
den=Denominator[list[[i]]];
list[[i]]=num/den//Expand,
{i,Length[list]}
];
Sum[list[[i]],{i,1,Length[list]}]//Expand
]


ADMVSimp[terms_]:=Module[
{l2,l3,j2,j3,k2,k3,list,term,termsnew},
  list=Level[terms+u0,1];
Do[
term=list[[i]];
l2=Exponent[term,t2];
l3=Exponent[term,t3];
j2=Exponent[term,u2];
j3=Exponent[term,u3];
k2=Exponent[term,s12];
k3=Exponent[term,s13];
list[[i]]=Which[
l3!=0&&j2!=0,
term/.{t3->t2,u2->u3},
l3!=0&&j3!=0,
term/.{t3->t2,u3->u2},
l3!=0&&k2!=0,
term/.{t3->t2,s12->s13},
l3!=0&&k3!=0,
term/.{t3->t2,s13->s12},
l3!=0,
term/.{t3->t2},
True,
term
],
(*Print[list[[i]],test1[list[[i]]]],*)
{i,Length[list]}
];
termsnew=Sum[list[[i]],{i,1,Length[list]}]/.u0->0
]


CoefSimp[terms_]:=Module[
{list,listpow={},elnew,listterms,termsnew,termsexp},
termsexp=terms//.{t1->t,u1->u,u->s23-Q^2-s-t}//Expand;
list=Level[termsexp,1];
Do[
elnew={Exponent[list[[i]],t2],Exponent[list[[i]],t3],Exponent[list[[i]],u2],Exponent[list[[i]],u3],Exponent[list[[i]],s12],Exponent[list[[i]],s13]};
If[FreeQ[listpow,elnew],AppendTo[listpow,elnew]],
{i,Length[list]}
];
listterms=Table[{0,0},{i,Length[listpow]}];
Do[
Do[
elnew={Exponent[list[[i]],t2],Exponent[list[[i]],t3],Exponent[list[[i]],u2],Exponent[list[[i]],u3],Exponent[list[[i]],s12],Exponent[list[[i]],s13]};
Which[
elnew==listpow[[j]]&&Exponent[list[[i]],s23]!=-1,
listterms[[j]][[1]]=listterms[[j]][[1]]+list[[i]],
elnew==listpow[[j]]&&Exponent[list[[i]],s23]==-1,
listterms[[j]][[2]]=listterms[[j]][[2]]+list[[i]]
],
{i,Length[list]}
];
listterms[[j]][[1]]=listterms[[j]][[1]]//Simplify;
listterms[[j]][[2]]=listterms[[j]][[2]]//Simplify,
{j,Length[listpow]}
];
listterms=Flatten[listterms];
Sum[Factor[listterms[[i]]],{i,1,Length[listterms]}]
]


CoefSimp2[terms_]:=Module[
{list,listpow={},elnew,listterms,termsnew,termsexp},
termsexp=terms//.{t1->t,u1->u,u->s23-Q^2-s-t}//Expand;
list=Level[termsexp,1];
Do[
elnew={Exponent[list[[i]],t2],Exponent[list[[i]],t3],Exponent[list[[i]],u2],Exponent[list[[i]],u3],Exponent[list[[i]],s12],Exponent[list[[i]],s13]};
If[FreeQ[listpow,elnew],AppendTo[listpow,elnew]],
{i,Length[list]}
];
listterms=Table[0,{i,Length[listpow]}];
Do[
Do[
elnew={Exponent[list[[i]],t2],Exponent[list[[i]],t3],Exponent[list[[i]],u2],Exponent[list[[i]],u3],Exponent[list[[i]],s12],Exponent[list[[i]],s13]};
If[elnew==listpow[[j]],
listterms[[j]]=listterms[[j]]+list[[i]]
],
{i,Length[list]}
];
listterms[[j]]=listterms[[j]]//Simplify,
{j,Length[listpow]}
];
Sum[Factor[listterms[[i]]],{i,1,Length[listterms]}]
]


procint[terms_]:=Module[
{list,termsnew},
  list=Level[terms+u0,1];
Do[
list[[i]]=psint[list[[i]]],
(*Print[list[[i]],test1[list[[i]]]],*)
{i,Length[list]}
];
termsnew=Sum[list[[i]],{i,1,Length[list]}]/.u0->0
]


procint2[terms_,prefac_]:=Module[
{list,term,termsnew},
  list=Level[terms+u0,1];
Do[
(*Print[psint[list[[i]]]];*)
term=prefac*psint[list[[i]]];
term=PowerExpand[term,{s23}];
(*term=Collect[term,s23,Simplify]*)
(*Print[term];*)
term=term//.{s23^(-1-\[Epsilon])->-1/\[Epsilon] \[Delta][s23](1-\[Epsilon]*Log[B]+1/2 \[Epsilon]^2 Log[B]^2)+Plus1B[B,s23]-\[Epsilon]*Plus2B[B,s23],s23^(-1-2\[Epsilon])->-1/(2\[Epsilon]) \[Delta][s23](1-2\[Epsilon]*Log[B]+1/2 (2\[Epsilon])^2 Log[B]^2)+Plus1B[B,s23]-2\[Epsilon]*Plus2B[B,s23],
             s23^(-2-\[Epsilon])->1/s23(-1/\[Epsilon] \[Delta][s23](1-\[Epsilon]*Log[B]+1/2 \[Epsilon]^2 Log[B]^2)+Plus1B[B,s23]-\[Epsilon]*Plus2B[B,s23]),s23^(-2-2\[Epsilon])->1/s23(-1/(2\[Epsilon]) \[Delta][s23](1-2\[Epsilon]*Log[B]+1/2 (2\[Epsilon])^2 Log[B]^2)+Plus1B[B,s23]-2\[Epsilon]*Plus2B[B,s23])};
list[[i]]=Normal[Series[term//.subut,{\[Epsilon],0,0}]](*//Simplify[#,Assumptions\[Rule]{s^2+4 Q^2 s23+2 s t+t^2>0,s23>0}]&*),
(*Print[list[[i]]],*)
(*Print[list[[i]],test1[list[[i]]]],*)
{i,Length[list]}
];
termsnew=Sum[list[[i]],{i,1,Length[list]}]/.u0->0
]


(*This function and the above one treats terms with s23^-2 in them*)
procint3[terms_,prefac_]:=Module[
{list,term,termsnew,termsnews23,listnew},
  list=Level[terms+u0,1];
Do[
list[[i]]=psint[list[[i]]],
{i,Length[list]}
];
termsnew=Sum[list[[i]],{i,1,Length[list]}]/.u0->0;
termsnew=Collect[termsnew,s23,Simplify];
If[FullSimplify[Select[termsnew,MatchQ[#,s23^-2__]&]]==0,
termsnew=termsnew/.s23^-2->0,
Print["s23^-2 term appears"];Abort[]
];
listnew=PowerExpand[prefac*Level[termsnew,1]];
termsnew=Sum[listnew[[i]],{i,1,Length[listnew]}];
(*termsnew=PowerExpand[termsnew,{s23}]*)
termsnew=termsnew//.{s23^(-1-\[Epsilon])->-1/\[Epsilon] \[Delta][s23](1-\[Epsilon]*Log[B]+1/2 \[Epsilon]^2 Log[B]^2)+Plus1B[B,s23]-\[Epsilon]*Plus2B[B,s23],s23^(-1-2\[Epsilon])->-1/(2\[Epsilon]) \[Delta][s23](1-2\[Epsilon]*Log[B]+1/2 (2\[Epsilon])^2 Log[B]^2)+Plus1B[B,s23]-2\[Epsilon]*Plus2B[B,s23]};
termsnew=Normal[Series[termsnew//.subut,{\[Epsilon],0,0}]]
]

