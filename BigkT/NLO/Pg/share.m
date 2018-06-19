(* ::Package:: *)

prefac222virt=2*Pi(*\[Mu]^(2\[Epsilon])**)(*((Q^2+s)^2/((-t-s-Q^2)(s*t)))^\[Epsilon]*)/(2Pi)^4;
prefac222subTR=2*Pi/(2Pi)^4;
prefac222subBW=prefac222virt;
(*The total makeup dimension is *\[Mu]^(4\[Epsilon]). For the virtual process, a factor *\[Mu]^(2\[Epsilon]) is included 
in the definition of the loop integral, so only a factor *\[Mu]^(2\[Epsilon]) is included in prefac222. *)
prefac223=2^(2\[Epsilon]-4)*Pi^(\[Epsilon]-2)*s23^-\[Epsilon]*(Gamma[1-\[Epsilon]]/Gamma[1-2\[Epsilon]])*\[Mu]^(2\[Epsilon])(*((Q^2+s)^2/(u*(s*t+Q^2*s23)))^\[Epsilon]*)/(2Pi)^4;(*(Eq.30)*)


avefac[processname_] = Which[
  processname == "Aqtoqg" || processname == "Aqtogq" , 1/2/2/3,
  processname == "Agtoqqb", 1/2/2/(1-\[Epsilon])/8
  ];


collectLog1[expr_]:=Module[
{rule1,rule2,a,b,c,x},
rule2=Log[a_^x_]->x*Log[a];
rule11=Log[a_*b_]->HoldForm[Log[a]+Log[b]];
rule12=Log[a_/b_]->Log[a]-Log[b];
rule13=Log[1/b_]->-Log[b];
(*expr=expr//Refine[#,Assumptions->{v<1,v>0,s>0,\[Mu]>0}]&;*)
(*expr//.rule11*)
ReleaseHold[expr/.Log[c_]->Log[HoldForm[Factor[c]]]]/.-1+v->-HoldForm[(1-v)]//.rule11//.rule12//.rule13//.rule2//ReleaseHold//ReleaseHold//ReleaseHold//ReleaseHold(*//Refine[#,Assumptions->{v<1,v>0,s>0,\[Mu]>0}]&*)
];
