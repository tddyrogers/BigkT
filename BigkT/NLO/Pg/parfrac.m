(* ::Package:: *)

(* ::Text::RGBColor[1, 0, 0]:: *)
(*Partial fraction algorithm*)


test2V[term_]:=Module[
{den,ttest=0,utest=0,stest=0},
den=Denominator[term];
If[Exponent[den,t2]!=0&&Exponent[den,t3]!=0,ttest=1];
If[Exponent[den,u2]!=0&&Exponent[den,u3]!=0,utest=1];
If[Exponent[den,s12]!=0&&Exponent[den,s13]!=0,stest=1];
{ttest,utest,stest}
]


twoVR[term_]:=Module[{},
term*If[test2V[term][[1]]==0,1,(t2+t3)/(-t1-s-2Q^2)]
*If[test2V[term][[2]]==0,1,(u2+u3)/(-u1-s-Q^2)]
*If[test2V[term][[3]]==0,1,(s12+s13)/(s-s23)]//Expand
]


procpf1[terms_]:=Module[
{list,diff,termstmp,termsnew},
termstmp=terms;
diff=1;
While[diff!=0,
  list=Level[termstmp,1];
Do[
list[[i]]=twoVR[list[[i]]],
(*Print[list[[i]],test1[list[[i]]]],*)
{i,Length[list]}
];
termsnew=Sum[list[[i]],{i,1,Length[list]}];
diff=If[termstmp===termsnew,0,1];
termstmp=termsnew
];
termsnew
]


test3V[term_]:=Module[
{den,t2u2s12test=0,t2u2s13test=0,t2u3s12test=0,t2u3s13test=0,t3u2s12test=0,t3u2s13test=0,t3u3s12test=0,t3u3s13test=0},
den=Denominator[term];
If[Exponent[den,t2]!=0&&Exponent[den,u2]!=0&&Exponent[den,s12]!=0,t2u2s12test=1];
If[Exponent[den,t2]!=0&&Exponent[den,u2]!=0&&Exponent[den,s13]!=0,t2u2s13test=1];
If[Exponent[den,t2]!=0&&Exponent[den,u3]!=0&&Exponent[den,s12]!=0,t2u3s12test=1];
If[Exponent[den,t2]!=0&&Exponent[den,u3]!=0&&Exponent[den,s13]!=0,t2u3s13test=1];
If[Exponent[den,t3]!=0&&Exponent[den,u2]!=0&&Exponent[den,s12]!=0,t3u2s12test=1];
If[Exponent[den,t3]!=0&&Exponent[den,u2]!=0&&Exponent[den,s13]!=0,t3u2s13test=1];
If[Exponent[den,t3]!=0&&Exponent[den,u3]!=0&&Exponent[den,s12]!=0,t3u3s12test=1];
If[Exponent[den,t3]!=0&&Exponent[den,u3]!=0&&Exponent[den,s13]!=0,t3u3s13test=1];
{t2u2s12test,t2u2s13test,t2u3s12test,t2u3s13test,t3u2s12test,t3u2s13test,t3u3s12test,t3u3s13test}
]


threeVR[term_]:=Module[{},
term*If[test3V[term][[1]]==0,1,(t2+u2+s12)/(-s23-Q^2)]
*If[test3V[term][[2]]==0,1,(s13-t2-u2)/(s+Q^2)]
*If[test3V[term][[3]]==0,1,(u3-s12-t2)/(Q^2+t1)]
*If[test3V[term][[4]]==0,1,(u3+s13-t2)/-u1]
*If[test3V[term][[5]]==0,1,(s12+u2-t3)/-u1]
*If[test3V[term][[6]]==0,1,(u2-t3-s13)/(Q^2+t1)]
*If[test3V[term][[7]]==0,1,(s12-t3-u3)/(Q^2+s)]
*If[test3V[term][[8]]==0,1,(s13+t3+u3)/(-Q^2-s23)]//Expand
]


procpf2[terms_]:=Module[
{list,diff,termstmp,termsnew},
termstmp=terms;
diff=1;
While[diff!=0,
  list=Level[termstmp,1];
Do[
list[[i]]=threeVR[list[[i]]],
(*Print[list[[i]],test1[list[[i]]]],*)
{i,Length[list]}
];
termsnew=Sum[list[[i]],{i,1,Length[list]}];
diff=If[termstmp===termsnew,0,1];
termstmp=termsnew
];
termsnew
]


testden[term_]:=Module[
{den,t2test=0,t3test=0,u2test=0,u3test=0,s12test=0,s13test=0},
den=Denominator[term];
If[Exponent[den,t2]!=0,t2test=1];
If[Exponent[den,t3]!=0,t3test=1];
If[Exponent[den,u2]!=0,u2test=1];
If[Exponent[den,u3]!=0,u3test=1];
If[Exponent[den,s12]!=0,s12test=1];
If[Exponent[den,s13]!=0,s13test=1];
{t2test,t3test,u2test,u3test,s12test,s13test}
]


NumR[term_]:=Module[
{den,num,numnew,termnew,eqns},
den=Denominator[term];
num=Numerator[term];
eqns={t1+t2+t3+s+2Q^2==0,u1+u2+u3+s+Q^2==0,s12+s13+s23==s,s+Q^2==s12-t3-u3};
numnew=Which[
(*two dens case testden={t2test,t3test,u2test,u3test,s12test,s13test}*)
testden[term]=={1,0,1,0,0,0},num//.Flatten[Solve[eqns,{t3,u3,s12,s13}]],
testden[term]=={1,0,0,1,0,0},num//.Flatten[Solve[eqns,{t3,u2,s12,s13}]],
testden[term]=={0,1,1,0,0,0},num//.Flatten[Solve[eqns,{t2,u3,s12,s13}]],
testden[term]=={0,1,0,1,0,0},num//.Flatten[Solve[eqns,{t2,u2,s12,s13}]],

testden[term]=={1,0,0,0,1,0},num//.Flatten[Solve[eqns,{t3,s13,u2,u3}]],
testden[term]=={1,0,0,0,0,1},num//.Flatten[Solve[eqns,{t3,s12,u2,u3}]],
testden[term]=={0,1,0,0,1,0},num//.Flatten[Solve[eqns,{t2,s13,u2,u3}]],
testden[term]=={0,1,0,0,0,1},num//.Flatten[Solve[eqns,{t2,s12,u2,u3}]],

testden[term]=={0,0,1,0,1,0},num//.Flatten[Solve[eqns,{u3,s13,t2,t3}]],
testden[term]=={0,0,1,0,0,1},num//.Flatten[Solve[eqns,{u3,s12,t2,t3}]],
testden[term]=={0,0,0,1,1,0},num//.Flatten[Solve[eqns,{u2,s13,t2,t3}]],
testden[term]=={0,0,0,1,0,1},num//.Flatten[Solve[eqns,{u2,s12,t2,t3}]],

(*one den case: convention for one extra ADMV:t2\[Rule]t2u2,t3\[Rule]t3u3,u2\[Rule]u2s12,u3\[Rule]u3s13,s12\[Rule]s12u2,s13\[Rule]s13u3*)
testden[term]=={1,0,0,0,0,0},num//.Flatten[Solve[eqns,{t3,u3,s12,s13}]],
testden[term]=={0,1,0,0,0,0},num//.Flatten[Solve[eqns,{t2,u2,s12,s13}]],

testden[term]=={0,0,1,0,0,0},num//.Flatten[Solve[eqns,{u3,s13,t2,t3}]],
testden[term]=={0,0,0,1,0,0},num//.Flatten[Solve[eqns,{u2,s12,t2,t3}]],

testden[term]=={0,0,0,0,1,0},num//.Flatten[Solve[eqns,{u3,s13,t2,t3}]],
testden[term]=={0,0,0,0,0,1},num//.Flatten[Solve[eqns,{u2,s12,t2,t3}]],

(*no den case: convention for two extra ADMVs: u2s12*)
testden[term]=={0,0,0,0,0,0},num//.Flatten[Solve[eqns,{u3,s13,t2,t3}]]
];
termnew=numnew/den//Expand
]


procpf3[terms_]:=Module[
{list,termsnew},
  list=Level[terms,1];
Do[
list[[i]]=NumR[list[[i]]],
(*Print[list[[i]],test1[list[[i]]]],*)
{i,Length[list]}
];
termsnew=Sum[list[[i]],{i,1,Length[list]}]
]
