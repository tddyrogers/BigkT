(* ::Package:: *)

setvars = {

p1.p2 -> (s+Q^2)/2, q1.q2 -> s/2,
p2.p1 -> (s+Q^2)/2, q2.q1 -> s/2,
p1.q1 -> -(t+Q^2)/2, p1.q2 -> -(u+Q^2)/2,
q1.p1 -> -(t+Q^2)/2, q2.p1 -> -(u+Q^2)/2,
p2.q1 -> -u/2, p2.q2 -> -t/2,
q1.p2 -> -u/2, q2.p2 -> -t/2,
p1.p1->-Q^2, p2.p2->0, q1.q1->0, q2.q2->0
,p1^2->-Q^2, p2^2->0, q1^2->0, q2^2->0
};
setparvars = {

u -> -Q^2-t-s,
t -> t

};

setvarsinv = {
s ->  2 p1.p2,
t1 -> -2 p1.q1, t2 -> -2 p1.q2,
u1 -> -2 p2.q1, u2 -> -2 p2.q2
};
(*setvars={1->1}*)
(*setvarsinv={1->1}*)





Protected[dim,\[Epsilon],\[Epsilon]feyn];


arguments[_[x__]]:={x};
(*den[a__]:=Module[{args=arguments[a],table={}},table=Table[If[i\[Equal]j,0,args[[i]] args[[j]]],{i,1,Length[args]},{j,1,Length[args]}]/.Times[-1,b_,c_]\[Rule]-b.c/.Times[b_,c_]\[Rule]b.c/.Dot[-1,c_,d_]\[Rule]-c.d;
1/(Plus@@Plus@@@table)//.setvars//.setparvars
]*)
den[a_]:=1/prod[a,a];


IntegralType[nummomentum_,momentumshift_]:=Which[

(momentumshift//Length)==1,tadpole[nummomentum,momentumshift],
(momentumshift//Length)==2,bubble[nummomentum,momentumshift],
(momentumshift//Length)==3,triangle[nummomentum,momentumshift],
(momentumshift//Length)==4,box[nummomentum,momentumshift]

]


expandnummomentum[nummomentum_]:=Module[
{powers=nummomentum[[All,1]],momenta=nummomentum[[All,2]]},
Join@@Table[Table[momenta[[i]],{j,1,powers[[i]]}],{i,1,nummomentum//Length}]
]







(*Evaluation of bubble integrals*)


(* all momenta must be written in terms of external momenta. This functions does not 
work with, for instance, k1, since it sets Power[mom,2]\[Rule]0, assuming on-shellness and massless particles *)
testlightlike[mom_]:=Module[{momsqr=prod[mom,mom]//Simplify},

If[momsqr==0,1,0,0]

]


testspacelike[mom_]:=If[
prod[mom,mom]//Simplify//Reduce[#>=0&&s>0&&t<0&&-Q^2-t-s<0&&Q>0]&//MatchQ[#,False]&,
1,0,0
]


testtimelike[mom_]:=If[
testspacelike[mom]==0&&testlightlike==0,
1,0,0
]


testspacetimelight[mom_]:=If[testlightlike[mom]==1,0,If[testspacelike[mom]==1,-1,1,1],If[testspacelike[mom]==1,-1,1,1]]


powersub={p1^2->-Q^2, p2^2->0, q1^2->0, q2^2->0};
(*returns the product of two momenta in terms of invariants. Works only in the same cases as 'testlightlike' above.*)
prod[a_,b_]:=Module[{prod=a*b},

(prod//Expand)/.powersub(*/.Power[c_,2]\[Rule]0*)/.Times[c_,mom1_,mom2_]->Dot[c,Dot[mom1,mom2]]/.Times[mom1_/;mom1=!=Q^2,mom2_/;mom2=!=Q^2]->Dot[mom1,mom2]/.Dot[c_,Dot[mom1_,mom2_]]->Times[c,Dot[mom1,mom2]]//.setvars//.setparvars

]



bubbleeval[nummomentum_,momentumshift_]:=Module[{num=nummomentum//expandnummomentum},
Which[

(num//Length)==0,B0v[momentumshift],
(num//Length)==1,bubble1[num,momentumshift]

]
]


(*num should have only one element*)
bubble1eval[num_,momentumshift_]:=B1v[momentumshift]*prod[num[[1]],momentumshift[[2]]-momentumshift[[1]]]








B1eval[momentumshift_]:=Module[
(*num should have one elements*)

{lightlike=testlightlike[momentumshift[[2]]],
R1=0,R2=0,
mom1=0,mom2=0,
momarray={},
formfactorarray={},
invertible=1,
output=0},


mom1=momentumshift[[2]]-momentumshift[[1]];

lightlike=testlightlike[mom1];
(*use solution only if mom1 is not lightlike*)

If[lightlike==0,

-(1/2)B0v[momentumshift],

(*B1v[momentumshift]*)0,

B1v[momentumshift]


]

]

(**eval functions that come form g\[Mu]\[Nu] contraction  (  separate case see reference    )*)
(*todo*)
C00eval[momentumshift_]:=Module[


(*num should have one elements*)

{
mom1=0,mom2=0,
output=0
},


mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];


output=prod[mom1,mom1]*C1v[momentumshift];
output+=(prod[mom2,mom2] + 2 prod[mom1,mom2])*C2v[momentumshift];
output+= B0v[{momentumshift[[2]],momentumshift[[3]]}]; 
output *=1/2/(dim-2)(*/.dim\[Rule]4*)(*BW:general dim? leave as it is for now since it has no effect on double pole*)


]

B1eval[{q1,p1}]

BUBBLEFF[term_]:=term/.bubble->bubbleeval/.bubble1->bubble1eval;







getGrahamMatrixTriangle[momentumshift_]:=Module[{momenta={},len=0},
(*extra shift -momentumshift[[1]] is the one from the reduction procedure. indices dont change, only need to shift mom in graham matrix *)
momenta={momentumshift[[2]]-momentumshift[[1]],momentumshift[[3]]-momentumshift[[2]]-0*momentumshift[[1]]}//Simplify;
(*BW: external momenta depends only on adjacent mom shifts*)
len=momenta//Length;

Table[Table[prod[momenta[[i]],momenta[[j]]],{j,1,len}],{i,1,len}]

]



triangleeval[nummomentum_,momentumshift_]:=Module[{num=nummomentum//expandnummomentum},
Which[

(num//Length)==0,C0v[momentumshift],
(num//Length)==1,triangle1[num,momentumshift],
(num//Length)==2,triangle2[num,momentumshift],
(num//Length)==3,triangle3[num,momentumshift]

]

]
(***write form factors***)


triangle1eval[num_,momentumshift_]:=Module[


(*num should have one elements*)

{Graham=getGrahamMatrixTriangle[momentumshift],
R1=0,R2=0,
mom1=0,mom2=0,
momarray={},
formfactorarray={},
invertible=1,
output=0},


mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];


(*use solution only if Graham matrix is invertible*)

invertible=If[Det[Graham]==0,0,1,1];

If[invertible==1,

momarray={
prod[num[[1]],mom1],
prod[num[[1]],mom2]
};

formfactorarray={

C1v[momentumshift],
C2v[momentumshift]

};

output = momarray.formfactorarray,

triangle1[num,momentumshift],(*if Graham matrix is singular, leave unevaluated, for now*)

triangle1[num,momentumshift]    (*if Graham matrix is singular, leave unevaluated, for now*)

]


]


(*Right hand side equations*)

Rc1[momentumshift_]:=Module[
{mom1=0,mom2=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

1/2 B0v[{momentumshift[[1]],momentumshift[[3]]}]
-1/2 B0v[{momentumshift[[2]],momentumshift[[3]]}]
-prod[mom1,mom1]/2 C0v[momentumshift]
]



Rc2[momentumshift_]:=Module[
{mom1=0,mom2=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

1/2 B0v[{momentumshift[[1]],momentumshift[[2]]}]
-1/2 B0v[{momentumshift[[1]],momentumshift[[3]]}]
-(prod[mom2,mom2]+2 prod[mom1,mom2])/2*C0v[momentumshift]
]


(**solutions for Ci**)

code1[i_]:="[momentumshift_]:=Module[\[IndentingNewLine]{\[IndentingNewLine]Graham=getGrahamMatrixTriangle[momentumshift],\[IndentingNewLine]Rarray={Rc1[momentumshift],Rc2[momentumshift]},\[IndentingNewLine]invertible=0,\[IndentingNewLine]sol={}\[IndentingNewLine]},\[IndentingNewLine]\[IndentingNewLine]invertible=If[Det[Graham]\[Equal]0,0,1,1];\[IndentingNewLine]\[IndentingNewLine]If[invertible\[Equal]1,\[IndentingNewLine]\[IndentingNewLine]sol=Inverse[Graham].(Rarray);\[IndentingNewLine]\[IndentingNewLine]sol[["<>ToString[i]<>"]],\[IndentingNewLine]\[IndentingNewLine]C"<>ToString[i]<>"v[momentumshift],(*if Graham matrix is singular, leave unevaluated, for now*)\[IndentingNewLine]\[IndentingNewLine]C"<>ToString[i]<>"v[momentumshift] (*if Graham matrix is singular, leave unevaluated, for now*)\[IndentingNewLine]\[IndentingNewLine]]\[IndentingNewLine]\[IndentingNewLine]\[IndentingNewLine]]"

stringfunct1=Table["C"<>ToString[i]<>"eval"<>code1[i],{i,1,2}];


stringfunct1//ToExpression;



code1[2]
\!\(TraditionalForm\`"\<[momentumshift_]:=Module[\[IndentingNewLine]{\[IndentingNewLine]Graham=getGrahamMatrixTriangle[momentumshift],\[IndentingNewLine]Rarray={Rc1[momentumshift],Rc2[momentumshift]},\[IndentingNewLine]invertible=0,\[IndentingNewLine]sol={}\[IndentingNewLine]},\[IndentingNewLine]\[IndentingNewLine]invertible=If[Det[Graham]\[Equal]0,0,1,1];\[IndentingNewLine]\[IndentingNewLine]If[invertible\[Equal]1,\[IndentingNewLine]\[IndentingNewLine]sol=Inverse[Graham].(Rarray);\[IndentingNewLine]\[IndentingNewLine]sol[[2]],\[IndentingNewLine]\[IndentingNewLine]C2v[momentumshift],(*if Graham matrix is singular, leave unevaluated, for now*)\[IndentingNewLine]\[IndentingNewLine]C2v[momentumshift] (*if Graham matrix is singular, leave unevaluated, for now*)\[IndentingNewLine]\[IndentingNewLine]]\[IndentingNewLine]\[IndentingNewLine]\[IndentingNewLine]]\>"\)

triangle2eval[num_,momentumshift_]:=Module[

{Graham=getGrahamMatrixTriangle[momentumshift],
R11=0,R12=0,
R21=0,R22=0,
mom1=0,mom2=0,
momarray={},
formfactorarray={},
invertible=1,
output=0},

(*num above should be a list with two momenta*)

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];



(*use solution only if Graham matrix is invertible*)

invertible=If[Det[Graham]==0,0,1,1];

If[invertible==1,

momarray={
prod[num[[1]],mom1]*prod[num[[2]],mom1],
prod[num[[1]],mom1]*prod[num[[2]],mom2]
};

formfactorarray={

C11[momentumshift],
C12[momentumshift]

};


output   =momarray.formfactorarray;

momarray={
prod[num[[1]],mom2]*prod[num[[2]],mom1],
prod[num[[1]],mom2]*prod[num[[2]],mom2]
};

formfactorarray={

C12[momentumshift],
C22[momentumshift]

};

output +=momarray.formfactorarray;


output += prod[num[[1]],num[[2]]] *C00[momentumshift]

,

triangle2[num,momentumshift],(*if Graham matrix is singular, leave unevaluated, for now*)

triangle2[num,momentumshift]    (*if Graham matrix is singular, leave unevaluated, for now*)

]


]

(*Right hand side equations*)


Rc11[momentumshift_]:=Module[
{mom1=0,mom2=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

1/2 B1v[{momentumshift[[1]],momentumshift[[3]]}]
+1/2 B0v[{momentumshift[[2]],momentumshift[[3]]}]

-C00[momentumshift]

-prod[mom1,mom1]/2 C1v[momentumshift]
]

Rc12[momentumshift_]:=Module[
{mom1=0,mom2=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];


1/2 B1v[{momentumshift[[1]],momentumshift[[2]]}]
-1/2 B1v[{momentumshift[[1]],momentumshift[[3]]}]
-(prod[mom2,mom2]+2 prod[mom1,mom2])/2*C1v[momentumshift]
]

Rc21[momentumshift_]:=Module[
{mom1=0,mom2=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];


1/2 B1v[{momentumshift[[1]],momentumshift[[3]]}]
-1/2 B1v[{momentumshift[[2]],momentumshift[[3]]}]
-prod[mom1,mom1]/2*C2v[momentumshift](*BW:should be C2*)
]



Rc22[momentumshift_]:=Module[
{mom1=0,mom2=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];


-(1/2)B1v[{momentumshift[[1]],momentumshift[[3]]}]

-C00[momentumshift]

-(prod[mom2,mom2]+2 prod[mom1,mom2])/2 C2v[momentumshift](*BW:should be C2*)
]




(*solutions for Cij*)

code2[i_,j_]:="[momentumshift_]:=Module[\[IndentingNewLine]{\[IndentingNewLine]Graham=getGrahamMatrixTriangle[momentumshift],\[IndentingNewLine]Rarray={\[IndentingNewLine]{Rc11[momentumshift],Rc12[momentumshift]},\[IndentingNewLine]{Rc21[momentumshift],Rc22[momentumshift]}\[IndentingNewLine]},\[IndentingNewLine]invertible=0,\[IndentingNewLine]sol={}\[IndentingNewLine]},\[IndentingNewLine]\[IndentingNewLine]invertible=If[Det[Graham]\[Equal]0,0,1,1];\[IndentingNewLine]\[IndentingNewLine]If[invertible\[Equal]1,\[IndentingNewLine]\[IndentingNewLine]sol=Inverse[Graham].(Rarray[["<>ToString[i]<>"]]);\[IndentingNewLine]\[IndentingNewLine]sol[["<>ToString[j]<>"]],\[IndentingNewLine]\[IndentingNewLine]C"<>ToString[i]<>ToString[j]<>"[momentumshift],(*if Graham matrix is singular, leave unevaluated, for now*)\[IndentingNewLine]\[IndentingNewLine]C"<>ToString[i]<>ToString[j]<>"[momentumshift] (*if Graham matrix is singular, leave unevaluated, for now*)\[IndentingNewLine]\[IndentingNewLine]]\[IndentingNewLine]]"

stringfunct2=Table["C"<>ToString[i]<>ToString[j]<>"eval"<>code2[i,j],{i,1,2},{j,1,2}];

stringfunct2//ToExpression;


(*replacements for form factors coming from g\[Mu]\[Nu]*)

(*todo    eq A.61 *)
D00eval[momentumshift_]:=Module[


(*num should have one elements*)

{
mom1=0,mom2=0,mom3=0,
output=0
},


mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];


output   =C0v[{momentumshift[[2]],momentumshift[[3]],momentumshift[[4]]}];
output+= prod[mom1,mom1]*D1v[momentumshift];
output+= (prod[mom2,mom2]+2*prod[mom1,mom2])*D2v[momentumshift];
output+= (prod[mom3,mom3]+2*prod[mom3,mom1]+2*prod[mom3,mom2])*D3v[momentumshift];

output *=  1/2/(dim-3)(*/.dim\[Rule]4*)(*BW:keep dim general? leave as it is for now since it has not effect on double pole*)



]





TRIANGLEFF[term_]:=term/.triangle->triangleeval/.triangle2->triangle2eval/.triangle1->triangle1eval;

C11eval[{a,b},{0,b,c}]
\!\(TraditionalForm\`C11eval({a, b}, {0, b, c})\)



(*Evaluation of box integrals*)
getGrahamMatrixBox[momentumshift_]:=Module[{momenta={},len=0},
(*extra shift -momentumshift[[1]] is the one from the reduction procedure. indices dont change, only need to shift mom in graham matrix *)
momenta={momentumshift[[2]]-momentumshift[[1]],momentumshift[[3]]-momentumshift[[2]]-0*momentumshift[[1]],momentumshift[[4]]-momentumshift[[3]]-0*momentumshift[[1]]};
(*BW: external momenta depends only on adjacent mom shifts*)
len=momenta//Length;

Table[Table[prod[momenta[[i]],momenta[[j]]],{j,1,len}],{i,1,len}]

]


boxeval[nummomentum_,momentumshift_]:=Module[{num=nummomentum//expandnummomentum},
Which[

(num//Length)==0,D0v[momentumshift],
(num//Length)==1,box1[num,momentumshift],
(num//Length)==2,box2[num,momentumshift],
(num//Length)==3,box3[num,momentumshift],
(num//Length)==4,box4[num,momentumshift]

]
]



box1eval[num_,momentumshift_]:=Module[


(*num should have one elements*)

{Graham=getGrahamMatrixBox[momentumshift],
R1=0,R2=0,R3=0,
mom1=0,mom2=0,mom3=0,
momarray={},formfactorarray={},
invertible=1,
output=0},


mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];



(*use solution only if Graham matrix is invertible*)

invertible=If[Det[Graham]==0,0,1,1];

If[invertible==1,

momarray={
prod[num[[1]],mom1],
prod[num[[1]],mom2],
prod[num[[1]],mom3]
};

formfactorarray={
D1v[momentumshift],
D2v[momentumshift],
D3v[momentumshift]
};

output=momarray.formfactorarray,

box1[num,momentumshift],(*if Graham matrix is singular, leave unevaluated, for now*)

box1[num,momentumshift]    (*if Graham matrix is singular, leave unevaluated, for now*)

]


]


(**right hand side of equations*)


Rd1[momentumshift_]:=Module[
{mom1=0,mom2=0,mom3=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];


1/2 C0v[{momentumshift[[1]],momentumshift[[3]],momentumshift[[4]]}]
-1/2 C0v[{momentumshift[[2]],momentumshift[[3]],momentumshift[[4]]}]
-prod[mom1,mom1]/2 D0v[momentumshift]

]



Rd2[momentumshift_]:=Module[
{mom1=0,mom2=0,mom3=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];


1/2 C0v[{momentumshift[[1]],momentumshift[[2]],momentumshift[[4]]}]
-1/2 C0v[{momentumshift[[1]],momentumshift[[3]],momentumshift[[4]]}]
-(prod[mom2,mom2]+2 prod[mom1,mom2])/2 D0v[momentumshift]

]



Rd3[momentumshift_]:=Module[
{mom1=0,mom2=0,mom3=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];


1/2 C0v[{momentumshift[[1]],momentumshift[[2]],momentumshift[[3]]}]
-1/2 C0v[{momentumshift[[1]],momentumshift[[2]],momentumshift[[4]]}]
-(prod[mom3,mom3]+2 prod[mom2,mom3]+2 prod[mom1,mom3])/2 D0v[momentumshift]

]


(**solutions for Di**)

code3[i_]:="[momentumshift_]:=Module[\[IndentingNewLine]{\[IndentingNewLine]Graham=getGrahamMatrixBox[momentumshift],\[IndentingNewLine]Rarray={Rd1[momentumshift],Rd2[momentumshift],Rd3[momentumshift]},\[IndentingNewLine]invertible=0,\[IndentingNewLine]sol={}\[IndentingNewLine]},\[IndentingNewLine]\[IndentingNewLine]invertible=If[Det[Graham]\[Equal]0,0,1,1];\[IndentingNewLine]\[IndentingNewLine]If[invertible\[Equal]1,\[IndentingNewLine]\[IndentingNewLine]sol=Inverse[Graham].(Rarray);\[IndentingNewLine]\[IndentingNewLine]sol[["<>ToString[i]<>"]],\[IndentingNewLine]\[IndentingNewLine]D"<>ToString[i]<>"v[momentumshift],(*if Graham matrix is singular, leave unevaluated, for now*)\[IndentingNewLine]\[IndentingNewLine]D"<>ToString[i]<>"v[momentumshift] (*if Graham matrix is singular, leave unevaluated, for now*)\[IndentingNewLine]\[IndentingNewLine]]\[IndentingNewLine]\[IndentingNewLine]\[IndentingNewLine]]"

stringfunct3=Table["D"<>ToString[i]<>"eval"<>code3[i],{i,1,3}];


stringfunct3//ToExpression;














box2eval[num_,momentumshift_]:=Module[

{Graham=getGrahamMatrixBox[momentumshift],
R11=0,R12=0,R13=0,
R21=0,R22=0,R23=0,
R31=0,R32=0,R33=0,
mom={0,0,0},
mom1=0,mom2=0,mom3=0,
momarray={},formfactorarray={},
invertible=1,
output=0},

(*num above should be a list with two momenta*)

mom[[1]]=momentumshift[[2]]-momentumshift[[1]];

mom[[2]]=momentumshift[[3]]-momentumshift[[2]];

mom[[3]]=momentumshift[[4]]-momentumshift[[3]];

(*mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];*)
(*use solution only if Graham matrix is invertible*)

invertible=If[Det[Graham]==0,0,1,1];

If[invertible==1,
output=0;
(*BW: Now loop over i,j*)
Do[
output+=prod[num[[1]],mom[[i]]]*
prod[num[[2]],mom[[j]]]*
ToExpression["D"<>ToString[i]<>ToString[j]<>"["<>ToString[momentumshift]<>"]"],
{i,1,3},{j,1,3}];

(*momarray={
prod[num[[1]],mom1]*prod[num[[2]],mom1],
prod[num[[1]],mom1]*prod[num[[2]],mom2],
prod[num[[1]],mom1]*prod[num[[2]],mom3]
};



formfactorarray={
D11[momentumshift],
D12[momentumshift],
D13[momentumshift]
};

output =momarray.formfactorarray;


momarray={
prod[num[[1]],mom1]*prod[num[[2]],mom2],
prod[num[[1]],mom2]*prod[num[[2]],mom2],
prod[num[[1]],mom2]*prod[num[[2]],mom3]
};


formfactorarray={
D12[momentumshift],
D22[momentumshift],
D23[momentumshift]
};


output+=momarray.formfactorarray;


momarray={
prod[num[[1]],mom1]*prod[num[[2]],mom3],
prod[num[[1]],mom2]*prod[num[[2]],mom3],
prod[num[[1]],mom3]*prod[num[[2]],mom3]
};


formfactorarray={
D13[momentumshift],
D23[momentumshift],
D33[momentumshift]
};


output+=momarray.formfactorarray;*)


output += prod[num[[1]],num[[2]]]*D00[momentumshift]


,

box2[num,momentumshift],(*if Graham matrix is singular, leave unevaluated, for now*)

box2[num,momentumshift]    (*if Graham matrix is singular, leave unevaluated, for now*)

]


]

(*right hand side if eqs*)
Rd11[momentumshift_]:=Module[
{mom1=0,mom2=0,mom3=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];


1/2 C0v[{momentumshift[[2]],momentumshift[[3]],momentumshift[[4]]}]
+1/2 C1v[{momentumshift[[1]],momentumshift[[3]],momentumshift[[4]]}]
-prod[mom1,mom1]/2 D1v[momentumshift]
-D00[momentumshift]

]

Rd12[momentumshift_]:=Module[
{mom1=0,mom2=0,mom3=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];

-(1/2)C1v[{momentumshift[[1]],momentumshift[[3]],momentumshift[[4]]}]
+1/2 C1v[{momentumshift[[1]],momentumshift[[2]],momentumshift[[4]]}]
-(prod[mom2,mom2]+2 prod[mom1,mom2])/2 D1v[momentumshift]

]


Rd13[momentumshift_]:=Module[
{mom1=0,mom2=0,mom3=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];

+(1/2)C1v[{momentumshift[[1]],momentumshift[[2]],momentumshift[[3]]}]
-1/2 C1v[{momentumshift[[1]],momentumshift[[2]],momentumshift[[4]]}]
-(prod[mom3,mom3]+2 prod[mom1,mom3]+2 prod[mom2,mom3])/2 D1v[momentumshift]

]


(********)

Rd21[momentumshift_]:=Module[
{mom1=0,mom2=0,mom3=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];



-(1/2)C1v[{momentumshift[[2]],momentumshift[[3]],momentumshift[[4]]}]
+1/2 C1v[{momentumshift[[1]],momentumshift[[3]],momentumshift[[4]]}]
-prod[mom1,mom1]/2 D2v[momentumshift]

]

Rd22[momentumshift_]:=Module[
{mom1=0,mom2=0,mom3=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];


-(1/2)C1v[{momentumshift[[1]],momentumshift[[3]],momentumshift[[4]]}]
+1/2 C2v[{momentumshift[[1]],momentumshift[[2]],momentumshift[[4]]}]
-(prod[mom2,mom2]+2 prod[mom2,mom1])/2 D2v[momentumshift]
-D00[momentumshift]

]


Rd23[momentumshift_]:=Module[
{mom1=0,mom2=0,mom3=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];

+(1/2)C2v[{momentumshift[[1]],momentumshift[[2]],momentumshift[[3]]}]
-1/2 C2v[{momentumshift[[1]],momentumshift[[2]],momentumshift[[4]]}]
-(prod[mom3,mom3]+2 prod[mom1,mom3]+2 prod[mom2,mom3])/2 D2v[momentumshift]
]


(********)




Rd31[momentumshift_]:=Module[
{mom1=0,mom2=0,mom3=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];



-(1/2)C2v[{momentumshift[[2]],momentumshift[[3]],momentumshift[[4]]}]
+1/2 C2v[{momentumshift[[1]],momentumshift[[3]],momentumshift[[4]]}]
-prod[mom1,mom1]/2 D3v[momentumshift]

]

Rd32[momentumshift_]:=Module[
{mom1=0,mom2=0,mom3=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];


-(1/2)C2v[{momentumshift[[1]],momentumshift[[3]],momentumshift[[4]]}]
+1/2 C2v[{momentumshift[[1]],momentumshift[[2]],momentumshift[[4]]}]
-(prod[mom2,mom2]+2 prod[mom2,mom1])/2 D3v[momentumshift]
]


Rd33[momentumshift_]:=Module[
{mom1=0,mom2=0,mom3=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];

-(1/2)C2v[{momentumshift[[1]],momentumshift[[2]],momentumshift[[4]]}]

-(prod[mom3,mom3]+2 prod[mom1,mom3]+2 prod[mom2,mom3])/2 D3v[momentumshift]

-D00[momentumshift]
]


(**solutions for Di**)

code4[i_,j_]:="[momentumshift_]:=Module[\[IndentingNewLine]{\[IndentingNewLine]Graham=getGrahamMatrixBox[momentumshift],\[IndentingNewLine]Rarray={\[IndentingNewLine]{Rd11[momentumshift],Rd12[momentumshift],Rd13[momentumshift]},\[IndentingNewLine]{Rd21[momentumshift],Rd22[momentumshift],Rd23[momentumshift]},\[IndentingNewLine]{Rd31[momentumshift],Rd32[momentumshift],Rd33[momentumshift]}\[IndentingNewLine]},\[IndentingNewLine]invertible=0,\[IndentingNewLine]sol={}\[IndentingNewLine]},\[IndentingNewLine]\[IndentingNewLine]invertible=If[Det[Graham]\[Equal]0,0,1,1];\[IndentingNewLine]\[IndentingNewLine]If[invertible\[Equal]1,\[IndentingNewLine]\[IndentingNewLine]sol=Inverse[Graham].(Rarray[["<>ToString[i]<>"]]);\[IndentingNewLine]\[IndentingNewLine]sol[["<>ToString[j]<>"]],\[IndentingNewLine]\[IndentingNewLine]D"<>ToString[i]<>""<>ToString[j]<>"[momentumshift],(*if Graham matrix is singular, leave unevaluated, for now*)\[IndentingNewLine]\[IndentingNewLine]D"<>ToString[i]<>""<>ToString[j]<>"[momentumshift] (*if Graham matrix is singular, leave unevaluated, for now*)\[IndentingNewLine]\[IndentingNewLine]]\[IndentingNewLine]\[IndentingNewLine]\[IndentingNewLine]]"

stringfunct4=Table["D"<>ToString[i]<>ToString[j]<>"eval"<>code4[i,j],{i,1,3},{j,1,3}];


stringfunct4//ToExpression;
(*BW:symmetrize Cij*)
Do[
"D"<>ToString[j]<>ToString[i]<>"eval=D"<>ToString[i]<>ToString[j]<>"eval"//ToExpression,
{i,1,3},{j,1,3}];



box3eval[num_,momentumshift_]:=Module[

{Graham=getGrahamMatrixBox[momentumshift]//Simplify,
R001=0,R002=0,R003=0,
R121=0,R122=0,R123=0,
R131=0,R132=0,R133=0,
R231=0,R232=0,R233=0,
R111=0,R112=0,R113=0,
R221=0,R222=0,R223=0,
R331=0,R332=0,R333=0,
mom={0,0,0},
mom1=0,mom2=0,mom3=0,
momarray={},formfactorarray={},
invertible=1,
output=0},

(*num above should be a list with three momenta*)

mom[[1]]=momentumshift[[2]]-momentumshift[[1]];

mom[[2]]=momentumshift[[3]]-momentumshift[[2]];

mom[[3]]=momentumshift[[4]]-momentumshift[[3]];

(*mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];*)

(*use solution only if Graham matrix is invertible*)

invertible=If[Det[Graham]==0,0,1,1];

If[invertible==1,
(*BW: in the following, one should not implement A.25 according to table A.2, since table A.2 does not include all terms in A.25. One should instead loop over the indices in A25 directory.*)

(******R00i*******)

momarray={

     prod[num[[1]],num[[2]]]*prod[num[[3]],mom[[1]]]
+prod[num[[1]],num[[3]]]*prod[num[[2]],mom[[1]]]
+prod[num[[2]],num[[3]]]*prod[num[[1]],mom[[1]]],
   
      prod[num[[1]],num[[2]]]*prod[num[[3]],mom[[2]]]
+prod[num[[1]],num[[3]]]*prod[num[[2]],mom[[2]]]
+prod[num[[2]],num[[3]]]*prod[num[[1]],mom[[2]]],

     prod[num[[1]],num[[2]]]*prod[num[[3]],mom[[3]]]
+prod[num[[1]],num[[3]]]*prod[num[[2]],mom[[3]]]
+prod[num[[2]],num[[3]]]*prod[num[[1]],mom[[3]]]

(*prod[num[[1]],num[[2]]]*prod[num[[3]],mom1]
+prod[num[[1]],num[[3]]]*prod[num[[2]],mom1]
+prod[num[[2]],num[[3]]]*prod[num[[1]],mom1],
   
      prod[num[[1]],num[[2]]]*prod[num[[3]],mom2]
+prod[num[[1]],num[[3]]]*prod[num[[2]],mom2]
+prod[num[[2]],num[[3]]]*prod[num[[1]],mom2],

     prod[num[[1]],num[[2]]]*prod[num[[3]],mom3]
+prod[num[[1]],num[[3]]]*prod[num[[2]],mom3]
+prod[num[[2]],num[[3]]]*prod[num[[1]],mom3]*)
};(*BW: Need to multiply by 2 since mu nu can be switched in A25 ? No, the origninal P-V paper does not have it. So no change is made here*)

formfactorarray={
D001[momentumshift],
D002[momentumshift],
D003[momentumshift]
};

output   =momarray.formfactorarray;
(*BW: now loop over i,j,k of 2nd term of A25*)
Do[
output+=prod[num[[1]],mom[[i]]]*
prod[num[[2]],mom[[j]]]*
prod[num[[3]],mom[[k]]]*
ToExpression["D"<>ToString[i]<>ToString[j]<>ToString[k]<>"["<>ToString[momentumshift]<>"]"],
{i,1,3},{j,1,3},{k,1,3}];
output

(*(******R12i*******)

momarray={
prod[num[[1]],mom1]*prod[num[[2]],mom1]*prod[num[[3]],mom2],
prod[num[[1]],mom1]*prod[num[[2]],mom2]*prod[num[[3]],mom2],
prod[num[[1]],mom1]*prod[num[[2]],mom2]*prod[num[[3]],mom3]
};

formfactorarray={
D112[momentumshift],
D122[momentumshift],
D123[momentumshift]
};

output+=momarray.formfactorarray;

(*******R13i******)

momarray={
prod[num[[1]],mom1]*prod[num[[2]],mom1]*prod[num[[3]],mom3],
prod[num[[1]],mom1]*prod[num[[2]],mom2]*prod[num[[3]],mom3],
prod[num[[1]],mom1]*prod[num[[2]],mom3]*prod[num[[3]],mom3]
};


formfactorarray={
D113[momentumshift],
D123[momentumshift],
D133[momentumshift]
};

output+=momarray.formfactorarray;


(*******R23i******)

momarray={
prod[num[[1]],mom1]*prod[num[[2]],mom2]*prod[num[[3]],mom3],
prod[num[[1]],mom2]*prod[num[[2]],mom2]*prod[num[[3]],mom3],
prod[num[[1]],mom2]*prod[num[[2]],mom3]*prod[num[[3]],mom3]
};


formfactorarray={
D123[momentumshift],
D223[momentumshift],
D233[momentumshift]
};

output+=momarray.formfactorarray;


(*******R11i******)

momarray={
prod[num[[1]],mom1]*prod[num[[2]],mom1]*prod[num[[3]],mom1],
prod[num[[1]],mom1]*prod[num[[2]],mom1]*prod[num[[3]],mom2],
prod[num[[1]],mom1]*prod[num[[2]],mom1]*prod[num[[3]],mom3]
};


formfactorarray={
D111[momentumshift],
D112[momentumshift],
D113[momentumshift]
};

output+=momarray.formfactorarray;

(*******R22i******)

momarray={
prod[num[[1]],mom1]*prod[num[[2]],mom2]*prod[num[[3]],mom2],
prod[num[[1]],mom2]*prod[num[[2]],mom2]*prod[num[[3]],mom2],
prod[num[[1]],mom2]*prod[num[[2]],mom2]*prod[num[[3]],mom3]
};


formfactorarray={
D122[momentumshift],
D222[momentumshift],
D223[momentumshift]
};

output+=momarray.formfactorarray;

(*******R33i******)

momarray={
prod[num[[1]],mom1]*prod[num[[2]],mom3]*prod[num[[3]],mom3],
prod[num[[1]],mom2]*prod[num[[2]],mom3]*prod[num[[3]],mom3],
prod[num[[1]],mom3]*prod[num[[2]],mom3]*prod[num[[3]],mom3]
};


formfactorarray={
D133[momentumshift],
D233[momentumshift],
D333[momentumshift]
};

output+=momarray.formfactorarray*)

,
Print["Singular Graham"];
symfac3*box3[num,momentumshift],(*if Graham matrix is singular, leave unevaluated, for now*)
Print["Singular Graham"];
symfac3*box3[num,momentumshift]    (*if Graham matrix is singular, leave unevaluated, for now*)

]

]

(*right hand of eqs**)

Rd001[momentumshift_]:=Module[
{mom1=0,mom2=0,mom3=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];



-(1/2)C00[{momentumshift[[2]],momentumshift[[3]],momentumshift[[4]]}]
+1/2 C00[{momentumshift[[1]],momentumshift[[3]],momentumshift[[4]]}]
-prod[mom1,mom1]/2 D00[momentumshift]



]

Rd002[momentumshift_]:=Module[
{mom1=0,mom2=0,mom3=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];



-(1/2)C00[{momentumshift[[1]],momentumshift[[3]],momentumshift[[4]]}]
+1/2 C00[{momentumshift[[1]],momentumshift[[2]],momentumshift[[4]]}]
-(prod[mom2,mom2] + 2 prod[mom2,mom1])/2 D00[momentumshift]



]


Rd003[momentumshift_]:=Module[
{mom1=0,mom2=0,mom3=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];



-(1/2)C00[{momentumshift[[1]],momentumshift[[2]],momentumshift[[4]]}]
+1/2 C00[{momentumshift[[1]],momentumshift[[2]],momentumshift[[3]]}]
-1/2 (prod[mom3,mom3] + 2 prod[mom3,mom1]+ 2 prod[mom3,mom2])D00[momentumshift]


]



(********)


Rd121[momentumshift_]:=Module[
{mom1=0,mom2=0,mom3=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];



+(1/2)C1v[{momentumshift[[2]],momentumshift[[3]],momentumshift[[4]]}]
+1/2 C11[{momentumshift[[1]],momentumshift[[3]],momentumshift[[4]]}]
-prod[mom1,mom1]/2 D12[momentumshift]
-D002[momentumshift]


]


Rd122[momentumshift_]:=Module[
{mom1=0,mom2=0,mom3=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];




+(1/2)C12[{momentumshift[[1]],momentumshift[[2]],momentumshift[[4]]}]
-1/2 C11[{momentumshift[[1]],momentumshift[[3]],momentumshift[[4]]}]
-(prod[mom2,mom2]+2 prod[mom1,mom2])/2 D12[momentumshift]
-D001[momentumshift]


]





Rd123[momentumshift_]:=Module[
{mom1=0,mom2=0,mom3=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];





1/2 C12[{momentumshift[[1]],momentumshift[[2]],momentumshift[[3]]}]
-1/2 C12[{momentumshift[[1]],momentumshift[[2]],momentumshift[[4]]}]
-1/2 (prod[mom3,mom3] + 2 prod[mom3,mom1]+ 2 prod[mom3,mom2])D12[momentumshift]


]




(***************)


Rd131[momentumshift_]:=Module[
{mom1=0,mom2=0,mom3=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];




   1/2 C2v[{momentumshift[[2]],momentumshift[[3]],momentumshift[[4]]}]
+1/2 C12[{momentumshift[[1]],momentumshift[[3]],momentumshift[[4]]}]
-prod[mom1,mom1]/2 D13[momentumshift]
- D003[momentumshift]



]




Rd132[momentumshift_]:=Module[
{mom1=0,mom2=0,mom3=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];





-(1/2)C12[{momentumshift[[1]],momentumshift[[3]],momentumshift[[4]]}]
+1/2 C12[{momentumshift[[1]],momentumshift[[2]],momentumshift[[4]]}]
-(prod[mom2,mom2]+2 prod[mom2,mom1])/2 D13[momentumshift]


]

Rd133[momentumshift_]:=Module[
{mom1=0,mom2=0,mom3=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];



-(1/2)C12[{momentumshift[[1]],momentumshift[[2]],momentumshift[[4]]}]
-1/2 (prod[mom3,mom3]+2 prod[mom1,mom3]+2 prod[mom2,mom3])D13[momentumshift]
-D001[momentumshift]



]






(********)


Rd231[momentumshift_]:=Module[
{mom1=0,mom2=0,mom3=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];



-(1/2)C12[{momentumshift[[2]],momentumshift[[3]],momentumshift[[4]]}]
+1/2 C12[{momentumshift[[1]],momentumshift[[3]],momentumshift[[4]]}]
-prod[mom1,mom1]/2 D23[momentumshift]



]


Rd232[momentumshift_]:=Module[
{mom1=0,mom2=0,mom3=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];




1/2 C22[{momentumshift[[1]],momentumshift[[2]],momentumshift[[4]]}]
-1/2 C12[{momentumshift[[1]],momentumshift[[3]],momentumshift[[4]]}]
-(prod[mom2,mom2]+2 prod[mom2,mom1])/2 D23[momentumshift]

-D003[momentumshift]

]



Rd233[momentumshift_]:=Module[
{mom1=0,mom2=0,mom3=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];





-(1/2)C22[{momentumshift[[1]],momentumshift[[2]],momentumshift[[4]]}]

-1/2 (prod[mom3,mom3]+2 prod[mom1,mom3]+2 prod[mom2,mom3])D23[momentumshift]

-D002[momentumshift]

]





(***************)


Rd111[momentumshift_]:=Module[
{mom1=0,mom2=0,mom3=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];





-(1/2)C0v[{momentumshift[[2]],momentumshift[[3]],momentumshift[[4]]}]
+1/2 C11[{momentumshift[[1]],momentumshift[[3]],momentumshift[[4]]}]

-prod[mom1,mom1]/2 D11[momentumshift]

-2 D001[momentumshift]
]



Rd112[momentumshift_]:=Module[
{mom1=0,mom2=0,mom3=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];






-(1/2)C11[{momentumshift[[1]],momentumshift[[3]],momentumshift[[4]]}]
+1/2 C11[{momentumshift[[1]],momentumshift[[2]],momentumshift[[4]]}]

-(prod[mom2,mom2]+2 prod[mom1,mom2])/2 D11[momentumshift]

]



Rd113[momentumshift_]:=Module[
{mom1=0,mom2=0,mom3=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];





  1/2 C11[{momentumshift[[1]],momentumshift[[2]],momentumshift[[3]]}]
-1/2 C11[{momentumshift[[1]],momentumshift[[2]],momentumshift[[4]]}]

-1/2 (prod[mom3,mom3]+2 prod[mom1,mom3]+2 prod[mom2,mom3])D11[momentumshift]



]






(***********)


Rd221[momentumshift_]:=Module[
{mom1=0,mom2=0,mom3=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];



-(1/2)C11[{momentumshift[[2]],momentumshift[[3]],momentumshift[[4]]}]
+1/2 C11[{momentumshift[[1]],momentumshift[[3]],momentumshift[[4]]}]

-prod[mom1,mom1]/2 D22[momentumshift]


]


Rd222[momentumshift_]:=Module[
{mom1=0,mom2=0,mom3=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];



-(1/2)C11[{momentumshift[[1]],momentumshift[[3]],momentumshift[[4]]}]
+1/2 C22[{momentumshift[[1]],momentumshift[[2]],momentumshift[[4]]}]

-(prod[mom2,mom2]+2 prod[mom1,mom2])/2 D22[momentumshift]
-2 D002[momentumshift]



]


Rd223[momentumshift_]:=Module[
{mom1=0,mom2=0,mom3=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];




  1/2 C22[{momentumshift[[1]],momentumshift[[2]],momentumshift[[3]]}]
-1/2 C22[{momentumshift[[1]],momentumshift[[2]],momentumshift[[4]]}]

-1/2 (prod[mom3,mom3]+2 prod[mom1,mom3] +2 prod[mom2,mom3])D22[momentumshift]


]






(***********)



Rd331[momentumshift_]:=Module[
{mom1=0,mom2=0,mom3=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];



-(1/2)C22[{momentumshift[[2]],momentumshift[[3]],momentumshift[[4]]}]
+1/2 C22[{momentumshift[[1]],momentumshift[[3]],momentumshift[[4]]}]

-prod[mom1,mom1]/2 D33[momentumshift]





]

Rd332[momentumshift_]:=Module[
{mom1=0,mom2=0,mom3=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];







-(1/2)C22[{momentumshift[[1]],momentumshift[[3]],momentumshift[[4]]}]
+1/2 C22[{momentumshift[[1]],momentumshift[[2]],momentumshift[[4]]}]

-(prod[mom2,mom2]+2 prod[mom1,mom2])/2 D33[momentumshift]




]




Rd333[momentumshift_]:=Module[
{mom1=0,mom2=0,mom3=0},

mom1=momentumshift[[2]]-momentumshift[[1]];

mom2=momentumshift[[3]]-momentumshift[[2]];

mom3=momentumshift[[4]]-momentumshift[[3]];




  -(1/2)C22[{momentumshift[[1]],momentumshift[[2]],momentumshift[[4]]}]

-1/2 (prod[mom3,mom3]+2 prod[mom1,mom3]+2 prod[mom2,mom3])D33[momentumshift]

-2 D003[momentumshift]

]




(**Dijk solutions **)


Dijkeqindex[i_,j_]:=Which[

i==1,Which[
j==1,4,
j==2,1,
j==3,2
],

i==2,Which[
j==2,5,
j==3,3
],

i==3,6,

i==0,7

]

code5[i_,j_,k_]:="[momentumshift_]:=Module[\[IndentingNewLine]\[IndentingNewLine]{Graham=getGrahamMatrixBox[momentumshift]//Simplify,\[IndentingNewLine]\[IndentingNewLine]\[IndentingNewLine]Rarray={\[IndentingNewLine]\[IndentingNewLine]{Rd121[momentumshift],Rd122[momentumshift],Rd123[momentumshift]},\[IndentingNewLine]{Rd131[momentumshift],Rd132[momentumshift],Rd133[momentumshift]},\[IndentingNewLine]{Rd231[momentumshift],Rd232[momentumshift],Rd233[momentumshift]},\[IndentingNewLine]{Rd111[momentumshift],Rd112[momentumshift],Rd113[momentumshift]},\[IndentingNewLine]{Rd221[momentumshift],Rd222[momentumshift],Rd223[momentumshift]},\[IndentingNewLine]{Rd331[momentumshift],Rd332[momentumshift],Rd333[momentumshift]},\[IndentingNewLine]{Rd001[momentumshift],Rd002[momentumshift],Rd003[momentumshift]}\[IndentingNewLine]\[IndentingNewLine]},\[IndentingNewLine]\[IndentingNewLine]index=1,\[IndentingNewLine]\[IndentingNewLine]invertible=1,\[IndentingNewLine]sol=0,\[IndentingNewLine]output=0},\[IndentingNewLine]\[IndentingNewLine]index=Dijkeqindex["<>ToString[i]<>","<>ToString[j]<>"];\[IndentingNewLine]\[IndentingNewLine](*use solution only if Graham matrix is invertible*)\[IndentingNewLine]\[IndentingNewLine]invertible=If[Det[Graham]\[Equal]0,0,1,1];\[IndentingNewLine]\[IndentingNewLine]If[invertible\[Equal]1,\[IndentingNewLine]\[IndentingNewLine]\[IndentingNewLine]sol=Inverse[Graham].(Rarray[[index]]);\[IndentingNewLine]\[IndentingNewLine]sol[["<>ToString[k]<>"]],\[IndentingNewLine]\[IndentingNewLine]\[IndentingNewLine]D"<>ToString[i]<>""<>ToString[j]<>""<>ToString[k]<>"[momentumshift],(*if Graham matrix is singular, leave unevaluated, for now*)\[IndentingNewLine]\[IndentingNewLine]D"<>ToString[i]<>""<>ToString[j]<>""<>ToString[k]<>"[momentumshift] (*if Graham matrix is singular, leave unevaluated, for now*)\[IndentingNewLine]\[IndentingNewLine]]\[IndentingNewLine]]"


temp=Table["D00"<>ToString[k]<>"eval"<>code5[0,0,k],{k,1,3}];

stringfunct5=Join[Table["D"<>ToString[i]<>ToString[j]<>ToString[k]<>"eval"<>code5[i,j,k],{i,1,3},{j,i,3},{k,j,3}],temp];


stringfunct5//ToExpression;
(*BW:symmetrize Dijk*)
Do[
"D"<>ToString[i]<>ToString[k]<>ToString[j]<>"eval=D"<>ToString[i]<>ToString[j]<>ToString[k]<>"eval"//ToExpression;
"D"<>ToString[j]<>ToString[i]<>ToString[k]<>"eval=D"<>ToString[i]<>ToString[j]<>ToString[k]<>"eval"//ToExpression;
"D"<>ToString[j]<>ToString[k]<>ToString[i]<>"eval=D"<>ToString[i]<>ToString[j]<>ToString[k]<>"eval"//ToExpression;
"D"<>ToString[k]<>ToString[i]<>ToString[j]<>"eval=D"<>ToString[i]<>ToString[j]<>ToString[k]<>"eval"//ToExpression;
"D"<>ToString[k]<>ToString[j]<>ToString[i]<>"eval=D"<>ToString[i]<>ToString[j]<>ToString[k]<>"eval"//ToExpression,
{i,1,3},{j,i,3},{k,j,3}];

BOXFF[term_]:=term/.box->boxeval/.box3->box3eval/.box2->box2eval/.box1->box1eval


(*Clear[D233eval]
D233eval[momentumshift_]:=Module[

{Graham=getGrahamMatrixBox[momentumshift]//Simplify,


Rarray={

{Rd121[momentumshift],Rd122[momentumshift],Rd123[momentumshift]},
{Rd131[momentumshift],Rd132[momentumshift],Rd133[momentumshift]},
{Rd231[momentumshift],Rd232[momentumshift],Rd233[momentumshift]},
{Rd111[momentumshift],Rd112[momentumshift],Rd113[momentumshift]},
{Rd221[momentumshift],Rd222[momentumshift],Rd223[momentumshift]},
{Rd331[momentumshift],Rd332[momentumshift],Rd333[momentumshift]},
{Rd001[momentumshift],Rd002[momentumshift],Rd003[momentumshift]}

},

index=1,

invertible=1,
sol=0,
output=0},

index=Dijkeqindex[2,3];

(*use solution only if Graham matrix is invertible*)

invertible=If[Det[Graham]\[Equal]0,0,1,1];

If[invertible\[Equal]1,


sol=Inverse[Graham].(Rarray[[index]]);

sol[[3]],


D233[momentumshift],(*if Graham matrix is singular, leave unevaluated, for now*)

D233[momentumshift] (*if Graham matrix is singular, leave unevaluated, for now*)

]
]*)






(*Replace expressions from FORM.py to different integral types.*)

(*prepare terms in table. *)

gettable[term_]:=Module[{expanded={}},

(*first replacement to simplify terms kdot[k1,k1], which reduce the rank of the integral *)

expanded=(term/.kden[-k1]->kden[k1]/.kdot[k1,k1]->1/kden[k1])//Expand//Cancel//Expand;

Table[expanded[[i]],{i,1,Length[expanded]}]

];

(*get ranks *)

getranks[termarray_]:=termarray/. kdot[_,_]->power/.__ power^p_->p/.__ power->1/.kden[__]->0;

(*get indices organized by rank*)

getrankindices[ranks_]:={Join@@Position[ranks,0],Join@@Position[ranks,1],Join@@Position[ranks,2],Join@@Position[ranks,3],Join@@Position[ranks,4]};

(*get number of propagatros*)

getnumberofprop[termarray_]:=termarray/. kden[__]->power/.__ power->1/.__ power^p_ ->p;

getmomentumshifts[termarray_]:=Module[{step1={}},
step1=(Cases[#,kden[a_]->a]&/@termarray);
step1/.Plus[Times[-1,k1],a_]->Plus[k1,-a]/.Times[-1,k1]->k1/.k1->0

]

(*getmomentumshifts[termarray_]:=(Cases[#,kden[a_]\[Rule]-a]&/@termarray)/.k1\[Rule]0;*)

(*momenta in the numerator*)

getnummomenta[termarray_]:=Module[{list1={},list2={}},
list1=Cases[#,kdot[k1,b_]->{1,b}]&/@termarray;
list2=Cases[#,kdot[k1,b_]^p_->{p,b}]&/@termarray;
Join[list1,list2,2]
]
(*BW:In some terms kdot[k1,k1] cancels kden[k1]. If there are other kdots in numerator, the momentum of this term must be shifted. This is done by the following function*)
MOMSHIFT[termarray_]:=
Module[
{table1={},momentumshifts={}},
table1=termarray//gettable;
momentumshifts=table1//getmomentumshifts;
Do[
If[
momentumshifts[[i]][[1]]=!=0,
table1[[i]]=table1[[i]]/.{kden[k1+k_]->kden[k1-momentumshifts[[i]][[1]]+k],kden[-k1+k_]->kden[-k1+momentumshifts[[i]][[1]]+k],kdot[k1,a_]->kdot[k1,a]-momentumshifts[[i]][[1]].a}/.(-a_).b_->-a.b/.setvars
],
{i,1,Length[table1]}
];
Sum[table1[[i]],{i,1,Length[table1]}]//Expand
]

(*generate table ready for replacement*)
TENSORINTEGRALS[termarray_]:=
Module[
{table1={},table2={},table3={},table4={},nummomenta={},momentumshifts={},termarray1},
termarray1=termarray//MOMSHIFT;
table1=termarray1//gettable;
nummomenta=table1//getnummomenta;
momentumshifts=table1//getmomentumshifts;
table2=(table1/.kden[__]->1/.kdot[_,_]->1);
table3=Table[IntegralType[nummomenta[[i]],momentumshifts[[i]]],{i,1,Length[nummomenta]}];
table4=Times@@{table3,table2};

Plus@@table4


];



(*Reduction chains*)
(*put together all form factor replacements*)

FORMFACTORS[term_]:=(term//BOXFF//TRIANGLEFF//BUBBLEFF)//.setparvars;


(*reduction of box form factors*)
boxstep1=Join@@Join@@Table["D"<>ToString[i]<>ToString[j]<>ToString[k]<>"->"<>"D"<>ToString[i]<>ToString[j]<>ToString[k]<>"eval",{i,1,3},{j,1,3},{k,1,3}];
boxstep2=Table["D00"<>ToString[k]<>"->"<>"D00"<>ToString[k]<>"eval",{k,1,3}];
boxstep3=Join@@Table["D"<>ToString[i]<>ToString[j]<>"->"<>"D"<>ToString[i]<>ToString[j]<>"eval",{i,1,3},{j,1,3}];
boxstep4="{D00->D00eval}";
boxstep5=Table["D"<>ToString[i]<>"v->"<>"D"<>ToString[i]<>"eval",{i,1,3}];

REDUCEBOX[term_]:=term/.ToExpression[boxstep1]/.ToExpression[boxstep2]/.ToExpression[boxstep3]/.ToExpression[boxstep4]/.ToExpression[boxstep5]


(*reduction of triangle form factors*)
trianglestep1=Join@@Table["C"<>ToString[i]<>ToString[j]<>"->"<>"C"<>ToString[i]<>ToString[j]<>"eval",{i,1,2},{j,i,2}];
trianglestep2="{C00->C00eval}";
trianglestep3=Table["C"<>ToString[i]<>"v->"<>"C"<>ToString[i]<>"eval",{i,1,2}];

REDUCETRIANGLE[term_]:=term/.ToExpression[trianglestep1]/.ToExpression[trianglestep2]/.ToExpression[trianglestep3];

(*reduction of bubble form factors*)

bubblestep1=Table["B"<>ToString[i]<>"v->"<>"B"<>ToString[i]<>"eval",{i,1,1}]

REDUCEBUBBLE[term_]:=term/.ToExpression[bubblestep1]



PASSARINOVELTMAN[term_]:=(term//TENSORINTEGRALS//FORMFACTORS//REDUCEBOX//REDUCETRIANGLE//REDUCEBUBBLE//Expand)/.dim->4-2\[Epsilon]


\!\(TraditionalForm\`{"\<B1->B1eval\>"}\)




(*Scalar loop integrals*)


r\[CapitalGamma]=Gamma[1-\[Epsilon]]^2*Gamma[1+\[Epsilon]]/Gamma[1-2 \[Epsilon]];

(*prefactor[\[Mu]_]:=I r\[CapitalGamma]/16/\[Pi]^2;*)(*BW:keep dim general for pi?*)
prefactor[\[Mu]_]:=I r\[CapitalGamma]/2^(4-2\[Epsilon])/\[Pi]^(2-\[Epsilon]); (*BW: Note here to keep the dimension of the integral, the factor \[Mu]^(2\[Epsilon]) is not removed from the integral. The net dimension factor for the squared amplitude is \[Mu]^(4\[Epsilon]) from gs^4. Do not forget to include the other \[Mu]^(2\[Epsilon]) factor in the two body phase space in the end*)

B0eval[momshifts_,\[Mu]_]:=Module[
{
mom=0,
mom2,mom2val,
output=0,
test=2
},

mom=momshifts[[2]]-momshifts[[1]]//Simplify;

mom2val=prod[mom,mom];

test=testspacetimelight[mom];



(*prefactor[\[Mu]]*(-4 \[Pi] \[Mu]^2/mom2)^\[Epsilon]*(1/\[Epsilon]+2)/.\[Mu]^2\[Rule]\[Mu]^2 Exp[EulerGamma]/4/\[Pi]//Series[#,{\[Epsilon],0,0}]&//Normal*)

(*output=\[ImaginaryI]/(16 \[Pi]^2 \[Epsilon])+(\[ImaginaryI] (2+Log[-( (\[Mu]^2)/mom2)]))/(16 \[Pi]^2);*)(*in MSbar*)

(*BW: prefactor[\[Mu]]*(-4 \[Pi] \[Mu]^2/mom2)^\[Epsilon]*(1/\[Epsilon]+2)//Series[#,{\[Epsilon],0,0}]&*)
(*output=\[ImaginaryI]/(16 \[Pi]^2 \[Epsilon])-(\[ImaginaryI] (-2+EulerGamma-Log[-((4 \[Pi] \[Mu]^2)/mom2)]))/(16 \[Pi]^2);*)
output=I/(16 \[Pi]^2 \[Epsilon])-(I (-2+E4Pi-Log[-(\[Mu]^2/mom2)]))/(16 \[Pi]^2);
Which[

test==-1,output/.mom2->mom2val(*/.\[Epsilon]\[Rule]\[Epsilon]UV*),

test==1, output/.Log[a_]->Log[-a]+I \[Pi]/.mom2->mom2val(*/.\[Epsilon]\[Rule]\[Epsilon]UV*),(*BW:not -mom2val? Log[-a]?*)

(*test\[Equal]0,B0v[{0,mom},\[Mu]] *)
test==0,0

]



]



C0type[mom1_,mom2_,mom3_]:=Module[{test1=mom1//testspacetimelight,test2=mom2//testspacetimelight,test3=mom3//testspacetimelight},


If[(test1==0&&test2==0)||(test1==0&&test3==0)||(test2==0&&test3==0),1,2,0]


]

C0type1nonvanishing[mom1_,mom2_,mom3_]:=Module[(*exactly 2 args arguments have to already be know lightlike momenta, third space or time-like*)

{test1=mom1//testspacetimelight,test2=mom2//testspacetimelight,test3=mom3//testspacetimelight},


Which[

(test1==0&&test2==0),mom3,
(test1==0&&test3==0),mom2,
(test2==0&&test3==0),mom1


]


]

C0type2nonvanishing[mom1_,mom2_,mom3_]:=Module[(*exactly 1 args arguments have to already be know lightlike momenta, rest space or time-like*)

{test1=mom1//testspacetimelight,test2=mom2//testspacetimelight,test3=mom3//testspacetimelight},


Which[

(test1==0),{mom2,mom3},
(test2==0),{mom1,mom3},
(test3==0),{mom1,mom2}


]


]


C0type1eval[momshifts_,\[Mu]_]:=Module[
{mom=0,
output=0,
mom1=0,mom2=0,mom3=0,
type=0,
s,sval=0,
test
},




mom1=momshifts[[2]]-momshifts[[1]];

mom2=momshifts[[3]]-momshifts[[2]];

mom3=momshifts[[1]]-momshifts[[3]];

mom=C0type1nonvanishing[mom1,mom2,mom3];

sval=prod[mom,mom];

test=mom//testspacetimelight;


(*prefactor[\[Mu]]*(-4 \[Pi] \[Mu]^2/s)^\[Epsilon]/\[Epsilon]^2 /s/.\[Mu]^2\[Rule]\[Mu]^2 Exp[EulerGamma]/4/\[Pi]//Series[#,{\[Epsilon],0,0}]&//Normal(*in MSBAR*)*)


(*output=\[ImaginaryI]/(16 \[Pi]^2 s \[Epsilon]^2)+(\[ImaginaryI] Log[-(\[Mu]^2/s)])/(16 \[Pi]^2 s \[Epsilon])-(\[ImaginaryI] (\[Pi]^2-6 Log[-(\[Mu]^2/s)]^2))/(192 \[Pi]^2 s);*)

(*output=\[ImaginaryI]/(16 \[Pi]^2 s \[Epsilon]^2)-(\[ImaginaryI] (EulerGamma-log[-((4 \[Pi] \[Mu]^2)/s)]))/(16 \[Pi]^2 s \[Epsilon])+(\[ImaginaryI] (-\[Pi]^2+6 (EulerGamma-Log[4 \[Pi]])^2+6 log[-(\[Mu]^2/s)] (-2 EulerGamma+log[-((16 \[Pi]^2 \[Mu]^2)/s)])))/(192 \[Pi]^2 s);*)
output=I/(16 \[Pi]^2 s \[Epsilon]^2)-(I (E4Pi-log[-(\[Mu]^2/s)]))/(16 \[Pi]^2 s \[Epsilon])-(I (-6 E4Pi^2+\[Pi]^2+12 E4Pi log[-(\[Mu]^2/s)]-6 log[-(\[Mu]^2/s)]^2))/(192 \[Pi]^2 s);
(*Which[

test\[Equal]-1,output/.s\[Rule]sval,

test\[Equal]1, output/.Log[a_]\[Rule]Log[-a]+I \[Pi]/.s\[Rule]sval(*BW:not -sval? Log[-a]?*),

test\[Equal]0,C0v[momshifts,\[Mu]] 

]*)

Which[

test==-1,output/.log[a_]->Log[a]/.s->sval,

test==1, output/.log[a_]->Log[-a]+I \[Pi]/.s->sval,(*as in BW's*)

test==0,C0v[momshifts,\[Mu]] 

]


]




C0type2eval[momshifts_,\[Mu]_]:=Module[
{mom=0,
output=0,
mom1=0,mom2=0,mom3=0,
type=0,
s1,s2,s3,
sval={},
test={}
},




mom1=momshifts[[2]]-momshifts[[1]];

mom2=momshifts[[3]]-momshifts[[2]];

mom3=momshifts[[1]]-momshifts[[3]];

mom=C0type2nonvanishing[mom1,mom2,mom3];

sval={  prod[mom[[1]],mom[[1]]],prod[mom[[2]],mom[[2]]]  };


test=testspacetimelight/@mom;



(*prefactor[\[Mu]]*((-4 \[Pi] \[Mu]^2/s1)^\[Epsilon]-(-4 \[Pi] \[Mu]^2/s2)^\[Epsilon])/\[Epsilon]^2 /(s1-s2)/.\[Mu]^2\[Rule]\[Mu]^2 Exp[EulerGamma]/4/\[Pi]//Series[#,{\[Epsilon],0,0}]&//Normal;*)

(*output=type2c*((\[ImaginaryI] (Log1[-(\[Mu]^2/s1)]-Log2[-(\[Mu]^2/s2)]))/(16 \[Pi]^2 (s1-s2) \[Epsilon])+(\[ImaginaryI] (Log1[-(\[Mu]^2/s1)]^2-Log2[-(\[Mu]^2/s2)]^2))/(32 \[Pi]^2 (s1-s2)));(*in MSBAR*)*)
(*output=(\[ImaginaryI] (Log[-(\[Mu]^2/s1)]-Log[-(\[Mu]^2/s2)]))/(16 \[Pi]^2 (s1-s2) \[Epsilon])+(\[ImaginaryI] (Log[-(\[Mu]^2/s1)]-Log[-(\[Mu]^2/s2)]) (-2 EulerGamma+Log[-((16 \[Pi]^2 \[Mu]^2)/s1)]+Log[-(\[Mu]^2/s2)]))/(32 \[Pi]^2 (s1-s2)); *)(*not needed for prompt photon*)
(*output=(\[ImaginaryI] (log1[-(\[Mu]^2/s1)]-log2[-(\[Mu]^2/s2)]))/(16 \[Pi]^2 (s1-s2) \[Epsilon])+(\[ImaginaryI] (log1[-(\[Mu]^2/s1)]-log2[-(\[Mu]^2/s2)]) (-2 EulerGamma+log1[-((16 \[Pi]^2 \[Mu]^2)/s1)]+log2[-(\[Mu]^2/s2)]))/(32 \[Pi]^2 (s1-s2)); *)(*as in BW's*)
output=(I (log1[-(\[Mu]^2/s1)]-log2[-(\[Mu]^2/s2)]))/(16 \[Pi]^2 (s1-s2) \[Epsilon])-(I (2 E4Pi log1[-(\[Mu]^2/s1)]-log1[-(\[Mu]^2/s1)]^2-2 E4Pi log2[-(\[Mu]^2/s2)]+log2[-(\[Mu]^2/s2)]^2))/(32 \[Pi]^2 (s1-s2));
(*Which[

(test[[1]]\[Equal]-1)&&(test[[2]]\[Equal]-1),output/.s1\[Rule]sval[[1]]/.s2\[Rule]sval[[2]],

(test[[1]]\[Equal]   1)&&(test[[2]]\[Equal]-1),output/.log1[a_]\[Rule]Log[a]+I \[Pi]/.s1\[Rule]sval[[1]]/.s2\[Rule]sval[[2]](*BW:not s1\[Rule]-sval[[1]]?*),

(test[[1]]\[Equal]-1)&&(test[[2]]\[Equal]   1),output/.log2[a_]\[Rule]Log[a]+I \[Pi]/.s2\[Rule]sval[[2]]/.s1\[Rule]sval[[1]](*BW:not s2\[Rule]-sval[[2]]?*),

(test[[1]]\[Equal]1)&&(test[[2]]\[Equal]1),output/.log1[a_]\[Rule]Log[a]+I \[Pi]/.s1\[Rule]sval[[1]]/.log2[a_]\[Rule]Log[a]+I \[Pi]/.s2\[Rule]sval[[2]]
(*BW:not s1\[Rule]-sval[[1]],s2\[Rule]-sval[[2]]?*)
]*)

Which[

(test[[1]]==-1)&&(test[[2]]==-1),output/.log1[a_]->Log[a]/.s1->sval[[1]]/.log2[a_]->Log[a]/.s2->sval[[2]],

(test[[1]]==   1)&&(test[[2]]==-1),output/.log1[a_]->Log[-a]+I \[Pi]/.s1->sval[[1]]/.log2[a_]->Log[a]/.s2->sval[[2]],

(test[[1]]==-1)&&(test[[2]]==   1),output/.log1[a_]->Log[a]/.s1->sval[[1]]/.log2[a_]->Log[-a]+I \[Pi]/.s2->sval[[2]],

(test[[1]]==1)&&(test[[2]]==1),output/.log1[a_]->Log[-a]+I \[Pi]/.s1->sval[[1]]/.log2[a_]->Log[-a]+I \[Pi]/.s2->sval[[2]]

]


]





C0eval[momshifts_,\[Mu]_]:=Module[
{mom=0,
mom1,mom2,mom3,
output1=0,
type=0
},




mom1=momshifts[[2]]-momshifts[[1]];

mom2=momshifts[[3]]-momshifts[[2]];

mom3=momshifts[[1]]-momshifts[[3]];

type=C0type[mom1,mom2,mom3];




Which[

type==1,C0type1eval[momshifts,\[Mu]],

type==2,C0type2eval[momshifts,\[Mu]]


]



]










D0type[mom1_,mom2_,mom3_,mom4_,mom5_,mom6_]:=Module[{test1=mom1//testspacetimelight,test2=mom2//testspacetimelight,test3=mom3//testspacetimelight,test4=mom4//testspacetimelight,test5=mom5//testspacetimelight,test6=mom6//testspacetimelight,type=0},


If[test1!=0,type++];
If[test2!=0,type++];
If[test3!=0,type++];
If[test4!=0,type++];
If[test5!=0,type++];
If[test6!=0,type++];

If[

type>3,3,

type-1

]

]

D0type1eval[momshifts_,\[Mu]_]:=Module[
{mom={},
mom1,mom2,mom3,mom4,mom5,mom6,
s12,s23,
sij={},
sijtype={},(*space- or time-like*)
output=0,
type=0
},


mom1=momshifts[[2]]-momshifts[[1]];

mom2=momshifts[[3]]-momshifts[[2]];

mom3=momshifts[[4]]-momshifts[[3]];

mom4=-(mom1+mom2+mom3);

mom5=mom1+mom2;

mom6=mom2+mom3;

mom={mom1,mom2,mom3,mom4,mom5,mom6};

(*keep only non light-like momenta*)

mom=Table[If[(mom[[i]]//testspacetimelight)!=0,mom[[i]],0],{i,1,mom//Length}]//DeleteCases[#,0]&;

(*mom shifts should contain only two non-zero invariants*)

sij=prod[#,#]&/@mom//Simplify;

(*check spacelike or timelike*)

sijtype=testspacetimelight/@mom;


(*prefactor[\[Mu]]*(2/\[Epsilon]^2(((4 \[Pi] \[Mu]^2)/(-s12-I \[Epsilon]F))^\[Epsilon]+((4 \[Pi] \[Mu]^2)/(-s23-I \[Epsilon]F))^\[Epsilon])-( Log[-s12-I \[Epsilon]F] - Log[-s23 - I \[Epsilon]F] )^2 -\[Pi]^2)/s12/s23/.\[Mu]^2\[Rule]\[Mu]^2 Exp[EulerGamma]/4/\[Pi]//Series[#,{\[Epsilon],0,0}]&//Normal;*)

(*output=\[ImaginaryI]/(4 \[Pi]^2 s12 s23 \[Epsilon]^2)+(\[ImaginaryI] (log1[\[Mu]^2/-s12]+log2[\[Mu]^2/-s23]))/(8 \[Pi]^2 s12 s23 \[Epsilon])-(\[ImaginaryI] (4 \[Pi]^2+3 (log1[-s12]-log2[-s23])^2-3 log1[\[Mu]^2/-s12]^2-3 log2[\[Mu]^2/-s23]^2))/(48 \[Pi]^2 s12 s23)  ;*)
(*output=\[ImaginaryI]/(4 \[Pi]^2 s12 s23 \[Epsilon]^2)-(\[ImaginaryI] (2 EulerGamma-log1[-((16 \[Pi]^2 \[Mu]^2)/s12)]-log2[-(\[Mu]^2/s23)])*1)/(8 \[Pi]^2 s12 s23 \[Epsilon])+(1*1)/(48 \[Pi]^2 s12 s23)\[ImaginaryI] (-4 \[Pi]^2+6 (EulerGamma-Log[4 \[Pi]])^2-3 (log1[-s12]-log2[-s23])^2+3 log1[-(\[Mu]^2/s12)] (-2 EulerGamma+log1[-((16 \[Pi]^2 \[Mu]^2)/s12)])+3 log2[-(\[Mu]^2/s23)] (-2 EulerGamma+log2[-((16 \[Pi]^2 \[Mu]^2)/s23)]));*)
output=I/(4 \[Pi]^2 s12 s23 \[Epsilon]^2)-(I (2 E4Pi-log1[-(\[Mu]^2/s12)]-log2[-(\[Mu]^2/s23)]))/(8 \[Pi]^2 s12 s23 \[Epsilon])+1/(48 \[Pi]^2 s12 s23) I (6 E4Pi^2-4 \[Pi]^2-3 (log1[-s12]-log2[-s23])^2-6 E4Pi log1[-(\[Mu]^2/s12)]+3 log1[-(\[Mu]^2/s12)]^2+3 log2[-(\[Mu]^2/s23)] (-2 E4Pi+log2[-(\[Mu]^2/s23)]));
If[
sijtype[[1]]==1,
output=output/.log1[a_]->Log[-a]+I \[Pi],
output=output/.log1[a_]->Log[a]          
];


If[
sijtype[[2]]==1,
output=output/.log2[a_]->Log[-a]+I \[Pi],
output=output/.log2[a_]->Log[a]          
];

output=output/.{s12->sij[[1]],s23->   sij[[2]]};

output//ReleaseHold



]


D0type2eval[momshifts_,\[Mu]_]:=Module[
{mom={},
mom1,mom2,mom3,mom4,mom5,mom6,
s12,s23,p42,
sij={},
sijtype={},(*space- or time-like*)
output=0,
type=0
},


mom1=momshifts[[2]]-momshifts[[1]];

mom2=momshifts[[3]]-momshifts[[2]];

mom3=momshifts[[4]]-momshifts[[3]];

mom4=-(mom1+mom2+mom3);

mom5=mom1+mom2;

mom6=mom2+mom3;

mom={mom1,mom2,mom3,mom4,mom5,mom6};

(*keep only non light-like momenta*)

mom=Table[If[(mom[[i]]//testspacetimelight)!=0,mom[[i]],0],{i,1,mom//Length}]//DeleteCases[#,0]&;

(*mom shifts should contain only two non-zero invariants*)
mom=DeleteCases[mom,Select[mom,prod[#,#]==-Q^2&][[1]]];
sij=prod[#,#]&/@mom//Simplify;
Print[mom];
Print[sij];
(*check spacelike or timelike*)

sijtype=testspacetimelight/@mom;


(*prefactor[\[Mu]]*(2/\[Epsilon]^2(((4 \[Pi] \[Mu]^2)/(-s12-I \[Epsilon]F))^\[Epsilon]+((4 \[Pi] \[Mu]^2)/(-s23-I \[Epsilon]F))^\[Epsilon])-( Log[-s12-I \[Epsilon]F] - Log[-s23 - I \[Epsilon]F] )^2 -\[Pi]^2)/s12/s23/.\[Mu]^2\[Rule]\[Mu]^2 Exp[EulerGamma]/4/\[Pi]//Series[#,{\[Epsilon],0,0}]&//Normal;*)

(*output=\[ImaginaryI]/(4 \[Pi]^2 s12 s23 \[Epsilon]^2)+(\[ImaginaryI] (log1[\[Mu]^2/-s12]+log2[\[Mu]^2/-s23]))/(8 \[Pi]^2 s12 s23 \[Epsilon])-(\[ImaginaryI] (4 \[Pi]^2+3 (log1[-s12]-log2[-s23])^2-3 log1[\[Mu]^2/-s12]^2-3 log2[\[Mu]^2/-s23]^2))/(48 \[Pi]^2 s12 s23)  ;*)
(*output=\[ImaginaryI]/(8 \[Pi]^2 s12 s23 \[Epsilon]^2)-(\[ImaginaryI] (EulerGamma+log1[-(1/(4 \[Pi] p42))]-log2[-(1/s12)]-log3[-(1/s23)]))/(8 \[Pi]^2 s12 s23 \[Epsilon])-1/(32 \[Pi]^2 s12 s23)\[ImaginaryI] (2 EulerGamma^2-\[Pi]^2-4 Li2[1-f1m p42]+4 Li2[1-f1m s12]+4 Li2[1-f1m s23]+8 Log[2]^2-2 log1[-(1/p42)]^2+4 log1[-(1/p42)] (EulerGamma-Log[4 \[Pi]])+2 Log[\[Pi]] Log[16 \[Pi]]-4 EulerGamma log2[-((4 \[Pi])/s12)]+2 log2[-(1/s12)] log2[-((16 \[Pi]^2)/s12)]+2 log3[-(1/s23)] (-2 EulerGamma+log3[-((16 \[Pi]^2)/s23)]));*)
output=I/(8 \[Pi]^2 s12 s23 \[Epsilon]^2)-(I (E4Pi+log1[-(\[Mu]^2/p42)]-log2[-(\[Mu]^2/s12)]-log3[-(\[Mu]^2/s23)]))/(8 \[Pi]^2 s12 s23 \[Epsilon])-1/(32 \[Pi]^2 s12 s23) I (-2 E4Pi^2+\[Pi]^2+4 Li21[1-p42/s12]+4 Li22[1-p42/s23]+2 log4[s12/s23]^2-4 E4Pi log1[-(\[Mu]^2/p42)]+2 log1[-(\[Mu]^2/p42)]^2+4 E4Pi log2[-(\[Mu]^2/s12)]-2 log2[-(\[Mu]^2/s12)]^2+4 E4Pi log3[-(\[Mu]^2/s23)]-2 log3[-(\[Mu]^2/s23)]^2);
(*If[
sijtype[[1]]\[Equal]1,
output=output/.log1[a_]\[Rule]Log[-a]+I \[Pi]/.p42\[Rule]sij[[1]](*BW:not -sij[[1]]?*),
output=output/.log1[a_]\[Rule]Log[a]          /.p42\[Rule]  sij[[1]]
];*)
output=output/.log1[a_]->Log[a] ;

If[
sijtype[[1]]==1,
output=output/.log2[a_]->Log[-a]+I \[Pi]/.Li21[a_]->1/3Pi^2-1/2Log[a]^2-I \[Pi]*Log[a]-PolyLog[2,1/a],


output=output/.log2[a_]->Log[a]/.Li21[a_]->PolyLog[2,a]
];

If[
sijtype[[2]]==1,
output=output/.log3[a_]->Log[-a]+I \[Pi]/.Li22[a_]->1/3Pi^2-1/2Log[a]^2-I \[Pi]*Log[a]-PolyLog[2,1/a],

output=output/.log3[a_]->Log[a] /.Li22[a_]->PolyLog[2,a]         
];

(*If[
(sijtype[[1]]==1&&sijtype[[2]]==-1)||(sijtype[[1]]==-1&&sijtype[[2]]==1),
output=output/.log4[a_]->Log[-a]+I \[Pi],
output=output/.log4[a_]->Log[a]          
];old implementation, cannot handle ln(s/u)*)
Which[
(sijtype[[1]]==1&&sijtype[[2]]==-1),
output=output/.log4[a_]->Log[-a]-I \[Pi],
(sijtype[[1]]==-1&&sijtype[[2]]==1),
output=output/.log4[a_]->Log[-a]+I \[Pi],
True,
output=output/.log4[a_]->Log[a]          
];
output=output/.{p42->  -Q^2,s12->sij[[1]],s23->sij[[2]]};
output//ReleaseHold



]



D0eval[momshifts_,\[Mu]_]:=Module[
{mom=0,
mom1,mom2,mom3,mom4,mom5,mom6,
output1=0,
type=0
},




mom1=momshifts[[2]]-momshifts[[1]];

mom2=momshifts[[3]]-momshifts[[2]];

mom3=momshifts[[1]]-momshifts[[3]];

mom4=-(mom1+mom2+mom3);

mom5=mom1+mom2;

mom6=mom2+mom3;

type=D0type[mom1,mom2,mom3,mom4,mom5,mom6];




Which[

type==1,D0type1eval[momshifts,\[Mu]],

type==2,D0type2eval[momshifts,\[Mu]],

type==3,Print["type 3 of D0 not implemented "];Abort[]
]

]





removeUVtemp[term_]:=term/.\[Epsilon]UV^-1->1;

ORGANIZEPOLES[term_,\[Epsilon]_]:=

Module[{test1=1,test2=1,termmod={} ,termmod2={} ,finite={},vanishing={},singlepole={},doublepole={}},
(*separatation before overall factor 'integralnorm'*)
termmod=Collect[term//Expand,{\[Epsilon]}];

doublepole=termmod//Select[#,MatchQ[#,\[Epsilon]^-2 __]&]&;

singlepole=termmod//Select[#,MatchQ[#,\[Epsilon]^-1 __]&]&;

finite=termmod//Select[#,FreeQ[#,\[Epsilon] __]&&FreeQ[#,\[Epsilon]^p_ __]&]&;

vanishing=termmod-doublepole-singlepole-finite;


test1=(vanishing/.\[Epsilon]->0);

test2=((termmod-doublepole-singlepole-finite-vanishing)//Expand);

If[test1==0&&test2==0,Print["ok"],Print["wrong"]];

termmod={doublepole,singlepole,finite,(vanishing/.\[Epsilon]->0)}//Simplify//Collect[#,{\[Epsilon]}]&




]

ORGANIZEPOLES2[termarray_,\[Epsilon]_]:=

Module[{test1=1,test2=1,termarraymod={} ,termarraymod2={} ,finite={},vanishing={},singlepole={},doublepole={}},
(*separatation before overall factor 'integralnorm'*)
termarraymod=Collect[termarray//Expand,{\[Epsilon]},Simplify];

doublepole=Select[#,MatchQ[#,\[Epsilon]^-2 __]&]&/@termarraymod;

singlepole=Select[#,MatchQ[#,\[Epsilon]^-1 __]&]&/@termarraymod;

finite=Select[#,FreeQ[#,\[Epsilon] __]&&FreeQ[#,\[Epsilon]^p_ __]&]&/@termarraymod;

vanishing=termarraymod-doublepole-singlepole-finite;


test1=(Plus@@vanishing/.\[Epsilon]->0);

test2=(Plus@@(termarraymod-doublepole-singlepole-finite-vanishing)//Expand);

If[test1==0&&test2==0,Print["ok"],Print["wrong"]];

termarraymod=Plus@@@{doublepole,singlepole,finite,(vanishing/.\[Epsilon]->0)}//Collect[#,{\[Epsilon]}]&


]



EVALB0[term_]:=term/.B0v[a_]->B0v[a,\[Mu]]/.B0v->B0eval;

EVALC0[term_]:=term/.C0v[a_]->C0v[a,\[Mu]]/.C0v->C0eval;

EVALD0[term_]:=term/.D0v[a_]->D0v[a,\[Mu]]/.D0v->D0eval;

EVALSCALARINT[term_]:=term//EVALB0//EVALC0//EVALD0;


virtual[processname_]:=Module[
{psfac,term1,term2,total,imaginary,real,final,SQAMPs,
v00,v01,v02,v03,v04,v05,v06,v07,v08,v09,v010,v10,v11,v12,v13,v14,v15,v16,v17,v18,v19,v110},
SetDirectory[NotebookDirectory[]];
SQAMPs = "interf_" <> processname <> ".m" ;
Import[SQAMPs];

v00=interf00//PASSARINOVELTMAN//EVALSCALARINT//Simplify;
v01=interf01//PASSARINOVELTMAN//EVALSCALARINT//Simplify;
v02=interf02//PASSARINOVELTMAN//EVALSCALARINT//Simplify;
v03=interf03//PASSARINOVELTMAN//EVALSCALARINT//Simplify;
v04=interf04//PASSARINOVELTMAN//EVALSCALARINT//Simplify;
v05=interf05//PASSARINOVELTMAN//Collect[#,{A0v[__],B0v[__],C0v[__],D0v[__],B1v[__]},Simplify]&//EVALSCALARINT;
v06=interf06//PASSARINOVELTMAN//Collect[#,{A0v[__],B0v[__],C0v[__],D0v[__],B1v[__]},Simplify]&//EVALSCALARINT;
v07=interf07//PASSARINOVELTMAN//Collect[#,{A0v[__],B0v[__],C0v[__],D0v[__],B1v[__]},Simplify]&//EVALSCALARINT;
v08=interf08//PASSARINOVELTMAN//Collect[#,{A0v[__],B0v[__],C0v[__],D0v[__],B1v[__]},Simplify]&//EVALSCALARINT;
v09=interf09//PASSARINOVELTMAN//Collect[#,{A0v[__],B0v[__],C0v[__],D0v[__],B1v[__]},Simplify]&//EVALSCALARINT;
v010=interf010//PASSARINOVELTMAN//Collect[#,{A0v[__],B0v[__],C0v[__],D0v[__],B1v[__]},Simplify]&//EVALSCALARINT;

v10=interf10//PASSARINOVELTMAN//EVALSCALARINT//Simplify;
v11=interf11//PASSARINOVELTMAN//EVALSCALARINT//Simplify;
v12=interf12//PASSARINOVELTMAN//EVALSCALARINT//Simplify;
v13=interf13//PASSARINOVELTMAN//Collect[#,{A0v[__],B0v[__],C0v[__],D0v[__],B1v[__]},Simplify]&//EVALSCALARINT//Simplify;
v14=interf14//PASSARINOVELTMAN//EVALSCALARINT//Simplify;
v15=interf15//PASSARINOVELTMAN//Collect[#,{A0v[__],B0v[__],C0v[__],D0v[__],B1v[__]},Simplify]&//EVALSCALARINT//Simplify;
v16=interf16//PASSARINOVELTMAN//Collect[#,{A0v[__],B0v[__],C0v[__],D0v[__],B1v[__]},Simplify]&//EVALSCALARINT//Simplify;
v17=interf17//PASSARINOVELTMAN//Collect[#,{A0v[__],B0v[__],C0v[__],D0v[__],B1v[__]},Simplify]&//EVALSCALARINT//Simplify;
v18=interf18//PASSARINOVELTMAN//Collect[#,{A0v[__],B0v[__],C0v[__],D0v[__],B1v[__]},Simplify]&//EVALSCALARINT;
v19=interf19//PASSARINOVELTMAN//Collect[#,{A0v[__],B0v[__],C0v[__],D0v[__],B1v[__]},Simplify]&//EVALSCALARINT;
v110=interf110//PASSARINOVELTMAN//Collect[#,{A0v[__],B0v[__],C0v[__],D0v[__],B1v[__]},Simplify]&//EVALSCALARINT;
(* term1:  terms v0i and v1i for i=0,..7  can be fully reduced. No integrals left to solve  *)
(* term2:  terms v0i and v1i for i=8,9,10 haveterms with scalar integrals that dont correspond to any in the main reference.  *)
(*psfac=1/2/2/(1-\[Epsilon])/8*1/(2(s+Q^2))*-1/(4Pi(s+Q^2))*1/Gamma[1-\[Epsilon]]*((4Pi(Q^2+s)^2*\[Mu]^2)/((-t-s-Q^2)(s*t)))^\[Epsilon];*)(*for the moment put 1/(1-\[Epsilon]) here for gluon average to be consisitent with real part, the whole factor 1/2/2/(1-\[Epsilon])/8 is the spin color average, which is not included in FORM output in interf.m. Here 1/2 is put for spin average of virtual photon, which is unnecessary for a intermediate particle, but just to be consistent with other parts of calculations*)


psfac=avefac[processname]*prefac222virt;
term1=
Series[psfac*Expand[(v00+v01+v02+v03+v04+v05+v06+v07+
v10+v11+v12+v13+v14+v15+v16+v17)/.\[Epsilon]UV^-1->0],{\[Epsilon],0,0}]//Normal;

term2=Series[psfac*((v08+v09+v010+v18+v19+v110)/.\[Epsilon]UV^-1->0),{\[Epsilon],0,0}]//Normal;




(*ORGANIZEPOLES separates in array the order {-2,-2,0,>0} in \[Epsilon], setting \[Epsilon]\[Rule]0 at the end for terms >0.*)
(*Note that the last one may give \[Epsilon] or \[Epsilon]UV as argument, but terms like \[Epsilon]/\[Epsilon]UV may be disregarded this way. Better to set \[Epsilon]UV -> \[Epsilon] before using ORGANIZEPOLES. *)
(*UV poles come only from B0 type integrals, they can be easily identified *)

total = -1*(term1+term2)/.\[Epsilon]UV->\[Epsilon];
(*totalUV=ORGANIZEPOLES[total,\[Epsilon]UV];
UVpiece=totalUV[[2]]/.\[Epsilon]\[Rule]0//Simplify*)

imaginary=total//Expand//Select[#,MatchQ[#,__ Complex[_,_]]&]&;
real           =total//Expand//Select[#,FreeQ[#,__ Complex[_,_]]&]&;
total -imaginary-real //Simplify; (*should output 0, if separation is good *)

(*'final' has no UV piece and is real. Here, 'final' is given in an array organized by powers of \[Epsilon], using 'ORGANIZEPOLES' *)
final = 2*(real)(*/.\[Epsilon]UV\[Rule]\[Epsilon]*)//Expand//ORGANIZEPOLES[#,\[Epsilon]]&//Collect[#,{Log[x__]},Simplify]&
]



virtualdoublepole[xsc_]:=xsc[[1]]//FullSimplify



virtualsinglepole[xsc_]:=Module[
{},
((xsc[[2]]//Simplify)/.s23->0/.t->-t//collectLog1//Simplify)/.t->-t/.E4Pi->EulerGamma-Log[4Pi]//FullSimplify
]
