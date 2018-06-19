(* ::Package:: *)

SetDirectory[NotebookDirectory[]]


<< parfrac.m;
<< psint.m;
<< sqamp.m;
<< cntxsc.m;


real[processname_]:=Module[
{SQAMPTOT,symmfac,prefac,SQAMP2,SQAMP3,test1,xsc},

SQAMPTOT = "SQAMP" <> processname <> "TOT" // ToExpression;
symmfac = Which[
  processname == "Aqtoqgg" || processname == "Aqbtoqbgg" || 
   processname == "Aqtoqbqq" || processname == "Aqbtoqqbqb", 1/2,
  True, 1
  ];

(*prefac223=-1/(2*(s + Q^2))*1/(2^7*(Pi)^4*(s + Q^2)*
      Gamma[1 - 2 \[Epsilon]])*((16 Pi^2 (Q^2 + s)^2*\[Mu]^4)/(s23*
        u*(s*t + Q^2*s23)))^\[Epsilon]*symmfac;*)(*symmetry factor 1/2 in the end when p2 and p3 are \
momenta of same particles, to avoid double counting (Peskin P108)*)
prefac=prefac223*symmfac;(*(Eq.30)*)
SQAMP2 = SQAMPTOT // procpf1 // procpf2 // procpf3;
tmpsqamp2=SQAMP2;
test1=(SQAMP2 - SQAMPTOT)/SQAMP2 //. {t3 -> -t1 - t2 - s - 2 Q^2, 
     u3 -> -u1 - u2 - s - Q^2, 
     s13 -> s - s12 - s23} //. {u2 -> -t2 - s12 - s23 - Q^2} //. 
  Q -> Sqrt[s23 - s - u1 - t1] /. {t2 -> -100.8, s12 -> 50.2, 
  s -> 320.5, t1 -> -30.1, s23 -> 20.3, u1 -> -42.38, D -> 4.2, 
  g -> 1.2, gp -> 1.5};
If[test1>10^-5,Print["wrong partial fraction"];Abort[]];
SQAMP3 = SQAMP2 // ADMVSimp // CoefSimp // CoefSimp;(*need to simplify twice to get correct exponent of s23 for \
each term*)
tmpsqamp3=SQAMP3;
xsc = procint2[SQAMP3, prefac]
]


(*prefacborn[s_, t_] = 
  1/(2 (s + Q^2))*-1/(4 Pi (s + Q^2))*1/
   Gamma[1 - \[Epsilon]]*((
    4 Pi (Q^2 + s)^2*\[Mu]^2)/((-t - s - Q^2) (s*t + 
       Q^2*s23)))^\[Epsilon];*)
prefacborn[s_,t_]=prefac222subTR;
(*The total makeup dimension is *\[Mu]^(4\[Epsilon]). For the virtual process, a factor *\[Mu]^(2\[Epsilon]) is included 
in the definition of the loop integral, so only a factor *\[Mu]^(2\[Epsilon]) is included in prefac222. *)


subtraction[processname_]:=Module[
{proccntxsc,symmfac,cntxsc},

proccntxsc = "proccntxsc" <> processname // ToExpression;
symmfac = Which[
  processname == "Aqtoqgg" || processname == "Aqbtoqbgg" || 
   processname == "Aqtoqbqq" || processname == "Aqbtoqqbqb", 1/2,
  True, 1
  ];
cntxsc = proccntxsc[prefacborn]
]


realdoublepole[xsc_]:=Expand[Select[Expand[xsc],MatchQ[#,\[Epsilon]^-2__]&]]/.s23->0//Simplify


realsinglepole[xsc_,cntxsc_]:=Module[
{xsc1,xsc2,xsc3,ndp},
ndp=checknonedeltapole[xsc,cntxsc];
(*Print[ndp];*)
If[ndp===0,Print["OK: all none-delta poles cancel"],Print["None-delta poles not all canceled. Check!!! "];Abort[]];
xsc1=Expand[Select[Expand[xsc-cntxsc],MatchQ[#,\[Epsilon]^-1__]&&MatchQ[#,\[Delta][s23]__]&&FreeQ[#,s23^-1__]&]+Simplify[Select[Expand[xsc-cntxsc],MatchQ[#,\[Epsilon]^-1__]&&MatchQ[#,\[Delta][s23]__]&&MatchQ[#,s23^-1__]&]]]/.s23->0//Collect[#,{Log[__]},Simplify]&//Simplify[#,s>0&&t<0&&Q>0]&;
(*ignore warnings  which is from terms like 1/(1+1/s23) when setting s23\[Rule]0. The terms are 0 in the limit but mathematica gives warnings about them*)
xsc2=xsc1/.Abs[a_]->a/.t->-t//collectLog1//FullSimplify;(*1.It is easy to see or check that Abs[a_]\[Rule]a and Abs[a_]\[Rule]-a gives same result due to the log structure.
2. mathematica gives wrong separation for Log[-ab] for logs if there is a minus sign in it. So here I use t\[Rule]-t to remove the sign, and in the end do t\[Rule]-t again. *)
xsc3=xsc2/.t->-t//Simplify
]


SimpPlus1[exp_]:=Module[
{Plus1terms,nonPlusterms,out},
Plus1terms=Select[exp//Expand,MatchQ[#,Plus1B[B,s23]__]&]/Plus1B[B,s23];
nonPlusterms=Select[exp//Expand,FreeQ[#,Plus1B[B,s23]__]&];
f[a_]:=Plus1terms/.s23->a;
(*out=f[0]*)
out=Simplify[(f[s23]-f[0])/s23+nonPlusterms]+Simplify[f[0]]*Plus1B[B,s23]
]


checknonedeltapole[xsc_,cntxsc_]:=Module[
{xsc1,xsc2},
xsc1=Expand[Select[Expand[xsc],FreeQ[#,\[Delta][s23]__]&&MatchQ[#,\[Epsilon]^-1__]&&FreeQ[#,\[Epsilon]^-2__]&]]//SimpPlus1;(*ignore warnings  which is from terms like 1/(1+1/s23) when setting s23\[Rule]0. The terms are 0 in the limit but mathematica gives warnings about them*)
xsc2=Expand[Select[Expand[cntxsc],FreeQ[#,\[Delta][s23]__]&&MatchQ[#,\[Epsilon]^-1__]&&FreeQ[#,\[Epsilon]^-2__]&]]/.Subscript[C, A]->3//SimpPlus1;
Expand[xsc1-xsc2]//Simplify
]


SimpPlus1PA[exp_](*PA stands for partially*):=Module[
{subplus,Plus1terms,nonPlusterms,out},
subplus={s23^(a_/;a>=1)*Plus1B[B,s23]->s23^(a-1),s23*Plus1B[B,s23]->1};
nonPlusterms=exp//Collect[#,Plus1B[B,s23]]&//Select[#,(*MatchQ[#,t3^1__]&&*)FreeQ[#,Plus1B[B,s23]^1__]&]&;
Plus1terms=Coefficient[exp,Plus1B[B,s23]]*Plus1B[B,s23]//Collect[#,{Log[__]},Simplify]&;
Plus1terms=Expand[Plus1terms]/.subplus//Collect[#,{Log[__]},Simplify]&;
Plus1terms=Plus1terms/.Log[\[Mu]]->Log[\[Mu]2]/2-Log[Pi]/2-Log[2]//Collect[#,{Log[__]},Simplify]&;
out=Plus1terms+nonPlusterms
]
