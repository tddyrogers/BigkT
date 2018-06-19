(* ::Package:: *)

(* ::Text:: *)
(*General Formulars*)


sqampAgtoqqb[s_,t_]:=(-((8 g^2 (-1+\[Epsilon]) (-2 s t-2 t^2+Q^4 (-1+\[Epsilon])+s^2 (-1+\[Epsilon])-2 Q^2 (t-s \[Epsilon])))/(3 t (Q^2+s+t))))*3/8*1/(1-\[Epsilon]);(*note 3/8 accounts for the wrong average for gluon, 1/(1-\[Epsilon]) for wrong gluon spin avergage*)
sqampAgtoqbq[s_,t_]:=-((8 g^2 (-1+\[Epsilon]) (-2 s t-2 t^2+Q^4 (-1+\[Epsilon])+s^2 (-1+\[Epsilon])-2 Q^2 (t-s \[Epsilon])))/(3 t (Q^2+s+t)))*3/8*1/(1-\[Epsilon]);
sqampAqtoqg[s_,t_]:=(-((8 g^2 (-1+\[Epsilon]) (-2 Q^4-2 Q^2 (s+t)+s^2 (-1+\[Epsilon])+t^2 (-1+\[Epsilon])+2 s t \[Epsilon]))/(3 s t)));
sqampAqbtoqbg[s_,t_]:=(-((8 g^2 (-1+\[Epsilon]) (-2 Q^4-2 Q^2 (s+t)+s^2 (-1+\[Epsilon])+t^2 (-1+\[Epsilon])+2 s t \[Epsilon]))/(3 s t)));

sqampAqtogq[s_,t_]:=(8 g^2 (-1+\[Epsilon]) (-2 s^2-2 s t+Q^4 (-1+\[Epsilon])+t^2 (-1+\[Epsilon])-2 Q^2 (s-t \[Epsilon])))/(3 s (Q^2+s+t));

Rgtog[z_,M_]=(2*Subscript[C, A]*(Plus1[1-z]*z+(1-z)/z+z(1-z))+\[Delta][1-z]*(11Subscript[C, A]-4Subscript[n, f] Subscript[T, R])/6)*(-1/\[Epsilon]-Log[4*Pi]+EulerGamma(*+Log[M^2/\[Mu]^2]*));
Rgtoq[z_,M_]=Subscript[T, R](z^2+(1-z)^2)*(-1/\[Epsilon]-Log[4*Pi]+EulerGamma(*+Log[M^2/\[Mu]^2]*));
Rqtoq[z_,M_]=Subscript[C, F]*(Plus1[1-z]*(1+z^2)+3/2*\[Delta][1-z])*(-1/\[Epsilon]-Log[4*Pi]+EulerGamma(*+Log[M^2/\[Mu]^2]*));
Rqtog[z_,M_]=Subscript[C, F]*(1+(1-z)^2)/z*(-1/\[Epsilon]-Log[4*Pi]+EulerGamma(*+Log[M^2/\[Mu]^2]*));

(*Convention: e.g. Rp2gtogs2 means parton pdf with initial momentum p2, for the process p1+p2\[Rule]q1+q2*)
Rp2gtogs2[s2_,M_]=Rgtog[1-s2/(s+u+Q^2),M]/.{Plus1[s2/(s+u+Q^2)]->(s2-t)*(Plus1B[B,s2]+\[Delta][s2]*Log[B/-t]),\[Delta][s2/(s+u+Q^2)]->(-t)*\[Delta][s2],u->s2-Q^2-s-t}//Simplify;
Rq1gtogs2[s2_,M_]=Rgtog[1-s2/s,M]/.{Plus1[s2/s]->s*(Plus1B[B,s2]+\[Delta][s2]*Log[B/s]),\[Delta][s2/s]->(s)*\[Delta][s2],u->s2-Q^2-s-t}//Simplify;

Rp2gtoqs2[s2_,M_]=Rgtoq[1-s2/(s+u+Q^2),M]/.{Plus1[s2/(s+u+Q^2)]->(s2-t)*(Plus1B[B,s2]+\[Delta][s2]*Log[B/-t]),\[Delta][s2/(s+u+Q^2)]->(-t)*\[Delta][s2],u->s2-Q^2-s-t}//Simplify;
Rq1gtoqs2[s2_,M_]=Rgtoq[1-s2/s,M]/.{Plus1[s2/s]->s*(Plus1B[B,s2]+\[Delta][s2]*Log[B/s]),\[Delta][s2/s]->(s)*\[Delta][s2],u->s2-Q^2-s-t}//Simplify;

Rp2qtoqs2[s2_,M_]=Rqtoq[1-s2/(s+u+Q^2),M]/.{Plus1[s2/(s+u+Q^2)]->(s2-t)*(Plus1B[B,s2]+\[Delta][s2]*Log[B/-t]),\[Delta][s2/(s+u+Q^2)]->(-t)*\[Delta][s2],u->s2-Q^2-s-t}//Simplify;
Rq1qtoqs2[s2_,M_]=Rqtoq[1-s2/s,M]/.{Plus1[s2/s]->s*(Plus1B[B,s2]+\[Delta][s2]*Log[B/s]),\[Delta][s2/s]->(s)*\[Delta][s2],u->s2-Q^2-s-t}//Simplify;

Rp2qtogs2[s2_,M_]=Rqtog[1-s2/(s+u+Q^2),M]/.{Plus1[s2/(s+u+Q^2)]->(s2-t)*(Plus1B[B,s2]+\[Delta][s2]*Log[B/-t]),\[Delta][s2/(s+u+Q^2)]->(-t)*\[Delta][s2],u->s2-Q^2-s-t}//Simplify;
Rq1qtogs2[s2_,M_]=Rqtog[1-s2/s,M]/.{Plus1[s2/s]->s*(Plus1B[B,s2]+\[Delta][s2]*Log[B/s]),\[Delta][s2/s]->(s)*\[Delta][s2],u->s2-Q^2-s-t}//Simplify;



(* ::Text:: *)
(*Cntxscs for various Processes*)
(*(Note:*)
(*1.Exchanging particles q2 and q3 gives the same contribution, e.g. Ag2qqbg=Ag2qgqb. Do NOT exchange q2 and q3 to avoid double counting, since they are integrated over.*)
(*2.Exchanging q and qb gives the same contribution, e.g. Ag2qqbg=Ag2qbqg.*)
(* )*)


proccntxscAgtoqqbg[prefac_]:=Module[
{counterxsec,counterxsec1,counterxsec2,counterxsec3,bornxsecAgtoqqb,bornxsecAqtoqg},
bornxsecAgtoqqb[s_,t_]:=prefac[s,t]*sqampAgtoqqb[s,t];
bornxsecAqtoqg[s_,t_]:=prefac[s,t]*sqampAqtoqg[s,t];
counterxsec1=Series[-1/t*bornxsecAgtoqqb[-((s*t+Q^2*s23)/(s23-t)),t]*Rp2gtogs2[s23,M]*1/(8Pi^2),{\[Epsilon],0,0}]/.Subscript[T, R]->1/2/.Subscript[C, A]->3//Simplify//Normal;
counterxsec2=Series[-1/t*bornxsecAqtoqg[-((s*t+Q^2*s23)/(s23-t)),t]*Rp2gtoqs2[s23,M]*1/(8Pi^2),{\[Epsilon],0,0}]/.Subscript[T, R]->1/2/.w->1-y//Simplify//Normal;
counterxsec3=Series[-1/(s23-s)*bornxsecAgtoqqb[s,-((s*t+Q^2*s23)/(s23-s))]*Rq1qtoqs2[s23,M]*1/(8Pi^2),{\[Epsilon],0,0}]/.Subscript[T, R]->1/2/.Subscript[C, F]->4/3/.w->1-y//Simplify//Normal;
counterxsec=1*counterxsec1+1*counterxsec2+1*counterxsec3//Normal//Collect[#,\[Epsilon]]&
];(*pole term checked*)


proccntxscAgtoqbqg[prefac_]:=proccntxscAg2qqbg[prefac];


proccntxscAgtogqqb[prefac_]:=Module[
{counterxsec,counterxsec1,counterxsec2,bornxsecAgtoqqb,bornxsecAqtogq},
bornxsecAgtoqqb[s_,t_]:=prefac[s,t]*sqampAgtoqqb[s,t];
bornxsecAqtogq[s_,t_]:=prefac[s,t]*sqampAqtogq[s,t];
counterxsec1=Series[-1/t*bornxsecAqtogq[-((s*t+Q^2*s23)/(s23-t)),t]*Rp2gtoqs2[s23,M]*1/(8Pi^2),{\[Epsilon],0,0}]/.Subscript[T, R]->1/2/.w->1-y//Simplify//Normal;
counterxsec2=Series[-1/(s23-s)*bornxsecAgtoqqb[s,-((s*t+Q^2*s23)/(s23-s))]*Rq1qtogs2[s23,M]*1/(8Pi^2),{\[Epsilon],0,0}]/.Subscript[T, R]->1/2/.Subscript[C, F]->4/3/.w->1-y//Simplify//Normal;
counterxsec=2*counterxsec1+2*counterxsec2//Normal//Collect[#,\[Epsilon]]&
];(*pole term checked*)


proccntxscAqtoqgg[prefac_]:=Module[
{counterxsec,counterxsec1,counterxsec2,bornxsecAqtoqg},
bornxsecAqtoqg[s_,t_]:=prefac[s,t]*sqampAqtoqg[s,t];
counterxsec1=Series[-1/t*bornxsecAqtoqg[-((s*t+Q^2*s23)/(s23-t)),t]*Rp2qtoqs2[s23,M]*1/(8Pi^2),{\[Epsilon],0,0}]/.Subscript[T, R]->1/2/.Subscript[C, F]->4/3/.w->1-y//Simplify//Normal;
counterxsec2=Series[-1/(s23-s)*bornxsecAqtoqg[s,-((s*t+Q^2*s23)/(s23-s))]*Rq1qtoqs2[s23,M]*1/(8Pi^2),{\[Epsilon],0,0}]/.Subscript[T, R]->1/2/.Subscript[C, F]->4/3/.w->1-y//Simplify//Normal;
counterxsec=1*counterxsec1+1*counterxsec2//Normal//Collect[#,\[Epsilon]]&
];(*pole term checked*)


proccntxscAqbtoqbgg[prefac_]:=proccntxscAqtoqgg[prefac];


proccntxscAqtogqg[prefac_]:=Module[
{counterxsec,counterxsec1,counterxsec2,counterxsec3,bornxsecAqtoqg,bornxsecAqtogq},
bornxsecAqtoqg[s_,t_]:=prefac[s,t]*sqampAqtoqg[s,t];
bornxsecAqtogq[s_,t_]:=prefac[s,t]*sqampAqtogq[s,t];
counterxsec1=Series[-1/t*bornxsecAqtogq[-((s*t+Q^2*s23)/(s23-t)),t]*Rp2qtoqs2[s23,M]*1/(8Pi^2),{\[Epsilon],0,0}]/.Subscript[C, A]->3/.Subscript[T, R]->1/2/.Subscript[C, F]->4/3/.w->1-y//Simplify//Normal;
counterxsec2=Series[-1/(s23-s)*bornxsecAqtoqg[s,-((s*t+Q^2*s23)/(s23-s))]*Rq1qtogs2[s23,M]*1/(8Pi^2),{\[Epsilon],0,0}]/.Subscript[C, A]->3/.Subscript[T, R]->1/2/.Subscript[C, F]->4/3/.w->1-y//Simplify//Normal;
counterxsec3=Series[-1/(s23-s)*bornxsecAqtogq[s,-((s*t+Q^2*s23)/(s23-s))]*Rq1gtogs2[s23,M]*1/(8Pi^2),{\[Epsilon],0,0}]/.Subscript[C, A]->3/.Subscript[T, R]->1/2/.Subscript[C, F]->4/3/.w->1-y//Simplify//Normal;

counterxsec=1*counterxsec1+1*counterxsec2+1*counterxsec3//Normal//Collect[#,\[Epsilon]]&
];(*pole term checked*)


proccntxscAqbtogqbg[prefac_]:=proccntxscAqtogqg[prefac];


proccntxscAqtoqqqb[prefac_]:=Module[
{counterxsec,counterxsec1,counterxsec2,bornxsecAgtoqqb,bornxsecAqtogq},
bornxsecAgtoqqb[s_,t_]:=prefac[s,t]*sqampAgtoqqb[s,t];
bornxsecAqtogq[s_,t_]:=prefac[s,t]*sqampAqtogq[s,t];
counterxsec1=Series[-1/t*bornxsecAgtoqqb[-((s*t+Q^2*s23)/(s23-t)),t]*Rp2qtogs2[s23,M]*1/(8Pi^2),{\[Epsilon],0,0}]/.Subscript[T, R]->1/2/.Subscript[C, F]->4/3/.w->1-y//Simplify//Normal;
counterxsec2=Series[-1/(s23-s)*bornxsecAqtogq[s,-((s*t+Q^2*s23)/(s23-s))]*Rq1gtoqs2[s23,M]*1/(8Pi^2),{\[Epsilon],0,0}]/.Subscript[T, R]->1/2/.Subscript[C, F]->4/3/.w->1-y//Simplify//Normal;
counterxsec=1*counterxsec1+1*counterxsec2//Normal//Collect[#,\[Epsilon]]&
];(*pole term checked*)


proccntxscAqbtoqbqbq[prefac_]:=proccntxscAqtoqqqb[prefac];


proccntxscAqtoqbqq[prefac_]:=proccntxscAqtoqqqb[prefac];(*The cross sections are not the same*)


proccntxscAqbtoqqbqb[prefac_]:=proccntxscAqtoqbqq[prefac];


proccntxscAqtoqQQb[prefac_]:=0;(*q and Q are two different quark flavors*)


proccntxscAqbtoqbQQb[prefac_]:=0;


proccntxscAqtoQqQb[prefac_]:=Module[
{counterxsec,counterxsec1,counterxsec2,bornxsecAgtoqqb,bornxsecAqtogq},
bornxsecAgtoqqb[s_,t_]:=prefac[s,t]*sqampAgtoqqb[s,t];
bornxsecAqtogq[s_,t_]:=prefac[s,t]*sqampAqtogq[s,t];
counterxsec1=Series[-1/t*(gp/g)^2*bornxsecAgtoqqb[-((s*t+Q^2*s23)/(s23-t)),t]*Rp2qtogs2[s23,M]*1/(8Pi^2),{\[Epsilon],0,0}]/.Subscript[T, R]->1/2/.Subscript[C, F]->4/3/.w->1-y//Simplify//Normal;
(*g and gp are electric charges for q and Q respectively*)
counterxsec2=Series[-1/(s23-s)*bornxsecAqtogq[s,-((s*t+Q^2*s23)/(s23-s))]*Rq1gtoqs2[s23,M]*1/(8Pi^2),{\[Epsilon],0,0}]/.Subscript[T, R]->1/2/.Subscript[C, F]->4/3/.w->1-y//Simplify//Normal;
counterxsec=1*counterxsec1+1*counterxsec2//Normal//Collect[#,\[Epsilon]]&
];(*pole term checked*)


proccntxscAqbtoQqbQb[prefac_]:=proccntxscAqtoQqQb[prefac];


proccntxscAqtoQbqQ[prefac_]:=proccntxscAqtoQqQb[prefac];


proccntxscAqbtoQbqbQ[prefac_]:=proccntxscAqtoQqQb[prefac];
