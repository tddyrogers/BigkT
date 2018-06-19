#!/apps/mathematica/bin/math -script
(*#!/home/nsato/apps/bin/math -script*)

<<share.m;
<<real.m;
<<virtual.m;
<<counter.m;

Channels={
{{"Agtoqqbg"},{"Agtoqqb"}},
{{"Aqtoqgg","Aqtoqqqb","AqtoqQQb"},{"Aqtoqg"}},
{{"Aqtogqg"},{"Aqtogq"}},
{{"Agtogqqb"},{"none"}},
{{"Aqtoqbqq"},{"none"}},
{{"AqtoQqQb"},{"none"}}};

(* Subscript[n, f]-1 from  identical contributions from Subscript[n, f]-1 
massless quarks QQb in AqtoqQQb *)

wt={
{{1},{1}},
{{1,1,Subscript[n, f]-1},{1}},
{{1},{1}},
{{1},{1}},
{{1},{1}},
{{1},{1}}};



CALC[ichannel_,type_]:=Module[
{out},
out= Which[
type == "real" ,Plus@@((real/@Channels[[ichannel]][[1]])*wt[[ichannel]][[1]]), 
type == "subtraction" ,Plus@@((subtraction/@Channels[[ichannel]][[1]])*wt[[ichannel]][[1]]),
type=="virtual",If[Channels[[ichannel]][[2]]=={"none"},{0,0,0},Plus@@virtual/@Channels[[ichannel]][[2]]], 
type=="counter",If[Channels[[ichannel]][[2]]=={"none"},0,Plus@@counter/@Channels[[ichannel]][[2]]]
]]


(*Pick a channel to calculate  from 1 to 6  *)
Do[
Print[ich];
realpiece=CALC[ich,"real"];
subtractionpiece=CALC[ich,"subtraction"];
virtualpiece=CALC[ich,"virtual"];
counterpiece=CALC[ich,"counter"];
Save["mdata/real-"<>ToString[ich]<> ".mx",realpiece];
Save["mdata/subtraction-"<>ToString[ich]<> ".mx",subtractionpiece];
Save["mdata/virtual-"<>ToString[ich]<> ".mx",virtualpiece];
Save["mdata/counter-"<>ToString[ich]<> ".mx",counterpiece],
{ich, 6}];






