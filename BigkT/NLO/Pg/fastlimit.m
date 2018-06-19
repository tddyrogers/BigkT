FastLimit[term_] := Module[
   {tablein = (term // Expand // Apply[List, #] &), tableout = {}, 
    bad = {}, limit = 0, limit2 = 0, i = 0, finaloutput = 0},
   For[i = 1, i <= (tablein // Length), i++,
    limit = 
     Quiet[Check[(tablein[[i]] // Together // 
         ReplaceAll[#, s23 -> 0] &), BAD]];
    If[limit === BAD,
     (*Print["bad at "<>ToString[i]<>"..."<>ToString[limit]];*)
     
     bad = Append[bad, tablein[[i]]],
     tableout = Append[tableout, limit]]
    ];
   tableout = tableout // DeleteCases[#, 0] &;
   limit2 = 
    Check[(bad // Apply[Plus, #] &) // Limit[#, s23 -> 0] &, BAD];
   Print["Test limit: " <> ToString[(limit2 // FullForm)]];
   finaloutput = (tableout // Apply[Plus, #] &) + limit2;
   finaloutput
   ];
      
   
Simple1[expr_] := 
  expr // Expand // 
   Collect[#, {Log[x__] Log[x__], Log[x__] Log[y__], Log[__], 
      PolyLog[__, __]}, Simplify] &;

