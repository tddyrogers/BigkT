id U(sc?,phi?,p?,s?,M?)*UB(sc?,phi?,p?,s?,M?) = gi_(sc,p)+gi_(sc)*M;
id V(sc?,phi?,p?,s?,M?)*VB(sc?,phi?,p?,s?,M?) = gi_(sc,p)-gi_(sc)*M;
id eps(phi?,p1?,s?,M?,o1?)*epsB(phi?,p1?,s?,M?,o2?)=d_(o1,o2);
id eps(phi?,q?,s?,M?,o1?,j1?)*epsB(phi?,q?,s?,M?,o2?,j2?)=-d_(o1,o2)*d_(j1,j2);
id den(p?,M?)=den(p.p-M*M);
id X(j1?)*XB(j2?)=d_(j1,j2);
Trace4,sc1;
Trace4,sc2;
Trace4,sc3;
Trace4,sc4;
id MA=0;
id p1.p1=-Q2;
id p1.p2=(s+Q2)/2;
id p1.q1=-(t1+Q2)/2;
id p1.q2=-(t2+Q2)/2;
id p1.q3=-(t3+Q2)/2;
id p2.p2=0;
id p2.q1=-u1/2;
id p2.q2=-u2/2;
id p2.q3=-u3/2;
id q1.q1=0;
id q1.q2=s12/2;
id q1.q3=s13/2;
id q2.q2=0;
id q2.q3=s23/2;
id q3.q3=0;
argument;
id MA=0;
id p1.p1=-Q2;
id p1.p2=(s+Q2)/2;
id p1.q1=-(t1+Q2)/2;
id p1.q2=-(t2+Q2)/2;
id p1.q3=-(t3+Q2)/2;
id p2.p2=0;
id p2.q1=-u1/2;
id p2.q2=-u2/2;
id p2.q3=-u3/2;
id q1.q1=0;
id q1.q2=s12/2;
id q1.q3=s13/2;
id q2.q2=0;
id q2.q3=s23/2;
id q3.q3=0;
endargument;
id den(x?)=1/x;
Multiply replace_ (Mu,0);
Multiply replace_ (MG,0);
#call SUn;
id NA=NF^2-1;
id NF=3;
id NF^-1=1/3;
*id D=4;
id a=1/2;
id nf=1;
*id ave=1/2/3;
*id ave=1/2/3;
