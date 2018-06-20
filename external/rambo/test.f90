program main
implicit none
integer :: LW
double precision :: ET
double precision :: XM(3)
double precision :: P(4,3)
double precision :: WT
ET=10d0
XM(1)=0d0
XM(2)=0d0
XM(3)=0d0
LW=1

call RAMBO(3,ET,XM,P,WT,LW)
print*,P 

end program 


