* gzvert_his5_qud5.F
* form factors for {MNE1, MNE1} -> {0, MZ}
* generated by FormCalc 12 Apr 2013 17:41

	subroutine gzvert_his5_qud5
	implicit none

#include "gzvars.h"

	double complex dup1571

        dup1571 = -(1/2.D0*(AAABR(1054)*MTR121(qud5,his5))) - 
     -   1/2.D0*(AAABR(1054)*MTR122(qud5,his5)) + 
     -   AAABR(1052)*(MTR121(qud5,his5) + MTR122(qud5,his5))

        Cloop(1) = Cloop(1) + 
     -   (B0(S,qudMass(qud5)**2,qudMass(qud5)**2)*
     -       (-(1/(96.D0*Pi**2)*
     -            (Abb5*dup1571*EE*MTR183(1,1,his5)*qudMass(qud5)))
     -           + 1/(96.D0*Pi**2)*
     -          (Abb7*dup1571*EE*MTR184(1,1,his5)*qudMass(qud5)))+
     -        Cval(cc00,iint80(qud5))*
     -       (1/(24.D0*Pi**2)*
     -          (Abb5*dup1571*EE*MTR183(1,1,his5)*qudMass(qud5)) - 
     -         1/(24.D0*Pi**2)*
     -          (Abb7*dup1571*EE*MTR184(1,1,his5)*qudMass(qud5)))+
     -        Cval(cc12,iint80(qud5))*
     -       (-(1/(24.D0*Pi**2)*
     -            (AbbSum30*dup1571*EE*MTR183(1,1,his5)*
     -              qudMass(qud5))) + 
     -         1/(24.D0*Pi**2)*
     -          (AbbSum31*dup1571*EE*MTR184(1,1,his5)*
     -            qudMass(qud5))) + 
     -      Cval(cc0,iint80(qud5))*
     -       (-(1/(96.D0*Pi**2)*
     -            (Abb5*dup1571*EE*MNE1**2*MTR183(1,1,his5)*
     -              qudMass(qud5))) + 
     -         1/(192.D0*Pi**2)*
     -          (Abb5*dup1571*EE*(T + U)*MTR183(1,1,his5)*
     -            qudMass(qud5)) + 
     -         1/(96.D0*Pi**2)*
     -          (AbbSum21*EE*(AAABR(1052) - 1/2.D0*AAABR(1054))*
     -            MTR121(qud5,his5)*MTR183(1,1,his5)*qudMass(qud5))
     -           + 1/(96.D0*Pi**2)*
     -          (AbbSum20*EE*(AAABR(1052) - 1/2.D0*AAABR(1054))*
     -            MTR122(qud5,his5)*MTR183(1,1,his5)*qudMass(qud5))
     -           + 1/(96.D0*Pi**2)*
     -          (Abb7*dup1571*EE*MNE1**2*MTR184(1,1,his5)*
     -            qudMass(qud5)) - 
     -         1/(192.D0*Pi**2)*
     -          (Abb7*dup1571*EE*(T + U)*MTR184(1,1,his5)*
     -            qudMass(qud5)) - 
     -         1/(96.D0*Pi**2)*
     -          (AbbSum29*EE*(AAABR(1052) - 1/2.D0*AAABR(1054))*
     -            MTR121(qud5,his5)*MTR184(1,1,his5)*qudMass(qud5))
     -           - 1/(96.D0*Pi**2)*
     -          (AbbSum28*EE*(AAABR(1052) - 1/2.D0*AAABR(1054))*
     -            MTR122(qud5,his5)*MTR184(1,1,his5)*qudMass(qud5))
     -         ))/(S - hisMass(his5)**2+(0D0,1D0)*hisMass(his5)*hisW(his5))

#ifdef DEBUG
	DEB 'vert_his5_qud5: Cloop(1) =', Cloop(1)
#endif

	end

