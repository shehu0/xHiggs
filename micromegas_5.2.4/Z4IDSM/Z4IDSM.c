/*====== Modules ===============
   Keys to switch on 
   various modules of micrOMEGAs  
================================*/
      
#define MASSES_INFO      
  /* Display information about mass spectrum  */

//#define CONSTRAINTS
//#define HIGGSBOUNDS
//#define HIGGSSIGNALS
//#define LILITH
//#define SMODELS

//#define OMEGA     /* Calculate relic density and display contribution of  individual channels */
//#define FREEZEIN  /*  Calculate relic density in Freeze-in scenario  */


//#define INDIRECT_DETECTION  
  /* Compute spectra of gamma/positron/antiprotons/neutrinos for DM annihilation; 
     Calculate <sigma*v>;
     Integrate gamma signal over DM galactic squared density for given line 
     of sight; 
     Calculate galactic propagation of positrons and antiprotons.      
  */
      
/*#define RESET_FORMFACTORS*/
  /* Modify default nucleus form factors, 
    DM velocity distribution,
    A-dependence of Fermi-dencity
  */     
//#define CDM_NUCLEON     
  /* Calculate amplitudes and cross-sections for  CDM-mucleon collisions */  

//#define CDM_NUCLEUS 
   //Calculate  exclusion rate for direct detection experiments Xenon1T and DarkSide50
//#define CROSS_SECTIONS 
//#define DECAYS   
#define SMODELS       
  
/*===== end of Modules  ======*/

/*===== Options ========*/
//#define SHOWPLOTS
     /* Display  graphical plots on the screen */ 

//#define CLEAN    to clean intermediate files

/*===== End of DEFINE  settings ===== */


#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"
#include"lib/pmodel.h"


int main(int argc,char** argv)
{  int err;
   char cdmName[10];
   int spin2, charge3,cdim;
   FILE *f1 = fopen("small-splitting-nearlydegen-csHCHC-ctauHCHC-general-scan-with-DD","w"); 
   FILE *f2 = fopen("Z4IDSM_2_WIMPs","r"); 
   //FILE *f2 = fopen("data1","r");    
  ForceUG=1;  /* to Force Unitary Gauge assign 1 */
  //useSLHAwidth=0;
  VZdecay=0; VWdecay=0;  

  if(argc==1)
  { 
      printf(" Correct usage:  ./main  <file with parameters> \n");
      printf("Example: ./main data1.par\n");
      exit(1);
  }
                               
  err=readVar(argv[1]);
//  toFeebleList("~sc");
  if(err==-1)     {printf("Can not open the file\n"); exit(1);}
  else if(err>0)  { printf("Wrong file contents at line %d\n",err);exit(1);}
           
  int i,j,k,l,m,n,o;
  double Omega, pval, sigmaV, width1, width2, ctau1;
  double csSIp,csSIn,csSDp,csSDn, csSIp_,csSIn_,csSDp_,csSDn_;  
  double MA0, MH0, MS, MHC, la3, la4, la5, laS2, csSIp1, omega1, omega2;
  double csAH, csHCHC, csAHC, csHHC;
  double Width1, Width2;
  double laS12, laS21, laS1;
  double ctau2;
/*  double la4, la5;
  la5=1e-10;
  la4=-1e-10;
  double MA0, MH0, MHC;
  MA0=120;
  assignValW("MH3", MA0);
  MH0=Sqrt(pow(MA0,2)+la5*pow(246,2));
  assignValW("MHX", MH0);
  MHC=Sqrt((pow(MH0,2)+ pow(MA0,2) - la4*pow(246,2))/2);
  assignValW("MHC", MHC);*/
  /*laS2 = 1e-01;
  assignValW("laS2",laS2);
  for(i=0;i<45;i++)
      {
  laS2=laS2+2e-02;
  assignValW("laS2",laS2);*/
/*  double laS12 = 1e-01; 
  double laS21 = -laS12;
  assignValW("laS12",laS12);
  assignValW("laS21",laS21);
  for(i=0;i<45;i++)
      {
  laS12=laS12+2e-02;
  assignValW("laS12",laS12);
  laS21=-laS12;
  assignValW("laS21",laS21);*/

//  while( fscanf(f2,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf  %lf  %lf ",&MA0,&MH0,&MS,&MHC,&la3,&la4,&la5,&laS2,&csSIp,&csSIp1,&omega1,&omega2,&Omega,&pval,&sigmaV,&Width1,&Width2) !=EOF)
// { 
//  while( fscanf(f2,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf  %lf  %lf",&MA0,&MH0,&MS,&MHC,&la3,&laS1,&laS2,&laS12,&laS21,&csSIp,&csSIp1,&omega1,&omega2,&Omega,&pval,&sigmaV) !=EOF)
// { 
//      if(MA0<MH0 && MA0<2*MS)
//{
//      if(MHC-MA0<5 && MH0-MA0<5)
//{
//      assignValW("la3", la3);   
//      assignValW("MH3", MA0);
//      assignValW("MHX", MH0);
//      assignValW("Msc", MS);
//      assignValW("MHC", MHC);
//      assignValW("laS2", laS2);
//     assignValW("laS1", laS1);
//     assignValW("laS12", laS12);
//     assignValW("laS21", laS21);

  err=sortOddParticles(cdmName);
  if(err) { printf("Can't calculate %s\n",cdmName); return 1;}
  
  if(CDM1) 
  { 
     qNumbers(CDM1, &spin2, &charge3, &cdim);
     printf("\nDark matter candidate is '%s' with spin=%d/2 mass=%.2E\n",CDM1,  spin2,Mcdm1); 
     if(charge3) printf("Dark Matter has electric charge %d/3\n",charge3);
     if(cdim!=1) printf("Dark Matter is a color particle\n");
  }
  if(CDM2) 
  { 
     qNumbers(CDM2, &spin2, &charge3, &cdim);
     printf("\nDark matter candidate is '%s' with spin=%d/2 mass=%.2E\n",CDM2,spin2,Mcdm2); 
     if(charge3) printf("Dark Matter has electric charge %d/3\n",charge3);
     if(cdim!=1) printf("Dark Matter is a color particle\n");
  }

  
#ifdef MASSES_INFO
{
  printf("\n=== MASSES OF HIGGS AND ODD PARTICLES: ===\n");
  printHiggs(stdout);
  printMasses(stdout,1);
}
#endif

#ifdef CONSTRAINTS
{ double csLim;
  if(Zinvisible()) printf("Excluded by Z->invisible\n");
  if(LspNlsp_LEP(&csLim)) printf("LEP excluded by e+,e- -> DM q qbar Cross Section= %.2E pb\n",csLim);
}
#endif


#if defined(HIGGSBOUNDS) || defined(HIGGSSIGNALS)
{  int NH0=3, NHch=1; // number of neutral and charged Higgs particles.
   int HB_id[3]={0,0,0},HB_result[3];
   double  HB_obsratio[3],HS_observ=-1,HS_chi2, HS_pval;
   char HB_chan[3][100]={""}, HB_version[50], HS_version[50]; 
   NH0=hbBlocksMO("HB.in",&NHch); 
//    NH0= hbBlocksMDL("HB.in",&NHch); 
   system("echo 'BLOCK DMASS\n 25  2  '>> HB.in");
#include "../include/hBandS.inc"
#ifdef HIGGSBOUNDS
   printf("HiggsBounds(%s)\n", HB_version);
   for(int i=0;i<3;i++) if(HB_id[i]) printf("  id= %d  result = %d  obsratio=%.2E  channel= %s \n", HB_id[i],HB_result[i],HB_obsratio[i],HB_chan[i]);
#endif 
#ifdef HIGGSSIGNALS
   if(HS_observ>=0)
   {
     printf("HiggsSignals(%s)\n",HS_version); 
     printf("  Nobservables=%.0f chi^2 = %.2E pval= %.2E\n",HS_observ,HS_chi2, HS_pval);
   }
#endif   
}
#endif

#ifdef LILITH
{  double m2logL, m2logL_reference=0,pvalue;
   int exp_ndf,n_par=0,ndf;
   char call_lilith[100], Lilith_version[20];

   if(LilithMDL("Lilith_in.xml"))
   {        
#include "../include/Lilith.inc"
      if(ndf)
      {
        printf("LILITH(DB%s):  -2*log(L): %.2f; -2*log(L_reference): %.2f; ndf: %d; p-value: %.2E \n", 
        Lilith_version,m2logL,m2logL_reference,ndf,pvalue);
      }  
   } else printf("LILITH: there is no Higgs candidate\n");
}   
#endif


#ifdef OMEGA
{ int fast=1;
  double Beps=1.E-4, cut=0.01;
  double Xf;  
  int i,err; 
  double Lmin,Lmax;
  printf("\n==== Calculation of relic density =====\n");   

//ExcludedFor2DM="1110 2220 1120  1210 1220 2210  1122 1112 1222";
//  ExcludedFor2DM="1122";
  Omega= darkOmega2(fast,Beps);


 /* displayPlot("vs11","T",Tend,Tstart,0,5
      ,"vs1100",0,vs1100F,NULL
      ,"vs1120",0,vs1120F,NULL
      ,"vs1122",0,vs1122F,NULL
      ,"vs1110",0,vs1110F,NULL
      ,"vs1112",0,vs1112F,NULL
      );*/
      
 /* displayPlot("vs12","T",Tend,Tstart,0,4
             ,"vs1210",0,vs1210F,NULL
             ,"vs1222",0,vs1222F,NULL
             ,"vs1120",0,vs1120F,NULL
             ,"vs1211",0,vs1211F,NULL);*/
                                
/*  displayPlot("vs22","T",Tend,Tstart,0,5
             ,"vs2200",0,vs2200F,NULL
             ,"vs2210",0,vs2210F,NULL
             ,"vs2211",0,vs2211F,NULL
             ,"vs2220",0,vs2220F,NULL
             ,"vs2221",0,vs2221F,NULL
             );*/


 // displayPlot("Y","T",   Tend,Tstart,0,4,"Y1" ,0,Y1F,NULL,"Y1eq",0,Yeq1,NULL,"Y2eq",0,Yeq2,NULL,"Y2",0,Y2F,NULL);

 
  // printf("Tend=%.3E\n",Tend);
  // printf("Tstart=%.3E\n",Tstart);  
/*   double j;                             
   displayPlot("Y","T",   Tend,Tstart,0,2,"Y1" ,0,Y1F,NULL,"Yeq1",0,Yeq1,NULL);     
   displayPlot("Y","T",   Tend,Tstart,0,2,"Y2" ,0,Y2F,NULL,"Yeq2",0,Yeq2,NULL);*/
/*  for(j=4.42; j<4.43; j=j+0.001)
  {
      printf("T=%.3E\n",j);
      printf("Y2F/Yeq2=%.3E\n", Y2F(j)/Yeq2(j));

  }*/
  printf("omega1=%.2E\n", Omega*(1-fracCDM2));
  printf("omega2=%.2E\n", Omega*fracCDM2);
  printf("fracCDM2=%E\n", fracCDM2);
 // fracCDM2=1;
 // printf("fracCDM2=%E\n", fracCDM2);
}

#endif

#ifdef FREEZEIN
{
  double TR=1E10;
  double omegaFi;  
  toFeebleList(CDM1);
  VWdecay=0; VZdecay=0;
  
  omegaFi=darkOmegaFi(TR,&err);
  printf("omega freeze-in=%.3E\n", omegaFi);
  printf("   omega1=%.3E omega2= %.3E fracCDM2=%.3E\n",omegaFi*(1-fracCDM2), omegaFi*fracCDM2,fracCDM2); 
  printChannelsFi(0,0,stdout);
}
#endif



#ifdef INDIRECT_DETECTION
{ 
  int err,i;
  double Emin=1,/* Energy cut  in GeV   */  sigmaV;
  double vcs_gz,vcs_gg;
  char txt[100];
  double SpA[NZ],SpE[NZ],SpP[NZ];
  double FluxA[NZ],FluxE[NZ],FluxP[NZ];
  double * SpNe=NULL,*SpNm=NULL,*SpNl=NULL;
  double Etest=Mcdm/2;
  
printf("\n==== Indirect detection =======\n");  

  sigmaV=calcSpectrum(4,SpA,SpE,SpP,SpNe,SpNm,SpNl ,&err);
    /* Returns sigma*v in cm^3/sec.     SpX - calculated spectra of annihilation.
       Use SpectdNdE(E, SpX) to calculate energy distribution in  1/GeV units.
       
       First parameter 1-includes W/Z polarization
                       2-includes gammas for 2->2+gamma
                       4-print cross sections             
    */
  printf("sigmav=%.2E[cm^3/s]\n",sigmaV);  


  { 
     double fi=0.1,dfi=0.05; /* angle of sight and 1/2 of cone angle in [rad] */ 

     gammaFluxTab(fi,dfi, sigmaV, SpA,  FluxA);     
     printf("Photon flux  for angle of sight f=%.2f[rad]\n"
     "and spherical region described by cone with angle %.2f[rad]\n",fi,2*dfi);
#ifdef SHOWPLOTS
     sprintf(txt,"Photon flux[cm^2 s GeV]^{1} at f=%.2f[rad], cone angle %.2f[rad]",fi,2*dfi);
     displayPlot(txt,"E[GeV]",Emin,Mcdm,0,1,"flux",0,SpectdNdE,FluxA);
#endif
     printf("Photon flux = %.2E[cm^2 s GeV]^{-1} for E=%.1f[GeV]\n",SpectdNdE(Etest, FluxA), Etest);       
  }

  { 
    posiFluxTab(Emin, sigmaV, SpE,  FluxE);
#ifdef SHOWPLOTS     
    displayPlot("positron flux [cm^2 s sr GeV]^{-1}","E[GeV]",Emin,Mcdm,0,1,"flux",0,SpectdNdE,FluxE);
#endif
    printf("Positron flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxE),  Etest);           
  }
  
  { 
    pbarFluxTab(Emin, sigmaV, SpP,  FluxP  ); 
#ifdef SHOWPLOTS    
     displayPlot("antiproton flux [cm^2 s sr GeV]^{-1}" ,"E[GeV]",Emin, Mcdm,0,1,"flux",0,SpectdNdE,FluxP);
#endif
    printf("Antiproton flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxP),  Etest);             
  }
}  
#endif

#ifdef RESET_FORMFACTORS
{
/* 
   The user has approach to form factors  which specifies quark contents 
   of  proton and nucleon via global parametes like
      <Type>FF<Nucleon><q>
   where <Type> can be "Scalar", "pVector", and "Sigma"; 
         <Nucleon>     "P" or "N" for proton and neutron
         <q>            "d", "u","s"

   calcScalarQuarkFF( Mu/Md, Ms/Md, sigmaPiN[MeV], sigmaS[MeV])  
   calculates and rewrites Scalar form factors
*/
  printf("\n======== RESET_FORMFACTORS ======\n");
 
  printf("protonFF (default) d %.2E, u %.2E, s %.2E\n",ScalarFFPd, ScalarFFPu,ScalarFFPs);                               
  printf("neutronFF(default) d %.2E, u %.2E, s %.2E\n",ScalarFFNd, ScalarFFNu,ScalarFFNs);

//                    To restore default form factors of  version 2  call 
     calcScalarQuarkFF(0.553,18.9,55.,243.5);

  printf("protonFF (new)     d %.2E, u %.2E, s %.2E\n",ScalarFFPd, ScalarFFPu,ScalarFFPs);                               
  printf("neutronFF(new)     d %.2E, u %.2E, s %.2E\n",ScalarFFNd, ScalarFFNu,ScalarFFNs);

//                    To restore default form factors  current version  call 
//  calcScalarQuarkFF(0.56,20.2,34,42);



/* Option to change parameters of DM velocity  distribution  */   
   SetfMaxwell(220.,600.);
/* 
    dN  ~  exp(-v^2/arg1^2)*Theta(v-arg2)  d^3v     
    Earth velocity with respect to Galaxy defined by 'Vearth' parameter.
    All parameters are  in [km/s] units.       
*/


}
#endif

#ifdef DECAYS
{ char*  pname1 = "~~X";
  txtList L1;

  if(pname1)
  {
    width1=pWidth(pname1,&L1);
    printf("\n%s :   total width=%E \n and Branchings:\n",pname1,width1);
    printTxtList(L1,stdout);
    ctau1=197e-18/width1;
  }
 char*  pname2 = "~~H+";
  txtList L2;

  if(pname2)
  {
    width2=pWidth(pname2,&L2);
    printf("\n%s :   total width=%E \n and Branchings:\n",pname2,width2);
    printTxtList(L2,stdout);
    ctau2=197e-18/width2;
  }

}
#endif



#ifdef CDM_NUCLEON
{ double pA0[2],pA5[2],nA0[2],nA5[2];
  double Nmass=0.939; /*nucleon mass*/
  double SCcoeff;        
//  double csSIp,csSIn,csSDp,csSDn, csSIp_,csSIn_,csSDp_,csSDn_;
  
printf("\n==== Calculation of CDM-nucleons amplitudes  =====\n");   

    nucleonAmplitudes(CDM1, pA0,pA5,nA0,nA5);
    printf("%s[%s]-nucleon micrOMEGAs amplitudes\n",CDM1,aCDM1? aCDM1:CDM1 );
    printf("proton:  SI  %.3E [%.3E]  SD  %.3E [%.3E]\n",pA0[0], pA0[1],  pA5[0], pA5[1] );
    printf("neutron: SI  %.3E [%.3E]  SD  %.3E [%.3E]\n",nA0[0], nA0[1],  nA5[0], nA5[1] ); 

    SCcoeff=4/M_PI*3.8937966E8*pow(Nmass*Mcdm1/(Nmass+ Mcdm1),2.);
    csSIp=  SCcoeff*pA0[0]*pA0[0];  csSIp_=  SCcoeff*pA0[1]*pA0[1];
    csSDp=3*SCcoeff*pA5[0]*pA5[0];  csSDp_=3*SCcoeff*pA5[1]*pA5[1];
    csSIn=  SCcoeff*nA0[0]*nA0[0];  csSIn_=  SCcoeff*nA0[1]*nA0[1];
    csSDn=3*SCcoeff*nA5[0]*nA5[0];  csSDn_=3*SCcoeff*nA5[1]*nA5[1];
    
    printf("%s[%s]-nucleon cross sections[pb] :\n",CDM1,aCDM1);
    printf(" proton  SI %.3E [%.3E] SD %.3E [%.3E]\n", csSIp,csSIp_,csSDp,csSDp_);
    printf(" neutron SI %.3E [%.3E] SD %.3E [%.3E]\n", csSIn,csSIn_,csSDn,csSDn_); 
                
    nucleonAmplitudes(CDM2, pA0,pA5,nA0,nA5);
    printf("%s[%s]-nucleon micrOMEGAs amplitudes\n",CDM2,aCDM2?aCDM2:CDM2);
    printf("proton:  SI  %.3E [%.3E]  SD  %.3E [%.3E]\n",pA0[0], pA0[1],  pA5[0], pA5[1] );
    printf("neutron: SI  %.3E [%.3E]  SD  %.3E [%.3E]\n",nA0[0], nA0[1],  nA5[0], nA5[1] ); 

    SCcoeff=4/M_PI*3.8937966E8*pow(Nmass*Mcdm2/(Nmass+ Mcdm2),2.);
    csSIp=  SCcoeff*pA0[0]*pA0[0];  csSIp_=  SCcoeff*pA0[1]*pA0[1];
    csSDp=3*SCcoeff*pA5[0]*pA5[0];  csSDp_=3*SCcoeff*pA5[1]*pA5[1];
    csSIn=  SCcoeff*nA0[0]*nA0[0];  csSIn_=  SCcoeff*nA0[1]*nA0[1];
    csSDn=3*SCcoeff*nA5[0]*nA5[0];  csSDn_=3*SCcoeff*nA5[1]*nA5[1];
                    
    printf("%s[%s]-nucleon cross sections[pb]:\n",CDM2,aCDM2);
    printf(" proton  SI %.3E [%.3E] SD %.3E [%.3E]\n", csSIp,csSIp_,csSDp,csSDp_);
    printf(" neutron SI %.3E [%.3E] SD %.3E [%.3E]\n", csSIn,csSIn_,csSDn,csSDn_);             
}
#endif
  
#ifdef CDM_NUCLEUS
{ char* expName; 
  printf("\n===== Direct detection exclusion:======\n");
  pval=DD_pval(AllDDexp, Maxwell, &expName);
  if(pval<0.1 )  printf("Excluded by %s  %.1f%%\n", expName, 100*(1-pval)); 
  else printf("Not excluded by DD experiments  at 90%% level \n");  
  //fprintf(f1, "%e %e %e %e %e\n", laS2, Omega*(1-fracCDM2), Omega*fracCDM2, Omega, pval);  
  //fprintf(f1, "%e %e %e %e %e\n", laS12, Omega*(1-fracCDM2), Omega*fracCDM2, Omega, pval);  
}

#endif 
//}
#ifdef CROSS_SECTIONS
{
  char* pName1= "~~X";
  char* pName2= "~~H3";

     double csAH, Pcm=6500, Qren, Qfact, pTmin=-1;
     int nf=5;

     Qren=Qfact= MA0+MH0; 
 
     printf("\npp > nextOdd  at sqrt(s)=%.2E GeV\n",2*Pcm);  
  
/*     csAH=hCollider(Pcm,1,nf,Qren, Qfact, pName1,pName2,pTmin,1);
     printf("Production of odd particle: cs(pp-> %s,%s)=%.2E[pb]\n",pName1,pName2, csAH);*/
  MHC=findValW("MHC");
  char*pName3 = "~~H+";
  char*pName4 = "~~H-"; 
  printf("MHC=%3E", MHC);
     double csHCHC, Tmin=-1;

     Qren=Qfact==2*MHC; 
 
     printf("\npp > nextOdd  at sqrt(s)=%.2E GeV\n",2*Pcm);  
  
     csHCHC=hCollider(Pcm,1,nf,Qren, Qfact, pName3,pName4,Tmin,1);
     printf("Production of 'pName3' and 'pName4': cs(pp-> %s,%s)=%.2E[pb]\n",pName3,pName4, csHCHC);
    
/*  char* pName5 = "~~X";
  char* pName6 = "~~H-";
     double csHHC;

  
     Qren=Qfact==MHC+MH0; 
 
     printf("\npp > nextOdd  at sqrt(s)=%.2E GeV\n",2*Pcm);  

     csHHC=hCollider(Pcm,1,nf,Qren, Qfact, pName5,pName6,Tmin,1);
     printf("Production of 'pName5' and 'pName6': cs(pp-> %s,%s)=%.2E[pb]\n",pName5,pName6, csHHC);
  
  char* pName7 = "~~H3";
  char* pName8 = "~~H-";

     double csAHC;

  
     Qren=Qfact=MHC+MA0; 
 
     printf("\npp > nextOdd  at sqrt(s)=%.2E GeV\n",2*Pcm);  
  
     csAHC=hCollider(Pcm,1,nf,Qren, Qfact, pName7,pName8,Tmin,1);
     printf("Production of 'pName7' and 'pName8': cs(pp-> %s,%s)=%.2E[pb]\n",pName7,pName8, csAHC);*/
     


if(pval>0.1)
{
  fprintf(f1, "%e %e %e %e %e %e %e %e %e %e %e %e %e\n", MA0, MH0, MHC, MS, la3, laS1, laS2, laS12, laS21, width2, ctau2, csHCHC, pval);
}
}
#endif 

#ifdef SMODELS
{  int result=0;
   double Rvalue=0;
   char analysis[30]={},topology[30]={};
   int LHCrun=LHC8|LHC13;  //  LHC8  - 8TeV; LHC13  - 13TeV; 
#include "../include/SMODELS.inc" 
}   
#endif 
//}
//}
//}      
#ifdef CLEAN
  system("rm -f HB.* HS.* hb.* hs.*  debug_channels.txt debug_predratio.txt  Key.dat");
  system("rm -f Lilith_*   particles.py*");
  system("rm -f  smodels.in  smodels.log  smodels.out  summary.*");  
#endif 


  killPlots();
  return 0;
}
