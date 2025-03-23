HumanPBPK.code <- '
$PARAM @annotated

//Physiological parameters		

BW                  :   63    : kg,                  GBZ/T200.1-2007
QCC                 :   15    : L/h/kg^0.75,         GBZ/T200.3-2014
QLC                 :   0.255 : Unitless,            GBZ/T200.3-2014
QKC                 :   0.19	: Unitless,            GBZ/T200.3-2014
QLuC                :   1   	: Unitless,            Pulmonary circulation equals systemic circulation
QFC                 :   0.05	: Unitless,            GBZ/T200.3-2014
Htc                 :   0.44  : Unitless,            Hematocrit for human; ICRP Publication 89 (2003)
VLC                 :   0.022 : Unitless,            GBZ/T200.2-2007
VKC                 :   0.0046: Unitless,            GBZ/T200.2-2007
VLuC                :   0.02  : Unitless,            GBZ/T200.2-2007
VFC                 :   0.143 : Unitless,            GBZ/T200.2-2007
VPlasC              :   0.0312: L/kg BW,             GBZ/T200.2-2007
VfilC               :  0.00046: L/kg BW,             Fraction vol. of filtrate; 10% of Kidney volume; (Worley and Fisher et al., 2015)
VPTCC               :  1.35e-4: L/kg kidney,         Volume of proximal tubule cells (60 million PTC cells/gram kidney, 1 PTC = 2250 um3)(Hsu et al., 2014)
FVBK                :    0.160: Unitless,            Blood volume fraction of kidney (Brown, 1997)
protein             :   2.0e-6: mg protein/PTCs,     Amount of protein in proximal tubule cells (Addis et al., 1936)

//Chemical-specific parameters 		

PL                  :   3.73  : Unitless,            Liver/ plasma PC; (In-house experiment) 
PK                  :   0.75  : Unitless,            Kidney/ plasma PC; (In-house experiment)
PLu                 :   0.63  : Unitless,            Lung/ plasma PC; (In-house experiment)
PF                  :   0.028 : Unitless,            Fat/ plasma PC; (Loccisano et al., 2012)
PRest               :   0.64  : Unitless,            Restofbody/ plasma PC; (Chou et al., 2019)
MW                  : 532.58  : g/mol,               6:2 Cl-PFESA  molecular mass  
MKC                 : 0.0046  : Unitless,            Fraction mass of kidney (percent of BW); Brown, 1997
Free                : 4.33e-3 : Unitless,            Free fraction; (In-house experiment)  
Vmax_baso_invitro   :   479   : pmol/mg Protein/min, Vmax of basolateral transporter; averaged in vitro value of rOAT1 and rOAT3 (Nakagawa, 2007); initial value asumsed the same as PFOA (Worley and Fisher, 2015)
Km_baso             :  64.4   : mg/L,                Km of basolateral transpoter, averaged in vitro value of rOAT1 and rOAT3 (Nakagawa, 2007); initial value asumsed the same as PFOA (Worley and Fisher, 2015)
Vmax_apical_invitro : 51803   : pmol/mg protein/min, Vmax of apical transporter; averaged invitro value of Oatp1a1 (Weaver, 2010);Fitting by Chou and Lin (2019)
Km_apical           :  20.1   : mg/L,                Km of apical transpoter, in vitro value for Oatp1a1 (Weaver, 2010);Fitting by Chou and Lin (2019)
RAFbaso             :   1     : Unitless             Relative activity factor, basolateral transpoters (male) (fit to data); initial value asumsed the same as PFOA (Worley and Fisher, 2015)
RAFapi              :  0.001  : Unitless             Relative acitivty factor, apical transpoters (fit to data); 0.001356 (female);initial value asumsed the same as PFOA (Worley and Fisher, 2015)
GFRC                :  24.19  : L/hr/kg kiney,       Glomerular filtration rate (male) (Corley, 2005)
Kdif                :  0.001  : L/h,                 Diffusion rate from proximal tubule cells to kidney serum; Fitting by Chou and Lin (2019)
Kabsc               :   0.015 : 1/(h*BW^-0.25),      Rate of absorption of 6:2 Cl-PFESA  from small intestine to liver; (In-house experiment)                         
KunabsC             : 1.49e-3 : 1/(h*BW^-0.25),      Rate of unabsorbed dose to appear in feces; (In-house experiment) 
GEC                 :  3.51   : 1/(h*BW^0.25),       Gastric emptying time  (Yang et al., 2013)
K0C                 :   0.13  : 1/(h*BW^-0.25),      Rate of uptake from the stomach into the liver (initial value assumed the same as PFOA (1) from Worley and Fisher, 2015 and then re-fitting)
KeffluxC            :   0.74  : 1/(h*BW^-0.25),      Rate of clearance of 6:2 Cl-PFESA  from proximal tubule cells into blood (initial value assumed the same as PFOA (2.49) from Worley and Fisher, 2015 and then re-fitting by Chou and Lin (2019))
KbileC              : 3.55e-4 : 1/(h*BW^-0.25),      Biliary elimination rate (male); liver to feces storage (initial value assumed the same as PFOA (0.004) from Worley and Fisher, 2015 and then re-fitting by Chou and Lin (2019))  
KurineC             : 5.68e-6 : 1/(h*BW^-0.25),      Rate of urine elimination from urine storage (male); (In-house experiment) 
HL                  :  5588   : day,                 Half-Live of rats (Shi et al.,2016)

$MAIN

// #+ Time varabiles: day and age
// #+ YEAR          

double YEAR           = TIME/(24*365); 

double QC = QCC*pow(BW, 0.75)*(1-Htc);               // L/h, Cardiac output (adjusted for plasma)
double QK = QKC*QC;                                  // L/h, Plasma flow to kidney
double QL = QLC*QC;                                  // L/h, Plasma flow to liver
double QLu= QLuC*QC;                                 // L/h, Plasma flow to lung
double QF = QFC*QC;                                  // L/h, Plasma flow to Fat
double QRest = QC-QK-QL-QF;                          // L/h, Plasma flow to the rest of body

double VL = VLC*BW;                                  // L,   Volume of liver 
double VLu= VLuC*BW;                                 // L,   Volume of lung
double VF = VFC*BW;                                  // L,   Volume of Fat
double VPlas = VPlasC*BW;                            // L,   Volume of plasma
double VK = VKC*BW;                                  // L,   Volume of kidney 
double Vfil = VfilC*BW;                              // L,   Volume of filtrate  
double VPTC = VK*VPTCC;                              // L,   Volume of proximal tubule cells 
double VKb = VK*FVBK;                                // L,   Volume of blood in the kidney; fraction blood volume of kidney (0.16) from Brown, 1997
double VRest = (0.93*BW)-VL-VLu-VF-VK-VPlas-VPTC-Vfil; // L,   Rest of body; volume of remaining tissue (L); Revised original equation (VR = (0.93*BW) - VPlas - VPTC - Vfil - VL) from Worley and Fisher, 2015  
double MK = VK*1000;                                 // g,   kidney weight in gram
double PTC = MKC*6e7*1000;                           // cells/kg BW, Number of PTC (cells/kg BW) (based on 60 million PTC/gram kidney); Revised original equation (PTC = MKC*6e7) from Worley and Fisher, 2015
double MPTC = VPTC*1000;                             // g,           mass of the proximal tubule cells (assuming density 1 kg/L)
double Vmax_basoC = (Vmax_baso_invitro*RAFbaso*PTC*protein*60*(MW/1e12)*1000);         // mg/h/kg BW^0.75, Vmax of basolateral transporters (average Oat1 and Oat3) equation from Worley and Fisher, 2015
double Vmax_apicalC = (Vmax_apical_invitro*RAFapi*PTC*protein*60*(MW/1e12)*1000);      // mg/h/kg BW^0.75, Vmax of apical transpoters in in vitro studies (Oatp1a1) equation from Worley and Fisher, 2015
double Vmax_baso = Vmax_basoC*pow(BW,0.75);          // mg/h, Vmax of basolateral transporters
double Vmax_apical = Vmax_apicalC*pow(BW,0.75);      // mg/h, Vmax of apical transpoters
double Kbile = KbileC*pow(BW,(-0.25));               // 1/h, Biliary elimination, liver to feces storage
double Kurine = KurineC*pow(BW,(-0.25));             // 1/h, Urinary elimination; 
double Kefflux = KeffluxC*pow(BW,(-0.25));           // 1/h, Efflux clearance rate from PTC to blood 
double GFR = GFRC*(MK/1000);                         // L/h, Glomerular filtration rate, scaled to mass of kidney 

//GI tract parameters
double Kabs = Kabsc*pow(BW,(-0.25));                 // 1/h, Rate of absorption of 6:2 Cl-PFESA  from small intestine to liver
double Kunabs = KunabsC*pow(BW,(-0.25));             // 1/h, Rate of unabsorbed dose to appear in feces
double GE = GEC*pow(BW,(-0.25));                     // 1/h, Gasric emptying time 
double K0 = K0C*pow(BW,(-0.25));                     // 1/h, Rate of uptake from the stomach into the liver

//Metabolic extraction rate of the liver, equation reference to (Verner et al.,2009) and (Zhang et al.,2022)
double Vd = (PL*VL+PK*VK+PLu*VLu+PF*VF+PRest*VRest)+VPlas;   //L, Volumes of the tissue
double Cle = log(2)/HL*Vd;                                   //L/d, Clearance rate of 6:2 Cl-PFESA 
double Clint = (Cle/(1-Cle/QL))/VL;                          //L/d/Kg, Intrinsic clearance value per kilogram of the liver 
double Eh = Clint*VL/(Clint*VL+QL);                          //Unitless, Metabolic extraction rate of the liver 

// Mass balance adjusted factor; avoding the negative occur at time = 0
double KDOSE        = (TIME==0)?0:1;

$CMT ADOSE A_baso A_apical Adif Aefflux ACI AUCCA_free APTC AUCCPTC AFil AUCCfil Aurine ARest AUCCRest AST 
AabsST ASI AabsSI Afeces AL AUCCL ALu AF AVPlas_free AAPlas_free AUCCV_free AUCCK AUCCF Ameta AUCCLu AKb

$ODE

// Concentrations in plasma
double CA_free  = AAPlas_free/(VPlas * 0.2);          // mg/L, Free 6:2 Cl-PFESA  concentration in the arterial plasma
double CV_free  = AVPlas_free/(VPlas * 0.8);          // mg/L, Free 6:2 Cl-PFESA  concentration in the venous plasma
double CA       = CA_free/Free;                       // mg/L, Concentration of total 6:2 Cl-PFESA  in the arterial plasma
double CV       = CV_free/Free;                       // mg/L, Concentration of total 6:2 Cl-PFESA  in the venous plasma
double Cplas    = ((AAPlas_free+AVPlas_free)/VPlas)/Free;// mg/L, Concentration of total 6:2 Cl-PFESA  in the plasma

// Concentrations in liver
double CL = AL/VL;                                   // mg/L, Concentration of 6:2 Cl-PFESA  in the liver compartment
double CVL = CL/PL;                                  // mg/L, Concentration of 6:2 Cl-PFESA  in venous plasma leaving liver
// Concentrations in kidney
double CKb = AKb/VKb;                                // mg/L, Concetraitons of 6:2 Cl-PFESA  in venous plasma leaving kidney 
double CVK = CKb;                                    // mg/L, Concentration of 6:2 Cl-PFESA  in plasma leaving kidney
double CK  = CVK*PK;                                 // mg/L, Concetraitons of 6:2 Cl-PFESA  in Kidney compartment
// Concentrations in lung
double CLu = ALu/VLu;                                // mg/L, Concentration of 6:2 Cl-PFESA  in the lung compartment
double CVLu= CLu/PLu;                                // mg/L, Concentration of 6:2 Cl-PFESA  in venous plasma leaving lung
// Concentrations in Fat
double CF  = AF/VF;                                  // mg/L, Concentration of 6:2 Cl-PFESA  in the fat compartment
double CVF = CF/PF;                                  // mg/L, Concentration of 6:2 Cl-PFESA  in venous plasma leaving fat
// Concentrations in Rest of body
double CRest = ARest/VRest;                          // mg/L, Concentration of 6:2 Cl-PFESA  in the rest of the body
double CVRest = ARest/(VRest*PRest);                 // mg/L, Concentration of 6:2 Cl-PFESA  in the venous palsma leaving the rest of body

// Kidney compartment plus 2 subcompartment (Proximal Tubule cells: PTCs, Filtrate: Fil)
// Concentration in kidney, PTCs and fil

double CPTC = APTC/VPTC;                             // mg/L, Concetraitons of 6:2 Cl-PFESA  in PTCs
double Cfil = AFil/Vfil;                             // mg/L, Concetraitons of 6:2 Cl-PFESA  in Fil

// Virtural compartment: 
// Basolateral (baso) 
// transport, Diffusion (dif), 
// Apical (apical) transport, and 
// efflux clearance (efflux)
// Clerance (CL) via glormerular filtration 

double RA_baso = (Vmax_baso*CKb)/(Km_baso+CKb);                                             // mg/h, Rate of basolateral transpoters 
dxdt_A_baso = RA_baso;                                                                      // mg,   Amount of basolateral transpoters
double RA_apical = (Vmax_apical*Cfil)/(Km_apical + Cfil);                                   // mg/h, Rate of apical transpoter 
dxdt_A_apical = RA_apical;                                                                  // mg,   Amount of apical transpoter
double Rdif = Kdif*(CKb - CPTC);                                                            // mg/h, Rate of diffusion from into the PTC
dxdt_Adif = Rdif;                                                                           // mg,   Amount moved via glomerular filtration
double RAefflux = Kefflux*APTC;                                                             // mg/h, Rate of efflux clearance rate from PTC to blood
dxdt_Aefflux = RAefflux;                                                                    // mg,   Amount of efflux clearance rate from PTC to blood
double RCI = CA*GFR*Free;                                                                   // mg/h, Rate of clerance (CL) to via glomerular filtration (GFR) 
dxdt_ACI = RCI;                                                                             // mg,   Amount of clearance via GFR

// {6:2 Cl-PFESA  distribution in each compartment}
// Free 6:2 Cl-PFESA  in plasma
double RV_free = (QRest*CVRest*Free)+(QK*CVK*Free)+(QL*CVL*Free)+(QF*CVF*Free)-(QLu*CV*Free)+RAefflux;  // mg/h,   Rate of change in the plasma
dxdt_AVPlas_free = RV_free;                                                               // mg,     Amount of free 6:2 Cl-PFESA  in the plasma
dxdt_AUCCV_free = CV_free;
double RA_free = (QLu*CVLu*Free)- (QC*CA*Free);  // mg/h,   Rate of change in the plasma
dxdt_AAPlas_free = RA_free;                                                               // mg,     Amount of free 6:2 Cl-PFESA  in the plasma
dxdt_AUCCA_free = CA_free;                                                                  // mg*h/L, Area under curve of free 6:2 Cl-PFESA  in plasma compartment

// Proximal Tubule Cells (PTCs)
double RPTC = Rdif + RA_apical + RA_baso - RAefflux;                                        // mg/h,   Rate of change in PTCs 
dxdt_APTC = RPTC;                                                                           // mg,     Amount in the PTCs
dxdt_AUCCPTC = CPTC;                                                                        // mg*h/L, Area under curve of 6:2 Cl-PFESA  in the compartment of PTCs

// Proximal Tubule Lumen/ Filtrate (Fil)
double Rfil = CA*GFR*Free - RA_apical - AFil*Kurine;                                        // mg/h,   Rate of change in the Fil
dxdt_AFil = Rfil;                                                                           // mg,     Amount in the Fil
dxdt_AUCCfil = Cfil;                                                                         // mg*h/L, Area under curve of 6:2 Cl-PFESA  in the compartment of Fil

// Urine elimination
double Rurine = Kurine*AFil;                                                                // mg/h,   Rate of change in urine
dxdt_Aurine = Rurine;                                                                       // mg,     Amount in urine

// Kidney compartment
double RKb = QK*(CA-CVK)*Free - CA*GFR*Free - Rdif - RA_baso;                               // mg/h,   Rate of change in Kidney compartment
dxdt_AKb = RKb;                                                                             // mg,     Amount in kidney compartment
dxdt_AUCCK= CK;                                                                            // mg*h/L, Area under curve of 6:2 Cl-PFESA  in the Kidney compartment

// 6:2 Cl-PFESA  in the compartment of rest of body, flow-limited model
double RRest = QRest*(CA-CVRest)*Free;                                                      // mg/h,   Rate of change in rest of body
dxdt_ARest = RRest;                                                                         // mg,     Amount in rest of body 
dxdt_AUCCRest = CRest;                                                                      // mg*h/L, Area under curve of 6:2 Cl-PFESA  in the compartment of rest of body

// 6:2 Cl-PFESA  in the compartment of Lung, flow-limited model
double RLu = QLu*(CV-CVLu)*Free;                                                            // mg/h,   Rate of change in Lung
dxdt_ALu   = RLu;                                                                           // mg,     Amount in rest of body 
dxdt_AUCCLu= CLu;                                                                           // mg*h/L, Area under curve of 6:2 Cl-PFESA  in the compartment of rest of body

// 6:2 Cl-PFESA  in the compartment of fat, flow-limited model
double RF  = QF*(CA-CVF)*Free;                                                              // mg/h,   Rate of change in rest of body
dxdt_AF = RF;                                                                         // mg,     Amount in rest of body 
dxdt_AUCCF = CF;                                                                      // mg*h/L, Area under curve of 6:2 Cl-PFESA  in the compartment of rest of body
    
// Gastrointestinal (GI) tract
// Stomach compartment
double RST = - K0*AST - GE*AST;                                                             // mg/h,   Rate of chagne in stomach caomprtment
dxdt_AST = RST;                                                                             // mg,     Amount in Stomach
double RabsST = K0*AST;                                                                     // mg/h,   Rate of absorption in the stomach
dxdt_AabsST = RabsST;                                                                       // mg,     Amount absorbed in the stomach

// Biliary excretion
double Rbile = Kbile*AL;                                                                    // mg,     Amount of 6:2 Cl-PFESA  in bile excretion

// Small intestine compartment
double RSI = GE*AST - Kabs*ASI - Kunabs*ASI + Kbile*AL;                                                // mg/h,   Rate of chagne in small intestine caomprtment
dxdt_ASI = RSI;                                                                             // mg,     Amount in small intestine
double RabsSI = Kabs*ASI;                                                                   // mg/h,   Rate of absorption in the small intestine
dxdt_AabsSI = RabsSI;                                                                       // mg,     Amount absorbed in the small intestine
double Total_oral_uptake = AabsSI + AabsST;                                                 // mg,     Total oral uptake in the GI

// Feces compartment
double Rfeces = Kunabs*ASI;                                                                 // mg/h,   Rate of change in feces compartment
dxdt_Afeces = Rfeces;                                                                       // mg,     Amount of the feces compartment

// 6:2 Cl-PFESA  in liver compartment, flow-limited model
double Rmeta = Eh*QL*CA*Free; 
dxdt_Ameta = Rmeta; 
double RL = QL*(CA-CVL)*Free - Kbile*AL - Eh*QL*CA*Free + Kabs*ASI + K0*AST;                                // mg/h,   Rate of chagne in liver caomprtment
dxdt_AL = RL;                                                                               // mg,     Amount in liver compartment
dxdt_AUCCL = CL;                                                                            // mg*h/L, Area under curve of 6:2 Cl-PFESA  in liver compartment

// #+ Virtural compartment for estmating input dose
dxdt_ADOSE        = 0;

// {Mass balance equations}
double Qbal = QC-QL-QK-QRest-QF;
double Tmass = AAPlas_free + AVPlas_free + ARest + AKb + AFil + APTC + AL + AST + ASI + ALu + AF;
double Loss = Aurine + Afeces + Ameta;
double ATotal = Tmass + Loss;
double Bal = ADOSE - ATotal*KDOSE;

$TABLE  
capture VP     = CV_free/Free;
capture AP     = CA_free/Free;
capture Plas   = ((AAPlas_free+AVPlas_free)/VPlas)/Free;
capture Liver  = AL/VL;
capture Kidney = CK;
capture Lung   = ALu/VLu;
capture Fat    = AF/VF;
capture Rest   = ARest/VRest;
capture AUC_CA = AUCCA_free;
capture AUC_CV = AUCCV_free;
capture AUC_CL = AUCCL;
capture AUC_CK = AUCCK;
capture AUC_CLu= AUCCLu;
capture AUC_CF = AUCCF;
capture Balance= Bal;
capture QB     = Qbal;
capture AT     = Tmass;'
