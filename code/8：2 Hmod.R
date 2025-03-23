HumanPBPK.code <- '
$PARAM @annotated

//Physiological parameters		

BW                  :   63    : kg,                  GBZ/T200.1-2007
QCC                 :   15    : L/h/kg^0.75,         GBZT 200.1-2007
QLC                 :   0.255 : Unitless,            FGBZ/T200.3-2014
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
FVBK                :    0.160: Unitless,            Blood volume fraction of kidney (Brown, 1997)

//Chemical-specific parameters 		

PL                  :   8.78  : Unitless,            Liver/ plasma PC; (In-house experiment) 
PK                  :   0.72  : Unitless,            Kidney/ plasma PC; (In-house experiment)
PLu                 :   1.49  : Unitless,            Lung/ plasma PC; (In-house experiment)
PF                  :   0.064 : Unitless,            Fat/ plasma PC; (Loccisano et al., 2012)
PRest               :   1.90  : Unitless,            Restofbody/ plasma PC; (Chou et al., 2019)
MW                  : 632.60  : g/mol,               8:2 Cl-PFESA molecular mass  
Free                :  0.009  : Unitless,            Free fraction; (In-house experiment)  
Tm                  :   2032  : mg/h,                8:2 Cl-PFESA Transport maximum; (original data from Kim et al., 2018)
Kt	                :   30    : mg/L,                8:2 Cl-PFESA Transport affinity constant; (original data from Kim et al., 2018)
Kabsc               :   0.010 : 1/(h*BW^-0.25),      Rate of absorption of 8:2 Cl-PFESA  from small intestine to liver; (In-house experiment)                         
KunabsC             : 7.93e-5 : 1/(h*BW^-0.25),      Rate of unabsorbed dose to appear in feces; (In-house experiment) 
GEC                 :  3.51   : 1/(h*BW^0.25),       Gastric emptying time  (Yang et al., 2013)
K0C                 :  0.14   : 1/(h*BW^-0.25),      Rate of uptake from the stomach into the liver (initial value assumed the same as PFOA (1) from Worley and Fisher, 2015 and then re-fitting)
KbileC              : 2.63e-4 : 1/(h*BW^-0.25),      Biliary elimination rate (male); liver to feces storage (initial value assumed the same as PFOA (0.004) from Worley and Fisher, 2015 and then re-fitting by Chou and Lin (2019))  
KurineC             : 4.54e-6 : 1/(h*BW^-0.25),      Rate of urine elimination from urine storage (male); (In-house experiment) 

$MAIN

// #+ Time varabiles: day and age
// #+ YEAR          

double YEAR           = TIME/(24*365); 

double QC = QCC*pow(BW, 0.75)*(1-Htc);               // L/h, Cardiac output (adjusted for plasma)
double QK = QKC*QC;                                  // L/h, Plasma flow to kidney
double QL = QLC*QC;                                  // L/h, Plasma flow to liver
double Qfil= 0.2*QK;  		                          // L/h, Plasma flow to filtrate compartment (L/h); 20% of QK (Kim et al., 2018)
double QLu= QLuC*QC;                                 // L/h, Plasma flow to lung
double QF = QFC*QC;                                  // L/h, Plasma flow to Fat
double QRest = QC-QK-QL-QF-Qfil;                     // L/h, Plasma flow to the rest of body

double VL = VLC*BW;                                  // L,   Volume of liver 
double VLu= VLuC*BW;                                 // L,   Volume of lung
double VF = VFC*BW;                                  // L,   Volume of Fat
double VPlas = VPlasC*BW;                            // L,   Volume of plasma
double VK = VKC*BW;                                  // L,   Volume of kidney 
double Vfil = VfilC*BW;                              // L,   Volume of filtrate  
double VKb = VK*FVBK;                                // L,   Volume of blood in the kidney; fraction blood volume of kidney (0.16) from Brown, 1997
double VRest = (0.93*BW)-VL-VLu-VF-VK-VPlas-Vfil;    // L,   Rest of body; volume of remaining tissue (L); Revised original equation (VR = (0.93*BW) - VPlas - VPTC - Vfil - VL) from Worley and Fisher, 2015  
double Kbile = KbileC*pow(BW,(-0.25));               // 1/h, Biliary elimination, liver to feces storage
double Kurine = KurineC*pow(BW,(-0.25));             // 1/h, Urinary elimination; 

//GI tract parameters
double Kabs = Kabsc*pow(BW,(-0.25));                 // 1/h, Rate of absorption of 8:2 Cl-PFESA  from small intestine to liver
double Kunabs = KunabsC*pow(BW,(-0.25));             // 1/h, Rate of unabsorbed dose to appear in feces
double GE = GEC*pow(BW,(-0.25));                     // 1/h, Gasric emptying time 
double K0 = K0C*pow(BW,(-0.25));                     // 1/h, Rate of uptake from the stomach into the liver


// Mass balance adjusted factor; avoding the negative occur at time = 0
double KDOSE        = (TIME==0)?0:1;

$CMT ADOSE AFil AUCCfil Aurine ARest AUCCRest AST AabsST ASI AabsSI Afeces 
AL AUCCL ALu AF AVPlas_free AAPlas_free AUCCV_free AUCCA_free AUCCK AUCCF AUCCLu AKb

$ODE

// Concentrations in plasma
double CA_free  = AAPlas_free/(VPlas * 0.2);          // mg/L, Free 8:2 Cl-PFESA  concentration in the arterial plasma
double CV_free  = AVPlas_free/(VPlas * 0.8);          // mg/L, Free 8:2 Cl-PFESA  concentration in the venous plasma
double CA       = CA_free/Free;                       // mg/L, Concentration of total 8:2 Cl-PFESA  in the arterial plasma
double CV       = CV_free/Free;                       // mg/L, Concentration of total 8:2 Cl-PFESA  in the venous plasma
double Cplas    = ((AAPlas_free+AVPlas_free)/VPlas)/Free;// mg/L, Concentration of total 8:2 Cl-PFESA  in the plasma

// Concentrations in liver
double CL = AL/VL;                                   // mg/L, Concentration of 8:2 Cl-PFESA  in the liver compartment
double CVL = CL/PL;                                  // mg/L, Concentration of 8:2 Cl-PFESA  in venous plasma leaving liver

// Concentrations in kidney
double CKb = AKb/VKb;                                // mg/L, Concetraitons of 8:2 Cl-PFESA  in venous plasma leaving kidney 
double CVK = CKb;                                    // mg/L, Concentration of 8:2 Cl-PFESA  in plasma leaving kidney
double CK  = CVK*PK;                                 // mg/L, Concetraitons of 8:2 Cl-PFESA  in Kidney compartment
double Cfil = AFil/Vfil;                             // mg/L, Concetraitons of 8:2 Cl-PFESA  in Fill

// Concentrations in lung
double CLu = ALu/VLu;                                // mg/L, Concentration of 8:2 Cl-PFESA  in the lung compartment
double CVLu= CLu/PLu;                                // mg/L, Concentration of 8:2 Cl-PFESA  in venous plasma leaving lung

// Concentrations in Fat
double CF  = AF/VF;                                  // mg/L, Concentration of 8:2 Cl-PFESA  in the fat compartment
double CVF = CF/PF;                                  // mg/L, Concentration of 8:2 Cl-PFESA  in venous plasma leaving fat

// Concentrations in Rest of body
double CRest = ARest/VRest;                          // mg/L, Concentration of 8:2 Cl-PFESA  in the rest of the body
double CVRest = CRest/PRest;                         // mg/L, Concentration of 8:2 Cl-PFESA  in the venous palsma leaving the rest of body

// {8:2 Cl-PFESA  distribution in each compartment}
// Free 8:2 Cl-PFESA in plasma
double RV_free = (QRest*CVRest*Free)+(QK*CVK*Free)+(QL*CVL*Free)+(QF*CVF*Free)-(QLu*CV*Free);  // mg/h,   Rate of change in the plasma
dxdt_AVPlas_free = RV_free;                                                               // mg,     Amount of free 8:2 Cl-PFESA  in the plasma
dxdt_AUCCV_free = CV_free;
double RA_free = (QLu*CVLu*Free)- (QC*CA*Free);  // mg/h,   Rate of change in the plasma
dxdt_AAPlas_free = RA_free;                                                               // mg,     Amount of free 8:2 Cl-PFESA  in the plasma
dxdt_AUCCA_free = CA_free;                                                                  // mg*h/L, Area under curve of free 8:2 Cl-PFESA  in plasma compartment

// Proximal Tubule Lumen/ Filtrate (Fil)
double Rfil = Qfil*(CA-Cfil)*Free- Tm*Cfil/(Kt+Cfil) - AFil*Kurine;                                        // mg/h,   Rate of change in the Fil
dxdt_AFil = Rfil;                                                                           // mg,     Amount in the Fil
dxdt_AUCCfil = Cfil;                                                                         // mg*h/L, Area under curve of 8:2 Cl-PFESA  in the compartment of Fil

// Urine elimination
double Rurine = (Qfil*Cfil*Free) - Kurine*AFil;                                                                // mg/h,   Rate of change in urine
dxdt_Aurine = Rurine;                                                                       // mg,     Amount in urine

// Kidney compartment
double RKb = QK*(CA-CVK)*Free + Tm*Cfil/(Kt+Cfil);                                         // mg/h,   Rate of change in Kidney compartment
dxdt_AKb = RKb;                                                                             // mg,     Amount in kidney compartment
dxdt_AUCCK= CK;                                                                            // mg*h/L, Area under curve of 8:2 Cl-PFESA  in the Kidney compartment

// 8:2 Cl-PFESA  in the compartment of rest of body, flow-limited model
double RRest = QRest*(CA-CVRest)*Free;                                                      // mg/h,   Rate of change in rest of body
dxdt_ARest = RRest;                                                                         // mg,     Amount in rest of body 
dxdt_AUCCRest = CRest;                                                                      // mg*h/L, Area under curve of 8:2 Cl-PFESA  in the compartment of rest of body

// 8:2 Cl-PFESA  in the compartment of Lung, flow-limited model
double RLu = QLu*(CV-CVLu)*Free;                                                            // mg/h,   Rate of change in Lung
dxdt_ALu   = RLu;                                                                           // mg,     Amount in rest of body 
dxdt_AUCCLu= CLu;                                                                           // mg*h/L, Area under curve of 8:2 Cl-PFESA  in the compartment of rest of body

// 8:2 Cl-PFESA in the compartment of fat
double RF  = QF*(CA-CVF)*Free;                                                              // mg/h,   Rate of change in rest of body
dxdt_AF = RF;                                                                         // mg,     Amount in rest of body 
dxdt_AUCCF = CF;                                                                      // mg*h/L, Area under curve of 8:2 Cl-PFESA  in the compartment of rest of body
    
// Gastrointestinal (GI) tract
// Stomach compartment
double RST = - K0*AST - GE*AST;                                                             // mg/h,   Rate of chagne in stomach caomprtment
dxdt_AST = RST;                                                                             // mg,     Amount in Stomach
double RabsST = K0*AST;                                                                     // mg/h,   Rate of absorption in the stomach
dxdt_AabsST = RabsST;                                                                       // mg,     Amount absorbed in the stomach

// Small intestine compartment
double RSI = GE*AST - Kabs*ASI - Kunabs*ASI + Kbile*AL;                                                // mg/h,   Rate of chagne in small intestine caomprtment
dxdt_ASI = RSI;                                                                             // mg,     Amount in small intestine
double RabsSI = Kabs*ASI;                                                                   // mg/h,   Rate of absorption in the small intestine
dxdt_AabsSI = RabsSI;                                                                       // mg,     Amount absorbed in the small intestine
double Total_oral_uptake = AabsSI + AabsST;                                                 // mg,     Total oral uptake in the GI

// Biliary excretion
double Rbile = Kbile*AL;                                                                    // mg,     Amount of 8:2 Cl-PFESA  in bile excretion

// Feces compartment
double Rfeces = Kunabs*ASI;                                                                 // mg/h,   Rate of change in feces compartment
dxdt_Afeces = Rfeces;                                                                       // mg,     Amount of the feces compartment

// 8:2 Cl-PFESA  in liver compartment, flow-limited model
double RL = QL*(CA-CVL)*Free - Kbile*AL + Kabs*ASI + K0*AST;                                // mg/h,   Rate of chagne in liver caomprtment
dxdt_AL = RL;                                                                               // mg,     Amount in liver compartment
dxdt_AUCCL = CL;                                                                        // mg*h/L, Area under curve of 8:2 Cl-PFESA  in liver compartment

// #+ Virtural compartment for estmating input dose
dxdt_ADOSE        = 0;

// {Mass balance equations}
double Qbal = QC-QL-QK-QRest-QF-Qfil;
double Tmass = AAPlas_free + AVPlas_free + ARest + AKb + AFil + AL + AST + ASI + ALu + AF;
double Loss = Aurine + Afeces;
double Input = ADOSE;
double Bal = Input- Tmass - Loss;

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
capture AT     = Tmass;
'
