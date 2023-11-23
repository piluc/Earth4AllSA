#using GlobalSensitivity, Statistics, OrdinaryDiffEq, QuasiMonteCarlo, Plots
using Interpolations
using DifferentialEquations
using ModelingToolkit
using DelimitedFiles
include("parameters.jl")
include("initialisations.jl")
include("system.jl")
_variables = getinitialisations()
_params = getparameters()


global text = "add_equation!(eqs, TRHGS ~ (TRSS1980 * ( (OW + 297) / 297) ) * (AL / AL1980) )
add_equation!(eqs, D(CO2A) ~ CO2E - CO2AB + 2 * CO2FCH4)
add_equation!(eqs, KN2OEKF ~ KN2OKF1980 * exp(-RDN2OKF * (t - 1980)) * IfElse.ifelse(t > 2022, exp(-ERDN2OKF2022 * (t - 2022)), 1))
add_equation!(eqs, MMN2OE ~ FEUS * KN2OEKF / 1000)
add_equation!(eqs, NN2OE ~ withlookup(t, [(1980.0, 0.009), (2020.0, 0.009), (2099.27, 0.0)]))
add_equation!(eqs, N2OE ~ NN2OE + MMN2OE)
add_equation!(eqs, N2OC1980 ~ N2OA1980 / MAT)
add_equation!(eqs, D(N2OA) ~ N2OE - N2OBD)
add_equation!(eqs, N2OBD ~ N2OA / LN2OA)
add_equation!(eqs, N2OCA ~ N2OA / GN2OPP)
add_equation!(eqs, N2OFPP ~ withlookup(t, [(1980.0, 0.43), (2000.0, 0.64), (2010.0, 0.73), (2020.0, 0.8), (2100.0, 1.0)]))
add_equation!(eqs, FN2O ~ N2OCA * N2OFPP)
add_equation!(eqs, KCH4EKC ~ KCH4KC1980 * exp(-RDCH4KC * (t - 1980)) * IfElse.ifelse(t > 2022, exp(-ERDCH4KC2022 * (t - 2022)), 1))
add_equation!(eqs, MMCH4E ~ CRSU * KCH4EKC / 1000)
add_equation!(eqs, NCH4E ~ withlookup(t, [(1980.0, 0.19), (2020.0, 0.19), (2100.0, 0.19)]))
add_equation!(eqs, CH4E ~ NCH4E + MMCH4E)
add_equation!(eqs, CH4C1980 ~ CH4A1980 / MAT)
add_equation!(eqs, D(CH4A) ~ CH4E - CH4BD)
add_equation!(eqs, CH4BD ~ CH4A / LCH4A)
add_equation!(eqs, CH4CA ~ CH4A / GCH4PP)
add_equation!(eqs, CH4FPP ~ withlookup(t, [(1980.0, 0.82), (2000.0, 0.94), (2020.0, 1.01), (2100.0, 1.1)]))
add_equation!(eqs, FCH4 ~ CH4CA * CH4FPP)
add_equation!(eqs, OWLCO2 ~ IfElse.ifelse(t > 2022, 1 + SOWLCO2 * (OW / OBWA2022 - 1), 1))
add_equation!(eqs, LECO2A ~ LECO2A1980 * OWLCO2)
add_equation!(eqs, CO2FCH4 ~ CH4BD * TCO2PTCH4)
add_equation!(eqs, CO2AB ~ (CO2A - CO2A1850) / LECO2A)
add_equation!(eqs, CO2E ~ CO2EI + CO2ELULUC - DACCO2)
add_equation!(eqs, CO2GDP ~ (CO2E / GDP) * 1000)
add_equation!(eqs, CAC ~ DACCO2 * CCCSt)
add_equation!(eqs, DACCO2 ~ IfElse.ifelse(t > 2022, ramp(t, (DACCO22100) / IPP, 2022, 2022 + IPP), 0))
add_equation!(eqs, CO2CA ~ CO2A / GCO2PP)
add_equation!(eqs, CO2FPP ~ withlookup(t, [(1980.0, 0.0032), (1990.0, 0.0041), (2000.0, 0.0046), (2020.0, 0.0051), (2100.0, 0.006)]))
add_equation!(eqs, FCO2 ~ CO2CA * CO2FPP)
add_equation!(eqs, FOG ~ withlookup(t, [(1980.0, 0.18), (2000.0, 0.36), (2020.0, 0.39), (2050.0, 0.37), (2100.0, 0.0)]))
add_equation!(eqs, MMF ~ FCO2 + FOG + FCH4 + FN2O)
add_equation!(eqs, GHGE ~ CO2E * TCO2ETCO2 + CH4E * TCO2ETCH4 + N2OE * TCO2ETN2O)
add_equation!(eqs, AL1980  ~ (ISCEGA1980 * ALIS + (GLSU - ISCEGA1980) * ALGAV) / GLSU)
add_equation!(eqs, AL ~ (ISCEGA * ALIS + (GLSU - ISCEGA) * ALGAV) / GLSU)
add_equation!(eqs, HTS ~ EHS * TRHGS )
add_equation!(eqs, MRS ~ MRS1980 * (OW / WA1980))
add_equation!(eqs, MRDI ~ MRS / SVDR)
add_equation!(eqs, ECIM ~ MRDI * AI1980 * TPM3I * HRMI)
add_equation!(eqs, MEL ~ ISCEGA * MRS)
add_equation!(eqs, D(ISCEGA) ~ - MEL)
add_equation!(eqs, ISC ~ ISCEGA * 100)
add_equation!(eqs, WVC ~ WVC1980 * (1 + OWWV * (OW / WA1980 - 1) ))
add_equation!(eqs, WVF ~ WVF1980 * (1 + WVWVF * (WVC / WVC1980 - 1) ))
add_equation!(eqs, TMMF ~ MMF + WVF)
add_equation!(eqs, EWFF ~ (TMMF * GLSU) * 31.5 / 1000)
add_equation!(eqs, OW ~ WA1980 + (EHS - EH1980) * WFEH)
add_equation!(eqs, REHE ~ withlookup( OW, [(0.0, 1.0), (1.2, 4.8), (2.0, 8.6), (2.9, 14.0), (5.2, 40.0)]))
add_equation!(eqs, TRHGA ~ TRSA1980 * ( (OW + 287) / 287 ) )
add_equation!(eqs, HDO ~ EHS * TRHGA)
add_equation!(eqs, D(EHS) ~ EWFF - ECIM - HDO - HTS)
smooth!(eqs, PWA, OW, PD)
add_equation!(eqs, BCIL ~ CFWB + CFGB)
add_equation!(eqs, BCISNI ~ BCIL / NI)
add_equation!(eqs, BITRO ~ min(1, ITRO1980) + ramp(t, (ITRO2022 - ITRO1980) / 42, 1980, 2022) + ramp(t, (GITRO - ITRO2022) / 78, 2022, 2100))
add_equation!(eqs, CANCD ~ pulse(t, 2022, 1) * GD * FGDC2022)
add_equation!(eqs, CD ~ WCD - STW + OC - STO)
add_equation!(eqs, CFGB ~ GIC + GP - GND)
add_equation!(eqs, CFWB ~ WIC + WP - WND)
add_equation!(eqs, CPP ~ CD / POP)
add_equation!(eqs, CSGDP ~ CD / NI)
add_equation!(eqs, CONTR ~ CSGDP + GSGDP + SSGDP)
add_equation!(eqs, EGTF2022 ~ IfElse.ifelse(t > 2022, EGTRF2022 + EETF2022 + EPTF2022, 0) * NI)
smooth!(eqs, ETF2022, GETF2022, TINT)
add_equation!(eqs, ETTAF2022 ~ IfElse.ifelse(t > 2022, ECTAF2022 * FETACPET, 0))
smooth!(eqs, FGBW, GFGBW, TINT)
add_equation!(eqs, GCIN ~ GNI - CFGB)
add_equation!(eqs, D(GD) ~ GND - CANCD - GP)
add_equation!(eqs, GDB ~ GD / NI)
add_equation!(eqs, GETF2022 ~ EGTF2022 + ETTAF2022)
add_equation!(eqs, GFGBW ~ FT1980 + IfElse.ifelse(t > 2022, ETGBW, 0))
add_equation!(eqs, GFSNI ~ (GIC + GP) / NI)
add_equation!(eqs, GGI ~ WT + OT + STO + STW + IC2022)
add_equation!(eqs, GGIS ~ GGI / NI)
add_equation!(eqs, GIC ~ GD * GBC)
add_equation!(eqs, GIPC ~ PGCIN - GPU)
add_equation!(eqs, GND ~ max(0, (MGD - GD) / GDDP) + step(t, GSF2022, 2022) * NI)
add_equation!(eqs, GNI ~ (WT + OT + STO + STW + IC2022) - TP + ST)
add_equation!(eqs, GNISNI ~ GNI / NI)
add_equation!(eqs, GP ~ GD / GPP)
add_equation!(eqs, GPU ~ PGCIN * GCF)
add_equation!(eqs, GS ~ GPU + GIPC)
add_equation!(eqs, GSGDP ~ GS / NI)
add_equation!(eqs, IC2022 ~ NI * IfElse.ifelse(t > 2022, ramp(t, GEIC / IPP, 2022, 2020 + IPP), 0))
add_equation!(eqs, INEQ ~ OOIAT / WIAT)
add_equation!(eqs, INEQI ~ INEQ / INEQ1980)
add_equation!(eqs, ITO ~ BITRO * NI * (1 - WSO))
add_equation!(eqs, ITW ~ BITRW * NI * WSO)
add_equation!(eqs, MGD ~ NI * MGDB)
add_equation!(eqs, MWD ~ WI * MWDB)
add_equation!(eqs, OC ~ POCI * OCF)
add_equation!(eqs, OCF ~ 1 - OSF)
add_equation!(eqs, OCIN ~ OOIAT)
add_equation!(eqs, OI ~ NI * (1 - WSO))
add_equation!(eqs, OOIAT ~ OI - OT)
add_equation!(eqs, OS ~ POCI - OC)
add_equation!(eqs, OSF ~ OSF1980 * (1 + GDPOSR * (EGDPP / GDPP1980 - 1)))
add_equation!(eqs, OT ~ ITO + ETF2022 * FETPO)
add_equation!(eqs, OTR ~ OT / OI)
smooth!(eqs, PGCIN, GCIN, TAB)
smooth!(eqs, POCI, OCIN, TAOC)
smooth!(eqs, PWCIN, WCIN, TAWC)
add_equation!(eqs, SSGDP ~ TS / NI)
add_equation!(eqs, ST ~ STW + STO)
add_equation!(eqs, STO ~ OC * STR)
add_equation!(eqs, STW ~ WCD * STR)
add_equation!(eqs, TP ~ (WT + OT + STO + STW + IC2022) * FGBW)
add_equation!(eqs, TPP ~ WCIN + GCIN + OCIN - ST)
add_equation!(eqs, TS ~ OS + WS)
add_equation!(eqs, WCD ~ PWCIN * WCF)
add_equation!(eqs, WCIN ~ WIAT - CFWB)
add_equation!(eqs, D(WD) ~ WND - WP)
add_equation!(eqs, WDI ~ PWCIN / WF)
add_equation!(eqs, WFCSI ~ CFWB / WIAT)
add_equation!(eqs, WI ~ NI * WSO)
add_equation!(eqs, WIAT ~ WI - WT + TP)
add_equation!(eqs, WIC ~ WD * WBC)
add_equation!(eqs, WT ~ ITW + ETF2022 * (1 - FETPO))
add_equation!(eqs, WTR ~ WT / WI)
add_equation!(eqs, WND ~ max(0, (MWD - WD) / WDP))
add_equation!(eqs, WP ~ WD / WPP)
add_equation!(eqs, WDB ~ WD / WIAT)
add_equation!(eqs, WS ~ PWCIN - WCD)
add_equation!(eqs, CNEL ~ NEP * CNED)
add_equation!(eqs, NFCO2PP ~ MNFCO2PP * (1 - exp(- ( GDPP / 10) ) ) )
add_equation!(eqs, CO2NFIP ~ (NFCO2PP / 1000) * POP * (1 - FCO2SCCS))
add_equation!(eqs, FCO2SCCS ~ FCO2SCCS2022 + ramp(t, (GFCO2SCCS - FCO2SCCS2022) / IPP, 2022, 2022 + IPP))
add_equation!(eqs, CO2EI ~ CO2EP + CO2NFIP)
add_equation!(eqs, CO2EP ~ UFF * (TCO2PT / 1000) * (1 - FCO2SCCS))
add_equation!(eqs, CO2EMPP ~ (CO2EI / POP) * 1000)
add_equation!(eqs, CCCSG ~ CCCSt * ICCSC)
add_equation!(eqs, ICCSC ~ FCO2SCCS * (CO2NFIP + CO2EP) / (1 - FCO2SCCS))
add_equation!(eqs, TCO2PT ~ 2.8 * exp(ROCTCO2PT * (t - 1980)))
add_equation!(eqs, D(EEPI2022) ~ IEEPI)
add_equation!(eqs, IEEPI ~ EROCEPA2022 * 0 + step(t, EROCEPA2022, 2022))
add_equation!(eqs, TPPUEBEE ~ withlookup( GDPP, [(0.0, 0.0), (10.0, 4.0), (20.0, 7.0), (30.0, 9.0), (50.0, 12.0), (65.0, 13.0)]))
add_equation!(eqs, TPPUFFNEUBEE ~ withlookup( GDPP, [(0.0, 0.3), (15.0, 2.0), (25.0, 3.1), (35.0, 4.0), (50.0, 5.0)]))
add_equation!(eqs, DEBNE ~ (POP * TPPUEBEE * exp(-NIEE * (t - 1980))) / EEPI2022)
add_equation!(eqs, DFFNEUBNE ~ (POP * TPPUFFNEUBEE * exp(-NIEE * (t - 1980))) / EEPI2022)
add_equation!(eqs, FNE ~ FNE1980 + ramp(t, (FNE2022 - FNE1980) / 42, 1980, 2022) + ramp(t, (GFNE - FNE2022) / IPP, 2022, 2022 + IPP))
add_equation!(eqs, ERDNEFFFNE ~ FNE * DFFNEUBNE)
add_equation!(eqs, CNE ~ (ECRUNEFF / 1000) * ERDNEFFFNE)
add_equation!(eqs, EIDEFNE ~ ERDNEFFFNE * EUEPRUNEFF)
add_equation!(eqs, DFFFNEU ~ DFFNEUBNE - ERDNEFFFNE)
add_equation!(eqs, UFF ~ DFFFNEU + FFE)
add_equation!(eqs, DE ~ DEBNE + EIDEFNE)
add_equation!(eqs, DRES ~ REFF1980 + ramp(t, (REFF2022 - REFF1980) / 42, 1980, 2022) + ramp(t, (GREF - REFF2022) / IPP, 2022, 2022 + IPP))
add_equation!(eqs, DSRE ~ DE * DRES)
add_equation!(eqs, DREC ~ DSRE / RCUT)
add_equation!(eqs, DRECC ~ DREC - REC)
add_equation!(eqs, D(REC) ~ AREC - DIREC)
add_equation!(eqs, AREC ~ max(0, (DRECC / RECT) + (DIREC)))
add_equation!(eqs, DIREC ~ REC / LREC)
add_equation!(eqs, ASWC ~ AREC)
add_equation!(eqs, D(ACSWCF1980) ~ ASWC)
add_equation!(eqs, NDSWC ~ log(2) + log( ACSWCF1980 / SWC1980))
add_equation!(eqs, CISWC ~ (1 - CRDSWC) ^ NDSWC)
add_equation!(eqs, CAPEXRED ~ CAPEXRE1980 * CISWC)
add_equation!(eqs, CAPEXREG ~ CAPEXRED * AREC)
add_equation!(eqs, OPEXREG ~ OPEXRED * REP)
add_equation!(eqs, CRE ~ CAPEXREG + OPEXREG)
add_equation!(eqs, CAPEXFEG ~ CAPEXFED * AFEC)
add_equation!(eqs, OPEXFEG ~ OPEXFED * FEP)
add_equation!(eqs, CFE ~ CAPEXFEG + OPEXFEG)
add_equation!(eqs, CEL ~ CFE + CRE + CNEL)
add_equation!(eqs, REP ~ REC * RCUT)
add_equation!(eqs, GHMH2 ~ REP * FREH / KWEPKGH2)
add_equation!(eqs, GHMt ~ GHMH2 * TPTH2)
add_equation!(eqs, RHP ~ BEM + GHMt)
add_equation!(eqs, TWEPEJEE ~ TWHPEJCE * EFPP)
add_equation!(eqs, IIASAREP ~ REP / TWEPEJEE + RHP / MTPEJCE)
add_equation!(eqs, FTWEPMt ~ TWEPEJEE / MTPEJCE)
add_equation!(eqs, IIASAFEP ~ UFF / MTPEJCE)
add_equation!(eqs, LCEP ~ REP + NEP)
add_equation!(eqs, DFE ~ max(0, DE - LCEP))
add_equation!(eqs, DFEC ~ DFE / EKHPY)
add_equation!(eqs, DFECC ~ (DFEC - FEC) / FECCT + DIFEC)
add_equation!(eqs, AFEC ~ max(0, DFECC))
add_equation!(eqs, D(FEC) ~ AFEC - DIFEC)
add_equation!(eqs, LFEC ~ NLFEC * FCUTLOFC)
add_equation!(eqs, DIFEC ~ FEC / LFEC)
add_equation!(eqs, FCUT ~ DFE / FEC)
add_equation!(eqs, FCUTLOFC ~ 1 + sFCUTLOFC * ((FCUT / EKHPY) - 1))
add_equation!(eqs, FEP ~ FEC * FCUT)
add_equation!(eqs, NC ~ withlookup(t, [(1980.0, 75.0), (2000.0, 310.0), (2020.0, 310.0), (2098.9, 310.0)]))
add_equation!(eqs, NEP ~ NC * NCUT)
add_equation!(eqs, EP ~ FEP + NEP + REP)
add_equation!(eqs, ELB ~ EP / DE)
add_equation!(eqs, FFPNE ~ (FEP + NEP) / EP)
add_equation!(eqs, EU ~ DFFFNEU + EP / FTWEPMt + RHP)
add_equation!(eqs, EUPP ~ EU / POP)
add_equation!(eqs, FFE ~ FEP / FTWEPMt)
add_equation!(eqs, TCEG ~ (DEBNE * TCE / 1000) * AFMCM)
add_equation!(eqs, TCFFFNEUG ~ (DFFNEUBNE * TCFFFNEU / 1000) * AFMCM)
add_equation!(eqs, CFFFNEU ~ (DFFFNEU * TCFFFNEU) / 1000)
add_equation!(eqs, CG ~ EP * TC)
add_equation!(eqs, TGC ~ DEBNE * TC)
add_equation!(eqs, TCEN ~ TCEG + TCFFFNEUG + TGC)
add_equation!(eqs, TCENSGDP ~ TCEN / GDP)
add_equation!(eqs, CE ~ CFFFNEU + CEL + CG + CNE + CCCSG + CAC)
add_equation!(eqs, RECTEC ~ CE / TCEN)
add_equation!(eqs, CESGDP ~ CE / GDP)
add_equation!(eqs, ECETSGDP ~ IfElse.ifelse(t > 2022, (CE - TCEN) / GDP, 0))
add_equation!(eqs, CBC ~ CCSD + NCCR)
add_equation!(eqs, CBC1980 ~ NSR + NBBM + NBOM + NCCR)
smooth!(eqs, CCSD, TIR + NBOM, FSRT)
add_equation!(eqs, D(CBSR) ~ CSR)
add_equation!(eqs, CSR ~ (ISR - CBSR) / SRAT)
smooth!(eqs, ELTI, PI, IEFT)
add_equation!(eqs, GBC ~ TIR)
add_equation!(eqs, ISR ~ NSR * (1 + INSR * (PI / IT - 1) + UNSR * (PU / UT - 1)))
add_equation!(eqs, NCCR ~ 0.02 * (1 + GRCR * (OGR / 0.03 - 1)))
smooth!(eqs, PI, IR, IPTCB)
smooth!(eqs, PU, UR, UPTCB)
add_equation!(eqs, TGIR ~ GBC + ELTI)
add_equation!(eqs, TIR ~ CBSR + NBBM)
add_equation!(eqs, WBC ~ CCSD)
add_equation!(eqs, ACY ~ (DCYCA * (1 - FRA) * SQICA + CYRA * FRA) * CO2ELY * WELY)
add_equation!(eqs, ALFL ~ 1 - exp(-FFLR / TFFLR) )
add_equation!(eqs, AFSRA ~ 268 - SFU)
add_equation!(eqs, D(BALA) ~ CRLO)
add_equation!(eqs, BIUS ~ withlookup(t, [(1980.0, 0.0), (1990.0, 0.0), (2000.0, 0.0), (2020.0, 0.0), (2100.0, 0.0)]))
add_equation!(eqs, CEM ~ IfElse.ifelse(t > 2022, 1 - SSP2LMA * ramp(t, (1 - 0) / 78, 2022, 2100), 1))
add_equation!(eqs, CIRA ~ (1 - CRDRA) ^ NDRA)
add_equation!(eqs, CO2AFL ~ FOLA * (CO2AFLH / 1000) * CO2ELY * WELY)
add_equation!(eqs, CO2AFLH ~ 1.6 * FAM)
add_equation!(eqs, CO2ELULUC ~ CO2RFC - CO2AFL - ECO2ARA)
add_equation!(eqs, CO2ELY ~ IfElse.ifelse(t > 2022, 1 + CO2CEACY * (CO2CA / CO2C2022 - 1), 1))
add_equation!(eqs, CO2RFC ~ ( (OGRE + CREX) * CO2RHFC) / 1000 )
add_equation!(eqs, COFE ~ FEUS * CTF / 1000)
add_equation!(eqs, COFO ~ AFGDP * GDP + CRA + COFE)
add_equation!(eqs, CRA ~ (ECRA * RAA) / 1000)
add_equation!(eqs, CRBA ~ CRUS / DCS)
add_equation!(eqs, CRBI ~ BIUS * TCTB)
add_equation!(eqs, CRDE ~ (TUCERM + FERM + CRBI))
add_equation!(eqs, CREX ~ IfElse.ifelse( FOLA > 0, CRLA * CREXR, 0) * ALFL * CEM)
add_equation!(eqs, CREXR ~ 1 / 200 + CBECLE * (PCB - 1))
add_equation!(eqs, D(CRLA) ~ CREX - CRLO - UREX)
add_equation!(eqs, CRLO ~ CRLA * LER)
add_equation!(eqs, CRSU ~ ACY * CRLA)
add_equation!(eqs, CRUS ~ CRSU * (1 + CWR))
add_equation!(eqs, CRUSP ~ CRUS / POP)
add_equation!(eqs, CSQCA ~ ROCSQCA * SQICA)
add_equation!(eqs, CSRA ~ CYRA * CRLA * FRA)
add_equation!(eqs, CWR ~ ramp(t, GCWR / IPP, 2022, 2022 + IPP))
add_equation!(eqs, DCS ~ CRDE)
add_equation!(eqs, DCSCA ~ DCS - CSRA)
add_equation!(eqs, DCYCA ~ DCSCA / (CRLA * (1 - FRA) ) )
add_equation!(eqs, DRM ~ ((POP * DRMP) / 1000) * (1 - FNRM))
add_equation!(eqs, DRMP ~ TURMP)
add_equation!(eqs, ECFT ~ CRA - FCR)
add_equation!(eqs, ECFTSGDP ~ ECFT / GDP)
add_equation!(eqs, ECO2ARA ~ RAA * CO2ARA / 1000)
add_equation!(eqs, ECRA ~ ECRA22 * CIRA)
add_equation!(eqs, FAM ~ IfElse.ifelse(t > 2022, 1 + SSP2LMA * ramp(t, (MFAM - 1) / 78, 2022, 2100), 1))
add_equation!(eqs, FCR ~ (AFSRA / 1000) * RAA * CTF)
add_equation!(eqs, FEER ~ 1 + FUELER * (FUCA / SFU - 1))
add_equation!(eqs, FERM ~ RMF * KCKRM)
add_equation!(eqs, FEUS ~ CRLA * (1 - FRA) * FUCA / 1000)
add_equation!(eqs, FFI ~ FOFO / FF80)
add_equation!(eqs, FFLR ~ max(0, FOLA / initFOLA) )
add_equation!(eqs, FFLREOGRR ~ max(1, 1 + FFLREOGRRM * (FFLR - TFFLR) ) )
add_equation!(eqs, FNRM ~ ramp(t, GFNRM / IPP, 2022, 2022 + IPP))
add_equation!(eqs, FOFO ~ CRLA * FEUS)
add_equation!(eqs, D(FOLA) ~ NFL - CREX)
add_equation!(eqs, FPI ~ exp(ROCFP * (t - 1980)))
add_equation!(eqs, FRA ~ ramp(t, GFRA / IPP, 2022, 2020 + IPP))
add_equation!(eqs, FSPI ~ exp(ROCFSP * (t - 1980)) * IfElse.ifelse(t > 2022, exp(EROCFSP * (t - 2022)), 1))
add_equation!(eqs, FUCA ~ TFUCA / FPI)
add_equation!(eqs, FUP ~ (FEUS / POP) * 1000)
add_equation!(eqs, GLY ~ GLY80 + 0 * CO2CA - 0 * OW)
add_equation!(eqs, GLY80 ~ 14 * CO2ELY * WELY)
add_equation!(eqs, D(GRLA) ~ NGL)
add_equation!(eqs, IUL ~ POP * ULP)
add_equation!(eqs, LERM ~ IfElse.ifelse(t > 2022, 1 - SSP2LMA * ramp(t, (1 - 0) / 78, 2022, 2100), 1))
add_equation!(eqs, LER ~ LER80 * FEER * LERM)
add_equation!(eqs, LFL ~ OGRE + CREX)
add_equation!(eqs, LOCR ~ CRLO + UREX)
add_equation!(eqs, NDRA ~ log( (RAA + EGB22) / EGB22) / 0.693)
add_equation!(eqs, NFL ~ OGRE * (1 - FCG))
add_equation!(eqs, NGL ~ OGRE * FCG)
add_equation!(eqs, D(OGFA) ~ - NFL - NGL)
add_equation!(eqs, OGRE ~ OGFA * OGRR * OGRRM)
add_equation!(eqs, OGRR ~ OGRR80 * FFLREOGRR)
add_equation!(eqs, OGRRM ~ IfElse.ifelse(t > 2022, 1 - SSP2LMA * ramp(t, (1 - 0) / 78, 2022, 2100), 1))
add_equation!(eqs, PCB ~ CRBA / (1 + DRC))
add_equation!(eqs, PRMGL ~ GRLA * GLY / 1000)
add_equation!(eqs, RAA ~ CRLA * FRA)
add_equation!(eqs, RMF ~ DRM - RMGL)
add_equation!(eqs, RMGL ~ min( DRM, PRMGL))
add_equation!(eqs, RMSP ~ (RMGL + RMF * min(1, CRBA )) * 1000 / POP)
add_equation!(eqs, ROCSQCA ~ 0 + FUESQ * (FUCA / SFU - 1))
add_equation!(eqs, D(SQICA) ~ CSQCA)
add_equation!(eqs, TFA ~ OGFA + FOLA)
add_equation!(eqs, TFUCA ~ withlookup( DCYCA, [(1.0, 0.0), (2.0, 40.0), (2.5, 50.0), (3.0, 60.0), (3.5, 70.0), (4.5, 100.0), (6.5, 200.0), (10.0, 600.0)]))
add_equation!(eqs, TUC ~ TUCP * POP / 1000)
add_equation!(eqs, TUCERM ~ (TUC - TUFRM) / FSPI)
add_equation!(eqs, TUCERMP ~ TUCERM * 1000 / POP)
add_equation!(eqs, TUCP ~ withlookup( GDPP, [(0.0, 400.0), (6.1, 680.0), (8.7, 780.0), (13.9, 950.0), (20.0, 1050.0), (30.0, 1150.0), (40.0, 1250.0), (60.0, 1350.0), (100.0, 1550.0)]))
add_equation!(eqs, TUFRM ~ ( ( (TURMP / 1000) * POP) - RMGL) * KCKRM)
add_equation!(eqs, TURMP ~ withlookup( GDPP, [(0.0, 0.0), (6.1, 6.0), (8.8, 8.5), (14.0, 13.0), (30.0, 27.0), (40.0, 32.0), (50.0, 33.0), (100.0, 25.0)]))
add_equation!(eqs, UREX ~ max(0, (IUL - URLA) / UDT))
add_equation!(eqs, D(URLA) ~ UREX)
add_equation!(eqs, WELY ~ IfElse.ifelse(t > 2022, 1 + OWEACY * (OW / OW2022 - 1), 1))
add_equation!(eqs, CDDI ~ ROCDDI * DELDI)
add_equation!(eqs, CPI ~ PRIN * IR)
add_equation!(eqs, DEL ~ ((EPP / PPU) / (DELDI / DDI1980) ) * IfElse.ifelse(t > 1984, PNIS, 1))
add_equation!(eqs, D(DELDI) ~ CDDI)
add_equation!(eqs, DEPU ~ 0 + PH * pulse(t, 2020, 5))
add_equation!(eqs, DSWI ~ 1 + INVEOSWI * (PRI / DRI - 1))
smooth!(eqs, EPP, TPP, DAT)
add_equation!(eqs, GDP ~ OUTP * PPU)
add_equation!(eqs, IC ~ INV / RS)
add_equation!(eqs, D(INV) ~ OUTP - DEL)
add_equation!(eqs, IR ~ INVEOIN * (PRI / MRIWI - 1))
add_equation!(eqs, NI ~ SA)
add_equation!(eqs, OUTP ~ ORO * SSWI / SWI1980)
add_equation!(eqs, D(PRIN) ~ CPI)
add_equation!(eqs, PNIS ~ 1)
smooth!(eqs, PRI, (IC / DIC), ICPT)
add_equation!(eqs, ROCDDI ~ 0 + INVEODDI * (PRI / SRI - 1))
smooth!(eqs, RS, DEL, SAT)
add_equation!(eqs, SA ~ DEL * PPU)
smooth!(eqs, SSWI, DSWI, TAS)
add_equation!(eqs, AGIW ~ WARA * AHW)
add_equation!(eqs, AHW ~ NHW / PFTJ)
add_equation!(eqs, AHW1980 ~ initNHW / PFTJ80)
add_equation!(eqs, AVWO ~ WAP * LPR)
add_equation!(eqs, CECLR ~ ROCECLR * ECLR)
add_equation!(eqs, CHWO ~ (OPWO - WF) / HFD)
add_equation!(eqs, CWSO ~ WSO * ROCWSO)
add_equation!(eqs, CWRA ~ WARA * ROCWSO)
add_equation!(eqs, D(ECLR) ~ CECLR)
add_equation!(eqs, GDPPEROCCLR ~ max(0, 1 + GDPPEROCCLRM * (GDPP / GDPP1980 - 1)))
add_equation!(eqs, ENLPR2022 ~ ramp(t, GENLPR / IPP, 2022, 2022 + IPP))
add_equation!(eqs, HFD ~ TYLD / 3)
add_equation!(eqs, HWMGDPP ~ 1 + TENHW * (GDPP / GDPP1980 - 1))
add_equation!(eqs, IWEOCLR ~ 1 + WSOECLR * (WSO / initWSO - 1))
add_equation!(eqs, ILPR ~ NLPR - PSW)
add_equation!(eqs, LAPR ~ (OUTP * PRUN) / LAUS)
add_equation!(eqs, LAUS ~ WF * AHW)
add_equation!(eqs, LAUS80 ~ initWF * AHW1980)
smooth!(eqs, LPR, ILPR, TELLM)
add_equation!(eqs, LTEWSO ~ WSO * RWER)
smooth!(eqs, NHW, initNHW * HWMGDPP, TAHW)
add_equation!(eqs, NLPR ~ NLPR80 * (1 + WSOELPR * (WSO / initWSO - 1)) + ENLPR2022)
add_equation!(eqs, OCLR ~ ECLR * WEOCLR)
add_equation!(eqs, OPWO ~ (CAP / OCLR) * PFTJ)
add_equation!(eqs, PART ~ LPR * (1 - PURA))
add_equation!(eqs, PSW ~ AUR * (1 + PUELPR * (PURA / AUR - 1)))
smooth!(eqs, PURA, UR, UPT)
add_equation!(eqs, ROCECLR ~ ROCECLR80 * GDPPEROCCLR)
add_equation!(eqs, ROCWSO ~ withlookup( PURA / AUR, [(0.0, 0.06), (0.5, 0.02), (1.0, 0.0), (1.5, -0.007), (2.0, -0.01)]))
add_equation!(eqs, TCT ~ TYLD / 3)
add_equation!(eqs, UNEM ~ max(0, AVWO - WF))
add_equation!(eqs, UPT ~ TYLD / 3)
add_equation!(eqs, UR ~ UNEM / AVWO)
add_equation!(eqs, WAP ~ A20PA)
add_equation!(eqs, D(WARA) ~ CWRA - WRE)
add_equation!(eqs, WASH ~ WARA / LAPR)
smooth!(eqs, WEOCLR, IWEOCLR, TCT)
add_equation!(eqs, D(WF) ~ CHWO)
add_equation!(eqs, WRE ~ WARA * WRER)
add_equation!(eqs, WRER ~ IR * (1 - FIC))
add_equation!(eqs, D(WSO) ~ CWSO - LTEWSO)
add_equation!(eqs, CFETA ~ COFO + CE)
add_equation!(eqs, CTA ~ CFETA)
add_equation!(eqs, FB15 ~ 1 - (1 / (1 + exp(-LK * (GDPP - 14)))))
add_equation!(eqs, IEL ~ 1 + INELOK * (INEQ / 0.5 - 1))
add_equation!(eqs, LK ~ NK * IEL)
add_equation!(eqs, PB15 ~ POP * FB15)
add_equation!(eqs, RGGDPP ~ ( (GDPP - PGDPP) / PGDPP) / TEGR)
smooth!(eqs, PGDPP, GDPP, TEGR)
add_equation!(eqs, AVCA ~ TS + FCI)
add_equation!(eqs, CAPIS ~ CUCPIS / CTPIS)
add_equation!(eqs, CAPUS ~ CUCPUS / CTPUS)
add_equation!(eqs, CIPIS ~ max( (INCPIS + OBSGIPIS * GDP) / COCA, 0))
add_equation!(eqs, CIPUS ~ max( (GIPC + OBSGIPUS * GDP) / COCA, 0))
add_equation!(eqs, CDPIS ~ CPIS / LCPIS)
add_equation!(eqs, CDPUS ~ CPUS / LCPUS)
add_equation!(eqs, COCA ~ CC1980 * OWECC)
add_equation!(eqs, D(CPIS) ~ CAPIS - CDPIS)
add_equation!(eqs, D(CPUS) ~ CAPUS - CDPUS)
add_equation!(eqs, CRR ~ CAPIS / CPIS)
add_equation!(eqs, D(CUCPIS) ~ CIPIS - CAPIS)
add_equation!(eqs, CUCPIS1980 ~ (CAPPIS1980 / LCPIS1980) * CTPIS * EMCUC)
add_equation!(eqs, D(CUCPUS) ~ CIPUS - CAPUS)
add_equation!(eqs, CUCPUS1980 ~ (CAPPUS1980 / LCPUS1980) * CTPUS * EMCUC)
add_equation!(eqs, ECR ~ (ITFP - ETFP) * CRR)
add_equation!(eqs, D(ETFP) ~ ECR)
add_equation!(eqs, INCPIS ~ AVCA * FACNC)
add_equation!(eqs, LCPUS ~ LCPUS1980)
add_equation!(eqs, LCPUS1980 ~ 15 * OWELC)
add_equation!(eqs, OBSGIPIS ~ IfElse.ifelse(t > 2022, USPIS2022, 0))
add_equation!(eqs, OBSGIPUS ~ IfElse.ifelse(t > 2022, 0.01 + USPUS2022, 0.01))
add_equation!(eqs, OWECC ~ IfElse.ifelse(t > 2022, 1 + OWECCM * (OW / OW2022 - 1), 1))
add_equation!(eqs, OWELC ~ IfElse.ifelse(t > 2022, 1 + OWELCM * (OW / OW2022 - 1), 1))
add_equation!(eqs, CAP ~ CPIS + CPUS)
add_equation!(eqs, CBCEFCA ~ 1 + CBCEFRA * (CBC / CBC1980 - 1))
add_equation!(eqs, EDE ~ TPP / OOV)
add_equation!(eqs, EDEFCA ~ 1 + EDEFRA * (PEDE / ED1980 - 1))
add_equation!(eqs, EDELC ~ 1 + EDELCM * (PEDE / ED1980 - 1))
smooth!(eqs, FACNC, FRA1980 * FRACAMGDPPL * (WSOEFCA + CBCEFCA + EDEFCA) / 3, IPT)
add_equation!(eqs, FRACAMGDPPL ~ max( FRACAM, 1 + GDPPEFRACA * (GDPP / GDPP1980 - 1)))
add_equation!(eqs, FRACAMGDPPT ~ withlookup( GDPP / GDPP1980, [(0.0, 1.0), (1.0, 1.0), (2.0, 0.85), (2.1, 0.84), (4.0, 0.65), (8.0, 0.55), (16.0, 0.5)]))
add_equation!(eqs, ISGDP ~ (INCPIS + GIPC) / GDP)
add_equation!(eqs, LCPIS ~ (LCPIS1980 * OWELC) / EDELC)
smooth!(eqs, OGR, (ORO - OLY) / OLY, 1)
smooth!(eqs, OLY, ORO, 1)
add_equation!(eqs, OOV ~ ORO * PRUN)
add_equation!(eqs, ORO ~ OO1980 * ( (CPIS + CPUS) / (CAPPIS1980 + CAPPUS1980) ) ^ KAPPA * (LAUS / LAUS1980) ^ LAMBDA * (ETFP))
smooth!(eqs, PEDE, EDE, TOED)
add_equation!(eqs, WSOEFCA ~ 1 + WSOEFRA * (WASH / initWSO - 1))
add_equation!(eqs, D(A0020) ~ BIRTHS - PASS20)
add_equation!(eqs, D(A2040) ~ PASS20 - PASS40)
add_equation!(eqs, A20PA ~ A2040 + A4060 + A60PL - OP)
add_equation!(eqs, D(A4060) ~ PASS40 - PASS60)
add_equation!(eqs, D(A60PL) ~ PASS60 - DEATHS)
add_equation!(eqs, BIRTHR ~ BIRTHS / POP)
add_equation!(eqs, BIRTHS ~ A2040 * FW * (OF / FP))
add_equation!(eqs, CEFR ~ CMFR * EFR)
add_equation!(eqs, DEATHR ~ DEATHS / POP)
delay_n!(eqs, PASS60, RT_DEATHS, LV_DEATHS, LE60, ORDER)
add_equation!(eqs, DEATHS ~ RT_DEATHS[10])
add_equation!(eqs, DNC ~ ( (DNCM + (DNC80 - DNCM) * exp(-DNCG * (EGDPP - initEGDPP) ) ) * (1 + DNCA * (EGDPP - initEGDPP) ) ) * (1 - EFR) * FM)
add_equation!(eqs, DR ~ (A0020 + A60PL) / (A2040 + A4060))
add_equation!(eqs, EFR ~ ramp(t, GEFR / IPP, 2022, 2022 + IPP))
add_equation!(eqs, EPA ~ ramp(t, (GEPA - initEPA) / IPP, 2022, 22022 + IPP))
add_equation!(eqs, D(EGDPP) ~ (GDPP - EGDPP) / TAHI)
add_equation!(eqs, FM ~ IfElse.ifelse( SSP2FA2022F > 0, IfElse.ifelse(t > 2022, 1 + ramp(t, (MFM - 1) / 78, 2022, 2100), 1), 1))
add_equation!(eqs, GDPP ~ GDP / POP)
add_equation!(eqs, LE ~ ((LEMAX - (LEMAX - initLE) * exp(-LEG * (EGDPP - initEGDPP) )) * (1 + LEA * (EGDPP - initEGDPP) )) * WELE * LEM)
add_equation!(eqs, LE60 ~ LE - 60)
add_equation!(eqs, LEM ~ IfElse.ifelse( SSP2FA2022F > 0, IfElse.ifelse(t > 2022, 1 + ramp(t, (MLEM - 1) / 78, 2022, 2100), 1), 1))
add_equation!(eqs, OF ~ DNC * FADFS)
add_equation!(eqs, OP ~ A60PL * (LE - PA) / (LE - 60))
add_equation!(eqs, PA ~ IfElse.ifelse( LE < initLE, initPA, initPA + LEEPA * (LE + EPA - initLE) ))
delay_n!(eqs, BIRTHS, RT_PASS20, LV_PASS20, 20, ORDER)
add_equation!(eqs, PASS20 ~ RT_PASS20[10])
delay_n!(eqs, PASS20, RT_PASS40, LV_PASS40, 20, ORDER)
add_equation!(eqs, PASS40 ~ RT_PASS40[10])
delay_n!(eqs, PASS40, RT_PASS60, LV_PASS60, 20, ORDER)
add_equation!(eqs, PASS60 ~ RT_PASS60[10])
add_equation!(eqs, PGR ~ BIRTHR - DEATHR)
add_equation!(eqs, POP ~ A0020 + A2040 + A4060 + A60PL)
add_equation!(eqs, PW ~ OP / A20PA)
add_equation!(eqs, WELE ~ IfElse.ifelse(t > 2022, max(0, 1 + OWELE * (OW / OW2022 - 1)), 1))
add_equation!(eqs, CTFP ~ RTA * TFPEE5TA)
add_equation!(eqs, DRTA ~ (DROTA1980 + IfElse.ifelse(t > 2022, EDROTA2022, 0) * (1 + SCROTA * ((SC / SC1980) - 1))))
add_equation!(eqs, ECTAF2022 ~ max(0, CTA - CTA2022))
add_equation!(eqs, ECTAGDP ~ ECTAF2022 / GDP)
add_equation!(eqs, GSSGDP ~ GS / GDP)
add_equation!(eqs, IPR ~ CPUS / GPU)
add_equation!(eqs, IROTA ~ IfElse.ifelse(t > 2022, max(0, MIROTA2022 * (1 - 1 * (GDPP / GDPTL - 1))), 0))
add_equation!(eqs, ITFP ~ TFPEE5TA * OWTFP)
add_equation!(eqs, OWTFP ~ IfElse.ifelse(t > 2022, 1 + OWETFP * (OW / OW2022 - 1), 1))
add_equation!(eqs, PLUA ~ ECTAGDP * FUATA)
add_equation!(eqs, PPP ~ max(0, 1 + IPRVPSS * log( IPR / IPR1980) ))
add_equation!(eqs, PSEP ~ VPSS / POP)
add_equation!(eqs, PSP ~ GS / POP)
add_equation!(eqs, D(TFPEE5TA) ~ CTFP)
add_equation!(eqs, RROTAI ~ min(1, 1 + IIEEROTA * (INEQI / 1 - 1)))
add_equation!(eqs, RTA ~ (DRTA + 0) * RROTAI + IROTA)
smooth!(eqs, RTFPUA, PLUA, IPT + CTPIS)
add_equation!(eqs, SC ~ VPSS / GDP)
add_equation!(eqs, TFPIE5TA ~ TFPEE5TA * (1 - RTFPUA))
add_equation!(eqs, VPSS ~ GPU * PPP)
add_equation!(eqs, XECTAGDP ~ XETAC2022 + ramp(t, (XETAC2100 - XETAC2022) / 78, 2022, 2022 + 78))
add_equation!(eqs, AWBDI ~ exp(DRDI + log( WDI / TDI) ))
add_equation!(eqs, AWBIN ~ 1 + IEAWBIF * (INEQ / TI - 1))
add_equation!(eqs, AWBGW ~ max( MWBGW, min(1, 1 + GWEAWBGWF * (PWA / TW - 1))))
add_equation!(eqs, AWBI ~ (0.5 * AWBDI + 0.5 * AWBPS) * AWBIN * AWBGW * AWBP)
add_equation!(eqs, AWBP ~ (1 + PREAWBF * (ORP - TPR )) * WBEP)
add_equation!(eqs, AWBPS ~ exp(DRPS + log( PSP / TPS) ))
add_equation!(eqs, IEST ~ WorldDynamics.interpolate( INEQ / AI, tables[:IEST], ranges[:IEST]))
add_equation!(eqs, IPP ~ IfElse.ifelse( EIPF > 0, EIP, RD))
add_equation!(eqs, IRD ~ NRD * STRERD * STEERD)
add_equation!(eqs, IST ~ PSESTR * IEST)
smooth!(eqs, ORP, ( (AWBI - PAWBI) / AWBI ) / AWBPD, AWBPD)
smooth!(eqs, PAWBI, AWBI, AWBPD)
add_equation!(eqs, PSSGDP ~ PSP / GDPP)
add_equation!(eqs, PSESTR ~ WorldDynamics.interpolate( PSSGDP / SPS, tables[:PSESTR], ranges[:PSESTR]))
smooth!(eqs, RD, IRD, TCRD)
add_equation!(eqs, STE ~ 1 + PESTF * (ORP - AP))
add_equation!(eqs, STEERD ~ 1 + STEERDF * (STE / initSTE - 1))
add_equation!(eqs, STRERD ~ 1 + STRERDF * (SOTR / initSOTR - 1))
smooth!(eqs, SOTR, IST, TEST)
add_equation!(eqs, WBEP ~ 1 + PAEAWBF * (LPR / THPA - 1))"

function variables_index_init_con(sol)
    return sol[1]
end

function variables_index(e4a)
    nsp_v = ModelingToolkit.namespace_variables(e4a)
    dictionaryVariablesIndex = Dict{String,Int64}()
    dictionaryIndexVariables = Dict{Int64,String}()
    index = 0
    for i in nsp_v
        c = string(i)
        c = split(c, "₊")
        c_2 = c[2]
        c_2 = split(c_2, "(t)")
        variable = c_2[1]
        if (c_2[2] != "")
            order = replace(c_2[2], ")" => "")
            variable = variable * order
        end
        index = index + 1
        dictionaryVariablesIndex[variable] = index
        dictionaryIndexVariables[index] = variable
    end
    return [dictionaryVariablesIndex, dictionaryIndexVariables]
end

"This function returns a file text with the dictionary with variables name as keys and values as index"
function write_dict_variables_index()
    variables_dictionary = variables_index(e4a)[1]

    open("dict_k_variables_v_index.txt", "w") do io
        for line in variables_dictionary
            println(io, line)
        end
    end
end

"This function returns a file text with the dictionary  with key = index, and value = name of parameters "
function write_dict_index_variables()
    variables_dictionary = variables_index(e4a)[2]

    open("dict_k_index_v_variables.txt", "w") do io
        for line in variables_dictionary
            println(io, line)
        end
    end
end

"This function returns a file text with the dictionary  with key = name of parameters, and value = index of parameters "
function write_dict_parameters_index()
    parameters_dictionary = parameters_index_namespace(e4a)[1]

    open("dict_k_parameters_v_index.txt", "w") do io
        for line in parameters_dictionary
            println(io, line)
        end
    end
end

"This function returns a file text with the dictionary  with key = index, and value = name of parameters "
function write_dict_index_parameters()
    parameters_dictionary = parameters_index_namespace(e4a)[2]

    open("dict_index_parameters.txt", "w") do io
        for line in parameters_dictionary
            println(io, line)
        end
    end
end


function return_initial_condition_by_index(sol, dict, index)
    dictionaryIndexVariables = dict[2]
    return sol[dictionaryIndexVariables[index]][1]
end

function return_initial_condition_by_name(sol, arg)
    initial_condition_variable = sol[arg][1]
    return initial_condition_variable
end


function dict_name_and_initial_conditions(sol, e4a)
    nsp_v = ModelingToolkit.namespace_variables(e4a)
    dictionaryNameInitialCondition = Dict{Any,Any}()
    for i in nsp_v
        dictionaryNameInitialCondition[i] = sol[i][1]
    end
    return dictionaryNameInitialCondition
end

"This function build the array of the variable initial conditions and save it in a file txt. 
 NOTE: the variable values are sorted according to the indices"
function array_initial_conditions_variable(sol, e4a)
    nsp_v = ModelingToolkit.namespace_variables(e4a)
    dict_variables_index = read_file_variables()
    dictionaryNameInitialCondition = Dict{Any,Any}()
    dictionaryIndexInitialCondition = Dict{Int64,Float64}()
    for i in nsp_v
        c = string(i)
        c = split(c, "₊")
        c_2 = c[2]
        c_2 = split(c_2, "(t)")
        variable = c_2[1]
        if (c_2[2] != "")
            order = replace(c_2[2], ")" => "")
            variable = variable * order
        end
        dictionaryNameInitialCondition[variable] = sol[i][1]
    end
    for (k, v) in dict_variables_index
        for (j, t) in dictionaryNameInitialCondition
            if (k == j)
                dictionaryIndexInitialCondition[v] = t
            end
        end
    end

    array_initial_conditions = []
    index = 0
    open("variable_initial_conditions.txt", "w") do io
        print(io, "u0 = [")
        for key in sort(collect(keys(dictionaryIndexInitialCondition)))
            index = index + 1
            push!(array_initial_conditions, dictionaryIndexInitialCondition[key])
            print(io, dictionaryIndexInitialCondition[key])

            if (index % 10 == 0)
                println(io, ";")
            else
                print(io, ";")
            end
            #println("$key => $(dictionaryNameInitialCondition[key])")
        end
        print(io, "]")
    end
    return [dictionaryIndexInitialCondition, dictionaryNameInitialCondition]
end



"This function build the array of the parameters and save it in a file txt. 
 NOTE: the variable are sorted according to the indices"
function array_parameter_values(parameters=_params)
    dict_parameters = read_file_parameters()
    #nsp_p = ModelingToolkit.namespace_parameters(e4a) 

    dictionaryIndexValues = Dict{Any,Any}()

    for (k, v) in dict_parameters

        for (j, l) in parameters

            if (k == string(j))

                dictionaryIndexValues[v] = l
            end
        end

    end

    p = []
    index = 0
    open("parameters_default_values.txt", "w") do io
        print(io, "p = [")
        for key in sort(collect(keys(dictionaryIndexValues)))
            push!(p, dictionaryIndexValues[key])
            print(io, dictionaryIndexValues[key])
            index = index + 1
            if (index % 10 == 0)
                println(io, ",")
            else
                print(io, ",")
            end
            #println("$key => $(dictionaryNameInitialCondition[key])")
        end
        print(io, "]")
    end
    return dictionaryIndexValues
end

function comparison(e4a, parameters=_params)
    par_namespace = ModelingToolkit.namespace_parameters(e4a)
    println("keys ", typeof(keys(parameters)))
    index = 0


    for i in par_namespace

        c = string(i)
        c = split(c, "₊")
        c_2 = c[2]

        if (Symbol(c_2) in (keys(parameters)))
            delete!(parameters, Symbol(c_2))

        else
            println("not found", c_2)
        end
    end
    println(parameters)


end

function variable_index_name_value(_variables)

    dictionaryVariablesIndex = Dict{Any,Any}()
    dictionaryIndexVariables = Dict{Any,Any}()
    dictionaryIndexValue = Dict{Any,Any}()

    index = 0
    for (k, v) in _variables
        index = index + 1
        dictionaryVariablesIndex[string(k)] = index
        dictionaryIndexVariables[index] = string(k)
        dictionaryIndexValue[index] = v
    end
    return [dictionaryVariablesIndex, dictionaryIndexVariables, dictionaryIndexValue]
end

function parameters_index_namespace(e4a, p=_params)
    nsp_p = ModelingToolkit.namespace_parameters(e4a)
    dictionaryParametersIndex = Dict{String,Int64}()
    dictionaryIndexParameters = Dict{Int64,String}()
    index = 0
    for i in nsp_p
        c = string(i)
        c = split(c, "₊")
        c_2 = c[2]
        index = index + 1
        dictionaryParametersIndex[c_2] = index
        dictionaryIndexParameters[index] = c_2
    end

    return [dictionaryParametersIndex, dictionaryIndexParameters]
end



function params_index_name_value(_params)

    dictionaryParametersIndex = Dict{Any,Any}()
    dictionaryIndexParameters = Dict{Any,Any}()
    dictionaryIndexValue = Dict{Any,Any}()
    index = 0
    for (k, v) in _params
        index = index + 1
        dictionaryParametersIndex[string(k)] = index
        dictionaryIndexParameters[index] = string(k)
        dictionaryIndexValue[index] = v
    end

    return [dictionaryParametersIndex, dictionaryIndexParameters, dictionaryIndexValue]
end


"This function splits equations based on their characteristics add or smooth equation."

function equation_one_by_one()
    equations = split(text, "\n")
    chunk_array = []

    for e in equations

        if (string(e[1]) == "a")
            add_split_equation(e, chunk_array)

        elseif (string(e[1]) == "s")
            smooth_translate_diff_eq(e, chunk_array)
        else
            println("not found", e)
        end

    end
    return chunk_array
end

"This function writes all the equations in a txt file"
function write_txt()
    translated_equations = equation_one_by_one()

    open("all_equations.txt", "w") do io
        for line in translated_equations
            println(io, line)
        end
    end
end

"This function splits the add_equations in differential equations and not differential equations"

function add_split_equation(e, chunk_array)

    all_components = []
    e_split = getindex.(split.(e, "eqs, "), 2)
    e_split = chop(e_split, tail=1)
    push!(all_components, e_split)

    for e_split in all_components
        if (string(e_split[1]) == "D" && string(e_split[2]) == "(")

            #println("Go to write differential equations")
            write_differential_eq(e_split, chunk_array)

        end
        if (string(e_split[1]) != "D" && string(e_split[2]) != "(")
            #println("Go to another method")
            write_not_differential_eq(e_split, chunk_array)

        end


    end

    return [all_components, chunk_array]
end

"this function read the file txt with the dictionary of variables name and index"
function read_file_variables()
    vars = Dict{String,Int64}()
    f = open("dict_k_variables_v_index.txt")
    for line in eachline(f)
        var, val = split(line, "=>")
        var = chop(var, head=1, tail=2)
        vars[string(var)] = parse(Int64, val)
    end
    return vars
end

"this function read the file txt with the dictionary of parameters name and index"

function read_file_parameters()
    vars = Dict{String,Int64}()
    f = open("dict_k_parameters_v_index.txt")
    for line in eachline(f)
        var, val = split(line, "=>")
        var = chop(var, head=1, tail=2)
        vars[string(var)] = parse(Int64, val)
    end
    return vars

end

"This function translates the add equations that are not differential equations"
function write_not_differential_eq(e_split, chunk_array)
    #variable_index = variable_index_name_value(_variables)[1]
    variable_index = read_file_variables()
    parameter_index = read_file_parameters()

    chunk = split.(e_split, " ")
    new_chunk = ""

    for i in chunk


        if (length(i) >= 2)

            if (i in (keys(variable_index)))
                new_chunk = new_chunk * "u[" * string(variable_index[i]) * "]"
            end
            if (i in (keys(parameter_index)))

                new_chunk = new_chunk * "p[" * string(parameter_index[i]) * "]"
            end
            if (!(i in (keys(parameter_index))) && !(i in (keys(variable_index))))

                if (string(i[1]) == "e")
                    i_d = split.(i, "(")
                    if (string(i_d[2]) in (keys(variable_index)))
                        new_chunk = new_chunk * "exp(u[" * string(variable_index[i_d[2]]) * "] "
                    end
                    if (string(i_d[2]) in (keys(parameter_index)))

                        new_chunk = new_chunk * "exp(p[" * string(parameter_index[i_d[2]]) * "] "
                    end

                    if (string(first(i_d[2])) == "-")
                        c = i_d[2]
                        c = replace(c, "-" => "")

                        if (c in (keys(variable_index)))
                            new_chunk = new_chunk * "exp(-u[" * string(variable_index[c]) * "] "

                        elseif (c in (keys(parameter_index)))
                            new_chunk = new_chunk * "exp(-p[" * string(parameter_index[c]) * "] "

                        else
                            (!(c in (keys(parameter_index))) && !(c in (keys(parameter_index))))

                            new_chunk = new_chunk * string(i)
                        end
                    end
                elseif (string(i[1]) == "(" && string(last(i)) != ")")
                    i_n = replace(i, "(" => "")
                    if (i_n in (keys(variable_index)))
                        new_chunk = new_chunk * "(u[" * string(variable_index[i_n]) * "]"

                    elseif (i_n in (keys(parameter_index)))
                        new_chunk = new_chunk * "(p[" * string(parameter_index[i_n]) * "]"
                    else
                        new_chunk = new_chunk * " " * string(i) * " "
                    end
                elseif (string(last(i)) == ")" && string(i[1]) != "(")
                    c = string(last(i))
                    i_c = string(chop(i, tail=1))

                    if (string(i_c) in (keys(variable_index)))
                        new_chunk = new_chunk * "u[" * string(variable_index[i_c]) * "]" * c * " "
                    elseif (i_c in (keys(parameter_index)))
                        new_chunk = new_chunk * "p[" * string(parameter_index[i_c]) * "]" * c * " "
                    else
                        new_chunk = new_chunk * " " * string(i) * " "
                    end

                elseif (string(last(i)) == ",")
                    if (string(last(i, 2)) == "),")
                        c = string(last(i, 2))
                        i_c = string(chop(i, tail=2))
                        if (string(i_c) in (keys(variable_index)))
                            new_chunk = new_chunk * "u[" * string(variable_index[i_c]) * "]" * c * " "

                        elseif (i_c in (keys(parameter_index)))
                            new_chunk = new_chunk * "p[" * string(parameter_index[i_c]) * "]" * c * " "
                        else
                            new_chunk = new_chunk * " " * string(i) * " "
                        end
                    else
                        c = string(last(i))
                        i_c = string(chop(i, tail=1))
                        if (string(i_c) in (keys(variable_index)))
                            new_chunk = new_chunk * "u[" * string(variable_index[i_c]) * "]" * c * " "
                        elseif (i_c in (keys(parameter_index)))
                            new_chunk = new_chunk * "p[" * string(parameter_index[i_c]) * "]" * c * " "

                        else
                            new_chunk = new_chunk * " " * string(i) * " "
                        end
                    end
                elseif (string(i[1]) == "(" && string(last(i)) == ")")
                    i_c = replace(i, "(" => "", ")" => "")
                    if (i_c in (keys(variable_index)))
                        new_chunk = new_chunk * "(u[" * string(variable_index[i_c]) * "]) "
                    end
                    if (i_c in (keys(parameter_index)))
                        new_chunk = new_chunk * "(p[" * string(parameter_index[i_c]) * "]) "
                    end

                else
                    new_chunk = new_chunk * " " * string(i) * " "
                end

            end


        elseif (length(i) == 1)
            if (i == "~")
                new_chunk = new_chunk * " = "
            else
                new_chunk = new_chunk * " " * string(i) * " "
            end
        end

    end
    push!(chunk_array, new_chunk)
end


"This function translates the add equations that are differential equations"
function write_differential_eq(e_split, chunk_array)

    #variable_index = variables_index(e4a)[1]
    #parameter_index = parameters_index_namespace(e4a)[1]
    variable_index = read_file_variables()
    parameter_index = read_file_parameters()

    chunk = split.(e_split, " ")
    new_chunk = ""

    for i in chunk

        if (length(i) >= 2)

            if (string(i[1]) == "D" && string(i[2]) == "(")
                i_n = chop(i, head=2, tail=1)
                new_chunk = new_chunk * "du[" * string(variable_index[i_n]) * "]"
            end

            if (i in (keys(variable_index)))
                new_chunk = new_chunk * "u[" * string(variable_index[i]) * "]"
            end
            if (i in (keys(parameter_index)))
                new_chunk = new_chunk * "p[" * string(parameter_index[i]) * "]"
            end
            if (string(i[1]) == "(" && string(last(i)) != ")" && string(i[1]) != "D")
                i_n = replace(i, "(" => "")
                if (i_n in (keys(variable_index)))
                    new_chunk = new_chunk * "(u[" * string(variable_index[i_n]) * "]"

                elseif (i_n in (keys(parameter_index)))
                    new_chunk = new_chunk * "(p[" * string(parameter_index[i_n]) * "]"
                else
                    new_chunk = new_chunk * " " * string(i) * " "
                end
            end

            if (string(last(i)) == ")" && string(i[1]) != "(" && string(i[1]) != "D")
                c = string(last(i))
                i_c = string(chop(i, tail=1))

                if (string(i_c) in (keys(variable_index)))
                    new_chunk = new_chunk * "u[" * string(variable_index[i_c]) * "]" * c * " "
                elseif (i_c in (keys(parameter_index)))
                    new_chunk = new_chunk * "p[" * string(parameter_index[i_c]) * "]" * c * " "
                else
                    new_chunk = new_chunk * " " * string(i) * " "
                end
            end


        elseif (length(i) == 1)
            if (i == "~")
                new_chunk = new_chunk * " = "
            else
                new_chunk = new_chunk * " " * string(i) * " "
            end
        end

    end
    push!(chunk_array, new_chunk)
    return chunk_array
end

"This function translates the the smooth equations"
function smooth_translate_diff_eq(e, chunk_array)

    #variable_index = variables_index(e4a)[1]
    #parameter_index = parameters_index_namespace(e4a)[1]
    variable_index = read_file_variables()
    parameter_index = read_file_parameters()
    new_chunk1 = ""
    new_chunk2 = ""
    new_chunk3 = ""
    chunk_n = split.(e, "eqs, ")
    chunk = chop(chunk_n[2], tail=1)
    chunk = split.(chunk, ", ")

    #Define e translate du[variable]
    if (string(chunk[1]) in (keys(variable_index)))
        new_chunk1 = new_chunk1 * "du[" * string(variable_index[string(chunk[1])]) * "]" * " = "
    end

    #define and translate chunk 2 

    s2 = chunk[2]
    s2 = split.(s2, " ")
    length_s2 = length(s2)

    if (length_s2 == 1)
        if (string(chunk[2]) in (keys(variable_index)))
            new_chunk2 = new_chunk2 * "u[" * string(variable_index[string(chunk[2])]) * "]"
        end
    end

    if (length_s2 > 1)
        for i in s2

            if (length(i) > 1)
                if (string(first(i)) == "(")
                    i_n = replace(i, "(" => "")
                    if (i_n in (keys(variable_index)))
                        new_chunk2 = new_chunk2 * "(" * "u[" * string(variable_index[i_n]) * "]"

                    elseif (i_n in (keys(parameter_index)))
                        new_chunk2 = new_chunk2 * "(" * "p[" * string(parameter_index[i_n]) * "]"
                    else
                        new_chunk2 = new_chunk2 * " " * string(i) * " "
                    end


                elseif (string(last(i)) == ")")
                    i_n = chop(i, tail=1)
                    if (i_n in (keys(variable_index)))
                        new_chunk2 = new_chunk2 * "u[" * string(variable_index[string(i_n)]) * "])"

                    elseif (string(i_n) in (keys(parameter_index)))
                        new_chunk2 = new_chunk2 * "p[" * string(parameter_index[string(i_n)]) * "])"
                    else
                        new_chunk2 = new_chunk2 * " " * string(i) * " "
                    end

                elseif (string(last(i)) != ")" && string(first(i)) != "(" && string(first(i)) != "e")
                    if (i in (keys(variable_index)))
                        new_chunk2 = new_chunk2 * "u[" * string(variable_index[i]) * "]"
                    end
                    if (i in (keys(parameter_index)))
                        new_chunk2 = new_chunk2 * "p[" * string(parameter_index[i]) * "]"
                    end

                elseif (string(first(i)) == "e")
                    i_d = split.(i, "(")
                    if (string(first(i_d[2])) in (keys(variable_index)))
                        new_chunk2 = new_chunk2 * "exp(u[" * string(variable_index[i_d[2]]) * "] "
                    end
                    if (string(first(i_d[2])) in (keys(parameter_index)))
                        new_chunk2 = new_chunk2 * "exp(p[" * string(parameter_index[i_d[2]]) * "] "
                    end
                    if (string(first(i_d[2])) == "-")
                        c = i_d[2]
                        c = replace(c, "-" => "")
                        if (c in (keys(variable_index)))
                            new_chunk2 = new_chunk2 * "exp(-u[" * string(variable_index[c]) * "] "
                        end
                        if (c in (keys(parameter_index)))
                            new_chunk2 = new_chunk2 * "exp(-p[" * string(parameter_index[c]) * "] "
                        end
                    end

                else
                    new_chunk2 = new_chunk2 * " " * string(i) * " "
                end
            end
            if (length(i) == 1)
                if (i == "~")
                    new_chunk2 = new_chunk2 * " = "
                else
                    new_chunk2 = new_chunk2 * " " * string(i) * " "
                end

            end
        end
    end



    #define and translate chunk 3
    if (string(chunk[3]) in (keys(parameter_index)))
        new_chunk3 = new_chunk3 * " / p[" * string(parameter_index[string(chunk[3])]) * "]"

    end


    new_chunk_final = new_chunk1 * "(" * new_chunk2 * ")" * new_chunk3
    push!(chunk_array, new_chunk_final)

    return chunk_array


end

