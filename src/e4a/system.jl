include("functions.jl")
@register ramp(x, slope, startx, endx)

@variables t
D = Differential(t)

function earth4all(; name, params=_params, inits=_inits, tables=_tables, ranges=_ranges)
    ORDER = Int64(params[:ORDER])

    # dummy parameter and variable to be used because of bug of GlobalSensitivity.jl
    @parameters DUMMYP = 1.0 [description = "Dummy parameter to be used as first because of bug of GlobalSensitivity.jl"]
    @variables DUMMYV(t) = 0.0 [description = "Dummy variable to be used in first equation because of bug of GlobalSensitivity.jl"]
    eqs = [D(DUMMYV) ~ DUMMYP - log10(10^DUMMYP)]

    # climate
    @parameters CCCSt = params[:CCCSt] [description = "Cost of CCS Dollar/tCO2"]
    @parameters ERDN2OKF2022 = params[:ERDN2OKF2022] [description = "Extra rate of decline in N2O per kg fertilizer from 2022"]
    @parameters RDN2OKF = params[:RDN2OKF] [description = "Rate of decline in N2O per kg fertilizer 1/y"]
    @parameters KN2OKF1980 = params[:KN2OKF1980] [description = "kg N2O per kg fertilizer in 1980 "]
    @parameters N2OA1980 = params[:N2OA1980] [description = "N2O in atm in 1980 GtN2O"]
    @parameters MAT = params[:MAT] [description = "Mass of ATmosphere Zt"]
    @parameters LN2OA = params[:LN2OA] [description = "Life of N2O in atmosphere y"]
    @parameters GN2OPP = params[:GN2OPP] [description = "GtN2O per ppm"]
    @parameters ERDCH4KC2022 = params[:ERDCH4KC2022] [description = "Extra rate of decline in CH4 per kg crop after 2022 1/y"]
    @parameters RDCH4KC = params[:RDCH4KC] [description = "Rate of decline in CH4 per kg crop 1/y"]
    @parameters KCH4KC1980 = params[:KCH4KC1980] [description = "kg CH4 per kg crop in 1980"]
    @parameters CH4A1980 = params[:CH4A1980] [description = "CH4 in atm in 1980 GtCH4"]
    @parameters LCH4A = params[:LCH4A] [description = "Life of CH4 in atmosphere y"]
    @parameters GCH4PP = params[:GCH4PP] [description = "GtCH4 per ppm"]
    @parameters OBWA2022 = params[:OBWA2022] [description = "OBserved WArming in 2022 deg C"]
    @parameters SOWLCO2 = params[:SOWLCO2] [description = "sOWeoLoCO2>0"]
    @parameters LECO2A1980 = params[:LECO2A1980] [description = "Life of extra CO2 in atmosphere in 1980 y"]
    @parameters CO2A1850 = params[:CO2A1850] [description = "CO2 in atmosphere in 1850 GtCO2"]
    @parameters TCO2PTCH4 = params[:TCO2PTCH4] [description = "tCO2 per tCH4"]
    @parameters DACCO22100 = params[:DACCO22100] [description = "Direct Air Capture of CO2 in 2100 GtCO2/y"]
    @parameters GCO2PP = params[:GCO2PP] [description = "GtCO2 per ppm"]
    @parameters TCO2ETN2O = params[:TCO2ETN2O] [description = "tCO2e/tN2O"]
    @parameters TCO2ETCH4 = params[:TCO2ETCH4] [description = "tCO2e/tCH4"]
    @parameters TCO2ETCO2 = params[:TCO2ETCO2] [description = "tCO2e/tCO2"]
    @parameters ALGAV = params[:ALGAV] [description = "ALbedo global average"]
    @parameters ALIS = params[:ALIS] [description = "ALbedo Ice and snow"]
    @parameters GLSU = params[:GLSU] [description = "GLobal SUrface Mkm2"]
    @parameters ISCEGA1980 = params[:ISCEGA1980] [description = "Ice and snow cover excluding Greenland and Antarctica in 1980 Mkm2"]
    @parameters TRSS1980 = params[:TRSS1980] [description = "Transfer rate surface-space in 1980 1/y"]
    @parameters MRS1980 = params[:MRS1980] [description = "Melting rate surface in 1980 1/y"]
    @parameters WA1980 = params[:WA1980] [description = "Warming in 1980 deg C"]
    @parameters SVDR = params[:SVDR] [description = "Surface vs deep rate"]
    @parameters AI1980 = params[:AI1980] [description = "Amount of ice in 1980 MKm3"]
    @parameters TPM3I = params[:TPM3I] [description = "Ton per m3 ice"]
    @parameters HRMI = params[:HRMI] [description = "Heat required to melt ice kJ/kg"]
    @parameters OWWV = params[:OWWV] [description = "sOWeoWV>0"]
    @parameters WVC1980 = params[:WVC1980] [description = "Water Vapour Concentration in 1980 g/kg"]
    @parameters WVF1980 = params[:WVF1980] [description = "Water Vapour Feedback in 1980 W/m2"]
    @parameters WVWVF = params[:WVWVF] [description = "sWVeoWVF>0"]
    @parameters WFEH = params[:WFEH] [description = "Warming from Extra Heat deg/ZJ"]
    @parameters EH1980 = params[:EH1980] [description = "Extra Heat in 1980 ZJ"]
    @parameters TRSA1980 = params[:TRSA1980] [description = "Transfer rate surface-abyss in 1980 1/y"]
    @parameters PD = params[:PD] [description = "Perception delay y"]

    @variables KN2OEKF(t) [description = "kg N2O emission per kg fertiliser"]
    @variables MMN2OE(t) [description = "Man-made N2O emissions GtN2O/y"]
    @variables NN2OE(t) [description = "Natural N2O emissions GtN2O/y"]
    @variables N2OE(t) [description = "N2O emissions GtN2O/y"]
    @variables N2OC1980(t) [description = "N2O conc in 1980 ppm"]
    @variables N2OA(t) = inits[:N2OA] [description = "N2O in Atmosphere GtN2O"]
    @variables N2OBD(t) [description = "N2O BreakDown GtN2O/y"]
    @variables N2OCA(t) [description = "N2O concentration in atm ppm"]
    @variables N2OFPP(t) [description = "N2O forcing per ppm W/m2/ppm"]
    @variables FN2O(t) [description = "Forcing from N2O W/m2"]
    @variables KCH4EKC(t) [description = "kg CH4 emission per kg crop"]
    @variables MMCH4E(t) [description = "Man-made CH4 emissions GtCH4/y"]
    @variables NCH4E(t) [description = "Natural CH4 emissions GtCH4/y"]
    @variables CH4E(t) [description = "CH4 emissions GtCH4/y"]
    @variables CH4C1980(t) [description = "CH4 conc in 1980 ppm"]
    @variables CH4A(t) = inits[:CH4A] [description = "CH4 in Atmosphere GtCH4"]
    @variables CH4BD(t) [description = "CH4 BreakDown GtCH4/y"]
    @variables CH4CA(t) [description = "CH4 concentration in atm ppm"]
    @variables CH4FPP(t) [description = "CH4 forcing per ppm W/m2/ppm"]
    @variables FCH4(t) [description = "Forcing from CH4 W/m2"]
    @variables OWLCO2(t) [description = "OWeoLoCO2"]
    @variables LECO2A(t) [description = "Life of extra CO2 in atm y"]
    @variables CO2FCH4(t) [description = "CO2 from CH4 GtCO2/y"]
    @variables CO2AB(t) [description = "CO2 absorption GtCO2/y"]
    @variables CO2E(t) [description = "CO2 emissions GtCO2/y"]
    @variables CO2GDP(t) [description = "CO2 per GDP (kgCO2/Dollar)"]
    @variables CAC(t) [description = "Cost of air capture GDollar/y"]
    @variables DACCO2(t) [description = "Direct Air Capture of CO2 GtCO2/y"]
    @variables CO2A(t) = inits[:CO2A] [description = "CO2 in Atmosphere GtCO2"]
    @variables CO2CA(t) [description = "CO2 concentration in atm ppm"]
    @variables CO2FPP(t) [description = "CO2 forcing per ppm W/m2/ppm"]
    @variables FCO2(t) [description = "Forcing from CO2 W/m2"]
    @variables FOG(t) [description = "Forcing from other gases W/m2"]
    @variables MMF(t) [description = "Man-made Forcing W/m2"]
    @variables GHGE(t) [description = "GHG emissions GtCO2e/y"]
    @variables AL1980(t) [description = "ALbedo in 1980 (1)"]
    @variables AL(t) [description = "ALbedo (1)"]
    @variables TRHGS(t) [description = "Transfer rate for heat going to space 1/y"]
    @variables HTS(t) [description = "Heat to space ZJ/y"]
    @variables MRS(t) [description = "Melting rate surface 1/y"]
    @variables MRDI(t) [description = "Melting rate deep ice 1/y"]
    @variables ECIM(t) [description = "Extra cooling from ice melt ZJ/y"]
    @variables MEL(t) [description = "Melting Mha/y"]
    @variables ISCEGA(t) = inits[:ISCEGA] [description = "Ice and snow cover excl G&A Mkm2"]
    @variables ISC(t) [description = "Ice and snow cover Mha"]
    @variables WVC(t) [description = "Water Vapor Concentration g/kg"]
    @variables WVF(t) [description = "Water Vapour Feedback W/m2"]
    @variables TMMF(t) [description = "Total man-made forcing W/m2"]
    @variables EWFF(t) [description = "Extra Warming from forcing ZJ/y"]
    @variables OW(t) [description = "OBserved WArming deg C"]
    @variables REHE(t) [description = "Risk of extreme heat event (1)"]
    @variables PWA(t) = inits[:PWA] [description = "Perceived WArming deg C"]
    @variables TRHGA(t) [description = "Transfer rate for heat going to abyss 1/y"]
    @variables HDO(t) [description = "Heat to deep ocean ZJ/y"]
    @variables EHS(t) = inits[:EHS] [description = "Extra heat in surface ZJ"]

    # demand
    @parameters BITRW = params[:BITRW] [description = "Basic Income Tax Rate Workers"]
    @parameters EETF2022 = params[:EETF2022] [description = "Extra Empowerment Tax From 2022 (share of NI)"]
    @parameters EGTRF2022 = params[:EGTRF2022] [description = "Extra General Tax Rate From 2022"]
    @parameters EPTF2022 = params[:EPTF2022] [description = "Extra Pension Tax From 2022 (share of NI)"]
    @parameters ETGBW = params[:ETGBW] [description = "Extra Transfer of Govmnt Budget to Workers"]
    @parameters FETACPET = params[:FETACPET] [description = "Fraction of Extra TA Cost Paid by Extra Taxes"]
    @parameters FETPO = params[:FETPO] [description = "Fraction of Extra Taxes Paid by Owners"]
    @parameters FGDC2022 = params[:FGDC2022] [description = "Fraction of Govmnt Debt Cancelled in 2022 1/y"]
    @parameters FT1980 = params[:FT1980] [description = "Fraction Transferred in 1980"]
    @parameters GCF = params[:GCF] [description = "Govmnt Consumption Fraction"]
    @parameters GDDP = params[:GDDP] [description = "Govmnt DrawDown Period y"]
    @parameters GDPOSR = params[:GDPOSR] [description = "sGDPeoOSR<0"]
    @parameters GDPP1980 = params[:GDPP1980] [description = "GDP per Person in 1980 kDollar/p/y"]
    @parameters GEIC = params[:GEIC] [description = "Goal for Extra Income from Commons (share of NI)"]
    @parameters GITRO = params[:GITRO] [description = "Goal for Income Tax Rate Owners"]
    @parameters GPP = params[:GPP] [description = "Govmnt Payback Period y"]
    @parameters GSF2022 = params[:GSF2022] [description = "Govmnt Stimulus From 2022 (share of NI)"]
    @parameters INEQ1980 = params[:INEQ1980] [description = "INEQuality in 1980"]
    @parameters ITRO1980 = params[:ITRO1980] [description = "Income Tax Rate Owners in 1980"]
    @parameters ITRO2022 = params[:ITRO2022] [description = "Income Tax Rate Owners in 2022"]
    @parameters MATGF = params[:MATGF] [description = "Mult to Avoid Transient in Govmnt Finance"]
    @parameters MATWF = params[:MATWF] [description = "Mult to Avoid Transient in Worker Finance"]
    @parameters MGDB = params[:MGDB] [description = "Max Govmnt Debt Burden y"]
    @parameters MWDB = params[:MWDB] [description = "Max Workers Debt Burden y"]
    @parameters OSF1980 = params[:OSF1980] [description = "Owner Savings Fraction in 1980"]
    @parameters STR = params[:STR] [description = "Sales Tax Rate"]
    @parameters TAB = params[:TAB] [description = "Time to Adjust Budget y"]
    @parameters TAOC = params[:TAOC] [description = "Time to Adjust Owner Consumption y"]
    @parameters TAWC = params[:TAWC] [description = "Time to Adjust Worker Consumption y"]
    @parameters TINT = params[:TINT] [description = "Time to Implement New Taxes y"]
    @parameters WCF = params[:WCF] [description = "Worker Consumption Fraction"]
    @parameters WDP = params[:WDP] [description = "Worker Drawdown Period y"]
    @parameters WPP = params[:WPP] [description = "Worker Payback Period y"]

    @variables BCIL(t) [description = "Bank Cash Inflow from Lending GDollar/y"]
    @variables BCISNI(t) [description = "Bank Cash Inflow as Share of NI (1)"]
    @variables BITRO(t) [description = "Basic Income Tax Rate Owners (1)"]
    @variables CANCD(t) [description = "CANCellation of Debt GDollar/y"]
    @variables CD(t) [description = "Consumption Demand GDollar/y"]
    @variables CFGB(t) [description = "Cash Flow from Govmnt to Banks GDollar/y"]
    @variables CFWB(t) [description = "Cash Flow from Workers to Banks GDollar/y"]
    @variables CONTR(t) [description = "CONTRol: (C+G+S)/NI = 1"]
    @variables CPP(t) [description = "Consumption Per Person GDollar/y"]
    @variables CSGDP(t) [description = "Consumption Share of GDP (1)"]
    @variables EGTF2022(t) [description = "Extra General Tax From 2022 Gdollar/y"]
    @variables ETF2022(t) = inits[:ETF2022] [description = "Extra Taxes From 2022 GDollar/y"]
    @variables ETTAF2022(t) [description = "Extra Taxes for TAs From 2022 GDollar/y"]
    @variables FGBW(t) = inits[:FGBW] [description = "Fraction of Govmnt Budget to Workers (1)"]
    @variables GCIN(t) [description = "Govmnt Cash INflow GDollar/y"]
    @variables GD(t) = 28087 * params[:MATGF] [description = "Govmnt Debt Gdollar"]
    @variables GDB(t) [description = "Govmnt Debt Burden y"]
    @variables GETF2022(t) [description = "Goal for Extra Taxes From 2022 GDollar/y"]
    @variables GFGBW(t) [description = "Goal for Fraction of Govmnt Budget to Workers (1)"]
    @variables GFSNI(t) [description = "Govmnt Finance as Share of NI (1)"]
    @variables GGI(t) [description = "Govmnt Gross Income GDollar/y"]
    @variables GGIS(t) [description = "Govmnt Gross Income (as Share of NI)"]
    @variables GIC(t) [description = "Govmnt Interest Cost GDollar/y"]
    @variables GIPC(t) [description = "Govmnt Investment in Public Capacity GDollar/y"]
    @variables GND(t) [description = "Govmnt New Debt GDollar/y"]
    @variables GNI(t) = inits[:GNI] [description = "Govmnt Net Income GDollar/y"]
    @variables GNISNI(t) [description = "Govmnt Net Income as Share of NI (1)"]
    @variables GP(t) [description = "Govmnt Payback GDollar/y"]
    @variables GPU(t) [description = "Govmnt PUrchases GDollar/y"]
    @variables GS(t) [description = "Govmnt Spending GDollar/y"]
    @variables GSGDP(t) [description = "Govmnt Share of GDP (1)"]
    @variables IC2022(t) [description = "Income from Commons from 2022 GDollar/y"]
    @variables INEQ(t) [description = "INEQuality (1)"]
    @variables INEQI(t) [description = "INEQuality Index (1980=1)"]
    @variables ITO(t) [description = "Income Tax Owners (1)"]
    @variables ITW(t) [description = "Income Tax Workers (1)"]
    @variables MGD(t) [description = "Max Govmnt Debt GDollar"]
    @variables MWD(t) [description = "Max Workers Debt GDollar"]
    @variables OC(t) [description = "Owner Consumption GDollar/y"]
    @variables OCF(t) [description = "Owner Consumptin Fraction (1)"]
    @variables OCIN(t) [description = "Owner Cash INflow GDollar/y"]
    @variables OI(t) [description = "Owner Income GDollar/y"]
    @variables OOIAT(t) [description = "Owner Operating Income After Tax GDollar/y"]
    @variables OS(t) [description = "Owner Savings GDollar/y"]
    @variables OSF(t) [description = "Owner Savings Fraction (1)"]
    @variables OT(t) [description = "Owner Taxes GDollar/y"]
    @variables OTR(t) [description = "Owner Tax Rate (1)"]
    @variables PGCIN(t) = inits[:PGCIN] [description = "Permanent Govmnt Cash INflow GDollar/y"]
    @variables POCI(t) = inits[:POCI] [description = "Permanent Owner Cash Inflow GDollar/y"]
    @variables PWCIN(t) = inits[:PWCIN] [description = "Permanent Worker Cash INflow GDollar/y"]
    @variables SSGDP(t) [description = "Savings Share of GDP (1)"]
    @variables ST(t) [description = "Sales Tax GDollar/y"]
    @variables STO(t) [description = "Sales Tax Owners GDollar/y"]
    @variables STW(t) [description = "Sales Tax Workers GDollar/y"]
    @variables TP(t) [description = "Transfer Payments GDollar/y"]
    @variables TPP(t) [description = "Total Purchasing Power GDollar/y"]
    @variables TS(t) [description = "Total Savings GDollar/y"]
    @variables WCD(t) [description = "Worker consumption demand GDollar/y"]
    @variables WCIN(t) [description = "Worker Cash INflow GDollar/y"]
    @variables WD(t) = 18992 * params[:MATWF] [description = "Workers Debt GDollar"]
    @variables WDB(t) [description = "Worker Debt Burden y"]
    @variables WDI(t) [description = "Worker Disposable Income kDollar/p/y"]
    @variables WFCSI(t) [description = "Worker Finance Cost as Share of Income (1)"]
    @variables WI(t) [description = "Worker Income GDollar/y"]
    @variables WIAT(t) [description = "Worker Income After Tax GDollar/y"]
    @variables WIC(t) [description = "Worker Interest Cost GDollar/y"]
    @variables WND(t) [description = "Workers New Debt GDollar/y"]
    @variables WP(t) [description = "Workers Payback GDollar/y"]
    @variables WS(t) [description = "Worker Savings GDollar/y"]
    @variables WT(t) [description = "Worker Taxes GDollar/y"]
    @variables WTR(t) [description = "Worker Tax Rate (1)"]

    # energy
    @parameters MNFCO2PP = params[:MNFCO2PP] [description = "Max non-fossil CO2 per person tCO2/p/y"]
    @parameters FCO2SCCS2022 = params[:FCO2SCCS2022] [description = "Fraction of CO2-sources with CCS in 2022"]
    @parameters GFCO2SCCS = params[:GFCO2SCCS] [description = "Goal for fraction of CO2-sources with CCS"]
    @parameters CCCSt = params[:CCCSt] [description = "Cost of CCS Dollar/tCO2"]
    @parameters ROCTCO2PT = params[:ROCTCO2PT] [description = "ROC in tCO2 per toe 1/y"]
    @parameters EROCEPA2022 = params[:EROCEPA2022] [description = "Extra ROC in energy productivity after 2022 1/y"]
    @parameters NIEE = params[:NIEE] [description = "Normal increase in energy efficiency 1/y"]
    @parameters GFNE = params[:GFNE] [description = "Goal for fraction new electrification"]
    @parameters FNE2022 = params[:FNE2022] [description = "Fraction new electrification in 2022"]
    @parameters FNE1980 = params[:FNE1980] [description = "Fraction new electrification in 1980"]
    @parameters EUEPRUNEFF = params[:EUEPRUNEFF] [description = "Extra use of electricity per reduced use of non-el FF MWh/toe"]
    @parameters ECRUNEFF = params[:ECRUNEFF] [description = "Extra cost per reduced use og non-el FF Dollar/toe"]
    @parameters GREF = params[:GREF] [description = "Goal for renewable el fraction"]
    @parameters REFF2022 = params[:REFF2022] [description = "Renewable el fraction in 2022"]
    @parameters REFF1980 = params[:REFF1980] [description = "Renewable el fraction in 1980"]
    @parameters RCUT = params[:RCUT] [description = "Renewable capacity up-time kh/y"]
    @parameters RECT = params[:RECT] [description = "Renewable el contruction time y"]
    @parameters LREC = params[:LREC] [description = "Life of renewable el capacity y"]
    @parameters SWC1980 = params[:SWC1980] [description = "Sun and wind capacity in 1980 GW"]
    @parameters CRDSWC = params[:CRDSWC] [description = "Cost reduction per dubling of sun and wind capacity"]
    @parameters CAPEXRE1980 = params[:CAPEXRE1980] [description = "CAPEX of renewable el in 1980 Dollar/W"]
    @parameters CAPEXFED = params[:CAPEXFED] [description = "CAPEX fossil el Dollar/W"]
    @parameters OPEXRED = params[:OPEXRED] [description = "OPEX renewable el Dollar/kWh"]
    @parameters OPEXFED = params[:OPEXFED] [description = "OPEX fossil el Dollar/kWh"]
    @parameters CNED = params[:CNED] [description = "Cost of Nuclear El Dollar/kWh"]
    @parameters FREH = params[:FREH] [description = "Fraction of renewable electricity to hydrogen"]
    @parameters KWEPKGH2 = params[:KWEPKGH2] [description = "kWh electricity per kg of hydrogen"]
    @parameters TPTH2 = params[:TPTH2] [description = "toe per tH2"]
    @parameters BEM = params[:BEM] [description = "Biomass energy Mtoe/y"]
    @parameters EFPP = params[:EFPP] [description = "Efficiency of fossil power plant TWh-el/TWh-heat"]
    @parameters TWHPEJCE = params[:TWHPEJCE] [description = "TWh-heat per EJ - calorific equivalent"]
    @parameters MTPEJCE = params[:MTPEJCE] [description = "Mtoe per EJ - calorific equivalent"]
    @parameters EKHPY = params[:EKHPY] [description = "8 khours per year"]
    @parameters FECCT = params[:FECCT] [description = "Fossil el cap construction time y"]
    @parameters NLFEC = params[:NLFEC] [description = "Normal life of fossil el capacity y"]
    @parameters sFCUTLOFC = params[:sFCUTLOFC] [description = "sFCUTeoLOFC>0"]
    @parameters NCUT = params[:NCUT] [description = "Nuclear capacity up-time kh/y"]
    @parameters TCE = params[:TCE] [description = "Traditional cost of electricity Dollar/kWh"]
    @parameters AFMCM = params[:AFMCM] [description = "Adjustment factor to make cost match 1980 - 2022"]
    @parameters TCFFFNEU = params[:TCFFFNEU] [description = "Traditional cost of fossil fuels for non-el use Dollar/toe"]
    @parameters TC = params[:TC] [description = "Transmission cost Dollar/kWh"]

    @variables CNEL(t) [description = "Cost of Nuclear ELectricity Gdollar/y"]
    @variables CO2NFIP(t) [description = "CO2 from non-fossil industrial processes GtCO2/y"]
    @variables FCO2SCCS(t) [description = "Fraction of CO2-sources with CCS (1)"]
    @variables NFCO2PP(t) [description = "Non-fossil CO2 per person tCO2/p/y"]
    @variables CO2EI(t) [description = "CO2 from energy and industry GtCO2/y"]
    @variables CO2EP(t) [description = "CO2 from energy production GtCO2/y"]
    @variables CO2EMPP(t) [description = "CO2 EMissions per person tCO2/y"]
    @variables CCCSG(t) [description = "Cost of CCS GDollar/y"]
    @variables ICCSC(t) [description = "Installed CCS capacity GtCO2/y"]
    @variables TCO2PT(t) [description = "tCO2 per toe"]
    @variables EEPI2022(t) = inits[:EEPI2022] [description = "Extra energy productivity index 2022=1"]
    @variables IEEPI(t) [description = "Increase in extra energy productivity index 1/y"]
    @variables TPPUEBEE(t) [description = "Traditional per person use of electricity before EE MWh/p/y"]
    @variables TPPUFFNEUBEE(t) [description = "Traditional per person use of fossil fuels for non-el-use before EE toe/p/y"]
    @variables DEBNE(t) [description = "Demand for electricity before NE TWh/y"]
    @variables DFFNEUBNE(t) [description = "Demand for fossil fuel for non-el use before NE Mtoe/y"]
    @variables FNE(t) [description = "Fraction new electrification (1)"]
    @variables ERDNEFFFNE(t) [description = "Extra reduction in demand for non-el fossil fuel from NE Mtoe/y"]
    @variables CNE(t) [description = "Cost of new electrification GDollar/y"]
    @variables EIDEFNE(t) [description = "Extra increase in demand for electricity from NE TWh/y"]
    @variables DFFFNEU(t) [description = "Demand for fossil fuel for non-el use Mtoe/y"]
    @variables UFF(t) [description = "Use of fossil fuels Mtoe/y"]
    @variables DE(t) [description = "Demand for electricity TWh/y"]
    @variables DRES(t) [description = "Desired renewable electricity share (1)"]
    @variables DSRE(t) [description = "Desired supply of renewable electricity TWh/y"]
    @variables DREC(t) [description = "Desired renewable el capacity GW"]
    @variables DRECC(t) [description = "Desired renewable el capacity change GW"]
    @variables REC(t) = inits[:REC] [description = "Renewable electricity capacity GW"]
    @variables AREC(t) [description = "Addition of renewable el capacity GW/y"]
    @variables DIREC(t) [description = "Discard of renewable el capacity GW/y"]
    @variables ASWC(t) [description = "Addition of sun and wind capacity GW/y"]
    @variables ACSWCF1980(t) = inits[:ACSWCF1980] [description = "ACcumulated sun and wind capacity from 1980 GW"]
    @variables NDSWC(t) [description = "Number of doublings in sun and wind capacity (1)"]
    @variables CISWC(t) [description = "Cost index for sun and wind capacity (1)"]
    @variables CAPEXRED(t) [description = "CAPEX renewable el dollar/W"]
    @variables CAPEXREG(t) [description = "CAPEX renewable el Gdollar/y"]
    @variables OPEXREG(t) [description = "OPEX renewable el Gdollar/y"]
    @variables CRE(t) [description = "Cost of renewable electricity GDollar/y"]
    @variables CAPEXFEG(t) [description = "CAPEX fossil el GDollar/y"]
    @variables OPEXFEG(t) [description = "OPEX fossil el GDollar/y"]
    @variables CFE(t) [description = "Cost of fossil electricity GDollar/y"]
    @variables CEL(t) [description = "Cost of electricity GDollar/y"]
    @variables REP(t) [description = "Renewable electricity production TWh/y"]
    @variables GHMH2(t) [description = "Green hydrogen MtH2/y"]
    @variables GHMt(t) [description = "Green hydrogen Mtoe/y"]
    @variables RHP(t) [description = "Renewable heat production Mtoe/y"]
    @variables TWEPEJEE(t) [description = "TWh-el per EJ - engineering equivalent"]
    @variables IIASAREP(t) [description = "IIASA Renewable energy production EJ/yr"]
    @variables FTWEPMt(t) [description = "4 TWh-el per Mtoe"]
    @variables IIASAFEP(t) [description = "IIASA Fossil energy production EJ/yr"]
    @variables LCEP(t) [description = "Low-carbon el production TWh/y"]
    @variables DFE(t) [description = "Demand for fossil electricity TWh/y"]
    @variables DFEC(t) [description = "Desired fossil el capacity GW"]
    @variables DFECC(t) [description = "Desired fossil el capacity change GW/y"]
    @variables AFEC(t) [description = "Addition of fossil el capacity GW/y"]
    @variables FEC(t) = inits[:FEC] [description = "Fossil electricity capacity GW"]
    @variables LFEC(t) [description = "Life of fossil el capacity y"]
    @variables DIFEC(t) [description = "DIscard of Fossil El Capacity GW/y"]
    @variables FCUT(t) [description = "Fossil capacity up-time kh/y"]
    @variables FCUTLOFC(t) [description = "FCUTeoLOFC (1)"]
    @variables FEP(t) [description = "Fossil Electricity Production TWh/y"]
    @variables NC(t) [description = "Nuclear capacity GW"]
    @variables NEP(t) [description = "Nuclear electricity production TWh/y"]
    @variables EP(t) [description = "Electricity production TWh/y"]
    @variables ELB(t) [description = "ELectricity balance (1)"]
    @variables FFPNE(t) [description = "Fraction fossil plus nuclear electricity (1)"]
    @variables EU(t) [description = "Energy use Mtoe/y"]
    @variables EUPP(t) [description = "Energy use per person toe/p/y"]
    @variables FFE(t) [description = "Fossil fuels for electricity Mtoe/y"]
    @variables TCEG(t) [description = "Traditional cost of electricity GDollar/y"]
    @variables TCFFFNEUG(t) [description = "Traditional cost of fossil fuel for non-el use Gdollar/y"]
    @variables CFFFNEU(t) [description = "Cost of Fossil Fuel For Non-El Use Gdollar/y"]
    @variables CG(t) [description = "Cost of grid GDollar/y"]
    @variables TGC(t) [description = "Traditional grid cost GDollar/y"]
    @variables TCEN(t) [description = "Traditional Cost of ENergy Gdollar/y"]
    @variables TCENSGDP(t) [description = "Traditional cost of energy as share of GDP (1)"]
    @variables CE(t) [description = "Cost of energy GDollar/y"]
    @variables RECTEC(t) [description = "Ratio of Energy cost to Trad Energy cost (1)"]
    @variables CESGDP(t) [description = "Cost of energy as share of GDP (1)"]
    @variables ECETSGDP(t) [description = "Extra cost of Energy Turnaround as share of GDP (1)"]

    # finance
    @parameters FSRT = params[:FSRT] [description = "Financial Sector Response Time y"]
    @parameters GRCR = params[:GRCR] [description = "sGReoCR<0: Growth Rate effect on Credit Risk"]
    @parameters IEFT = params[:IEFT] [description = "Inflation Expectation Formation Time y"]
    @parameters INSR = params[:INSR] [description = "sINeoSR>0: INflation effect on Signal Rate"]
    @parameters IPTCB = params[:IPTCB] [description = "Inflation Perception Time CB y"]
    @parameters IT = params[:IT] [description = "Inflation Target 1/y"]
    @parameters NBOM = params[:NBOM] [description = "Normal Bank Operating Margin 1/y"]
    @parameters NBBM = params[:NBBM] [description = "Normal Basic Bank Margin 1/y"]
    @parameters NSR = params[:NSR] [description = "Normal Signal Rate 1/y"]
    @parameters SRAT = params[:SRAT] [description = "Signal Rate Adjustment Time y"]
    @parameters UPTCB = params[:UPTCB] [description = "Unemployment Perception Time CB y"]
    @parameters UNSR = params[:UNSR] [description = "sUNeoSR<0: UNemployment effect on Signal Rate"]
    @parameters UT = params[:UT] [description = "Unemployment Target"]

    @variables CBC(t) [description = "Corporate Borrowing Cost 1/y"]
    @variables CBSR(t) = inits[:CBSR] [description = "Central Bank Signal Rate 1/y"]
    @variables CCSD(t) = inits[:CCSD] [description = "Cost of Capital for Secured Debt 1/y"]
    @variables CSR(t) [description = "Change in Signal Rate 1/yy"]
    @variables GBC(t) [description = "Govmnt Borrowing Cost 1/y"]
    @variables CBC1980(t) [description = "Corporate Borrowing Cost in 1980 1/y"]
    @variables ELTI(t) = inits[:ELTI] [description = "Expected Long Term Inflation 1/y"]
    @variables ISR(t) [description = "Indicated Signal Rate 1/y"]
    @variables NCCR(t) [description = "Normal Corporate Credit Risk 1/y"]
    @variables PI(t) = inits[:PI] [description = "Perceived Inflation CB 1/y"]
    @variables PU(t) = inits[:PU] [description = "Perceived Unemployment CB (1)"]
    @variables TGIR(t) [description = "10-yr Govmnt Interest Rate 1/y"]
    @variables TIR(t) [description = "3m Interest Rate 1/y"]
    @variables WBC(t) [description = "Worker Borrowing Cost 1/y"]

    # foodland
    @parameters AFGDP = params[:AFGDP] [description = "Agriculture as Fraction of GDP"]
    @parameters CBECLE = params[:CBECLE] [description = "sFBeoCLE<0: Crob Balance Effect on CropLand Expansion"]
    @parameters CO2ARA = params[:CO2ARA] [description = "CO2 Absorbed in Reg Ag tCO2/ha/y"]
    @parameters CO2C2022 = params[:CO2C2022] [description = "Climate.CO2 concentration in 2022 ppm"]
    @parameters CO2CEACY = params[:CO2CEACY] [description = "sCO2CeoACY>0: CO2 Concentration Effect on Average Crop Yeld"]
    @parameters CO2RHFC = params[:CO2RHFC] [description = "CO2 Release per Ha of Forest Cut tCO2/ha"]
    @parameters CRDRA = params[:CRDRA] [description = "Cost Reduction per Doubling in Regenerative Agriculture"]
    @parameters CTF = params[:CTF] [description = "Cost per Ton Fertilizer dollar/t"]
    @parameters CYRA = params[:CYRA] [description = "Crop Yield in Reg Ag t-crop/ha/y"]
    @parameters DRC = params[:DRC] [description = "Desired Reserve Capacity"]
    @parameters ECRA22 = params[:ECRA22] [description = "Extra Cost of Reg Ag in 2022 dollar/ha/y"]
    @parameters EGB22 = params[:EGB22] [description = "Experience Gained Before 2022 Mha"]
    @parameters EROCFSP = params[:EROCFSP] [description = "Extra ROC in Food Sector Productivity from 2022 1/y"]
    @parameters FCG = params[:FCG] [description = "Fraction Cleared for Grazing"]
    @parameters FF80 = inits[:CRLA] * params[:FU80] [description = "Food Footprint in 1980"]
    @parameters FFLREOGRRM = params[:FFLREOGRRM] [description = "sFFLReoOGRR<0: Fraction Forestry Land Remaining Effect on Old Growth Removal Rate Multiplier"]
    @parameters FU80 = params[:FU80] [description = "Fertilizer Use in 1980 Mt/y"]
    @parameters FUELER = params[:FUELER] [description = "sFUeoLER>0: Fertilizer Use Effect on Land Erosion Rate"]
    @parameters FUESQ = params[:FUESQ] [description = "sFUeoSQ<0: Fertilizer Use Effect on Soil Quality"]
    @parameters GCWR = params[:GCWR] [description = "Goal for Crop Waste Reduction"]
    @parameters GFNRM = params[:GFNRM] [description = "Goal for Fraction New Red Meat"]
    @parameters GFRA = params[:GFRA] [description = "Goal for fraction regenerative agriculture"]
    @parameters KCKRM = params[:KCKRM] [description = "Kg-Crop per Kg-Red-Meat"]
    @parameters LER80 = params[:LER80] [description = "Land Erosion Rate in 1980 1/y"]
    @parameters MFAM = params[:MFAM] [description = "Max Forest Absorption Multiplier"]
    @parameters OGRR80 = params[:OGRR80] [description = "OGRR in 1980 1/y"]
    @parameters OW2022 = params[:OW2022] [description = "Climate.Observed Warming in 2022 deg C"]
    @parameters OWEACY = params[:OWEACY] [description = "sOWeoACY<0: Observed Warming Effect on Average Crop Yeld"]
    @parameters ROCFP = params[:ROCFP] [description = "ROC in Fertilizer Productivity 1/y"]
    @parameters ROCFSP = params[:ROCFSP] [description = "ROC in Food Sector Productivity 1/y"]
    @parameters SFU = params[:SFU] [description = "Sustainable Fertiliser Use kgN/ha/y"]
    @parameters SSP2LMA = params[:SSP2LMA] [description = "SSP2 land management action from 2022?"]
    @parameters TCTB = params[:TCTB] [description = "Ton Crops per Toe Biofuel"]
    @parameters TFFLR = params[:TFFLR] [description = "Threshold FFLR"]
    @parameters UDT = params[:UDT] [description = "Urban Development Time y"]
    @parameters ULP = params[:ULP] [description = "Urban Land per Population ha/p"]

    @variables ACY(t) [description = "Average Crop Yield t-crop/ha/y"]
    @variables ALFL(t) [description = "Acceptable Loss of Forestry Land (1)"]
    @variables AFSRA(t) [description = "Amount of Fertilizer Saved in Reg Ag kgN/ha/y"]
    @variables BALA(t) = inits[:BALA] [description = "BArren LAnd Mha"]
    @variables BIUS(t) [description = "BIofuels USe Mtoe/y"]
    @variables CEM(t) [description = "Cropland Expansion Multiplier (1)"]
    @variables CIRA(t) [description = "Cost Index for Regenerative Agriculture (1)"]
    @variables CO2AFL(t) [description = "CO2 Absorption in Forestry Land GtCO2/y"]
    @variables CO2AFLH(t) [description = "CO2 Absorption in Forestry Land tCO2/ha/y"]
    @variables CO2ELULUC(t) [description = "CO2 Emissions from LULUC GtCO2/y"]
    @variables CO2ELY(t) [description = "CO2 Effect on Land Yield (1)"]
    @variables CO2RFC(t) [description = "CO2 Release from Forest Cut GtCO2/y"]
    @variables COFE(t) [description = "COst of FErtilizer Gdollar/y"]
    @variables COFO(t) [description = "COst of FOod Gdollar/y"]
    @variables CRA(t) [description = "Cost of Regenerative Agriculture Gdollar/y"]
    @variables CRBA(t) [description = "CRop BAlance (1)"]
    @variables CRBI(t) [description = "CRops for BIofuel Mt-crop/y"]
    @variables CRDE(t) [description = "CRop DEmand Mt-crop/y"]
    @variables CREX(t) [description = "CRopland EXpansion Mha/y"]
    @variables CREXR(t) [description = "CRopland EXpansion Rate 1/y"]
    @variables CRLA(t) = inits[:CRLA] [description = "CRopLAnd Mha"]
    @variables CRLO(t) [description = "CRopland LOss Mha/y"]
    @variables CRSU(t) [description = "CRop SUpply (after 20 % waste) Mt-crop/y"]
    @variables CRUS(t) [description = "CRop USe Mt/y"]
    @variables CRUSP(t) [description = "CRop USe per Person t-crop/p/y"]
    @variables CSQCA(t) [description = "Change in Soil Quality in Conv Ag t-crop/ha/y/y"]
    @variables CSRA(t) [description = "Crop Supply Reg Ag Mt-crop/y"]
    @variables CWR(t) [description = "Crop Waste Reduction (1)"]
    @variables DCS(t) [description = "Desired Crop Supply Mt-crop/y"]
    @variables DCSCA(t) [description = "Desired Crop Supply Conv Ag Mt-crop/y"]
    @variables DCYCA(t) [description = "Desired Crop Yield in Conv Ag t-crop/ha/y"]
    @variables DRM(t) [description = "Demand for Red Meat Mt-red-meat/y"]
    @variables DRMP(t) [description = "Demand for Red Meat per Person kg-red-meat/p/y"]
    @variables ECFT(t) [description = "Extra Cost of Food Turnaround Gdollar/y"]
    @variables ECFTSGDP(t) [description = "Extra Cost of Food Turnaround as Share of GDP (1)"]
    @variables ECO2ARA(t) [description = "Extra CO2 Absorption in Reg Ag GtCO2/y"]
    @variables ECRA(t) [description = "Extra Cost of Reg Ag dollar/ha/y"]
    @variables FAM(t) [description = "Forest Absorption Multipler (1)"]
    @variables FCR(t) [description = "Fertilizer Cost Reduction Gdollar/y"]
    @variables FEER(t) [description = "Fertilizer Effect on Erosion Rate (1)"]
    @variables FERM(t) [description = "FEed for Red Meat Mt-crop/y"]
    @variables FEUS(t) [description = "FErtilizer USe Mt/y"]
    @variables FFI(t) [description = "Food Footprint Index (1980=1)"]
    @variables FFLR(t) [description = "Fraction Forestry Land Remaining (1)"]
    @variables FFLREOGRR(t) [description = "FFLReoOGRR"]
    @variables FNRM(t) [description = "Fraction New Red Meat (1)"]
    @variables FOFO(t) [description = "FOod FOotprint"]
    @variables FOLA(t) = params[:INITFOLA] [description = "FOrestry LAnd Mha"]
    @variables FPI(t) [description = "Fertilizer Productivity Index (1980=1)"]
    @variables FRA(t) [description = "Fraction Regenerative Agriculture (1)"]
    @variables FSPI(t) [description = "Food Sector Productivity Index (1980=1)"]
    @variables FUCA(t) [description = "Fertilizer Use in Conv Ag kgN/ha/y"]
    @variables FUP(t) [description = "FErtilizer USe per Person kg/p/y"]
    @variables GLY(t) [description = "Grazing Land Yield kg-red-meat/ha/y"]
    @variables GLY80(t) [description = "Grazing Land Yied in 1980 kg-red-meat/ha/y"]
    @variables GRLA(t) = inits[:GRLA] [description = "GRazing LAnd Mha"]
    @variables IUL(t) [description = "Indicated Urban Land Mha"]
    @variables LERM(t) [description = "Land ERosion Multiplier (1)"]
    @variables LER(t) [description = "Land Erosion Rate 1/y"]
    @variables LFL(t) [description = "Loss of Forest Land Mha/y"]
    @variables LOCR(t) [description = "LOss of CRopland Mha/y"]
    @variables NDRA(t) [description = "Number of Doublings in Reg Ag (1)"]
    @variables NFL(t) [description = "New Forestry Land Mha/y"]
    @variables NGL(t) [description = "New Grazing Land Mha/y"]
    @variables OGFA(t) = inits[:OGFA] [description = "Old Growth Forest Area Mha 1"]
    @variables OGRE(t) [description = "Old Growth Removal Mha/y"]
    @variables OGRR(t) [description = "Old Growth Removal Rate 1/y"]
    @variables OGRRM(t) [description = "Old Growth Removal Rate Multiplier (1)"]
    @variables PCB(t) [description = "Perceived Crop Balance (1)"]
    @variables PRMGL(t) [description = "Potential Red Meat from Grazing Land Mt-red-meat/y"]
    @variables RAA(t) [description = "Regenerative Agriculture Area Mha"]
    @variables RMF(t) [description = "Red Meat from Feedlots Mt-red-meat/y"]
    @variables RMGL(t) [description = "Red Meat from Grazing Land Mt-red-meat/y"]
    @variables RMSP(t) [description = "Red meat Supply per Person kg-red-meat/p/y"]
    @variables ROCSQCA(t) [description = "ROC in Soil Quality in Conv Ag 1/y"]
    @variables SQICA(t) = inits[:SQICA] [description = "Soil Quality Index in Conv Ag (1980=1)"]
    @variables TFA(t) [description = "Total Forest Area Mha"]
    @variables TFUCA(t) [description = "Traditional Fertilizer Use in Conv Ag kgN/ha/y"]
    @variables TUC(t) [description = "Traditional Use of Crops Mt/y"]
    @variables TUCERM(t) [description = "Traditional Use of Crops Ex Red Meat Mt/y"]
    @variables TUCERMP(t) [description = "Traditional Use of Crops Ex Red Meat per Person kg-crop/p/y"]
    @variables TUCP(t) [description = "Traditional Use of Crops per Person kg-crop/p/y"]
    @variables TUFRM(t) [description = "Traditional Use of Feed for Red Meat Mt-crop/y"]
    @variables TURMP(t) [description = "Traditional Use of Red Meat per Person kg-red-meat/p/y"]
    @variables UREX(t) [description = "URban EXpansion Mha/y"]
    @variables URLA(t) = inits[:URLA] [description = "URban LAnd Mha"]
    @variables WELY(t) [description = "Warming Effect on Land Yield (1)"]

    # inventory
    @parameters DAT = params[:DAT] [description = "Demand Adjustment Time y"]
    @parameters DDI1980 = params[:DDI1980] [description = "DDI in 1980 y"]
    @parameters DIC = params[:DIC] [description = "Desired Inventory Coverage y"]
    @parameters DRI = params[:DRI] [description = "Desired Relative Inventory"]
    @parameters ICPT = params[:ICPT] [description = "Inventory Coverage Perception Time y"]
    @parameters MRIWI = params[:MRIWI] [description = "Minimum Relative Inventory Without Inflation"]
    @parameters OO = params[:OO] [description = "Optimal Output in 1980 Gu/y"]
    @parameters PH = params[:PH] [description = "Pulse Height"]
    @parameters PPU = params[:PPU] [description = "Price Per Unit /u"]
    @parameters SAT = params[:SAT] [description = "Sales Averaging Time y"]
    @parameters INVEODDI = params[:INVEODDI] [description = "sINVeoDDI < 0: INVentory Effect On Delivery Delay Index"]
    @parameters INVEOIN = params[:INVEOIN] [description = "sINVeoIN < 0: INVentory Effect On INflation"]
    @parameters INVEOSWI = params[:INVEOSWI] [description = "sINVeoSWI < 0: INVentory Effect On Shifts Worked Index"]
    @parameters SRI = params[:SRI] [description = "Sufficient Relative Inventory"]
    @parameters SWI1980 = params[:SWI1980] [description = "SWI in 1980"]
    @parameters TAS = params[:TAS] [description = "Time to Adjust Shifts y"]

    @variables DELDI(t) = inits[:DELDI] [description = "DELivery Delay - Index (1)"]
    @variables EPP(t) = inits[:EPP] [description = "Effective Purchasing Power Gdollar/y"]
    @variables INV(t) = inits[:INV] [description = "INVentory Gu"]
    @variables PRIN(t) = inits[:PRIN] [description = "PRice INdex (1980=1)"]
    @variables PRI(t) = inits[:PRI] [description = "Perceived Relative Inventory (1)"]
    @variables RS(t) = inits[:RS] [description = "Recent Sales Gu/y"]
    @variables SSWI(t) = inits[:SSWI] [description = "ShiftS Worked - Index (1)"]
    @variables CDDI(t) [description = "Change in DDI 1/y"]
    @variables CPI(t) [description = "Change in Price Index 1/y"]
    @variables DEL(t) [description = "DELiveries Gu/y"]
    @variables DEPU(t) [description = "DEmand PUlse 2020-25 (1)"]
    @variables DSWI(t) [description = "Desired Shifts Worked - Index (1)"]
    @variables GDP(t) [description = "GDP Gdollar/y"]
    @variables IC(t) [description = "Inventory Coverage y"]
    @variables IR(t) [description = "Inflation Rate 1/y"]
    @variables NI(t) [description = "National Income Gdollar/y"]
    @variables OUTP(t) [description = "OUTPut Gu/y"]
    @variables PNIS(t) [description = "Pink Noise In Sales (1)"]
    @variables ROCDDI(t) [description = "ROC in DDI 1/y"]
    @variables SA(t) [description = "SAles Gdollar/y"]

    # labourmarket
    @parameters AUR = params[:AUR] [description = "Acceptable Unemployment Rate (1)"]
    @parameters FIC = params[:FIC] [description = "Fraction of Inflation Compensated (1)"]
    @parameters GDPP1980 = params[:GDPP1980] [description = "GDP per Person in 1980"]
    @parameters GDPPEROCCLRM = params[:GDPPEROCCLRM] [description = "sGDPppeoROCCLR<0: GDPP Effect on Rate Of Change in Capital Labour Ratio"]
    @parameters GENLPR = params[:GENLPR] [description = "Goal for Extra Normal LPR (1)"]
    @parameters NLPR80 = params[:NLPR80] [description = "Normal LPR in 1980 (1)"]
    @parameters PFTJ = params[:PFTJ] [description = "Persons per Full-Time Job p/ftj"]
    @parameters PFTJ80 = params[:PFTJ80] [description = "Persons per Full-Time Job in 1980 p/ftj"]
    @parameters PRUN = params[:PRUN] [description = "PRice per UNit dollar/u"]
    @parameters PUELPR = params[:PUELPR] [description = "sPUNeoLPR>0: Perceived Unemployment Effect on Labour Participation Rate"]
    @parameters ROCECLR80 = params[:ROCECLR80] [description = "ROC in ECLR in 1980 1/y"]
    @parameters RWER = params[:RWER] [description = "Real Wage Erosion Rate 1/y"]
    @parameters TAHW = params[:TAHW] [description = "Time to Adjust Hours Worked y"]
    @parameters TELLM = params[:TELLM] [description = "Time to Enter/Leave Labor Market y"]
    @parameters TENHW = params[:TENHW] [description = "sTIeoNHW<0: Time Effect on Number Hours Worked"]
    @parameters TYLD = params[:TYLD] [description = "10-Yr Loop Delay y"]
    @parameters WSOECLR = params[:WSOECLR] [description = "sWSOeoCLR>0: Worker Share of Output Effect on Capital Labour Ratio "]
    @parameters WSOELPR = params[:WSOELPR] [description = "sWSOeoLPR>0: Worker Share of Output Effect on Labour Participation Rate"]

    @variables AGIW(t) [description = "Average Gross Income per Worker kdollar/p/y"]
    @variables AHW(t) [description = "Average Hours Worked kh/y"]
    @variables AHW1980(t) [description = "Average Hours Worked in 1980 kh/y"]
    @variables AVWO(t) [description = "AVailable WOrkforce Mp"]
    @variables CECLR(t) [description = "Change in Embedded CLR kcu/ftj/y"]
    @variables CHWO(t) [description = "CHange in WOrkforce Mp/y"]
    @variables CWRA(t) [description = "Change in Wage RAte dollar/ph/y"]
    @variables CWSO(t) [description = "Change in WSO 1/y"]
    @variables ECLR(t) = inits[:ECLR] [description = "Embedded CLR kcu/ftj"]
    @variables ENLPR2022(t) [description = "Extra Normal LPR from 2022 (1)"]
    @variables GDPPEROCCLR(t) [description = "GDPppeoROCCLR"]
    @variables HFD(t) [description = "Hiring/Firing Delay y"]
    @variables HWMGDPP(t) [description = "Hours Worked Mult from GDPpP (1)"]
    @variables IWEOCLR(t) [description = "Indicated Wage Effect on Optimal CLR (1)"]
    @variables ILPR(t) [description = "Indicated Labour Participation Rate (1)"]
    @variables LAPR(t) [description = "LAbour PRoductivity dollar/ph"]
    @variables LAUS(t) = inits[:LAUS] [description = "LAbour USe Gph/y"]
    @variables LAUS80(t) [description = "LAbour USe in 1980 Gph/y"]
    @variables LPR(t) = inits[:ILPR] [description = "Labour Participation Rate (1)"]
    @variables LTEWSO(t) [description = "Long-Term Erosion of WSO 1/y"]
    @variables NHW(t) = inits[:NHW] [description = "Normal Hours Worked kh/ftj/y"]
    @variables NLPR(t) [description = "Normal LPR (1)"]
    @variables OCLR(t) [description = "Optimal Capital Labour Ratio kcu/ftj"]
    @variables OPWO(t) [description = "OPtimal WOrkforce Mp"]
    @variables PART(t) [description = "Participation (1)"]
    @variables PSW(t) [description = "Perceived Surplus Workforce (1)"]
    @variables PURA(t) = inits[:PURA] [description = "Perceived Unemployment RAte (1)"]
    @variables ROCECLR(t) [description = "ROC in ECLR 1/y"]
    @variables ROCWSO(t) [description = "ROC in WSO - Table 1/y"]
    @variables TCT(t) [description = "Time to Change Tooling y"]
    @variables UNEM(t) [description = "UNEMployed Mp"]
    @variables UR(t) [description = "UNemployment RAte (1)"]
    @variables UPT(t) [description = "Unemployment perception time y"]
    @variables WAP(t) [description = "Working Age Population Mp"]
    @variables WARA(t) = inits[:WARA] [description = "WAge RAte dollar/ph"]
    @variables WASH(t) [description = "WAge Share (1)"]
    @variables WEOCLR(t) = inits[:WEOCLR] [description = "Wage Effect on Optimal CLR (1)"]
    @variables WF(t) = inits[:WF] [description = "WorkForce Mp"]
    @variables WRER(t) [description = "Wage Rate Erosion Rate 1/y"]
    @variables WRE(t) [description = "Wage Rate Erosion dollar/ph/y"]
    @variables WSO(t) = inits[:WSO] [description = "Worker Share of Output (1)"]

    # other
    @parameters INELOK = params[:INELOK] [description = "sINEeoLOK<0: INequity Effect on LOgistic K"]
    @parameters NK = params[:NK] [description = "Normal K"]
    @parameters TEGR = params[:TEGR] [description = "Time to Establish Growth Rate y"]

    @variables CFETA(t) [description = "Cost of Food and Energy TAs GDollar/y"]
    @variables CTA(t) [description = "Cost of TAs GDollar/y"]
    @variables FB15(t) [description = "Fraction Below 15 kDollar/p/y (1)"]
    @variables IEL(t) [description = "Inequity Effect on Logistic k (1)"]
    @variables LK(t) [description = "Logistic K (1)"]
    @variables PB15(t) [description = "Population Below 15 kDollar/p/y Mp"]
    @variables PGDPP(t) = inits[:PGDPP] [description = "Past GDP per Person kDollar/y"]
    @variables RGGDPP(t) [description = "Rate of Growth in GDP per Person 1/y"]

    # output
    @parameters CAPPIS1980 = params[:CAPPIS1980] [description = "CAP PIS in 1980 Gcu"]
    @parameters CAPPUS1980 = params[:CAPPUS1980] [description = "CAP PUS in 1980 Gcu"]
    @parameters CC1980 = params[:CC1980] [description = "Cost of Capacity in 1980 dollar/cu"]
    @parameters CTPIS = params[:CTPIS] [description = "Construction Time PIS y"]
    @parameters CTPUS = params[:CTPUS] [description = "Construction Time PUS y"]
    @parameters EMCUC = params[:EMCUC] [description = "Extra Mult on CUC, to avoid initial transient in investment share of GDP"]
    @parameters FCI = params[:FCI] [description = "Foreign Capital Inflow Gdollar/y"]
    @parameters LCPIS1980 = params[:LCPIS1980] [description = "Life of Capacity PIS in 1980 y"]
    @parameters OW2022 = params[:OW2022] [description = "OBserved WArming in 2022 deg C"]
    @parameters OWECCM = params[:OWECCM] [description = "sOWeoCOC>0: Observed Warming Effect on Cost of Capacity Multiplier"]
    @parameters OWELCM = params[:OWELCM] [description = "sOWeoLOC<0: Observed Warming Effect on Life of Capacity Multiplier"]
    @parameters USPIS2022 = params[:USPIS2022] [description = "Unconventional Stimulus in PIS from 2022 (share of GDP)"]
    @parameters USPUS2022 = params[:USPUS2022] [description = "Unconventional Stimulus in PUS from 2022 (share of GDP)"]
    @parameters CBCEFRA = params[:CBCEFRA] [description = "sCBCeoFRA<0: Corporate Borrowing Cost Effect on FRA"]
    @parameters CUCPIS1980 = (params[:CAPPIS1980] / params[:LCPIS1980]) * params[:CTPIS] * params[:EMCUC] [description = "CUC PIS in 1980 Gcu"]
    @parameters ED1980 = params[:ED1980] [description = "Excess Demand in 1980"]
    @parameters EDEFRA = params[:EDEFRA] [description = "sEDeoFRA>0: Excess Demand Effect on FRA"]
    @parameters EDELCM = params[:EDELCM] [description = "sEDeoLOC>0: Excess Demand Effect on Life of Capacity"]
    @parameters FRA1980 = params[:FRA1980] [description = "FRA in 1980"]
    @parameters FRACAM = params[:FRACAM] [description = "FRACA Min"]
    @parameters GDPP1980 = params[:GDPP1980] [description = "GDP per Person in 1980"]
    @parameters GDPPEFRACA = params[:GDPPEFRACA] [description = "sGDPppeoFRACA<0: GDP per Person Effect on FRACA"]
    @parameters IPT = params[:IPT] [description = "Investment Planning Time y"]
    @parameters JOBS1980 = params[:JOBS1980] [description = "JOBS in 1980 M-ftj"]
    @parameters KAPPA = params[:KAPPA] [description = "Kappa"]
    @parameters LAMBDA = params[:LAMBDA] [description = "Lambda"]
    @parameters LAUS1980 = params[:LAUS1980] [description = "Labour Use in 1980 Gph/y"]
    @parameters OG1980 = params[:OG1980] [description = "Output Growth in 1980 1/y (to avoid transient)"]
    @parameters OO1980 = params[:OO1980] [description = "Optimal output in 1980 Gu/y"]
    @parameters PRUN = params[:CU1980] * (1 + params[:MA1980]) [description = "Price per Unit dollar/u"]
    @parameters TOED = params[:TOED] [description = "Time to observe excess demand y"]
    @parameters WSOEFRA = params[:WSOEFRA] [description = "sWSOeoFRA<0: Worker Share of Output Effect on FRA"]

    @variables AVCA(t) [description = "AVailable CApital Gdollar/y"]
    @variables CAPIS(t) [description = "Capacity Addition PIS Gcu/y"]
    @variables CAPUS(t) [description = "Capacity Addition PUS Gcu/y"]
    @variables CDPIS(t) [description = "Capacity Discard PIS Gcu/y"]
    @variables CDPUS(t) [description = "Capacity Discard PUS Gcu/y"]
    @variables CIPIS(t) [description = "Capacity Initiation PIS Gcu/y"]
    @variables CIPUS(t) [description = "Capacity Initiation PUS Gcu/y"]
    @variables COCA(t) [description = "COst of CApacity dollar/cu"]
    @variables CPIS(t) = params[:CAPPIS1980] [description = "Capacity PIS Gcu"]
    @variables CPUS(t) = params[:CAPPUS1980] [description = "Capacity PUS Gcu"]
    @variables CRR(t) [description = "Capacity Renewal Rate 1/y"]
    @variables CUCPIS(t) = inits[:CUCPIS] [description = "Capacity Under Construction PIS Gcu"]
    @variables CUCPIS1980(t) [description = "CUC PIS in 1980 Gcu"]
    @variables CUCPUS(t) = inits[:CUCPUS] [description = "Capacity Under Construction PUS Gcu"]
    @variables CUCPUS1980(t) [description = "CUC PUS in 1980 Gcu"]
    @variables ECR(t) [description = "Effect of Capacity Renewal 1/y"]
    @variables ETFP(t) = inits[:ETFP] [description = "Embedded TFP (1)"]
    @variables INCPIS(t) [description = "Investment in New Capacity PIS Gdollar/y"]
    @variables LCPUS(t) [description = "Life of Capacity PUS y"]
    @variables LCPUS1980(t) [description = "Life of Capacity PUS in 1980 y"]
    @variables OBSGIPIS(t) [description = "Off-Balance Sheet Govmnt Inv in PIS (share of GDP)"]
    @variables OBSGIPUS(t) [description = "Off-Balance-Sheet Govmnt Inv in PUS (share of GDP)"]
    @variables OWECC(t) [description = "OWeoCOC (1)"]
    @variables OWELC(t) [description = "OWeoLOC (1)"]
    @variables CBCEFCA(t) [description = "CBC Effect on Flow to Capacity Addion (1)"]
    @variables CAP(t) [description = "CAPacity Gcu"]
    @variables EDE(t) [description = "Excess DEmand (1)"]
    @variables EDEFCA(t) [description = "ED Effect on Flow to Capacity Addition (1)"]
    @variables EDELC(t) [description = "Excess Demand Effect on Life of Capacity (1)"]
    @variables FACNC(t) = inits[:FACNC] [description = "Fraction of Available Capital to New Capacity (1)"]
    @variables FRACAMGDPPL(t) [description = "FRACA Mult from GDPpP - Line (1)"]
    @variables FRACAMGDPPT(t) [description = "FRACA Mult from GDPpP - Table (1)"]
    @variables ISGDP(t) [description = "Investment Share of GDP (1)"]
    @variables LCPIS(t) [description = "Life of Capacity PIS y"]
    @variables OGR(t) = params[:OG1980] [description = "Output Growth Rate 1/y"]
    @variables OLY(t) = inits[:OLY] [description = "Output Last Year Gdollar/y"]
    @variables OOV(t) [description = "Optimal Ouput - Value Gdollar/y"]
    @variables ORO(t) [description = "Optimal Real Output Gu/y"]
    @variables PEDE(t) = 1 [description = "Perceived Excess DEmand (1)"]
    @variables WSOEFCA(t) [description = "WSO Effect on Flow to Capacity Addition (1)"]

    # population
    @parameters CMFR = params[:CMFR] [description = "Cost of Max Fertility Reduction (share of GDP)"]
    @parameters DNC80 = params[:DNC80] [description = "DNC in 1980"]
    @parameters DNCA = params[:DNCA] [description = "DNCalfa<0"]
    @parameters DNCG = params[:DNCG] [description = "DNCgamma"]
    @parameters DNCM = params[:DNCM] [description = "DNCmin"]
    @parameters EIP = params[:EIP] [description = "Exogenous Introduction Period y"]
    @parameters EIPF = params[:EIPF] [description = "Exogenous Introduction Period Flag"]
    @parameters FADFS = params[:FADFS] [description = "Fraction Achieving Desired Family Size"]
    @parameters FP = params[:FP] [description = "Fertile Period"]
    @parameters FW = params[:FW] [description = "Fraction Women"]
    @parameters GEFR = params[:GEFR] [description = "Goal for Extra Fertility Reduction"]
    @parameters GEPA = params[:GEPA] [description = "Goal for Extra Pension Age y"]
    @parameters LEA = params[:LEA] [description = "LEalfa"]
    @parameters LEEPA = params[:LEEPA] [description = "sLEeoPa>0: Life Expectancy Effect on Pension Age"]
    @parameters LEG = params[:LEG] [description = "LEgamma"]
    @parameters LEMAX = params[:LEMAX] [description = "LEmax"]
    @parameters MFM = params[:MFM] [description = "Max Fertility Multiplier"]
    @parameters MLEM = params[:MLEM] [description = "Max Life Expectancy Multiplier"]
    @parameters OW2022 = params[:OW2022] [description = "Observed Warming in 2022 deg C"]
    @parameters OWELE = params[:OWELE] [description = "sOWeoLE<0: Observed Warming Effect on Life Expectancy"]
    @parameters SSP2FA2022F = params[:SSP2FA2022F] [description = "SSP2 Family Action from 2022 Flag"]
    @parameters TAHI = params[:TAHI] [description = "Time to adapt to higher income y"]

    @variables A0020(t) = inits[:A0020] [description = "Aged 0-20 years Mp"]
    @variables A2040(t) = inits[:A2040] [description = "Aged 20-40 years Mp"]
    @variables A20PA(t) [description = "Aged 20-Pension Age Mp"]
    @variables A4060(t) = inits[:A4060] [description = "Aged 40-60 Mp"]
    @variables A60PL(t) = inits[:A60PL] [description = "Aged 60 + Mp"]
    @variables BIRTHR(t) [description = "Birth Rate 1/y"]
    @variables BIRTHS(t) [description = "Births Mp/y"]
    @variables CEFR(t) [description = "Cost of Extra Fertility Reduction (share of GDP)"]
    @variables DEATHR(t) [description = "Death Rate 1/y"]
    @variables DEATHS(t) [description = "Deaths Mp/y"]
    @variables DNC(t) [description = "Desired No of Children 1"]
    @variables DR(t) [description = "Dependency Ratio p/p"]
    @variables EFR(t) [description = "Extra Fertility Reduction (1)"]
    @variables EGDPP(t) = inits[:EGDPP] [description = "Effective GDP per Person kDollar/p/y"]
    @variables EPA(t) = inits[:EPA] [description = "Extra Pension Age y"]
    @variables FM(t) [description = "Fertility Multiplier (1)"]
    @variables GDPP(t) [description = "GDP per Person kDollar/p/y"]
    @variables LE(t) = inits[:LE] [description = "Life Expectancy y"]
    @variables LE60(t) [description = "LE at 60 y"]
    @variables LEM(t) [description = "Life Expectancy Multipler (1)"]
    @variables (LV_DEATHS(t))[1:ORDER] = fill(inits[:DEATHS] * (inits[:LE] - 60) / ORDER, ORDER) [description = "LV functions for deaths Mp/y"]
    @variables (LV_PASS20(t))[1:ORDER] = fill(inits[:PASS20] * 20 / ORDER, ORDER) [description = "LV functions for passing 20 Mp/y"]
    @variables (LV_PASS40(t))[1:ORDER] = fill(inits[:PASS40] * 20 / ORDER, ORDER) [description = "LV functions for passing 40 Mp/y"]
    @variables (LV_PASS60(t))[1:ORDER] = fill(inits[:PASS60] * 20 / ORDER, ORDER) [description = "LV functions for passing 60 Mp/y"]
    @variables OF(t) [description = "Observed Fertility 1"]
    @variables OP(t) [description = "On Pension Mp"]
    @variables PA(t) = inits[:PA] [description = "Pension Age y"]
    @variables PASS20(t) [description = "Passing 20 Mp/y"]
    @variables PASS40(t) = inits[:PASS40] [description = "Passing 40 Mp/y"]
    @variables PASS60(t) = inits[:PASS60] [description = "Passing 60 Mp/y"]
    @variables PGR(t) [description = "Population Growth Rate 1/y"]
    @variables POP(t) [description = "Population Mp"]
    @variables PW(t) [description = "Pensioners per Worker p/p"]
    @variables (RT_DEATHS(t))[1:ORDER] = fill(inits[:DEATHS] * (inits[:LE] - 60) / ORDER, ORDER) [description = "RT functions for deaths Mp/y"]
    @variables (RT_PASS20(t))[1:ORDER] = fill(inits[:PASS20] * 20 / ORDER, ORDER) [description = "RT functions for passing 20 Mp/y"]
    @variables (RT_PASS40(t))[1:ORDER] = fill(inits[:PASS40] * 20 / ORDER, ORDER) [description = "RT functions for passing 40 Mp/y"]
    @variables (RT_PASS60(t))[1:ORDER] = fill(inits[:PASS60] * 20 / ORDER, ORDER) [description = "RT functions for passing 60 Mp/y"]
    @variables WELE(t) [description = "Warming Effect on Life Expectancy (1)"]

    # public
    @parameters OW2022 = params[:OW2022] [description = "Observed Warming in 2022 deg C"]
    @parameters CTA2022 = params[:CTA2022] [description = "Cost of TAs in 2022 GDollar/y"]
    @parameters CTPIS = params[:CTPIS] [description = "Construction time PIS y"]
    @parameters DROTA1980 = params[:DROTA1980] [description = "Domestic ROTA in 1980 1/y"]
    @parameters EDROTA2022 = params[:EDROTA2022] [description = "Extra Domestic ROTA in 2022 1/y"]
    @parameters FUATA = params[:FUATA] [description = "Fraction Unprofitable Activity in TAs"]
    @parameters GDPTL = params[:GDPTL] [description = "GDPpp of Technology Leader kDollar/p/k"]
    @parameters IIEEROTA = params[:IIEEROTA] [description = "sIIEeoROTA<0: Inequality Index Effect on ROTA"]
    @parameters IPR1980 = params[:IPR1980] [description = "Infrastructure Purchases Ratio in 1980 y"]
    @parameters IPRVPSS = params[:IPRVPSS] [description = "sIPReoVPSS>0: Infrastructure Purchase Ratio effect on Value of Pubblic Services Supplied"]
    @parameters IPT = params[:IPT] [description = "Investment planning time y"]
    @parameters MIROTA2022 = params[:MIROTA2022] [description = "Max Imported ROTA from 2022 1/y"]
    @parameters OWETFP = params[:OWETFP] [description = "sOWeoTFP<0: Observed Warming Effect on Total Factor Productivity"]
    @parameters SC1980 = params[:SC1980] [description = "State Capacity in 1980 (fraction of GDP)"]
    @parameters SCROTA = params[:SCROTA] [description = "sSCeoROTA>0: State Capacity effect on Rate Of Technological Advance"]
    @parameters XETAC2022 = params[:XETAC2022] [description = "XExtra TA Cost in 2022 (share of GDP)"]
    @parameters XETAC2100 = params[:XETAC2100] [description = "XExtra TA Cost in 2100 (share of GDP)"]

    @variables CTFP(t) [description = "Change in TFP 1/y"]
    @variables DRTA(t) [description = "Domestic Rate of Technological Advance 1/y"]
    @variables ECTAF2022(t) [description = "Extra Cost of TAs From 2022 GDollar/y"]
    @variables ECTAGDP(t) [description = "Extra Cost of TAs as share of GDP (1)"]
    @variables GSSGDP(t) [description = "Govmnt Spending as Share of GDP"]
    @variables IPR(t) [description = "Infrastructure Purchases Ratio y"]
    @variables IROTA(t) [description = "Imported ROTA 1/y"]
    @variables ITFP(t) [description = "Indicated TFP (1)"]
    @variables OWTFP(t) [description = "OWeoTFP"]
    @variables PLUA(t) [description = "Productivity Loss from Unprofitable Activity (1)"]
    @variables PPP(t) [description = "Productivity of Public Purchases (1)"]
    @variables PSEP(t) [description = "Public SErvices per Person kdollar/p/y"]
    @variables PSP(t) [description = "Public Spending per Person kDollar/p/y"]
    @variables RROTAI(t) [description = "Reduction in ROTA from Inequality 1/y"]
    @variables RTA(t) [description = "Rate of Technological Advance 1/y"]
    @variables RTFPUA(t) = inits[:RTFPUA] [description = "Reduction in TFP from Unprofitable Activity (1)"]
    @variables SC(t) [description = "State Capacity (fraction of GDP)"]
    @variables TFPEE5TA(t) = inits[:TFPEE5TA] [description = "TFP Excluding Effect of 5TAs (1)"]
    @variables TFPIE5TA(t) [description = "TFP Including Effect of 5TAs (1)"]
    @variables VPSS(t) [description = "Value of Public Services Supplied GDollar/y"]
    @variables XECTAGDP(t) [description = "XExtra Cost of TAs as share of GDP (1)"]

    # wellbeing
    @parameters AI = params[:AI] [description = "Acceptable Inequality"]
    @parameters AP = params[:AP] [description = "Acceptable Progress 1/y"]
    @parameters AWBPD = params[:AWBPD] [description = "Average WellBeing Perception Delay y"]
    @parameters DRDI = params[:DRDI] [description = "Diminishing Return Disposable Income"]
    @parameters DRPS = params[:DRPS] [description = "Diminishing Return Public Spending"]
    @parameters EIP = params[:EIP] [description = "Exogenous Introduction Period y"]
    @parameters EIPF = params[:EIPF] [description = "Exogenous Introduction Period Flag"]
    @parameters GWEAWBGWF = params[:GWEAWBGWF] [description = "sGWeoAW<0: Global Warming Effect on Average WellBeing from Global Warming Flag"]
    @parameters IEAWBIF = params[:IEAWBIF] [description = "sIIeoAW<0: Inequality Effect on Average WellBeing from Inequality Flag"]
    @parameters MWBGW = params[:MWBGW] [description = "Min WellBeing from Global Warming"]
    @parameters NRD = params[:NRD] [description = "Normal Reform Delay y"]
    @parameters PAEAWBF = params[:PAEAWBF] [description = "sLPeoAWP>0: PArtecipation Effect on Average WellBeing Flag"]
    @parameters PREAWBF = params[:PREAWBF] [description = "sROPeoAW>0: PRogress Effect on Average WellBeing Flag"]
    @parameters PESTF = params[:PESTF] [description = "sPPReoSTE<0: Progress Effect on Social Tension Flag"]
    @parameters SPS = params[:SPS] [description = "Satisfactory Public Spending"]
    @parameters STEERDF = params[:STEERDF] [description = "sSTEeoRD>0: Social Tension Effect on Reform Delay Flag"]
    @parameters STRERDF = params[:STRERDF] [description = "sSTReoRD<0: Social Trust Effect on Reform Delay Flag"]
    @parameters TCRD = params[:TCRD] [description = "Time to Change Reform Delay y"]
    @parameters TDI = params[:TDI] [description = "Threshold Disposable Income kdollar/p/y"]
    @parameters TEST = params[:TEST] [description = "Time to Establish Social Trust y"]
    @parameters TI = params[:TI] [description = "Threshold Inequality"]
    @parameters THPA = params[:THPA] [description = "Threshold Participation"]
    @parameters TPR = params[:TPR] [description = "Threshold Progress Rate 1/y"]
    @parameters TPS = params[:TPS] [description = "Threshold Public Spending kdollar/p/y"]
    @parameters TW = params[:TW] [description = "Threshold Warming deg C"]

    @variables AWBDI(t) [description = "Average WellBeing from Disposable Income (1)"]
    @variables AWBGW(t) [description = "Average WellBeing from Global Warming (1)"]
    @variables AWBI(t) [description = "Average WellBeing Index (1)"]
    @variables AWBIN(t) [description = "Average WellBeing from INequality (1)"]
    @variables AWBPS(t) [description = "Average WellBeing from Public Spending (1)"]
    @variables AWBP(t) [description = "Average WellBeing from Progress (1)"]
    @variables IEST(t) [description = "Inequity Effect on Social Trust (1)"]
    @variables IPP(t) [description = "Introduction Period for Policy y"]
    @variables IRD(t) [description = "Indicated Reform Delay y"]
    @variables IST(t) [description = "Indicated Social Trust (1)"]
    @variables ORP(t) = inits[:ORP] [description = "Observed Rate of Progress 1/y"]
    @variables PAWBI(t) = inits[:PAWBI] [description = "Past AWI (1)"]
    @variables PSESTR(t) [description = "Public Spending Effect on Social TRust (1)"]
    @variables PSSGDP(t) [description = "Public Spending as Share of GDP"]
    @variables RD(t) = inits[:RD] [description = "Reform Delay y"]
    @variables SOTR(t) = inits[:SOTR] [description = "Social Trust (1)"]
    @variables STE(t) = inits[:STE] [description = "Social TEnsion (1)"]
    @variables STEERD(t) [description = "Social TEnsion Effect on Reform Delay (1)"]
    @variables STRERD(t) [description = "Social TRust Effect on Reform Delay (1)"]
    @variables WBEP(t) [description = "WellBeing Effect of Participation (1)"]

    #other parameters 
    @parameters INITFOLA = params[:INITFOLA]
    @parameters INITNHW = params[:INITNHW]
    @parameters INITWSO = params[:INITWSO]
    @parameters INITWF = params[:INITWF]
    @parameters INITEGDPP = params[:INITEGDPP]
    @parameters INITEPA = params[:INITEPA]
    @parameters INITLE = params[:INITLE]
    @parameters INITPA = params[:INITPA]
    @parameters INITSTE = params[:INITSTE]
    @parameters INITSOTR = params[:INITSOTR]


    # climate
    add_equation!(eqs, KN2OEKF ~ KN2OKF1980 * exp(-(RDN2OKF) * (t - 1980)) * IfElse.ifelse(t > 2022, exp(-(ERDN2OKF2022) * (t - 2022)), 1))
    add_equation!(eqs, MMN2OE ~ FEUS * KN2OEKF / 1000)
    add_equation!(eqs, NN2OE ~ withlookup(t, [(1980.0, 0.009), (2020.0, 0.009), (2099.27, 0.0)]))
    add_equation!(eqs, N2OE ~ NN2OE + MMN2OE)
    add_equation!(eqs, N2OC1980 ~ N2OA1980 / MAT)
    add_equation!(eqs, D(N2OA) ~ N2OE - N2OBD)
    add_equation!(eqs, N2OBD ~ N2OA / LN2OA)
    add_equation!(eqs, N2OCA ~ N2OA / GN2OPP)
    add_equation!(eqs, N2OFPP ~ withlookup(t, [(1980.0, 0.43), (2000.0, 0.64), (2010.0, 0.73), (2020.0, 0.8), (2100.0, 1.0)]))
    add_equation!(eqs, FN2O ~ N2OCA * N2OFPP)
    add_equation!(eqs, KCH4EKC ~ KCH4KC1980 * exp(-(RDCH4KC) * (t - 1980)) * IfElse.ifelse(t > 2022, exp(-(ERDCH4KC2022) * (t - 2022)), 1))
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
    add_equation!(eqs, D(CO2A) ~ CO2E - CO2AB + 2 * CO2FCH4) #STRANGE EQUATION!
    add_equation!(eqs, CO2CA ~ CO2A / GCO2PP)
    add_equation!(eqs, CO2FPP ~ withlookup(t, [(1980.0, 0.0032), (1990.0, 0.0041), (2000.0, 0.0046), (2020.0, 0.0051), (2100.0, 0.006)]))
    add_equation!(eqs, FCO2 ~ CO2CA * CO2FPP)
    add_equation!(eqs, FOG ~ withlookup(t, [(1980.0, 0.18), (2000.0, 0.36), (2020.0, 0.39), (2050.0, 0.37), (2100.0, 0.0)]))
    add_equation!(eqs, MMF ~ FCO2 + FOG + FCH4 + FN2O)
    add_equation!(eqs, GHGE ~ CO2E * TCO2ETCO2 + CH4E * TCO2ETCH4 + N2OE * TCO2ETN2O)
    add_equation!(eqs, AL1980 ~ (ISCEGA1980 * ALIS + (GLSU - ISCEGA1980) * ALGAV) / GLSU)
    add_equation!(eqs, AL ~ (ISCEGA * ALIS + (GLSU - ISCEGA) * ALGAV) / GLSU)
    add_equation!(eqs, TRHGS ~ (TRSS1980 * ((OW + 297) / 297)) * (AL / AL1980))
    add_equation!(eqs, HTS ~ EHS * TRHGS) # EHS
    add_equation!(eqs, MRS ~ MRS1980 * (OW / WA1980))
    add_equation!(eqs, MRDI ~ MRS / SVDR)
    add_equation!(eqs, ECIM ~ MRDI * AI1980 * TPM3I * HRMI)
    add_equation!(eqs, MEL ~ ISCEGA * MRS)
    add_equation!(eqs, D(ISCEGA) ~ -MEL)
    add_equation!(eqs, ISC ~ ISCEGA * 100)
    add_equation!(eqs, WVC ~ WVC1980 * (1 + OWWV * (OW / WA1980 - 1)))
    add_equation!(eqs, WVF ~ WVF1980 * (1 + WVWVF * (WVC / WVC1980 - 1)))
    add_equation!(eqs, TMMF ~ MMF + WVF)
    add_equation!(eqs, EWFF ~ (TMMF * GLSU) * 31.5 / 1000)
    add_equation!(eqs, OW ~ WA1980 + (EHS - EH1980) * WFEH)
    add_equation!(eqs, REHE ~ withlookup(OW, [(0.0, 1.0), (1.2, 4.8), (2.0, 8.6), (2.9, 14.0), (5.2, 40.0)]))
    add_equation!(eqs, TRHGA ~ TRSA1980 * ((OW + 287) / 287))
    add_equation!(eqs, HDO ~ EHS * TRHGA)
    add_equation!(eqs, D(EHS) ~ EWFF - ECIM - HDO - HTS)
    smooth!(eqs, PWA, OW, PD)

    # demand
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

    # energy
    add_equation!(eqs, CNEL ~ NEP * CNED)
    add_equation!(eqs, NFCO2PP ~ MNFCO2PP * (1 - exp(-(GDPP / 10))))
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
    add_equation!(eqs, TPPUEBEE ~ withlookup(GDPP, [(0.0, 0.0), (10.0, 4.0), (20.0, 7.0), (30.0, 9.0), (50.0, 12.0), (65.0, 13.0)]))
    add_equation!(eqs, TPPUFFNEUBEE ~ withlookup(GDPP, [(0.0, 0.3), (15.0, 2.0), (25.0, 3.1), (35.0, 4.0), (50.0, 5.0)]))
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
    add_equation!(eqs, NDSWC ~ log(2) + log(ACSWCF1980 / SWC1980))
    add_equation!(eqs, CISWC ~ (1 - CRDSWC)^NDSWC)
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

    # finance
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

    # foodland
    add_equation!(eqs, ACY ~ (DCYCA * (1 - FRA) * SQICA + CYRA * FRA) * CO2ELY * WELY)
    add_equation!(eqs, ALFL ~ 1 - exp(-FFLR / TFFLR))
    add_equation!(eqs, AFSRA ~ 268 - SFU)
    add_equation!(eqs, D(BALA) ~ CRLO)
    add_equation!(eqs, BIUS ~ withlookup(t, [(1980.0, 0.0), (1990.0, 0.0), (2000.0, 0.0), (2020.0, 0.0), (2100.0, 0.0)]))
    add_equation!(eqs, CEM ~ IfElse.ifelse(t > 2022, 1 - SSP2LMA * ramp(t, (1 - 0) / 78, 2022, 2100), 1))
    add_equation!(eqs, CIRA ~ (1 - CRDRA)^NDRA)
    add_equation!(eqs, CO2AFL ~ FOLA * (CO2AFLH / 1000) * CO2ELY * WELY)
    add_equation!(eqs, CO2AFLH ~ 1.6 * FAM)
    add_equation!(eqs, CO2ELULUC ~ CO2RFC - CO2AFL - ECO2ARA)
    add_equation!(eqs, CO2ELY ~ IfElse.ifelse(t > 2022, 1 + CO2CEACY * (CO2CA / CO2C2022 - 1), 1))
    add_equation!(eqs, CO2RFC ~ ((OGRE + CREX) * CO2RHFC) / 1000)
    add_equation!(eqs, COFE ~ FEUS * CTF / 1000)
    add_equation!(eqs, COFO ~ AFGDP * GDP + CRA + COFE)
    add_equation!(eqs, CRA ~ (ECRA * RAA) / 1000)
    add_equation!(eqs, CRBA ~ CRUS / DCS)
    add_equation!(eqs, CRBI ~ BIUS * TCTB)
    add_equation!(eqs, CRDE ~ (TUCERM + FERM + CRBI))
    add_equation!(eqs, CREX ~ IfElse.ifelse(FOLA > 0, CRLA * CREXR, 0) * ALFL * CEM)
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
    add_equation!(eqs, DCYCA ~ DCSCA / (CRLA * (1 - FRA)))
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
    add_equation!(eqs, FFLR ~ max(0, FOLA / INITFOLA))
    add_equation!(eqs, FFLREOGRR ~ max(1, 1 + FFLREOGRRM * (FFLR - TFFLR)))
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
    add_equation!(eqs, NDRA ~ log((RAA + EGB22) / EGB22) / 0.693)
    add_equation!(eqs, NFL ~ OGRE * (1 - FCG))
    add_equation!(eqs, NGL ~ OGRE * FCG)
    add_equation!(eqs, D(OGFA) ~ -NFL - NGL)
    add_equation!(eqs, OGRE ~ OGFA * OGRR * OGRRM)
    add_equation!(eqs, OGRR ~ OGRR80 * FFLREOGRR)
    add_equation!(eqs, OGRRM ~ IfElse.ifelse(t > 2022, 1 - SSP2LMA * ramp(t, (1 - 0) / 78, 2022, 2100), 1))
    add_equation!(eqs, PCB ~ CRBA / (1 + DRC))
    add_equation!(eqs, PRMGL ~ GRLA * GLY / 1000)
    add_equation!(eqs, RAA ~ CRLA * FRA)
    add_equation!(eqs, RMF ~ DRM - RMGL)
    add_equation!(eqs, RMGL ~ min(DRM, PRMGL))
    add_equation!(eqs, RMSP ~ (RMGL + RMF * min(1, CRBA)) * 1000 / POP)
    add_equation!(eqs, ROCSQCA ~ 0 + FUESQ * (FUCA / SFU - 1))
    add_equation!(eqs, D(SQICA) ~ CSQCA)
    add_equation!(eqs, TFA ~ OGFA + FOLA)
    add_equation!(eqs, TFUCA ~ withlookup(DCYCA, [(1.0, 0.0), (2.0, 40.0), (2.5, 50.0), (3.0, 60.0), (3.5, 70.0), (4.5, 100.0), (6.5, 200.0), (10.0, 600.0)]))
    add_equation!(eqs, TUC ~ TUCP * POP / 1000)
    add_equation!(eqs, TUCERM ~ (TUC - TUFRM) / FSPI)
    add_equation!(eqs, TUCERMP ~ TUCERM * 1000 / POP)
    add_equation!(eqs, TUCP ~ withlookup(GDPP, [(0.0, 400.0), (6.1, 680.0), (8.7, 780.0), (13.9, 950.0), (20.0, 1050.0), (30.0, 1150.0), (40.0, 1250.0), (60.0, 1350.0), (100.0, 1550.0)]))
    add_equation!(eqs, TUFRM ~ (((TURMP / 1000) * POP) - RMGL) * KCKRM)
    add_equation!(eqs, TURMP ~ withlookup(GDPP, [(0.0, 0.0), (6.1, 6.0), (8.8, 8.5), (14.0, 13.0), (30.0, 27.0), (40.0, 32.0), (50.0, 33.0), (100.0, 25.0)]))
    add_equation!(eqs, UREX ~ max(0, (IUL - URLA) / UDT))
    add_equation!(eqs, D(URLA) ~ UREX)
    add_equation!(eqs, WELY ~ IfElse.ifelse(t > 2022, 1 + OWEACY * (OW / OW2022 - 1), 1))

    # inventory
    add_equation!(eqs, CDDI ~ ROCDDI * DELDI)
    add_equation!(eqs, CPI ~ PRIN * IR)
    add_equation!(eqs, DEL ~ ((EPP / PPU) / (DELDI / DDI1980)) * IfElse.ifelse(t > 1984, PNIS, 1))
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

    # labourmarket
    add_equation!(eqs, AGIW ~ WARA * AHW)
    add_equation!(eqs, AHW ~ NHW / PFTJ)
    add_equation!(eqs, AHW1980 ~ INITNHW / PFTJ80)
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
    add_equation!(eqs, IWEOCLR ~ 1 + WSOECLR * (WSO / INITWSO - 1))
    add_equation!(eqs, ILPR ~ NLPR - PSW)
    add_equation!(eqs, LAPR ~ (OUTP * PRUN) / LAUS)
    add_equation!(eqs, LAUS ~ WF * AHW)
    add_equation!(eqs, LAUS80 ~ INITWF * AHW1980)
    smooth!(eqs, LPR, ILPR, TELLM)
    add_equation!(eqs, LTEWSO ~ WSO * RWER)
    smooth!(eqs, NHW, INITNHW * HWMGDPP, TAHW)
    add_equation!(eqs, NLPR ~ NLPR80 * (1 + WSOELPR * (WSO / INITWSO - 1)) + ENLPR2022)
    add_equation!(eqs, OCLR ~ ECLR * WEOCLR)
    add_equation!(eqs, OPWO ~ (CAP / OCLR) * PFTJ)
    add_equation!(eqs, PART ~ LPR * (1 - PURA))
    add_equation!(eqs, PSW ~ AUR * (1 + PUELPR * (PURA / AUR - 1)))
    smooth!(eqs, PURA, UR, UPT)
    add_equation!(eqs, ROCECLR ~ ROCECLR80 * GDPPEROCCLR)
    add_equation!(eqs, ROCWSO ~ withlookup(PURA / AUR, [(0.0, 0.06), (0.5, 0.02), (1.0, 0.0), (1.5, -0.007), (2.0, -0.01)]))
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

    # other
    add_equation!(eqs, CFETA ~ COFO + CE)
    add_equation!(eqs, CTA ~ CFETA)
    add_equation!(eqs, FB15 ~ 1 - (1 / (1 + exp(-LK * (GDPP - 14)))))
    add_equation!(eqs, IEL ~ 1 + INELOK * (INEQ / 0.5 - 1))
    add_equation!(eqs, LK ~ NK * IEL)
    add_equation!(eqs, PB15 ~ POP * FB15)
    add_equation!(eqs, RGGDPP ~ ((GDPP - PGDPP) / PGDPP) / TEGR)
    smooth!(eqs, PGDPP, GDPP, TEGR)

    # output
    add_equation!(eqs, AVCA ~ TS + FCI)
    add_equation!(eqs, CAPIS ~ CUCPIS / CTPIS)
    add_equation!(eqs, CAPUS ~ CUCPUS / CTPUS)
    add_equation!(eqs, CIPIS ~ max((INCPIS + OBSGIPIS * GDP) / COCA, 0))
    add_equation!(eqs, CIPUS ~ max((GIPC + OBSGIPUS * GDP) / COCA, 0))
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
    add_equation!(eqs, FRACAMGDPPL ~ max(FRACAM, 1 + GDPPEFRACA * (GDPP / GDPP1980 - 1)))
    add_equation!(eqs, FRACAMGDPPT ~ withlookup(GDPP / GDPP1980, [(0.0, 1.0), (1.0, 1.0), (2.0, 0.85), (2.1, 0.84), (4.0, 0.65), (8.0, 0.55), (16.0, 0.5)]))
    add_equation!(eqs, ISGDP ~ (INCPIS + GIPC) / GDP)
    add_equation!(eqs, LCPIS ~ (LCPIS1980 * OWELC) / EDELC)
    smooth!(eqs, OGR, (ORO - OLY) / OLY, 1)
    smooth!(eqs, OLY, ORO, 1)
    add_equation!(eqs, OOV ~ ORO * PRUN)
    add_equation!(eqs, ORO ~ OO1980 * ((CPIS + CPUS) / (CAPPIS1980 + CAPPUS1980))^KAPPA * (LAUS / LAUS1980)^LAMBDA * (ETFP))
    smooth!(eqs, PEDE, EDE, TOED)
    add_equation!(eqs, WSOEFCA ~ 1 + WSOEFRA * (WASH / INITWSO - 1))

    # population
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
    add_equation!(eqs, DEATHS ~ RT_DEATHS[ORDER])
    add_equation!(eqs, DNC ~ ((DNCM + (DNC80 - DNCM) * exp(-DNCG * (EGDPP - INITEGDPP))) * (1 + DNCA * (EGDPP - INITEGDPP))) * (1 - EFR) * FM)
    add_equation!(eqs, DR ~ (A0020 + A60PL) / (A2040 + A4060))
    add_equation!(eqs, EFR ~ ramp(t, GEFR / IPP, 2022, 2022 + IPP))
    add_equation!(eqs, EPA ~ ramp(t, (GEPA - INITEPA) / IPP, 2022, 22022 + IPP))
    add_equation!(eqs, D(EGDPP) ~ (GDPP - EGDPP) / TAHI)
    add_equation!(eqs, FM ~ IfElse.ifelse(SSP2FA2022F > 0, IfElse.ifelse(t > 2022, 1 + ramp(t, (MFM - 1) / 78, 2022, 2100), 1), 1))
    add_equation!(eqs, GDPP ~ GDP / POP)
    add_equation!(eqs, LE ~ ((LEMAX - (LEMAX - INITLE) * exp(-LEG * (EGDPP - INITEGDPP))) * (1 + LEA * (EGDPP - INITEGDPP))) * WELE * LEM)
    add_equation!(eqs, LE60 ~ LE - 60)
    add_equation!(eqs, LEM ~ IfElse.ifelse(SSP2FA2022F > 0, IfElse.ifelse(t > 2022, 1 + ramp(t, (MLEM - 1) / 78, 2022, 2100), 1), 1))
    add_equation!(eqs, OF ~ DNC * FADFS)
    add_equation!(eqs, OP ~ A60PL * (LE - PA) / (LE - 60))
    add_equation!(eqs, PA ~ IfElse.ifelse(LE < INITLE, INITPA, INITPA + LEEPA * (LE + EPA - INITLE)))
    delay_n!(eqs, BIRTHS, RT_PASS20, LV_PASS20, 20, ORDER)
    add_equation!(eqs, PASS20 ~ RT_PASS20[ORDER])
    delay_n!(eqs, PASS20, RT_PASS40, LV_PASS40, 20, ORDER)
    add_equation!(eqs, PASS40 ~ RT_PASS40[ORDER])
    delay_n!(eqs, PASS40, RT_PASS60, LV_PASS60, 20, ORDER)
    add_equation!(eqs, PASS60 ~ RT_PASS60[ORDER])
    add_equation!(eqs, PGR ~ BIRTHR - DEATHR)
    add_equation!(eqs, POP ~ A0020 + A2040 + A4060 + A60PL)
    add_equation!(eqs, PW ~ OP / A20PA)
    add_equation!(eqs, WELE ~ IfElse.ifelse(t > 2022, max(0, 1 + OWELE * (OW / OW2022 - 1)), 1))

    # public
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
    add_equation!(eqs, PPP ~ max(0, 1 + IPRVPSS * log(IPR / IPR1980)))
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

    # wellbeing
    add_equation!(eqs, AWBDI ~ exp(DRDI + log(WDI / TDI)))
    add_equation!(eqs, AWBIN ~ 1 + IEAWBIF * (INEQ / TI - 1))
    add_equation!(eqs, AWBGW ~ max(MWBGW, min(1, 1 + GWEAWBGWF * (PWA / TW - 1))))
    add_equation!(eqs, AWBI ~ (0.5 * AWBDI + 0.5 * AWBPS) * AWBIN * AWBGW * AWBP)
    add_equation!(eqs, AWBP ~ (1 + PREAWBF * (ORP - TPR)) * WBEP)
    add_equation!(eqs, AWBPS ~ exp(DRPS + log(PSP / TPS)))
    add_equation!(eqs, IEST ~ WorldDynamics.interpolate(INEQ / AI, tables[:IEST], ranges[:IEST]))
    add_equation!(eqs, IPP ~ IfElse.ifelse(EIPF > 0, EIP, RD))
    add_equation!(eqs, IRD ~ NRD * STRERD * STEERD)
    add_equation!(eqs, IST ~ PSESTR * IEST)
    smooth!(eqs, ORP, ((AWBI - PAWBI) / AWBI) / AWBPD, AWBPD)
    smooth!(eqs, PAWBI, AWBI, AWBPD)
    add_equation!(eqs, PSSGDP ~ PSP / GDPP)
    add_equation!(eqs, PSESTR ~ WorldDynamics.interpolate(PSSGDP / SPS, tables[:PSESTR], ranges[:PSESTR]))
    smooth!(eqs, RD, IRD, TCRD)
    add_equation!(eqs, STE ~ 1 + PESTF * (ORP - AP))
    add_equation!(eqs, STEERD ~ 1 + STEERDF * (STE / INITSTE - 1))
    add_equation!(eqs, STRERD ~ 1 + STRERDF * (SOTR / INITSOTR - 1))
    smooth!(eqs, SOTR, IST, TEST)
    add_equation!(eqs, WBEP ~ 1 + PAEAWBF * (LPR / THPA - 1))

    return ODESystem(eqs; name=name)
end
