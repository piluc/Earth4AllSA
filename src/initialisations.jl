
#variables: finance - climate - demand - energy - foodland - inventory - labourmarket - other 
_inits = Dict{Symbol,Float64}(
    # climate
    :N2OA => 1.052,
    :CH4A => 2.5,
    :CO2A => 2600,
    :ISCEGA => 12,
    :PWA => 0.4,
    :EHS => 0,
    # demand
    :ETF2022 => 0, # No extra taxes before 2022
    :FGBW => 0.3, # Equal to fraction transfer in 1980
    :POCI => 7081, # It was OCI in 1980
    :WD => 7406.88,
    :PWCIN => 13000,
    :PGCIN => 5400,
    :GNI => 6531.07,
    # energy
    :EEPI2022 => 1,
    :REC => 300,
    :ACSWCF1980 => 10,
    :FEC => 980,
    # finance
    :CCSD => 0.04, # Taken from Vensim table
    :ELTI => 0.02,
    :PI => 0.02,
    :PU => 0.0326951,
    :CBSR => 0.02,
    # foodland
    :BALA => 3000,
    :CRLA => 1450,
    :GRLA => 3300,
    :FOLA => 1100,
    :OGFA => 2600,
    :SQICA => 1,
    :URLA => 215,
    # inventory
    :DELDI => 1,
    :EPP => 28087, # Taken from Vensim table
    :INV => 11234.8, # Taken from Vensim table
    :PRIN => 1, # Taken from Vensim table
    :PRI => 1,
    :RS => 28087,
    :SSWI => 1,
    # labourmarket
    :ECLR => 41,
    :ILPR => 0.8, # Taken from Vensim table
    :LAUS => 3060, # Taken from Vensim table
    :NHW => 2,
    :PURA => 0.05,
    :WARA => 3.6715, # Taken from Vensim table
    :WEOCLR => 1, # Taken from Vensim table
    :WF => 1530,
    :WSO => 0.5,
    # other
    :PGDPP => 6.4 * 0.93,
    # output
    :CUCPIS => 10072.5, # Taken from Vensim table
    :CUCPUS => 909.5, # Taken from Vensim table
    :ETFP => 1,
    :FACNC => 1.05149, # Taken from Vensim table
    :LAUS => 3060, # Taken from Labour and market sector
    :OLY => 26497.1, # Taken from Vensim table
    :WSO => 0.5, # Taken from Labour and market sector
    # population
    :A0020 => 2170,
    :A2040 => 1100,
    :A4060 => 768,
    :A60PL => 382,
    :DEATHS => 30,
    :EGDPP => 6.4,
    :EPA => 0,
    :LE => 67,
    :PA => 62,
    :PASS20 => 100,
    :PASS40 => 64,
    :PASS60 => 38,
    # public
    :RTFPUA => 0, # Taken from Vensim table
    :TFPEE5TA => 1,
    # wellbeing
    :ORP => 0,
    :PAWBI => 0.65,
    :RD => 30,
    :STE => 1.3,
    :SOTR => 0.6,
)

getinitialisations() = copy(_inits)
