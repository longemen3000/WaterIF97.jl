@test mysignif(SpecificG(3.0, 0.307_842_258E3), TestDigits) ≈ mysignif(SpecificG_Ps(3.0, 0.5), TestDigits)
@test mysignif(SpecificG(3.0u"MPa", 0.307_842_258E3u"K"), TestDigits) ≈ mysignif(SpecificG_Ps(3.0u"MPa", 0.5u"kJ/kg/K"), TestDigits)

@test mysignif(SpecificG(80.0, 0.309_979_785E3), TestDigits) ≈ mysignif(SpecificG_Ps(80.0, 0.5), TestDigits)
@test mysignif(SpecificG(80.0u"MPa", 0.309_979_785E3u"K"), TestDigits) ≈ mysignif(SpecificG_Ps(80.0u"MPa", 0.5u"kJ/kg/K"), TestDigits)

@test mysignif(SpecificG(80.0, 0.565_899_909E3), TestDigits) ≈ mysignif(SpecificG_Ps(80.0, 3.0), TestDigits)
@test mysignif(SpecificG(80.0u"MPa", 0.565_899_909E3u"K"), TestDigits) ≈ mysignif(SpecificG_Ps(80.0u"MPa", 3.0u"kJ/kg/K"), TestDigits)

@test mysignif(SpecificF(3.0, 0.307_842_258E3), TestDigits) ≈ mysignif(SpecificF_Ps(3.0, 0.5), TestDigits)
@test mysignif(SpecificF(3.0u"MPa", 0.307_842_258E3u"K"), TestDigits) ≈ mysignif(SpecificF_Ps(3.0u"MPa", 0.5u"kJ/kg/K"), TestDigits)

@test mysignif(SpecificF(80.0, 0.309_979_785E3), TestDigits) ≈ mysignif(SpecificF_Ps(80.0, 0.5), TestDigits)
@test mysignif(SpecificF(80.0u"MPa", 0.309_979_785E3u"K"), TestDigits) ≈ mysignif(SpecificF_Ps(80.0u"MPa", 0.5u"kJ/kg/K"), TestDigits)

@test mysignif(SpecificF(80.0, 0.565_899_909E3), TestDigits) ≈ mysignif(SpecificF_Ps(80.0, 3.0), TestDigits)
@test mysignif(SpecificF(80.0u"MPa", 0.565_899_909E3u"K"), TestDigits) ≈ mysignif(SpecificF_Ps(80.0u"MPa", 3.0u"kJ/kg/K"), TestDigits)

@test SpecificV(3.0, 0.307_842_258E3) ≈ SpecificV_Ps(3.0, 0.5)
@test SpecificV(3.0u"MPa", 0.307_842_258E3u"K") ≈ SpecificV_Ps(3.0u"MPa", 0.5u"kJ/kg/K")

@test SpecificV(80.0, 0.309_979_785E3) ≈ SpecificV_Ps(80.0, 0.5)
@test SpecificV(80.0u"MPa", 0.309_979_785E3u"K") ≈ SpecificV_Ps(80.0u"MPa", 0.5u"kJ/kg/K")

@test SpecificV(80.0, 0.565_899_909E3) ≈ SpecificV_Ps(80.0, 3.0)
@test SpecificV(80.0u"MPa", 0.565_899_909E3u"K") ≈ SpecificV_Ps(80.0u"MPa", 3.0u"kJ/kg/K")

@test SpecificU(3.0, 0.307_842_258E3) ≈ SpecificU_Ps(3.0, 0.5)
@test SpecificU(3.0u"MPa", 0.307_842_258E3u"K") ≈ SpecificU_Ps(3.0u"MPa", 0.5u"kJ/kg/K")

@test SpecificU(80.0, 0.309_979_785E3) ≈ SpecificU_Ps(80.0, 0.5)
@test SpecificU(80.0u"MPa", 0.309_979_785E3u"K") ≈ SpecificU_Ps(80.0u"MPa", 0.5u"kJ/kg/K")

@test SpecificU(80.0, 0.565_899_909E3) ≈ SpecificU_Ps(80.0, 3.0)
@test SpecificU(80.0u"MPa", 0.565_899_909E3u"K") ≈ SpecificU_Ps(80.0u"MPa", 3.0u"kJ/kg/K")

@test SpecificH(3.0, 0.307_842_258E3) ≈ SpecificH_Ps(3.0, 0.5)
@test SpecificH(3.0u"MPa", 0.307_842_258E3u"K") ≈ SpecificH_Ps(3.0u"MPa", 0.5u"kJ/kg/K")

@test SpecificH(80.0, 0.309_979_785E3) ≈ SpecificH_Ps(80.0, 0.5)
@test SpecificH(80.0u"MPa", 0.309_979_785E3u"K") ≈ SpecificH_Ps(80.0u"MPa", 0.5u"kJ/kg/K")

@test SpecificH(80.0, 0.565_899_909E3) ≈ SpecificH_Ps(80.0, 3.0)
@test SpecificH(80.0u"MPa", 0.565_899_909E3u"K") ≈ SpecificH_Ps(80.0u"MPa", 3.0u"kJ/kg/K")

@test SpecificCP(3.0, 0.307_842_258E3) ≈ SpecificCP_Ps(3.0, 0.5)
@test SpecificCP(3.0u"MPa", 0.307_842_258E3u"K") ≈ SpecificCP_Ps(3.0u"MPa", 0.5u"kJ/kg/K")

@test SpecificCP(80.0, 0.309_979_785E3) ≈ SpecificCP_Ps(80.0, 0.5)
@test SpecificCP(80.0u"MPa", 0.309_979_785E3u"K") ≈ SpecificCP_Ps(80.0u"MPa", 0.5u"kJ/kg/K")

@test SpecificCP(80.0, 0.565_899_909E3) ≈ SpecificCP_Ps(80.0, 3.0)
@test SpecificCP(80.0u"MPa", 0.565_899_909E3u"K") ≈ SpecificCP_Ps(80.0u"MPa", 3.0u"kJ/kg/K")

@test SpecificCV(3.0, 0.307_842_258E3) ≈ SpecificCV_Ps(3.0, 0.5)
@test SpecificCV(3.0u"MPa", 0.307_842_258E3u"K") ≈ SpecificCV_Ps(3.0u"MPa", 0.5u"kJ/kg/K")

@test SpecificCV(80.0, 0.309_979_785E3) ≈ SpecificCV_Ps(80.0, 0.5)
@test SpecificCV(80.0u"MPa", 0.309_979_785E3u"K") ≈ SpecificCV_Ps(80.0u"MPa", 0.5u"kJ/kg/K")

@test SpecificCV(80.0,0.565_899_909E3) ≈ SpecificCV_Ps(80.0, 3.0)
@test SpecificCV(80.0u"MPa", 0.565_899_909E3u"K") ≈ SpecificCV_Ps(80.0u"MPa", 3.0u"kJ/kg/K")

@test SpeedOfSound(3.0, 0.307_842_258E3) ≈ SpeedOfSound_Ps(3.0, 0.5)
@test SpeedOfSound(3.0u"MPa", 0.307_842_258E3u"K") ≈ SpeedOfSound_Ps(3.0u"MPa", 0.5u"kJ/kg/K")

@test SpeedOfSound(80.0, 0.309_979_785E3) ≈ SpeedOfSound_Ps(80.0, 0.5)
@test SpeedOfSound(80.0u"MPa", 0.309_979_785E3u"K") ≈ SpeedOfSound_Ps(80.0u"MPa", 0.5u"kJ/kg/K")

@test SpeedOfSound(80.0, 0.565_899_909E3) ≈ SpeedOfSound_Ps(80.0, 3.0)
@test SpeedOfSound(80.0u"MPa", 0.565_899_909E3u"K") ≈ SpeedOfSound_Ps(80.0u"MPa", 3.0u"kJ/kg/K")
