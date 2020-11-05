@test mysignif(SpecificG(40.0, 0.743_056_411E3), TestDigits) ≈ mysignif(SpecificG_Ph(40.0, 2700.0), TestDigits)
@test mysignif(SpecificG(40.0u"MPa", 0.743_056_411E3u"K"), TestDigits) ≈ mysignif(SpecificG_Ph(40.0u"MPa", 2700.0u"kJ/kg"), TestDigits)

@test mysignif(SpecificG(60.0, 0.791_137_067E3), TestDigits) ≈ mysignif(SpecificG_Ph(60.0, 2700.0), TestDigits)
@test mysignif(SpecificG(60.0u"MPa", 0.791_137_067E3u"K"), TestDigits) ≈ mysignif(SpecificG_Ph(60.0u"MPa", 2700.0u"kJ/kg"), TestDigits)

@test mysignif(SpecificG(60.0, 0.882_756_860E3), TestDigits) ≈ mysignif(SpecificG_Ph(60.0, 3200.0), TestDigits)
@test mysignif(SpecificG(60.0u"MPa", 0.882_756_860E3u"K"), TestDigits) ≈ mysignif(SpecificG_Ph(60.0u"MPa", 3200.0u"kJ/kg"), TestDigits)

@test mysignif(SpecificF(40.0, 0.743_056_411E3), TestDigits) ≈ mysignif(SpecificF_Ph(40.0, 2700.0), TestDigits)
@test mysignif(SpecificF(40.0u"MPa", 0.743_056_411E3u"K"), TestDigits) ≈ mysignif(SpecificF_Ph(40.0u"MPa", 2700.0u"kJ/kg"), TestDigits)

@test mysignif(SpecificF(60.0, 0.791_137_067E3), TestDigits) ≈ mysignif(SpecificF_Ph(60.0, 2700.0), TestDigits)
@test mysignif(SpecificF(60.0u"MPa", 0.791_137_067E3u"K"), TestDigits) ≈ mysignif(SpecificF_Ph(60.0u"MPa", 2700.0u"kJ/kg"), TestDigits)

@test mysignif(SpecificF(60.0, 0.882_756_860E3), TestDigits) ≈ mysignif(SpecificF_Ph(60.0, 3200.0), TestDigits)
@test mysignif(SpecificF(60.0u"MPa", 0.882_756_860E3u"K"), TestDigits) ≈ mysignif(SpecificF_Ph(60.0u"MPa", 3200.0u"kJ/kg"), TestDigits)

@test SpecificV(40.0, 0.743_056_411E3) ≈ SpecificV_Ph(40.0, 2700.0)
@test SpecificV(40.0u"MPa", 0.743_056_411E3u"K") ≈ SpecificV_Ph(40.0u"MPa", 2700.0u"kJ/kg")

@test SpecificV(60.0, 0.791_137_067E3) ≈ SpecificV_Ph(60.0, 2700.0)
@test SpecificV(60.0u"MPa", 0.791_137_067E3u"K") ≈ SpecificV_Ph(60.0u"MPa", 2700.0u"kJ/kg")

@test SpecificV(60.0, 0.882_756_860E3) ≈ SpecificV_Ph(60.0, 3200.0)
@test SpecificV(60.0u"MPa", 0.882_756_860E3u"K") ≈ SpecificV_Ph(60.0u"MPa", 3200.0u"kJ/kg")

@test SpecificU(40.0, 0.743_056_411E3) ≈ SpecificU_Ph(40.0, 2700.0)
@test SpecificU(40.0u"MPa", 0.743_056_411E3u"K") ≈ SpecificU_Ph(40.0u"MPa", 2700.0u"kJ/kg")

@test SpecificU(60.0, 0.791_137_067E3) ≈ SpecificU_Ph(60.0, 2700.0)
@test SpecificU(60.0u"MPa", 0.791_137_067E3u"K") ≈ SpecificU_Ph(60.0u"MPa", 2700.0u"kJ/kg")

@test SpecificU(60.0, 0.882_756_860E3) ≈ SpecificU_Ph(60.0, 3200.0)
@test SpecificU(60.0u"MPa", 0.882_756_860E3u"K") ≈ SpecificU_Ph(60.0u"MPa", 3200.0u"kJ/kg")

@test SpecificS(40.0, 0.743_056_411E3) ≈ SpecificS_Ph(40.0, 2700.0)
@test SpecificS(40.0u"MPa", 0.743_056_411E3u"K") ≈ SpecificS_Ph(40.0u"MPa", 2700.0u"kJ/kg")

@test SpecificS(60.0, 0.791_137_067E3) ≈ SpecificS_Ph(60.0, 2700.0)
@test SpecificS(60.0u"MPa", 0.791_137_067E3u"K") ≈ SpecificS_Ph(60.0u"MPa", 2700.0u"kJ/kg")

@test SpecificS(60.0, 0.882_756_860E3) ≈ SpecificS_Ph(60.0, 3200.0)
@test SpecificS(60.0u"MPa", 0.882_756_860E3u"K") ≈ SpecificS_Ph(60.0u"MPa", 3200.0u"kJ/kg")

@test SpecificCP(40.0, 0.743_056_411E3) ≈ SpecificCP_Ph(40.0, 2700.0)
@test SpecificCP(40.0u"MPa", 0.743_056_411E3u"K") ≈ SpecificCP_Ph(40.0u"MPa", 2700.0u"kJ/kg")

@test SpecificCP(60.0, 0.791_137_067E3) ≈ SpecificCP_Ph(60.0, 2700.0)
@test SpecificCP(60.0u"MPa", 0.791_137_067E3u"K") ≈ SpecificCP_Ph(60.0u"MPa", 2700.0u"kJ/kg")

@test SpecificCP(60.0, 0.882_756_860E3) ≈ SpecificCP_Ph(60.0, 3200.0)
@test SpecificCP(60.0u"MPa", 0.882_756_860E3u"K") ≈ SpecificCP_Ph(60.0u"MPa", 3200.0u"kJ/kg")

@test SpecificCV(40.0, 0.743_056_411E3) ≈ SpecificCV_Ph(40.0, 2700.0)
@test SpecificCV(40.0u"MPa", 0.743_056_411E3u"K") ≈ SpecificCV_Ph(40.0u"MPa", 2700.0u"kJ/kg")

@test SpecificCV(60.0, 0.791_137_067E3) ≈ SpecificCV_Ph(60.0, 2700.0)
@test SpecificCV(60.0u"MPa", 0.791_137_067E3u"K") ≈ SpecificCV_Ph(60.0u"MPa", 2700.0u"kJ/kg")

@test SpecificCV(60.0, 0.882_756_860E3) ≈ SpecificCV_Ph(60.0, 3200.0)
@test SpecificCV(60.0u"MPa", 0.882_756_860E3u"K") ≈ SpecificCV_Ph(60.0u"MPa", 3200.0u"kJ/kg")

@test SpeedOfSound(40.0, 0.743_056_411E3) ≈ SpeedOfSound_Ph(40.0, 2700.0)
@test SpeedOfSound(40.0u"MPa", 0.743_056_411E3u"K") ≈ SpeedOfSound_Ph(40.0u"MPa", 2700.0u"kJ/kg")

@test SpeedOfSound(60.0, 0.791_137_067E3) ≈ SpeedOfSound_Ph(60.0, 2700.0)
@test SpeedOfSound(60.0u"MPa", 0.791_137_067E3u"K") ≈ SpeedOfSound_Ph(60.0u"MPa", 2700.0u"kJ/kg")

@test SpeedOfSound(60.0, 0.882_756_860E3) ≈ SpeedOfSound_Ph(60.0, 3200.0)
@test SpeedOfSound(60.0u"MPa", 0.882_756_860E3u"K") ≈ SpeedOfSound_Ph(60.0u"MPa", 3200.0u"kJ/kg")
