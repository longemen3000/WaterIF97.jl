@test mysignif(SpecificG(0.1, 0.399_517_097E3), TestDigits) ≈ mysignif(SpecificG_Ps(0.1, 7.5), TestDigits)
@test mysignif(SpecificG(0.1u"MPa", 0.399_517_097E3u"K"), TestDigits) ≈ mysignif(SpecificG_Ps(0.1u"MPa", 7.5u"kJ/kg/K"), TestDigits)

@test mysignif(SpecificG(0.1, 0.514_127_081E3), TestDigits) ≈ mysignif(SpecificG_Ps(0.1, 8.0), TestDigits)
@test mysignif(SpecificG(0.1u"MPa", 0.514_127_081E3u"K"), TestDigits) ≈ mysignif(SpecificG_Ps(0.1u"MPa", 8.0u"kJ/kg/K"), TestDigits)

@test mysignif(SpecificG(2.5, 0.103_984_917E4), TestDigits) ≈ mysignif(SpecificG_Ps(2.5, 8.0), TestDigits)
@test mysignif(SpecificG(2.5u"MPa", 0.103_984_917E4u"K"), TestDigits) ≈ mysignif(SpecificG_Ps(2.5u"MPa", 8.0u"kJ/kg/K"), TestDigits)

@test mysignif(SpecificF(0.1, 0.399_517_097E3), TestDigits) ≈ mysignif(SpecificF_Ps(0.1, 7.5), TestDigits)
@test mysignif(SpecificF(0.1u"MPa", 0.399_517_097E3u"K"), TestDigits) ≈ mysignif(SpecificF_Ps(0.1u"MPa", 7.5u"kJ/kg/K"), TestDigits)

@test mysignif(SpecificF(0.1, 0.514_127_081E3), TestDigits) ≈ mysignif(SpecificF_Ps(0.1, 8.0), TestDigits)
@test mysignif(SpecificF(0.1u"MPa", 0.514_127_081E3u"K"), TestDigits) ≈ mysignif(SpecificF_Ps(0.1u"MPa", 8.0u"kJ/kg/K"), TestDigits)

@test mysignif(SpecificF(2.5, 0.103_984_917E4), TestDigits) ≈ mysignif(SpecificF_Ps(2.5, 8.0), TestDigits)
@test mysignif(SpecificF(2.5u"MPa", 0.103_984_917E4u"K"), TestDigits) ≈ mysignif(SpecificF_Ps(2.5u"MPa", 8.0u"kJ/kg/K"), TestDigits)

@test SpecificV(0.1, 0.399_517_097E3) ≈ SpecificV_Ps(0.1, 7.5)
@test SpecificV(0.1u"MPa", 0.399_517_097E3u"K") ≈ SpecificV_Ps(0.1u"MPa", 7.5u"kJ/kg/K")

@test SpecificV(0.1, 0.514_127_081E3) ≈ SpecificV_Ps(0.1, 8.0)
@test SpecificV(0.1u"MPa", 0.514_127_081E3u"K") ≈ SpecificV_Ps(0.1u"MPa", 8.0u"kJ/kg/K")

@test SpecificV(2.5, 0.103_984_917E4) ≈ SpecificV_Ps(2.5, 8.0)
@test SpecificV(2.5u"MPa", 0.103_984_917E4u"K") ≈ SpecificV_Ps(2.5u"MPa", 8.0u"kJ/kg/K")

@test SpecificU(0.1, 0.399_517_097E3) ≈ SpecificU_Ps(0.1, 7.5)
@test SpecificU(0.1u"MPa", 0.399_517_097E3u"K") ≈ SpecificU_Ps(0.1u"MPa", 7.5u"kJ/kg/K")

@test SpecificU(0.1, 0.514_127_081E3) ≈ SpecificU_Ps(0.1, 8.0)
@test SpecificU(0.1u"MPa", 0.514_127_081E3u"K") ≈ SpecificU_Ps(0.1u"MPa", 8.0u"kJ/kg/K")

@test SpecificU(2.5, 0.103_984_917E4) ≈ SpecificU_Ps(2.5, 8.0)
@test SpecificU(2.5u"MPa", 0.103_984_917E4u"K") ≈ SpecificU_Ps(2.5u"MPa", 8.0u"kJ/kg/K")

@test SpecificH(0.1, 0.399_517_097E3) ≈ SpecificH_Ps(0.1, 7.5)
@test SpecificH(0.1u"MPa", 0.399_517_097E3u"K") ≈ SpecificH_Ps(0.1u"MPa", 7.5u"kJ/kg/K")

@test SpecificH(0.1, 0.514_127_081E3) ≈ SpecificH_Ps(0.1, 8.0)
@test SpecificH(0.1u"MPa", 0.514_127_081E3u"K") ≈ SpecificH_Ps(0.1u"MPa", 8.0u"kJ/kg/K")

@test SpecificH(2.5, 0.103_984_917E4) ≈ SpecificH_Ps(2.5, 8.0)
@test SpecificH(2.5u"MPa", 0.103_984_917E4u"K") ≈ SpecificH_Ps(2.5u"MPa", 8.0u"kJ/kg/K")

@test SpecificCP(0.1, 0.399_517_097E3) ≈ SpecificCP_Ps(0.1, 7.5)
@test SpecificCP(0.1u"MPa", 0.399_517_097E3u"K") ≈ SpecificCP_Ps(0.1u"MPa", 7.5u"kJ/kg/K")

@test SpecificCP(0.1, 0.514_127_081E3) ≈ SpecificCP_Ps(0.1, 8.0)
@test SpecificCP(0.1u"MPa", 0.514_127_081E3u"K") ≈ SpecificCP_Ps(0.1u"MPa", 8.0u"kJ/kg/K")

@test SpecificCP(2.5, 0.103_984_917E4) ≈ SpecificCP_Ps(2.5, 8.0)
@test SpecificCP(2.5u"MPa", 0.103_984_917E4u"K") ≈ SpecificCP_Ps(2.5u"MPa", 8.0u"kJ/kg/K")

@test SpecificCV(0.1, 0.399_517_097E3) ≈ SpecificCV_Ps(0.1, 7.5)
@test SpecificCV(0.1u"MPa", 0.399_517_097E3u"K") ≈ SpecificCV_Ps(0.1u"MPa", 7.5u"kJ/kg/K")

@test SpecificCV(0.1, 0.514_127_081E3) ≈ SpecificCV_Ps(0.1, 8.0)
@test SpecificCV(0.1u"MPa", 0.514_127_081E3u"K") ≈ SpecificCV_Ps(0.1u"MPa", 8.0u"kJ/kg/K")

@test SpecificCV(2.5, 0.103_984_917E4) ≈ SpecificCV_Ps(2.5, 8.0)
@test SpecificCV(2.5u"MPa", 0.103_984_917E4u"K") ≈ SpecificCV_Ps(2.5u"MPa", 8.0u"kJ/kg/K")

@test SpeedOfSound(0.1, 0.399_517_097E3) ≈ SpeedOfSound_Ps(0.1, 7.5)
@test SpeedOfSound(0.1u"MPa", 0.399_517_097E3u"K") ≈ SpeedOfSound_Ps(0.1u"MPa", 7.5u"kJ/kg/K")

@test SpeedOfSound(0.1, 0.514_127_081E3) ≈ SpeedOfSound_Ps(0.1, 8.0)
@test SpeedOfSound(0.1u"MPa", 0.514_127_081E3u"K") ≈ SpeedOfSound_Ps(0.1u"MPa", 8.0u"kJ/kg/K")

@test SpeedOfSound(2.5, 0.103_984_917E4) ≈ SpeedOfSound_Ps(2.5, 8.0)
@test SpeedOfSound(2.5u"MPa", 0.103_984_917E4u"K") ≈ SpeedOfSound_Ps(2.5u"MPa", 8.0u"kJ/kg/K")
