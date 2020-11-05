@test Psat(300.0) ≈ 0.353_658_941E-2
@test Psat(300.0u"K") ≈ 0.353_658_941E-2u"MPa"

@test Psat(500.0) ≈ 0.263_889_776E1
@test Psat(500.0u"K") ≈ 0.263_889_776E1u"MPa"

@test Psat(600.0) ≈ 0.123_443_146E2
@test Psat(600.0u"K") ≈ 0.123_443_146E2u"MPa"
