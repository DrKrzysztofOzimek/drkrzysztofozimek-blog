from statsmodels.stats.power import TTestPower

analysis = TTestPower()
sample_size = analysis.solve_power(
    effect_size=0.6, power=0.9, alpha=0.05, alternative="larger"
)
sample_size
