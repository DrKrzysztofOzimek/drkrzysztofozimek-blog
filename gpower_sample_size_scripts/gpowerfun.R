library(pwr)

sample_size <- pwr::pwr.t.test(d = 0.6, power = 0.9, sig.level = 0.05, type = "paired", alternative = "greater")
sample_size$n
