
def helium_content(burnup):
    # Assume helium production is proportional to burnup
    return 1.5e-3 * burnup  # Moles of He per ton of fuel

def yield_stress_embrittlement(yield_stress_original, helium_content):
    k = 1e-5  # Linear reduction factor for helium embrittlement
    return yield_stress_original * (1 - k * helium_content)