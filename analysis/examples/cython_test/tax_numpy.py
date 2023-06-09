def calculate_pop_tax(income):
    # Personfradrag, http://www.skatteetaten.no/no/Tabeller-og-satser/Personfradrag/
    # Minstefradrag, http://www.skatteetaten.no/no/Tabeller-og-satser/Minstefradrag/
    tax_income = income - 53150 - 31800
    tax_income[tax_income < 0] = 0

    # Alminnelig inntekt, http://www.skatteetaten.no/no/Tabeller-og-satser/Alminnelig-inntekt/
    tax = 0.24 * tax_income

    # Trinnskatt, http://www.skatteetaten.no/no/Tabeller-og-satser/trinnskatt/
    idx = tax_income >= 934_050
    tax[idx] += 0.1452 * (tax_income[idx] - 934_050) + 49761.155
    idx = (tax_income < 934_050) & (tax_income >= 580_650)
    tax[idx] += 0.1152 * (tax_income[idx] - 580_650) + 9049.475
    idx = (tax_income < 580_650) & (tax_income >= 230_950)
    tax[idx] += 0.0241 * (tax_income[idx] - 230_950) + 621.705
    idx = (tax_income < 230_950) & (tax_income >= 164_100)
    tax[idx] += 0.0093 * (tax_income[idx] - 164_100)

    return tax
