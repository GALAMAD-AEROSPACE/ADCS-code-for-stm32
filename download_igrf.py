import requests
import json

# Download and save all IGRF coefficients for all epochs and SV

def download_all_igrf_coefficients():
    url = "https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf14coeffs.txt"
    response = requests.get(url)
    response.raise_for_status()
    lines = response.text.split('\n')

    # Find the header line with years
    header_line = None
    for line in lines:
        if 'g/h' in line:
            header_line = line
            break
    if header_line is None:
        raise ValueError("Could not find coefficient header line")

    print("Header line:", repr(header_line))

    header_parts = header_line.split()
    # Find all year columns and their indices
    year_indices = []
    for i, part in enumerate(header_parts):
        try:
            year = float(part)
            if 1900 <= year <= 2100:
                year_indices.append((year, i))
        except ValueError:
            continue
    # Find SV column (look for '2025-30' or similar)
    sv_index = None
    for i, part in enumerate(header_parts):
        if part.strip().replace(' ', '') == '2025-30':
            sv_index = i
            break
    if sv_index is None:
        raise ValueError("Could not find SV (2025-30) column")

    years = [y for y, _ in year_indices]
    coeffs = {'g': {}, 'h': {}}
    sv_year = years[-1]  # SV is for the last year

    # Parse coefficients
    start_parsing = False
    for line in lines:
        if line.startswith('#'):
            continue
        if not start_parsing:
            if 'g/h' in line:
                start_parsing = True
            continue
        if not line.strip():
            continue
        parts = line.split()
        if len(parts) < sv_index + 1:
            continue
        kind = parts[0]  # 'g' or 'h'
        n = int(parts[1])
        m = int(parts[2])
        key = f"{n},{m}"
        values = []
        for _, idx in year_indices:
            try:
                values.append(float(parts[idx]))
            except Exception:
                values.append(0.0)
        # Add SV as last value
        try:
            sv = float(parts[sv_index])
        except Exception:
            sv = 0.0
        values.append(sv)
        coeffs[kind][key] = values

    out = {
        "years": years,
        "SV_year": sv_year,
        "coeffs": coeffs
    }
    with open("igrf_coefficients_full.json", "w") as f:
        json.dump(out, f)
    print(f"Saved all IGRF coefficients for years {years[0]} to {years[-1]} and SV.")

if __name__ == "__main__":
    download_all_igrf_coefficients() 