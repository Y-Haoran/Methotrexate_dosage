# Family Recommendation Summary: methotrexate

## Input Context

- formulation_context: `semi_solid_printing`
- route_of_administration: `topical`
- processing_context: `drying_after_extrusion`
- pH_context: `topical_or_semi_solid_unknown`
- top_k_families_per_state: `4`

## API State Classification

### neutral

- label: `neutral`
- formal_charge: `0`
- primary_api_class: `hydrophilic_neutral`
- assigned_api_classes: `acidic_or_anionic,amphiphilic,highly_flexible_polar,hydrophilic_neutral`
- water_affinity_proxy: `0.85`
- self_aggregation_risk_proxy: `0.38`
- crystallization_risk_proxy: `0.32`
- descriptor_source_note: `parent neutral smiles`

### monoanion

- label: `monoanionic`
- formal_charge: `-1`
- primary_api_class: `acidic_or_anionic`
- assigned_api_classes: `acidic_or_anionic,amphiphilic,highly_flexible_polar`
- water_affinity_proxy: `0.95`
- self_aggregation_risk_proxy: `0.37`
- crystallization_risk_proxy: `0.32`
- descriptor_source_note: `parent smiles with metadata-only charge override`

### dianion

- label: `dianionic`
- formal_charge: `-2`
- primary_api_class: `acidic_or_anionic`
- assigned_api_classes: `acidic_or_anionic,amphiphilic,highly_flexible_polar`
- water_affinity_proxy: `1.00`
- self_aggregation_risk_proxy: `0.36`
- crystallization_risk_proxy: `0.32`
- descriptor_source_note: `parent smiles with metadata-only charge override`

## Top Family Recommendations

### neutral

- rank 1: `pvp + water` (score `4.06`, classes `amphiphilic,highly_flexible_polar,hydrophilic_neutral`)
  rationale: Amphiphilic APIs often benefit from a balanced polar-apolar local environment and strong local binder contacts.
- rank 2: `hpc + water` (score `3.95`, classes `hydrophilic_neutral`)
  rationale: Cellulose-ether aqueous route that usually tolerates polar APIs while preserving a printable semi-solid path.
- rank 3: `hpc + ethanol_water` (score `3.93`, classes `amphiphilic,highly_flexible_polar`)
  rationale: Flexible polar APIs often sit well in cellulose-ether routes with a controlled volatile co-solvent component.
- rank 4: `gelatin + water` (score `3.35`, classes `highly_flexible_polar`)
  rationale: Hydrogen-bond-rich peptide-like polymer comparator for flexible polar APIs.

### monoanion

- rank 1: `hpc + ethanol_water` (score `4.61`, classes `acidic_or_anionic,amphiphilic,highly_flexible_polar`)
  rationale: Cellulose-ether route is a strong first screen for acidic or anionic APIs in printable semi-solid contexts.
- rank 2: `pvp + water` (score `4.16`, classes `acidic_or_anionic,amphiphilic,highly_flexible_polar`)
  rationale: PVP remains a useful comparator for acidic APIs because it can still offer strong local H-bonding and dispersion support.
- rank 3: `xanthan + water` (score `3.60`, classes `acidic_or_anionic`)
  rationale: An anionic polysaccharide gel route is worth screening when the goal includes robust gel-network behavior.
- rank 4: `hpc + water` (score `3.40`, classes `acidic_or_anionic`)
  rationale: General fallback cellulose-ether route for acidic or anionic APIs.

### dianion

- rank 1: `hpc + ethanol_water` (score `4.61`, classes `acidic_or_anionic,amphiphilic,highly_flexible_polar`)
  rationale: Cellulose-ether route is a strong first screen for acidic or anionic APIs in printable semi-solid contexts.
- rank 2: `pvp + water` (score `4.16`, classes `acidic_or_anionic,amphiphilic,highly_flexible_polar`)
  rationale: PVP remains a useful comparator for acidic APIs because it can still offer strong local H-bonding and dispersion support.
- rank 3: `xanthan + water` (score `3.60`, classes `acidic_or_anionic`)
  rationale: An anionic polysaccharide gel route is worth screening when the goal includes robust gel-network behavior.
- rank 4: `hpc + water` (score `3.40`, classes `acidic_or_anionic`)
  rationale: General fallback cellulose-ether route for acidic or anionic APIs.

## Downstream Mechanistic Screen Matrix

- generated rows: `8`
- neutral and ionic states are both supported in the recommendation layer, but their mechanistic rankings should still be interpreted separately.
- descriptor-level charge-state handling is metadata-based unless an explicit `state_smiles` is supplied for the charged form.

