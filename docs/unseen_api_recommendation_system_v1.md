# Unseen API Recommendation System v1

## Goal

Build a practical recommendation workflow that starts from an unseen API and returns a ranked shortlist of formulation families for lab screening.

This system is intended to support:

- family-level excipient prioritization
- mechanistic preformulation screening
- DFT shortlist definition
- purchase and wet-lab planning

It is not intended to output a final formulation recipe in v1.

## Core Principle

The system should not jump directly from `API SMILES` to a final formulation.

The correct order is:

1. define the API chemistry
2. compute API-level descriptors and charge-state options
3. map the API into formulation-relevant categories
4. assign family-level formulation priors
5. run replicate mechanistic screening on shortlisted families
6. aggregate the screening outputs into a stable ranking
7. send the shortlist to DFT and then experiment

## End-To-End Workflow

```text
API SMILES
  -> descriptor extraction
  -> charge/salt state expansion
  -> API class assignment
  -> family-level prior recommendation
  -> mechanistic cluster generation
  -> replicate relaxation and scoring
  -> aggregated ranking
  -> DFT shortlist
  -> experimental shortlist
```

## System Inputs

### Required

- `api_name`
- `api_smiles`

### Optional but strongly recommended

- `known_charge_state`
- `salt_form`
- `pH_context`
- `formulation_context`
- `route_of_administration`
- `processing_context`

### Example v1 input

```yaml
api_name: methotrexate
api_smiles: CN(C)C1=NC(=NC2=CC(=CC=C21)C(=O)NCCC(=O)O)NCC3=CN=CC=C3
known_charge_state: [neutral, anionic]
salt_form: sodium_optional
pH_context: topical_or_semi_solid_unknown
formulation_context: semi_solid_printing
route_of_administration: topical
processing_context: drying_after_extrusion
```

## Layer 1: API Feature Extraction

Descriptors should be computed automatically from SMILES and charge-state candidates.

### Core molecular descriptors

- molecular weight
- heavy atom count
- cLogP
- topological polar surface area
- H-bond donor count
- H-bond acceptor count
- rotatable bond count
- aromatic ring count
- formal charge
- fraction of hetero atoms

### Formulation-relevant proxy features

- ionizable group count
- candidate protonation states
- acidic/basic/amphoteric flag
- rigidity or flexibility proxy
- polarity balance proxy
- water-affinity proxy
- self-aggregation risk proxy
- crystallization risk proxy

### Important boundary

`dissolution` should not be a direct input feature in v1.

Reason:

- dissolution is a higher-level experimental endpoint
- it depends on polymer, solvent, microstructure, and processing
- it cannot be treated as a reliable direct function of API SMILES alone

Instead, v1 should use descriptor-level proxies that influence dissolution indirectly.

## Layer 2: Charge-State and Salt-State Expansion

Each new API should expand into a small set of chemically plausible screening states before family recommendation.

### Required outputs

- neutral state if chemically plausible
- major ionic state under the intended pH context
- common salt-associated state if relevant

### Why this matters

- polymer compatibility can change strongly with ionization
- solvent competition can change strongly with charge
- Carbopol-like or polyelectrolyte systems are especially sensitive to ionic state

The output of this layer is a set of `API states`, not yet formulation candidates.

## Layer 3: API Class Assignment

The system should assign the API into formulation-relevant classes based on descriptors and charge-state behavior.

### Example v1 classes

- hydrophilic neutral
- poorly soluble neutral
- acidic or anionic
- basic or cationic
- amphiphilic
- highly flexible polar
- aromatic aggregation-prone

This layer can begin as a rule-based classifier before moving to a trained model.

### Example rule logic

- high `cLogP` and low `TPSA` -> poor-solubility route
- strong acidity and dominant anionic state -> anionic route
- high aromaticity and low flexibility -> aggregation-risk route
- high H-bonding capacity and moderate polarity -> polymer-binding candidate

## Layer 4: Family-Level Recommendation

This is the first recommendation layer.

The system should recommend families, not exact recipes.

### Output examples

- `PVP + water-like family`
- `HPMC/HPC + ethanol/water-like family`
- `Carbopol-like + glycerol/water-like family`
- `PEG-containing co-solvent family`

### Inputs to family recommendation

- API descriptors
- API class assignment
- charge-state context
- formulation context
- literature or internal co-occurrence priors

### Priority rule

The system should recommend `2-5` families for mechanistic screening, not the whole space.

### Why this layer exists

It keeps the mechanistic screen focused on chemically plausible formulation routes instead of brute-force enumeration.

## Layer 5: Literature and Data Priors

Family recommendation should be informed by an explicit prior table.

### Recommended prior table fields

- `api_class`
- `formulation_context`
- `polymer_family`
- `solvent_family`
- `co_solvent_family`
- `prior_strength`
- `evidence_type`
- `evidence_source`
- `notes`

### Example uses

- acidic APIs may receive higher prior weight for cellulose or gel-network routes
- neutral poorly soluble APIs may receive higher prior weight for PVP-like amorphous stabilization routes
- semi-solid printing context may boost high-viscosity cellulose or carbomer-like families

In v1, this can be a manually curated CSV or markdown-backed ruleset.

## Layer 6: Mechanistic Screening

This is where the current repository already has working infrastructure.

For each shortlisted family, build local candidate systems such as:

- `API + polymer fragment`
- `API + polymer fragment + solvent`
- `API + polymer fragment + co-solvent cluster`
- `API + API + polymer`
- `polymer + polymer`
- `API(ionized) + polymer + solvent` when relevant

### Mechanistic questions to score

- does the API prefer polymer association or self-association
- does solvent preserve or disrupt API-polymer contacts
- does the local polymer environment keep useful cohesion
- does ionic state change the ranking
- are there local hotspot or instability signals

### Current v1 metrics

- interaction energy proxy
- API-polymer contact count
- API-solvent contact count
- API-API contact count
- polymer-polymer contact count
- minimum API-polymer distance
- relaxation score shift

## Layer 7: Replicate Screening

Every shortlisted candidate should be run with multiple initial placements and random seeds.

### Required v1 settings

- `n_replicates = 5` as default
- `n_replicates = 10` for a higher-confidence follow-up
- random rotation
- random initial placement offsets
- same relaxation protocol across replicates

### Required outputs

- mean final score
- score standard deviation
- mean post interaction energy
- contact statistics
- mean rank
- rank standard deviation
- top-1 frequency
- top-2 frequency

This is the minimum needed to convert a single-shot scoring script into a stable recommendation layer.

## Layer 8: Aggregation Logic

The final recommendation should not rely on one opaque number.

### Recommended v1 ranking view

- `family_prior_score`
- `mechanistic_mean_score`
- `mechanistic_score_sd`
- `rank_stability`
- `top1_frequency`
- `top2_frequency`
- `screening_confidence_tier`

### Proposed interpretation rule

- high priority: good mean score and stable rank across replicates
- medium priority: good mean score but unstable rank
- low priority: weak mean score or collapse after relaxation

### Confidence tier examples

- `A`: stable across replicates and consistent with DFT spot check
- `B`: stable across replicates but not yet DFT checked
- `C`: single-screen or high-variance candidate
- `X`: ionic or reference-uncertain branch, not cross-comparable yet

## Layer 9: Final Outputs

The v1 system should return a human-readable shortlist plus machine-readable tables.

### Human-readable output

For unseen API `X`, recommend:

1. `PVP + water-like family`
   reason: strongest local stabilization and low aggregation tendency in replicate screen
2. `HPMC/HPC + ethanol/water-like family`
   reason: good compatibility and plausible printability route
3. `Carbopol-like + glycerol/water-like family`
   reason: useful gel-network route worth experimental testing

### Machine-readable outputs

- ranked family table
- per-candidate replicate table
- aggregated candidate table
- DFT shortlist table
- suggested purchase list
- suggested lab validation plan

## Suggested Output Schema

### `family_recommendations.csv`

- `api_name`
- `api_state`
- `formulation_context`
- `family_id`
- `polymer_family`
- `solvent_family`
- `co_solvent_family`
- `family_prior_score`
- `mechanistic_mean_score`
- `mechanistic_score_sd`
- `mean_rank`
- `rank_sd`
- `top1_frequency`
- `top2_frequency`
- `confidence_tier`
- `recommended_action`

### `dft_shortlist.csv`

- `api_name`
- `candidate_id`
- `family_id`
- `selection_reason`
- `mlip_rank`
- `confidence_tier`

## Validation Ladder

The correct validation order is:

1. replicate mechanistic screen
2. DFT spot checks on top, middle, and weak candidates
3. wet-lab shortlist confirmation

### DFT questions

- is the relative ordering preserved
- are the relaxed contact motifs chemically reasonable
- does the API remain associated with polymer
- does the aggregation control behave as expected

### Wet-lab questions

- miscibility and amorphous stabilization
- rheology and printability
- drying behavior
- release behavior
- chemical stability

## What v1 Should Not Claim

The system should not claim that it can directly predict:

- exact formulation percentages
- final printability window
- final dissolution curve
- long-term storage success

The proper v1 framing is:

`unseen API family recommendation and mechanistic shortlist generation`

## Immediate Build Targets For This Repo

To turn this design into an implemented system, the next concrete repo tasks should be:

1. add an `API descriptor extraction` script driven by SMILES
2. add a `family prior` table under `data/preformulation/`
3. add a rule-based `API class assignment` layer
4. connect family recommendation to `run_preformulation_mechanistic_screen.py`
5. output a family-level recommendation summary for neutral candidates first
6. validate on `2-3` real unseen APIs, including methotrexate

## Recommended v1 Milestone

The first complete milestone should be:

`Input unseen API SMILES -> generate descriptors and charge-state candidates -> recommend 2-5 formulation families -> run replicate neutral mechanistic screen -> aggregate the ranking -> output DFT and lab shortlist`
