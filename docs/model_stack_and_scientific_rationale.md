# Model Stack And Scientific Rationale

## Why This Platform Makes Sense

This repository is built around a simple scientific idea:

preformulation triage for a new API is first a `local interaction ranking` problem, not a final bulk-property prediction problem.

That means the platform does not need to predict the full final formulation recipe in one step. It needs to do three things well:

- recognize what kind of API it is
- propose chemically plausible formulation families
- rank those families using local mechanistic evidence

That is exactly what the current model stack is designed to do.

## The Model Stack

The current platform uses two model layers plus one execution layer.

### 1. Descriptor And Prior Layer

Used for:

- fast unseen-API intake
- API class assignment
- first-pass family recommendation

Current components:

- RDKit descriptor extraction
- rule-based API classification
- curated formulation-family prior table

What it can evaluate:

- molecular weight
- cLogP
- topological polar surface area
- H-bond donor and acceptor counts
- rotatable bond count
- aromaticity
- formal charge and charge-state context
- water-affinity, aggregation-risk, and crystallization-risk proxies

Why this layer supports the platform:

- formulation family recommendation depends strongly on polarity, flexibility, ionization, and hydrogen-bonding capacity
- these are features that can be estimated directly from structure
- for early triage, interpretable descriptors are often more useful than a black-box end-to-end predictor

This layer is implemented in:

- `scripts/extract_api_descriptors.py`
- `scripts/classify_api.py`
- `scripts/generate_candidate_matrix.py`

## 2. Mechanistic Screening Layer

Used for:

- local cluster relaxation
- interaction ranking
- shortlist generation for DFT and wet-lab follow-up

Current model:

- checkpoint-based FAIRChem molecular interatomic potential through `FAIRChemCalculator`
- current archived relaxed GPU run uses:
  - `task_name = "omol"`
  - `device = "cuda"`
  - `relax_steps = 25`

The checkpoint path is intentionally configurable:

- `checkpoints/best_inference_ckpt.pt`

What this model can do in the repo:

- predict energies and forces for local molecular/polymer clusters
- relax initial structures
- score local API-polymer, API-API, polymer-polymer, and solvent competition patterns
- support replicate screening with randomized initial placements

Why this layer supports the platform:

- amorphous stabilization, dispersion, and local compatibility are controlled by short-range interactions
- these are exactly the types of local geometric and energetic questions a molecular MLIP can probe
- cluster relaxation reduces dependence on arbitrary initial placements
- replicate screening makes ranking more robust

This layer is implemented in:

- `scripts/run_preformulation_mechanistic_screen.py`

## 3. Execution And Aggregation Layer

Used for:

- user-facing API input
- end-to-end orchestration
- output explanation

Current components:

- `run_unseen_api_recommendation.py`
- `scripts/run_unseen_api_recommendation.py`
- `scripts/explain_results.py`

What this layer does:

- converts user input into a structured API payload
- runs descriptor extraction and family recommendation
- writes the candidate matrix
- optionally launches mechanistic screening
- packages outputs into:
  - `summary.md`
  - `ranking.csv`
  - `candidate_matrix_generated.csv`
  - `buy_list.txt`
  - `lab_test_plan.txt`

## What The Platform Is Actually Modeling

The platform is strongest when the decision is governed by:

- API-polymer association versus API self-association
- solvent competition around the API and polymer
- local polymer cohesion
- charge-state sensitivity
- shortlist-level ranking before expensive validation

The platform is not directly modeling:

- full bulk rheology
- final printability window
- exact dissolution profile
- long-time physical stability
- final formulation percentages

That is why the platform should be described as:

- recommendation and screening platform
- preformulation triage platform
- shortlist generator for DFT and lab evaluation

not as:

- final formulation predictor
- printability predictor
- full formulation design engine

## Why The Existing Results Support The Platform

### Evidence 1: archived mechanistic screening is nontrivial

In the relaxed GPU archive:

- the paracetamol neutral candidates changed strongly after relaxation
- the top candidate remained `PVP + water`
- the middle ranking changed after relaxation

This is important because it shows:

- the mechanistic layer is doing meaningful structural work
- the platform is not just echoing input ordering
- relaxation affects decision quality

### Evidence 2: methotrexate recommendation is chemically structured

In the methotrexate family recommendation archive:

- the neutral branch and ionic branch do not collapse to the same family ordering
- the neutral shortlist is led by `PVP + water`
- the monoanionic shortlist is led by `HPC + ethanol/water`

This is important because it shows:

- the recommendation layer is sensitive to charge-state context
- the platform is producing mechanistically interpretable families rather than random combinations

## Current Honest Claim

The honest claim is:

this platform can intake a new API, generate plausible formulation families, rank them using structure-based priors and mechanistic local screening, and produce a shortlist for downstream DFT and lab evaluation.

That claim is already supported by the tracked code, archived outputs, and the current evidence figures.
