# **SAFE-TUS: Signal Artifact Filtering for EEG with FUS**

## Authors:
** **Xavier Corominas-Terue**<sup>^</sup>,...  & **Lasse Christiansen**<sup>#</sup>

Danish Research Center For Magnetic Resonance, DRCMR 2025.

<sup>^</sup> First authorship 
<sup>#</sup> Senior Authorship

Please see for more details: 

**Reproducibility Instructions**
This repository provides an example preprocessing script for cleaning EEG recorded during concurrent transcranial focused ultrasound (tFUS–EEG). The script offers a simple yet effective approach to suppress short-latency artifacts that occur at the onset and offset of each tFUS pulse.


**Artifact profiles**

In systems with embedded hardware high-pass filtering (e.g., ~0.016 Hz), slow trends are removed and tFUS artifacts typically appear as sharp transitions at pulse onset/offset. This is the scenario assumed in the example.

**When you may need to adapt the code**

Depending on your acquisition hardware and built-in filters, tFUS pulse artifacts may instead present as DC step-like shift.

If your recording system lacks hardware high-pass filtering (or uses a very low cutoff), adding a detrend or high-pass filter before the artifact-suppression step to handle DC steps appropriately may be recommended. The provided script is a starting point optimized for systems with a ~0.016 Hz hardware high-pass; please adjust the preprocessing (e.g., high-pass cutoff and artifact handling) to match your device’s filtering and the artifact morphology you observe.



