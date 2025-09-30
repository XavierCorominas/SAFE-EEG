# **SAFE-TUS: Signal Artifact Filtering for concurrent tFUS-EEG**

**Authors:**
** **Xavier Corominas-Teruel**<sup>^</sup>,...  & **Lasse Christiansen**<sup>#</sup>

Danish Research Center For Magnetic Resonance, DRCMR 2025.

<sup>^</sup> First authorship (xavi@drcmr.dk) 
<sup>#</sup> Senior Authorship


<img width="1301" height="904" alt="PictureEG" src="https://github.com/user-attachments/assets/725cae0d-cd09-41dc-a968-27c636908341" />


**Dependencies**

The current processing script employ functions from EEGlab (https://eeglab.org/), and the TESA toolbox (https://github.com/nigelrogasch/TESA).  Please install them locally before executing the SAFE-EEG processing.

  

**Reproducibility Instructions**

This repository provides an example preprocessing script for cleaning EEG signals recorded during concurrent low-intensity transcranial focused ultrasound and EEG (tFUS–EEG). The script offers a simple yet effective approach to suppress short-latency artifacts that occur at the onset and offset of each tFUS pulse.


  
**Artifact profiles**

In systems with embedded hardware high-pass filtering (e.g., ~0.016 Hz), slow trends are removed and tFUS artifacts typically appear as sharp transitions at pulse onset/offset. This is the scenario assumed in the example.

  
**When you may need to adapt the code:**

Depending on your acquisition hardware and built-in filters, tFUS pulse artifacts may instead present as DC step-like shift. If your recording system lacks hardware high-pass filtering (or uses a very low cutoff), adding a detrend or high-pass filter before the artifact-suppression step to handle DC steps appropriately may be recommended. The provided script is a starting point optimized for systems with a ~0.016 Hz hardware high-pass; please adjust the preprocessing (e.g., high-pass cutoff and artifact handling) to match your device’s filtering and the artifact morphology you observe.



**Further considerations**
The current script focuses on removing pulse-like artifacts. However, in realistic human experiments, additional artifacts are common. It is recommended extending the pipeline—if needed—to address ocular activity (blinks/saccades), muscular contamination, line noise (50/60 Hz and harmonics), motion/electrode drift, and other physiological or non-physiological sources that could confound the EEG.



**Processing EEG steps**
The pipeline includes minimal processing steps: 
(1) Epoching the data centered at the onset of the first pulse of the pulse train.
(2) De-meaning.
(3) Removal and cubic interpolation of segments containing sharp pulse artifacts (from -1.5ms to +1.5ms around each artifact, depending on EEG sampling frequency).
(4) Band-pass filtering with a 2nd order zero-phase Butterworth filter, baseline correction (-210, -10ms) and re-referencing to the average.



