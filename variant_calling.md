# Variant calling

Variant calling on ONT reads Illumina-corrected with Ratatosk can be performed with the Pepper-MARGIN-DeepVariant pipeline. Pepper-MARGIN-DeepVariant models are now provided for ONT reads from R9.4 flowcells, basecalled with Guppy 5+ SUP model and corrected with Ratatosk 0.7.5. Using these models on reads from a different type of flowcell or basecalled with a different version of Guppy will work but the call accuracy will be suboptimal.

## Requirements

* [Pepper-MARGIN-DeepVariant](https://github.com/kishwarshafin/pepper)
* [Pepper-MARGIN-DeepVariant models for Ratatosk](https://drive.google.com/file/d/1AbkKIGY19xbnvVI6PUF_R4YhVOLeXiZw/view?usp=sharing)

## Installation

Decompress the Pepper-MARGIN-DeepVariant models:
```
tar -xvzf R9_GUPPY_SUP_MODELS.tar.gz
```

The output should be 5 files:
```
ls -lh
# R9_GUPPY_SUP_DEEPVARIANT.data-00000-of-00001
# R9_GUPPY_SUP_DEEPVARIANT.index
# R9_GUPPY_SUP_DEEPVARIANT.meta
# R9_GUPPY_SUP_PEPPER_HP.pkl
# R9_GUPPY_SUP_PEPPER_SNP.pkl
```

## Input

* Corrected long reads mapping (BAM file)
* Pepper-MARGIN-DeepVariant models for Ratatosk

The corrected long reads must be mapped to a reference genome. We recommend [Winnowmap2](https://github.com/marbl/Winnowmap) with the `-x map-pb` preset or [Minimap2](https://github.com/lh3/minimap2) with the `-x map-hifi` preset for the mapping. The BAM file must be sorted (`samtools sort`) and indexed (`samtools index`).

## Usage

```
singularity exec --bind /usr/lib/locale/ \
pepper_deepvariant_r0.7.sif \
run_pepper_margin_deepvariant call_variant \
-b "${INPUT_DIR}/${BAM}" \
-f "${INPUT_DIR}/${REF}" \
-o "${OUTPUT_DIR}" \
-p "${OUTPUT_PREFIX}" \
-s "${SAMPLE}"
-t ${THREADS} \
--ont_r9_guppy5_sup \
--dv_model "R9_GUPPY_SUP_DEEPVARIANT" \
--pepper_model R9_GUPPY_SUP_PEPPER_SNP.pkl \
--pepper_hp_model R9_GUPPY_SUP_PEPPER_HP.pkl \
```

The Docker command line should be similar to the Singularity one. See [Pepper-MARGIN-DeepVariant](https://github.com/kishwarshafin/pepper) for more information.
