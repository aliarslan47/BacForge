"""Modül kayıt defteri (registry). Sequence = M00 -- M18 pipeline flow.
"""
from __future__ import annotations

from .m00_input import InputDetectionModule
from .m01_read_qc import ReadQCModule
from .m02_taxonomic_qc import TaxonomicQCModule
from .m03_assembly import AssemblyModule
from .m04_polishing_qc import PolishingGenomeQCModule
from .m05_species_ref import SpeciesReferenceIdentificationModule
from .m06_annotation import GenomeAnnotationModule
from .m07_typing import StrainTypingModule
from .m08_amr import AMRModule
from .m09_virulence import VirulenceModule
from .m10_plasmid import PlasmidModule
from .m11_mge import MobileGeneticElementsModule
from .m12_phage_crispr import PhageCRISPRDefenseModule
from .m13_variants import VariantMutationModule
from .m14_context import GenomicContextModule
from .m15_comparative import ComparativeGenomicsModule
from .m16_phylogenomics import PhylogenomicsModule
from .m17_statistics import StatisticsVisualizationModule
from .m18_report import FinalReportExportModule

REGISTRY = [
    InputDetectionModule,                   # M00
    ReadQCModule,                           # M01
    TaxonomicQCModule,                      # M02
    AssemblyModule,                         # M03
    PolishingGenomeQCModule,                # M04
    SpeciesReferenceIdentificationModule,   # M05
    GenomeAnnotationModule,                 # M06
    StrainTypingModule,                     # M07
    AMRModule,                              # M08
    VirulenceModule,                        # M09
    PlasmidModule,                          # M10
    MobileGeneticElementsModule,            # M11
    PhageCRISPRDefenseModule,               # M12
    VariantMutationModule,                  # M13
    GenomicContextModule,                   # M14
    ComparativeGenomicsModule,              # M15
    PhylogenomicsModule,                    # M16
    StatisticsVisualizationModule,          # M17
    FinalReportExportModule,                # M18
]
