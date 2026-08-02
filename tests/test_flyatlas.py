"""Tests for FlyAtlas expression parsing."""

from dataset_finder.clients.flyatlas import FlyAtlasClient

SAMPLE_TABLE = """FlyBase ID\tFBgn0000114
Annotation Symbol\tCG31762
Symbol\tbru1
Name\tbruno 1

\tAdult Male\t\t\tAdult Female\t\t\tMale v. Female\t\tLarval
Tissue\tFPKM\tSD\tEnrichment\tFPKM\tSD\tEnrichment\tM/F\tp value\tFPKM\tSD\tEnrichment
Head\t14.33\t2.2\t0.43\t14.8\t0.52\t0.12\t0.97\tn.s.\t-\t-\t-
Brain / CNS\t9.89\t0.78\t0.3\t8.54\t1.12\t0.07\t1.16\tn.s.\t4.56\t1.07\t1.14
Hindgut\t37.64\t4.21\t1.14\t48.37\t4.78\t0.4\t0.78\tp > 0.05\t71.39\t19.93\t17.8
Ovary\t-\t-\t-\t180.06\t8.06\t1.48\t-\t-\t-\t-\t-
Testis\t210.17\t66.27\t6.37\t-\t-\t-\t-\t-\t-\t-\t-
Whole body\t32.97\t3.71\t1\t121.92\t10.58\t1\t0.27\tp > 0.01\t4.01\t0.22\t1
"""


def test_parse_flyatlas_gene_table() -> None:
    expression = FlyAtlasClient.parse_table(SAMPLE_TABLE)

    assert expression.flybase_id == "FBgn0000114"
    assert expression.symbol == "bru1"
    assert expression.brain_male_fpkm == 9.89
    assert expression.brain_female_fpkm == 8.54
    assert expression.brain_larval_fpkm == 4.56
    assert expression.head_male_fpkm == 14.33
    assert expression.head_female_fpkm == 14.8
    assert expression.top_male_tissue == "Testis"
    assert expression.top_male_fpkm == 210.17
    assert expression.top_female_tissue == "Ovary"
    assert expression.top_female_fpkm == 180.06
    assert expression.top_larval_tissue == "Hindgut"
    assert expression.top_larval_fpkm == 71.39
