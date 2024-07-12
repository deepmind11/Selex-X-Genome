import math
import sqlite3

# def test_get_probe_count(my_countTable):
#     counts = my_countTable.get_probe_count()
#     assert counts[0] == 1330386


# def test_update_database(my_countTable):
#     my_countTable.update_database()

#     with sqlite3.connect(
#         "/burg/home/hg2604/hblab/Projects/Selex-X-Genome/database/Selex_X_Genome.db"
#     ) as conn:

#         cursor = conn.cursor()

#         cursor.execute(
#             """
#             SELECT "r1_file", "probe_count", "r0_count", "r1_count" FROM "count_tables" WHERE "r1_file" = ?
#             """,
#             (1,),
#         )
#         results = cursor.fetchone()
#         r1_file, probe_count, r0_count, r1_count = results

#     assert probe_count == 1330386
#     assert r0_count > 0
#     assert r1_count > 0


def test_score(my_countTable_sample, my_mono_motif):
    my_countTable_sample.score(my_mono_motif)

    with sqlite3.connect(
        "/burg/home/hg2604/hblab/Projects/Selex-X-Genome/database/Selex_X_Genome.db"
    ) as conn:

        cursor = conn.cursor()

        cursor.execute(
            """
            SELECT "type", "tf", "organism", "count_table_id", "score", "r0_count", "r1_count", "enrichment" FROM "motif" WHERE "count_table_id" = ?
            """,
            (1,),
        )
        results = cursor.fetchone()
        # print(results)
        type1, tf, organism, count_table_id, score, r0_count, r1_count, enrichment = (
            results
        )

    assert type1 == "Mononucleotide"
    assert tf == "CTCF"
    assert organism == "9606"
    assert count_table_id == 1
    assert math.isclose(score, -7010.382055, rel_tol=1e-4, abs_tol=1e-4)
    assert r0_count == 44
    assert r1_count == 977
    assert math.isclose(enrichment, 22.204545, rel_tol=1e-4, abs_tol=1e-4)
