CREATE TABLE IF NOT EXISTS "experiments" (
    "id" INTEGER,
    "accession" TEXT UNIQUE NOT NULL,
    "control" TEXT CHECK("response" IN ('y', 'n', NULL)),
    "life_stage_age" TEXT,
    "perturbed" TEXT CHECK("response" IN ('y', 'n', NULL)),
    "lab" TEXT,
    "biosample_class" TEXT,
    "developmental_slims" TEXT,
    "system_slims" TEXT,
    "organ_slims" TEXT,
    "cell_slims" TEXT,
    PRIMARY KEY ("id")
);


CREATE TABLE IF NOT EXISTS "libraries" (
    "id" INTEGER,
    "accession" TEXT UNIQUE NOT NULL,
    "antibody" TEXT,
    "biosample" TEXT,
    "technical_rep_number" INTEGER,
    "biological_rep_number" INTEGER,
    "experiment_id" INTEGER NOT NULL,
    FOREIGN KEY ("experiment_id") REFERENCES "experiments" ("id"),
    PRIMARY KEY ("id")
);

CREATE TABLE IF NOT EXISTS "files" (
    "id" INTEGER,
    "accession" TEXT UNIQUE NOT NULL,
    "read_length" INTEGER,
    "experiment_id" INTEGER NOT NULL,
    "library_id" INTEGER NOT NULL,
    "paired_end" TEXT CHECK("response" IN ('y', 'n', NULL)),
    FOREIGN KEY ("experiment_id") REFERENCES "experiments" ("id"),
    FOREIGN KEY ("library_id") REFERENCES "library" ("id"),
    PRIMARY KEY ("id")
);

CREATE TABLE IF NOT EXISTS "count_tables" (
    "id" INTEGER,
    "r1_file" TEXT NOT NULL UNIQUE,
    "probe_count" INTEGER NOT NULL,
    "r0_count" INTEGER NOT NULL,
    "r1_count" INTEGER NOT NULL,
    FOREIGN KEY ("r1_file") REFERENCES "files" ("id"),
    PRIMARY KEY ("id")
);


CREATE TABLE IF NOT EXISTS "motif" (
    "id" INTEGER,
    "type" TEXT NOT NULL,
    "tf" TEXT NOT NULL,
    "organism" TEXT NOT NULL,
    "count_table_id" INTEGER NOT NULL,
    "score" NUMERIC,
    "r0_count" INTEGER,
    "r1_count" INTEGER,
    "enrichment" NUMERIC,
    "correlation" NUMERIC,
    FOREIGN KEY ("count_table_id") REFERENCES "count_tables" ("id"),
    PRIMARY KEY ("id"),
    UNIQUE ("tf", "organism", "count_table_id")
);