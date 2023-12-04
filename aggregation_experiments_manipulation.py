from pdb_manipulation import (average_files, flatten_pdb, randomize_pdb,
                              rotate_pdb)

if __name__ == "__main__":
    # rotate
    rotate_pdb(
        "aggregation_test/up_template/2upa.pdb",
        "aggregation_test/up_template_rotated/2upa.pdb",
        "1",
    )
    rotate_pdb(
        "aggregation_test/low_template/2low.pdb",
        "aggregation_test/low_template_rotated/2low.pdb",
        "2",
    )
    average_files(
        "aggregation_test/up_template_rotated/2upa.pdb",
        "aggregation_test/low_template_rotated/2low.pdb",
        "aggregation_test/aggregated_rotated/2agg.pdb",
    )
    # randomize
    randomize_pdb(
        "aggregation_test/low_template/2low.pdb",
        "aggregation_test/low_template_randomized/2low.pdb",
        20,
        None,
        100,
    )
    randomize_pdb(
        "aggregation_test/up_template/2upa.pdb",
        "aggregation_test/up_template_randomized/2upa.pdb",
        20,
        100,
        None,
    )
    average_files(
        "aggregation_test/up_template_randomized/2upa.pdb",
        "aggregation_test/low_template_randomized/2low.pdb",
        "aggregation_test/aggregated_randomized/2agg.pdb",
    )
    rotate_pdb(
        "aggregation_test/up_template_randomized/2upa.pdb",
        "aggregation_test/up_template_randomized/2upa.pdb",
        "1",
    )
    rotate_pdb(
        "aggregation_test/low_template_randomized/2low.pdb",
        "aggregation_test/low_template_randomized/2low.pdb",
        "2",
    )
    # flatten
    ## flatten up
    flatten_pdb(
        "aggregation_test/full_template/2ubq.pdb",
        "aggregation_test/up_template_flattened/2upa.pdb",
        2,
        100,
        None,
    )
    flatten_pdb(
        "aggregation_test/full_template/2ubq.pdb",
        "aggregation_test/low_template_flattened/2low.pdb",
        2,
        None,
        100,
    )
    average_files(
        "aggregation_test/low_template_flattened/2low.pdb",
        "aggregation_test/up_template_flattened/2upa.pdb",
        "aggregation_test/aggregated_flattened/2agg.pdb",
    )
    rotate_pdb(
        "aggregation_test/up_template_flattened/2upa.pdb",
        "aggregation_test/up_template_flattened/2upa.pdb",
        "1",
    )
    rotate_pdb(
        "aggregation_test/low_template_flattened/2low.pdb",
        "aggregation_test/low_template_flattened/2low.pdb",
        "2",
    )
    # RFDiffusion
    average_files(
        "aggregation_test/low_template_rfdiff/2low.pdb",
        "aggregation_test/up_template_rfdiff/2upa.pdb",
        "aggregation_test/aggregated_rfdiff/2agg.pdb",
    )
    rotate_pdb(
        "aggregation_test/up_template_rfdiff/2upa.pdb",
        "aggregation_test/up_template_rfdiff/2upa.pdb",
        "1",
    )
    rotate_pdb(
        "aggregation_test/low_template_rfdiff/2low.pdb",
        "aggregation_test/low_template_rfdiff/2low.pdb",
        "2",
    )