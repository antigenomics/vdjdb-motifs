from __future__ import annotations

import importlib
import shutil
import sys
import tempfile
import types
import unittest
from pathlib import Path

import pandas as pd


def _load_io_module():
    stub_compat = types.ModuleType("redcea.vdjdb_redcea.compat")

    class _BackgroundTransform:
        pass

    stub_compat.BackgroundTransform = _BackgroundTransform
    sys.modules["redcea.vdjdb_redcea.compat"] = stub_compat
    return importlib.import_module("redcea.vdjdb_redcea.io")


IO_MODULE = _load_io_module()
build_processed_airr_from_tcremp_representations = IO_MODULE.build_processed_airr_from_tcremp_representations
prepare_output_dirs = IO_MODULE.prepare_output_dirs


class ProcessedAirrTests(unittest.TestCase):
    def test_prepare_output_dirs_creates_processed_airr_directory(self):
        temp_dir = Path(tempfile.mkdtemp(prefix="vdjdb_motifs_test_"))
        try:
            paths = prepare_output_dirs(temp_dir)
            self.assertTrue(paths.airr_dir.is_dir())
            self.assertTrue(paths.airr_processed_dir.is_dir())
            self.assertTrue(paths.tcremp_dir.is_dir())
            self.assertEqual(paths.airr_processed_dir, temp_dir / "airr_format_processed")
        finally:
            shutil.rmtree(temp_dir)

    def test_build_processed_airr_uses_processed_representations_and_filters_invalid_rows(self):
        representations = pd.DataFrame(
            {
                "cdr3aa_beta": ["CASSLGQETQYF", "BADSEQ", "CASSQETQYF", None],
                "v_beta": ["TRBV7-9", "TRBV5-1", "TRBV6-5", "TRBV2"],
                "j_beta": ["TRBJ2-5", "TRBJ1-2", None, "TRBJ2-7"],
                "clone_id": [101, 102, 103, 104],
            }
        )

        result = build_processed_airr_from_tcremp_representations(representations, "TRB")

        expected = pd.DataFrame(
            {
                "junction_aa": ["CASSLGQETQYF"],
                "v_call": ["TRBV7-9"],
                "j_call": ["TRBJ2-5"],
                "locus": ["beta"],
            }
        )
        pd.testing.assert_frame_equal(result.reset_index(drop=True), expected)

    def test_build_processed_airr_falls_back_to_generic_column_names(self):
        representations = pd.DataFrame(
            {
                "junction_aa": ["CAVRDTDKLIF", "NOTCANONICAL"],
                "v_call": ["TRAV1-2", "TRAV12-2"],
                "j_call": ["TRAJ34", "TRAJ24"],
            }
        )

        result = build_processed_airr_from_tcremp_representations(representations, "TRA")

        self.assertEqual(list(result.columns), ["junction_aa", "v_call", "j_call", "locus"])
        self.assertEqual(result.to_dict(orient="records"), [{"junction_aa": "CAVRDTDKLIF", "v_call": "TRAV1-2", "j_call": "TRAJ34", "locus": "alpha"}])


if __name__ == "__main__":
    unittest.main()
