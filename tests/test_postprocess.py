"""Tests for the post-processing results writer module."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from opensees_pipeline.postprocess import write_results


class TestWriteResults:
    """Tests for write_results() output file generation."""

    def _make_dummy_data(self, n_nodes: int = 5) -> tuple[
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        np.ndarray,
        dict,
        dict,
    ]:
        """Create minimal valid inputs for write_results."""
        n_elem = n_nodes - 1
        x = np.linspace(0, 10, n_nodes)
        ux = np.zeros(n_nodes)
        uy = np.random.default_rng(42).uniform(-0.01, 0.01, n_nodes)
        uz = np.zeros(n_nodes)
        my_i = np.random.default_rng(42).uniform(-100, 100, n_elem)
        my_j = np.random.default_rng(42).uniform(-100, 100, n_elem)
        n_force = np.random.default_rng(42).uniform(-50, 50, n_elem)
        s_ax = np.random.default_rng(42).uniform(-10, 10, n_nodes)
        s_vt = np.random.default_rng(42).uniform(-20, 20, n_nodes)
        trib = np.full(n_nodes, 10.0 / n_nodes)
        pipe = {"outer_R": 0.5, "t_wall": 0.025, "E": 210e9, "nu": 0.3}
        meta = {"unit_system": "SI", "analysis_type": "test"}
        return x, ux, uy, uz, my_i, my_j, n_force, s_ax, s_vt, trib, pipe, meta

    def test_creates_all_files(self, tmp_path: Path) -> None:
        """write_results creates all four expected output files."""
        args = self._make_dummy_data()
        write_results(str(tmp_path), *args)

        assert (tmp_path / "pipe_nodes.csv").exists()
        assert (tmp_path / "pipe_elements.csv").exists()
        assert (tmp_path / "spring_forces.csv").exists()
        assert (tmp_path / "model_metadata.json").exists()

    def test_nodes_csv_shape(self, tmp_path: Path) -> None:
        """pipe_nodes.csv has correct number of rows and columns."""
        n_nodes = 7
        args = self._make_dummy_data(n_nodes)
        write_results(str(tmp_path), *args)

        df = pd.read_csv(tmp_path / "pipe_nodes.csv")
        assert len(df) == n_nodes
        assert list(df.columns) == ["node_index", "x_m", "UX_m", "UY_m", "UZ_m"]

    def test_elements_csv_shape(self, tmp_path: Path) -> None:
        """pipe_elements.csv has correct number of rows and columns."""
        n_nodes = 7
        args = self._make_dummy_data(n_nodes)
        write_results(str(tmp_path), *args)

        df = pd.read_csv(tmp_path / "pipe_elements.csv")
        assert len(df) == n_nodes - 1
        expected = ["elem_index", "x_i_m", "x_j_m", "N_N", "My_i_Nm", "My_j_Nm"]
        assert list(df.columns) == expected

    def test_springs_csv_shape(self, tmp_path: Path) -> None:
        """spring_forces.csv has correct shape."""
        n_nodes = 7
        args = self._make_dummy_data(n_nodes)
        write_results(str(tmp_path), *args)

        df = pd.read_csv(tmp_path / "spring_forces.csv")
        assert len(df) == n_nodes

    def test_metadata_json_keys(self, tmp_path: Path) -> None:
        """model_metadata.json contains required keys."""
        args = self._make_dummy_data()
        write_results(str(tmp_path), *args)

        with open(tmp_path / "model_metadata.json") as fh:
            meta = json.load(fh)

        assert "timestamp" in meta
        assert meta["pipe_od_m"] == 1.0  # 2 * 0.5
        assert meta["pipe_t_m"] == 0.025
        assert meta["pipe_E_Pa"] == 210e9
        assert meta["unit_system"] == "SI"
        assert meta["analysis_type"] == "test"

    def test_creates_output_dir(self, tmp_path: Path) -> None:
        """write_results creates the output directory if it does not exist."""
        out = tmp_path / "nested" / "deep"
        args = self._make_dummy_data()
        write_results(str(out), *args)

        assert out.exists()
        assert (out / "pipe_nodes.csv").exists()

    def test_extra_metadata_preserved(self, tmp_path: Path) -> None:
        """Extra keys in model_metadata are written to JSON."""
        args = list(self._make_dummy_data())
        args[-1] = {"unit_system": "SI", "analysis_type": "test", "custom_key": 42}
        write_results(str(tmp_path), *args)

        with open(tmp_path / "model_metadata.json") as fh:
            meta = json.load(fh)

        assert meta["custom_key"] == 42
