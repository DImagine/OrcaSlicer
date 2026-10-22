"""Check seam-exchange argument rejection using the built CLI, without loading a model."""

from pathlib import Path
import subprocess
import sys
import tempfile


def main():
    # This suite is skipped by unit-only CI installations that do not contain the app.
    binary = Path(sys.argv[1]).resolve()
    if not binary.is_file():
        return 77
    with tempfile.TemporaryDirectory(prefix="orca-ai-seam-cli-") as temporary:
        root = Path(temporary)
        cases = [
            ["--slice", "0", "--export-seam-loops", str(root / "export.json")],
            ["--slice", "1", "--import-seam-points", str(root / "response.json"),
             "--load-slicedata", str(root)],
        ]
        for index, arguments in enumerate(cases):
            output = root / str(index)
            output.mkdir()
            # The targeted rejection must precede project loading and G-code generation.
            result = subprocess.run(
                [str(binary), "--datadir", str(root / "data"), "--outputdir", str(output), *arguments],
                capture_output=True, text=True, timeout=60, check=False,
            )
            assert result.returncode != 0, result.stdout
            assert "AI seam exchange requires --slice N" in result.stdout + result.stderr, result.stderr
            assert not (root / "export.json").exists()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
