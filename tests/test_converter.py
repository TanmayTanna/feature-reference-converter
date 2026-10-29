import csv
import gzip
import tempfile
import unittest
from pathlib import Path

from update_feature_csv import convert_row, parse_gtf, update_csv


GTF = (
    'chr1\tsynthetic\tgene\t1\t10\t.\t+\t.\tgene_name "GeneA"; gene_id "GENE001";\n'
    'chr1\tsynthetic\tgene\t20\t30\t.\t+\t.\tgene_id "GENE002"; gene_name "GeneB";\n'
    'chr1\tsynthetic\tgene\t40\t50\t.\t+\t.\tgene_id "GENE003"; gene_name "Shared";\n'
    'chr1\tsynthetic\tgene\t60\t70\t.\t+\t.\tgene_id "GENE004"; gene_name "Shared";\n'
)


class ConverterTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.gtf = self.root / "reference.gtf"
        self.gtf.write_text(GTF)
        self.ids, self.names = parse_gtf(self.gtf)

    def row(self, gene_id, gene_name):
        return {
            "id": "feature_1", "sequence": "ACGTACGT",
            "target_gene_id": gene_id, "target_gene_name": gene_name,
        }

    def test_exact_and_name_update(self):
        row, status = convert_row(self.row("GENE001", "GeneA"), self.ids, self.names)
        self.assertEqual(status, "unchanged")
        row, status = convert_row(self.row("GENE001", "OldName"), self.ids, self.names)
        self.assertEqual(row["target_gene_name"], "GeneA")
        self.assertEqual(status, "name_updated")

    def test_unique_symbol_updates_id(self):
        row, status = convert_row(self.row("OLD_ID", "GeneB"), self.ids, self.names)
        self.assertEqual(row["target_gene_id"], "GENE002")
        self.assertEqual(status, "id_updated")

    def test_control_and_ignore(self):
        for label in ["Non-Targeting", "intergenic", "negative_control"]:
            row, status = convert_row(self.row(label, label), self.ids, self.names)
            self.assertEqual(status, "control")
            self.assertEqual(row["target_gene_id"], "Non-Targeting")
        row, status = convert_row(self.row("Ignore", "Ignore"), self.ids, self.names)
        self.assertEqual(status, "ignored")
        self.assertEqual(row["target_gene_name"], "Ignore")

    def test_no_substring_control_matching_or_symbol_case_guessing(self):
        for gene_id, gene_name in [("OLD", "ControlLikeGene"), ("OLD", "genea")]:
            row, status = convert_row(self.row(gene_id, gene_name), self.ids, self.names)
            self.assertEqual(status, "unresolved")

    def test_ambiguous_symbol_not_guessed(self):
        row, status = convert_row(self.row("OLD", "Shared"), self.ids, self.names)
        self.assertEqual(status, "ambiguous")
        self.assertEqual(row["target_gene_id"], "NA")

    def test_gzip_and_empty_gtf(self):
        compressed = self.root / "reference.gtf.gz"
        with gzip.open(compressed, "wt") as handle:
            handle.write(GTF)
        self.assertEqual(parse_gtf(compressed)[0], self.ids)
        empty = self.root / "empty.gtf"
        empty.write_text("")
        with self.assertRaises(ValueError):
            parse_gtf(empty)

    def test_preserves_other_fields_and_refuses_overwrite(self):
        source, output = self.root / "input.csv", self.root / "output.csv"
        original = self.row("OLD", "GeneB")
        with source.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(original))
            writer.writeheader()
            writer.writerow(original)
        counts = update_csv(source, output, self.ids, self.names, strict=True)
        with output.open(newline="") as handle:
            converted = list(csv.DictReader(handle))[0]
        self.assertEqual(counts["id_updated"], 1)
        self.assertEqual(list(converted), list(original))
        self.assertEqual(converted["sequence"], original["sequence"])
        self.assertEqual(converted["id"], original["id"])
        with self.assertRaises(ValueError):
            update_csv(source, output, self.ids, self.names)
        with self.assertRaises(ValueError):
            update_csv(source, source, self.ids, self.names)

    def test_strict_mode_does_not_write_unresolved_output(self):
        source, output = self.root / "input.csv", self.root / "output.csv"
        source.write_text("target_gene_id,target_gene_name\nUNKNOWN,Unknown\n")
        with self.assertRaises(ValueError):
            update_csv(source, output, self.ids, self.names, strict=True)
        self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
