#!/usr/bin/env python3
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Auto-generate Quarto markdown documentation for the vstol CLI.

Builds the same argparse parser tree that ``vstolib.cli.cli_main`` constructs
(``init_arg_parser`` plus the ``run`` subparser registrations), walks the
registered subparsers, and emits one ``.qmd`` per subcommand plus an
``index.qmd`` listing them all.

Usage:
    python scripts/docs/generate_docs.py

Output is written to ``<repo_root>/docs/cli`` regardless of the current
working directory. The script never deletes anything.

Per-subcommand pages (``<subcommand>.qmd``) are regenerated on every run, so
any manual edits to those files will be lost — edit the argparse parsers in
``src/vstolib/cli/`` instead.

The index page (``docs/cli/index.qmd``) is treated as user-editable: it is only
created the first time, and subsequent runs leave it alone. Edit it freely.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Iterable

# Prefer the in-repo source over any installed copy so the generated docs always
# reflect the working tree. vstol uses a src-layout and is typically installed
# non-editably (``pip install .``); without this, a stale site-packages copy
# would be documented instead. In CI this is a no-op (the source was just
# installed from this same tree).
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

import vstolib  # noqa: E402
from vstolib.constants import VariantCallingMethod  # noqa: E402
from vstolib.cli.cli_annotate import add_cli_annotate_arg_parser  # noqa: E402
from vstolib.cli.cli_collapse import add_cli_collapse_arg_parser
from vstolib.cli.cli_diff import add_cli_diff_arg_parser
from vstolib.cli.cli_intersect import add_cli_intersect_arg_parser
from vstolib.cli.cli_merge import add_cli_merge_arg_parser
from vstolib.cli.cli_overlap import add_cli_overlap_arg_parser
from vstolib.cli.cli_vcf2tsv import add_cli_vcf2tsv_arg_parser


# Match the order in cli_main.run so generated docs reflect the live CLI.
SUBPARSER_REGISTRARS = [
    add_cli_annotate_arg_parser,
    add_cli_collapse_arg_parser,
    add_cli_diff_arg_parser,
    add_cli_intersect_arg_parser,
    add_cli_merge_arg_parser,
    add_cli_overlap_arg_parser,
    add_cli_vcf2tsv_arg_parser,
]


# Resolve paths from the script's own location so the script works regardless
# of the current working directory. Script lives at
# ``<repo_root>/scripts/docs/generate_docs.py``, so the repo root is two
# directories up.
REPO_ROOT = Path(__file__).resolve().parents[2]
DOCS_DIR = REPO_ROOT / "docs"
CLI_DOCS_DIR = DOCS_DIR / "cli"


# Concrete usage example per subcommand. Used to populate the "Example" block
# on each page alongside the auto-generated "Template" block. Mirrors the
# runnable scripts in ``examples/`` (the COLO829 / HCC1395 SAVANA callsets).
EXAMPLES: dict[str, str] = {
    "vcf2tsv": """\
vstol vcf2tsv \\
    --vcf-file colo829-dna-pacbio_savana.somatic.vcf.gz \\
    --method savana \\
    --platform pacbio \\
    --source-id colo829 \\
    --output-tsv-file colo829-dna-pacbio_savana.somatic.tsv""",
    "annotate": """\
vstol annotate \\
    --tsv-file colo829-dna-pacbio_savana.somatic.tsv \\
    --annotator gencode \\
    --gencode-gtf-file gencode.v41.annotation.gtf.gz \\
    --gencode-version v41 \\
    --gencode-species human \\
    --gencode-levels 1 2 \\
    --num-processes 2 \\
    --output-tsv-file colo829-dna-pacbio_savana.somatic.annotated.tsv""",
    "diff": """\
vstol diff \\
    --target-tsv-file colo829-dna-pacbio_savana.somatic.tsv \\
    --query-tsv-file hcc1395-dna-pacbio_savana.somatic.tsv \\
    --match-both-positions yes \\
    --max-breakpoint-distance 1000 \\
    --match-operation-types yes \\
    --output-tsv-file colo829-specific.somatic.tsv""",
    "intersect": """\
vstol intersect \\
    --tsv-files colo829-dna-pacbio_savana.somatic.tsv hcc1395-dna-pacbio_savana.somatic.tsv \\
    --output-tsv-file colo829-hcc1395_intersect.tsv""",
    "merge": """\
vstol merge \\
    --tsv-files colo829-dna-pacbio_savana.somatic.tsv hcc1395-dna-pacbio_savana.somatic.tsv \\
    --output-tsv-file colo829-hcc1395_merge.tsv""",
    "collapse": """\
vstol collapse \\
    --tsv-file colo829-hcc1395_merge.tsv \\
    --tumor-sample-ids colo829-dna-pacbio hcc1395-dna-pacbio \\
    --output-tsv-file colo829-hcc1395_merge.collapsed.tsv""",
    "overlap": """\
vstol overlap \\
    --tsv-file colo829-dna-pacbio_savana.somatic.tsv \\
    --regions-tsv-file hg38_ucsc_gap_table.tsv \\
    --output-tsv-file colo829-dna-pacbio_savana.gap_overlap.tsv""",
}


# I/O contract per subcommand, rendered as an "At a glance" callout under the
# Description on each page. "inputs" lists the file types the command consumes;
# "outputs" lists what it produces; "next_steps" links to the command(s) that
# typically run after this one in a variant-standardization workflow.
IO_CONTRACT: dict[str, dict[str, str]] = {
    "vcf2tsv": {
        "inputs": "`*.vcf` / `*.vcf.gz` (variant-caller output)",
        "outputs": "`*.tsv` (Occam's Variant Grammar)",
        "next_steps": "[`annotate`](annotate.qmd), [`merge`](merge.qmd), [`intersect`](intersect.qmd), [`diff`](diff.qmd), [`overlap`](overlap.qmd)",
    },
    "annotate": {
        "inputs": "`*.tsv` (variants), `*.gtf.gz` (GENCODE gene annotation)",
        "outputs": "`*.tsv` (variants with gene-level annotation)",
        "next_steps": "Endpoint — annotated callset",
    },
    "diff": {
        "inputs": "`*.tsv` (target variants), one or more `*.tsv` (query variants)",
        "outputs": "`*.tsv` (variants specific to the target list)",
        "next_steps": "[`annotate`](annotate.qmd)",
    },
    "intersect": {
        "inputs": "two or more variant `*.tsv`",
        "outputs": "`*.tsv` (variants shared across all input lists)",
        "next_steps": "[`annotate`](annotate.qmd)",
    },
    "merge": {
        "inputs": "two or more variant `*.tsv`",
        "outputs": "`*.tsv` (union of variants, with a `group_id` per cluster)",
        "next_steps": "[`collapse`](collapse.qmd)",
    },
    "collapse": {
        "inputs": "`*.tsv` (a `vstol merge` output with a `group_id` column)",
        "outputs": "`*.tsv` (one row per `group_id`, somatic-filtered)",
        "next_steps": "[`annotate`](annotate.qmd)",
    },
    "overlap": {
        "inputs": "`*.tsv` (variants), `*.tsv` (genomic regions: chromosome/start/end)",
        "outputs": "`*.tsv` (variants overlapping the regions)",
        "next_steps": "Endpoint — region-filtered callset",
    },
}


# Curated metadata for the "Variant Callers" reference page, keyed by the
# canonical ``--method`` value (the same identifiers in
# ``VariantCallingMethod.ALL``). Each value is (display name, homepage URL,
# version tested). The homepages and tested versions are not derivable from the
# source, so they live here; ``render_variant_callers_page`` cross-checks the
# keys against ``VariantCallingMethod.ALL`` so a newly-added caller can't ship
# without a documented version. Rows are emitted in ``ALL`` order.
CALLER_METADATA: dict[str, tuple[str, str, str]] = {
    "clair3": ("Clair3", "https://github.com/HKU-BAL/Clair3", "1.2.0"),
    "clairs": ("ClairS", "https://github.com/HKU-BAL/ClairS", "0.3.0"),
    "cutesv": ("CuteSV", "https://github.com/tjiangHIT/cuteSV", "2.1.1"),
    "dbsnp": ("dbSNP", "https://www.ncbi.nlm.nih.gov/snp/", "146"),
    "deepsomatic": ("DeepSomatic", "https://github.com/google/deepsomatic", "1.9.0"),
    "deepvariant": ("DeepVariant", "https://github.com/google/deepvariant", "1.9.0"),
    "delly2": ("Delly2", "https://github.com/dellytools/delly", "1.2.6"),
    "dysgu": ("Dysgu", "https://github.com/kcleal/dysgu", "1.8.5"),
    "gridss2": ("GRIDSS2", "https://github.com/PapenfussLab/gridss", "2.13.2"),
    "haplotypecaller": (
        "HaplotypeCaller",
        "https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller",
        "4.4.0.0 (GATK4)",
    ),
    "longshot": ("Longshot", "https://github.com/pjedge/longshot", "1.0.0"),
    "lumpy": ("Lumpy", "https://github.com/arq5x/lumpy-sv", "0.3.1"),
    "manta": ("Manta", "https://github.com/Illumina/manta", "1.6.0"),
    "mutect2": (
        "Mutect2",
        "https://gatk.broadinstitute.org/hc/en-us/articles/360035531132--How-to-Call-somatic-mutations-using-GATK4-Mutect2",
        "4.4.0.0 (GATK4)",
    ),
    "nanocaller": ("NanoCaller", "https://github.com/WGLab/NanoCaller", "3.6.2"),
    "nanomonsv": ("Nanomonsv", "https://github.com/friend1ws/nanomonsv", "0.8.0"),
    "nanovar": ("NanoVar", "https://github.com/benoukraflab/NanoVar", "1.8.3"),
    "octopus": ("Octopus", "https://github.com/luntergroup/octopus", "0.7.4"),
    "pbsv": ("PBSV", "https://github.com/PacificBiosciences/pbsv", "2.9.0"),
    "savana": ("Savana", "https://github.com/cortes-ciriano-lab/savana", "1.3.5"),
    "severus": ("Severus", "https://github.com/KolmogorovLab/Severus", "1.1"),
    "sniffles2": ("Sniffles2", "https://github.com/fritzsedlazeck/Sniffles", "2.3.3"),
    "strelka2-germline": ("Strelka2", "https://github.com/Illumina/strelka", "2.9.10"),
    "strelka2-somatic": ("Strelka2", "https://github.com/Illumina/strelka", "2.9.10"),
    "svaba": ("SvABA", "https://github.com/walaj/svaba", "1.2.0"),
    "svim": ("Svim", "https://github.com/eldariont/svim", "2.0.0"),
    "svisionpro": ("SVision-pro", "https://github.com/songbowang125/SVision-pro", "1.8"),
}


def render_variant_callers_page() -> str:
    """Render the standalone "Variant Callers" reference page.

    One row per ``--method`` value (in ``VariantCallingMethod.ALL`` order) with
    the caller's homepage and the version VSTOL was tested against. Warns if any
    method lacks metadata or any metadata key is no longer a valid method, so the
    page stays in lock-step with the CLI.
    """
    documented = set(CALLER_METADATA)
    supported = set(VariantCallingMethod.ALL)
    missing = supported - documented
    extra = documented - supported
    if missing:
        print(f"  WARNING: --method value(s) without caller metadata: {sorted(missing)}")
    if extra:
        print(f"  WARNING: caller metadata key(s) not in VariantCallingMethod.ALL: {sorted(extra)}")

    lines = [
        "---",
        'title: "Variant Callers"',
        'description: "Variant callers supported by vstol vcf2tsv and the versions tested."',
        "---",
        "",
        "`vstol vcf2tsv` converts the VCF output of the variant callers below into the",
        "[Occam's Variant Grammar](https://pirlblog.substack.com/p/variant-standardization)",
        "TSV. Select the caller with the `--method` value shown in the table.",
        "",
        "::: {.callout-note}",
        'The "Version Tested" column records the release each converter was validated',
        "against; other versions of the same caller usually work but are not guaranteed.",
        "Run `vstol vcf2tsv --help` for the authoritative list of accepted `--method`",
        "values.",
        ":::",
        "",
        "## DNA Variant Callers",
        "",
        "| Variant Caller | `--method` | Version Tested |",
        "|:---|:---|:---|",
    ]
    for method in VariantCallingMethod.ALL:
        meta = CALLER_METADATA.get(method)
        if not meta:
            # Surfaced via the warning above; still emit a row so the gap is
            # visible in the rendered page rather than silently dropped.
            lines.append(f"| (undocumented) | `{method}` | — |")
            continue
        display, url, version = meta
        lines.append(f"| [{display}]({url}) | `{method}` | {version} |")
    lines += [
        "",
        ': {.striped .hover tbl-colwidths="[45, 30, 25]"}',
        "",
    ]
    return "\n".join(lines).rstrip() + "\n"


def build_root_parser() -> argparse.ArgumentParser:
    """Mirror ``cli_main.init_arg_parser`` + ``cli_main.run`` registration."""
    parser = argparse.ArgumentParser(
        description="VSTOL: Variant Selection, Tabulation, and Operations Library."
    )
    parser.add_argument(
        "--version", "-v",
        action="version",
        version=f"%(prog)s version {vstolib.__version__}",
    )
    sub_parsers = parser.add_subparsers(help="VSTOL sub-commands.")
    for register in SUBPARSER_REGISTRARS:
        register(sub_parsers=sub_parsers)
    return parser


def iter_subparsers(
    parser: argparse.ArgumentParser,
) -> Iterable[tuple[str, argparse.ArgumentParser, str]]:
    """Yield (name, subparser, top-level help) for every registered subcommand."""
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            help_by_name = {
                choice.dest: (choice.help or "")
                for choice in action._choices_actions
            }
            for name, subparser in action.choices.items():
                yield name, subparser, help_by_name.get(name, "")


def render_type(action: argparse.Action) -> str:
    if isinstance(action, (argparse._StoreTrueAction, argparse._StoreFalseAction)):
        return "flag"
    t = action.type
    if t is None:
        return "str"
    name = getattr(t, "__name__", str(t))
    # vstol uses a `str2bool` converter for yes/no arguments; show the
    # user-facing spelling rather than leaking the internal converter name.
    # Use a slash, not a pipe: the type is rendered inside a backtick code span
    # in a Markdown table, where a literal `|` breaks the cell and an escaped
    # `\|` would show the backslash verbatim (code spans are not unescaped).
    if name == "str2bool":
        return "yes/no"
    return name


def yaml_escape(text: str) -> str:
    """Make a string safe to embed inside a double-quoted YAML scalar."""
    return text.replace("\\", "\\\\").replace('"', '\\"').replace("\n", " ").strip()


_DEFAULT_PAREN_RE = re.compile(r"\s*\(default[^)]*\)", re.IGNORECASE)
_CHOICES_PAREN_RE = re.compile(r"\s*\(choices[^)]*\)", re.IGNORECASE)
_ALLOWED_RE = re.compile(r"\s*Allowed options:.*", re.IGNORECASE)

# Above this many choices, list them in a callout rather than cramming them into
# the table's Type column (e.g. vcf2tsv --method has ~27 supported callers).
MAX_INLINE_CHOICES = 8


def _format_default(value: object) -> str:
    """Compact default-value rendering (lists rendered as comma-joined)."""
    if value is None or value is argparse.SUPPRESS:
        return ""
    if isinstance(value, list):
        if not value:
            return ""
        return "`" + ", ".join(str(v) for v in value) + "`"
    if value == "":
        return ""
    return f"`{value}`"


def _clean_help(help_text: str | None) -> str:
    """Strip redundant '(default: ...)' / '(choices: ...)' / 'Allowed options: ...'."""
    if not help_text:
        return ""
    text = help_text.strip().replace("\n", " ")
    text = _DEFAULT_PAREN_RE.sub("", text)
    text = _CHOICES_PAREN_RE.sub("", text)
    text = _ALLOWED_RE.sub("", text)
    # Heal artifacts left by stripping the parenthetical: " ." → "." and any
    # collapsed double-spaces.
    text = re.sub(r"\s+([.,;:])", r"\1", text)
    text = re.sub(r"\s+", " ", text)
    text = text.strip().rstrip(".").rstrip()
    return text + "." if text else ""


def render_args_table(
    actions: list[argparse.Action], include_default: bool
) -> str:
    """Render argparse actions as a Quarto markdown table.

    Choices are folded into the Type column as ``str (a|b)`` so the table stays
    compact. Required-arg tables omit Default; optional-arg tables include it.
    """
    headers = ["Flag", "Type"]
    if include_default:
        headers.append("Default")
    headers.append("Description")

    rows = ["| " + " | ".join(headers) + " |"]
    rows.append("|" + "|".join([":---"] * len(headers)) + "|")

    for action in actions:
        flag = ", ".join(f"`{s}`" for s in action.option_strings) or f"`{action.dest}`"
        # `type_str` is rendered inside a backtick code span, so it must contain
        # no `|` (render_type guarantees this). Choices below carry their own
        # escaped pipes in plain text, where Pandoc does unescape them.
        type_str = render_type(action)
        if action.choices and len(action.choices) <= MAX_INLINE_CHOICES:
            choices_str = "\\|".join(str(c) for c in action.choices)
            type_cell = f"`{type_str}` ({choices_str})"
        else:
            type_cell = f"`{type_str}`"
        desc_cell = _clean_help(action.help).replace("|", "\\|")

        row = [flag, type_cell]
        if include_default:
            row.append(_format_default(action.default))
        row.append(desc_cell)
        rows.append("| " + " | ".join(row) + " |")

    # Column widths sum to 100. Required tables (no Default) give Description
    # more room.
    if include_default:
        colwidths = "[28, 14, 12, 46]"
    else:
        colwidths = "[30, 16, 54]"

    rows.append("")
    rows.append(f': {{.striped .hover tbl-colwidths="{colwidths}"}}')
    return "\n".join(rows)


def render_io_contract(name: str) -> str:
    """Render the per-subcommand I/O contract as a Quarto callout block."""
    contract = IO_CONTRACT.get(name)
    if not contract:
        return ""
    return "\n".join([
        '::: {.callout-note title="At a glance"}',
        f"**Inputs:** {contract['inputs']}",
        "",
        f"**Outputs:** {contract['outputs']}",
        "",
        f"**Typical next step:** {contract['next_steps']}",
        ":::",
        "",
    ])


def render_usage(name: str, subparser: argparse.ArgumentParser) -> str:
    """Build a multi-line ``vstol <name> ...`` snippet listing required args
    (as ``--flag <placeholder>``) and optional args (as ``[--flag PLACEHOLDER]``,
    matching argparse's default usage style)."""
    actions = [
        a for a in subparser._actions
        if not isinstance(a, argparse._HelpAction)
    ]
    if not actions:
        return f"vstol {name}"

    parts = [f"vstol {name}"]

    # Required first, then optional — same convention argparse uses.
    for action in actions:
        if not action.required:
            continue
        flag = action.option_strings[0] if action.option_strings else action.dest
        if isinstance(action, (argparse._StoreTrueAction, argparse._StoreFalseAction)):
            placeholder = ""
        elif action.choices and len(action.choices) <= MAX_INLINE_CHOICES:
            # Short choice sets (e.g. yes|no) stay inline; large ones (e.g.
            # vcf2tsv --method's ~27 callers) would bloat the template, so fall
            # back to a generic placeholder — the full list lives in the callout.
            placeholder = "<" + "|".join(str(c) for c in action.choices) + ">"
        else:
            placeholder = f"<{action.dest}>"
        parts.append(f"    {flag} {placeholder}".rstrip())

    for action in actions:
        if action.required:
            continue
        flag = action.option_strings[0] if action.option_strings else action.dest
        if isinstance(action, (argparse._StoreTrueAction, argparse._StoreFalseAction)):
            chunk = f"[{flag}]"
        else:
            placeholder = action.dest.upper()
            if action.nargs in ("*", "+"):
                placeholder = f"{placeholder} [{placeholder} ...]"
            chunk = f"[{flag} {placeholder}]"
        parts.append(f"    {chunk}")

    return " \\\n".join(parts)


def render_subcommand(
    name: str, subparser: argparse.ArgumentParser, help_text: str
) -> str:
    description = (subparser.description or help_text or "").strip()
    lines = [
        "---",
        f'title: "{name}"',
        f'description: "{yaml_escape(description)}"',
        "---",
        "",
        "## Usage",
        "",
        "Template:",
        "",
        "```bash",
        render_usage(name, subparser),
        "```",
        "",
    ]

    example = EXAMPLES.get(name)
    if example:
        lines += [
            "Example:",
            "",
            "```bash",
            example,
            "```",
            "",
        ]

    if description:
        lines += ["## Description", "", description, ""]

    contract_block = render_io_contract(name)
    if contract_block:
        lines.append(contract_block)

    # Args with many choices (e.g. vcf2tsv --method) would overflow the Type
    # column, so surface their allowed values in a callout — keeping the full
    # list readily visible.
    for action in subparser._actions:
        choices = getattr(action, "choices", None)
        if choices and len(choices) > MAX_INLINE_CHOICES:
            flag = action.option_strings[0] if action.option_strings else action.dest
            chips = ", ".join(f"`{c}`" for c in choices)
            lines += [
                f'::: {{.callout-tip title="Supported {flag} values"}}',
                chips,
                ":::",
                "",
            ]

    for group in subparser._action_groups:
        actions = [
            a for a in group._group_actions
            if not isinstance(a, argparse._HelpAction)
        ]
        if not actions:
            continue
        title = (group.title or "Arguments").strip()
        title = title[:1].upper() + title[1:]
        all_required = all(a.required for a in actions)
        lines += [
            f"## {title}",
            "",
            render_args_table(actions, include_default=not all_required),
            "",
        ]

    return "\n".join(lines).rstrip() + "\n"


def render_index(
    subcommands: list[tuple[str, argparse.ArgumentParser, str]],
) -> str:
    lines = [
        "---",
        'title: "Commands"',
        "---",
        "",
        f"`vstol` version {vstolib.__version__}.",
        "",
        "## Subcommands",
        "",
    ]
    for name, _subparser, help_text in sorted(subcommands, key=lambda t: t[0]):
        suffix = f" — {help_text.strip()}" if help_text.strip() else ""
        lines.append(f"- [`{name}`]({name}.qmd){suffix}")
    return "\n".join(lines) + "\n"


def main() -> int:
    print("Generating Quarto CLI documentation...")
    CLI_DOCS_DIR.mkdir(parents=True, exist_ok=True)

    parser = build_root_parser()
    subcommands = list(iter_subparsers(parser))

    for name, subparser, help_text in subcommands:
        out_path = CLI_DOCS_DIR / f"{name}.qmd"
        out_path.write_text(render_subcommand(name, subparser, help_text))
        print(f"  wrote {out_path.relative_to(REPO_ROOT)}")

    index_path = CLI_DOCS_DIR / "index.qmd"
    if index_path.exists():
        print(f"  kept {index_path.relative_to(REPO_ROOT)} (user-editable; not overwritten)")
    else:
        index_path.write_text(render_index(subcommands))
        print(f"  wrote {index_path.relative_to(REPO_ROOT)} (initial scaffold; safe to edit)")

    # Standalone reference page (curated versions); regenerated every run.
    callers_path = DOCS_DIR / "variant-callers.qmd"
    callers_path.write_text(render_variant_callers_page())
    print(f"  wrote {callers_path.relative_to(REPO_ROOT)}")

    print(f"Done! {len(subcommands)} subcommand pages generated.")
    print(f"  Preview with: cd {DOCS_DIR.relative_to(REPO_ROOT)} && quarto preview")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
