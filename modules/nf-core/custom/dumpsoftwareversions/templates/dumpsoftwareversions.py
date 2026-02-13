#!/usr/bin/env python

"""Provide functions to merge multiple versions.yml files."""

import platform
from textwrap import dedent

import yaml


def _normalize_versions(value):
    """
    Ensure each process' versions value is a dict of {tool: version}.
    Accepts dict, string, list, None, etc and returns a dict.
    """
    if value is None:
        return {}

    if isinstance(value, dict):
        return value

    if isinstance(value, list):
        merged = {}
        for item in value:
            if isinstance(item, dict):
                merged.update(item)
            else:
                merged[str(item)] = "unknown"
        return merged

    if isinstance(value, str):
        try:
            parsed = yaml.safe_load(value)
        except Exception:
            parsed = None

        if isinstance(parsed, dict):
            return parsed
        if isinstance(parsed, list):
            merged = {}
            for item in parsed:
                if isinstance(item, dict):
                    merged.update(item)
                else:
                    merged[str(item)] = "unknown"
            return merged

        return {"_raw": value}

    return {"_raw": str(value)}


def _make_versions_html(versions):
    """Generate a tabular HTML output of all versions for MultiQC."""
    html = [
        dedent(
            """\
            <style>
            #nf-core-versions tbody:nth-child(even) {
                background-color: #f2f2f2;
            }
            </style>
            <table class="table" style="width:100%" id="nf-core-versions">
                <thead>
                    <tr>
                        <th> Process Name </th>
                        <th> Software </th>
                        <th> Version  </th>
                    </tr>
                </thead>
            """
        )
    ]

    for process, tmp_versions in sorted(versions.items()):
        tmp_versions = _normalize_versions(tmp_versions)

        html.append("<tbody>")
        for i, (tool, version) in enumerate(sorted(tmp_versions.items())):
            html.append(
                dedent(
                    f"""\
                    <tr>
                        <td><samp>{process if (i == 0) else ''}</samp></td>
                        <td><samp>{tool}</samp></td>
                        <td><samp>{version}</samp></td>
                    </tr>
                    """
                )
            )
        html.append("</tbody>")

    html.append("</table>")
    return chr(10).join(html)


def main():
    """Load all version files and generate merged output."""
    versions_this_module = {}
    versions_this_module["${task.process}"] = {
        "python": platform.python_version(),
        "yaml": yaml.__version__,
    }

    with open("$versions") as f:
        loaded = yaml.load(f, Loader=yaml.BaseLoader) or {}

    versions_by_process = {k: _normalize_versions(v) for k, v in loaded.items()}
    versions_by_process.update(versions_this_module)

    versions_by_module = {}
    for process, process_versions in versions_by_process.items():
        module = process.split(":")[-1]
        process_versions = _normalize_versions(process_versions)

        if module not in versions_by_module:
            versions_by_module[module] = dict(process_versions)
        else:
            merged = dict(versions_by_module[module])
            for tool, ver in process_versions.items():
                if tool not in merged:
                    merged[tool] = ver
                elif merged[tool] != ver:
                    merged[f"{tool} (alt)"] = ver
            versions_by_module[module] = merged

    versions_by_module["Workflow"] = {
        "Nextflow": "$workflow.nextflow.version",
        "$workflow.manifest.name": "$workflow.manifest.version",
    }

    versions_mqc = {
        "id": "software_versions",
        "section_name": "${workflow.manifest.name} Software Versions",
        "section_href": "https://github.com/${workflow.manifest.name}",
        "plot_type": "html",
        "description": "are collected at run time from the software output.",
        "data": _make_versions_html(versions_by_module),
    }

    with open("software_versions.yml", "w") as f:
        yaml.dump(versions_by_module, f, default_flow_style=False)

    with open("software_versions_mqc.yml", "w") as f:
        yaml.dump(versions_mqc, f, default_flow_style=False)

    with open("versions.yml", "w") as f:
        yaml.dump(versions_this_module, f, default_flow_style=False)


if __name__ == "__main__":
    main()

